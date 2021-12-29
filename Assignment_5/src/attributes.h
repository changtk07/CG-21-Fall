#pragma once

#include <fstream>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Dense>

struct Light {
	Eigen::Vector4f position;
	Eigen::Vector4f intensity;
};

struct Material {
	Eigen::Vector4f ambient_color;
	Eigen::Vector4f diffuse_color;
	Eigen::Vector4f specular_color;
	float specular_exponent;
};

struct Camera {
	bool is_perspective;
	float field_of_view;
	Eigen::Vector4f position;
	Eigen::Vector4f gaze;
	Eigen::Vector4f min;
	Eigen::Vector4f max;
};

class VertexAttributes {
public:
	VertexAttributes(float x = 0, float y = 0, float z = 0, float w = 1, int a = 1) : normal(0,0,0,0) {
		position << x,y,z,w;
	}

	// Interpolates the vertex attributes
	static VertexAttributes interpolate(
		const VertexAttributes& a,
		const VertexAttributes& b,
		const VertexAttributes& c,
		const float alpha,
		const float beta,
		const float gamma
	) 
	{
		VertexAttributes r;
		r.position = alpha*(a.position/a.position[3]) + beta*(b.position/b.position[3]) + gamma*(c.position/c.position[3]);
		r.color = alpha*a.color + beta*b.color + gamma*c.color;
		return r;
	}

	int a;
	Eigen::Vector4f position;
	Eigen::Vector4f normal;
	Eigen::Vector4f color;
};

class FragmentAttributes {
public:
	FragmentAttributes(float r = 0, float g = 0, float b = 0, float a = 1) {
		color << r,g,b,a;
	}

	float depth;
	Eigen::Vector4f color;
};

class FrameBufferAttributes {
public:
	FrameBufferAttributes(uint8_t r = 0, uint8_t g = 0, uint8_t b = 0, uint8_t a = 255) : depth(-999) {
		color << r,g,b,a;
	}

	Eigen::Matrix<uint8_t,4,1> color;
	float depth;
};

class UniformAttributes {
public:
	Eigen::Matrix4f M1, M2;

	Eigen::Vector4f ambient_light;
	Light light;
	Material mat;
	Camera camera;
};

class Mesh {
public:
	Eigen::MatrixXf V; // n x 4 matrix (n points)
	Eigen::MatrixXi F; // m x 3 matrix (m triangles)

	void load_off(const std::string filename) {
		std::ifstream in(filename);
		std::string token;
		in >> token;
		int nv, nf, ne;
		in >> nv >> nf >> ne;
		V.resize(nv, 4);
		F.resize(nf, 3);
		for (int i = 0; i < nv; ++i) {
			in >> V(i, 0) >> V(i, 1) >> V(i, 2);
			V(i, 3) = 1;
		}
		for (int i = 0; i < nf; ++i) {
			int s;
			in >> s >> F(i, 0) >> F(i, 1) >> F(i, 2);
			assert(s == 3);
		}
	}

	void build_vas_triangles(std::vector<VertexAttributes> &vas) {
		vas.resize(F.rows()*3);
		for (int i = 0; i < F.rows(); ++i) {
			vas[i*3+0] = VertexAttributes(V(F(i,0), 0), V(F(i,0), 1), V(F(i,0), 2));
			vas[i*3+1] = VertexAttributes(V(F(i,1), 0), V(F(i,1), 1), V(F(i,1), 2));
			vas[i*3+2] = VertexAttributes(V(F(i,2), 0), V(F(i,2), 1), V(F(i,2), 2));
		}
	}

	void build_vas_lines(std::vector<VertexAttributes>& vas) {
		vas.resize(F.rows()*6);
		for (int i = 0; i < F.rows(); ++i) {
			vas[i*6+0] = VertexAttributes(V(F(i,0), 0), V(F(i,0), 1), V(F(i,0), 2));
			vas[i*6+1] = VertexAttributes(V(F(i,1), 0), V(F(i,1), 1), V(F(i,1), 2));
			vas[i*6+2] = VertexAttributes(V(F(i,1), 0), V(F(i,1), 1), V(F(i,1), 2));
			vas[i*6+3] = VertexAttributes(V(F(i,2), 0), V(F(i,2), 1), V(F(i,2), 2));
			vas[i*6+4] = VertexAttributes(V(F(i,2), 0), V(F(i,2), 1), V(F(i,2), 2));
			vas[i*6+5] = VertexAttributes(V(F(i,0), 0), V(F(i,0), 1), V(F(i,0), 2));
		}
	}

	void per_face_normal(const Camera &camera, std::vector<VertexAttributes> &vas) {
		for (int i = 0; i < F.rows(); ++i) {
			VertexAttributes &a = vas[i*3+0];
			VertexAttributes &b = vas[i*3+1];
			VertexAttributes &c = vas[i*3+2];

			Eigen::Vector3f apos(a.position[0], a.position[1], a.position[2]);
			Eigen::Vector3f bpos(b.position[0], b.position[1], b.position[2]);
			Eigen::Vector3f cpos(c.position[0], c.position[1], c.position[2]);
			Eigen::Vector3f n = (bpos-apos).cross(cpos-apos).normalized();

			Eigen::Vector4f face_normal(n[0], n[1], n[2], 0);
			face_normal = face_normal.dot(-camera.gaze) > 0 ? face_normal : -face_normal;

			a.normal = face_normal;
			b.normal = face_normal;
			c.normal = face_normal;
		}
	}

	void per_vertex_normal(const Camera &camera, std::vector<VertexAttributes> &vas) {
		std::vector<int> count(V.rows(), 0);
		std::vector<Eigen::Vector4f> normals(V.rows(), Eigen::Vector4f(0,0,0,0));

		for (int i = 0; i < F.rows(); ++i) {
			Eigen::Vector3f apos(V(F(i,0), 0), V(F(i,0), 1), V(F(i,0), 2));
			Eigen::Vector3f bpos(V(F(i,1), 0), V(F(i,1), 1), V(F(i,1), 2));
			Eigen::Vector3f cpos(V(F(i,2), 0), V(F(i,2), 1), V(F(i,2), 2));
			Eigen::Vector3f n = (bpos-apos).cross(cpos-apos).normalized();

			Eigen::Vector4f face_normal(n[0], n[1], n[2], 0);
			Eigen::Vector4f view = camera.is_perspective ? (camera.position-vas[i*3].position).normalized() : -camera.gaze;
			face_normal = face_normal.dot(view) > 0 ? face_normal : -face_normal;

			normals[F(i,0)] += face_normal;
			normals[F(i,1)] += face_normal;
			normals[F(i,2)] += face_normal;
			++count[F(i,0)];
			++count[F(i,1)];
			++count[F(i,2)];
		}

		for (int i = 0; i < F.rows(); ++i) {
			vas[i*3+0].normal = normals[F(i,0)] / count[F(i,0)];
			vas[i*3+1].normal = normals[F(i,1)] / count[F(i,1)];
			vas[i*3+2].normal = normals[F(i,2)] / count[F(i,2)];
		}
	}

	Eigen::AlignedBox4f bbox(const Eigen::Matrix4f &M) const {
		Eigen::AlignedBox4f bbox1;
		for (int i = 0; i < V.rows(); ++i) {
			Eigen::Vector4f v = V.row(i);
			bbox1.extend(v);
		}
		Eigen::AlignedBox4f bbox2;
		Eigen::Vector4f max = bbox1.max();
		Eigen::Vector4f min = bbox1.min();
		bbox2.extend(M*max);
		bbox2.extend(M*min);
		return bbox2;
	}

	Eigen::Vector4f barycenter() const {
		Eigen::Vector4f barycenter(0,0,0,0);
		for (int i = 0; i < V.rows(); ++i) {
			barycenter += V.row(i);
		}
		barycenter /= V.rows();
		return barycenter;
	}
};