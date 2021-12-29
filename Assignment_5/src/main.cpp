// C++ include
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <climits>

// Utilities for the Assignment
#include <gif.h>
#include "raster.h"

// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"

using namespace std;

Eigen::Matrix4f camera_matrix(const Camera &camera) {
	Eigen::Vector4f w = -camera.gaze.normalized();
	Eigen::Vector3f w3d(w[0], w[1], w[2]);
	Eigen::Vector3f t3d(0, 1, 0);
	Eigen::Vector3f u3d = t3d.cross(w3d).normalized();
	Eigen::Vector3f v3d = w3d.cross(u3d);

	Eigen::Matrix4f M;
	M.col(0) << u3d[0], u3d[1], u3d[2], 0;
	M.col(1) << v3d[0], v3d[1], v3d[2], 0;
	M.col(2) = w;
	M.col(3) = camera.position;
	return M.inverse();
}

Eigen::Matrix4f project_matrix(const Camera &camera) {
	Eigen::Matrix4f M, P;
	float l = camera.min[0], b = camera.min[1], f = camera.min[2];
	float r = camera.max[0], t = camera.max[1], n = camera.max[2];
	M <<
	2/(r-l), 0, 0, -(r+l)/(r-l),
	0, 2/(t-b), 0, -(t+b)/(t-b),
	0, 0, 2/(n-f), -(n+f)/(n-f),
	0, 0, 0, 1;
	if (camera.is_perspective) {
		P <<
		n, 0, 0, 0,
		0, n, 0,  0,
		0, 0, n+f, -n*f,
		0, 0, 1, 0;
	}
	else {
		P = Eigen::Matrix4f::Identity();
	}
	return M*P;
}

Eigen::Vector4f raytrace(const VertexAttributes& va, const UniformAttributes& uniform) {
	Eigen::Vector4f lights_color(0,0,0,1);
	Eigen::Vector4f ambient = uniform.mat.ambient_color.array() * uniform.ambient_light.array();
	lights_color += ambient;

	Eigen::Vector4f l = uniform.light.position - va.position;
	Eigen::Vector4f v = uniform.camera.is_perspective ? (uniform.camera.position-va.position).normalized() : -uniform.camera.gaze;
	Eigen::Vector4f h = (v+l).normalized();
	Eigen::Vector4f diffuse = uniform.mat.diffuse_color * std::max(l.dot(va.normal), (float)0.0);
	Eigen::Vector4f specular = uniform.mat.specular_color * std::pow(std::max(h.dot(va.normal), (float)0.0), uniform.mat.specular_exponent);
	lights_color += (diffuse+specular).cwiseProduct(uniform.light.intensity) / l.squaredNorm();

	// clamp to 0 and 1 to avoid underflow/overflow
	lights_color = lights_color.cwiseMax(0);
	lights_color = lights_color.cwiseMin(1);
	return lights_color;
}

int main(int argc, char* argv[]) {
	std::string render_mode = "perface"; // wireframe, perface, pervertex
	bool render_gif = true;

	// The Framebuffer storing the image rendered by the rasterizer
	Eigen::Matrix<FrameBufferAttributes,Eigen::Dynamic,Eigen::Dynamic> frameBuffer(800,450);

	// Basic rasterization program
	Program program;

	program.VertexShader = [](const VertexAttributes& va, const UniformAttributes& uniform) {
		VertexAttributes out;
		out.position = uniform.M1 * va.position;
		out.normal = uniform.M1.transpose().inverse() * va.normal;
		if (va.normal == Eigen::Vector4f(0,0,0,0)) // lines
			out.color = Eigen::Vector4f(0.6, 0.2, 0.2, 1);
		else // triangles
			out.color = raytrace(out, uniform);
		out.position = uniform.M2 * out.position;
		return out;
	};

	program.FragmentShader = [](const VertexAttributes& va, const UniformAttributes& uniform) {
		FragmentAttributes out;
		out.color = va.color;
		out.depth = va.position[2];
		return out;
	};

	program.BlendingShader = [](const FragmentAttributes& fa, const FrameBufferAttributes& previous) {
		if (fa.depth + 1e-4 > previous.depth) {
			FrameBufferAttributes out(fa.color[0]*255, fa.color[1]*255, fa.color[2]*255, fa.color[3]*255);
			out.depth = fa.depth;
			return out;
		}
		else
			return previous;
	};

	Mesh mesh;
	mesh.load_off(std::string(DATA_DIR)+"bunny.off");

	// Uniform attributes
	UniformAttributes uniform;
	uniform.ambient_light << 0.2, 0.2, 0.2, 0.0;

	// Lights
	uniform.light.position << 0, 7, 7, 0;
	uniform.light.intensity << 16, 16, 16, 0;

	// Material
	uniform.mat.ambient_color << 0.0, 0.5, 0.0, 0.0;
	uniform.mat.diffuse_color << 0.5, 0.5, 0.5, 0.0;
	uniform.mat.specular_color << 0.2, 0.2, 0.2, 0.0;
	uniform.mat.specular_exponent = 256.0;

	// Camera
	uniform.camera.is_perspective = true;
	uniform.camera.field_of_view = 0.2;
	uniform.camera.position << -0.02, 0.1, 1, 1;
	uniform.camera.gaze << 0, 0, -1, 0;

	// Matrix
	// to camera space
	Eigen::Matrix4f M_cam = camera_matrix(uniform.camera);

	// to canonical
	if (uniform.camera.is_perspective) {
		float n = -0.1, f = -2.5;
		float t = std::tan(uniform.camera.field_of_view/2) * std::abs(n);
		uniform.camera.max << t, t, n, 1;
		uniform.camera.min << -t, -t, f, 1;
	}
	else {
		Eigen::AlignedBox4f bbox = mesh.bbox(M_cam);
		uniform.camera.max = bbox.max(); // r, t, n
		uniform.camera.min = bbox.min(); // l, b, f
		uniform.camera.max = uniform.camera.max.array() + 0.01;
		uniform.camera.min = uniform.camera.min.array() - 0.01;
	}
	Eigen::Matrix4f M_proj = project_matrix(uniform.camera);

	// view
	Eigen::Matrix4f M_view = Eigen::Matrix4f::Identity();
	float aspect_ratio = (float)frameBuffer.cols() / frameBuffer.rows();
	if (aspect_ratio < 1)
		M_view(0,0) = aspect_ratio;
	else
		M_view(1,1) = 1/aspect_ratio;

	// set parameters for translation and rotation
	int frames = 1, delay = 10;
	float translate_start = 0, translate_step = 0;
	float rotation_start = 0, rotation_step = 0;
	if (render_gif) {
		frames = 128;
		translate_start = -1;
		translate_step = 0.01;
		rotation_start = 0;
		rotation_step = -2 * M_PI / frames;
	}

	// prepare vertex attributes
	vector<VertexAttributes> vas_triangles;
	vector<VertexAttributes> vas_lines;
	if (render_mode == "perface") {
		mesh.build_vas_triangles(vas_triangles);
		mesh.per_face_normal(uniform.camera, vas_triangles);
	}
	if (render_mode == "perface" || render_mode == "wireframe") {
		mesh.build_vas_lines(vas_lines);
	}
	if (render_mode == "pervertex") {
		mesh.build_vas_triangles(vas_triangles);
		mesh.per_vertex_normal(uniform.camera, vas_triangles);
	}

	// Compute Barycenter
	Eigen::Vector4f barycenter = mesh.barycenter();

	// Render
	GifWriter g;
	if (render_gif) {
		GifBegin(&g, "bunny.gif", frameBuffer.rows(), frameBuffer.cols(), delay);
	}

	for (int i = 0; i < frames; ++i) {
		frameBuffer.setConstant(FrameBufferAttributes());
		float t = translate_start + i*translate_step;
		float r = rotation_start + i*rotation_step;

		// Translate & Rotate
		Eigen::Matrix4f M_t1, M_r1, M_t2;
		M_t1 <<
		1, 0, 0, -barycenter[0],
		0, 1, 0, -barycenter[1],
		0, 0, 1, -barycenter[2],
		0, 0, 0, 1;
		M_r1 <<
		std::cos(r), -std::sin(r), 0, 0,
		std::sin(r), std::cos(r), 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1;
		M_t2 <<
		1, 0, 0, barycenter[0],
		0, 1, 0, barycenter[1],
		0, 0, 1, barycenter[2]+t,
		0, 0, 0, 1;


		uniform.M1 = M_view * M_cam * M_t2 * M_r1 * M_t1;
		uniform.M2 = M_proj;

		rasterize_triangles(program, uniform, vas_triangles, frameBuffer);
		rasterize_lines(program, uniform, vas_lines, 0.5, frameBuffer);

		vector<uint8_t> image;
		framebuffer_to_uint8(frameBuffer,image);
		if (render_gif) {
			GifWriteFrame(&g, image.data(), frameBuffer.rows(), frameBuffer.cols(), delay);
		}
		else {
			stbi_write_png("bunny.png", frameBuffer.rows(), frameBuffer.cols(), 4, image.data(), frameBuffer.rows()*4);
		}
	}
	
	return 0;
}
