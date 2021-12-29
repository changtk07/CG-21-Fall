// C++ include
#include <iostream>
#include <string>
#include <vector>
#include <cmath>

// Utilities for the Assignment
#include "utils.h"

// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"

// Shortcut to avoid Eigen:: everywhere, DO NOT USE IN .h
using namespace Eigen;

void raytrace_sphere() {
	std::cout << "Simple ray tracer, one sphere with orthographic projection" << std::endl;

	const std::string filename("sphere_orthographic.png");
	MatrixXd C = MatrixXd::Zero(800,800); // Store the color
	MatrixXd A = MatrixXd::Zero(800,800); // Store the alpha mask

	// The camera is orthographic, pointing in the direction -z and covering the unit square (-1,1) in x and y
	Vector3d origin(-1,1,1);
	Vector3d x_displacement(2.0/C.cols(),0,0);
	Vector3d y_displacement(0,-2.0/C.rows(),0);

	// Single light source
	const Vector3d light_position(-1,1,1);

	for (unsigned i=0; i < C.cols(); ++i) {
		for (unsigned j=0; j < C.rows(); ++j) {
			// Prepare the ray
			Vector3d ray_origin = origin + double(i)*x_displacement + double(j)*y_displacement;
			Vector3d ray_direction = RowVector3d(0,0,-1);

			// Intersect with the sphere
			// NOTE: this is a special case of a sphere centered in the origin and for orthographic rays aligned with the z axis
			Vector2d ray_on_xy(ray_origin(0),ray_origin(1));
			const double sphere_radius = 0.9;

			if (ray_on_xy.norm() < sphere_radius) {
				// The ray hit the sphere, compute the exact intersection point
				Vector3d ray_intersection(ray_on_xy(0),ray_on_xy(1),sqrt(sphere_radius*sphere_radius - ray_on_xy.squaredNorm()));

				// Compute normal at the intersection point
				Vector3d ray_normal = ray_intersection.normalized();

				// Simple diffuse model
				C(i,j) = (light_position-ray_intersection).normalized().transpose() * ray_normal;

				// Clamp to zero
				C(i,j) = std::max(C(i,j),0.);

				// Disable the alpha mask for this pixel
				A(i,j) = 1;
			}
		}
	}

	// Save to png
	write_matrix_to_png(C,C,C,A,filename);

}

void raytrace_parallelogram() {
	std::cout << "Simple ray tracer, one parallelogram with orthographic projection" << std::endl;

	const std::string filename("plane_orthographic.png");
	MatrixXd C = MatrixXd::Zero(800,800); // Store the color
	MatrixXd A = MatrixXd::Zero(800,800); // Store the alpha mask

	// The camera is orthographic, pointing in the direction -z and covering the unit square (-1,1) in x and y
	Vector3d origin(-1,1,1);
	Vector3d x_displacement(2.0/C.cols(),0,0);
	Vector3d y_displacement(0,-2.0/C.rows(),0);

	// Parameters of the parallelogram (position of the lower-left corner + two sides)
	Vector3d pgram_origin(-0.5, -0.5, 0);
 	Vector3d pgram_u(-0.3, 1, -0.5);
 	Vector3d pgram_v(1, -0.3, -1);

	// Single light source
	const Vector3d light_position(-1,1,1);

	for (unsigned i=0; i < C.cols(); ++i) {
		for (unsigned j=0; j < C.rows(); ++j) {
			// Prepare the ray
			Vector3d ray_origin = origin + double(i)*x_displacement + double(j)*y_displacement;
			Vector3d ray_direction = RowVector3d(0,0,-1);

			// Check if the ray intersects with the parallelogram
			// project everything to xy-plane, this is a special case where ray directon is always [0 0 -1]
			Vector2d u(pgram_u(0), pgram_u(1));
			Vector2d v(pgram_v(0), pgram_v(1));
			Vector2d p(ray_origin(0)-pgram_origin(0), ray_origin(1)-pgram_origin(1));
			Matrix2d uv, pu, pv;
			uv.col(0) = u; uv.col(1) = v;
			pu.col(0) = p; pu.col(1) = u;
			pv.col(0) = p; pv.col(1) = v;
			double m = pv.determinant() / uv.determinant();
			double n = -pu.determinant() / uv.determinant();

			if (0 <= m && m <= 1 && 0 <= n && n <= 1) {
				// The ray hit the parallelogram, compute the exact intersection point
				Vector3d ray_intersection = pgram_origin + m*pgram_u + n*pgram_v;

				// Compute normal at the intersection point
				Vector3d plane_normal = pgram_v.cross(pgram_u);
				Vector3d ray_normal = (ray_intersection + plane_normal).normalized();

				// Simple diffuse model
				C(i,j) = (light_position-ray_intersection).normalized().transpose() * ray_normal;

				// Clamp to zero
				C(i,j) = std::max(C(i,j),0.);

				// Disable the alpha mask for this pixel
				A(i,j) = 1;
			}
		}
	}

	// Save to png
	write_matrix_to_png(C,C,C,A,filename);
}

void raytrace_perspective() {
	std::cout << "Simple ray tracer, one parallelogram with perspective projection" << std::endl;

	const std::string filename("plane_perspective.png");
	MatrixXd C = MatrixXd::Zero(800,800); // Store the color
	MatrixXd A = MatrixXd::Zero(800,800); // Store the alpha mask

	// The camera is perspective, pointing in the direction -z and covering the unit square (-1,1) in x and y
	Vector3d camera(0,0,3);
	Vector3d origin(-1,1,1);
	Vector3d x_displacement(2.0/C.cols(),0,0);
	Vector3d y_displacement(0,-2.0/C.rows(),0);

	// Parameters of the parallelogram (position of the lower-left corner + two sides)
	Vector3d pgram_origin(-0.5,-1,0);
	Vector3d pgram_u(-0.5,2,0);
	Vector3d pgram_v(1.5,0,0);

	// Single light source
	const Vector3d light_position(-1,1,1);

	for (unsigned i=0; i < C.cols(); ++i) {
		for (unsigned j=0; j < C.rows(); ++j) {
			// Prepare the ray (origin point and direction)
			Vector3d ray_origin = origin + double(i)*x_displacement + double(j)*y_displacement;
			Vector3d ray_direction = (ray_origin - camera).normalized();

			// Check if the ray intersects with the parallelogram
			Vector3d plane_normal = pgram_v.cross(pgram_u);
			double t = (pgram_origin - ray_origin).dot(plane_normal) / ray_direction.dot(plane_normal);
			Vector3d plane_intersection = ray_origin + t*ray_direction;

			// project everything to xy-plane
			Vector2d u(pgram_u(0), pgram_u(1));
			Vector2d v(pgram_v(0), pgram_v(1));
			Vector2d p(plane_intersection(0)-pgram_origin(0), plane_intersection(1)-pgram_origin(1));
			Matrix2d uv, pu, pv;
			uv.col(0) = u; uv.col(1) = v;
			pu.col(0) = p; pu.col(1) = u;
			pv.col(0) = p; pv.col(1) = v;
			double m = pv.determinant() / uv.determinant();
			double n = -pu.determinant() / uv.determinant();

			if (0 <= m && m <= 1 && 0 <= n && n <= 1) {
				// The ray hit the parallelogram, compute the exact intersection point
				Vector3d ray_intersection = plane_intersection;

				// Compute normal at the intersection point
				Vector3d ray_normal = (ray_intersection + plane_normal).normalized();

				// Simple diffuse model
				C(i,j) = (light_position-ray_intersection).normalized().transpose() * ray_normal;

				// Clamp to zero
				C(i,j) = std::max(C(i,j),0.);

				// Disable the alpha mask for this pixel
				A(i,j) = 1;
			}
		}
	}

	// Save to png
	write_matrix_to_png(C,C,C,A,filename);
}

void raytrace_shading(){
	std::cout << "Simple ray tracer, one sphere with different shading" << std::endl;

	const std::string filename("shading");
	MatrixXd Z = MatrixXd::Zero(800,800);
	MatrixXd C = MatrixXd::Zero(800,800); // Store the color
	MatrixXd A = MatrixXd::Zero(800,800); // Store the alpha mask

	// The camera is perspective, pointing in the direction -z and covering the unit square (-1,1) in x and y
	Vector3d camera(0,0,10);
	Vector3d origin(-1,1,1);
	Vector3d x_displacement(2.0/C.cols(),0,0);
	Vector3d y_displacement(0,-2.0/C.rows(),0);

	// Single light source
	const Vector3d light_position(-1,1,1);
	double ambient = 0.05;
	MatrixXd diffuse = MatrixXd::Zero(800, 800);
	MatrixXd specular = MatrixXd::Zero(800, 800);

	for (unsigned i=0; i < C.cols(); ++i) {
		for (unsigned j=0; j < C.rows(); ++j) {
			// Prepare the ray
			Vector3d ray_origin = origin + double(i)*x_displacement + double(j)*y_displacement;
			Vector3d ray_direction = (ray_origin - camera).normalized();

			// Intersect with the sphere
			const double sphere_radius = 0.9;
			double ed = ray_origin.dot(ray_direction);
			double ee = ray_origin.dot(ray_origin);
			double delta = ed*ed - ee + sphere_radius*sphere_radius;

			if (delta >= 0) {
				// The ray hit the sphere, compute the exact intersection point
				double t = -ed - std::sqrt(delta);
				Vector3d ray_intersection = ray_origin + t * ray_direction;

				// Compute normal at the intersection point
				Vector3d ray_normal = ray_intersection.normalized();

				// Add shading parameter here
				double kd = 0.5, ks = 0.5, phong = 15;
				Vector3d l = (light_position-ray_intersection).normalized();
				Vector3d v = (camera-ray_intersection).normalized();
				Vector3d h = (v+l).normalized();
				diffuse(i,j) = l.transpose() * ray_normal;
				specular(i,j) = h.transpose() * ray_normal;
				diffuse(i,j) = std::max(diffuse(i,j), 0.) * kd;
				specular(i,j) = std::pow(std::max(specular(i,j), 0.), phong) * ks;

				C(i,j) = ambient + diffuse(i,j) + specular(i,j);
				Z(i,j) = C(i,j) - 0.5;

				// Disable the alpha mask for this pixel
				A(i,j) = 1;
			}
		}
	}

	// Save to png
	write_matrix_to_png(C,Z,Z,A,filename+"_red.png");
	write_matrix_to_png(Z,C,Z,A,filename+"_green.png");
	write_matrix_to_png(Z,Z,C,A,filename+"_blue.png");
}

int main() {
	raytrace_sphere();
	raytrace_parallelogram();
	raytrace_perspective();
	raytrace_shading();

	return 0;
}
