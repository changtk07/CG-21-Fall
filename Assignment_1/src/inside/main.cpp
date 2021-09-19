////////////////////////////////////////////////////////////////////////////////
#include <algorithm>
#include <complex>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>
#include <cmath>
////////////////////////////////////////////////////////////////////////////////

typedef std::complex<double> Point;
typedef std::vector<Point> Polygon;

double inline det(const Point &u, const Point &v) {
	return u.real()*v.imag() - v.real()*u.imag();
}

// Return true iff [a,b] intersects [c,d], and store the intersection in ans
bool intersect_segment(const Point &a, const Point &b, const Point &c, const Point &d, Point &ans) {
	auto t = det(a-c, c-d) / det(a-b, c-d);
	auto u = det(b-a, a-c) / det(a-b, c-d);

	auto x1 = a.real(), y1 = a.imag();
	auto x2 = b.real(), y2 = b.imag();
	ans = { x1 + t * (x2-x1), y1 + t * (y2-y1) };

	return 0 <= t && t <= 1 && 0 <= u && u <= 1;
}

////////////////////////////////////////////////////////////////////////////////

bool is_inside(const Polygon &poly, const Point &query) {
	// 1. Compute bounding box and set coordinate of a point outside the polygon
	double min_x = poly[0].real(), min_y = poly[0].imag();
	for (auto &p : poly) {
		min_x = std::min(min_x, p.real());
		min_y = std::min(min_y, p.imag());
	}
	Point outside(min_x-10, min_y-10);
	// 2. Cast a ray from the query point to the 'outside' point, count number of intersections
	int count = 0;
	for (size_t i = 0; i < poly.size(); ++i) {
		auto &p1 = poly[i], &p2 = poly[ (i+1)%poly.size() ];
		Point inter;
		if (intersect_segment(p1, p2, query, outside, inter)) {
			if (inter == p1)
				count += det(query-outside, p2-outside) > 0;
			else if (inter == p2)
				count += det(query-outside, p1-outside) > 0;
			else
				count += 1;
		}
	}
	return count % 2;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Point> load_xyz(const std::string &filename) {
	std::vector<Point> points;
	std::ifstream in(filename);
	int N;
	in >> N;
	points.resize(N);
	for (auto &v : points) {
		double x, y, z;
		in >> x >> y >> z;
		v = {x, y};
	}
	return points;
}

Polygon load_obj(const std::string &filename) {
	Polygon poly;
	std::ifstream in(filename);
	char c;
	double x, y, z;
	while (in >> c && c == 'v') {
		in >> x >> y >> z;
		poly.push_back({x, y});
	}
	return poly;
}

void save_xyz(const std::string &filename, const std::vector<Point> &points) {
	std::ofstream out(filename);
	out << points.size() << '\n';
	for (auto &p : points) {
		out << p.real() << ' ' << p.imag() << " 0\n";
	}
	out << std::endl;
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char * argv[]) {
	if (argc <= 3) {
		std::cerr << "Usage: " << argv[0] << " points.xyz poly.obj result.xyz" << std::endl;
	}
	std::vector<Point> points = load_xyz(argv[1]);
	Polygon poly = load_obj(argv[2]);
	std::vector<Point> result, diff;
	for (size_t i = 0; i < points.size(); ++i) {
		if (is_inside(poly, points[i])) {
			result.push_back(points[i]);
		}
		else {
			diff.push_back(points[i]);
		}
	}
	// for (Point &p : result) std::cout << p << std::endl;
	save_xyz(argv[3], result);
	save_xyz("diff.xyz", diff);

	return 0;
}
