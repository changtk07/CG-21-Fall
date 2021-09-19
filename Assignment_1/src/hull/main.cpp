////////////////////////////////////////////////////////////////////////////////
#include <algorithm>
#include <complex>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>
////////////////////////////////////////////////////////////////////////////////

typedef std::complex<double> Point;
typedef std::vector<Point> Polygon;

double inline det(const Point &u, const Point &v) {
	return u.real()*v.imag() - v.real()*u.imag();
}

double inline dist(const Point &u, const Point &v) {
	return (u.real()-v.real())*(u.real()-v.real()) +
				 (u.imag()-v.imag())*(u.imag()-v.imag());
}

struct Compare {
	Point p0; // Leftmost point of the poly
	bool operator ()(const Point &p1, const Point &p2) {
		double r = det(p1-p0, p2-p0);
		return r > 0 || (r == 0 && dist(p0, p1) > dist(p0, p2));
	}
};

bool inline salientAngle(Point &a, Point &b, Point &c) {
	return det(b-a, c-a) > 0;
}

////////////////////////////////////////////////////////////////////////////////

Polygon convex_hull(std::vector<Point> &points) {
  Compare order;
  order.p0 = *min_element(points.begin(), points.end(), [](Point &p, Point &q) {
      return p.imag() < q.imag() || (p.imag() == q.imag() && p.real() < q.real());
  });

  std::sort(points.begin(), points.end(), order);
  auto uiq = std::unique(points.begin(), points.end(), [&](Point &p, Point &q) {
      return det(p-order.p0, q-order.p0) == 0;
  });
  points.erase(uiq, points.end());

  Polygon hull{order.p0};
  for (Point &v : points) {
      size_t n;
      while ((n = hull.size()) > 1 && !salientAngle(hull[n-2], hull[n-1], v)) {
          hull.pop_back();
      }
      hull.push_back(v);
  }
  return hull;
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

void save_obj(const std::string &filename, Polygon &poly) {
	std::ofstream out(filename);
	if (!out.is_open()) {
		throw std::runtime_error("failed to open file " + filename);
	}
	out << std::fixed;
	for (const auto &v : poly) {
		out << "v " << v.real() << ' ' << v.imag() << " 0\n";
	}
	for (size_t i = 0; i < poly.size(); ++i) {
		out << "l " << i+1 << ' ' << 1+(i+1)%poly.size() << "\n";
	}
	out << std::endl;
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char * argv[]) {
	if (argc <= 2) {
		std::cerr << "Usage: " << argv[0] << " points.xyz output.obj" << std::endl;
	}
	std::vector<Point> points = load_xyz(argv[1]);
	Polygon hull = convex_hull(points);
	// for (Point &p : hull) std::cout << p << std::endl;
	save_obj(argv[2], hull);
	return 0;
}
