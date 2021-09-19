# Assignment 1 Report
This is a short report of what I've implemented in the first assignment as well as comments and explanations.


## System Info
- OS: macOs Big Sur 11.6
- Compiler: AppleClang 12.0.5.12050022


## Convex Hull
In this section, I implemented Graham's scan algorithm to construct the convex hull of a set of 2D points. I completed the `TODO` part in the file `src/hull/main.cpp`, specifically:

- `double inline det(const Point &u, const Point &v)`
- `struct Compare`
- `bool inline salientAngle(Point &a, Point &b, Point &c)`
- `Polygon convex_hull(std::vector<Point> &points)`
- `std::vector<Point> load_xyz(const std::string &filename)`

and added a helper function:

- `double inline dist(const Point &u, const Point &v)`

#### `double inline det(const Point &u, const Point &v)`
This function returns the determinant of two 2d-vectors, `u` and `v`.

By the definition `det(u, v) = ux * vy - vx * uy`


#### `double inline dist(const Point &u, const Point &v)`
This function returns the distance of two points `u` and `v`.

#### `struct Compare`
This structure provided the comparison that used to sort points in counter-clockwise order w.r.t `p0`. In case of tie, we sort the points by distance to `p0` in descending order. To compare the angle of two vectors `vec1 (p0, p1)` and `vec2 (p0, p2)`, we compute their determinant, if it's greater than 0, then `vec1` to `vec2` is a counter-clockwise turn; if it's less than 0, it's a clockwise turn; if it's 0, then `vec1` and `vec2` are collinear.


#### `bool inline salientAngle(Point &a, Point &b, Point &c)`
This function returns true if the vector `(a, b)` to `(a, c)` makes a counter-clockwise turn. Again, we make use of their determinant.


#### `Polygon convex_hull(std::vector<Point> &points)`
This function is the main part of Graham's algorithm.

First we find the bottom-left-most point and assign it to `p0`. This is simple, just compare their `x` and `y` coordinate.

Then sort all the points using `struct Compare`.

Next remove all collinear points w.r.t `p0`, only keep the furthest one. Here I use `std::unique` and `std::erase` to achieve this.

And finally traverse all the points and construct the convex hull. If the last two points in convex hull and the new point make a counter-clockwise turn, then add it to the convex hull. Otherwise, first remove the last point in convex hull before adding the new point.


#### `std::vector<Point> load_xyz(const std::string &filename)`
This function reads points from a file to a `std::vector` and returns it.


## Point In Polygon
In this section, I implemented an algorithm that filter all the points that are inside a given polygon. I completed the `TODO` part in the file `/src/inside/main.cpp`, specifically:

- `double inline det(const Point &u, const Point &v)`
- `bool intersect_segment(const Point &a, const Point &b, const Point &c, const Point &d, Point &ans)`
- `bool is_inside(const Polygon &poly, const Point &query)`
- `std::vector<Point> load_xyz(const std::string &filename)`
- `Polygon load_obj(const std::string &filename)`
- `void save_xyz(const std::string &filename, const std::vector<Point> &points)`

#### `double inline det(const Point &u, const Point &v)`
This function is identical to the one in previous section.


#### `bool intersect_segment(const Point &a, const Point &b, const Point &c, const Point &d, Point &ans)`
This function returns true iff two line segment `[a, b]` and `[c, d]` intersects, and store the intersection in `ans`.

The approach I used is described here: https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection#Given_two_points_on_each_line_segment


#### `bool is_inside(const Polygon &poly, const Point &query)`
This function returns true if the point `query` is inside the polygon `poly`.

First we find a point `outside` that is outside the polygon. The point I used is `(min_x - 10, min_y - 10)`, where `min_x` is the minimum of `x` coordinate, and `min_y` is the minimum of `y` coordinate.

Then traverse each edge in `poly` and count the number of intersections with `[query, outside]`. If the number is odd then `query` is inside. I also handled special cases when an intersection is one of then vertex of polygon and when two line segments are overlapping.

If the intersection point is a vertex of a polygon edge, then the intersection counts only if the second vertex of the edge lies below the ray. This is achieved by calculating their determinants.

If one of the edge and the ray are overlapping, then don't count the intersection.


#### `std::vector<Point> load_xyz(const std::string &filename)`
This function is identical to the one in previous section.


#### `Polygon load_obj(const std::string &filename)`
This function reads a polygon from `filename` to a `Polygon` and returns it.


#### `void save_xyz(const std::string &filename, const std::vector<Point> &points)`
This function writes a list of points to `filename`.


## Results
Here's the result of convex hull:

![](img/convex_hull.png?raw=true)

And the result of points in side:

![](img/inside.png?raw=true)


## Authors
Kevin Chang: tc3149@nyu.edu
