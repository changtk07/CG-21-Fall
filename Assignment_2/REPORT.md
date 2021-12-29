# Assignment 2 Report
This is a short report of what I've implemented in the second assignment as well as comments and explanations.


## System Info
- OS: macOs Big Sur 11.6
- Compiler: AppleClang 12.0.5.12050022


## Parallelogram Orthographic Projection
We know the representation of a ray `p = e + td`. And the reprentation of a parallelogram, an origin `o` and two sides vector `u` and `v`.

First calculate the plane of parallelogram. Because this is a special case that rays are all in the same direction and are all `[0, 0, -1]`, so to make it easier, I project everything to the xy-plane. It can be done by converting a `Vector3d` to a `Vector2d` keeping only `x` and `y` coordinates.

Now we have the ray-plane intersetion, the intersection is in the parallelogram when `0 <= m <= 1` and `0 <= n <= 1`. `m` and `n` can be easily calculate by

```
m = det(p,v) / det(u,v)
n = det(p,u) / det(v,u) = -det(p,u) / det(u,v)
```

Finally, the normal at intersection is `o + u × v`


## Parallelogram Perspective Projection
This is almost the same as orthpgraphic, except for the ray. I choose a point as the camera `c`, and the ray direction is now `e - c` rather than `[0, 0, -1]`.

And then the intersection point. The idea is basicaly the same, except this time I calculate `t' = ((o - e) * n) / (d * n)`, where `n` is the plane normal. And the itersection is `p = e + t'd`.

Everything else is exactly the same.


## Sphere Perspective And Shading
For perspective ray-sphere intersection,
```
ray: p = e + td
sphere: (p-c) * (p-c) - r^2 = 0
c = [0, 0, 0]

dd t^2 + 2ed t + ee - r^2 = 0

delta (b^2 - 4ac) = (ed)^2 - (ee - r^2)
t = -ed - sqrt(delta)
```
so the intersection is `e +td`.

For shading, I define three parameters for diffuse and specular, `kd`, `ks` and `phong`, and by using the formula presented in class.

For coloring, I simply reduce the light in other two channels by 0.5.


## Screenshots
Parallelogram orthographic, tesing differnet pgram angle 
1 | 2
:-------------------------:|:-------------------------:
![](img/plane_orthographic_1.png) |  ![](img/plane_orthographic_2.png)

Parallelogram perspective, tesing differnet pgram angle
1 | 2
:-------------------------:|:-------------------------:
![](img/plane_perspective_1.png) |  ![](img/plane_perspective_2.png)

Sphere shading, testing different parameter values
Groups | R | G | B
:-------------------------:|:-------------------------:|:-------------------------:|:-------------------------:
control group | ![](img/shading_red_1.png) | ![](img/shading_green_1.png) | ![](img/shading_blue_1.png)
higher ambient | ![](img/shading_red_2.png) | ![](img/shading_green_2.png) | ![](img/shading_blue_2.png)
higher diffuse | ![](img/shading_red_3.png) | ![](img/shading_green_3.png) | ![](img/shading_blue_3.png)
higher specular | ![](img/shading_red_4.png) | ![](img/shading_green_4.png) | ![](img/shading_blue_4.png)
higher phong | ![](img/shading_red_5.png) | ![](img/shading_green_5.png) | ![](img/shading_blue_5.png)
