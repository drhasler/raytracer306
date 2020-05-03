# raytracer 306
a simple C++ raytracer  

author: Matthias Hasler  
CSE306 - Ecole polytechnique, Bachelor of Science - SS20  

## files
- `vec.cpp` vector and ray classes
- `random.cpp` random sampling of rays
- `mat.cpp` material and UV classes
- `sphere.cpp` self explanatory
- `mesh.cpp` triangle, mesh, BVH tree and tiny `.obj` parser
- `main.cpp` scene, camera classes; rendering

(for grading)
`pandoc <(echo "# Report"; cat README.md) -o report.pdf`

## features
this raytracer can render:

- diffuse, refractive, reflective and emitting spheres
- meshes, with UV mapping

it does not handle phong reflections nor participating media (smoke, etc.)

## details and tricks
### parallelism
OpenMP, dynamic scheduling over the lines rendered

### monte carlo
- effects like soft shadows are achieved by shooting
several rays per pixel and sampling bouncing rays from random distributions
- depth of field is achieved by shooting rays from various positions
in the aperture aiming at the focal plane
- for antialiasing, we introduce normal noise when aiming at the focal plane

### refraction
We check whether the ray coming in or out.
Since we only care about n1/n2
[(Snell's law)](https://en.wikipedia.org/wiki/Snell%27s_law),
we can simply invert the refractive index depending on the case.
To emulate refractive indices of nested objects,
one should divide the real refractive index
by the ambient one.

### parser
expects `f v/vt/vn` and `vn x y z`

## external res
- `Models_F0202A090/*`
[cat model and texture](http://www.cadnav.com/3d-models/model-47556.html)
- `stb_image_write.h` and `stb_image.h`
[STB library](https://github.com/nothings/stb)

## renders
recursion depth limit set to 5, resolution 600 x 400

![little cat and hollow sphere (5 samples)](renders/cat2.png)

![big boi (20 samples. 2m53s)](renders/bigg.png)

![from left to right: spheres showing diffraction, refraction, reflection and
a hollow refractive sphere (90 samples, 8.6s)](renders/balls.png)
