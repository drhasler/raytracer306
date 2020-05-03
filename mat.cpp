enum MaterialType { DIFFUSE, LIGHT, MIRROR, REFR };

struct Material {
    Material(const vec& albedo) : albedo(albedo), T(DIFFUSE) {}
    Material(double f, MaterialType T) : T(T), I(f) {} // REFR or LIGHT
    Material() : T(MIRROR) {}
    union { vec albedo; double I; double n; };
    MaterialType T;
};

struct UV {
    std::vector<vec> img;
    int w,h;

    vec col(double x, double y) const {
        x -= floor(x);
        y -= floor(y);
        int X = x*w;
        int Y = (1-y)*h;
        return img[Y*w+X];
    }
};
