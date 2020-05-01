enum MaterialType { DIFFUSE, LIGHT, MIRROR, REFR };

struct Material {
    Material(const vec& albedo) : albedo(albedo), T(DIFFUSE) {}
    Material(double f, MaterialType T) : T(T), I(f) {} // REFR or LIGHT
    Material() : T(MIRROR) {}
    union { vec albedo; double I; double n; };
    MaterialType T;
};

struct UV {
    unsigned char* img;
    int w,h;

    vec col(double x, double y) const {
        if ((x -= round(x)) < 0) x += 1;
        if ((y -= round(y)) < 0) y += 1;
        int X = x*w;
        int Y = (1-y)*h;
        return vec{
            (double)img[3*(Y*w+X) + 0],
            (double)img[3*(Y*w+X) + 1],
            (double)img[3*(Y*w+X) + 2],
        } / 255;
    }
};
