struct vec {
    double x=0,y=0,z=0;
    constexpr vec() {}
    vec(double x, double y, double z)
        : x(x), y(y), z(z) {}
    vec& operator+=(const vec& o) {
        x += o.x; y += o.y; z += o.z;
        return *this;
    }
    vec operator+(const vec& o) const {
        return {x+o.x, y+o.y, z+o.z};
    }
    vec operator-(const vec& o) const {
        return {x-o.x, y-o.y, z-o.z};
    }
    vec operator*(double f) const {
        return {f*x, f*y, f*z};
    }
    vec operator/(double f) const {
        return *this * (1/f);
    }
    vec operator*(const vec& o) const {
        return {x*o.x, y*o.y, z*o.z};
    }
};

vec operator*(double f, const vec& v)
{ return {f*v.x, f*v.y, f*v.z}; }

double dot(const vec& a, const vec& b)
{ return a.x*b.x + a.y*b.y + a.z*b.z; }

vec cross(const vec& a, const vec& b) {
    return { a.y*b.z - a.z*b.y,
             a.z*b.x - a.x*b.z,
             a.x*b.y - a.y*b.x };
}

double abs(const vec& v)
{ return sqrt(dot(v,v)); }
vec norm(const vec& v)
{ return v / abs(v); }
vec norm(const vec& v, double& d)
{ return v / (d = abs(v)); }

vec refl(const vec& wi, const vec& N)
{ return wi - N * (2*dot(wi,N)); }

vec refr(const vec& wi, const vec& N, double n) { // n = n1/n2
    double f = dot(wi,N);
    vec wT = n * (wi - f * N);
    f = 1 - n*n * (1-f*f);
    if (f < 0) return refl(wi,N);
    return wT - sqrt(f) * N;
}

struct Ray { vec O,u; };

std::ostream& operator<<(std::ostream& os, const vec& v)
{ return os << '(' << v.x << ',' << v.y << ',' << v.z << ')'; }
std::ostream& operator<<(std::ostream& os, const Ray& r)
{ return os << "Ray{" << r.O << r.u << '}'; }
