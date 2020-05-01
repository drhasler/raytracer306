struct Sphere {
    Sphere(const vec& C, double R, const Material& M)
        : C(C), R(R), M(M) {}
    double intersect(const Ray& r, vec& P, vec& N) const {
        vec OC = C - r.O;
        double uOC = dot(r.u,OC);
        double delta = uOC*uOC - (dot(OC,OC) - R*R);
        if (delta < 0) return inf;
        delta = sqrt(delta);
        double t;
        if ((t = uOC - delta) > eps) {
            P = r.O + r.u * t;
            N = norm(P-C);
        } else if ((t = uOC + delta) > eps) {
            P = r.O + r.u * t;
            N = norm(C-P);
        } else return inf;
        return t;
    }
    vec C;
    double R;
    Material M;
};
