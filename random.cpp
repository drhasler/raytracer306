thread_local static std::random_device rd;
thread_local static std::mt19937 rng(rd());
thread_local static std::uniform_real_distribution udis(0.,1.);

double uni() { return udis(rng); }
void box_muller(double stdev, double& x, double& y) {
    double r = sqrt(-2*log(uni()))*stdev;
    double ang = 2*M_PI * uni();
    x = r*cos(ang), y = r*sin(ang);
}
vec random_cos(const vec& N) {
    double z = uni(), ang = 2*M_PI*uni();
    double r = sqrt(1-z);
    double x = r * cos(ang);
    double y = r * sin(ang);
    z = sqrt(z);
    vec u = abs(N.x) > sqrt(.5) ? vec{0,1,0} : vec{1,0,0};
    u = norm(u - dot(N,u)*N);
    vec v = cross(N,u);
    return x*u + y*v + z*N;
}
vec random_pow(const vec& N, double alpha) {
    double z = pow(uni(),1/(alpha+1)), ang = 2*M_PI * uni();
//    double r = sqrt()
}

vec uni_sphere() {
    double r1 = uni(), r2 = uni();
    return {
        cos(2*M_PI*r1) * sqrt(r2*(1-r2)),
        sin(2*M_PI*r1) * sqrt(r2*(1-r2)),
        1 - 2*r2
    };
}
