#include <bits/stdc++.h>

struct vec2 {
	double x,y;

	vec2 operator+(const vec2& o) const { return { x+o.x, y+o.y }; }
	vec2 operator-(const vec2& o) const { return { x-o.x, y-o.y }; }
	vec2 operator*(double d) const { return { d*x, d*y }; }

	vec2& operator+=(const vec2& o) { return *this = *this + o; }
	vec2& operator-=(const vec2& o) { return *this = *this - o; }
	vec2& operator*=(double d) { return *this = *this * d; }

	vec2 orth() const { return {-y, x }; }
};

double dot(const vec2& a, const vec2& b) { return a.x*b.x + a.y*b.y; }
double cross(const vec2& a, const vec2& b) { return a.x*b.y - a.y*b.x; }
double norm(const vec2& v) { return sqrt(dot(v,v)); }

std::ostream& operator<<(std::ostream& os, const vec2& p) {
	os << '(' << p.x << ',' << p.y << ')';
	return os;
}

std::pair< std::vector<vec2>, std::vector<double> >
sample(const size_t N) {
    std::mt19937 rng(time(0));
    std::uniform_real_distribution<double> uni(0,1);
    vec2 C = {.5,.5};

    std::vector<vec2> pts(N);
    std::generate(pts.begin(), pts.end(),
            [&](){ return vec2{uni(rng), uni(rng)}; });
    std::vector<double> lam(N);
    std::transform(pts.begin(), pts.end(), lam.begin(),
            [&](vec2& v){ return exp(-6*norm(v-C)); });
    double s = std::accumulate(lam.begin(), lam.end(), 0.);
    std::for_each(lam.begin(), lam.end(), [&](double& d){ d/=s; });
    return { pts, lam };
}

struct Cell {
    vec2 kern; double wk;
    std::vector<vec2> pts = {{0,0},{1,0},{1,1},{0,1}};

    Cell(const vec2& k, double w) : kern(k), wk(w) { // center points in kern
        std::for_each(pts.begin(), pts.end(), [&](vec2& v) { v -= k; });
    }
    void add_pt(vec2 p, double w);
    double area() const {
        const int n = pts.size();
        double ans = 0;
        for (int i = 0; i < n; i++)
            ans += cross(pts[i], pts[(i+1) % n]);
        return ans / 2;
    }
    double inertia() const {
        const int n = pts.size();
        double ans = 0;
        for (int i = 0; i < n; i++) {
            const vec2 &a = pts[i], &b = pts[(i+1) % n];
            ans += cross(a,b) * (dot(a,a) + dot(a,b) + dot(b,b));
        }
        return ans / 12;
    }
};

typedef std::vector<Cell> Diagram;

struct args_t {
    const std::vector<vec2>& seeds;
    const std::vector<double>& lam;
};

