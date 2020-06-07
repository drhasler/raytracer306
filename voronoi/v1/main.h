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

std::ostream& operator<<(std::ostream& os, const vec2& p) {
	os << '(' << p.x << ',' << p.y << ')';
	return os;
}

std::vector<vec2> sample_uni(const size_t N) {
    std::mt19937 rng(time(0));
    std::uniform_real_distribution<double> uni(0,1);

    std::vector<vec2> pts(N);
    std::generate(pts.begin(), pts.end(), [&]()
            { return vec2{uni(rng), uni(rng)}; });
    return pts;
}

struct Cell {
    vec2 kern;
    std::vector<vec2> pts = {{0,0},{1,0},{1,1},{0,1}};

    Cell(const vec2& k) : kern(k) { // center points in kern
        std::transform(pts.begin(), pts.end(), pts.begin(),
                [&](vec2& v) { return v-k; });
    }

    void add_pt(vec2 p) {
        p = (p-kern) * .5;
        double th = dot(p,p);
        auto inside = [&](const vec2& o) { return dot(o,p) < th; };

        int l,r, n = pts.size();
        for (l = 0; l < n && inside(pts[l]); l++);
        if (l == n) return; // all inside
        for (l++; !inside(pts[l%n]); l++); // first inside
        for (r = l+1; inside(pts[r%n]); r++); // first outside

        auto find = [&](const vec2& p1, const vec2& p2) -> vec2 {
            double c1 = dot(p, p1),
                   c2 = dot(p, p2);
            // a + b = 1, a c1 + b c2 = th
            double a = (th-c2) / (c1-c2),
                   b = 1 - a;
            return p1*a + p2*b;
        };
        vec2 pl = find(pts[(l-1) % n], pts[l % n]),
             pr = find(pts[(r-1) % n], pts[r % n]);

        std::vector<vec2> nxt = { pr, pl };
        for (int i = l; i < r; i++) nxt.push_back(pts[i % n]);
        std::swap(pts, nxt);
    }
};

typedef std::vector<Cell> Diagram;

