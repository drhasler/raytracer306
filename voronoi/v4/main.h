#pragma once
#include <bits/stdc++.h>
#include <lbfgs.c>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

// 2D vector

struct vec2 {
	double x,y;

	vec2 operator+(const vec2& o) const { return { x+o.x, y+o.y }; }
	vec2 operator-(const vec2& o) const { return { x-o.x, y-o.y }; }
	vec2 operator*(double d) const { return { d*x, d*y }; }
	vec2 operator/(double d) const { return *this * (1/d); }

	vec2& operator+=(const vec2& o) { return *this = *this + o; }
	vec2& operator-=(const vec2& o) { return *this = *this - o; }
	vec2& operator*=(double d) { return *this = *this * d; }
	vec2& operator/=(double d) { return *this = *this / d; }

	vec2 orth() const { return {-y, x }; }
};

double dot(const vec2& a, const vec2& b) { return a.x*b.x + a.y*b.y; }
double cross(const vec2& a, const vec2& b) { return a.x*b.y - a.y*b.x; }
double norm(const vec2& v) { return sqrt(dot(v,v)); }

std::ostream& operator<<(std::ostream& os, const vec2& p)
{ return os << '(' << p.x << ',' << p.y << ')'; }

struct Bbox { int xlo, xhi, ylo, yhi; };

// convex power cell

struct Cell {
    vec2 kern; double wk;
    std::vector<vec2> pts = {{0,0},{1,0},{1,1},{0,1}};

    Cell() {}
    Cell(const vec2& k, double w) : kern(k), wk(w) // offset pts by kern
    { std::for_each(pts.begin(), pts.end(), [&](vec2& v){ v -= k; }); }
    void add_pt(vec2 p, double w); // consider other cell

    double area() const;
    double inertia() const; // wrt kern
    vec2 centroid() const;
    Bbox bbox(int res) const;
};

typedef std::vector<Cell> Diagram;
Diagram get_diagram(const std::vector<vec2>& seeds, const double* w);

// random stuff

std::mt19937 rng(time(0));

vec2 random_pos() { // vector in unit square
    std::uniform_real_distribution<double> uni(0,1);
    return { uni(rng), uni(rng) };
}

vec2 random_disk() { // vector in unit disk
    std::uniform_real_distribution<double> uni(0,1);
    double r = uni(rng) + uni(rng);
    r = r > 1 ? 2-r : r;
    double t = 2*M_PI*uni(rng);
    return { r*cos(t), r*sin(t) };
}

std::pair< std::vector<vec2>, std::vector<double> >
petri_dish(int N) { // food in the middle, ppl everywhere
    vec2 C = {.5,.5};
    std::vector<vec2> pts(N);
    std::generate(pts.begin(), pts.end(), random_pos);

    std::vector<double> lam(N);
    std::transform(pts.begin(), pts.end(), lam.begin(),
            [&](vec2& v){ return exp(-6*norm(v-C)); });

    double s = std::accumulate(lam.begin(), lam.end(), 0.);
    std::for_each(lam.begin(), lam.end(), [&](double& d){ d/=s; });

    return { pts, lam };
}

std::vector<vec2> social_distancing(int N) {
    std::vector<vec2> seeds(N);
    std::generate(seeds.begin(), seeds.end(), random_pos);

    double w[N];
    std::fill(w, w+N, 0);

    for (int i = 0; i < 10; i++) { // lloyd's algorithm
        auto cells = get_diagram(seeds, w);
        std::transform(cells.begin(), cells.end(), seeds.begin(),
                [](const Cell& c){ return c.centroid(); });
    }
    return seeds;
}

std::vector<vec2> crop_circle(int N) {
    vec2 C = {.5,.5};
    std::vector<vec2> pts(N);
    std::generate(pts.begin(), pts.end(),
            [&](){ return random_disk() * .3 + C; });
    return pts;
}

// dell optiplex

struct args_t {
    const std::vector<vec2>& seeds;
    const std::vector<double>& lam;
};

static double evaluate(void *instance, const double *w, double *g, const int n, const double step) {
    auto& [seeds, lam] = *(args_t*) instance;
    auto cells = get_diagram(seeds, w);
    double f = 0;
#pragma omp parallel for
    for (int i = 0; i < n; i++) {
        double A = cells[i].area(), 
               I = cells[i].inertia();
        g[i] = A - lam[i];
#pragma omp atomic
        f += w[i] * g[i] - I;
    }
    return f;
}

static int progress(void *instance, const double *x, const double *g, const double fx, const double xnorm, const double gnorm, const double step, int n, int k, int ls) {
    printf("Iteration %d:\n", k);
    printf("  fx = %f\n", fx);
    printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
    printf("\n");
    return 0;
}

void optimize(args_t args, double* x) {
    const int n = args.seeds.size();

    double fx;
    lbfgs_parameter_t param;
    lbfgs_parameter_init(&param);
    int ret = lbfgs(n, x, &fx, evaluate, progress, &args, &param);
    std::cout << "status code: " << ret << std::endl;
}

// enter the simulation

std::pair< std::vector<vec2>, std::vector<double> >
init_sim(int N, int M, double f) {
    auto pts = crop_circle(N); // watr
    auto air = social_distancing(M);
    pts.insert(pts.end(), air.begin(), air.end());

    std::vector<double> lam(N+M);
    std::fill(lam.begin(), lam.begin()+N, f/N);
    std::fill(lam.begin()+N, lam.end(), (1-f)/M);

    return { pts, lam };
}

void save_frame(const Diagram &cells, std::string filename, int frame, int N);
void simulate(int N, int M, double f) {
    auto [pts, lam] = init_sim(N,M,f);
    std::vector<vec2> vel(N+M);
    double k = 30;
    vec2 G = { 0, -9.8 };

    int n_steps = 10;
    double dt = .03;

    double w[N+M];
    std::fill(w, w+N+M, 0);

    for (int i = 0; i < n_steps; i++) {
        optimize( {pts, lam} , w);
        auto cells = get_diagram(pts, w);
        save_frame(cells, "out", i, N);
        for (int i = 0; i < N+M; i++) {
            vec2 f = (cells[i].centroid() - pts[i]) * k;
            if (i < N) f += G; // watr
            vel[i] += f * dt;
            pts[i] += vel[i] * dt;
        }
    }
}

// messy details (warning: includes math)

void Cell::add_pt(vec2 p, double w) {
    p -= kern;
    vec2 p2 = p * ( .5  + (wk-w) / dot(p,p) );
    double th = dot(p2,p);
    auto inside = [&](const vec2& o) { return dot(o,p) < th; };

    int l,r, n = pts.size();
    for (l = 0; l < n && inside(pts[l]); l++);
    if (l == n) return; // all inside
    if (!l) {
        for (l++; l<n && !inside(pts[l]); l++);
        if (l == n) { pts = {}; return; } // all outside
    } else for (l++; !inside(pts[l%n]); l++); // first inside
    for (r = l+1; inside(pts[r%n]); r++); // first outside

    auto find = [&](const vec2& p1, const vec2& p2) -> vec2 {
        double c1 = dot(p, p1), c2 = dot(p, p2);
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
double Cell::area() const {
    const int n = pts.size();
    double ans = 0;
    for (int i = 0; i < n; i++)
        ans += cross(pts[i], pts[(i+1) % n]);
    return ans / 2;
}
double Cell::inertia() const {
    const int n = pts.size();
    double ans = 0;
    for (int i = 0; i < n; i++) {
        const vec2 &a = pts[i], &b = pts[(i+1) % n];
        ans += cross(a,b) * (dot(a,a) + dot(a,b) + dot(b,b));
    }
    return ans / 12;
}
vec2 Cell::centroid() const {
    const int n = pts.size();
    if (!n) return kern;
    vec2 C = {0,0};
    for (int i = 0; i < n; i++) {
        const vec2 &a = pts[i], &b = pts[(i+1) % n];
        C += (a+b) * cross(a,b);
    }
    return kern + C / (6 * this->area());
}
Diagram get_diagram(const std::vector<vec2>& seeds, const double* w) {
    const int n = seeds.size();
    Diagram ans(n);
#pragma omp parallel for
    for (int i = 0; i < n; i++) {
        ans[i] = Cell(seeds[i], w[i]);
        for (int j = 0; j < n; j++)
            if (i != j) ans[i].add_pt(seeds[j], w[j]);
    }
    return ans;
}

// taking pics

Bbox Cell::bbox(int res) const {
    double x0 = 1e9, x1 = -1e9,
           y0 = 1e9, y1 = -1e9;
    for (vec2 v: pts) {
        x0 = std::min(x0, v.x);
        x1 = std::max(x1, v.x);
        y0 = std::min(y0, v.y);
        y1 = std::max(y1, v.y);
    }
    return { (x0 + kern.x)*res, (x1 + kern.x)*res,
             (y0 + kern.y)*res, (y1 + kern.y)*res };
}

void save_frame(const Diagram &cells, std::string filename, int frame, int N) {
    int W = 1000, H = 1000;
    std::vector<unsigned char> image(W*H * 3, 255);
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < cells.size(); i++) {
        const auto& cell = cells[i];
        int n = cell.pts.size();
        if (!n) continue;

        auto [bx0,bx1,by0,by1] = cell.bbox(W);
        for (int y = by0; y < by1; y++)
        for (int x = bx0; x < bx1; x++) {
            vec2 v = vec2{x,y} / W - cell.kern;
            bool in = 1;
            double dst = 1e9;
            for (int j = 0; j < n; j++) {
                vec2 a = cell.pts[j],
                b = cell.pts[(j+1) % n];
                b = b-a;
                a = v-a;

                double alt = cross(b,a) / norm(b);
                if (alt < -1e-3) { in = 0; break; }
                dst = std::min(dst, alt);
            }
            if (!in) continue;

            int idx = ((H - y - 1)*W + x) * 3;
            if (dst <= 2e-3) {
                image[idx + 0] = 0;
                image[idx + 1] = 0;
                image[idx + 2] = 0;
            } else if (i < N) { // is watr
                image[idx + 0] = 0;
                image[idx + 1] = 0;
                image[idx + 2] = 255;
            }
        }
    }
    std::ostringstream os;
    os << filename << frame << ".png";
    stbi_write_png(os.str().c_str(), W, H, 3, &image[0], 0);
}

