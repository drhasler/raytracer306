#include <lbfgs.c>
#include "main.h"
#include "svg_writer.cpp"

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

Diagram get_diagram(const std::vector<vec2>& seeds, const double* w) {
    Diagram ans;
    const int n = seeds.size();
    for (int i = 0; i < n; i++) {
        Cell cell(seeds[i], w[i]);
        for (int j = 0; j < n; j++)
            if (i != j) cell.add_pt(seeds[j], w[j]);
        ans.push_back(cell);
    }
    return ans;
}

static lbfgsfloatval_t evaluate(void *instance, const lbfgsfloatval_t *w, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step) {
    auto& [seeds, lam] = *(args_t*) instance;
    auto cells = get_diagram(seeds, w);
    double f = 0;
    for (int i = 0; i < n; i++) {
        double A = cells[i].area(), 
               I = cells[i].inertia();
        g[i] = A - lam[i];
        f += w[i] * g[i] - I;
        // std::cout << w[i] << ' ' << A << ' ' << I << ' ' << lam[i] << '\n';
    }
    std::cout << f << '\n';
    /*
    save_svg("out.svg", cells);
    getchar();
    */
    return f;
}

static int progress(void *instance, const lbfgsfloatval_t *x, const lbfgsfloatval_t *g, const lbfgsfloatval_t fx, const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step, int n, int k, int ls) {
    printf("Iteration %d:\n", k);
    printf("  fx = %f\n", fx);
    printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
    printf("\n");
    return 0;
}

void get_weights(const std::vector<vec2>& seeds,
        const std::vector<double>& lam) {
    const int n = seeds.size();
    args_t args { seeds, lam };

    int i, ret = 0;
    lbfgsfloatval_t fx;
    lbfgsfloatval_t *x = lbfgs_malloc(n);
    lbfgs_parameter_t param;
    if (x == NULL) {
        printf("ERROR: Failed to allocate a memory block for variables.\n");
        exit(1);
    }
    /* Initialize the variables. */
    std::fill(x, x+n, 0);
    /* Initialize the parameters for the L-BFGS optimization. */
    lbfgs_parameter_init(&param);
    /*param.linesearch = LBFGS_LINESEARCH_BACKTRACKING;*/
    /* Start the L-BFGS optimization; this will invoke the callback functions
       evaluate() and progress() when necessary.  */
    ret = lbfgs(n, x, &fx, evaluate, progress, &args, &param);
    /* Report the result. */
    printf("L-BFGS optimization terminated with status code = %d\n", ret);

    auto partition = get_diagram(seeds, x);
    save_svg("out.svg", partition);
    lbfgs_free(x);
}

int main() {
    const int N = 100;
    auto [pts,lam] = sample(N);

    get_weights(pts, lam);
}
