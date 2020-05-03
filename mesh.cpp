struct tri_idx {
    int vtx[3];
    int norm[3];
    int uv[3];
};

struct bbox {
    double lx=inf,hx=-inf, ly=inf,hy=-inf, lz=inf,hz=-inf;
    bbox() {}
    bbox(int l, int r, const std::vector<tri_idx>& tri,
         const std::vector<vec>& vtx) {
        for ( ; l<r; l++) for (int i: tri[l].vtx) {
            auto [x,y,z] = vtx[i];
            if (x < lx) lx = x;
            if (x > hx) hx = x;
            if (y < ly) ly = y;
            if (y > hy) hy = y;
            if (z < lz) lz = z;
            if (z > hz) hz = z;
        }
    }
    bool intersect(const Ray& r) const {
        double tl = -inf, th = inf;
        for (auto [d,p,l,h]: { std::tuple
                {r.u.x,r.O.x,lx,hx},
                {r.u.y,r.O.y,ly,hy},
                {r.u.z,r.O.z,lz,hz}, })
        {
            if (abs(d) < eps) continue;
            if (d > 0) {
                tl = std::max(tl,(l-p)/d);
                th = std::min(th,(h-p)/d);
            } else {
                tl = std::max(tl,(h-p)/d);
                th = std::max(th,(l-p)/d);
            }
        }
        return th > 0 & tl < th;
    }
};

struct BVH {
    BVH(int l, int r,
        std::vector<tri_idx>& trii,
        const std::vector<vec>& vtx) : bb(l,r,trii,vtx) {
        if (r-l <= 5) {
            leaf = 1;
            tri = {l,r};
            return;
        }
        leaf = 0;
        double dx = bb.hx - bb.lx;
        double dy = bb.hy - bb.ly;
        double dz = bb.hz - bb.lz;

        auto center = [&](const tri_idx& i)
        { return vtx[i.vtx[0]] + vtx[i.vtx[1]] + vtx[i.vtx[2]]; };

        if (dx > dy && dx > dz) { // sort along x
            std::sort(&trii[l],&trii[r],[&](const tri_idx& a, const tri_idx& b)
                    { return center(a).x < center(b).x; });
        } else if (dy > dz) { // sort along y
            std::sort(&trii[l],&trii[r],[&](const tri_idx& a, const tri_idx& b)
                    { return center(a).y < center(b).y; });
        } else { // sort along z
            std::sort(&trii[l],&trii[r],[&](const tri_idx& a, const tri_idx& b)
                    { return center(a).z < center(b).z; });
        }
        int m = (l+r) / 2;
        child = { new BVH(l,m,trii,vtx), new BVH(m,r,trii,vtx) };
    }
    union { std::pair<int,int> tri;
            std::pair<BVH*,BVH*> child; };
    bbox bb;
    bool leaf;
};

struct Mesh {
    void shift(vec o) { for (auto& v: vtx) v += o; }
    void scale(double f) { for (auto& v: vtx) v = v*f; }
    void rot(double yaw, double pit, double rol) {
        for (auto* vs: { &vtx, &norm })
        for (auto& v: *vs) {
            v = { v.x*cos(yaw) + v.z*sin(yaw),
                  v.y,
                 -v.x*sin(yaw) + v.z*cos(yaw) };
            v = { v.x,
                  v.y*cos(pit) + v.z*sin(pit),
                 -v.y*sin(pit) + v.z*cos(pit) };
            v = { v.x*cos(rol) - v.y*sin(rol),
                  v.x*sin(rol) + v.y*cos(rol),
                  v.z };
        }
    }
    void make_bvh() { bvh = new BVH(0,tri.size(), tri, vtx); }
    double intersect(const Ray& r, vec& P, vec& N, vec& alb) const {
        double t = inf;
        std::vector<BVH*> dfs = { bvh };
        while (dfs.size()) {
            auto cur = dfs.back(); dfs.pop_back();
            if (!cur->bb.intersect(r)) continue;
            if (!cur->leaf) {
                dfs.push_back(cur->child.first);
                dfs.push_back(cur->child.second);
                continue;
            }
            auto [bi,ei] = cur->tri;
            while (bi < ei) {
                auto [ia,ib,ic] = tri[bi].vtx;
                vec A = vtx[ia], B = vtx[ib], C = vtx[ic];
                vec e1 = B - A, e2 = C - A;
                vec OA = r.O - A;
                double den = dot(e1, cross(r.u,e2));
                double b = dot(OA, cross(r.u,e2)) / den;
                double c = dot(e1,cross(r.u,OA)) / den;
                double a = 1-b-c;
                if ( a >= 0 & a <= 1
                   & b >= 0 & b <= 1
                   & c >= 0 & c <= 1 ) {
                    double alt = dot(e1,cross(e2,OA)) / den;
                    if (alt < t) {
                        t = alt;
                        P = r.O + t*r.u;
                        auto [na,nb,nc] = tri[bi].norm;
                        N = a*norm[na] + b*norm[nb] + c*norm[nc];
                        auto [ta,tb,tc] = tri[bi].uv;
                        alb = a*uv[ta]+b*uv[tb]+c*uv[tc];
                        // std::cout << alb << '\n';
                        alb = mat.col(alb.x, alb.y);
                    }
                }
                bi++;
            }
        }
        return t;
    }
    UV mat;
    std::vector<vec> vtx,norm,uv;
    std::vector<tri_idx> tri;
    BVH* bvh;
};

Mesh load_cat() {
    std::ifstream is("Models_F0202A090/cat.obj");
    if (!is.is_open()) { std::cerr << "couldnt open it\n"; exit(1); }

    Mesh m; std::string s;
    while (is >> s) {
        if (s[0] == 'v') {
            double x,y,z;
            is >> x >> y >> z;
            ( s.size()==1 ? m.vtx : s[1]=='n' ? m.norm : m.uv )
                .push_back(vec{x,y,z});
        } else if (s[0] == 'f') {
            char c; tri_idx t;
            for (int i=0; i<3; i++) {
                is >> t.vtx[i] >> c >> t.uv[i] >> c >> t.norm[i];
                t.vtx[i]--; t.uv[i]--; t.norm[i]--;
            }
            m.tri.push_back(t);
        }
        is.ignore(256,'\n');
    }
    int nc;
    m.mat.img = stbi_load("Models_F0202A090/cat_diff.png",
            &m.mat.w, &m.mat.h, &nc, 3);
    m.shift(vec{0,-20,0});
    m.scale(1.2);
    m.rot(-2.5,0,0);
    m.make_bvh();
    return m;
}

std::ostream& operator<<(std::ostream& os, const bbox& b)
{   return os << "bbox{"
    << b.lx << ',' << b.hx << ','
    << b.ly << ',' << b.hy << ','
    << b.lz << ',' << b.hz << '}'; }
