#include <bits/stdc++.h>
#include <omp.h>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

const double inf = 1e18;
const double eps = 1e-5;

#include "vec.cpp"
#include "mat.cpp"
#include "sphere.cpp"
#include "random.cpp"
#include "mesh.cpp"
#include "object.cpp"

struct Scene {
    Scene() {}
    std::vector<Sphere> spheres;
    std::vector<Mesh> meshes;
    double intersect(const Ray& r, vec& P, vec& N, const Sphere*& sp) const {
        double t = inf;
        vec P2,N2;
        for (const auto& s: spheres) {
            double alt = s.intersect(r,P2,N2);
            if (alt < t) {
                t = alt; sp = &s;
                P = P2; N = N2;
            }
        }
        return t;
    }
    vec get_color(const Ray& r, int depth, bool last_diffuse=0) const {
        if (!depth) return vec();
        const Sphere* sp = 0; const Material* mp; vec P,N;
        double ts = intersect(r,P,N,sp);
        if (sp) mp = &sp->M;
        double tm = inf;
        Material mm{vec()};
        for (const auto& m: meshes) {
            tri_idx ii;
            vec P2,N2, alb;
            double alt = m.intersect(r,P2,N2,alb);
            if (alt < tm) {
                tm = alt; mp = &mm;
                P = P2; N = N2; mm.albedo = alb;
            }
        }
        if (ts == inf) return vec();

        if (mp->T == DIFFUSE) {
            vec Lo{};
            for (const auto& l: spheres) if (l.M.T == LIGHT) {
                double d;
                vec D = norm(P-l.C);
                vec Np = random_cos(D);
                vec xp = l.R*Np + l.C;
                vec wi = norm(P-xp, d);

                Ray r { xp, wi }; vec N2,P2; const Sphere* tmp;
                if (intersect(r,N2,P2,tmp) < d - eps) continue;
                double pdf = dot(Np,D) / (M_PI*l.R*l.R);
                Lo += - 1 / (4*M_PI*M_PI*l.R*l.R) / M_PI
                    * dot(N,wi) * dot(Np,wi) / (d*d*pdf)
                    * l.M.I * mp->albedo;
            }
            return Lo + mp->albedo * get_color({P,random_cos(N)}, depth-1, 1);
        } if (mp->T == LIGHT) {
            if (last_diffuse) return vec(); // dont count twice
            vec OC = r.O - sp->C;
            double R2 = dot(OC,OC);
            return mp->I / (4*M_PI*M_PI*R2) * vec(1,1,1);
        } if (mp->T == MIRROR) { return get_color({P,refl(r.u,N)}, depth-1);
        } if (mp->T == REFR) {
            double n = dot(r.u, P-sp->C) > 0 ? mp->n : 1/mp->n;
            double k0 = (1-n)*(1-n) / ((1+n)*(1+n));
            double R = k0 + (1-k0) * pow(1 + dot(N,r.u),5);
            if (uni() < R) return get_color({P,refl(r.u,N)}, depth-1);
            return get_color({P,refr(r.u,N,n)}, depth-1);
        }
    }
};

void init(Scene& scene) {

    scene.spheres = {
        {{5,10,0},4,{}},
        {{-15,10,0},6,{vec{255,113,206}/255}},
        {{-4,-2,-5},8,{1.5,REFR}},
        {{15,0,0},10,{1.5,REFR}},
        {{15,0,0},9.5,{1/1.5,REFR}},
        {{0, 1000,0},940,{vec{1,205,254}/255}},
        {{0,- 955,0},940,{vec{255,155,131}/255}},
        {{0,0, 1000},940,{vec{185,103,255}/255}},
        {{0,0,-1000},940,{vec{5,255,161}/255}},
        {{-10,25,20},10,{1.5e10,LIGHT}},
    };
    
    scene.meshes.push_back( load_cat() );
}

struct Camera {
    Camera(const vec& Q, double alpha) : Q(Q), f(tan(alpha/2)) {}
    vec Q; double f;
};

int main() {
    Scene scene;
    init(scene);
    Camera cam({0,0,55},60*M_PI/180);

    const int W = 600, H = 400;
    const int n_sampl = 2, max_depth = 5;
    double DoF = .5;
    double AA = .2;
    uint8_t img[H][W][3];
#pragma omp parallel for schedule(dynamic,1)
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            vec u(j + .5 - W/2, (H - i - 1) + .5 - H/2, -W/(2*cam.f));
            vec col{0,0,0}; 
            for (int k = 0; k < n_sampl; k++) {
                double x,y; box_muller(AA,x,y);
                vec v = u; v.x += x; v.y += y;
                v = v * (-55/u.z) + cam.Q; // focal plane
                vec O = cam.Q;
                box_muller(DoF,x,y); O.x += x; O.y += y;
                col += scene.get_color({O,norm(v-O)}, max_depth);
            }
            col = col * (1./n_sampl);
            auto M_PImg = img[i][j];
            for (double c: {col.x,col.y,col.z}) { // r,g,b
                c = pow(c,1/2.2); // gamma correct
                *M_PImg++ = c < 0 ? 0 : c > 255 ? 255 : c; // clip
            }
        }
    }
    stbi_write_png("out.png",W,H,3,img,W*3);
}


