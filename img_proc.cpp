#include <bits/stdc++.h>
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include "vec.cpp"
#include "random.cpp"

vec from_char(unsigned char* a) { return {a[0], a[1], a[2]}; }
void from_vec(const vec& v, unsigned char* a) {
    auto clamp = [](double x)->unsigned char
    { return x < 0 ? 0 : x > 255 ? 255 : x; };
    a[0] = clamp(v.x);
    a[1] = clamp(v.y);
    a[2] = clamp(v.z);
}

void SOTCT(std::vector<vec>& I, const std::vector<vec>& M) {
    int n = I.size();
    if (M.size() != n) { std::cerr << "not same size!\n"; exit(1); }
    for (int i = 0; i < 100; i++) {
        vec v = uni_sphere();
        std::vector<std::pair<double,int>> vi(n), vm(n);
        for (int i = 0; i < n; i++)
            vi[i] = { dot(v,I[i]), i };
        for (int i = 0; i < n; i++)
            vm[i] = { dot(v,M[i]), i };
        sort(vi.begin(), vi.end());
        sort(vm.begin(), vm.end());
        for (int i = 0; i < n; i++)
            I[vi[i].second] += (vm[i].first - vi[i].first) * v;
    }
}

int main() {
    int w,h,nc;
    auto img = stbi_load("renders/bigg.png", &w, &h, &nc, 3);
    auto img2 = stbi_load("renders/balls.png", &w, &h, &nc, 3);
    std::vector<vec> im1(w*h), im2(w*h);
    for (int i = 0; i < w*h; i++) im1[i] = from_char(img+3*i);
    for (int i = 0; i < w*h; i++) im2[i] = from_char(img2+3*i);

    SOTCT(im1, im2);

    for (int i = 0; i < w*h; i++) from_vec(im1[i], img+3*i);
    stbi_write_png("out.png", w, h, 3, img, w*3);
}
