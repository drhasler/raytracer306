#include "main.h"
#include "svg_writer.cpp"

int main() {
    const int N = 100;
    auto pts = sample_uni2(N);
    auto lam = sample_uni1(N);

    Diagram partition;
    for (int i = 0; i < N; i++) {
        Cell cell(pts[i], lam[i]);
        for (int j = 0; j < N; j++)
            if (i != j) cell.add_pt(pts[j], lam[j]);
        partition.push_back(cell);
    }

    save_svg("out.svg", partition);
}
