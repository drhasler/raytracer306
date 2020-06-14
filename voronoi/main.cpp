#include "main.h"

int main() {
    int N = 1000; // watr molecules
    int M = 1000; // air molecules
    double f = .3; // watr area
    simulate(N,M,f); // fluid simulation

    // test_powercell(N); // food in the center
}
