////////////////////////////////////////////////////////////////////////////////
// benchmark.cpp                                                              //
// ------                                                                     //
// benchmark                                                                  //
////////////////////////////////////////////////////////////////////////////////

#include <mlinterp>

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>

using namespace mlinterp;
using namespace std;

double func(double x, double y, double z) { return sin(x) * cos(y) * sin(exp(0.5 * z)) + 10; }

int main() {

    // Boundaries of the interval [-pi, pi]
    constexpr double b = M_PI, a = -b;

    constexpr int nxd = 15, nyd = 15, nzd = 15, nd[] = {nxd, nyd, nzd};
    double xd[nxd];
    for (int i = 0; i < nxd; ++i) {
        xd[i] = a + (b - a) / (nxd - 1) * i;
    }
    double yd[nyd];
    for (int j = 0; j < nyd; ++j) {
        yd[j] = a + (b - a) / (nyd - 1) * j;
    }
    double zd[nzd];
    for (int j = 0; j < nzd; ++j) {
        zd[j] = a + (b - a) / (nzd - 1) * j;
    }

    double fd[nxd * nyd * nzd];
    for (int i = 0; i < nxd; ++i) {
        for (int j = 0; j < nyd; ++j) {
            for (int k = 0; k < nzd; ++k) {
                const int n = i * nyd * nzd + j * nzd + k;
                fd[n] = func(xd[i], yd[j], zd[k]);
            }
        }
    }

    // interpolate points
    constexpr int ni = 500000;
    std::mt19937 gen;
    std::uniform_real_distribution<> dis(a, b);
    std::vector<double> xi(ni);
    std::vector<double> yi(ni);
    std::vector<double> zi(ni);
    for (int i = 0; i < ni; i++) {
        xi[i] = dis(gen);
        yi[i] = dis(gen);
        zi[i] = dis(gen);
    }
    std::vector<double> fi(ni);

    // Perform the interpolation
    auto t1 = std::chrono::high_resolution_clock::now();
    interp(nd, ni,        // Number of points
           fd, fi.data(), // Output axis (z)
           xd, xi.data(), //
           yd, yi.data(), //
           zd, zi.data()  // Input axes (x,y,z)
    );
    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = t2 - t1;
    std::cout << "Time:" << diff.count() << " s\n";

    // Print the interpolated values
    cout << scientific << setprecision(8) << showpos;
    for (int n = 0; n < ni; ++n) {
        const double val = func(xi[n], yi[n], zi[n]);
        cout << fi[n] << "\t" << (fi[n] - val) / val << endl;
    }

    return 0;
}
