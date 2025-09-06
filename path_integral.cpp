#include <iostream>
#include <random>
#include <cmath>
int main() {
    double a = -3.0;
    double b = 3.0;
    double dy = 0.125;
    double m = 1.0;
    double w = 1.0;
    double dt = 1.0;
    double lambda = 0.0;
    int nt = 4;
    int ntraj = 10000000;
    int nbin = (b-a)/dy;
    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<> uni(a,b);
    double x[nt], k[nbin];
    for (int bin = 0; bin < nbin; ++bin) {
        k[bin] = 0;
    }
    for (int traj = 0; traj < ntraj; ++traj) {
        for (int t = 0; t < nt; ++t) {
            x[t] = uni(rng);
        }
        double S = 0;
        for (int t = 0; t < nt; ++t) {
            S += m*pow((x[(t+1)%nt]-x[t])/dt,2)/2+m*w*w*x[t]*x[t]/2+lambda*x[t]*x[t]*x[t]*x[t];
        }
        int i = 0;
        for (double y = a; y <= b; y += dy) {
            if (x[0] > y && x[0] < y+dy) {
                k[i] += pow((b-a),nt)*exp(-dt*S)/ntraj;
            }
            ++i;
        }
    }
    double sum = 0;
    for (int bin = 0; bin < nbin; ++bin) {
        sum += k[bin]*dy;
    }
    int i = 0;
    for (double y = a; y < 0; y += dy) {
        std::cout << y << '\t' << k[i]/sum << '\n';
        ++i;
    }
    for (double y = dy; y <= b; y += dy) {
        std::cout << y << '\t' << k[i]/sum << '\n';
        ++i;
    }
}
