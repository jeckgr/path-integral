#include <iostream>
#include <random>
#include <cmath>
int main() {
    int a = -6;
    int b = 6;
    int nt = 4;
    int npath = 1000000000;
    double dt = 0.25;
    double m2 = 1.0;
    double lambda = 0.0;
    std::mt19937 rng(rand());
    std::uniform_real_distribution<> uni(a,b);
    double x[nt];
    double c[nt];
    for (int i = 0; i < nt; ++i) {
        c[i] = 0;
    }
    double z = 0;
    for (int i = 0; i < npath; ++i) {
        for (int j = 0; j < nt; ++j) {
            x[j] = uni(rng);
        }
        double S = 0;
        for (int j = 0; j < nt; ++j) {
            S += pow((x[(j+1)%nt]-x[j])/dt,2)/2+m2*x[j]*x[j]/2+lambda*x[j]*x[j]*x[j]*x[j]/24;
        }
        for (int j = 0; j < nt; ++j) {
            c[j] += pow((b-a),nt)*x[j]*x[0]*exp(-dt*S)/npath;
        }
        z += pow((b-a),nt)*exp(-dt*S)/npath;
    }
    for (int i = 1; i < nt; ++i) {
        std::cout << acosh((c[(i+1)%nt]/z+c[(i-1+nt)%nt]/z)/(2*c[i]/z))/dt << '\n';
    }
}
