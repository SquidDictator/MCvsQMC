#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include "sobol.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2 0.70710678118654752440 // 1/sqrt(2)
#endif


// ----------- Black-Scholes European Call -----------
double bs_call(double S, double K, double r, double sigma, double T) {
    double d1 = (std::log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    double d2 = d1 - sigma * std::sqrt(T);
    auto N = [](double x) { return 0.5 * std::erfc(-x * M_SQRT1_2); };
    return S * N(d1) - K * std::exp(-r * T) * N(d2);
}

// ----------- Approximate Inverse Error Function (erfinv) -----------
double erfinv(double x) {
    double tt1, tt2, lnx, sgn;
    sgn = (x < 0) ? -1.0 : 1.0;
    x = (1 - x) * (1 + x);
    lnx = std::log(x);
    tt1 = 2 / (M_PI * 0.147) + 0.5 * lnx;
    tt2 = 1 / 0.147 * lnx;
    return (sgn * std::sqrt(-tt1 + std::sqrt(tt1 * tt1 - tt2)));
}

namespace std {
    double erfinv(double x) { return ::erfinv(x); }
}

// ----------- Metrics Functions -----------

// Pi estimation: unit quarter circle in [0,1]^2
double estimatePi(const std::vector<std::pair<double, double>>& pts) {
    int inside = 0;
    for (const auto& p : pts)
        if (p.first * p.first + p.second * p.second <= 1.0) inside++;
    return 4.0 * inside / pts.size();
}

// 2D integral of f(x, y) = x*y over [0,1]^2
double estimateXY(const std::vector<std::pair<double, double>>& pts) {
    double sum = 0.0;
    for (const auto& p : pts) sum += p.first * p.second;
    return sum / pts.size();
}

// 1D mean of exp(x) over [0,1]
double estimateExpX(const std::vector<double>& xs) {
    double sum = 0.0;
    for (double x : xs) sum += std::exp(x);
    return sum / xs.size();
}

// European option pricing by simulation
double estimate_euro_call(const std::vector<double>& zs,
                          double S0, double K, double r, double sigma, double T)
{
    double sum = 0.0;
    for (double Z : zs) {
        double ST = S0 * std::exp((r - 0.5 * sigma * sigma) * T + sigma * std::sqrt(T) * Z);
        double payoff = std::max(ST - K, 0.0);
        sum += std::exp(-r * T) * payoff;
    }
    return sum / zs.size();
}

// Generate standard normal samples from uniform [0,1] via inverse CDF
std::vector<double> uniform_to_normal(const std::vector<double>& u) {
    std::vector<double> zs;
    zs.reserve(u.size());
    for (double x : u)
        zs.push_back(std::sqrt(2) * std::erfinv(2 * x - 1));
    return zs;
}

int main() {
    // Option Parameters (can be changed)
    double S0 = 100.0, K = 100.0, r = 0.05, sigma = 0.5, T = 1.0; //sigma 0.5 - market far more volatile
    double analytic_call = bs_call(S0, K, r, sigma, T);

    // For error checking
    const double XY_exact = 0.25;
    const double expx_exact = std::exp(1.0) - 1.0; // mean = (e^1 - 1)/1 = e-1 ≈ 1.7182818

    // Sample Sizes
    std::vector<int> ns = {1000, 5000, 10000, 50000, 100000};

    // Output CSV
    std::ofstream fout("mc_vs_qmc_allmetrics.csv");
    fout << "N,Method,Pi,PiError,XY,XYError,ExpX,ExpXError,CallOption,CallOptionError,AnalyticCall\n";
    std::cout << std::fixed << std::setprecision(8);

    for (int n : ns) {
        // Monte Carlo: random
        std::mt19937 rng(42 + n); // different seed for each n
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        std::normal_distribution<double> norm(0.0, 1.0);

        // For 2D
        std::vector<std::pair<double, double>> pts2d;
        pts2d.reserve(n);
        for (int i = 0; i < n; ++i)
            pts2d.emplace_back(dist(rng), dist(rng));

        // For 1D
        std::vector<double> pts1d;
        pts1d.reserve(n);
        for (int i = 0; i < n; ++i)
            pts1d.push_back(dist(rng));

        // For option pricing (normal)
        std::vector<double> zs;
        zs.reserve(n);
        for (int i = 0; i < n; ++i)
            zs.push_back(norm(rng));

        double pi_mc = estimatePi(pts2d);
        double xy_mc = estimateXY(pts2d);
        double expx_mc = estimateExpX(pts1d);
        double call_mc = estimate_euro_call(zs, S0, K, r, sigma, T);

        fout << n << ",MC,"
             << pi_mc << "," << std::abs(pi_mc - M_PI) << ","
             << xy_mc << "," << std::abs(xy_mc - XY_exact) << ","
             << expx_mc << "," << std::abs(expx_mc - expx_exact) << ","
             << call_mc << "," << std::abs(call_mc - analytic_call) << ","
             << analytic_call << "\n";

        // Quasi Monte Carlo: Sobol
        // 2D
        unsigned long dim2d = 2;
        sobol::SobolGenerator sobolGen2(dim2d);
        std::vector<std::pair<double, double>> sobol2d;
        std::vector<double> sobol_pt(dim2d);
        for (int i = 0; i < n; ++i) {
            sobolGen2.nextPoint(sobol_pt);
            sobol2d.emplace_back(sobol_pt[0], sobol_pt[1]);
        }
        // 1D
        unsigned long dim1d = 1;
        sobol::SobolGenerator sobolGen1(dim1d);
        std::vector<double> sobol1d;
        std::vector<double> pt(dim1d);
        for (int i = 0; i < n; ++i) {
            sobolGen1.nextPoint(pt);
            sobol1d.push_back(pt[0]);
        }
        // For option pricing: need normal draws from Sobol
        std::vector<double> sobol_normal = uniform_to_normal(sobol1d);

        double pi_qmc = estimatePi(sobol2d);
        double xy_qmc = estimateXY(sobol2d);
        double expx_qmc = estimateExpX(sobol1d);
        double call_qmc = estimate_euro_call(sobol_normal, S0, K, r, sigma, T);

        fout << n << ",QMC,"
             << pi_qmc << "," << std::abs(pi_qmc - M_PI) << ","
             << xy_qmc << "," << std::abs(xy_qmc - XY_exact) << ","
             << expx_qmc << "," << std::abs(expx_qmc - expx_exact) << ","
             << call_qmc << "," << std::abs(call_qmc - analytic_call) << ","
             << analytic_call << "\n";
    }

    fout.close();
    std::cout << "Results saved to mc_vs_qmc_allmetrics.csv\n";
    return 0;
}
