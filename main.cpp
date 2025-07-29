#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <cmath>
#include <iomanip>

// Simple Sobol sequence generator for 2D (enough for this demo)
void sobol2D(int n, std::vector<std::pair<double, double>>& points) {
    // Van der Corput for base 2 and 3
    auto vdc = [](int n, int base) {
        double v = 0, denom = 1;
        while (n) {
            denom *= base;
            v += (n % base) / denom;
            n /= base;
        }
        return v;
    };
    points.clear();
    for (int i = 0; i < n; ++i)
        points.emplace_back(vdc(i+1, 2), vdc(i+1, 3));
}

// Monte Carlo: uniform random samples in 2D
void monteCarlo2D(int n, std::vector<std::pair<double, double>>& points) {
    std::mt19937 rng(42); // fixed seed for reproducibility
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    points.clear();
    for (int i = 0; i < n; ++i)
        points.emplace_back(dist(rng), dist(rng));
}

// Pi estimation (unit circle in unit square)
double estimatePi(const std::vector<std::pair<double, double>>& pts) {
    int inside = 0;
    for (const auto& p : pts) {
        if (p.first*p.first + p.second*p.second <= 1.0) inside++;
    }
    return 4.0 * inside / pts.size();
}

// 2D integral of f(x, y) = x*y over [0,1]^2, exact = 0.25
double estimateXY(const std::vector<std::pair<double, double>>& pts) {
    double sum = 0.0;
    for (const auto& p : pts) sum += p.first * p.second;
    return sum / pts.size();
}

// 1D mean of f(x) = exp(x) over [0,1], exact = (e - 1)/1 ≈ 1.71828
double estimateExpX(const std::vector<double>& xs) {
    double sum = 0.0;
    for (double x : xs) sum += std::exp(x);
    return sum / xs.size();
}

// For QMC 1D, just use first coordinate of Sobol/Van der Corput
void sobol1D(int n, std::vector<double>& xs) {
    auto vdc = [](int n, int base) {
        double v = 0, denom = 1;
        while (n) {
            denom *= base;
            v += (n % base) / denom;
            n /= base;
        }
        return v;
    };
    xs.clear();
    for (int i = 0; i < n; ++i) xs.push_back(vdc(i+1, 2));
}
void monteCarlo1D(int n, std::vector<double>& xs) {
    std::mt19937 rng(123); // different seed
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    xs.clear();
    for (int i = 0; i < n; ++i) xs.push_back(dist(rng));
}

int main() {
    std::ofstream fout("mc_vs_qmc_results.csv");
    fout << "N,Method,PiEstimate,XYIntegral,ExpXMean\n";
    std::vector<int> ns = {100, 500, 1000, 5000, 10000};
    for (int n : ns) {
        // Monte Carlo
        std::vector<std::pair<double,double>> pts2d;
        std::vector<double> pts1d;
        monteCarlo2D(n, pts2d);
        monteCarlo1D(n, pts1d);
        double pi_mc = estimatePi(pts2d);
        double xy_mc = estimateXY(pts2d);
        double expx_mc = estimateExpX(pts1d);
        fout << n << ",MC," << std::setprecision(8) << pi_mc << "," << xy_mc << "," << expx_mc << "\n";

        // Quasi Monte Carlo
        sobol2D(n, pts2d);
        sobol1D(n, pts1d);
        double pi_qmc = estimatePi(pts2d);
        double xy_qmc = estimateXY(pts2d);
        double expx_qmc = estimateExpX(pts1d);
        fout << n << ",QMC," << std::setprecision(8) << pi_qmc << "," << xy_qmc << "," << expx_qmc << "\n";
    }
    fout.close();
    std::cout << "Results saved to mc_vs_qmc_results.csv\n";
    return 0;
}
