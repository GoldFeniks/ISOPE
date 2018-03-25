#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <cstring>
#include <iomanip>
#include "model.hpp"
#include "fft.hpp"
#include "utils.hpp"

using namespace std::complex_literals;

int main() {
    const double k0 = 2 * M_PI / 1e-6 * 2, sigma = std::sqrt(2) * 1 / 1.0049e-6, eps = 2 * 2e-17 / 2, v = std::sqrt(2) * 0.01 / 1.0049e-6;
    const auto pm = sigma / k0 * std::sqrt(2 / eps);
    isope::model m(k0, sigma, eps, v, 256, -4e-05, 3.96875e-05,
                   static_cast<size_t>(std::round(2e-3 / 0.5e-8)), 0, 1);
    auto a = m.solve([=](double x) -> isope::fft::basic_fft::complex {
        return pm / std::cosh(sigma * x + sigma * v / k0 * 0.) * std::exp(-1.0i * v * x + (1.0i * (sigma * sigma - v * v) / k0 / 2.0) * 0.);
    }, 1, static_cast<size_t>(std::round(2e-3 / 0.5e-8)) / 1000);
    auto size = a.size();
    std::ofstream out("test.txt");
    for (const auto& i : a) {
        for (const auto& j : i)
            out << std::setprecision(6) << j.real() << ' ';
        out << std::endl;
    }
    return 0;
}