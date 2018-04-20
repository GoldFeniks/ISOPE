#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <cstring>
#include <iomanip>
#include "model.hpp"
#include "fft.hpp"
#include "utils.hpp"
#include <chrono>

using namespace std::complex_literals;

int main() {
    const double k0 = 2 * M_PI / 1e-6 * 2, sigma = std::sqrt(2) * 1 / 1.0049e-6, eps = 2 * 2e-17 / 2, v = std::sqrt(2) * 0. / 1.0049e-6;
    const auto pm = sigma / k0 * std::sqrt(2 / eps);
    isope::model m(k0, sigma, eps, v, 2048, -32e-05, 32e-05, 400000
                   /*static_cast<size_t>(std::round(2e-3 / 0.5e-8))*/, 0, 0.002);
    auto t1 = std::chrono::system_clock::now();
    auto a = m.solve([=](double x, double z = 0) -> isope::fft::basic_fft::complex {
        return pm / std::cosh(sigma * x + sigma * v / k0 * z) * std::exp(-1.0i * v * x + (1.0i * (sigma * sigma - v * v) / k0 / 2.0) * z);
    }, 0, static_cast<size_t>(std::round(2e-3 / 0.5e-8)) / 1000);
    auto t2 = std::chrono::system_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << std::endl;
    auto size = a.size();
    std::ofstream out("testBIN.bin", std::ios_base::binary);
    for (const auto& i : a)
        out.write(reinterpret_cast<const char*>(i.data()), sizeof(isope::fft::basic_fft::complex) * i.size());
    out.close();
    out.open("testBINABS.bin", std::ios_base::binary);
    for (const auto& i : a)
        for (const auto& j : i) {
            auto buff = std::abs(j);
            out.write(reinterpret_cast<char*>(&buff), sizeof(double));
        }
    out.close();
//    out.open("testINFO.txt");
//    out << m.dx() << ' ' << m.dz() * 400 << ' ' << a[0].size() << ' ' << a.size() << std::endl;
    out.close();
    out.open("testA1C.txt");
    for (const auto& i : a) {
        for (const auto& j : i)
            out << std::setprecision(6) << j.real() << ' ';
        out << std::endl;
    }
    return 0;
}