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

template<typename T, typename S, typename It = T*, typename F = std::function<T(const T&)>>
void export_to_stream(It begin, const It end, S& stream, const F& func = [](const T& value) { return value; }) {
    while(begin != end)
        stream << func(*begin++);
}

int main() {
    const double k0 = 2 * M_PI / 1e-6 * 2, sigma = std::sqrt(2) * 1 / 1.0049e-6, eps = 2 * 2e-17 / 2, v = std::sqrt(2) * 0. / 1.0049e-6;
    const auto pm = sigma / k0 * std::sqrt(2 / eps);
    isope::cap_model m(k0, sigma, eps, v, 256, -4e-05, 4e-05, 1200000
                   /*static_cast<size_t>(std::round(2e-3 / 0.5e-8))*/, 0, 0.006, 2e-11);
    auto t1 = std::chrono::system_clock::now();
    auto a = m.solve([=](double x, double z = 0) -> isope::fft::basic_fft::complex {
        return pm / std::cosh(sigma * x + sigma * v / k0 * z) * std::exp(-1.0i * v * x + (1.0i * (sigma * sigma - v * v) / k0 / 2.0) * z);
    }, 1, static_cast<size_t>(std::round(2e-3 / 0.5e-8)) / 1000);
    auto t2 = std::chrono::system_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << std::endl;
    auto size = a.size();
    std::ofstream out("tests/test3.bin", std::ios_base::binary);
    export_to_stream<isope::fft::basic_fft::complex>(a[0].data(), a[0].data() + a.size() * a[0].size() * sizeof(a[0][0]), out);
//    for (const auto& i : a)
//        out.write(reinterpret_cast<const char*>(i.data() + 128), sizeof(isope::fft::basic_fft::complex) * i.size() / 2);
    out.close();
    out.open("tests/test3ABS.bin", std::ios_base::binary);
    export_to_stream<isope::types::complex>(a[0].data(), a[0].data() + a.size() * a[0].size() * sizeof(a[0][0]), out, [](const isope::fft::basic_fft::complex& value) { return std::abs(value); });
//    for (const auto& i : a)
//        for (const auto& j : i) {
//            auto buff = std::abs(j);
//            out.write(reinterpret_cast<char*>(&buff), sizeof(double));
//        }
    out.close();
//    out.open("testINFO.txt");
//    out << m.dx() << ' ' << m.dz() * 400 << ' ' << a[0].size() << ' ' << a.size() << std::endl;
//    out.close();
    out.open("testA1C.txt");
    for (const auto& i : a) {
        for (const auto& j : i)
            out << std::setprecision(6) << j.real() << ' ';
        out << std::endl;
    }
    return 0;
}