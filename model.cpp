#include <iostream>
#include <cstring>
#include "model.hpp"
#include "utils.hpp"

using namespace std::complex_literals;

isope::model::model(const double k0, const double sigma, const double eps, const double v,
                    const size_t nx, const double x0, const double x1,
                    const size_t nz, const double z0, const double z1) :
    a_(1i / (2.0 * k0)), b_(1i * eps * k0 / 2.0), k0_(k0), sigma_(sigma), eps_(eps), v_(v),
    xs_(utils::mesh<rvector1d_t>(x0, x1, nx)), zs_(utils::mesh<rvector1d_t>(z0, z1, nz)),
    k_(vector_k_(x1 - x0, nx))
{
    cvector1d_t cA(nx);
    for (size_t i = 0; i < cA.size(); ++i)
        cA[i] = -k_[i] * a_;
    expA_ = cvector1d_t(nx);
    for (size_t i = 0; i < expA_.size(); ++i)
        expA_[i] = std::exp(dt_ * cA[i]);
    nlfacA_ = cvector1d_t(nx);
    for (size_t i = 0; i < nlfacA_.size(); ++i)
        nlfacA_[i] = (expA_[i] * (1. + 1. / cA[i] / dt_) - 1. / cA[i] / dt_ - 2.) / cA[i];
    nlfacAp_ = cvector1d_t(nx);
    for (size_t i = 0; i < nlfacAp_.size(); ++i)
        nlfacAp_[i] = (expA_[i] * (-1. / cA[i] / dt_) + 1. / cA[i] / dt_ + 1.) / cA[i];
    nlfacA_[0] = nlfacAp_[0] = dt_ / 2.;
}

isope::model::cvector2d_t isope::model::solve(std::function<complex(double)> init_cond, const size_t n, const size_t m) const {
    cvector2d_t result(zs_.size() / m + 1, cvector1d_t(xs_.size()));
    cvector2d_t as(n, cvector1d_t(xs_.size(), 0.));
    cvector2d_t ash(n, cvector1d_t(xs_.size(), 0.));
    cvector2d_t nl(n, cvector1d_t(xs_.size()));
    fft::fftw fft(static_cast<int>(xs_.size()), true);

    for (size_t i = 0; i < xs_.size(); ++i)
        result[0][i] = as[0][i] = init_cond(xs_[i]);
    std::memcpy(fft.forward_data(), as[0].data(), xs_.size() * sizeof(complex));
    fft.execute_forward();
    std::memcpy(ash[0].data(), fft.forward_data(), xs_.size() * sizeof(complex));
    size_t index = 1;
    for (size_t i = 0; i < zs_.size(); ++i) {
        non_linear_coeff(as, 0, fft.forward_data());
        fft.execute_forward();
        utils::multiply_by(fft.forward_data(), fft.forward_data_end(), b_);
        if (i == 0)
            std::memcpy(nl[0].data(), fft.backward_data(), xs_.size() * sizeof(complex));
        for (size_t i = 0; i < as[0].size(); ++i)
            ash[0][i] = ash[0][i] * expA_[i] + nlfacA_[i] * fft.forward_data()[i] + nlfacAp_[i] * nl[0][i];
        std::memcpy(nl[0].data(), fft.backward_data(), xs_.size() * sizeof(complex));
        std::memcpy(fft.forward_data(), ash[0].data(), xs_.size() * sizeof(complex));
        fft.execute_backward().normalize_backward();
        std::memcpy(as[0].data(), fft.forward_data(), xs_.size() * sizeof(complex));
        if (i % 40000 == 0)
            std::cout << i << std::endl;
        if (i % m == 0)
            std::memcpy(result[index++].data(), fft.forward_data(), xs_.size() * sizeof(complex));
    }
    return result;
}

void isope::model::non_linear_coeff(const cvector2d_t& values, const size_t n, complex* data) const {
    for (size_t i = 0; i < values[0].size(); ++i)
        data[i] = values[0][i] * std::abs(values[0][i]) * std::abs(values[0][i]);
}

isope::model::rvector1d_t isope::model::vector_k_(const double l, const size_t n) {
    rvector1d_t result(n);
    for (size_t i = 0; i < n / 2; ++i)
        result[i] = i * 2 * M_PI / l;
    result[n / 2] = 0;
    auto t = -static_cast<int64_t>(n) / 2 + 1;
    for (auto i = n / 2 + 1; i < n; ++i, ++t)
        result[i] = t * 2 * M_PI / l;
    for (auto &it : result)
        it *= it;
    result[n / 2] = (2 * M_PI / l) * n / 2 * (2 * M_PI / l) * n / 2;
    return result;
}
