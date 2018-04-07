#include <iostream>
#include <cstring>
#include "model.hpp"
#include "utils.hpp"

using namespace std::complex_literals;

isope::model::model(const double k0, const double sigma, const double eps, const double v,
                    const size_t nx, const double x0, const double x1,
                    const size_t nz, const double z0, const double z1) :
    a_(1i / (2. * k0)), b_(1i * eps * k0 / 2.), k0_(k0), sigma_(sigma), eps_(eps), v_(v), dx_((x1 - x0) / nx), dz_((z1 - z0) / nz),
    xs_(utils::mesh<rvector1d_t>(x0, x1, nx)), zs_(utils::mesh<rvector1d_t>(z0, z1, nz)),
    k_(vector_k_(x1 - x0, nx))
{
    cA_ = cvector1d_t(nx);
    for (size_t i = 0; i < cA_.size(); ++i)
        cA_[i] = -k_[i] * a_;
    expA_ = cvector1d_t(nx);
    for (size_t i = 0; i < expA_.size(); ++i)
        expA_[i] = std::exp(dz_ * cA_[i]);
    nlfacA_ = cvector1d_t(nx);
    for (size_t i = 0; i < nlfacA_.size(); ++i)
        nlfacA_[i] = (expA_[i] * (1. + 1. / cA_[i] / dz_) - 1. / cA_[i] / dz_ - 2.) / cA_[i];
    nlfacAp_ = cvector1d_t(nx);
    for (size_t i = 0; i < nlfacAp_.size(); ++i)
        nlfacAp_[i] = (expA_[i] * (-1. / cA_[i] / dz_) + 1. / cA_[i] / dz_ + 1.) / cA_[i];
    nlfacA_[0] = nlfacAp_[0] = dz_ / 2.;
}

isope::model::cvector2d_t isope::model::solve(std::function<complex(double, double)> init_cond, const size_t n, const size_t m) const {
    cvector2d_t result(zs_.size() / m + 1, cvector1d_t(xs_.size()));
    cvector2d_t as(n, cvector1d_t(xs_.size(), 0.));
    cvector2d_t ash(n, cvector1d_t(xs_.size(), 0.));
    cvector2d_t nl(n, cvector1d_t(xs_.size()));
    cvector2d_t sp(n, cvector1d_t(xs_.size(), 0.));
    cvector1d_t spn(xs_.size(), 0.);
    cvector1d_t buff(xs_.size(), 0.), buff1(xs_.size(), 0.);

    fft::fftw fft(static_cast<int>(xs_.size()), true);


//    for (size_t i = 0; i < xs_.size(); ++i)
//        fft.forward_data(i) = der(xs_[i], zs_[0]);
//
//    fft.execute_forward();
//    std::memcpy(sp[0].data(), fft.backward_data(), xs_.size() * sizeof(complex));

    for (size_t i = 0; i < xs_.size(); ++i)
        result[0][i] = as[0][i] = init_cond(xs_[i], 0.);
    std::memcpy(fft.forward_data(), as[0].data(), xs_.size() * sizeof(complex));
    fft.execute_forward();
    std::memcpy(ash[0].data(), fft.forward_data(), xs_.size() * sizeof(complex));

    size_t index = 1;

    for (size_t _ = 0; _ < zs_.size(); ++_) {
        for (size_t i = 0; i < n; ++i) {
            non_linear_coeff(as, i, fft.forward_data());
            fft.execute_forward();
            utils::multiply_by(fft.forward_data(), fft.forward_data_end(), b_);
            if (_ == 0) std::memcpy(nl[i].data(), fft.backward_data(), xs_.size() * sizeof(complex));
            calculate(ash, fft.forward_data(), nl, buff, i, buff1.data());
//            if (_ == 0 && i != 0) std::memcpy(sp[i].data(), buff1.data(), xs_.size() * sizeof(complex));
            std::memcpy(nl[i].data(), fft.backward_data(), xs_.size() * sizeof(complex));
            for (size_t j = 0; j < xs_.size(); ++j)
                spn[j] = (buff1[j] - ash[i][j]) / dz_;
            for (size_t j = 0; j < xs_.size(); ++j)
                buff[j] = spn[j] - expA_[j] * sp[i][j] + cA_[i] * (
                        buff1[j] * (1. + cA_[i] * dz_ / 2.) + ash[i][j] * expA_[j] * (cA_[i] * dz_ / 2. - 1.));
            std::swap(sp[i], spn);
            std::memcpy(ash[i].data(), buff1.data(), xs_.size() * sizeof(complex));
            std::memcpy(fft.forward_data(), buff1.data(), xs_.size() * sizeof(complex));
            fft.execute_backward().normalize_backward();
            std::memcpy(as[i].data(), fft.forward_data(), xs_.size() * sizeof(complex));
        }
        if (_ % 40000 == 0)
            std::cout << _ << std::endl;
        if (_ % m == 0) {
            for (size_t i = 2; i <= n; ++i)
                for (size_t j = 0; j < xs_.size(); ++j)
                    fft.forward_data(j) += as[n - i][j];
            for (size_t i = 0; i < xs_.size(); ++i)
                fft.forward_data(i) *= std::exp(1i * k0_ * zs_[_ + 1]);
            std::memcpy(result[index++].data(), fft.forward_data(), xs_.size() * sizeof(complex));
        }
    }
    return result;
}

void isope::model::non_linear_coeff(const cvector2d_t& values, const size_t n, complex* data) {
    for (size_t i = 0; i <= n; ++i)
        for (size_t j = 0; j <= n - i; ++j)
            for (size_t k = 0; k < values[0].size(); ++k)
                data[k] = values[i][k] * values[j][k] * std::conj(values[n - i - j][k]);
//    switch(n) {
//        case 0:
//            for (size_t i = 0; i < values[0].size(); ++i)
//                data[i] = values[0][i] * std::pow(std::abs(values[0][i]), 2);
//            return;
//        case 1:
//            for (size_t i = 0; i < values[0].size(); ++i)
//                data[i] = 2 * std::pow(std::abs(values[0][i]), 2) * values[1][i] + std::pow(values[0][i], 2) * std::conj(values[1][i]);
//            return;
//        default: return;
//    }
}

void isope::model::calculate(const cvector2d_t& values, const complex* nl, const cvector2d_t& nlp,
                             const cvector1d_t& buff, const size_t n, complex* data) const {
    if (n == 0)
        for (size_t i = 0; i < values[0].size(); ++i)
            data[i] = values[0][i] * expA_[i] + nlfacA_[i] * nl[i] + nlfacAp_[i] * nlp[0][i];
    else
        for (size_t i = 0; i < values[0].size(); ++i)
            data[i] = values[n][i] * expA_[i] + dz_ / 2. * (nl[i] + expA_[i] * nlp[n][i]) + a_ * buff[i];
}

isope::model::rvector1d_t isope::model::vector_k_(const double l, const size_t n) {
    rvector1d_t result(n);
    for (size_t i = 0; i < n / 2; ++i)
        result[i] = i * 2 * M_PI / l;
    auto t = -static_cast<int64_t>(n) / 2 + 1;
    for (auto i = n / 2 + 1; i < n; ++i, ++t)
        result[i] = t * 2 * M_PI / l;
    for (auto &it : result)
        it *= it;
    result[n / 2] = (2 * M_PI / l) * n / 2 * (2 * M_PI / l) * n / 2;
    return result;
}

isope::model::complex isope::model::a0sol(const double x, const double z) const {
    return sigma_ / k0_ * std::sqrt(2 / eps_) / std::cosh(sigma_ * x) * std::exp(1i * std::pow(sigma_ , 2) / 2. / k0_ * z);
}

isope::model::complex isope::model::a1sol(const double x, const double z) const {
    return -1i * std::pow(sigma_, 5) / 8. / std::pow(k0_, 4) * z * std::sqrt(2 / eps_) / std::cosh(sigma_ * x) * std::exp(1i * std::pow(sigma_ , 2) / 2. / k0_ * z);
}

isope::model::complex isope::model::der(const double x, const double z) const {
    return std::exp(1i * std::pow(sigma_, 2) * z / 2. / k0_) * std::pow(eps_, -0.5) * std::pow(sigma_, 5) * (std::pow(sigma_, 2) * z - 2i * k0_) / std::cosh(sigma_ * x) / 8. / std::sqrt(2) / std::pow(k0_, 5);
}