#include <iostream>
#include <cstring>
#include "model.hpp"
#include "utils.hpp"

using namespace std::complex_literals;

isope::model::model(const double k0, const double sigma, const double eps, const double v,
                    const size_t nx, const double x0, const double x1,
                    const size_t nz, const double z0, const double z1) :
    a_(1i / (2.0 * k0)), b_(1i * eps * k0 / 2.0), k0_(k0), sigma_(sigma), eps_(eps), v_(v), dx_((x1 - x0) / nx), dz_((z1 - z0) / nz),
    xs_(utils::mesh<rvector1d_t>(x0, x1, nx)), zs_(utils::mesh<rvector1d_t>(z0, z1, nz)),
    k_(vector_k_(x1 - x0, nx))
{
    cA = cvector1d_t(nx);
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
    cvector2d_t sp(n, cvector1d_t(xs_.size(), 0.));
    cvector1d_t spn(xs_.size(), 0.);
    cvector1d_t buff(xs_.size(), 0.);

    fft::fftw fft(static_cast<int>(xs_.size()), true);

    for (size_t i = 0; i < xs_.size(); ++i)
        result[0][i] = as[0][i] = init_cond(xs_[i]);
    std::memcpy(fft.forward_data(), as[0].data(), xs_.size() * sizeof(complex));
    fft.execute_forward();
    std::memcpy(ash[0].data(), fft.forward_data(), xs_.size() * sizeof(complex));

    size_t index = 1;

    for (size_t _ = 0; _ < zs_.size(); ++_) {
        for (size_t i = 0; i < n; ++i) {
            non_linear_coeff(as, 0, fft.forward_data());
            fft.execute_forward();
            utils::multiply_by(fft.forward_data(), fft.forward_data_end(), b_);
            if (_ == 0) std::memcpy(nl[i].data(), fft.backward_data(), xs_.size() * sizeof(complex));
            calculate(ash, fft.forward_data(), nl, buff, i, ash[i].data());
            std::memcpy(nl[i].data(), fft.backward_data(), xs_.size() * sizeof(complex));
            std::memcpy(fft.forward_data(), ash[i].data(), xs_.size() * sizeof(complex));
            fft.execute_backward().normalize_backward();
            for (size_t j = 0; j < xs_.size(); ++j)
                spn[j] = (fft.forward_data(j) - as[i][j]) / dz_;
            for (size_t j = 0; j < xs_.size(); ++j)
                buff[j] = a_ * (spn[j] - expA_[j] * sp[i][j] +
                    fft.forward_data(j)  * (cA[j] * cA[j] * dz_ + cA[j] + 1. / dz_) +
                    as[i][j] * expA_[j] * (cA[j] * cA[j] * dz_ - cA[j] - 1. / dz_ / expA_[j]));
            std::swap(sp[i], spn);
            std::memcpy(as[i].data(), fft.forward_data(), xs_.size() * sizeof(complex));
        }
        if (_ % 40000 == 0)
            std::cout << _ << std::endl;
        if (_ % m == 0) {
            for (size_t i = 1; i < n; ++i)
                for (size_t j = 0; j < xs_.size(); ++j)
                    fft.forward_data(j) += as[n - 1][j];
            for (size_t i = 0; i < xs_.size(); ++i)
                fft.forward_data(i) *= std::exp(1i * k0_ * zs_[_]);
            std::memcpy(result[index++].data(), fft.forward_data(), xs_.size() * sizeof(complex));
        }
    }
    return result;
}

void isope::model::non_linear_coeff(const cvector2d_t& values, const size_t n, complex* data) {
    switch(n) {
        case 0:
            for (size_t i = 0; i < values[0].size(); ++i)
                data[i] = values[0][i] * std::pow(std::abs(values[0][i]), 2);
            return;
        case 1:
            for (size_t i = 0; i < values[0].size(); ++i)
                data[i] = 2 * std::pow(std::abs(values[0][i]), 2) * values[1][i] + std::pow(values[0][i], 2) * std::conj(values[1][i]);
            return;
        default: return;
    }
}

void isope::model::calculate(const cvector2d_t& values, const complex* nl, const cvector2d_t& nlp,
                             const cvector1d_t& buff, const size_t n, complex* data) const {
    switch(n) {
        case 0:
            for (size_t i = 0; i < values[0].size(); ++i)
                data[i] = values[0][i] * expA_[i] + nlfacA_[i] * nl[i] + nlfacAp_[i] * nlp[0][i];
            return;
        case 1:
            for (size_t i = 0; i < values[0].size(); ++i)
                data[i] = values[1][i] * expA_[i] + dz_ * (nl[i] + expA_[i] * nlp[1][i]) + buff[i];
            return;
        default: return;
    }
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
