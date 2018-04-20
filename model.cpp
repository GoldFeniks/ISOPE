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
    k_(vector_k_(x1 - x0, nx)), x_size_(xs_.size()), z_size_(zs_.size()), bytes_size_(xs_.size() * sizeof(complex))
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

void isope::model::diff(const cvector1d_t& ash0, const cvector1d_t& ash1, cvector1d_t& sp) const {
    for (size_t i = 0; i < ash0.size(); ++i)
        sp[i] = (ash1[i] - ash0[i]) / 2. / dz_;
}


void isope::model::store(cvector2d_t& result, std::vector<size_t>& indexes, const cvector1d_t& values,
                         const size_t index, const size_t z_index, const bool do_store) const {
    if (!do_store) return;
    auto& ind = indexes[index];
    for (size_t i = 0; i < values.size(); ++i)
        result[ind][i] += values[i] * std::exp(1i * k0_ * zs_[z_index]);
    ++ind;
}

isope::model::cvector2d_t
isope::model::new_solve(std::function<complex(double, double)> init_cond, const size_t n, const size_t m) const {

    size_t index = 1;

    cvector2d_t result(zs_.size() / m + 1, cvector1d_t(x_size_, 0.)),
                as(n + 1, cvector1d_t(x_size_, 0.)),
                ash0(n + 1, cvector1d_t(x_size_, 0.)),
                ash1(n + 1, cvector1d_t(x_size_, 0.)),
                nl(n + 1, cvector1d_t(x_size_, 0.)),
                sp(n + 1, cvector1d_t(x_size_, 0.)),
                dd0(x_size_, cvector1d_t(x_size_, 0.)),
                dd1(x_size_, cvector1d_t(x_size_, 0.));

    fft::fftw fft(static_cast<int>(x_size_), true);

    std::vector<size_t> indexes(n + 1, 1);

    for (size_t i = 0; i < x_size_; ++i)
        result[0][i] = as[0][i] = init_cond(xs_[i], zs_[0]);

    std::memcpy(fft.forward_data(), as[0].data(), bytes_size_);
    fft.execute_forward();
    std::memcpy(ash1[0].data(), fft.backward_data(), bytes_size_);

    step<flags::first_call | flags::first_equation>(0, as, ash0[0], ash1[0], sp[0], dd0[0], dd1[0], nl[0], fft);
    for (size_t i = 1; i <= n; ++i) {
        step(0, as, ash0[0], ash1[0], sp[0], dd0[0], dd1[0], nl[0], fft);
        store(result, indexes, as[0], 0, i, i % m == 0);
        for (size_t j = 1; j < i; ++j) {
            step(j, as, ash0[j], ash1[j], sp[j], dd0[j - 1], dd1[j], nl[j], fft);
            store(result, indexes, as[j], j, i - j, (i - j) % m == 0);
        }
        step<flags::first_call>(i, as, ash0[i], ash1[i], sp[i], dd0[i], dd1[i], nl[i], fft);
        std::swap(dd0, dd1);
    }

    for (size_t i = n + 1; i < z_size_ - n; ++i) {
        step(0, as, ash0[0], ash1[0], sp[0], dd0[0], dd1[0], nl[0], fft);
        store(result, indexes, as[0], 0, i, i % m == 0);
        for (size_t j = 1; j <= n; ++j) {
            step(j, as, ash0[j], ash1[j], sp[j], dd0[j - 1], dd1[j], nl[j], fft);
            store(result, indexes, as[j], j, i - j, (i - j) % m == 0);
        }
        std::swap(dd0, dd1);
        if (i % 40000 == 0)
            std::cout << i << std::endl;
    }

    for (size_t i = 1; i <= n; ++i) {
        for (size_t j = n; j >= i; --j) {
            step(j, as, ash0[j], ash1[j], sp[j], dd0[j - 1], dd1[j], nl[j], fft);
            store(result, indexes, as[j], j, i - j, (i - j) % m == 0);
        }
        std::swap(dd0, dd1);
    }
    return result;
}



isope::model::cvector2d_t isope::model::solve(std::function<complex(double, double)> init_cond, const size_t n, const size_t m) const {
    cvector2d_t result(zs_.size() / m + 1, cvector1d_t(xs_.size()));
    cvector2d_t as(n, cvector1d_t(xs_.size(), 0.));
    cvector2d_t ash(n, cvector1d_t(xs_.size(), 0.));
    cvector2d_t ashp(n, cvector1d_t(xs_.size(), 0.));
    cvector2d_t nl(n, cvector1d_t(xs_.size()));
    cvector2d_t sp(n, cvector1d_t(xs_.size(), 0.));
    cvector1d_t spn(xs_.size(), 0.);
    cvector1d_t buff(xs_.size(), 0.), buff1(xs_.size(), 0.);

    fft::fftw fft(static_cast<int>(xs_.size()), true);

    for (size_t i = 0; i < xs_.size(); ++i)
        result[0][i] = as[0][i] = init_cond(xs_[i], 0.);
    std::memcpy(fft.forward_data(), as[0].data(), xs_.size() * sizeof(complex));
    fft.execute_forward();
    std::memcpy(ash[0].data(), fft.forward_data(), xs_.size() * sizeof(complex));

    size_t index = 1;

    for (size_t _ = 0; _ < zs_.size(); ++_) {
        for (size_t i = 0; i < n; ++i) {
//            step(as, ash[i], ashp[i], nl[i], sp[i], buff, fft, i, _);
            non_linear_coeff(as, i, fft.forward_data());
            fft.execute_forward();
            utils::multiply_by(fft.backward_data(), fft.backward_data_end(), b_);
            if (_ == 0) std::memcpy(nl[i].data(), fft.backward_data(), xs_.size() * sizeof(complex));
            calculate(ash[i], fft.backward_data(), nl[i], buff, i, buff1.data());
            std::memcpy(nl[i].data(), fft.backward_data(), xs_.size() * sizeof(complex));
            for (size_t j = 0; j < xs_.size(); ++j)
                spn[j] = (buff1[j] - ash[i][j]) / dz_;
            for (size_t j = 0; j < xs_.size(); ++j)
                buff[j] = spn[j] - expA_[j] * sp[i][j] + cA_[j] * (
                        buff1[j] * (1. + cA_[j] * dz_ / 2.) + ash[i][j] * expA_[j] * (cA_[j] * dz_ / 2. - 1.));
            std::swap(sp[i], spn);
            std::memcpy(ash[i].data(), buff1.data(), xs_.size() * sizeof(complex));
            std::memcpy(fft.backward_data(), buff1.data(), xs_.size() * sizeof(complex));
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

void isope::model::non_linear_coeff(const cvector2d_t& values, const size_t n, complex* data) const {
    for (size_t i = 0; i <= n; ++i)
        for (size_t j = 0; j <= n - i; ++j)
            for (size_t k = 0; k < values[0].size(); ++k)
                data[k] = values[i][k] * values[j][k] * std::conj(values[n - i - j][k]) * b_;
}

void isope::model::calculate(const cvector1d_t& values, const complex* nl, const cvector1d_t& nlp,
                             const cvector1d_t& buff, const size_t n, complex* data) const {
    if (n == 0)
        for (size_t i = 0; i < values.size(); ++i)
            data[i] = values[i] * expA_[i] + nlfacA_[i] * nl[i] + nlfacAp_[i] * nlp[i];// + a_ * buff[i];
    else
        for (size_t i = 0; i < values.size(); ++i)
            data[i] = values[i] * expA_[i] + nlfacA_[i] * nl[i] + nlfacAp_[i] * nlp[i] + a_ * buff[i];
//            data[i] = values[n][i] * expA_[i] + dz_ / 2. * (nl[i] + expA_[i] * nlp[n][i]) + a_ * buff[i];
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
