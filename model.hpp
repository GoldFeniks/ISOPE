#pragma once

#include <vector>
#include <functional>
#include <cstring>
#include <iostream>
#include "fft.hpp"
#include "types.hpp"

namespace isope {

    using namespace std::complex_literals;

    template<bool UseCAP>
    class model_ {

    public:

        using complex = types::complex;
        using rvector1d_t = types::vector1d_t<double>;
        using cvector1d_t = types::vector1d_t<complex>;
        using cvector2d_t = types::vector2d_t<complex>;

        model_(const double k0, const double eps,
              const size_t nx, const double x0, const double x1,
              const size_t nz, const double z0, const double z1, const double delta) :
        a_(1i / (2. * k0)), b_(1i * eps * k0 / 2.), k0_(k0),
                dx_((x1 - x0) / (nx - 1)), dz_((z1 - z0) / (nz - 1)),
        xs_(UseCAP ? utils::expanded_mesh<rvector1d_t>(x0, x1, nx) : utils::mesh<rvector1d_t>(x0, x1, nx)),
        zs_(utils::mesh<rvector1d_t>(z0, z1, nz)),
        k_(vector_k_(xs_.back() - xs_[0], xs_.size())),
        x_size_(xs_.size()), z_size_(zs_.size()), bytes_size_(xs_.size() * sizeof(complex))
        {
            if (UseCAP)
                utils::init_vector(s_, x_size_,
                                   [this, &x0, &x1, &delta](const size_t& i) {
                                       if (i < xs_.size() / 4)
                                           return std::pow((xs_[i] - x0) / delta, 2) / 2. / k0_;
                                       if (i > xs_.size() / 4 * 3)
                                           return std::pow((xs_[i] - x1) / delta, 2) / 2. / k0_;
                                       return 0.;
                });
            utils::init_vector(c_coeff_, x_size_, [this](const size_t& i) { return -k_[i] * a_; });
            utils::init_vector(exp_coeff_, x_size_, [this](const size_t& i) { return std::exp(dz_ * c_coeff_[i]); });
            utils::init_vector(nl_coeff_, x_size_,
                               [this](const size_t& i) {
                                    return i == 0
                                           ? dz_ / 2
                                           : (exp_coeff_[i] * (1. + 1. / c_coeff_[i] / dz_) - 1. / c_coeff_[i] / dz_ - 2.) / c_coeff_[i];
            });
            utils::init_vector(nlp_coeff_, x_size_,
                               [this](const size_t& i) {
                                   return i == 0
                                          ? dz_ / 2
                                          : (exp_coeff_[i] * (-1. / c_coeff_[i] / dz_) + 1. / c_coeff_[i] / dz_ + 1.) / c_coeff_[i];
                               });
            utils::init_vector(integral_coeff1_, x_size_, [this](const size_t& i) { return 1. + c_coeff_[i] * dz_ / 2.; });
            utils::init_vector(integral_coeff2_, x_size_, [this](const size_t& i) { return exp_coeff_[i] * (c_coeff_[i] * dz_ / 2. - 1.); });
        }

        template<typename F>
        cvector2d_t solve(const F& init_cond, const size_t n, const size_t m, const size_t verb) const {
            cvector2d_t result(zs_.size() / m, cvector1d_t(x_size_, 0.)),
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
                result[0][i] = as[0][i] = init_cond(xs_[i]);

            std::memcpy(fft.forward_data(), as[0].data(), bytes_size_);
            fft.execute_forward();
            std::memcpy(ash1[0].data(), fft.backward_data(), bytes_size_);

            step<flags::first_call | flags::first_equation>(0, as, ash0[0], ash1[0], sp[0], dd0[0], dd1[0], nl[0], fft);

            if (n > 0) {

                for (size_t i = 1; i <= n; ++i) {
                    step<flags::first_equation>(0, as, ash0[0], ash1[0], sp[0], dd0[0], dd1[0], nl[0], fft);
                    store(result, indexes, as[0], 0, i, i % m == 0);
                    for (size_t j = 1; j < i; ++j) {
                        step(j, as, ash0[j], ash1[j], sp[j], dd0[j - 1], dd1[j], nl[j], fft);
                        store(result, indexes, as[j], j, i - j, (i - j) % m == 0);
                    }
                    step<flags::first_call>(i, as, ash0[i], ash1[i], sp[i], dd0[i - 1], dd1[i], nl[i], fft);
                    std::swap(dd0, dd1);
                    if (verb && (i + 1) % verb == 0)
                        std::cout << (i + 1) << std::endl;
                }

                for (size_t i = n + 1; i < z_size_ - n; ++i) {
                    step<flags::first_equation>(0, as, ash0[0], ash1[0], sp[0], dd0[0], dd1[0], nl[0], fft);
                    store(result, indexes, as[0], 0, i, i % m == 0);
                    for (size_t j = 1; j < n; ++j) {
                        step(j, as, ash0[j], ash1[j], sp[j], dd0[j - 1], dd1[j], nl[j], fft);
                        store(result, indexes, as[j], j, i - j, (i - j) % m == 0);
                    }
                    step<flags::last_equation>(n, as, ash0[n], ash1[n], sp[n], dd0[n - 1], dd1[n], nl[n], fft);
                    store(result, indexes, as[n], n, i - n, (i - n) % m == 0);
                    std::swap(dd0, dd1);
                    if (verb && (i + 1) % verb == 0)
                        std::cout << (i + 1) << std::endl;
                }

                for (size_t i = 0; i < n; ++i) {
                    step<flags::last_equation>(n, as, ash0[n], ash1[n], sp[n], dd0[n - 1], dd1[n], nl[n], fft);
                    store(result, indexes, as[n], n, i - n, (i - n) % m == 0);
                    for (size_t j = n - 1; j > i; --j) {
                        step(j, as, ash0[j], ash1[j], sp[j], dd0[j - 1], dd1[j], nl[j], fft);
                        store(result, indexes, as[j], j, i - j, (i - j) % m == 0);
                    }
                    std::swap(dd0, dd1);
                    if (verb && (z_size_ - n + i + 1) % verb == 0)
                        std::cout << (z_size_ - n + i + 1) << std::endl;
                }
            }
            else {

                for (size_t i = 1; i < z_size_; ++i) {
                    step<flags::first_equation | flags::last_equation>(0, as, ash0[0], ash1[0], sp[0], dd0[0], dd1[0], nl[0], fft);
                    store(result, indexes, as[0], 0, i, i % m == 0);
                    if (verb && (i + 1) % verb == 0)
                        std::cout << (i + 1) << std::endl;
                }

            }

            return result;
        }

        auto dx() const {
            return dx_;
        }

        auto dz() const {
            return dz_;
        }

        const rvector1d_t& x_coords() const {
            return xs_;
        }

        const rvector1d_t& z_coords() const {
            return zs_;
        }

    private:

        const complex a_, b_;
        const double k0_, dx_, dz_;
        const rvector1d_t xs_, zs_, k_;
        rvector1d_t s_;
        cvector1d_t c_coeff_, exp_coeff_, nl_coeff_, nlp_coeff_, integral_coeff1_, integral_coeff2_;
        const size_t x_size_, z_size_, bytes_size_;

        struct flags {

            static constexpr size_t none = 0u;
            static constexpr size_t first_call = 1u << 0u;
            static constexpr size_t first_equation = 1u << 1u;
            static constexpr size_t last_equation = 1u << 2u;

        };

        static rvector1d_t vector_k_(const double l, const size_t n) {
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

        void non_linear_coeff(const cvector2d_t& values, const size_t n, complex* data) const {
            for (size_t i = 0; i <= n; ++i)
                for (size_t j = 0; j <= n - i; ++j)
                    for (size_t k = 0; k < values[0].size(); ++k)
                        data[k] = values[i][k] * values[j][k] * std::conj(values[n - i - j][k]) * b_;
            if (UseCAP) {
                for (size_t i = 0; i < x_size_ / 4; ++i)
                    data[i] -= values[n][i] * s_[i];
                for (size_t i = x_size_ / 4 * 3; i < x_size_; ++i)
                    data[i] -= values[n][i] * s_[i];
            }
        }

        inline void store(cvector2d_t& result, std::vector<size_t>& indexes, const cvector1d_t& values,
                          const size_t index, const size_t z_index, const bool do_store) const {
            if (!do_store) return;
            auto& ind = indexes[index];
            for (size_t i = 0; i < values.size(); ++i)
                result[ind][i] += values[i] * std::exp(1i * k0_ * zs_[z_index]);
            ++ind;
        }

        inline void diff(const cvector1d_t& ash0, const cvector1d_t& ash1, cvector1d_t& sp) const {
            for (size_t i = 0; i < x_size_; ++i)
                sp[i] = (ash1[i] - ash0[i]) / 2. / dz_;
        }

        template<size_t Flags>
        inline void calculate(const cvector1d_t& values, const complex* nl, const cvector1d_t& nlp,
                                     const cvector1d_t& buff, const size_t n, complex* data) const {
            for (size_t i = 0; i < values.size(); ++i)
                data[i] = values[i] * exp_coeff_[i] + nl_coeff_[i] * nl[i] + nlp_coeff_[i] * nlp[i] +
                          (utils::has_flag(Flags, flags::first_equation) ? 0. : a_ * buff[i]);
        }

        template<size_t Flags>
        inline void diff0(const cvector1d_t& ash, const cvector1d_t& nl, cvector1d_t& sp) const {
            for (size_t i = 0; i < ash.size(); ++i)
                sp[i] = utils::has_flag(Flags, flags::first_equation)
                        ? c_coeff_[i] * ash[i] + nl[i]
                        : ash[i] / dz_;
        }

        template<size_t Flags = flags::none>
        inline void step(const size_t n, cvector2d_t& as, cvector1d_t& ash0, cvector1d_t& ash1,
                         cvector1d_t& sp, cvector1d_t& dd0, cvector1d_t& dd1, cvector1d_t& nl, fft::fftw& fft) const {

            static cvector1d_t buff(x_size_, 0.);
            static cvector1d_t spn(x_size_, 0.);

            non_linear_coeff(as, n, fft.forward_data());
            fft.execute_forward();
            if (utils::has_flag(Flags, flags::first_call)) std::memcpy(nl.data(), fft.backward_data(), bytes_size_);
            calculate<Flags>(ash1, fft.backward_data(), nl, dd0, n, buff.data());
            std::memcpy(nl.data(), fft.backward_data(), bytes_size_);

            if (utils::has_flag(Flags, flags::first_call))
                diff0<Flags>(buff, nl, sp);
            else if (!utils::has_flag(Flags, flags::last_equation)) {
                if (utils::has_flag(Flags, flags::first_equation))
                    diff0<Flags>(ash1, nl, spn);
                else
                    diff(ash0, buff, spn);
                for (size_t i = 0; i < x_size_; ++i)
                    dd1[i] = spn[i] - exp_coeff_[i] * sp[i] + c_coeff_[i] * (ash1[i] * integral_coeff1_[i] + ash0[i] * integral_coeff2_[i]);
                std::swap(sp, spn);
            }
            std::swap(ash0, ash1);
            std::memcpy(ash1.data(), buff.data(), bytes_size_);
            std::memcpy(fft.backward_data(), buff.data(), bytes_size_);
            fft.execute_backward().normalize_backward();
            std::memcpy(as[n].data(), fft.forward_data(), bytes_size_);
        }

    };

    using model = model_<false>;
    using cap_model = model_<true>;

}

