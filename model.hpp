#pragma once

#include <vector>
#include <functional>
#include <cstring>
#include "fft.hpp"
#include "types.hpp"

namespace isope {

    class model {

    public:

        using complex = fft::basic_fft::complex;
        using rvector1d_t = types::vector1d_t<double>;
        using cvector1d_t = types::vector1d_t<complex>;
        using cvector2d_t = types::vector2d_t<complex>;

        model(double k0, double sigma, double eps, double v,
              size_t n, double x0, double x1, size_t nz, double z0, double z1);

        cvector2d_t solve(std::function<complex(double, double)> init_cond, size_t n, size_t m) const;

    private:

        const complex a_, b_;
        const double k0_, sigma_, eps_, v_, dx_, dz_;
        const rvector1d_t xs_, zs_, k_;
        cvector1d_t expA_, nlfacA_, nlfacAp_, cA_;
        const size_t x_size_, z_size_, bytes_size_;

        struct flags {

            static constexpr size_t none = 0;
            static constexpr size_t first_call = 1 << 0;
            static constexpr size_t first_equation = 1 << 1;
            static constexpr size_t last_equation = 1 << 2;

        };

        static rvector1d_t vector_k_(double l, size_t n);
        void non_linear_coeff(const cvector2d_t& values, size_t n, complex* data) const;
        inline void store(cvector2d_t& result, std::vector<size_t>& indexes, const cvector1d_t& values,
                          size_t index, size_t z_index, bool do_store) const;

        inline void diff(const cvector1d_t& ash0, const cvector1d_t& ash1, cvector1d_t& sp) const;

        template<size_t Flags>
        void calculate(const cvector1d_t& values, const complex* nl, const cvector1d_t& nlp,
                                     const cvector1d_t& buff, const size_t n, complex* data) const {
            for (size_t i = 0; i < values.size(); ++i)
                data[i] = values[i] * expA_[i] + nlfacA_[i] * nl[i] + nlfacAp_[i] * nlp[i] +
                          (utils::has_flag(Flags, flags::first_equation) ? 0. : a_ * buff[i]);
        }

        template<size_t Flags>
        inline void diff0(const cvector1d_t& ash, const cvector1d_t& nl, cvector1d_t& sp) const {
            for (size_t i = 0; i < ash.size(); ++i)
                sp[i] = utils::has_flag(Flags, flags::first_equation)
                        ? cA_[i] * ash[i] + nl[i]
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
                    dd1[i] = spn[i] - expA_[i] * sp[i] + cA_[i] * (
                            ash1[i] * (1. + cA_[i] * dz_ / 2.) + ash0[i] * expA_[i] * (cA_[i] * dz_ / 2. - 1.));
                std::swap(sp, spn);
            }
            std::swap(ash0, ash1);
            std::memcpy(ash1.data(), buff.data(), bytes_size_);
            std::memcpy(fft.backward_data(), buff.data(), bytes_size_);
            fft.execute_backward().normalize_backward();
            std::memcpy(as[n].data(), fft.forward_data(), bytes_size_);
        }

        complex a0sol(double x, double z) const;
        complex a1sol(double x, double z) const;
        complex der(double x, double z) const;

    };

}

