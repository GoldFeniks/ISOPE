#pragma once

#include <vector>
#include <functional>
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

        static rvector1d_t vector_k_(double l, size_t n);

        static void non_linear_coeff(const cvector2d_t& values, size_t n, complex* data);
        void calculate(const cvector2d_t& values, const complex* nl, const cvector2d_t& nlp,
                       const cvector1d_t& buff, size_t n, complex* data) const;

        complex a0sol(double x, double z) const;
        complex a1sol(double x, double z) const;
        complex der(double x, double z) const;

    };

}

