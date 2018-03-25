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

        cvector2d_t solve(std::function<complex(double)> init_cond, size_t n, size_t m) const;

    private:

        const complex a_, b_;
        const double k0_, sigma_, eps_, v_;
        const rvector1d_t xs_, zs_, k_;
        cvector1d_t expA_, nlfacA_, nlfacAp_;

        const double dt_ = 0.5e-8;

        static rvector1d_t vector_k_(double l, size_t n);

        void non_linear_coeff(const cvector2d_t& values, size_t n, complex* data) const;

    };

}

