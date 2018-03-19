#pragma once

#include <vector>
#include <functional>
#include <memory>
#include "fft.hpp"
#include "types.hpp"

namespace isope {

    class model {

    public:

        using complex = fft::basic_fft::complex;
        using rvector1d_t = types::vector1d_t<double>;
        using cvector1d_t = types::vector1d_t<complex>;
        using cvector2d_t = types::vector2d_t<complex>;

        model(complex k0, double sigma, double eps, double v,
              size_t n, double x0, double x1, size_t nz, double z0, double z1);

        cvector2d_t solve(std::function<typename complex(double)> init_cond) const;

    private:

        const complex k0_;
        const double sigma_, eps_, v_;
        const std::shared_ptr<fft::basic_fft> fft_;
        const rvector1d_t xs_, zs_;

    };

}

