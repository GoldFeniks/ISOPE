#include <iostream>
#include "model.hpp"
#include "utils.hpp"

isope::model::model(const complex k0, const double sigma, const double eps, const double v,
                    const size_t nx, const double x0, const double x1,
                    const size_t nz, const double z0, const double z1) :
    k0_(k0), sigma_(sigma), eps_(eps), v_(v), fft_(std::make_shared<fft::fftw>(nx, true)),
    xs_(utils::mesh<rvector1d_t>(x0, x1, nx)), zs_(utils::mesh<rvector1d_t>(z0, z1, nz))
{
    for (auto it : xs_)
        std::cout << it << std::endl;
}

isope::model::cvector2d_t isope::model::solve(std::function<typename complex(double)> init_cond) const {
    cvector2d_t result(zs_.size(), cvector1d_t(xs_.size()));
    for (size_t i = 0; i < xs_.size(); ++i)
        result[0][i] = init_cond(xs_[i]);
    for (size_t i = 1; i < zs_.size(); ++i) {

    }
}
