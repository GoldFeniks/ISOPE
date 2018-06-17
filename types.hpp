#pragma once

#include <vector>

namespace isope {

    namespace types {

        template<typename T>
        using vector1d_t = std::vector<T>;

        template<typename T>
        using vector2d_t = std::vector<std::vector<T>>;

        using real_t = fft::basic_fft::real_t;

        using complex_t = fft::basic_fft::complex_t;

    }

}
