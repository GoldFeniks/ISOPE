#pragma once

#include <vector>
#include <array>

namespace isope {

    namespace types {

        template<typename T>
        using vector1d_t = std::vector<T>;

        template<typename T>
        using vector2d_t = std::vector<std::vector<T>>;

        using real_t = fft::basic_fft::real_t;

        using complex_t = fft::basic_fft::complex_t;

        template<typename T, size_t N>
        using array1d_t = std::array<T, N>;

        template<typename T, size_t N>
        using array2d_t = std::array<std::array<T, N>, N>;

    }

}
