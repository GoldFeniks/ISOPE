#pragma once

#include <cstddef>

namespace isope {

    namespace utils {

        template<typename T>
        T mesh(const typename T::value_type& a, const typename T::value_type& b, const size_t n) {
            T result(n);
            typename T::value_type d = (b - a) / (n - 1);
            for (size_t i = 0; i < n; ++i)
                result[i] = a + i * d;
            return result;
        }

    }

}