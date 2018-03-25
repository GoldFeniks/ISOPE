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

        template<typename It1, typename It2, typename T>
        void multiply_by(It1 from, It2 to, const T& value) {
            while (from != to) {
                *from *= value;
                ++from;
            }
        };

        template<typename It1, typename It2, typename T>
        void divide_by(It1 from, It2 to, const T& value) {
            while (from != to) {
                *from /= value;
                ++from;
            }
        };


    }

}