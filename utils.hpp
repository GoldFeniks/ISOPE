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

        template<typename T>
        T expanded_mesh(const typename T::value_type& a, const typename T::value_type& b, const size_t n) {
            const auto d = (b - a) / 2;
            return mesh<T>(a - d, b + d, 2 * n);
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

        template<typename T>
        constexpr bool has_flag(const T& a, const T& b) {
            return a & b;
        }

        template<typename T, typename F>
        void init_vector(T& vector, size_t size, const F& init_fun) {
            vector.resize(size);
            for (size_t i = 0; i < size; ++i)
                vector[i] = init_fun(i);
        };

    }

}