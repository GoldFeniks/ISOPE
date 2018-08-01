#pragma once

#include <vector>
#include <functional>
#include <tuple>
#include "types.hpp"

namespace isope {

    template<typename T, size_t N>
    class coefficients {

    public:

        static constexpr auto order = N;

        coefficients() = default;

        coefficients(const types::array2d_t<T, N>& a_table, const types::array1d_t<T, N>& b_table, const types::array1d_t<T, N>& c_table)
                : _a_table(a_table), _b_table(b_table), _c_table(c_table) {}

        const T& get_a(const size_t i, const size_t j) const {
            return _a_table[i][j];
        };

        const std::vector<T>& get_a(const size_t i) const {
            return _a_table[i];
        }

        const T& get_b(const size_t i) const {
            return _b_table[i];
        }

        const T& get_c(const size_t i) const  {
            return _c_table[i];
        }

        std::tuple<std::vector<T>, T ,T> get(size_t i) const {
            return { get_a(i), get_b(i), get_c(i) };
        };

    private:

        const types::array2d_t<T, N> _a_table;
        const types::array1d_t<T, N> _b_table, _c_table;

    };

    template<typename T, size_t N>
    class runge_coefficients : public coefficients<T, N> {};

    template<typename T>
    class runge_coefficients<T, 1> : public coefficients<T, 1> {

    public:

        static constexpr types::array2d_t<T, 1> a_table = { T(1) };
        static constexpr types::array1d_t<T, 1> b_table = { T(1) };
        static constexpr types::array1d_t<T, 1> c_table = { T(0) };

        runge_coefficients() : coefficients<T, 1>(a_table, b_table, c_table) {};

    };

    template<typename T>
    class runge_coefficients<T, 2> : public coefficients<T, 2> {

    public:

        static constexpr types::array2d_t<T, 2> a_table = { T(0), T(0),  T(1) / T(2), T(0)  };
        static constexpr types::array1d_t<T, 2> b_table = { T(0), T(1) };
        static constexpr types::array1d_t<T, 2> c_table = { T(0), T(1) / T(2) };

        runge_coefficients() : coefficients<T, 2>(a_table, b_table, c_table) {};

    };

    template<typename T>
    class runge_coefficients<T, 4> : public coefficients<T, 4> {

    public:

        static constexpr types::array2d_t<T, 4> a_table = {
                                                            T(0),        T(0),        T(0), T(0),
                                                            T(1) / T(2), T(0),        T(0), T(0),
                                                            T(0),        T(1) / T(2), T(0), T(0),
                                                            T(0),        T(0),        T(1), T(0)
                                                          };
        static constexpr types::array1d_t<T, 4> b_table = {
                                                            T(1) / T(6),
                                                            T(1) / T(3),
                                                            T(1) / T(3),
                                                            T(1) / T(6)
                                                          };
        static constexpr types::array1d_t<T, 4> c_table = {
                                                            T(0),
                                                            T(1) / T(2),
                                                            T(1) / T(2),
                                                            T(1)
                                                          };

        runge_coefficients() : coefficients<T, 4>(a_table, b_table, c_table) {}
    };

    template<typename ValueType, typename ArgumentType, typename Coeffs = runge_coefficients<ValueType, 4>>
    class runge_kutta {};

}
