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

    template<typename T>
    class dormand_prince_coefficients1 : public coefficients<T, 7> {

    public:

        static constexpr types::array2d_t<T, 7> a_table =
                {
                        T(0),                T(0),               T(0),                T(0),             T(0),               T(0),          T(0),
                        T(1)     / T(5),     T(0),               T(0),                T(0),             T(0),               T(0),          T(0),
                        T(3)     / T(40),    T(3)     / T(40),   T(0),                T(0),             T(0),               T(0),          T(0),
                        T(44)    / T(45),   -T(56)    / T(15),   T(32)    / T(9),     T(0),             T(0),               T(0),          T(0),
                        T(19372) / T(6561), -T(25360) / T(2187), T(64448) / T(6561), -T(212) / T(729),  T(0),               T(0),          T(0),
                        T(9017)  / T(3168), -T(355)   / T(33),   T(46732) / T(5247),  T(49)  / T(176), -T(5103) / T(18656), T(0),          T(0),
                        T(35)    / T(384),   T(0),               T(500)   / T(1113),  T(125) / T(192), -T(2187) / T(6784),  T(11) / T(84), T(0),
                };
        static constexpr types::array1d_t<T, 7> b_table =
                {
                         T(35)   / T(384),
                         T(0),
                         T(500)  / T(1113),
                         T(125)  / T(192),
                        -T(2187) / T(6784),
                         T(11)   / T(84),
                         T(0)
                };
        static constexpr types::array1d_t<T, 7> c_table =
                {
                        T(0),
                        T(1) / T(5),
                        T(3) / T(10),
                        T(4) / T(5),
                        T(8) / T(9),
                        T(1),
                        T(1)
                };

        dormand_prince_coefficients1() : coefficients<T, 7>(a_table, b_table, c_table) {}
    };

    template<typename T>
    class dormand_prince_coefficients2 : public coefficients<T, 7> {

    public:

        static constexpr types::array2d_t<T, 7> a_table = dormand_prince_coefficients1<T>::a_table;
        static constexpr types::array1d_t<T, 7> b_table =
                {
                         T(5179)  / T(57600),
                         T(0),
                         T(7571)  / T(16695),
                         T(393)   / T(640),
                        -T(92097) / T(339200),
                         T(187)   / T(2100),
                         T(1)     / T(40)
                };
        static constexpr types::array1d_t<T, 7> c_table = dormand_prince_coefficients1<T>::c_table;

        dormand_prince_coefficients2() : coefficients<T, 7>(a_table, b_table, c_table) {}

    };

    template<size_t Order, typename ValueType, typename ArgumentType, typename Coeffs = runge_coefficients<ValueType, Order>>
    class runge_kutta {

    public:

        static auto create() {
            return std::bind(solve, std::placeholders::_1, std::placeholders::_2,
                             std::placeholders::_3, std::placeholders::_4, std::placeholders::_5, Coeffs());
        }

        static types::vector1d_t<ValueType> solve(const ArgumentType& a, const ArgumentType& b, const ValueType& y0,
                                                  const size_t n, const std::function<ValueType(ArgumentType, ValueType)>& f,
                                                  const Coeffs& coeffs) {
            const auto d = (b - a) / (n - 1);
            types::vector1d_t<ValueType> result(n, ValueType(0));
            result[0] = y0;
            types::array1d_t<ValueType, Order> k;
            auto x = a;
            for (size_t i = 1; i < n; ++i) {
                for (size_t j = 0; j < Order; ++j) {
                    auto buff = ValueType(0);
                    for (size_t m = 0; m < j; ++m)
                        buff += coeffs.get_a(j, m) * k[m];
                    k[j] = f(x + d * coeffs.get_c(j), result[i - 1] + buff * d);
                    result[i] += coeffs.get_b(j) * k[j];
                }
                result[i] *= d;
                result[i] += result[i - 1];
                x += d;
            }
            return result;
        }

    };

    template<size_t Order, typename ValueType, typename ArgumentType, typename Coeffs = runge_coefficients<ValueType, Order>>
    class lazy_runge_kutta {

    public:

        lazy_runge_kutta(const ArgumentType& a, const ArgumentType& b, const ValueType& y0,
                         const size_t n, const std::function<ValueType(ArgumentType, ValueType, size_t)>& f)
                         : _d((b - a) / (n - 1)), _lastValue(y0), _lastArgument(a), _f(f), _max(n) {}

        ValueType operator()(types::array1d_t<ValueType, Order>& k) {
            if (!*this) return ValueType(0);
            ++_current;
            auto result = ValueType(0);
            for (size_t i = 0; i < Order; ++i) {
                auto buff = ValueType(0);
                for (size_t j = 0; j < i; ++j)
                    buff += _coeffs.get_a(i, j) * k[j];
                result += _coeffs.get_b(i) * (k[i] = _f(_lastArgument + _d * _coeffs.get_c(i), _lastValue + buff * _d, i));
            }
            _lastArgument += _d;
            return _lastValue += result * _d;
        }

        ValueType operator()() {
            static types::array1d_t<ValueType, Order> k;
            return (*this)(k);
        }

        operator bool() const {
            return _max > _current;
        }

    private:

        const ArgumentType _d;
        ArgumentType _lastArgument;
        ValueType _lastValue;
        const std::function<ValueType(ArgumentType, ValueType, size_t)> _f;
        const Coeffs _coeffs;
        const size_t _max;
        size_t _current = 1;

    };

    template<size_t Order, typename ValueType, typename ArgumentType, typename Coeffs = runge_coefficients<ValueType, Order>>
    class lazy_2d_runge_kutta {

    public:

        lazy_2d_runge_kutta(const ArgumentType& a, const ArgumentType& b,
                const types::vector1d_t<ArgumentType>& x0, const types::vector1d_t<ValueType>& y0,
                const size_t n, const std::function<ValueType(ArgumentType, ArgumentType, ValueType, size_t)>& f)
                : _d((b - a) / (n - 1)), _lastValue(y0), _args(x0), _lastArgument(a), _f(f), _max(n) {}

        types::vector1d_t<ValueType> operator()(types::vector1d_t<types::array1d_t<ValueType, Order>>& k) {
            if (!*this) return types::vector1d_t<ValueType>(_lastValue.size(), 0);
            ++_current;
            for (size_t m = 0; m < _lastValue.size(); ++m) {
                auto lv = _lastValue[m];
                for (size_t i = 0; i < Order; ++i) {
                    auto buff = ValueType(0);
                    for (size_t j = 0; j < i; ++j)
                        buff += _coeffs.get_a(i, j) * k[m][j];
                    _lastValue[m] += _coeffs.get_b(i) * (k[m][i] = _f(_args[m], _lastArgument + _d * _coeffs.get_c(i), lv + buff * _d, i)) * _d;
                }
            }
            _lastArgument += _d;
            return _lastValue;
        }

        types::vector1d_t<ValueType> operator()() {
            static types::vector1d_t<types::array1d_t<ValueType, Order>> k(_lastValue.size());
            return (*this)(k);
        }

        operator bool() const {
            return _max > _current;
        }

    private:

        const ArgumentType _d;
        ArgumentType _lastArgument;
        types::vector1d_t<ValueType> _lastValue, _args;
        const std::function<ValueType(ArgumentType, ArgumentType, ValueType, size_t)> _f;
        const Coeffs _coeffs;
        const size_t _max;
        size_t _current = 1;
    };

}
