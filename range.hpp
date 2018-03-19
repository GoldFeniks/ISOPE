#pragma once
#include <cstdint>
#include <cstddef>

namespace range {

    struct bound {

        int64_t value;

        explicit bound(const int64_t& value) : value(value) {}
        bound& operator-();

    };

    struct range {

        const bound min, max;

        range(const bound& min, const bound& max) : min(min), max(max) {}

        template<typename T>
        explicit operator T() const {
            T result(max.value - min.value + 1);
            auto value = min.value;
            for (size_t i = 0; i < result.size(); ++i, ++value)
                result[i] = value;
            return result;
        }

    };

    template<typename T>
    struct scaled_range {

        T scale = 1;
        const bound min, max;

        scaled_range(const bound& min, const bound& max) : min(min), max(max) {}

        template<typename C>
        explicit operator C() const {
            C result(max.value - min.value + 1);
            T value = min.value * scale;
            for (size_t i = 0; i < result.size(); ++i, value += scale)
                result[i] = value;
            return result;
        }

        scaled_range<T>& operator*(const T& value) {
            scale *= value;
            return *this;
        }

        scaled_range<T>& operator/(const T& value) {
            scale /= value;
            return *this;
        }

    };

    template<typename T>
    scaled_range<T>& operator*(const T& value, scaled_range<T>&& range) {
        range.scale *= value;
        return range;
    }

    template<typename T>
    scaled_range<T>& operator/(const T& value, scaled_range<T>&& range) {
        range.scale /= value;
        return range;
    }

    template<typename T>
    scaled_range<T>& operator*(scaled_range<T>&& range, const T& value) {
        range.scale *= value;
        return range;
    }

    template<typename T>
    scaled_range<T>& operator/(scaled_range<T>&& range, const T& value) {
        range.scale /= value;
        return range;
    }

    template<typename T>
    scaled_range<T> operator*(const T& value, const range& range) {
        return scaled_range<T>(range.min, range.max) * value;
    }

    template<typename T>
    scaled_range<T> operator*(const range& range, const T& value) {
        return scaled_range<T>(range.min, range.max) * value;
    }

    template<typename T>
    scaled_range<T> operator/(const T& value, const range& range) {
        return scaled_range<T>(range.min, range.max) * value;
    }

    template<typename T>
    scaled_range<T> operator/(const range& range, const T& value) {
        return scaled_range<T>(range.min, range.max) * value;
    }

    inline range operator,(const bound& min, const bound& max) { return range(min, max); }

    namespace range_literals {

        inline bound operator "" b(const size_t value) { return bound(value); }
        inline range operator "" r(const size_t value) { return range(bound(0), bound(value)); }

    }

}
