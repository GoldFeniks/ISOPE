#pragma once

#include <complex>
#include "utils.hpp"
#include "fftw3.h"

namespace isope {

    namespace fft {

        class basic_fft {

        public:

            typedef double real_t;
            typedef std::complex<real_t> complex_t;

            virtual const basic_fft& execute_forward() const = 0;
            virtual const basic_fft& execute_backward() const = 0;
            virtual const basic_fft& normalize_forward() const = 0;
            virtual const basic_fft& normalize_backward() const = 0;
            virtual complex_t* forward_data() = 0;
            virtual const complex_t* forward_data() const = 0;
            virtual complex_t& forward_data(const size_t&) = 0;
            virtual const complex_t& forward_data(const size_t&) const = 0;
            virtual complex_t* backward_data() = 0;
            virtual const complex_t* backward_data() const = 0;
            virtual complex_t& backward_data(const size_t&) = 0;
            virtual const complex_t& backward_data(const size_t&) const = 0;
            virtual int size() const = 0;

        };


        class fftw : public basic_fft {

        public:

            fftw() = delete;
            fftw(const fftw&) = delete;
            fftw(fftw&& other) noexcept;

            explicit fftw(int size, bool same_data = false);
            ~fftw();

            fftw& operator=(fftw&& other) noexcept;

            const basic_fft& execute_forward() const override;
            const basic_fft& execute_backward() const override;
            const basic_fft& normalize_forward() const override;
            const basic_fft& normalize_backward() const override;
            complex_t* forward_data() override;
            const complex_t* forward_data() const override;
            complex_t& forward_data(const size_t& index) override;
            const complex_t& forward_data(const size_t& index) const override;
            complex_t* backward_data() override;
            const complex_t* backward_data() const override;
            complex_t& backward_data(const size_t& index) override;
            const complex_t& backward_data(const size_t& index) const override;
            const complex_t* forward_data_end() const;
            const complex_t* backward_data_end() const;
            int size() const override;

        private:

            int size_ = 0;
            real_t dsize_ = 0;
            bool same_data_ = false;
            complex_t* forward_data_ = nullptr, *backward_data_ = nullptr,
                *forward_data_end_ = nullptr, *backward_data_end_ = nullptr;
            fftw_plan plan_forward_ = fftw_plan(), plan_backward_ = fftw_plan();

        };

    }

}
