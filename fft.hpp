#pragma once

#include <complex>
#include "utils.hpp"
#include "fftw3.h"

namespace isope {

    namespace fft {

        class basic_fft {

        public:

            typedef std::complex<double> complex;

            virtual const basic_fft& execute_forward() const = 0;
            virtual const basic_fft& execute_backward() const = 0;
            virtual const basic_fft& normalize_forward() const = 0;
            virtual const basic_fft& normalize_backward() const = 0;
            virtual complex* forward_data() = 0;
            virtual const complex* forward_data() const = 0;
            virtual complex& forward_data(const size_t&) = 0;
            virtual const complex& forward_data(const size_t&) const = 0;
            virtual complex* backward_data() = 0;
            virtual const complex* backward_data() const = 0;
            virtual complex& backward_data(const size_t&) = 0;
            virtual const complex& backward_data(const size_t&) const = 0;
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
            complex* forward_data() override;
            const complex* forward_data() const override;
            complex& forward_data(const size_t& index) override;
            const complex& forward_data(const size_t& index) const override;
            complex* backward_data() override;
            const complex* backward_data() const override;
            complex& backward_data(const size_t& index) override;
            const complex& backward_data(const size_t& index) const override;
            const complex* forward_data_end() const;
            const complex* backward_data_end() const;
            int size() const override;

        private:

            int size_ = 0;
            double dsize_ = 0;
            bool same_data_ = false;
            complex* forward_data_ = nullptr, *backward_data_ = nullptr,
                *forward_data_end_ = nullptr, *backward_data_end_ = nullptr;
            fftw_plan plan_forward_ = fftw_plan(), plan_backward_ = fftw_plan();

        };

    }

}