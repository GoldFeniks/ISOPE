#pragma once

#include <complex>
#include "fftw3.h"

namespace isope {

    namespace fft {

        enum class normalization { none, once, twice };

        class basic_fft {

        public:

            typedef std::complex<double> complex;

            virtual void execute_forward(normalization) const = 0;
            virtual void execute_backward(normalization) const = 0;
            virtual complex* forward_data() const = 0;
            virtual complex* backward_data() const = 0;
            virtual int size() const = 0;

        };


        class fftw : public basic_fft {

        public:

            fftw() = delete;
            fftw(const fftw&) = delete;
            fftw(fftw&&) = delete;

            explicit fftw(int size, bool same_data = false);
            ~fftw();

            void execute_forward(normalization norm) const override;
            void execute_backward(normalization norm) const override;
            complex* forward_data() const override;
            complex* backward_data() const override;
            int size() const override;

        private:

            int size_ = 0;
            double dsize_ = 0;
            bool same_data_ = false;
            complex* forward_data_ = nullptr, *backward_data_ = nullptr;
            fftw_plan plan_forward_ = fftw_plan(), plan_backward_ = fftw_plan();

            void normalize(complex* data, normalization norm) const;

        };

    }

}