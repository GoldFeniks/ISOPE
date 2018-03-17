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

            explicit fftw(const int size, const bool same_data = false) : size_(size),
                                                                         dsize_(size),
                                                                         same_data_(same_data),
                                                                         forward_data_(new complex[size])
            {
                backward_data_ = same_data ? forward_data_ : new complex[size];
                plan_forward_ = fftw_plan_dft_1d(size, reinterpret_cast<fftw_complex*>(forward_data_),
                                                 reinterpret_cast<fftw_complex*>(backward_data_),
                                                 FFTW_FORWARD, FFTW_MEASURE);
                plan_backward_ = fftw_plan_dft_1d(size, reinterpret_cast<fftw_complex*>(backward_data_),
                                                 reinterpret_cast<fftw_complex*>(forward_data_),
                                                 FFTW_BACKWARD, FFTW_MEASURE);
            }

            ~fftw() {
                fftw_destroy_plan(plan_forward_);
                fftw_destroy_plan(plan_backward_);
                delete[] forward_data_;
                if (!same_data_)
                    delete[] backward_data_;
            }

            void execute_forward(const normalization norm) const override {
                fftw_execute(plan_forward_);
                normalize(backward_data_, norm);
            }

            void execute_backward(const normalization norm) const override {
                fftw_execute(plan_backward_);
                normalize(forward_data_, norm);
            }

            complex* forward_data() const override { return forward_data_; }
            complex* backward_data() const override { return backward_data_; }
            int size() const override { return size_; }

        private:

            int size_ = 0;
            double dsize_ = 0;
            bool same_data_ = false;
            complex* forward_data_ = nullptr, *backward_data_ = nullptr;
            fftw_plan plan_forward_ = fftw_plan(), plan_backward_ = fftw_plan();

            void normalize(complex* data, const normalization norm) const {
                if (norm == normalization::none) return;
                const double c = norm == normalization::twice ? dsize_ : std::sqrt(dsize_);
                for (auto i = 0; i < size(); ++i)
                    *(data++) /= c;
            }

        };

    }

}