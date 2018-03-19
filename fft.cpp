#include "fft.hpp"

isope::fft::fftw::fftw(const int size, const bool same_data) : size_(size),
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

isope::fft::fftw::~fftw() {
    fftw_destroy_plan(plan_forward_);
    fftw_destroy_plan(plan_backward_);
    delete[] forward_data_;
    if (!same_data_)
        delete[] backward_data_;
}

const isope::fft::basic_fft & isope::fft::fftw::execute_forward() const {
    fftw_execute(plan_forward_);
    return *this;
}

const isope::fft::basic_fft & isope::fft::fftw::execute_backward() const {
    fftw_execute(plan_backward_);
    return *this;
}

const isope::fft::basic_fft & isope::fft::fftw::normalize_forward() const {
    multiply_by(backward_data_, std::sqrt(dsize_));
    return *this;
}

const isope::fft::basic_fft & isope::fft::fftw::normalize_backward() const {
    multiply_by(forward_data_, dsize_);
    return *this;
}

isope::fft::basic_fft::complex *isope::fft::fftw::forward_data() const { return forward_data_; }

isope::fft::basic_fft::complex *isope::fft::fftw::backward_data() const { return backward_data_; }

int isope::fft::fftw::size() const { return size_; }

void isope::fft::fftw::multiply_by(complex* data, const double value) const {
    for (auto i = 0; i < size(); ++i)
        *(data++) /= value;
}
