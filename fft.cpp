#include <utility>
#include "fft.hpp"

isope::fft::fftw::fftw(const int size, const bool same_data) : size_(size),
                                                               dsize_(size),
                                                               same_data_(same_data)
{
    forward_data_ = static_cast<complex*>(fftw_malloc(sizeof(complex) * size));
    forward_data_end_ = forward_data_ + size;
    backward_data_ = same_data ? forward_data_ : static_cast<complex*>(fftw_malloc(sizeof(complex) * size));
    backward_data_end_ = backward_data_ + size;
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
    fftw_free(forward_data_);
    if (!same_data_)
        fftw_free(backward_data_);
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
    utils::divide_by(backward_data_, backward_data_end_, std::sqrt(dsize_));
    return *this;
}

const isope::fft::basic_fft & isope::fft::fftw::normalize_backward() const {
    utils::divide_by(forward_data_, forward_data_end_, dsize_);
    return *this;
}

isope::fft::basic_fft::complex* isope::fft::fftw::forward_data() const { return forward_data_; }

isope::fft::basic_fft::complex* isope::fft::fftw::backward_data() const { return backward_data_; }

isope::fft::basic_fft::complex* isope::fft::fftw::forward_data_end() const { return forward_data_end_; }

isope::fft::basic_fft::complex* isope::fft::fftw::backward_data_end() const { return backward_data_end_; }

int isope::fft::fftw::size() const { return size_; }

isope::fft::fftw::fftw(isope::fft::fftw&& other) noexcept {
    *this = std::move(other);
}

isope::fft::fftw& isope::fft::fftw::operator=(isope::fft::fftw &&other) noexcept {
    size_ = other.size_;
    dsize_ = other.dsize_;
    same_data_ = other.same_data_;
    std::swap(forward_data_, other.forward_data_);
    std::swap(backward_data_, other.backward_data_);
    plan_forward_ = other.plan_forward_;
    plan_backward_ = other.plan_backward_;
    return *this;
}
