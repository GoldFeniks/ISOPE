#pragma once

#include <vector>
#include <stdexcept>
#include <tuple>

namespace isope {

    namespace interpolation {

        template<typename CoordT, typename DataT = CoordT>
        class interpolation {

        public:

            template<typename CoordIt = CoordT*, typename DataIt = DataT*>
            interpolation(CoordIt c_begin, CoordIt c_end, DataIt d_begin, DataIt d_end)
                    : coords_(c_begin, c_end), data_(d_begin, d_end)
            {
                if (coords_.size() != data_.size())
                    throw std::runtime_error("Data count and coordinates count must be the same");
            }

            virtual DataT operator()(const CoordT&) const = 0;

        protected:

            const std::vector<CoordT> coords_;
            const std::vector<DataT> data_;

        };


        template<typename CoordT, typename DataT = CoordT>
        class linear : public interpolation<CoordT, DataT> {

        public:

            template<typename CoordIt = CoordT*, typename DataIt = DataT*>
            linear(CoordIt c_begin, CoordIt c_end, DataIt d_begin, DataIt d_end) :
                    interpolation<CoordT, DataT>(c_begin, c_end, d_begin, d_end) {}

            DataT operator()(const CoordT& x) const override {
                size_t a, b;
                std::tie(a, b) = find_interval(x);
                return this->data_[a] + (this->data_[b] - this->data_[a]) *
                                        (x - this->coords_[a]) / (this->coords_[b] - this->coords_[a]);
            }

        private:

            auto find_interval(const CoordT& x) const {
                auto result = std::make_pair(0ul, this->coords_.size() - 1);
                while (result.second - result.first > 1) {
                    const auto m = (result.second + result.first) / 2;
                    if (this->coords_[m] >= x)
                        result.second = m;
                    else
                        result.first = m;
                }
                return result;
            }

        };

        template<typename CoordT, typename DataT = CoordT>
        class polynomial : public interpolation<CoordT, DataT> {

        public:

            template<typename CoordIt = CoordT*, typename DataIt = DataT*>
            polynomial(CoordIt c_begin, CoordIt c_end, DataIt d_begin, DataIt d_end) :
                interpolation<CoordT, DataT>(c_begin, c_end, d_begin, d_end)
            {
                coeffs_.resize(this->coords_.size());
                for (size_t i = 0; i < coeffs_.size(); ++i) {
                    coeffs_[i].resize(coeffs_.size());
                    for (size_t j = 0; j < coeffs_.size(); ++j)
                        coeffs_[i][j] = i == j ? 0 : this->coords_[i] - this->coords_[j];
                }
            }

            DataT operator()(const CoordT& x) const override {
                DataT result = 0;
                for (size_t i = 0; i < coeffs_.size(); ++i) {
                    auto temp = this->data_[i];
                    for (size_t j = 0; j < coeffs_.size(); ++j)
                        temp *= i == j ? 1 : (x - this->coords_[j]) / coeffs_[i][j];
                    result += temp;
                }
                return result;
            }

        private:

            std::vector<std::vector<CoordT>> coeffs_;

        };

    }

}