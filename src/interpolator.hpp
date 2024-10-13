#pragma once
#include <array>
#include <iostream>

template<typename xType, typename yType, unsigned int N>
class NewtonInterpolator {
private:
    std::array<xType, N> points_;
    std::array<yType, N> dividedDifferences_;
public:
    NewtonInterpolator(const std::array<xType, N> &points, const std::array<yType, N>& values) noexcept : points_{points}, dividedDifferences_(values){
        for (std::size_t i = 0; i < N - 1; i++) {
            for (std::size_t j = N - 1; j > i; j--)
                dividedDifferences_[j] = (dividedDifferences_[j] - dividedDifferences_[j - 1]) / (points[j] - points[j - i - 1]); 
        }
    }
    yType interpolate(const xType& x) const noexcept{
        yType ans = dividedDifferences_[N-1];
        for (std::size_t i = N - 1; i > 0; i--){
            ans = dividedDifferences_[i - 1] + ans * (x - points_[i - 1]);
        }
        return ans;
    }
};
