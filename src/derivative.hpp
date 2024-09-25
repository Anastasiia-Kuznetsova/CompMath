#pragma once
#include <array>
#include <eigen3/Eigen/Dense>

unsigned int factorial(unsigned int n){
    return (n==1 || n==0) ? 1: n * factorial(n - 1);
}

template<typename RealType, unsigned int N>
struct DerivativeCoef {
    RealType centralCoef;
    std::array<RealType, N> otherCoefs;
};

template<typename RealType, unsigned int N, unsigned int L>
DerivativeCoef<RealType, N> calcDerivativeCoef(const std::array<RealType, N>& points) noexcept{
    Eigen::Matrix<RealType, N, N> matrix;
    Eigen::Vector<RealType, N> b; 
    Eigen::Vector<RealType, N> ans;

    for (std::size_t i = 0; i < N; i++){
        b(i) = 0;
        matrix(0,i) = points[i];
    }
    b(L-1) = factorial(L);
    for(std::size_t i = 1; i < N; i++){
        for(std::size_t j = 0; j < N; j++)
            matrix(i, j) = matrix(i - 1, j) * points[j];
    }

    ans = matrix.householderQr().solve(b);
    RealType central = -ans.sum();
    std::array<RealType, N> otherArray;
    std::copy(ans.data(), ans.data() + ans.size(), otherArray.begin());

    return {central, otherArray};
}
