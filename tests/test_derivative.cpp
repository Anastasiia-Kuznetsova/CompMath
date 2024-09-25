#include "gtest/gtest.h"
#include "derivative.hpp"


TEST(DERIVATIVE, TEST_1){
    const unsigned int N = 2;
    std::array<double, 2> points{-1, 1};
    const double centralCoeff = 0.0;
    const std::array<double, 2> otherCoeff {-0.5, 0.5};
    DerivativeCoef<double, N> answer = calcDerivativeCoef<double, N, 1>(points);
    for (std::size_t i = 0; i < N; i++)
        ASSERT_NEAR(otherCoeff[i], answer.otherCoefs[i], 1e-13);
    ASSERT_NEAR(centralCoeff, answer.centralCoef, 1e-13);
}

TEST(DERIVATIVE, TEST_2){
    const unsigned int N = 2;
    std::array<double, 2> points{1, 2};
    const double centralCoeff = -1.5;
    const std::array<double, 2> otherCoeff {2, -0.5};
    DerivativeCoef<double, N> answer = calcDerivativeCoef<double, N, 1>(points);
    for (std::size_t i = 0; i < N; i++)
        ASSERT_NEAR(otherCoeff[i], answer.otherCoefs[i], 1e-13);
    ASSERT_NEAR(centralCoeff, answer.centralCoef, 1e-13);
}

TEST(DERIVATIVE, TEST_3){
    const unsigned int N = 2;
    const unsigned int L = 2;
    std::array<double, 2> points{-1, 1};
    const double centralCoeff =-2.0;
    const std::array<double, 2> otherCoeff {1, 1};
    DerivativeCoef<double, N> answer = calcDerivativeCoef<double, N, L>(points);
    for (std::size_t i = 0; i < N; i++)
        ASSERT_NEAR(otherCoeff[i], answer.otherCoefs[i], 1e-13);
    ASSERT_NEAR(centralCoeff, answer.centralCoef, 1e-13);
}
