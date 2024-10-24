#pragma once
#include "ThreeDiagonalMatrix.hpp"


template<typename Type>
struct SplinePolynomCoeffs {
    Type a, b, c, d;
};

template<typename xType, typename yType>
class CubicSpline {
private:
    std::vector<xType> points_;
    std::vector<SplinePolynomCoeffs<DivisType<xType, yType>>> coeffs_;

public:
    CubicSpline( const std::vector<xType> &points,  const std::vector<yType>& values): points_(points){
        std::size_t n = points_.size();
        coeffs_.resize(n - 1);
        std::vector<xType> h(n - 1);
        std::vector<yType> u(n - 1);
        std::vector<DivisType<xType, yType>> column(n - 2);
        
        for(std::size_t i = 0; i < n - 1; i++){
            h[i] = points[i + 1] - points[i];
            u[i] = (values[i + 1] - values[i]) / h[i];
        }
        
        for (std::size_t i = 0; i < n - 2; i++) 
            column[i] = 6 * (u[i + 1] - u[i]) / (h[i + 1] + h[i]);

        std::vector<xType> b(n - 2, 2);
        std::vector<xType> a(n - 3), c(n - 3);
        for (std::size_t i = 0; i < n - 3; i++) {
            c[i] = h[i + 1] / (h[i + 1] + h[i]);
            a[i] = h[i] / (h[i + 1] + h[i]);
        }
       
        ThreeDiagonalMatrix<xType> matrix(a, b, c);
        std::vector<DivisType<xType, yType>> C_coeffs = solve(matrix, column);

        for (std::size_t i = 0; i < n - 1; i++) 
            coeffs_[i].a = values[i + 1];
         
        coeffs_[n - 2].c = 0;
        for (std::size_t i = 0; i < n - 2; i++)
            coeffs_[i].c = C_coeffs[i];
      
        coeffs_[0].d = coeffs_[0].c / h[0];
        for (std::size_t i = 1; i < n-1; i++) 
            coeffs_[i].d = (coeffs_[i].c - coeffs_[i - 1].c) / h[i];

        for (std::size_t i = 0; i < n - 1; i++) 
            coeffs_[i].b = coeffs_[i].c / 2 * h[i] - coeffs_[i].d/6 * h[i] * h[i]+ u[i];
    }
                        
    yType interpolate(const xType& x) const noexcept{
        for (std::size_t i = 0; i < coeffs_.size(); ++i) {
            if (points_[i] <= x && x <= points_[i + 1]) 
                return coeffs_[i].a + (x - points_[i + 1]) * (coeffs_[i].b  +  (x - points_[i + 1]) * (coeffs_[i].c / 2 + (x - points_[i + 1]) * coeffs_[i].d / 6));
        }
    }
};