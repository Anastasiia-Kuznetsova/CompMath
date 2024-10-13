#include <vector>
#include <type_traits>

template<typename Type>
class ThreeDiagonalMatrix {
private:
    std::vector<Type> a_;
    std::vector<Type> b_;
    std::vector<Type> c_;
    std::size_t matrix_size_;  
public:
    ThreeDiagonalMatrix(const std::vector<Type> &a,
                      const std::vector<Type> &b,
                      const std::vector<Type> &c) : a_(a), b_(b), c_(c), matrix_size_(b.size()){}
    
    std::size_t size() const {return matrix_size_;}
    
    Type a(std::size_t i) const{return a_[i];}
    Type b(std::size_t  i) const{return b_[i];}
    Type c(std::size_t i) const{return c_[i];} 
};

template<typename numeratorType, typename denominatorType>
using DivisType = decltype(std::declval<numeratorType>() / std::declval<denominatorType>());


template<typename mType, typename cType>
std::vector<DivisType<cType, mType>> solve( const ThreeDiagonalMatrix<mType>& matrix,
                                            const std::vector<cType>& column){

    std::size_t n = matrix.size();
    std::vector<mType> x(n);
    std::vector<mType> p(n);
    std::vector<mType> q(n);

    p[1] = -matrix.c(0) / matrix.b(0);
    q[1] = column[0] /  matrix.b(0);
    for (std::size_t i = 1; i + 1  < (n); i++) {
       p[i+1] = -matrix.c(i) / (matrix.a(i-1) * p[i] + matrix.b(i));
       q[i+1]  = (column[i] - matrix.a(i - 1) * q[i]) / (matrix.a(i - 1) * p[i] + matrix.b(i));
    }

    x[n - 1]  = (column[n - 1] - matrix.a(n - 2) * q[n-1]) / (matrix.a(n - 2) * p[n-1] + matrix.b(n - 1));
    for(std::size_t i = n - 2; i >= 1; i--){
        x[i] = p[i+1] * x[i + 1] + q[i+1];
    }
    x[0] = p[1] * x[1] + q[1];

    return x;
}