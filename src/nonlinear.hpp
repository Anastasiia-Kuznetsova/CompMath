#pragma once
#include <type_traits>
#include <cmath>
#include <exception>


double keplerSolver(double ecc, double meanAnomaly, unsigned int maxIter, double tol){
    double prev = meanAnomaly;
    double cur = 0;
    for(std::size_t i = 0; i < maxIter; i++){
        cur = prev - (prev - ecc * std::sin(prev) - meanAnomaly) / (1 - ecc * std::cos(prev));
        if (std::abs(cur - prev) < tol) 
            return cur;
        prev = cur;
    }
    throw "Can't find solution";
}

template<typename A>
struct ArgumentGetter;

template<typename R, typename Arg>
struct ArgumentGetter<R(Arg)> {
    using Argument = Arg;
};

template<typename Callable, typename RealType>
decltype(auto) solve(   
    const Callable& func,                                             // функция F
    const RealType& tau,                                              // шаг тау
    const typename ArgumentGetter<Callable>::Argument& initialGuess,  // начальное приближение
    const unsigned int nIteration                                     // количество итераций
                    ){
    typename std::decay_t<typename ArgumentGetter<Callable>::Argument> cur = initialGuess;
    for(std::size_t i = 0; i < nIteration; i++){
        cur = cur + tau * func(cur);
    }
    return cur;
 }
