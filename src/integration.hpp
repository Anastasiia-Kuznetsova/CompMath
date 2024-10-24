#include <array>
#include <type_traits>

template<typename A>
struct ArgumentGetter;

template<typename R, typename Arg>
struct ArgumentGetter<R(Arg)> {
    using Argument = Arg;
};

template<typename T>
using Dif = decltype(std::declval<T>() - std::declval<T>());


template<std::size_t N>
constexpr std::array<double, N> points();

template<std::size_t N>
constexpr std::array<double, N> weights();

template<>
constexpr std::array<double, 3> points<3>(){
    return {-0.77459666924148338, 0, 0.77459666924148338};
}

template<>
constexpr std::array<double, 3> weights<3>(){
    return  {0.5555555555555556, 0.8888888888888889,0.5555555555555556};
}


template<>
constexpr std::array<double, 4> points<4>(){
    return {-0.8611363115940526, -0.3399810435848562, 0.3399810435848562, 0.8611363115940526};
}

template<>
constexpr std::array<double, 4> weights<4>(){
    return  {0.3478548451374538, 0.8611363115940526, 0.8611363115940526, 0.3478548451374538};
}


/* Функция производит интегрирование на одном отрезке */
template<typename Callable, std::size_t N>
decltype(auto) integrate(   
    const Callable& func,  // Интегрируемая функция
    const typename ArgumentGetter<Callable>::Argument& start,  // начало отрезка
    const typename ArgumentGetter<Callable>::Argument& end  // конец отрезка
    ){                    
    constexpr auto x = points<N>();
    constexpr auto w = weights<N>();
    typename ArgumentGetter<Callable>::Argument sum =  w[0] * func(x[0] * (end - start) / 2 + (end + start) / 2);
	for(size_t i = 1; i < N; i++)
		sum += w[i] * func(x[i] * (end - start) / 2 + (end + start) / 2);
    return (end - start) / 2 * sum;
    }
      

/* Функция производит интегрирование, разбивая отрезок на подотрезки длиной не более dx */
template<typename Callable, std::size_t N>
decltype(auto) integrate(   
    const Callable& func,  // Интегрируемая функция
    const typename ArgumentGetter<Callable>::Argument& start,  // начало отрезка
    const typename ArgumentGetter<Callable>::Argument& end,  // конец отрезка
    const Dif<typename ArgumentGetter<Callable>::Argument>& dx  // Длина подотрезка
                        ){
    typename ArgumentGetter<Callable>::Argument ans = 0;
	for(size_t i = 0; i < std::round((end - start) / dx); i++){
		ans += integrate<Callable, N>(func, start + i * dx, start + (i + 1) * dx);
	}
	ans += integrate<Callable, N>(func, start + round((end - start) / dx) * dx, end);
	return ans;
    
    }