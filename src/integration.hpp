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





/* Функция производит интегрирование на одном отрезке */
template<typename Callable, std::size_t N>
decltype(auto) integrate(   
    const Callable& func,  // Интегрируемая функция
    const typename ArgumentGetter<Callable>::Argument& start,  // начало отрезка
    const typename ArgumentGetter<Callable>::Argument& end  // конец отрезка
    ){
    std::array<typename ArgumentGetter<Callable>::Argument, 4> x;
	std::array<typename ArgumentGetter<Callable>::Argument, 4> weights;                       
    if (N == 3){
        x = {-0.77459666924148338, 0, 0.77459666924148338, 0};
        weights = {0.5555555555555556, 0.8888888888888889,0.5555555555555556, 0};
    } else if (N == 4){
        x = {-0.8611363115940526, -0.3399810435848562, 0.3399810435848562, 0.8611363115940526};
        weights = {0.3478548451374538, 0.8611363115940526, 0.8611363115940526, 0.3478548451374538};
    }
    typename ArgumentGetter<Callable>::Argument sum = 0;
	for(size_t i = 0; i < N; i++)
		sum += weights[i] * func(x[i] * (end - start) / 2 + (end + start) / 2);
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