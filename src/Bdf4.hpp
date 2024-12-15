#include "Rk4.hpp"
#include <eigen3/Eigen/Dense>

/* Это коэффициенты для метода на 4 шагах. Я их не заполнил */
struct BDF4{
    static constexpr unsigned int size = 4;
    static constexpr std::array<double, size + 1> alpha = {-3./25, 16./25, -36./25, 48./25, 12./25};
};


class Oscillator_bdf {
public:
    
    static constexpr unsigned int dim = 2;  // размерность задачи
    
    using Argument = double;  // тип аргумента, тип t
    
    using State = Eigen::Vector<double, dim>;  // состояние
    
    struct StateAndArg{
        State state;
        Argument arg;
    };

    /*** Вычисляет разницу двух состояний (решения нелиненого уравнения) ***/
    State calcDif(const State& first, const State& second) const{
        return second - first;
    }
    
    /*** Вычисляет правую часть ДУ - функцию f***/
    Eigen::Vector<double, dim> calc(const StateAndArg& stateAndArg) const {
        return Eigen::Vector<double, dim>{stateAndArg.state(1), -stateAndArg.state(0)};
    } 
};


class Task_1_bdf {
public:
    static constexpr unsigned int dim = 1; // Размерность задачи
    using Argument = double; // тип аргумента, тип t
    using State = Eigen::Vector<double, dim>; // состояние
    
    struct StateAndArg {
        State state;
        Argument arg;
    };

    /*** Вычисляет разницу двух состояний (решения нелиненого уравнения) ***/
    State calcDif(const State& first, const State& second) const{
        return second - first;
    }
     /*** Вычисляет правую часть ДУ - функцию f***/
    Eigen::Vector<double, dim> calc(const StateAndArg& stateAndArg) const {
        return Eigen::Vector<double, dim>{stateAndArg.arg * stateAndArg.arg *stateAndArg.arg};
    }
};

struct IntegrationParameters {
    double step; // шаг интегрирования
    double epsilon; // точность решения нелинейного уравнения
    unsigned int maxIter; // максимальное количество итераций для
                          // решения нелинейного уравнения
};

const double mu = 3.9860044189e14;
const double epsilon = 3/2 * mu * 6378137 * 1.082e-3;

class Orbit {
public:
    static constexpr unsigned int dim = 6; // Размерность задачи
    using Argument = double; // тип аргумента, тип t
    using State = Eigen::Vector<double, dim>; // состояние
    
    struct StateAndArg {
        State state;
        Argument arg;
    };

    /*** Вычисляет разницу двух состояний (решения нелиненого уравнения) ***/
    State calcDif(const State& first, const State& second) const{
        return second - first;
    }
     /*** Вычисляет правую часть ДУ - функцию f***/
    Eigen::Vector<double, dim> calc(const StateAndArg& stateAndArg) const {
        return Eigen::Vector<double, dim>{stateAndArg.state[3], stateAndArg.state[4], stateAndArg.state[5],
                     -mu / (std::pow(stateAndArg.state[0],2)) - epsilon / (3 * std::pow(stateAndArg.state[0],4)) * (3 * std::pow(std::sin(stateAndArg.state[2]), 2) - 1),
                     0, 
                     -epsilon / (std::pow(stateAndArg.state[0], 4)) * std::sin(2*stateAndArg.state[2])};
    }
};


/***
BDF - структура с коэффициентами метода
RHS - правая часть Д.У.
RKTable -  таблица Бутчера метода Рунге-Кутты для разгона.
***/

template <typename BDF, typename RHS,
          typename RKTable> // таблица бутчера и класс правой части f
std::vector<typename RHS::StateAndArg>
integrate(const typename RHS::StateAndArg &initialState,
          const typename RHS::Argument &endTime,
          const IntegrationParameters &parameters, const RHS &rhs) {

    unsigned int n = std::ceil((endTime - initialState.arg) / parameters.step);
    std::vector<typename RHS::StateAndArg> res = integrate<RKTable, RHS>(initialState, initialState.arg + parameters.step * (BDF::size - 1), parameters.step, rhs);
    double t = initialState.arg + parameters.step * 4;
    Eigen::Vector<double, RHS::dim> tmp;
    for (std::size_t i = BDF::size; i <= n; i++) {
        for (std::size_t j = 0; j < parameters.maxIter; j++) {
                const typename RHS::StateAndArg a = {res.back().state, t};
                tmp = parameters.step * rhs.calc(a) * BDF::alpha.back();
                for (std::size_t k = 0; k < BDF::size; ++k) {
                    tmp += res[res.size() - 4 + k].state * BDF::alpha[k];
                }
            if (rhs.calcDif(tmp, res.back().state).norm() < parameters.epsilon) 
                 break;
         }
       res.push_back({tmp, t});
       t += parameters.step;
    }
    return res;
}