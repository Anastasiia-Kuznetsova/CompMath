#include <array>
#include <vector>
#include <eigen3/Eigen/Dense>


/* Это таблица Бутчера для метода Рунге-Кутты 4 порядка. Я ее не заполнил */
struct RK4Table{
    static constexpr unsigned int stages = 4;
    static constexpr std::array<std::array<double, stages>, stages> table = {{{0., 0., 0., 0.}, {0.5, 0., 0., 0.}, {0., 0.5, 0., 0.}, {0., 0., 1., 0.}}};
    static constexpr std::array<double, stages> cColumn = {0., 0.5, 0.5, 1};
    static constexpr std::array<double, stages> bString = {1./6, 1./3, 1./3, 1./6};
};


class Oscillator {
public:
    static constexpr unsigned int dim = 2;  // размерность задачи
    using Argument = double;  // тип аргумента, тип t
    using State = Eigen::Vector<double, dim>;  // состояние
    
    struct StateAndArg{
        State state;
        Argument arg;
    };
    /*** Вычисляет правую часть ДУ - функцию f***/
    Eigen::Vector<double, dim> calc(const StateAndArg& stateAndArg) const {
        return Eigen::Vector<double, dim>{stateAndArg.state(1), -stateAndArg.state(0)};
    } 
};

class Task_1 {
public:
    static constexpr unsigned int dim = 1; // Размерность задачи
    using Argument = double; // тип аргумента, тип t
    using State = Eigen::Vector<double, dim>; // состояние
    
    struct StateAndArg {
        State state;
        Argument arg;
    };
     /*** Вычисляет правую часть ДУ - функцию f***/
    Eigen::Vector<double, dim> calc(const StateAndArg& stateAndArg) const {
        return Eigen::Vector<double, dim>{stateAndArg.arg * stateAndArg.arg *stateAndArg.arg};
    }
};



template<typename Table, typename RHS>  // таблица бутчера и класс правой части f
std::vector<typename RHS::StateAndArg> integrate(
    const typename RHS::StateAndArg& initialState,
    const typename RHS::Argument& endTime,
    double step,
    const RHS& rhs
){
    unsigned int n = std::ceil((endTime - initialState.arg) / step) + 1;
    std::vector<typename RHS::StateAndArg> result;
    result.reserve(n);
    result.push_back(initialState);
    
    std::array<typename RHS::State, Table::stages> k_coef;

    typename RHS::State tmp;
    for(std::size_t i = 0; i < n; i++) {
        k_coef[0] = rhs.calc(result[i]);
        for (std::size_t j = 1; j < Table::stages; j++) {
            tmp = result[i].state;
            for (std::size_t k = 0; k < j; k++) 
                tmp += Table::table[j][k] * k_coef[k] * step;
            k_coef[j] = rhs.calc(typename RHS::StateAndArg{tmp , result[i].arg + step * Table::cColumn[j]});
        }
        
        tmp = Table::bString[0] * k_coef[0];
        for (std::size_t j = 1; j < Table::stages; j++) 
            tmp += Table::bString[j] * k_coef[j];

        result.push_back({result[i].state + tmp * step, result[i].arg + step});
    }
    return result;
}
