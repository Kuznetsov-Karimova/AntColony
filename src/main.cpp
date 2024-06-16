#include "algorithm.hpp"

auto main() -> int {
    // Параметры алгоритма
    double alpha = 1.0;
    double beta = 5.0;
    double evaporation = 0.5;
    int num_iterations = 1000;
    int num_ants = 10;
    
    // Создание объекта данных
    VRP data;
    data.read_file("./benchmarks/A/A-n32-k5.vrp");
    data.intialize_matricies();

    // Создание объекта муравьиной колонии
    AntColonyOptimization aco(data, num_ants, num_iterations, alpha, beta, evaporation);

    // Запуск алгоритма
    aco.run();

    return 0;
}