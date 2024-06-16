#include "algorithm.hpp"

auto main() -> int {
    double alpha = 1.0;
    double beta = 5.0;
    double evaporation = 0.5;
    int num_iterations = 1000;

    VRP data;
    data.read_file("./benchmarks/A/A-n32-k5.vrp");
    data.intialize_matricies();

    AntColony solver_vrp;
    auto best_solution = solver_vrp.find_solution(data, true);

    std::cout << "---" << std::endl;
    std::cout << "Final best solution:" << std::endl;
    std::cout << best_solution.second << std::endl;
    for (const auto& route : best_solution.first) {
        for (int node : route) {
            std::cout << node << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}