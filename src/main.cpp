#include "algorithm.hpp"
#include <filesystem>
#include <vector>

namespace fs = std::filesystem;

void save_best_solution(const std::string& filename, const std::pair<std::vector<std::vector<int>>, double>& best_solution) {
    std::ofstream outfile("./solutions/" + filename + ".sol");

    int route_number = 1;
    for (const auto& route : best_solution.first) {
        if (route.size() > 1) {
            outfile << "Route #" << route_number << ": ";
            for (size_t i = 1; i < route.size() - 1; ++i) { // Начинаем с 1, чтобы не включать депо (1), и по этой же причине не доходим до самого конца
                outfile << route[i] << " ";
            }
            outfile << std::endl;
            route_number++;
        }
    }

    outfile << "cost " << best_solution.second << std::endl;
    outfile.close();
}

void run_vrp_on_file(const std::string& filepath, const std::vector<std::tuple<double, double, double>>& params, int num_iterations) {
    double best_cost = std::numeric_limits<double>::max();
    std::pair<std::vector<std::vector<int>>, double> best_solution;

    double best_alpha;
    double best_beta;
    double best_evaporation;

    std::cout << "File: " << filepath << std::endl;

    std::cout << "Result without 2-opt on 10000 iterations, alpha = beta = 0.5, evaporation = 0.2" << std::endl;
    VRP data_without_2_opt(0.5, 0.5, 0.25, num_iterations, false);
    data_without_2_opt.read_file(filepath);
    data_without_2_opt.initialize_matrices();

    AntColony solver_vrp_without_2_opt;
    auto best_solution_without_2_opt = solver_vrp_without_2_opt.find_solution(data_without_2_opt);
    std::cout << "Cost:" << best_solution_without_2_opt.second << std::endl;
    for (const auto& route : best_solution.first) {
        for (int node : route) {
            std::cout << node << " ";
        }
        std::cout << std::endl;
    }

    for (const auto& param: params) {
        double alpha;
        double beta;
        double evaporation;
        std::tie(alpha, beta, evaporation) = param;
        VRP data(alpha, beta, evaporation, num_iterations, true);
        data.read_file(filepath);
        data.initialize_matrices();

        std::cout << "alpha = " << alpha << "beta = " << beta << "evaporation = " << evaporation << std::endl;

        AntColony solver_vrp;
        auto solution = solver_vrp.find_solution(data);

        if (solution.second < best_cost) {
            best_cost = solution.second;
            best_solution = solution;
            best_alpha = alpha;
            best_beta = beta;
            best_evaporation = evaporation;
        }
    }
    std::cout << "Best cost: " << best_cost << std::endl;
    std::cout << "Best params: alpha=" << best_alpha << ", beta=" << best_beta << ", evaporation=" << best_evaporation << std::endl;
    for (const auto& route : best_solution.first) {
        for (int node : route) {
            std::cout << node << " ";
        }
        std::cout << std::endl;
    }

    std::string filename = fs::path(filepath).filename().string();
    save_best_solution(filename, best_solution);
}

auto main() -> int {
    std::vector<std::tuple<double, double, double>> params {{1.0, 1.5, 0.5},
                                                     {1.5, 1.0, 0.25},
                                                     {1.5, 2.0, 0.25},
                                                     {1.5, 0.5, 0.5},
                                                     {1.5, 1.0, 0.5},
                                                     {1.5, 1.5, 0.5},
                                                     {2.0, 0.5, 0.25},
                                                     {2.0, 1.0, 0.25}};
    int num_iterations = 10000;

    std::vector<std::string> directories = {"./benchmarks/A", "./benchmarks/B"};

    for (const auto& dir : directories) {
        for (const auto& entry : fs::directory_iterator(dir)) {
            if (entry.is_regular_file() && entry.path().extension() == ".vrp") {
                run_vrp_on_file(entry.path().string(), params, num_iterations);
            }
        }
    }

    return 0;
}