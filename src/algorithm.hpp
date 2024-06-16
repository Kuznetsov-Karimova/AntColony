#pragma once

#include <iostream>
#include <fstream>
#include <map>
#include <random>
#include <set>
#include <string>
#include <vector>
#include <sstream>
#include <cmath>

struct Node {
    int id;
    int x;
    int y;
};


auto euclidean_distance(const Node& a, const Node& b) -> double {
    return std::sqrt(std::pow(b.x - a.x, 2) + std::pow(b.y - a.y, 2));
}


class VRP {
public:
    explicit VRP(double alpha = 0.5, double beta = 0.5, double evaporation = 0.25, double iterations_num = 1000) {
        ALPHA = alpha;
        BETA = beta;
        EVAPORATION = evaporation;
        ITERATIONS_NUM = iterations_num;
    }

    void read_file(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Unable to open file" << std::endl;
            return;
        }

        std::string line;

        while (std::getline(file, line)) {
            std::istringstream iss(line);
            std::string key;

            if (line.find("NAME") != std::string::npos) {
                continue;
            } else if (line.find("COMMENT") != std::string::npos) {
                size_t pos1 = line.find("No of trucks: ");
                size_t pos2 = line.find(", Optimal value:");
                if (pos1 != std::string::npos && pos2 != std::string::npos) {
                    NUMBER_OF_TRUCKS = std::stoi(line.substr(pos1 + 13, pos2 - pos1 - 13));
                }
            } else if (line.find("TYPE") != std::string::npos) {
                continue;
            } else if (line.find("DIMENSION") != std::string::npos) {
                iss >> key >> key >> DIMENSION;
                demands.resize(DIMENSION+1);
            } else if (line.find("CAPACITY") != std::string::npos) {
                iss >> key >> key >> CAPACITY;
            } else if (line.find("NODE_COORD_SECTION") != std::string::npos) {
                while (std::getline(file, line) && line.find("DEMAND_SECTION") == std::string::npos) {
                    Node node;
                    iss.clear();
                    iss.str(line);
                    iss >> node.id >> node.x >> node.y;
                    nodes.push_back(node);
                }
            }

            if (line.find("DEMAND_SECTION") != std::string::npos) {
                while (std::getline(file, line) && line.find("DEPOT_SECTION") == std::string::npos) {
                    int id;
                    int demand;
                    iss.clear();
                    iss.str(line);
                    iss >> id >> demand;
                    demands[id] = demand;
                }
            }

            if (line.find("DEPOT_SECTION") != std::string::npos) {
                break;
            }
        }

        file.close();
    }

    void intialize_matricies() {
        adjacency_matrix.resize(DIMENSION, std::vector<double>(DIMENSION, 0.0));
        pheromone_matrix.resize(DIMENSION, std::vector<double>(DIMENSION, 1.0));
        for (int i = 1; i < DIMENSION + 1; ++i) {
            for (int j = i + 1; j < DIMENSION + 1; ++j) {
                double distance = euclidean_distance(nodes[i], nodes[j]);
                adjacency_matrix[i][j] = distance;
                adjacency_matrix[j][i] = distance;
            }
        }
    }

    [[nodiscard]] auto get_distance(int from, int to) const -> double {
        return adjacency_matrix[from][to];
    }

    [[nodiscard]] auto get_pheromone(int from, int to) const -> double {
        return pheromone_matrix[from][to];
    }

    void add_pheromone(int from, int to, double amount) {
        pheromone_matrix[from][to] += amount;
    }

    void evaporate_pheromone(int from, int to, double evaporation_rate) {
        pheromone_matrix[from][to] *= (1.0 - evaporation_rate);
    }

    [[nodiscard]] auto get_demand(int node) const -> int {
        return demands[node];
    }

    [[nodiscard]] auto get_dimension() const -> int {
        return DIMENSION;
    }

    [[nodiscard]] auto get_capacity() const -> int {
        return CAPACITY;
    }

    void initialize_pheromone_matrix() {
        pheromone_matrix.resize(DIMENSION + 1, std::vector<double>(DIMENSION + 1, 1.0));
    }

private:
    double ALPHA;
    double BETA;
    double EVAPORATION;
    double ITERATIONS_NUM;

    int DIMENSION;
    int NUMBER_OF_TRUCKS;
    int CAPACITY;

    std::vector<Node> nodes;
    std::vector<int> demands;

    std::vector<std::vector<double>> adjacency_matrix;
    std::vector<std::vector<double>> pheromone_matrix;
};


class Ant {
public:
    Ant(VRP& graph, int capacity, double alpha, double beta)
        : graph(graph), ant_capacity(capacity), alpha(alpha), beta(beta) {
        reset_state();
    }

    auto find_solution() -> std::vector<std::vector<int>> {
        start_new_route();
        while (!cities_left.empty()) {
            int current_city = routes.back().back();
            int next_city = select_next_city(current_city);
            if (next_city == -1) {
                move_to_city(current_city, 1); // move to depot
                start_new_route();
            } else {
                move_to_city(current_city, next_city);
            }
        }
        move_to_city(routes.back().back(), 1); // complete last route
        return routes;
    }

    void reset_state() {
        capacity = ant_capacity;
        cities_left.clear();
        for (int i = 2; i <= graph.get_dimension(); ++i) {
            cities_left.insert(i);
        }
        routes.clear();
        total_path_cost = 0;
    }

    [[nodiscard]] auto get_total_path_cost() const -> double {
        return total_path_cost;
    }

private:
    VRP& graph;
    int ant_capacity;
    double alpha;
    double beta;

    int capacity;
    std::set<int> cities_left;
    std::vector<std::vector<int>> routes;
    double total_path_cost;

    void start_new_route() {
        capacity = ant_capacity;
        routes.emplace_back(1, 1);
        int first_city = select_random_city();
        move_to_city(1, first_city);
    }

    auto select_random_city() -> int {
        std::vector<int> available_cities;
        for (int city : cities_left) {
            if (capacity >= graph.get_demand(city)) {
                available_cities.push_back(city);
            }
        }
        if (available_cities.empty()) { return -1; }

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(0, available_cities.size() - 1);
        return available_cities[dis(gen)];
    }

    auto select_next_city(int current_city) -> int {
        std::vector<int> available_cities;
        std::vector<double> probabilities;
        double total_pheromone = 0.0;

        for (int city : cities_left) {
            if (capacity >= graph.get_demand(city)) {
                double pheromone = std::pow(graph.get_pheromone(current_city, city), alpha);
                double heuristic = std::pow(1.0 / (graph.get_distance(current_city, city) + 1e-10), beta);
                double score = pheromone * heuristic;
                available_cities.push_back(city);
                probabilities.push_back(score);
                total_pheromone += score;
            }
        }

        if (available_cities.empty()) { return -1; }

        std::random_device rd;
        std::mt19937 gen(rd());
        std::discrete_distribution<> dis(probabilities.begin(), probabilities.end());

        return available_cities[dis(gen)];
    }

    void move_to_city(int current_city, int next_city) {
        routes.back().push_back(next_city);
        if (next_city != 1) { cities_left.erase(next_city); }
        capacity -= graph.get_demand(next_city);
        total_path_cost += graph.get_distance(current_city, next_city);
    }
};


class AntColonyOptimization {
public:
    AntColonyOptimization(VRP& data, int num_ants, int num_iterations, double alpha, double beta, double evaporation)
        : data(data), num_ants(num_ants), num_iterations(num_iterations), alpha(alpha), beta(beta), evaporation(evaporation) {
        initialize_pheromone_matrix();
    }

    void run() {
        for (int iteration = 0; iteration < num_iterations; ++iteration) {
            std::vector<std::vector<std::vector<int>>> all_solutions;
            std::vector<double> all_costs;

            for (int i = 0; i < num_ants; ++i) {
                Ant ant(data, data.get_capacity(), alpha, beta);
                std::vector<std::vector<int>> solution = ant.find_solution();
                double cost = ant.get_total_path_cost();
                all_solutions.push_back(solution);
                all_costs.push_back(cost);
                update_pheromones(solution, cost);
            }

            evaporate_pheromones();
            auto min_cost_it = std::min_element(all_costs.begin(), all_costs.end());
            double min_cost = *min_cost_it;
            std::vector<std::vector<int>> best_solution = all_solutions[std::distance(all_costs.begin(), min_cost_it)];

            if (min_cost < best_tour_length) {
                best_tour_length = min_cost;
                best_tour = best_solution;
            }

            std::cout << "Iteration " << iteration + 1 << ": Best tour length = " << best_tour_length << std::endl;
        }

        print_best_tour();
    }

private:
    VRP& data;
    int num_ants;
    int num_iterations;
    double alpha;
    double beta;
    double evaporation;
    double best_tour_length = std::numeric_limits<double>::max();
    std::vector<std::vector<int>> best_tour;

    void initialize_pheromone_matrix() {
        data.initialize_pheromone_matrix();
    }

    void update_pheromones(const std::vector<std::vector<int>>& solution, double cost) {
        double pheromone_deposit = 1.0 / cost;
        for (const auto& route : solution) {
            for (size_t i = 0; i < route.size() - 1; ++i) {
                int from = route[i];
                int to = route[i + 1];
                data.add_pheromone(from, to, pheromone_deposit);
            }
        }
    }

    void evaporate_pheromones() {
        for (int i = 1; i <= data.get_dimension(); ++i) {
            for (int j = 1; j <= data.get_dimension(); ++j) {
                data.evaporate_pheromone(i, j, evaporation);
            }
        }
    }

    void print_best_tour() {
        std::cout << "Best tour length: " << best_tour_length << std::endl;
        for (const auto& route : best_tour) {
            for (int node : route) {
                std::cout << node << " ";
            }
            std::cout << std::endl;
        }
    }
};