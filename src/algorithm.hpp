#pragma once

#include <iostream>
#include <fstream>
#include <map>
#include <numeric>
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
            } else if (line.find("COMMENT") != std::string::npos) { //NOLINT
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
                nodes.resize(DIMENSION + 1);
            } else if (line.find("CAPACITY") != std::string::npos) {
                iss >> key >> key >> CAPACITY;
            } else if (line.find("NODE_COORD_SECTION") != std::string::npos) {
                while (std::getline(file, line) && line.find("DEMAND_SECTION") == std::string::npos) {
                    Node node;
                    iss.clear();
                    iss.str(line);
                    iss >> node.id >> node.x >> node.y;
                    nodes[node.id] = node;
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
        std::cout << "a" << std::endl;
        adjacency_matrix.resize(DIMENSION + 1, std::vector<double>(DIMENSION + 1, 0.0));
        pheromone_matrix.resize(DIMENSION + 1, std::vector<double>(DIMENSION + 1, 1.0));
        for (int i = 1; i < DIMENSION + 1; ++i) {
            for (int j = i + 1; j < DIMENSION + 1; ++j) {
                double distance = euclidean_distance(nodes[i], nodes[j]);
                adjacency_matrix[i][j] = distance;
                adjacency_matrix[j][i] = distance;
            }
        }
        std::cout << "Sucses intialize matricies" << std::endl;
    }

    void update_pheromone_map(const std::vector<std::pair<std::vector<std::vector<int>>, double>>& solutions) {
        for (size_t i = 1; i < pheromone_matrix.size(); ++i) {
            for (size_t j = i + 1; j < pheromone_matrix.size(); ++j) {
                pheromone_matrix[i][j] = std::max((1 - EVAPORATION) * pheromone_matrix[i][j], 1e-10);
                pheromone_matrix[j][i] = pheromone_matrix[i][j];
            }
        }

        for (const auto& solution : solutions) {
            double pheromone_increase = 1.0 / solution.second;
            for (const auto& route : solution.first) {
                for (size_t i = 0; i < route.size() - 1; ++i) {
                    int from = route[i];
                    int to = route[i + 1];
                    pheromone_matrix[from][to] += pheromone_increase;
                    pheromone_matrix[to][from] += pheromone_increase;
                }
            }
        }
    }

    void global_update_pheromone_map(const std::pair<std::vector<std::vector<int>>, double>& ant_solution, double capacity) {
        for (size_t i = 1; i < pheromone_matrix.size(); ++i) {
            for (size_t j = i + 1; j < pheromone_matrix.size(); ++j) {
                pheromone_matrix[i][j] = std::max((1 - EVAPORATION) * pheromone_matrix[i][j], 1e-10);
                pheromone_matrix[j][i] = pheromone_matrix[i][j];
            }
        }

        double pheromone_increase = capacity / ant_solution.second;
        for (const auto& route : ant_solution.first) {
            for (size_t i = 0; i < route.size() - 1; ++i) {
                int from = route[i];
                int to = route[i + 1];
                pheromone_matrix[from][to] += pheromone_increase;
                pheromone_matrix[to][from] += pheromone_increase;
            }
        }
    }

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
    explicit Ant(VRP& data)
        : p_data(data){
        reset_state();
    }

    auto select_next_point(int current_point) -> int {
        std::vector<int> available_points;
        for (int point : points_left) {
            if (capacity >= p_data.demands[point]) {
                available_points.push_back(point);
            }
        }

        if (available_points.empty()) {
            return -1;
        }

        std::vector<double> scores;
        for (int point : available_points) {
            double score = std::pow(p_data.pheromone_matrix[current_point][point], p_data.ALPHA) *
                           std::pow(1.0 / (p_data.adjacency_matrix[current_point][point] + 1e-10), p_data.BETA);
            scores.push_back(score);
        }

        double denominator = std::accumulate(scores.begin(), scores.end(), 0.0);
        std::vector<double> probabilities;
        for (double score : scores) {
            probabilities.push_back(score / denominator);
        }

        std::discrete_distribution<int> distribution(probabilities.begin(), probabilities.end());
        std::random_device rd;
        std::mt19937 gen(rd());
        int next_point = available_points[distribution(gen)];

        return next_point;
    }

    void move_to_point(int current_point, int next_point) {
        routes.back().push_back(current_point);
        if (next_point != 1) {
            points_left.erase(next_point);
        }
        capacity -= p_data.demands[next_point];
        total_path += p_data.adjacency_matrix[current_point][next_point];
    }

    void start_new_route() {
        capacity = p_data.CAPACITY;
        routes.push_back({1});

        std::vector<int> possible_points;
        for (int point : points_left) {
            if (capacity >= p_data.demands[point]) {
                possible_points.push_back(point);
            }
        }

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(0, static_cast<int>(possible_points.size()) - 1);
        int first_point = possible_points[dis(gen)];
        move_to_point(1, first_point);
    }

    auto find_solution() -> std::pair<std::vector<std::vector<int>>, double> {
        start_new_route();

        while (!points_left.empty()) {
            int current_point = routes.back().back();
            int next_point = select_next_point(current_point);
            if (next_point == -1) {
                move_to_point(current_point, 1);
                start_new_route();
            } else {
                move_to_point(current_point, next_point);
            }
        }

        move_to_point(routes.back().back(), 1);

        return {routes, total_path};
    }

    void reset_state() {
        capacity = p_data.CAPACITY;
        points_left.clear();
        for (size_t i = 1; i < p_data.DIMENSION + 1; ++i) {
            points_left.insert(static_cast<int>(i));
        }
        points_left.erase(1);
        routes.clear();
        total_path = 0.0;
    }

private:
    VRP& p_data;
    int capacity;
    std::set<int> points_left;
    std::vector<std::vector<int>> routes;
    double total_path;
};


class AntColony {
public:
    AntColony() = default;

    static auto get_route_cost(const std::vector<int>& route, const VRP& data) -> double {
        double total_cost = 0.0;
        for (size_t i = 0; i < route.size() - 1; ++i) {
            total_cost += data.adjacency_matrix[route[i]][route[i + 1]];
        }
        return total_cost;
    }

    auto get_route_cost_with_(const std::vector<int>& route, const VRP& data) -> double {
        double depot_costs = data.adjacency_matrix[1][route[0]] + data.adjacency_matrix[route.back()][1];
        return depot_costs + get_route_cost(route, data);
    }

    auto find_solution(VRP& vrp, bool verbose = true) -> std::pair<std::vector<std::vector<int>>, double> {
        std::vector<Ant> ants;
        for (int i = 0; i < vrp.NUMBER_OF_TRUCKS; ++i) {
            ants.emplace_back(vrp);
        }

        std::pair<std::vector<std::vector<int>>, double> best_solution;
        double best_cost = std::numeric_limits<double>::max();
        int tolerance = 0;

        for (int iter = 1; iter <= vrp.ITERATIONS_NUM; ++iter) {
            for (auto& ant : ants) {
                ant.reset_state();
            }

            std::vector<std::pair<std::vector<std::vector<int>>, double>> solutions;
            for (auto& ant : ants) {
                solutions.push_back(ant.find_solution());
            }

            auto candidate_best_solution = *std::min_element(solutions.begin(), solutions.end(), 
                [](const auto& lhs, const auto& rhs) {
                    return lhs.second < rhs.second;
                });

            vrp.update_pheromone_map(solutions);
            vrp.global_update_pheromone_map(candidate_best_solution, vrp.CAPACITY);

            if (candidate_best_solution.second < best_cost) {
                best_solution = candidate_best_solution;
                best_cost = candidate_best_solution.second;
                tolerance = 0;
            } else {
                tolerance += 1;
            }

            if (verbose && iter % 100 == 0) {
                std::cout << "Best solution in iteration " << iter << "/" << vrp.ITERATIONS_NUM << " = " << best_cost << std::endl;
            }

            if (tolerance >= static_cast<int>(0.3 * vrp.ITERATIONS_NUM)) {
                if (verbose) {
                    std::cout << "---" << std::endl;
                    std::cout << "Final best solution:" << std::endl;
                    std::cout << best_cost << std::endl;
                    for (const auto& route : best_solution.first) {
                        for (int node : route) {
                            std::cout << node << " ";
                        }
                        std::cout << std::endl;
                    }
                }
                return best_solution;
            }
        }

        if (verbose) {
            std::cout << "---" << std::endl;
            std::cout << "Final best solution:" << std::endl;
            std::cout << best_cost << std::endl;
            for (const auto& route : best_solution.first) {
                for (int node : route) {
                    std::cout << node << " ";
                }
                std::cout << std::endl;
            }
        }

        return best_solution;
    }

private:
    std::vector<std::vector<int>> routes;
    double cost = 0.0;
};