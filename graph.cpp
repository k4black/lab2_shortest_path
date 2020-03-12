#ifndef MIN_CUT_GRAPH
#define MIN_CUT_GRAPH

#include <vector>
#include <random>
#include <utility>
#include <algorithm>
#include <unordered_map>
#include <functional>
#include <memory>
#include <set>
#include <iostream>


enum Colour {
    BLACK, GREY, WHITE
};


class Edge {
public:
    size_t first;
    size_t second;
    int32_t weight;

    Edge(size_t first, size_t second, int32_t weight) {
        this->first = first;
        this->second = second;
        this->weight = weight;
    }

    bool operator< (const Edge &a) const {
        return this->first < a.first || this->first == a.first && this->second < a.second;
    }
};



class Graph {
private:
    std::vector<std::vector<size_t>> matrix;  // matrix
    std::vector<std::vector<std::pair<size_t, int32_t>>> adjacency;  // adjacency list with weight

    size_t num_of_nodes = 0;

public:
    Graph() = default;

    explicit Graph(size_t num_of_nodes) {
        reserve(num_of_nodes);
    }

    void reserve(size_t num_of_nodes) {
        this->num_of_nodes = num_of_nodes;

        matrix.clear();
        matrix.resize(num_of_nodes);

        for (auto &line: matrix) {
            line.resize(num_of_nodes);
            std::fill(line.begin(), line.end(), INT64_MAX);
        }

        adjacency.clear();
        adjacency.resize(num_of_nodes);
    }

    void build(const std::vector<Edge> &edges) {
        for (auto &edge : edges) {
            set_neighbour(edge.first, edge.second, edge.weight);
        }
    }

    void build(const std::set<Edge> &edges) {
        for (auto &edge : edges) {
            set_neighbour(edge.first, edge.second, edge.weight);
        }
    }

    void build_directed(const std::vector<Edge> &edges) {
        for (auto &edge : edges) {
            set_neighbour_directed(edge.first, edge.second, edge.weight);
        }
    }

    [[nodiscard]] const std::vector<std::pair<size_t, int32_t>>& operator[](size_t node_index) const {
        return adjacency[node_index];
    }

    [[nodiscard]] const std::vector<std::pair<size_t, int32_t>>&get_neighbours(size_t node_index) const {
        return adjacency[node_index];
    }


    void remove_neighbour(size_t node_index, size_t neighbour_index) {
        remove_neighbour_directed(node_index, neighbour_index);
        remove_neighbour_directed(neighbour_index, node_index);
    }

    void remove_neighbour_directed(size_t node_index, size_t neighbour_index) {
        matrix[node_index][neighbour_index] = INT64_MAX;

        auto position = adjacency[node_index].begin();
        for (; position != adjacency[node_index].end(); ++position) {
            if (position->first == neighbour_index) {
                break;
            }
        }
        if (position != adjacency[node_index].end()) {  // == adjacency.end() means the element was not found
            adjacency[node_index].erase(position);
        }
    }

    void set_neighbour(size_t node_index, size_t neighbour_index, int32_t weight) {
        set_neighbour_directed(node_index, neighbour_index, weight);
        set_neighbour_directed(neighbour_index, node_index, weight);
    }

    void set_neighbour_directed(size_t node_index, size_t neighbour_index, int32_t weight) {
        matrix[node_index][neighbour_index] = weight;

        auto position = adjacency[node_index].begin();
        for (; position != adjacency[node_index].end(); ++position) {
            if (position->first == neighbour_index) {
                break;
            }
        }
        if (position == adjacency[node_index].end()) {  // == adjacency.end() means the element was not found
            adjacency[node_index].emplace_back(neighbour_index, weight);
        }
    }

    void add_neighbour(size_t node_index, size_t neighbour_index, int32_t weight) {
        set_neighbour(node_index, neighbour_index, weight);
    }

    void add_neighbour_directed(size_t node_index, size_t neighbour_index, size_t weight) {
        set_neighbour_directed(node_index, neighbour_index, weight);
    }

    [[nodiscard]] bool check_neighbour(size_t node_index, size_t neighbour_index) const {
        return matrix[node_index][neighbour_index] != INT64_MAX;
    }

    void revert_edge(size_t node_index, size_t neighbour_index) {
        std::swap(matrix[node_index][neighbour_index], matrix[neighbour_index][node_index]);

        remove_neighbour(node_index, neighbour_index);
        set_neighbour_directed(node_index, neighbour_index, matrix[node_index][neighbour_index]);
        set_neighbour_directed(neighbour_index, node_index, matrix[neighbour_index][node_index]);
    }

    [[nodiscard]] size_t node_degree(size_t node_index) const {
        return adjacency[node_index].size();
    }


    class iterator {
    public:
        size_t current_it = 0;

        iterator(size_t current) {
            current_it = current;
        }

        bool operator==(const iterator &other) {
            return this->current_it == other.current_it;
        }

        void operator=(const iterator &other) {
            current_it = other.current_it;
        }

        bool operator!=(const iterator &other) {
            return !this->operator==(other);
        }

        iterator& operator++() {
            ++current_it;
            return *this;
        }

        size_t operator*() const {
            return current_it;
        }
    };

    [[nodiscard]] iterator begin() noexcept {
        return iterator(0);
    }

    [[nodiscard]] iterator end() noexcept {
        return iterator(num_of_nodes);
    }

    [[nodiscard]] iterator begin() const noexcept {
        return iterator(0);
    }

    [[nodiscard]] iterator end() const noexcept {
        return iterator(num_of_nodes);
    }
//
//    // available since C++ 2011.03
//    [[nodiscard]] const_iterator cbegin() const noexcept {
//        return iterator(0);
//    }
//
//    [[nodiscard]] const_iterator cend() const noexcept {
//        return iterator(num_of_nodes);
//    }

    [[nodiscard]] size_t size() const {
        return num_of_nodes;
    }
};


class GraphGenerator {
public:

//    static void build(Graph &graph, const std::vector<std::pair<size_t, size_t>> &edges) {
//        graph.build(edges);
//    }

    GraphGenerator() = default;

    // Generate tree (2 dir edges)
    static void generate_tree_edges(const Graph &graph, std::set<Edge> &edges) {
//        std::mt19937 gen;
//        edges.clear();
//
//        size_t last_taken = 0;
//
//        for (size_t i = 1; i < graph.size(); ++i) {
//            int32_t weight = 1;
//
//            size_t chosen = gen() % ++last_taken;
//            edges.emplace(chosen, last_taken, weight);
//        }
    }

    // Generate tree (2 dir edges)
    static void generate_tree_edges(const Graph &graph, std::vector<Edge> &edges) {
//        std::mt19937 gen;
//        edges.clear();
//
//        size_t last_taken = 0;
//
//        for (size_t i = 1; i < graph.size(); ++i) {
//            int32_t weight = 1;  // TODO
//
//            size_t chosen = gen() % ++last_taken;
//            edges.emplace_back(chosen, last_taken, weight);
//        }
    }

    static void build_by_prob(Graph &graph, float prob_f) {
        size_t prob = size_t(1 / prob_f);
        build_by_prob_inv(graph, prob);
    }

    // Generating random graph with one connected component
    // probability = 1 / prob
    static void build_by_prob_inv(Graph &graph, size_t prob) {
        std::mt19937 gen;
        std::set<Edge> edges;

        for (size_t i = 0; i < graph.size(); ++i) {
            for (size_t j = i + 1; j < graph.size(); ++j) {
//                std::pair<size_t, size_t> temp_edge(i, j);

                if (gen() % prob == 0) {
                    int32_t weight = 1;

                    edges.emplace(i, j, weight);
                }
            }
        }

        graph.build(edges);
    }

    // Generate graph with one connected component by num_of_edges
    static void build_by_number(Graph &graph, size_t num_of_edges) {
        std::mt19937 gen;
        std::set<Edge> edges;


        std::vector<std::pair<size_t , size_t >> total_edges_vector;
        for (size_t i = 0; i < graph.size(); ++i) {
            for (size_t j = 0; j < graph.size(); ++j) {
                total_edges_vector.emplace_back(i, j);
            }
        }
        std::shuffle(total_edges_vector.begin(), total_edges_vector.end(), gen);

        while (num_of_edges > 0 && !total_edges_vector.empty()) {
            std::pair<size_t, size_t> edge = total_edges_vector.back();
            int32_t weight = 1;

            edges.emplace(edge.first, edge.second, weight);
            total_edges_vector.pop_back();

            --num_of_edges;
        }


        graph.build(edges);
    }
};

#endif // MIN_CUT_GRAPH
