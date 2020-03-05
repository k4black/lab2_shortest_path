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


enum Colour {
    BLACK, GREY, WHITE
};


class Graph {
private:
    std::vector<std::vector<size_t>> neighbours;  // matrix
    size_t num_of_nodes = 0;

public:
    Graph() = default;

    explicit Graph(size_t num_of_nodes) {
        reserve(num_of_nodes);
    }

    void reserve(size_t num_of_nodes) {
        this->num_of_nodes = num_of_nodes;

        neighbours.clear();
        neighbours.resize(num_of_nodes);

        for (auto &line: neighbours) {
            l
        }
    }

    void build(const std::vector<std::pair<size_t, size_t>> &edges) {
        for (auto &edge : edges) {
            set_neighbour(edge.first, edge.second);
        }
    }

    void build(const std::set<std::pair<size_t, size_t>> &edges) {
        for (auto &edge : edges) {
            set_neighbour(edge.first, edge.second);
        }
    }

    void build_directed(const std::vector<std::pair<size_t, size_t>> &edges) {
        for (auto &edge : edges) {
            set_neighbour_directed(edge.first, edge.second);
        }
    }

    [[nodiscard]] const std::vector<size_t>& operator[](size_t node_index) const {
        return neighbours[node_index];
    }

    [[nodiscard]] const std::vector<size_t>&get_neighbours(size_t node_index) const {
        return neighbours[node_index];
    }

    void remove_neighbour(size_t node_index, size_t neighbour_index) {
        for(typename std::vector<size_t>::iterator it = neighbours[node_index].begin(); it != neighbours[node_index].end(); ++it) {
            if ((*it) == neighbour_index) {
                neighbours[node_index].erase(it);
                break;
            }
        }

        for(typename std::vector<size_t>::iterator it = neighbours[neighbour_index].begin(); it != neighbours[neighbour_index].end(); ++it) {
            if ((*it) == node_index) {
                neighbours[neighbour_index].erase(it);
                break;
            }
        }
    }

    void remove_neighbour_directed(size_t node_index, size_t neighbour_index) {
        for(typename std::vector<size_t>::iterator it = neighbours[node_index].begin(); it != neighbours[node_index].end(); ++it) {
            if ((*it) == neighbour_index) {
                neighbours[node_index].erase(it);
                break;
            }
        }
    }

    void set_neighbour(size_t node_index, size_t neighbour_index) {
        neighbours[node_index].push_back(neighbour_index);
        neighbours[neighbour_index].push_back(node_index);
    }

    void set_neighbour_directed(size_t node_index, size_t neighbour_index) {
        neighbours[node_index].push_back(neighbour_index);
    }

    void add_neighbour(size_t node_index, size_t neighbour_index) {
        if (!check_neighbour(node_index, neighbour_index)) {
            neighbours[node_index].push_back(neighbour_index);
        }
        if (!check_neighbour(neighbour_index, node_index)) {
            neighbours[neighbour_index].push_back(node_index);
        }
    }

    void add_neighbour_directed(size_t node_index, size_t neighbour_index) {
        if (!check_neighbour(node_index, neighbour_index)) {
            neighbours[node_index].push_back(neighbour_index);
        }
    }

    [[nodiscard]] bool check_neighbour(size_t node_index, size_t neighbour_index) const {
        for (size_t node_child : neighbours[node_index]) {
            if (node_child == neighbour_index) {
                return true;
            }
        }

        return false;
    }

    void revert_edge(size_t node_index, size_t neighbour_index) {
        if (!check_neighbour(node_index, neighbour_index)) {
            remove_neighbour_directed(node_index, neighbour_index);
            add_neighbour_directed(neighbour_index, node_index);
        }
    }

    [[nodiscard]] size_t node_degree(size_t node_index) const {
        return neighbours[node_index].size();
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
    static void generate_tree_edges(const Graph &graph, std::set<std::pair<size_t, size_t>> &edges) {
        std::mt19937 gen;
        edges.clear();

        size_t last_taken = 0;

        for (size_t i = 1; i < graph.size(); ++i) {
            size_t chosen = gen() % ++last_taken;
            edges.emplace(chosen, last_taken);
        }
    }

    // Generate tree (2 dir edges)
    static void generate_tree_edges(const Graph &graph, std::vector<std::pair<size_t, size_t>> &edges) {
        std::mt19937 gen;
        edges.clear();

        size_t last_taken = 0;

        for (size_t i = 1; i < graph.size(); ++i) {
            size_t chosen = gen() % ++last_taken;
            edges.emplace_back(chosen, last_taken);
        }
    }

    static void build_by_prob(Graph &graph, float prob_f) {
        size_t prob = size_t(1 / prob_f);
        build_by_prob_inv(graph, prob);
    }

    // Generating random graph with one connected component
    // probability = 1 / prob
    static void build_by_prob_inv(Graph &graph, size_t prob) {
        std::mt19937 gen;
        std::set<std::pair<size_t, size_t>> edges;

        // Generate tree (for one куомпонента связности)
        generate_tree_edges(graph, edges);

        // Generate others edges
        // by prob
        // TODO count prob accounting already taken edges
        for (size_t i = 0; i < graph.size(); ++i) {
            for (size_t j = i + 1; j < graph.size(); ++j) {
//                std::pair<size_t, size_t> temp_edge(i, j);

                if (gen() % prob == 0) {
                    edges.emplace(i, j);
                }
            }
        }

        graph.build(edges);
    }

    // Generate graph with one connected component by num_of_edges
    static void build_by_number(Graph &graph, size_t num_of_edges) {
        std::mt19937 gen;
        std::set<std::pair<size_t, size_t>> edges;

        // Generate tree (for one куомпонента связности)
        generate_tree_edges(graph, edges);

        std::vector<std::pair<size_t, size_t>> all_others_edges;
        for (size_t i = 0; i < graph.size(); ++i) {
            for (size_t j = i + 1; j < graph.size(); ++j) {
                std::pair<size_t, size_t> temp_edge(i, j);

                if (edges.find(temp_edge) == edges.end()) {
                    all_others_edges.push_back(temp_edge);
                }
            }
        }

        size_t left_edges = num_of_edges - (graph.size() - 1);
        std::shuffle(all_others_edges.begin(), all_others_edges.end(), gen);

        while (left_edges > 0 && !all_others_edges.empty()) {
            edges.insert(all_others_edges.back());
            all_others_edges.pop_back();
        }

        graph.build(edges);
    }
};

#endif // MIN_CUT_GRAPH
