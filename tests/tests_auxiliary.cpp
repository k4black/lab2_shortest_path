//#pragma clang diagnostic push
//#pragma clang diagnostic ignored "-Wunknown-pragmas"
#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch2.hpp"

#include <vector>
#include <utility>
#include <queue>
#include <algorithm>
#include <random>
#include <unordered_map>
#include <chrono>
#include <iostream>
#include <stack>
#include <chrono>
#include <memory>
#include <memory>
#include <map>
#include <limits>
#include <stack>
#include <iostream>
#include <queue>

namespace task {
    #include "../shortest.cpp"
}



TEST_CASE("Testing Dijkstra", "[dijkstra]") {

    SECTION("Simple case") {
        task::Graph graph(3);
        std::vector<task::Edge> edges{{0, 1, 1},
                                      {1, 2, 5},
                                      {0, 2, 100}};
        graph.build_directed(edges);

        std::vector<int64_t> shortest;
        task::Dijkstra(graph, 0, shortest);

        std::vector<int64_t> shortest_true{0, 1, 6};
        REQUIRE(shortest == shortest_true);
    }

    SECTION("Empty case") {
        task::Graph graph(100);
        std::vector<task::Edge> edges{};
        graph.build_directed(edges);

        std::vector<int64_t> shortest;
        task::Dijkstra(graph, 0, shortest);

        std::vector<int64_t> shortest_true(100, INT64_MAX);
        shortest_true[0] = 0;
        REQUIRE(shortest == shortest_true);
    }

    SECTION("One edge") {
        task::Graph graph(100);
        std::vector<task::Edge> edges{{10, 1, 0}};
        graph.build_directed(edges);

        std::vector<int64_t> shortest;
        task::Dijkstra(graph, 10, shortest);

        std::vector<int64_t> shortest_true(100, INT64_MAX);
        shortest_true[10] = 0;
        shortest_true[1] = 0;
        REQUIRE(shortest == shortest_true);
    }

    SECTION("Medium-Simple case") {
        task::Graph graph(6);
        std::vector<task::Edge> edges{{0, 1, 10},
                                      {0, 4, 3},
                                      {1, 4, 3},
                                      {1, 2, 4},
                                      {4, 2, 7},
                                      {2, 3, 1},
                                      {5, 3, 3},
                                      {5, 0, 3}};
        graph.build_directed(edges);

        std::vector<int64_t> shortest;
        task::Dijkstra(graph, 0, shortest);

        std::vector<int64_t> shortest_true{0, 10, 10, 11, 3, INT64_MAX};
        REQUIRE(shortest == shortest_true);
    }

    SECTION("Medium case") {
        task::Graph graph(9);
        // http://www.geeksforgeeks.org/wp-content/uploads/gat2012.png
        std::vector<task::Edge> edges{{0, 2, 1},
                                      {8, 0, 4},
                                      {2, 4, 1},
                                      {2, 3, 3},
                                      {4, 6, 2},
                                      {6, 4, 2},
                                      {4, 7, 4},
                                      {6, 7, 3},
                                      {3, 4, 1},
                                      {3, 7, 3},
                                      {3, 5, 5},
                                      {7, 5, 5},
                                      {8, 3, 7},
                                      {1, 3, 4},
                                      {8, 1, 3},
                                      {1, 8, 3}};
        graph.build_directed(edges);

        std::vector<int64_t> shortest;
        task::Dijkstra(graph, 0, shortest);

        std::vector<int64_t> shortest_true{0, INT64_MAX, 1, 4, 2, 9, 4, 6, INT64_MAX};
        REQUIRE(shortest == shortest_true);
    }
}


TEST_CASE("Testing FloydWarshall", "[floyd_warshall]") {

    SECTION("Simple case") {
        task::Graph graph(3);
        std::vector<task::Edge> edges{{0, 1, 1},
                                      {1, 2, 5},
                                      {0, 2, 100}};
        graph.build_directed(edges);

        std::vector<std::vector<int64_t>> shortest;
        task::FloydWarshall(graph, shortest);

        std::vector<std::vector<int64_t>> shortest_true{{0,         1,         6},
                                                        {INT64_MAX, 0,         5},
                                                        {INT64_MAX, INT64_MAX, 0}};
        REQUIRE(shortest == shortest_true);
    }

    SECTION("Medium-Simple case") {
        task::Graph graph(5);
        std::vector<task::Edge> edges{{0, 1, 10},
                                      {0, 4, 3},
                                      {1, 4, 3},
                                      {1, 2, 4},
                                      {4, 2, 7},
                                      {2, 3, 1}};
        graph.build_directed(edges);

        std::vector<std::vector<int64_t>> shortest;
        task::FloydWarshall(graph, shortest);

        for (auto node : graph) {
            std::vector<int64_t> shortest_true;

            task::Dijkstra(graph, node, shortest_true);
            REQUIRE(shortest[node] == shortest_true);
        }
    }

    SECTION("Medium case") {
        task::Graph graph(9);
        // http://www.geeksforgeeks.org/wp-content/uploads/gat2012.png
        std::vector<task::Edge> edges{{0, 2, 1},
                                      {8, 0, 4},
                                      {2, 4, 1},
                                      {2, 3, 3},
                                      {4, 6, 2},
                                      {6, 4, 2},
                                      {4, 7, 4},
                                      {6, 7, 3},
                                      {3, 4, 1},
                                      {3, 7, 3},
                                      {3, 5, 5},
                                      {7, 5, 5},
                                      {8, 3, 7},
                                      {1, 3, 4},
                                      {8, 1, 3},
                                      {1, 8, 3}};
        graph.build_directed(edges);

        std::vector<std::vector<int64_t>> shortest;
        task::FloydWarshall(graph, shortest);

        for (auto node : graph) {
            std::vector<int64_t> shortest_true;

            task::Dijkstra(graph, node, shortest_true);
            REQUIRE(shortest[node] == shortest_true);
        }
    }

    SECTION("Random case") {
        size_t graph_size = 100;
        size_t num_of_rep = 10;
        float prob = 0.2;
        task::Graph graph(graph_size);

        for (size_t i = 0; i < num_of_rep; ++i) {
            bool build_correct = false;
            while (!build_correct) {
                task::GraphGenerator::build_by_prob(graph, prob);

                std::vector<int64_t> shortest;
                bool negative = task::BellmanFord(graph, 0, shortest);

                if (negative == false) {
                    build_correct = false;
                } else {
                    build_correct = true;
                }
            }

            std::vector<std::vector<int64_t>> shortest;
            task::FloydWarshall(graph, shortest);

            for (auto node : graph) {
                std::vector<int64_t> shortest_true;

                task::Dijkstra(graph, node, shortest_true);
                REQUIRE(shortest[node] == shortest_true);
            }
        }
    }

}




TEST_CASE("Testing BellmanFord", "[bellman_ford]") {


    SECTION("Simple case") {
        task::Graph graph(3);
        std::vector<task::Edge> edges{{0, 1, 1},
                                      {1, 2, 5},
                                      {0, 2, 100}};
        graph.build_directed(edges);

        std::vector<int64_t> shortest;
        task::BellmanFord(graph, 0, shortest);

        std::vector<int64_t> shortest_true{0, 1, 6};
        REQUIRE(shortest == shortest_true);
    }


    SECTION("Simple Negative case") {
        task::Graph graph(3);
        std::vector<task::Edge> edges{{0, 1, -3},
                                      {1, 2, 1},
                                      {2, 0, 1}};
        graph.build_directed(edges);

        std::vector<int64_t> shortest;
        bool negative = task::BellmanFord(graph, 0, shortest);

        REQUIRE(negative == false);
    }


    SECTION("Medium Negative case") {
        task::Graph graph(5);
        std::vector<task::Edge> edges{{0, 1, 3},
                                      {1, 4, 1},
                                      {4, 0, 7},
                                      {4, 3, 2},
                                      {3, 2, 1},
                                      {2, 1, -5},
                                      {1, 3, 10}};
        graph.build_directed(edges);

        std::vector<int64_t> shortest;
        bool negative = task::BellmanFord(graph, 0, shortest);

        REQUIRE(negative == false);
    }

    SECTION("Empty case") {
        task::Graph graph(100);
        std::vector<task::Edge> edges{};
        graph.build_directed(edges);

        std::vector<int64_t> shortest;
        task::BellmanFord(graph, 0, shortest);

        std::vector<int64_t> shortest_true(100, INT64_MAX);
        shortest_true[0] = 0;
        REQUIRE(shortest == shortest_true);
    }

    SECTION("One edge") {
        task::Graph graph(100);
        std::vector<task::Edge> edges{{10, 1, 0}};
        graph.build_directed(edges);

        std::vector<int64_t> shortest;
        task::BellmanFord(graph, 10, shortest);

        std::vector<int64_t> shortest_true(100, INT64_MAX);
        shortest_true[10] = 0;
        shortest_true[1] = 0;
        REQUIRE(shortest == shortest_true);
    }

    SECTION("Medium-Simple case") {
        task::Graph graph(6);
        std::vector<task::Edge> edges{{0, 1, 10},
                                      {0, 4, 3},
                                      {1, 4, 3},
                                      {1, 2, 4},
                                      {4, 2, 7},
                                      {2, 3, 1},
                                      {5, 3, 3},
                                      {5, 0, 3}};
        graph.build_directed(edges);

        std::vector<int64_t> shortest;
        task::BellmanFord(graph, 0, shortest);

        std::vector<int64_t> shortest_true{0, 10, 10, 11, 3, INT64_MAX};
        REQUIRE(shortest == shortest_true);
    }

    SECTION("Medium case") {
        task::Graph graph(9);
        // http://www.geeksforgeeks.org/wp-content/uploads/gat2012.png
        std::vector<task::Edge> edges{{0, 2, 1},
                                      {8, 0, 4},
                                      {2, 4, 1},
                                      {2, 3, 3},
                                      {4, 6, 2},
                                      {6, 4, 2},
                                      {4, 7, 4},
                                      {6, 7, 3},
                                      {3, 4, 1},
                                      {3, 7, 3},
                                      {3, 5, 5},
                                      {7, 5, 5},
                                      {8, 3, 7},
                                      {1, 3, 4},
                                      {8, 1, 3},
                                      {1, 8, 3}};
        graph.build_directed(edges);

        std::vector<int64_t> shortest;
        task::BellmanFord(graph, 0, shortest);

        std::vector<int64_t> shortest_true{0, INT64_MAX, 1, 4, 2, 9, 4, 6, INT64_MAX};
        REQUIRE(shortest == shortest_true);
    }

}





TEST_CASE("Testing Johnson", "[johnson]") {

    SECTION("Simple case") {
        task::Graph graph(3);
        std::vector<task::Edge> edges{{0, 1, 1},
                                      {1, 2, 5},
                                      {0, 2, 100}};
        graph.build_directed(edges);

        std::vector<std::vector<int64_t>> shortest;
        task::Johnson(graph, shortest);

        std::vector<std::vector<int64_t>> shortest_true{{0,         1,         6},
                                                        {INT64_MAX, 0,         5},
                                                        {INT64_MAX, INT64_MAX, 0}};
        REQUIRE(shortest == shortest_true);
    }

    SECTION("Medium-Simple case") {
        task::Graph graph(5);
        std::vector<task::Edge> edges{{0, 1, 10},
                                      {0, 4, 3},
                                      {1, 4, 3},
                                      {1, 2, 4},
                                      {4, 2, 7},
                                      {2, 3, 1}};
        graph.build_directed(edges);

        std::vector<std::vector<int64_t>> shortest;
        task::Johnson(graph, shortest);

        for (auto node : graph) {
            std::vector<int64_t> shortest_true;

            task::Dijkstra(graph, node, shortest_true);
            REQUIRE(shortest[node] == shortest_true);
        }
    }

    SECTION("Random case") {
        size_t graph_size = 100;
        size_t num_of_rep = 10;
        float prob = 0.2;
        task::Graph graph(graph_size);

        for (size_t i = 0; i < num_of_rep; ++i) {
            bool build_correct = false;
            while (!build_correct) {
                task::GraphGenerator::build_by_prob(graph, prob);

                std::vector<int64_t> shortest;
                bool negative = task::BellmanFord(graph, 0, shortest);

                if (negative == false) {
                    build_correct = false;
                } else {
                    build_correct = true;
                }
            }

            std::vector<std::vector<int64_t>> shortest;
            task::Johnson(graph, shortest);

            std::vector<std::vector<int64_t>> shortest_true;
            task::FloydWarshall(graph, shortest_true);

            REQUIRE(shortest == shortest_true);
        }
    }

}