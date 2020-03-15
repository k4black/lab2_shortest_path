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


namespace graph {
#include "../graph.cpp"
}


TEST_CASE("Testing Graph", "[graph]") {

    SECTION("Default_constructor") {
        graph::Graph graph;

        REQUIRE(graph.size() == 0);
    }

    SECTION("Constructor") {
        graph::Graph graph(10);

        REQUIRE(graph.size() == 10);
        REQUIRE(!graph.check_neighbour(0, 1));
        REQUIRE(!graph.check_neighbour(0, 5));
    }

    SECTION("Copy Constructor") {
        graph::Graph graph(4);
        std::vector<graph::Edge> edges = { { 0, 3, 0 }, { 2, 3, 1 },
                                    { 1, 2, 3 }, { 0, 1, 5 } };
        graph.build(edges);


        graph::Graph new_graph(graph);

        REQUIRE(new_graph.size() == 4);

        new_graph.add_node();

        REQUIRE(graph.size() == 4);
        REQUIRE(new_graph.size() == 5);

        new_graph.set_neighbour_directed(3, 4, 4);
        graph.set_neighbour_directed(3, 3, 3);  // Test not influenced

        REQUIRE(graph.get_weight(0, 3) == new_graph.get_weight(0, 3));
        REQUIRE(graph.get_weight(3, 3) != new_graph.get_weight(3, 3));  // Test not influenced
        REQUIRE(new_graph.get_weight(3, 4) == 4);
    }

    SECTION("Iterator") {
        graph::Graph graph(10);

        REQUIRE(*graph.begin() == 0);
        REQUIRE(*graph.end() == 10);

        size_t i = 0;
        for (graph::Graph::iterator it = graph.begin(); it != graph.end(); ++it) {
            REQUIRE((*it) == i);
            ++i;
        }
    }

    SECTION("Iterator for") {
        graph::Graph graph(15);

        size_t i = 0;
        for (size_t it : graph) {
            REQUIRE(it == i);
            ++i;
        }
    }

    SECTION("Build_from_edges") {
        graph::Graph graph(4);

        std::vector<graph::Edge> edges{{0, 1, 0},
                                       {1, 2, 0},
                                       {2, 3, 0},
                                       {3, 1, 0}};
        graph.build(edges);

        REQUIRE(graph.check_neighbour(0, 1));
        REQUIRE(graph.check_neighbour(1, 0));

        REQUIRE(graph.check_neighbour(1, 2));
        REQUIRE(graph.check_neighbour(2, 1));

        REQUIRE(!graph.check_neighbour(0, 3));
    }

    SECTION("Build_from_edges_directed") {
        graph::Graph graph(4);

        std::vector<graph::Edge> edges{{0, 1, 0},
                                       {1, 2, 0},
                                       {2, 3, 0},
                                       {3, 1, 0}};
        graph.build_directed(edges);

        REQUIRE(graph.check_neighbour(0, 1));
        REQUIRE(!graph.check_neighbour(1, 0));

        REQUIRE(graph.check_neighbour(1, 2));
        REQUIRE(!graph.check_neighbour(2, 1));

        REQUIRE(!graph.check_neighbour(0, 3));
    }

    SECTION("Test_get_neighbours") {
        graph::Graph graph(4);

        std::vector<graph::Edge> edges{{0, 1, 0},
                                       {1, 2, 0},
                                       {2, 3, 1},
                                       {3, 1, 5}};
        graph.build(edges);

        REQUIRE(graph.get_neighbours(0)[0].first == 1);
        REQUIRE(graph.get_neighbours(0)[0].second == 0);

        std::vector<size_t> neighbours;
        std::vector<int64_t> neighbours_weight;
        for (auto node : graph.get_neighbours(3)) {
            neighbours.push_back(node.first);
            neighbours_weight.push_back(node.second);
        }

        std::sort(neighbours.begin(), neighbours.end());
        REQUIRE((std::vector<size_t>{1, 2}) == neighbours);

        std::sort(neighbours_weight.begin(), neighbours_weight.end());
        REQUIRE((std::vector<int64_t>{1, 5}) == neighbours_weight);
    }

    SECTION("Test_get_neighbours_directed") {
        graph::Graph graph(4);

        std::vector<graph::Edge> edges{{0, 1, 0},
                                       {1, 2, 0},
                                       {2, 3, 1},
                                       {3, 1, 5}};
        graph.build_directed(edges);

        REQUIRE(graph.get_neighbours(0)[0].first == 1);
        REQUIRE(graph.get_neighbours(0)[0].second == 0);

        REQUIRE(graph.get_neighbours(3).size() == 1);
        REQUIRE(graph.get_neighbours(3)[0].first == 1);
        REQUIRE(graph.get_neighbours(3)[0].second == 5);
    }

    SECTION("Test_remove_neighbours") {
        graph::Graph graph(4);

        std::vector<graph::Edge> edges{{0, 1, 0},
                                       {1, 2, 0},
                                       {2, 3, 1},
                                       {3, 1, 5}};
        graph.build(edges);

        graph.remove_neighbour(1, 0);
        REQUIRE(graph.get_neighbours(0).empty());

        graph.remove_neighbour(2, 3);
        REQUIRE(graph.get_neighbours(3).size() == 1);
        REQUIRE(graph.get_neighbours(3)[0].first == 1);
        REQUIRE(graph.get_neighbours(3)[0].second == 5);
    }

    SECTION("Test_remove_neighbours_directed") {
        graph::Graph graph(4);

        std::vector<graph::Edge> edges{{0, 1, 0},
                                       {1, 2, 0},
                                       {2, 3, 1},
                                       {3, 1, 5}};
        graph.build_directed(edges);

        graph.remove_neighbour(1, 0);
        REQUIRE(graph.get_neighbours(0)[0].first == 1);
        REQUIRE(graph.get_neighbours(0)[0].second == 0);

        graph.remove_neighbour(2, 3);
        REQUIRE(graph.get_neighbours(2).empty());
        REQUIRE(graph.get_neighbours(3).size() == 1);
        REQUIRE(graph.get_neighbours(3)[0].first == 1);
        REQUIRE(graph.get_neighbours(3)[0].second == 5);
    }

}