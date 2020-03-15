#include <iostream>


#include "graph.cpp"
#include "shortest.cpp"


int main() {

    Graph graph(4);
    std::vector<Edge> edges = { { 0, 3, 0 }, { 2, 3, 1 },
                                { 1, 2, 3 }, { 0, 1, 5 } };
    graph.build(edges);

    for (size_t i : graph) {
        std::cout << i << ":    ";
        for (auto edge : graph[i]) {
            std::cout << edge.first << " $" << edge.second << "; ";
        }
        std::cout << '\n';
    }

    for (size_t i : graph) {
        for (size_t j : graph) {
            int32_t weight = graph.get_weight(i, j);
            if (weight == INT32_MAX) {
                std::cout << "- ";
            } else {
                std::cout << weight << " ";
            }
        }
        std::cout << '\n';
    }

    Graph new_graph(graph);
    new_graph.add_node();
    new_graph.set_neighbour_directed(3, 4, 4);
    graph.set_neighbour_directed(3, 3, 3);  // Test not influenced

    std::cout << "------------------\n";

    for (size_t i : new_graph) {
        std::cout << i << ":    ";
        for (auto edge : new_graph[i]) {
            std::cout << edge.first << " $" << edge.second << "; ";
        }
        std::cout << '\n';
    }

    for (size_t i : new_graph) {
        for (size_t j : new_graph) {
            int32_t weight = new_graph.get_weight(i, j);
            if (weight == INT32_MAX) {
                std::cout << "- ";
            } else {
                std::cout << weight << " ";
            }
        }
        std::cout << '\n';
    }



//    Graph graph(5);
//    std::vector<Edge> edges = { { 0, 1, -1 }, { 0, 2, 4 },
//                              { 1, 2, 3 }, { 1, 3, 2 },
//                              { 1, 4, 2 }, { 3, 2, 5 },
//                              { 3, 1, 1 }, { 4, 3, -3 } };
//    graph.build_directed(edges);
//
//    BellmanFord(graph, 0);



//    Graph graph(4);
//    std::vector<Edge> edges = { { 0, 3, 10 }, { 2, 3, 1 },
//                                { 1, 2, 3 }, { 0, 1, 5 } };
//    graph.build_directed(edges);
//
//    std::vector<std::vector<int64_t>> output;
//
//    FloydWarshall(graph, output);


    return 0;
}