#include <iostream>


#include "graph.cpp"
#include "shortest.cpp"


int main() {

//    Graph graph(5);
//    std::vector<Edge> edges = { { 0, 1, -1 }, { 0, 2, 4 },
//                              { 1, 2, 3 }, { 1, 3, 2 },
//                              { 1, 4, 2 }, { 3, 2, 5 },
//                              { 3, 1, 1 }, { 4, 3, -3 } };
//    graph.build_directed(edges);
//
//    BellmanFord(graph, 0);



    Graph graph(4);
    std::vector<Edge> edges = { { 0, 3, 10 }, { 2, 3, 1 },
                                { 1, 2, 3 }, { 0, 1, 5 } };
    graph.build_directed(edges);

    std::vector<std::vector<int64_t>> output;

    Johnson(graph, output);

    int64_t result{};
    A_Star(graph, 0, 2, result, &simple_heuristic);


    return 0;
}