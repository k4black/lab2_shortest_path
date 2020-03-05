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

    FloydWarshall(graph);


    return 0;

//
//    int graph[V][V] = { {0, 5, INF, 10},
//                        {INF, 0, 3, INF},
//                        {INF, INF, 0, 1},
//                        {INF, INF, INF, 0}
//    };
//
//    // Print the solution
//    FloydWarshall(graph);


//    Graph graph(10);
//
//    GraphGenerator::build_by_prob(graph, 0.3);
//
//    for (auto i : graph) {
////        std::cout << "=" << i << "\n";
//
//        for (auto weight : graph[i]) {
//            std::cout << weight << " ";
//        }
//
//
//        std::cout << "\n";
//    }
//
//
//    return 0;
}