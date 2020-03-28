#include <iostream>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <chrono>

#include "graph.cpp"
#include "shortest.cpp"

using namespace boost::numeric::ublas;


int main() {

//    Graph graph(5);
//    std::vector<Edge> edges = { { 0, 1, -1 }, { 0, 2, 4 },
//                              { 1, 2, 3 }, { 1, 3, 2 },
//                              { 1, 4, 2 }, { 3, 2, 5 },
//                              { 3, 1, 1 }, { 4, 3, -3 } };
//    graph.build_directed(edges);
//
//    BellmanFord(graph, 0);

    matrix<int32_t> test_matrix {3, 3};
    test_matrix(1, 2) = 1;
    test_matrix(2, 0) = 5;
    std::cout << test_matrix << std::endl;
//    auto start = std::chrono::steady_clock::now();
//    auto test_squared = boost::numeric::ublas::element_prod(test_matrix, test_matrix);
//    auto end = std::chrono::steady_clock::now();
//    auto duration = end - start;
//    std::cout << std::chrono::duration<double, std::milli>(duration).count() << " ms" << std::endl;
//    auto test_squared = boost::numeric::ublas::prod(test_matrix, test_matrix);
//    std::cout << test_squared << std::endl;

    return 0;

    Graph graph(4);
    std::vector<Edge> edges = { { 0, 3, 10 }, { 2, 3, 1 },
                                { 1, 2, 3 }, { 0, 1, 5 } };
    graph.build_directed(edges);
    std::vector<std::vector<int64_t>> dist;

    Johnson(graph, dist);

    int64_t result{};
    A_Star(graph, 0, 2, result, &simple_heuristic);


    return 0;
}