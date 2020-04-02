#include <iostream>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <chrono>

#include "graph.cpp"
#include "shortest.cpp"

using namespace boost::numeric::ublas;


int main() {

//    boost::numeric::ublas::matrix<int32_t> test{};
//    std::cout << test.size1() <<std::endl;
//    std::cout << test.size2() <<std::endl;
//    return 0;

//    Graph graph(5);
//    std::vector<Edge> edges = { { 0, 1, -1 }, { 0, 2, 4 },
//                              { 1, 2, 3 }, { 1, 3, 2 },
//                              { 1, 4, 2 }, { 3, 2, 5 },
//                              { 3, 1, 1 }, { 4, 3, -3 } };
//    graph.build_directed(edges);
//
//    BellmanFord(graph, 0);

//    matrix<int32_t> test_matrix {3, 3};
//    test_matrix(1, 2) = 1;
//    test_matrix(2, 0) = 5;
//    std::cout << test_matrix << std::endl;
//    auto start = std::chrono::steady_clock::now();
//    auto test_squared = boost::numeric::ublas::element_prod(test_matrix, test_matrix);
//    auto end = std::chrono::steady_clock::now();
//    auto duration = end - start;
//    std::cout << std::chrono::duration<double, std::milli>(duration).count() << " ms" << std::endl;
//    auto test_squared = boost::numeric::ublas::prod(test_matrix, test_matrix);
//    std::cout << test_squared << std::endl;

//    return 0;

    Graph graph(4);
    std::vector<Edge> edges = { { 0, 3, 1 }, { 2, 3, 1 },
                                { 1, 2, 1 }, { 0, 1, 1 } };
    graph.build_directed(edges);
    std::vector<std::vector<int64_t>> dist;

    Johnson(graph, dist);

    std::vector<int32_t> heur{1, 1, 1, 1};
    std::vector<size_t> result{};

//    A_Star(graph, 0, 2, heur, result);

    Graph graph1{6};
    std::vector<Edge> edges1{{0, 1, 1},
                             {0, 2, 1},
                             {1, 2, 1},
                             {2, 4, 1},
                             {1, 3, 1},
                             {3, 5, 1}};
    graph1.build_unweighted(edges1);

    for (auto &elem : graph1.get_matrix()) {
        for (auto i : elem) std::cout << i << " ";
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::vector<std::vector<int64_t>> lengths{};
    std::vector<std::vector<size_t>> preds{};

    Seidel(graph1, lengths, preds, true);

    for (std::vector<int64_t> &i : lengths) {
        for (int64_t j : i) std::cout << j << " ";
        std::cout << std::endl;
    }

    std::cout << std::endl;

    for (std::vector<size_t> &i : preds) {
        for (int64_t j : i) std::cout << j << " ";
        std::cout << std::endl;
    }

    return 0;
}