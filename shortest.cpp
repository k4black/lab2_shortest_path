
#include "graph.cpp"

#include <set>
#include <algorithm>
#include <utility>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>


int64_t MinDistance(const Graph& graph, std::vector<int64_t> &dist, std::vector<bool> &sptSet)  {
    /// Auxiliary function for easily getting
    /// the node for the current iteration
    /// in the Dijkstra algorithm
    int64_t min = INT64_MAX;  //TODO: maybe some other stuff, less explicit?
    size_t min_index = 0;

    for (size_t node : graph) {
        if (sptSet[node] == false && dist[node] <= min) {
            min = dist[node];
            min_index = node;
        }
    }

    return min_index;
}


void Dijkstra(const Graph& graph, size_t src, std::vector<int64_t> &dist) {
    /// Solves single-source shortest-path problem for graphs with nonnegative edges
    // prerequisites
    dist.clear();
    dist.resize(graph.size());
    std::fill(dist.begin(), dist.end(), INT64_MAX);
    std::vector<bool> sptSet(graph.size(), false);  // TODO: Heap

    // initializing the starting vertex
    dist[src] = 0;

    // iterating over all vertices
    for (size_t count = 0; count < graph.size() - 1; count++) {
        // initializing vertex for the current iteration
        size_t current_node = MinDistance(graph, dist, sptSet);

        // marking vertex as visited (each vertex is visited one time)
        sptSet[current_node] = true;

        // iterating over neighbourhood of the current vertex
        for (auto neib : graph[current_node]) {
            // change vertex distance mark according to the algorithm's condition
            if (dist[current_node] != INT64_MAX && dist[current_node] + neib.second < dist[neib.first]) {
                dist[neib.first] = dist[current_node] + neib.second;
            }
        }
    }
}


bool BellmanFord(Graph &graph, size_t src, std::vector<int64_t> &dist) {
    /// Solves single-source shortest-path problem for graphs with nonnegative cycles
    /// Returns true if can be solved, false otherwise
    dist.clear();
    dist.resize(graph.size());
    std::fill(dist.begin(), dist.end(), INT64_MAX);
    dist[src] = 0;

    // getting graph as list of edges
    std::vector<Edge> edges;
    for (size_t node : graph) {
        for (auto &neib : graph[node]) {
            edges.emplace_back(node, neib.first, neib.second);
        }
    }

    // relaxation of all edges, N-1 times
    for (size_t i = 0; i < graph.size() - 1; ++i) {
        for (auto &edge : edges) {
            if (dist[edge.first] != INT64_MAX && edge.weight != INT32_MAX && dist[edge.first] + edge.weight < dist[edge.second]) {
                dist[edge.second] = dist[edge.first] + edge.weight;
            }
        }
    }

    // check for negative-weight cycles
    // algorithmically, checking whether next relaxation would change any distances
    for (auto &edge : edges) {
        size_t x = edge.first;
        size_t y = edge.second;
        int32_t weight = edge.weight;

        if (dist[x] != INT64_MAX && weight != INT32_MAX && dist[x] + weight < dist[y]) {
            return false;  // TODO: Raise error
        }
    }

    return true;
}


void FloydWarshall(Graph &graph, std::vector<std::vector<int64_t>> &dist) {
    /// Solves all-pairs shortest-path problem
    dist.clear();
    dist.resize(graph.size());
    for (size_t node : graph) {
        dist[node].resize(graph.size());
        std::fill(dist[node].begin(), dist[node].end(), INT64_MAX);

        for (auto &neib : graph[node]) {
            dist[node][neib.first] = neib.second;
        }

        dist[node][node] = 0;
    }

    for (size_t k : graph) {
        for (size_t i : graph) {
            for (size_t j : graph) {
                if (dist[i][k] == INT64_MAX || dist[k][j] == INT64_MAX) {
                    continue;
                }

                if (dist[i][k] + dist[k][j] < dist[i][j]) {
                    dist[i][j] = dist[i][k] + dist[k][j];
                }
            }
        }
    }
}


bool Johnson(Graph &graph, std::vector<std::vector<int64_t>> &dist) {
    /// Solves all-pairs shortest-path problem
    dist.clear();
    dist.resize(graph.size());

    // Running Bellman-Ford algorithm on the augmented graph
    Graph graph_prime(graph);
    graph_prime.add_node();  // Augmenting node
    // Augmenting zero-weight edges
    for (size_t node = 0; node < graph_prime.size() - 1; ++node) {
        graph_prime.set_neighbour_directed(graph.size(), node, 0);
    }

    std::vector<int64_t> belman_dist;

    if (!BellmanFord(graph_prime, graph.size(), belman_dist)) {
        std::cout << "Graph contains negative weight cycle" << std::endl;  // TODO: raise error
        return false;
    }

    // TODO: Optimize to copy constructor
    // Updating weights of the initial graph with Bellman-Ford results
    Graph graph_updated(graph.size());
    std::vector<Edge> edges_updated;
    for (size_t node : graph) {
        for (auto &neib : graph[node]) {
            if (neib.second != 0) {
                edges_updated.emplace_back(node, neib.first, neib.second + belman_dist[node] - belman_dist[neib.first]);
            } else {
                edges_updated.emplace_back(node, neib.first, 0);
            }
        }
    }

    graph_updated.build_directed(edges_updated);

    // Running Dijkstra on the updated graph
    for (size_t node : graph) {
        Dijkstra(graph_updated, node, dist[node]);
    }

    return true;
}


int64_t simple_heuristic(size_t vertex) {
    return 1;
}


void reconstruct_path(std::unordered_map<size_t, size_t> &came_from, size_t dest, std::vector<size_t> &path) {
    path.clear();
//    path.resize(came_from.size());
    path.push_back(dest);
    size_t current = dest;

    while (came_from.find(current) != came_from.end() && came_from[current] != current) {
        current = came_from[current];
        std::cout << "    " << current << std::endl;
        path.push_back(current);
    }
}


// TODO: modify for to-all paths
void A_Star(Graph &graph, size_t src, size_t dest, int64_t &result, std::function<int64_t (size_t)> heuristics = simple_heuristic) {
    // heap containing the nodes to be visited
    std::vector<std::pair<int32_t , size_t >> open_set {{0, src}};
    std::make_heap(open_set.begin(), open_set.end(), std::greater<>{});

    // TODO: for reconstructing the path
//    std::vector<size_t> came_from {};
    std::unordered_map<size_t, size_t> came_from{{src, src}};

    // g_score[n] is the cost of the cheapest path from start to n currently known
    std::unordered_map<size_t, int64_t> g_score {};
    for (size_t i = 0; i < graph.size(); ++i) {
        g_score[i] = INT64_MAX;
    }
    g_score[src] = 0;

    std::unordered_map<size_t, int64_t> f_score {};
    for (size_t i = 0; i < graph.size(); ++i) {
        f_score[i] = INT64_MAX;
    }
    f_score[src] = heuristics(src);

    while (!open_set.empty()) {
        std::pop_heap(open_set.begin(), open_set.end());
        size_t current = open_set[open_set.size()-1].second;
        open_set.resize(open_set.size()-1);

        std::cout << "current vertex is " << current << std::endl;

        // TODO: maybe if this condition is removed, then we get to-all variant?
        if (current == dest) {
            // TODO: reconstruct path and its price
            std::cout << "Algorithm finished" << std::endl;
            std::vector<size_t> path;
            reconstruct_path(came_from, dest, path);
            for (size_t &vertex : path) {
                std::cout << vertex << " ";
            }
            std::cout << std::endl;
            return;
        }

        for (auto &neib : graph[current]) {
            std::cout << neib.first << " " << neib.second << std::endl;
            std::cout << "    looking at neighbour " << neib.first << std::endl;
            auto relax_score = g_score[current] + neib.second;
            if (relax_score < g_score[neib.first]) {
                came_from[neib.first] = current;
                g_score[neib.first] = relax_score;
                f_score[neib.first] = g_score[neib.first] + heuristics(neib.first);
                // if neib not in open_set
                if (std::find_if(open_set.begin(), open_set.end(), [&neib](std::pair<int32_t, size_t> elem){return elem.second == neib.first;}) == open_set.end()) {
                    open_set.emplace_back(neib.second, neib.first);  // TODO: this way?
                    std::make_heap(open_set.begin(), open_set.end(), std::greater<>{});
                }
            }
        }
    }
}

//TODO: std::functional to pass custom matrix multiplication functions
//TODO: int32_t -> {0; 1}, as naive Seidel is for unweighted graph only
boost::numeric::ublas::matrix<int64_t> Seidel_internal(boost::numeric::ublas::matrix<int64_t> &A) {

    bool ret = true;
    for (size_t i = 0; i < A.size1(); ++i) {
        for (size_t j = 0; j < A.size1(); ++j) {
            if (A(i, j) == 0 && i != j) {
                ret = false;
                break;
            }
        }
        if (!ret) {break;}
    }
    if (ret) {return A;}

    boost::numeric::ublas::matrix<int32_t> Z = boost::numeric::ublas::prod(A, A);  //TODO: dedicated function for matrix exponentiation in uBLAS/in general?
    boost::numeric::ublas::matrix<int64_t> B {A.size1(), A.size1()};  //TODO: dedicated class for square matrices in uBLAS?
    for (size_t i = 0; i < A.size1(); ++i) {
        for (size_t j = 0; j < A.size1(); ++j) {
            if ((A(i, j) == 1 || Z(i, j) > 0) && i != j) {B(i, j) = 1;} else {B(i, j) = 0;}
        }
    }
    auto T = Seidel_internal(B);
    std::vector<int64_t> degree{};
    degree.reserve(A.size1());
    for (size_t i = 0; i < A.size1(); ++i) {
        //TODO: better sum of all elements in a row?
        for (size_t j = 0; j < A.size1(); ++j) {
            degree[i] += A(i, j);
        }
    }
    auto X = boost::numeric::ublas::prod(T, A);
    boost::numeric::ublas::matrix<int32_t> D {A.size1(), A.size1()};
    for (size_t i = 0; i < A.size1(); ++i) {
        for (size_t j = 0; j < A.size1(); ++j) {
            if (X(i, j) >= T(i, j)*degree[j]) {D(i, j) = 2*T(i, j);} else {D(i, j) = 2*T(i, j)-1;}
        }
    }
    return D;
}

//TODO: std::functional to pass custom matrix multiplication functions
//TODO: int32_t -> {0; 1}, as naive Seidel is for unweighted graph only
void Seidel(const Graph &graph, std::vector<std::vector<int64_t>> &output) {
    boost::numeric::ublas::matrix<int64_t> A {graph.size(), graph.size()};
    for (size_t i = 0; i < graph.size(); ++i) {
        for (size_t j = 0; j < graph.size(); ++j) {
            A(i, j) = graph.get_matrix()[i][j];
        }
    }
    auto ret = Seidel_internal(A);
    for (size_t i = 0; i < A.size1(); ++i) {
        for (size_t j = 0; j < A.size1(); ++j) {
            output[i][j] = ret(i, j);
        }
    }
}
