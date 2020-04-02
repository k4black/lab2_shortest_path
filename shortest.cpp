
#include "graph.cpp"

#include <set>
#include <algorithm>
#include <utility>
#include <boost/numeric/ublas/matrix.hpp>
#include <cmath>
#include <random>
#include <vector>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/operation_blocked.hpp>


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



void DijkstraSet(const Graph &graph, size_t src, std::vector<int64_t> &dist) {
    dist.clear();
    dist.resize(graph.size());
    std::fill(dist.begin(), dist.end(), INT64_MAX);
    dist[src] = 0;

    std::set< std::pair<int64_t, size_t>> setds;
    setds.emplace(0, src);


    // Find shortest path for all vertices
    while (!setds.empty()) {
        auto current_state = *(setds.begin());
        setds.erase(setds.begin());
        size_t current_node = current_state.second;
        int64_t current_dist = current_state.first;

        for (const auto &i : graph[current_node]){
            // Get vertex label and weight of current adjacent of u.
            size_t v = i.first;
            int64_t weight = i.second;

            //  If there is shorter path to v through u.
            if (current_dist != INT64_MAX && dist[v] > current_dist + weight && weight != INT32_MAX) {
                if (dist[v] != INT64_MAX)
                    setds.erase(setds.find(std::make_pair(dist[v], v)));

                // Updating distance of v
                dist[v] = current_dist + weight;
                setds.insert(std::make_pair(dist[v], v));
            }
        }
    }
}



int64_t BiDijkstra(const Graph &graph, size_t v, size_t w) {
    std::vector<bool> visited_v(graph.size(), false);
    std::vector<bool> visited_w(graph.size(), false);
    std::vector<int64_t> dist_v(graph.size(), INT64_MAX);
    std::vector<int64_t> dist_w(graph.size(), INT64_MAX);

    dist_v[v] = 0;
    dist_w[w] = 0;

    // Find shortest path for all vertices
    for (size_t count = 0; count < graph.size() - 1; count++) {
        size_t current_node_v = MinDistance(graph, dist_v, visited_v);
        size_t current_node_w = MinDistance(graph, dist_w, visited_w);

        if (dist_v[current_node_v] != INT64_MAX)
            for (auto node : graph[current_node_v]) {
                if (dist_v[current_node_v] != INT64_MAX && node.second != INT32_MAX && dist_v[current_node_v] + node.second < dist_v[node.first]) {
                    dist_v[node.first] = dist_v[current_node_v] + node.second;
                }
            }

        if (dist_w[current_node_w] != INT64_MAX)
            for (auto node : graph[current_node_w]) {
                if (dist_w[current_node_w] != INT64_MAX && node.second != INT32_MAX && dist_w[current_node_w] + node.second < dist_w[node.first]) {
                    dist_w[node.first] = dist_w[current_node_w] + node.second;
                }
            }

        if (visited_w[current_node_v] && visited_v[current_node_v]) {
            return dist_w[current_node_v] + dist_v[current_node_v];
        }
        if (visited_w[current_node_w] && visited_v[current_node_w]) {
            return dist_w[current_node_w] + dist_v[current_node_w];
        }

    }

    return dist_w[v];
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


bool Johnson(Graph &graph, std::vector<std::vector<int64_t>> &dist, bool useDijkstraSet=true) {
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
        if (useDijkstraSet) {
            DijkstraSet(graph_updated, node, dist[node]);
        } else {
            Dijkstra(graph_updated, node, dist[node]);
        }
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
        path.push_back(current);
    }
}


// TODO: modify for to-all paths
void A_Star(Graph &graph, size_t src, size_t dest, std::vector<int32_t> &heur, int64_t &length, std::vector<size_t> &output, bool reconstruct = false) {
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
    f_score[src] = heur[src];

    while (!open_set.empty()) {
        std::pop_heap(open_set.begin(), open_set.end());
        size_t current = open_set[open_set.size()-1].second;
        open_set.resize(open_set.size()-1);

//        std::cout << "current vertex is " << current << std::endl;

        // TODO: maybe if this condition is removed, then we get to-all variant?
        if (current == dest) {
            // TODO: reconstruct path and its price
//            std::cout << "Algorithm finished" << std::endl;
            if (reconstruct) reconstruct_path(came_from, dest, output);
//            for (size_t &vertex : path) {
//                std::cout << vertex << " ";
//            }
//            std::cout << std::endl;
            length = g_score[dest];
        }

        for (auto &neib : graph[current]) {
//            std::cout << neib.first << " " << neib.second << std::endl;
//            std::cout << "    looking at neighbour " << neib.first << std::endl;
            auto relax_score = g_score[current] + neib.second;
            if (relax_score < g_score[neib.first]) {
                came_from[neib.first] = current;
                g_score[neib.first] = relax_score;
                f_score[neib.first] = g_score[neib.first] + heur[neib.first];
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

    boost::numeric::ublas::matrix<int32_t> Z{A.size1(), A.size1()};
//    boost::numeric::ublas::axpy_prod(A, A, Z);
    Z = boost::numeric::ublas::block_prod<boost::numeric::ublas::matrix<int32_t>, 64>(A, A);
//    boost::numeric::ublas::prod(A, A);  //TODO: dedicated function for matrix exponentiation in uBLAS/in general?
    boost::numeric::ublas::matrix<int64_t> B {A.size1(), A.size1()};  //TODO: dedicated class for square matrices in uBLAS?
    for (size_t i = 0; i < A.size1(); ++i) {
        for (size_t j = 0; j < A.size1(); ++j) {
            if ((A(i, j) == 1 || Z(i, j) > 0) && i != j) {B(i, j) = 1;} else {B(i, j) = 0;}
        }
    }

    bool ret = true;
    for (size_t i = 0; i < A.size1(); ++i) {
        for (size_t j = 0; j < A.size1(); ++j) {
            if (B(i, j) == 0 && i != j) {
                ret = false;
                break;
            }
        }
        if (!ret) {break;}
    }
    if (ret) {return 2*B-A;}

    auto T = Seidel_internal(B);
    std::vector<int64_t> degree{};
    degree.reserve(A.size1());
    for (size_t i = 0; i < A.size1(); ++i) {
        //TODO: better sum of all elements in a row?
        for (size_t j = 0; j < A.size1(); ++j) {
            degree[i] += A(i, j);
        }
    }
    boost::numeric::ublas::matrix<int32_t> X {A.size1(), A.size1()};
//    auto X = boost::numeric::ublas::prod(T, A);
//    boost::numeric::ublas::axpy_prod(T, A, X);
    X = boost::numeric::ublas::block_prod<boost::numeric::ublas::matrix<int32_t>, 64>(T, A);
    boost::numeric::ublas::matrix<int32_t> D {A.size1(), A.size1()};
    for (size_t i = 0; i < A.size1(); ++i) {
        for (size_t j = 0; j < A.size1(); ++j) {
            if (X(i, j) >= T(i, j)*degree[j]) {D(i, j) = 2*T(i, j);} else {D(i, j) = 2*T(i, j)-1;}
        }
    }
    return D;
}

boost::numeric::ublas::matrix<int32_t> BPWM(boost::numeric::ublas::matrix<int64_t> &A, boost::numeric::ublas::matrix<int32_t> &Dr) {
//    std::cout << "A:" << A << std::endl;
//    std::cout << "Dr:" << Dr << std::endl;
    boost::numeric::ublas::matrix<int64_t> W{A.size1(), A.size2()};
    boost::numeric::ublas::axpy_prod(-A, Dr, W);
//    boost::numeric::ublas::matrix<int64_t> W = -boost::numeric::ublas::prod(A, Dr);
    for (size_t l = 0; l <= std::ceil(std::log2(W.size1())) - 1; ++l) {
        int64_t d = std::pow(2, l);
        for (size_t i = 0; i < std::ceil(3.42 * std::log2(W.size1())); ++i) {
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_int_distribution<size_t> dis(0, A.size1() - 1);
            std::vector<size_t> k{};
            for (size_t j = 0; j < d; ++j) {
                k.push_back(dis(gen));
            }
            boost::numeric::ublas::matrix<int32_t> X{A.size1(), (size_t) d};
            for (size_t n = 0; n < d; ++n) {
                for (size_t m = 0; m < A.size1(); ++m) {
                    X(m, n) = k[n] * A(m, k[n]);
                }
            }
            boost::numeric::ublas::matrix<int32_t> Y{(size_t) d, A.size1()};
            for (size_t m = 0; m < d; ++m) {
                for (size_t n = 0; n < A.size1(); ++n) {
                    Y(m, n) = Dr(k[m], n);
                }
            }
            boost::numeric::ublas::matrix<int64_t> C{A.size1(), A.size2()};
//            boost::numeric::ublas::prod(X, Y);
            boost::numeric::ublas::axpy_prod(X, Y, C);
//            std::cout << "C:" << C << std::endl;
            for (size_t m = 0; m < W.size1(); ++m) {
                for (size_t n = 0; n < W.size2(); ++n) {
//                    if (C(m, n) < 0 || C(m, n) >= A.size2() || C(m, n) >= Dr.size1() || m >= A.size1() || n >= Dr.size2()) {
//                        std::cout << "Will fail, C(m, n) = " << C(m, n) << std::endl;
//                        return boost::numeric::ublas::matrix<int32_t>{};
//                    }
//                    std::cout << C(m, n) << " ";
//                    std::cout << W(m, n) << " ";
//                    std::cout << A(m, C(m, n)) << " ";
//                    std::cout << Dr(C(m, n), n) << " ";
//                    std::cout << std::endl;
                    if ((W(m, n) < 0) && (A(m, C(m, n)) == 1) && (Dr(C(m, n), n) == 1)) {
                        W(m, n) = C(m, n);
                    }
                }
            }
        }
    }
    for (size_t m = 0; m < W.size1(); ++m) {
        for (size_t n = 0; n < W.size2(); ++n) {
            if (W(m, n) < 0) {
                for (size_t p = 0; p < W.size1(); ++p) {
                    if ((A(m, p) == 1) && (Dr(p, n) == 1)) {
                        W(m, n) = p;
                        break;
                    }
                }
            }
        }
    }
    return W;
}

//TODO: std::functional to pass custom matrix multiplication functions
//TODO: int32_t -> {0; 1}, as naive Seidel is for unweighted graph only
void Seidel(const Graph &graph, std::vector<std::vector<int64_t>> &lengths, std::vector<std::vector<size_t>> &preds, bool reconstruct = false) {
    lengths.clear();
    for (size_t i = 0; i < graph.size(); ++i) {
        lengths.emplace_back(graph.size());
    }
    preds.clear();
    for (size_t i = 0; i < graph.size(); ++i) {
        preds.emplace_back(graph.size());
    }
    boost::numeric::ublas::matrix<int64_t> A {graph.size(), graph.size()};
    for (size_t i = 0; i < graph.size(); ++i) {
        for (size_t j = 0; j < graph.size(); ++j) {
            A(i, j) = graph.get_matrix()[i][j];
        }
    }
    auto D = Seidel_internal(A);
//    std::cout << D << std::endl;
    for (size_t i = 0; i < A.size1(); ++i) {
        for (size_t j = 0; j < A.size1(); ++j) {
            lengths[i][j] = D(i, j);
        }
    }

    if (reconstruct) {
        std::vector<boost::numeric::ublas::matrix<int32_t>> Wr{};
        Wr.emplace_back(graph.size(), graph.size());
        Wr.emplace_back(graph.size(), graph.size());
        Wr.emplace_back(graph.size(), graph.size());
        for (int r = 0; r <= 2; ++r) {
            boost::numeric::ublas::matrix<int32_t> Dr{graph.size(), graph.size()};
            for (size_t i = 0; i < graph.size(); ++i) {
                for (size_t j = 0; j < graph.size(); ++j) {
                    if ((D(i, j) + 1) % 3 == r) { Dr(i, j) = 1; } else { Dr(i, j) = 0; }
                }
            }
            Wr[r] = BPWM(A, Dr);
//        if (Wr[r].size1() == 0 && Wr[r].size2() == 0) {
//            for (size_t i = 0; i < graph.size(); ++i) {
//                for (size_t j = 0; j < graph.size(); ++j) {
//                    preds[i][j] = -1;
//                }
//            }
//            return;
//        }
        }
        for (size_t i = 0; i < graph.size(); ++i) {
            for (size_t j = 0; j < graph.size(); ++j) {
                preds[i][j] = Wr[D(i, j) % 3](i, j);
            }
        }
    }
}
