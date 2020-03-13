
#include "graph.cpp"

#include <set>
#include <algorithm>
#include <utility>


int64_t MinDistance(const Graph& graph, std::vector<int64_t> &dist, std::vector<bool> &sptSet)  {
    // Initialize min value
    int64_t min = INT64_MAX;
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
    dist.clear();
    dist.resize(graph.size());
    std::fill(dist.begin(), dist.end(), INT64_MAX);
    dist[src] = 0;

    std::vector<bool> sptSet(graph.size(), false);  // TODO: Heap


    // Find shortest path for all vertices
    for (size_t count = 0; count < graph.size() - 1; count++) {
        size_t current_node = MinDistance(graph, dist, sptSet);

        sptSet[current_node] = true;

        for (auto neib : graph[current_node]) {
            if (dist[current_node] != INT64_MAX && dist[current_node] + neib.second < dist[neib.first]) {
//                std::cout << dist[neib.first] << " " << dist[current_node] << "+" << neib.second << '\n';
                dist[neib.first] = dist[current_node] + neib.second;
            }
        }
    }
}



bool BellmanFord(Graph &graph, size_t src, std::vector<int64_t> &dist) {
    // Return true if can be solved
    dist.clear();
    dist.resize(graph.size());
    std::fill(dist.begin(), dist.end(), INT64_MAX);
    dist[src] = 0;


    std::vector<Edge> edges;
    for (size_t node : graph) {
        for (auto &neib : graph[node]) {
            edges.emplace_back(node, neib.first, neib.second);
        }
    }


    for (int i = 0; i < graph.size() - 1; i++) {
        for (auto &edge : edges) {
            if (dist[edge.first] + edge.weight < dist[edge.second]) {
                dist[edge.second] = dist[edge.first] + edge.weight;
            }
        }
    }

    // check for negative-weight cycles.
    for (auto &edge : edges) {
        size_t x = edge.first;
        size_t y = edge.second;
        int32_t weight = edge.weight;

        if (dist[x] != INT64_MAX && dist[x] + weight < dist[y]) {
            return false;  // TODO Raise error
        }
    }

    return true;
}





void FloydWarshall(Graph &graph, std::vector<std::vector<int64_t>> &dist) {
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
    dist.clear();
    dist.resize(graph.size());

    // Run belman on dummy graph
    Graph graph_prime(graph.size() + 1);

    std::vector<Edge> edges;
    for (size_t node : graph) {
        for (auto &neib : graph[node]) {
            edges.emplace_back(node, neib.first, neib.second);
        }

        edges.emplace_back(graph.size(), node, 0);  // Dummy node
    }

    graph_prime.build_directed(edges);


    std::vector<int64_t> belman_dist;

    if (!BellmanFord(graph_prime, graph.size(), belman_dist)) {
        std::cout << "Graph contains negative weight cycle" << std::endl;  // TODO: raise error
        return false;
    }

    // Update weights
    Graph graph_updated(graph.size());
    std::vector<Edge> edges_updated;
    for (size_t node : graph) {
        for (auto &neib : graph[node]) {
            if (neib.second != 0) {
                edges.emplace_back(node, neib.first, neib.second + belman_dist[node] - belman_dist[neib.first]);
            } else {
                edges.emplace_back(node, neib.first, 0);
            }
        }
    }

    graph_updated.build_directed(edges);


    // Run Dijkstra
    for (size_t node : graph) {
        Dijkstra(graph_updated, node, dist[node]);
    }


    return true;
}