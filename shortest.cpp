
#include "graph.cpp"


void BellmanFord(Graph &graph, size_t src, std::vector<int64_t> &output) {
    std::vector<Edge> edges;

    // TODO: Think do better
    for (size_t node : graph) {
        size_t j = 0;
        for (auto weight : graph[node]) {
            if (weight != 0) {
                edges.emplace_back(node, j, weight);
            }
            ++j;
        }
    }

    for (auto &edge : edges) {
        std::cout << edge.first << "->" << edge.second << " $" << edge.weight << std::endl;
    }

    // Initialize distance of all vertices as 0.
    std::vector<int64_t> dis(graph.size(), INT64_MAX);

    // initialize distance of source as 0
    dis[src] = 0;


    for (int i = 0; i < graph.size() - 1; i++) {
        for (auto &edge : edges) {
            if (dis[edge.first] + edge.weight < dis[edge.second]) {
                dis[edge.second] = dis[edge.first] + edge.weight;
            }
        }
    }

    // check for negative-weight cycles.
    for (auto &edge : edges) {
        size_t x = edge.first;
        size_t y = edge.second;
        int32_t weight = edge.weight;

        if (dis[x] != INT64_MAX && dis[x] + weight < dis[y])
            std::cout << "Graph contains negative weight cycle" << std::endl;
    }

    std::cout << "Vertex Distance from Source" << std::endl;
    for (int i = 0; i < graph.size(); i++) {
        std::cout << i << "\t\t" << dis[i] << std::endl;
    }
}





void FloydWarshall(Graph &graph, std::vector<int64_t> &output) {
    std::vector<std::vector<int64_t>> dist(graph.size(), std::vector<int64_t>(graph.size(), 0));


    // TODO: Think do better
    for (size_t node : graph) {
        size_t j = 0;
        for (auto weight : graph[node]) {
            if (weight == 0) {
                dist[node][j] = INT64_MAX;
            } else {
                dist[node][j] = weight;
            }
            ++j;
        }
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


    // Print the shortest distance matrix
    std::cout<<"The following matrix shows the shortest distances between every pair of vertices \n";

    for (int i = 0; i < graph.size(); i++) {
        for (int j = 0; j < graph.size(); j++) {
            if (dist[i][j] == INT64_MAX) {
                std::cout << "INF" << "     ";
            } else {
                std::cout << dist[i][j] << "     ";
            }
        }
        std::cout << std::endl;
    }
}

