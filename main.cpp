#include <iostream>


#include "graph.cpp"


int main() {
    Graph graph(10);

    for (auto i : graph) {
        std::cout << i << " ";
    }

    return 0;
}