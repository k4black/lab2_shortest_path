from libc.stddef cimport size_t
from libcpp cimport bool
from libcpp.pair cimport pair
from libcpp.vector cimport vector
from libcpp.set cimport set
from libcpp.memory cimport shared_ptr
from libc.stdint cimport uint32_t, uint64_t, int32_t, int64_t


cdef extern from "graph.cpp":

    cdef enum Colour:
        pass

    cdef cppclass Edge:
        Edge() except +
        Edge(size_t first, size_t second, int32_t weight) except +

        bool operator<(const Edge& a) const

    cdef cppclass Graph:
        Graph() except +
        Graph(size_t num_of_nodes) except +

        void reserve(size_t num_of_nodes)
        void build(const vector[Edge] &edges)
        void build(const set[Edge] &edges)
        void build_unweighted(const vector[Edge] &edges)
        void build_directed(const vector[Edge] &edges)
        const vector[size_t] &operator[](size_t node_index) const
        const vector[size_t] &get_neighbours(size_t node_index) const
        const vector[vector[int32_t]] &get_matrix()
        int32_t get_weight(size_t node_index, size_t neighbour_index)
        # shared_ptr get_node(size_t node_index) const
        void remove_neighbour(size_t node_index, size_t neighbour_index)
        void remove_neighbour_directed(size_t node_index, size_t neighbour_index)
        # TODO: weight is size_t?
        void set_neighbour(size_t node_index, size_t neighbour_index, size_t weight)
        # TODO: weight is size_t?
        void set_neighbour_directed(size_t node_index, size_t neighbour_index, size_t weight)
        # TODO: weight is size_t?
        void add_neighbour(size_t node_index, size_t neighbour_index, size_t weight)
        # TODO: weight is size_t?
        void add_neighbour_directed(size_t node_index, size_t neighbour_index, size_t weight)
        bool check_neighbour(size_t node_index, size_t neighbour_index) const
        void revert_edge(size_t node_index, size_t neighbour_index)
        size_t node_degree(size_t node_index) const

        cppclass iterator:
            size_t current_it

            iterator(size_t current) except +
            bool operator==(const iterator &other)
            void operator=(const iterator &other)
            bool operator!=(const iterator &other)
            iterator &operator++()
            size_t operator*()

        # TODO: noexcept here?
        iterator begin()
        iterator end()
        iterator begin() const
        iterator end() const
        size_t size() const

    cdef cppclass GraphGenerator:
        GraphGenerator() except +
        # TODO: didnt' implement, strange place to construct Edges, same as for Graph.build()
#        @staticmethod
#        void generate_tree_edges(const Graph &graph, set[Edge] &edges)
#        @staticmethod
#        void generate_tree_edges(const Graph &graph, vector[Edge] &edges)
        @staticmethod
        void build_by_prob(Graph &graph, float prob_f)
        @staticmethod
        void build_by_prob_inv(Graph &graph, size_t prob)
        @staticmethod
        void build_by_number(Graph &graph, size_t num_of_edges)


cdef extern from "shortest.cpp":

    cdef void Dijkstra(Graph &graph, size_t src, vector[int64_t] &output)
    cdef void DijkstraHeap(Graph &graph, size_t src, vector[int64_t] &output)
    cdef void BellmanFord(Graph &graph, size_t src, vector[int64_t] &output)
    cdef void FloydWarshall(Graph &graph, vector[vector[int64_t]] &output)
    cdef void Johnson(Graph &graph, vector[vector[int64_t]] &output)
    cdef void JohnsonHeap(Graph &graph, vector[vector[int64_t]] &output)
    cdef void A_Star(Graph &graph, size_t src, size_t dest, vector[int32_t] &heur, vector[size_t] &output)
    cdef void Seidel(Graph &graph, vector[vector[int64_t]] &lengths, vector[vector[size_t]] &preds, bool reconstruct)
