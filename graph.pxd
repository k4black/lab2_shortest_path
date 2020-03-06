from libc.stddef cimport size_t, int
from libcpp cimport bool
from libcpp.pair cimport pair
from libcpp.vector cimport vector
from libcpp.set cimport set
from libcpp.memory cimport shared_ptr
from libc.stdint cimport uint32_t, uint64_t


cdef extern from "graph.cpp":

    cdef enum Colour:
        pass

    cdef cppclass Graph:
        Graph() except +
        Graph(size_t num_of_nodes) except +

        void reserve(size_t num_of_nodes)
        void build(const vector[pair[size_t, size_t]] &edges)
        void build(const set[pair[size_t, size_t]] &edges)
        void build_directed(const vector[pair[size_t, size_t]] &edges)
        const vector[size_t] &operator[](size_t node_index) const
        const vector[size_t] &get_neighbours(size_t node_index) const
        # shared_ptr get_node(size_t node_index) const
        void remove_neighbour(size_t node_index, size_t neighbour_index)
        void remove_neighbour_directed(size_t node_index, size_t neighbour_index)
        void set_neighbour(size_t node_index, size_t neighbour_index)
        void set_neighbour_directed(size_t node_index, size_t neighbour_index)
        void add_neighbour(size_t node_index, size_t neighbour_index)
        void add_neighbour_directed(size_t node_index, size_t neighbour_index)
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

        iterator begin()
        iterator end()
        iterator begin() const
        iterator end() const
        size_t size() const

    cdef cppclass GraphGenerator:
        GraphGenerator() except +
        @staticmethod
        void generate_tree_edges(const Graph &graph, set[pair[size_t, size_t]] &edges)
        @staticmethod
        void generate_tree_edges(const Graph &graph, vector[pair[size_t, size_t]] &edges)
        @staticmethod
        void build_by_prob(Graph &graph, float prob_f)
        @staticmethod
        void build_by_prob_inv(Graph &graph, size_t prob)
        @staticmethod
        void build_by_number(Graph &graph, size_t num_of_edges)


cdef extern from "shortest.cpp":

    cdef void BellmanFord(Graph &graph, size_t src, vector[int64_t] &output)
    cdef void FloydWarshall(Graph &graph, vector[vector[int64_t]] &output)
