# distutils: sources = graph.cpp
# distutils: language = c++

import typing
from libcpp.pair cimport pair
from libcpp.vector cimport vector
from libc.stddef cimport size_t
from libcpp.set cimport set
from libcpp.memory cimport shared_ptr
from cython.operator cimport dereference
from libc.stdint cimport uint32_t, uint64_t, int32_t, int64_t

from graph cimport Graph, GraphGenerator
from graph cimport BellmanFord as BellmanFord_cpp
from graph cimport FloydWarshall as FloydWarshall_cpp
from graph cimport Johnson as Johnson_cpp
#from graph cimport A_Star as A_Star_cpp
from graph cimport Seidel as Seidel_cpp


cdef class PyGraph:
    cdef Graph graph

    # TODO: what about default constructor?
    def __cinit__(self, num: int = 0):
        self.graph = Graph(num)

    def reserve(self, num_of_nodes: int):
        self.graph.reserve(num_of_nodes)

    def build(self, edges: typing.List[typing.List[int, int, int]]):
        # TODO current C++ implementation needs vector of Edges
        # To avoid this (unneeded) overhead, building here via set_neighbour directly
#        self.graph.build(<vector[pair[size_t, size_t]]> edges)
        for edge in edges:
            self.graph.set_neighbour(edge[0], edge[1], edge[2])

    def build(self, edges: typing.Set[typing.List[int, int, int]]):
        # TODO current C++ implementation needs set of Edges
        # To avoid this (unneeded) overhead, building here via set_neighbour directly
#        self.graph.build(<set[pair[size_t, size_t]]> edges)
        for edge in edges:
            self.graph.set_neighbour(edge[0], edge[1], edge[2])

    def build_directed(self, edges: typing.List[typing.List[int, int, int]]):
        # TODO current C++ implementation needs vector of Edges
        # To avoid this (unneeded) overhead, building here via set_neighbour directly
#        self.graph.build_directed(edges)vector
        for edge in edges:
            self.graph.set_neighbour_directed(edge[0], edge[1], edge[2])

    # TODO: accessors

    def remove_neighbour(self, node_index: int, neighbour_index: int):
        self.graph.remove_neighbour(node_index, neighbour_index)

    def remove_neighbour_directed(self, node_index: int, neighbour_index: int):
        self.graph.remove_neighbour_directed(node_index, neighbour_index)

    def set_neighbour(self, node_index: int, neighbour_index: int, weight: int):
        self.graph.set_neighbour(node_index, neighbour_index, weight)

    def set_neighbour_directed(self, node_index: int, neighbour_index: int, weight: int):
        self.graph.set_neighbour_directed(node_index, neighbour_index, weight)

    def add_neighbour(self, node_index: int, neighbour_index: int, weight: int):
        self.graph.add_neighbour(node_index, neighbour_index, weight)

    def add_neighbour_directed(self, node_index: int, neighbour_index: int, weight: int):
        self.graph.add_neighbour_directed(node_index, neighbour_index, weight)

    def check_neighbour(self, node_index: int, neighbour_index: int) -> bool:
        return self.graph.check_neighbour(node_index, neighbour_index)

    def revert_edge(self, node_index: int, neighbour_index: int):
        self.graph.revert_edge(node_index, neighbour_index)

    def node_degree(self, num: int) -> int:
        return self.graph.node_degree(num)

    # TODO: iterators

    def __len__(self):
        return self.graph.size()

cdef class PyGraphGenerator:
    cdef GraphGenerator graph_generator

    def __cinit__(self):
        self.graph_generator = GraphGenerator()

    # TODO: didnt' implement, strange place to construct Edges, same as for Graph.build()
#    @staticmethod
#    def generate_tree_edges(graph: PyGraph, edges: typing.Dict[typing.Tuple[int, int]]):
#        GraphGenerator.generate_tree_edges(graph.graph, edges)

#    @staticmethod
#    def generate_tree_edges(graph: PyGraph, edges: typing.List[typing.Tuple[int, int]]):
#        GraphGenerator.generate_tree_edges(graph.graph, edges)

    @staticmethod
    def build_by_prob(graph: PyGraph, prob: float):
        GraphGenerator.build_by_prob(graph.graph, prob)

    @staticmethod
    def build_by_prob_inv(graph: PyGraph, prob_inv: int):
        GraphGenerator.build_by_prob_inv(graph.graph, prob_inv)

    @staticmethod
    def build_by_number(graph: PyGraph, num_of_egdes: int):
        GraphGenerator.build_by_number(graph.graph, num_of_egdes)



def Dijkstra(graph: PyGraph, src: int) -> typing.List[int]:
    cdef vector[int64_t] output
    BellmanFord_cpp(graph.graph, src, output)
    out: typing.List[int] = []
    for dist in output:
        out.append(dist)
    return out


def BellmanFord(graph: PyGraph, src: int) -> typing.List[int]:
    cdef vector[int64_t] output
    BellmanFord_cpp(graph.graph, src, output)
    out: typing.List[int] = []
    for dist in output:
        out.append(dist)
    return out


def FloydWarshall(graph: PyGraph) -> typing.List[typing.List[int]]:
    cdef vector[vector[int64_t]] output
    FloydWarshall_cpp(graph.graph, output)
    out: typing.List[typing.List[int]] = []
    for line in output:
        out.append([])
        for dist in line:
            out[-1].append(dist)
    return out


def Johnson(graph: PyGraph) -> typing.List[typing.List[int]]:
    cdef vector[vector[int64_t]] output
    Johnson_cpp(graph.graph, output)
    out: typing.List[typing.List[int]] = []
    for line in output:
        out.append([])
        for dist in line:
            out[-1].append(dist)
    return out

def Seidel(graph: PyGraph) -> typing.List[typing.List[int]]:
    cdef vector[vector[int64_t]] output
    Seidel_cpp(graph.graph, output)
    out: typing.List[typing.List[int]] = []
    for line in output:
        out.append([])
        for dist in line:
            out[-1].append(dist)
    return out
