# distutils: sources = graph.cpp
# distutils: language = c++

import typing
from libcpp.pair cimport pair
from libcpp.vector cimport vector
from libc.stddef cimport size_t
from libcpp.set cimport set
from libcpp.memory cimport shared_ptr
from cython.operator cimport dereference
from libc.stdint cimport uint32_t, uint64_t

from graph cimport Graph, BucketSorter, RadixSorter, StdSorter
from graph cimport cut_classes_multimap as cut_classes_multimap_cpp, cut_classes_sort as cut_classes_sort_cpp, cut_edges as cut_edges_cpp, find_bridge as find_bridge_cpp
from graph cimport GraphGenerator


cdef class PyGraph:
    cdef Graph graph

    # TODO: what about default constructor?
    def __cinit__(self, num: int = 0):
        self.graph = Graph(num)

    def reserve(self, num_of_nodes: int):
        self.graph.reserve(num_of_nodes)

    def build(self, edges: typing.List[typing.Tuple[int, int]]):
        self.graph.build(<vector[pair[size_t, size_t]]> edges)

    def build(self, edges: typing.Set[typing.Tuple[int, int]]):
        self.graph.build(<set[pair[size_t, size_t]]> edges)

    def build_directed(self, edges: typing.List[typing.Tuple[int, int]]):
        self.graph.build_directed(edges)

    # TODO: accessors

    def remove_neighbour(self, node_index: int, neighbour_index: int):
        self.graph.remove_neighbour(node_index, neighbour_index)

    def remove_neighbour_directed(self, node_index: int, neighbour_index: int):
        self.graph.remove_neighbour_directed(node_index, neighbour_index)

    def set_neighbour(self, node_index: int, neighbour_index: int):
        self.graph.set_neighbour(node_index, neighbour_index)

    def set_neighbour_directed(self, node_index: int, neighbour_index: int):
        self.graph.set_neighbour_directed(node_index, neighbour_index)

    def add_neighbour(self, node_index: int, neighbour_index: int):
        self.graph.add_neighbour(node_index, neighbour_index)

    def add_neighbour_directed(self, node_index: int, neighbour_index: int):
        self.graph.add_neighbour_directed(node_index, neighbour_index)

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

    @staticmethod
    def generate_tree_edges(graph: PyGraph, edges: typing.Dict[typing.Tuple[int, int]]):
        GraphGenerator.generate_tree_edges(graph.graph, edges)

    @staticmethod
    def generate_tree_edges(graph: PyGraph, edges: typing.List[typing.Tuple[int, int]]):
        GraphGenerator.generate_tree_edges(graph.graph, edges)

    @staticmethod
    def build_by_prob(graph: PyGraph, prob: float):
        GraphGenerator.build_by_prob(graph.graph, prob)

    @staticmethod
    def build_by_prob_inv(graph: PyGraph, prob_inv: int):
        GraphGenerator.build_by_prob_inv(graph.graph, prob_inv)

    @staticmethod
    def build_by_number(graph: PyGraph, num_of_egdes: int):
        GraphGenerator.build_by_number(graph.graph, num_of_egdes)


def cut_classes_multimap(graph: PyGraph, precision: str = '32') -> typing.List[typing.List[typing.Tuple[int, int]]]:
    cdef vector[vector[pair[size_t, size_t]]] output
    if precision == '32':
        cut_classes_multimap_cpp[uint32_t](graph.graph, output)
    elif precision == '64':
        cut_classes_multimap_cpp[uint64_t](graph.graph, output)
    else:
        raise Exception(f'{precision} is not a supported precision for algorithm')
    out: typing.List[typing.List[typing.Tuple[int, int]]] = []
    for vec in output:
        out.append([])
        for pair in vec:
            out[-1].append((pair.first, pair.second))
    return out


def cut_classes_sort(graph: PyGraph, precision: str = '32', sort: str = 'std') -> typing.List[typing.List[typing.Tuple[int, int]]]:
    cdef vector[vector[pair[size_t, size_t]]] output
    if sort == 'bucket':
        if precision == '32':
            cut_classes_sort_cpp[uint32_t, BucketSorter](graph.graph, output)
        elif precision == '64':
            cut_classes_sort_cpp[uint64_t, BucketSorter](graph.graph, output)
        else:
            raise Exception(f'{precision} is not a supported precision for algorithm')
    elif sort == 'radix':
        if precision == '32':
            cut_classes_sort_cpp[uint32_t, RadixSorter](graph.graph, output)
        elif precision == '64':
            cut_classes_sort_cpp[uint64_t, RadixSorter](graph.graph, output)
        else:
            raise Exception(f'{precision} is not a supported precision for algorithm')
    elif sort == 'std':
        if precision == '32':
            cut_classes_sort_cpp[uint32_t, StdSorter](graph.graph, output)
        elif precision == '64':
            cut_classes_sort_cpp[uint64_t, StdSorter](graph.graph, output)
        else:
            raise Exception(f'{precision} is not a supported precision for algorithm')
    else:
        raise Exception(f'{sort} is not a supported sorting functor')
    out: typing.List[typing.List[typing.Tuple[int, int]]] = []
    for vec in output:
        out.append([])
        for pair in vec:
            out[-1].append((pair.first, pair.second))
    return out

def cut_edges(graph: PyGraph, precision: str = '32') -> typing.List[typing.Tuple[int, int]]:
    cdef vector[pair[size_t, size_t]] output
    if precision == '32':
        cut_edges_cpp[uint32_t](graph.graph, output)
    elif precision == '64':
        cut_edges_cpp[uint64_t](graph.graph, output)
    else:
        raise Exception(f'{precision} is not a supported precision for algorithm')
    out: typing.List[typing.Tuple[int, int]] = []
    for pair in output:
        out.append((pair.first, pair.second))
    return out

def find_bridges(graph: PyGraph) -> typing.List[typing.Tuple[int, int]]:
    cdef vector[pair[size_t, size_t]] output
    find_bridge_cpp(graph.graph, output)
    out: typing.List[typing.Tuple[int, int]] = []
    for pair in output:
        out.append((pair.first, pair.second))
    return out
