import random
import typing as tp
import matplotlib.pyplot as plt
import time
import numpy as np
import random
import math
from graph_python import PyGraph, PyEdge, A_star, Seidel

if __name__ == '__main__':
    graph1 = PyGraph(6)
    graph1.build_unweighted([[0, 1, 1],
                             [0, 2, 1],
                             [1, 2, 1],
                             [2, 4, 1],
                             [1, 3, 1],
                             [3, 5, 1]])

    lengths, preds = Seidel(graph1)
    #
    # for line in lengths:
    #     for elem in line:
    #         print(elem, end=' ')
    #     print()
    # print()
    #
    # for line in preds:
    #     for elem in line:
    #         print(elem, end=' ')
    #     print()
    # print()
