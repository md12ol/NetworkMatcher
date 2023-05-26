import copy
import random
import os
from math import exp, log
import matplotlib.pyplot as plt
from typing import List

alpha = 0.3
outp = "./NetMatchFigs/"


def main():
    print("HELLO!")

    file = "./NetMatchIn/raw_network80.dat"
    num_nodes = 80
    adj = [[0 for _ in range(num_nodes)] for _ in range(num_nodes)]
    edges = 0
    tot_weight = 0
    weight_cnt = [0 for _ in range(5)]
    with open(file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line != '\n':
                line = line.rstrip('\n')
                li = line.split('\t')
                edges += 1
                tot_weight += int(li[2])
                weight_cnt[int(li[2])-1] += 1
                adj[int(li[0])][int(li[1])] = int(li[2])
                adj[int(li[1])][int(li[0])] = int(li[2])
                pass
            pass
        pass

    lists = []
    with open(outp + "dublin_graph80_shuffle.dat", "w") as f:
        f.write(str(num_nodes) + " " + str(edges) + " " + str(tot_weight) + "\n")
        for val in weight_cnt:
            f.write(str(val) + " ")
            pass
        f.write("\n")
        for rowIdx in range(num_nodes):
            li = []
            for colIdx in range(num_nodes):
                for _ in range(int(adj[rowIdx][colIdx])):
                    li.append(int(colIdx))
                    f.write(str(colIdx) + " ")
                    pass
                pass
            lists.append(li)
            f.write("\n")
            pass
        pass

    print("DONE!")
    pass

main()