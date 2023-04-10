import copy
import random
import os
from math import exp, log
import matplotlib.pyplot as plt
from typing import List

alpha = 0.3
outp = "./BitProcessed/"


def infected(sick: int):
    beta = 1 - exp(sick * log(1 - alpha))
    if random.random() < beta:
        return True
    return False


def fitness_bare(adj_lists: List[List[int]], nodes: int, p0):
    temp_list = copy.deepcopy(adj_lists)
    n_state = [0 for _ in range(nodes)]  # susceptible
    n_state[p0] = 1
    epi_log = [[p0]]
    num_infected = 1
    ttl_infected = 0
    time_step = 0
    length = 0
    while num_infected > 0 and time_step < nodes:
        current_infected = num_infected / nodes
        inf_neighbours = [0 for _ in range(nodes)]
        current_infected = num_infected / nodes
        for n in range(nodes):
            if n_state[n] == 1:
                for nei in temp_list[n]:
                    inf_neighbours[nei] += 1
        for n in range(nodes):
            if n_state[n] == 0 and inf_neighbours[n] > 0:
                if infected(inf_neighbours[n]):
                    n_state[n] = 3
        ttl_infected += num_infected
        num_infected = 0
        new_inf = []
        for n in range(nodes):
            if n_state[n] == 1:  # infected -> removed
                n_state[n] = 2
            elif n_state[n] == 3:
                n_state[n] = 1
                num_infected += 1
                new_inf.append(n)
        epi_log.append(new_inf)
        length += 1
        time_step += 1
    return epi_log, ttl_infected, length


def make_prof(adj_lists: List[List[int]], nodes: int, p0):
    logs = []
    most = -1
    for _ in range(500):
        prof = fitness_bare(adj_lists, nodes, p0)[0]
        logs.append(prof)
        if len(prof) > most:
            most = len(prof)
            pass
        pass

    avg_prof = [0.0 for _ in range(most)]
    for idx in range(most):
        for pr in logs:
            if idx < len(pr):
                avg_prof[idx] += len(pr[idx])
                pass
            pass
        pass

    for idx in range(most):
        avg_prof[idx] = avg_prof[idx]/500
        pass

    return avg_prof


def get_vals(filename: str):
    dat = []
    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            dat.append(float(line))
            pass
        pass

    if dat[0] != 1:
        dat.insert(0, 1)
        pass

    if dat[-1] != 0.0:
        dat.append(0)
        pass

    return dat


def make_profiles_in_dir(dir: str):
    files = os.listdir(dir)
    for fidx, file in enumerate(files):
        if file.__contains__(".dat"):
            plt.rc('xtick', labelsize=6)
            plt.rc('ytick', labelsize=6)

            f = plt.figure()
            f.set_figheight(5)
            f.set_figwidth(8)
            plot = f.add_subplot(111)

            vals = get_vals(dir + file)
            xs = [i for i in range(len(vals))]
            plot.plot(vals, color="#87cefa", marker='o', mfc='r', mec='black')

            if file.__contains__(str(0)):
                f.suptitle("Dublin Epidemic Profile", fontsize=12)
                # for x,y in zip(xs, vals):
                #     if x %2 == 0:
                #         plot.annotate(str(y), (x, y), textcoords="offset points", xytext=(0, -13), ha='center',
                #                       fontsize=8)
                #     else:
                #         plot.annotate(str(y), (x, y), textcoords="offset points", xytext=(0, 5), ha='center',
                #                       fontsize=8)
                #     pass
                pass
            elif file.__contains__(str(1)):
                f.suptitle("Unimodal Epidemic Profile", fontsize=12)
            elif file.__contains__(str(7)):
                f.suptitle("Bimodal Epidemic Profile", fontsize = 12)
                pass

            plot.set_xlabel("Time Step", fontsize=10)
            plot.set_ylabel("New Infections", fontsize=10)
            plot.grid(visible="True", axis="y", which='both', color="darkgray", linewidth=0.75)
            f.tight_layout()

            if file.__contains__(str(0)):
                f.savefig(outp + "ProfileDublin.png", dpi=300)
            elif file.__contains__(str(1)):
                f.savefig(outp + "Profile1.png", dpi=300)
            elif file.__contains__(str(7)):
                f.savefig(outp + "Profile7.png", dpi=300)
                pass
            pass
        pass


def same(one, two):
    # print(str(len(one)) + "\t" + str(len(two)))
    # if len(one) != len(two):
    #     return False
    for idx in range(max(len(one), len(two))):
        val1 = round(one[idx], 1)
        val2 = round(two[idx], 1)
        if val1 != val2:
            return False
        pass
    return True


def main():
    print("HELLO!")
    make_profiles_in_dir("./Profiles/")

    file = "./raw_network.txt"
    adj = [[0 for _ in range(200)] for _ in range(200)]
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
    with open("dublin_graph_NEW.dat", "w") as f:
        f.write(str(200) + " " + str(edges) + " " + str(tot_weight) + "\n")
        for val in weight_cnt:
            f.write(str(val) + " ")
            pass
        f.write("\n")
        for rowIdx in range(200):
            li = []
            for colIdx in range(200):
                for _ in range(int(adj[rowIdx][colIdx])):
                    li.append(int(colIdx))
                    f.write(str(colIdx) + " ")
                    pass
                pass
            lists.append(li)
            f.write("\n")
            pass
        pass

    with open("./Profiles/Profile0.dat") as f:
        lines = f.readlines()
        other = []
        for line in lines:
            line=line.rstrip()
            other.append(float(line))
            pass
        pass

    print("DONE!")
    pass

main()