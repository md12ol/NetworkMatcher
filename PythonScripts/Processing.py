import math
import os
import sys
from math import cos, pi, sin
from operator import itemgetter

import matplotlib.pyplot as plt
from graphviz import Graph
import numpy as np

inp = "../../Conferences and Papers/2023 CIBCB/NetworkMatching/NetMatchOut/"
outp = "./NetMatchFigs/"
james_inp = "./JamesData/"
finame = "best.dat"
graph_in = "./NetMatchIn/"
graph_check_dir = "./NetMatchOut/"
samps = 30
precision = 6
col_width = 8 + precision


def make_nets_in_fold(parent_dir: str):
    file_names = os.listdir(parent_dir)
    nets = []
    names = []

    for fname in file_names:
        if fname.__contains__("graphs"):
            with open(parent_dir + fname) as f:
                lines = f.readlines()
                lines = lines[1:]
                nets.append(lines)
                names.append(fname)
                pass
            pass
        pass

    for idx, dat in enumerate(nets):
        graph_path = "JamesGraphs/"
        if not os.path.exists(graph_path):
            os.makedirs(graph_path)
            pass
        if names[idx].__contains__("large"):
            make_graph(edge_list(dat, 200), [], [], graph_path + names[idx], 200, -1)
        elif names[idx].__contains__("small"):
            make_graph(edge_list(dat, 80), [], [], graph_path + names[idx], 80, -1)
            pass
        pass
    pass


def get_base_data(file_path: str):
    data = []
    with open(file_path) as f:
        vals = []
        lines = f.readlines()
        for line in lines:
            if line == "\n":
                data.append(vals)
                vals = []
                pass
            elif not line.__contains__("EE"):
                vals.append(float(line))
                pass
            pass
        data.append(vals)
        pass
    mins = []
    for dat in data:
        mins.append(min(dat))
        pass

    to_return = [[i + 1, mins[i], data[i]] for i in range(len(data))]
    # to_return.sort(key=itemgetter(1))  # Ascending
    return to_return


def get_james_data():
    data = []

    # Large
    size_data = []
    exp_dat = get_base_data(james_inp + "rawlargeham1.txt")
    size_data.append(exp_dat)
    exp_dat = get_base_data(james_inp + "rawlargeham5.txt")
    size_data.append(exp_dat)
    exp_dat = get_base_data(james_inp + "rawlargeham50.txt")
    size_data.append(exp_dat)
    data.append(size_data)

    # Small
    size_data = []
    exp_dat = get_base_data(james_inp + "rawsmallham1.txt")
    size_data.append(exp_dat)
    exp_dat = get_base_data(james_inp + "rawsmallham5.txt")
    size_data.append(exp_dat)
    exp_dat = get_base_data(james_inp + "rawsmallham50.txt")
    size_data.append(exp_dat)
    data.append(size_data)
    return data


def writeStat(data: [], out, lower_better):
    mean = float(np.mean(data))
    mean = round(mean, precision)
    std = float(np.std(data, ddof=0))
    std = round(std, precision)  # Population standard deviation
    diff = 1.96 * std / math.sqrt(30)  # 95% CI
    diff = round(diff, precision)
    if lower_better:
        maxima = float(min(data))
        pass
    else:
        maxima = float(max(data))
        pass
    maxima = round(maxima, precision)
    out.write(str(maxima).ljust(col_width))
    out.write(str(mean).ljust(col_width))
    out.write(str(std).ljust(col_width))
    out.write(str(diff).ljust(col_width))
    return mean, maxima


def make_table(many_data: [], best_run, exp_info: [], fname: str, minimizing: bool):
    with open(fname, "w") as f:
        f.write("EXP".ljust(col_width))
        f.write("Parameters".ljust(3 * col_width))
        f.write("Best Run".ljust(col_width))
        f.write("Best Fit".ljust(col_width))
        f.write("Mean".ljust(col_width))
        f.write("SD".ljust(col_width))
        f.write("95% CI".ljust(col_width))
        f.write("\n")
        for di, data in enumerate(many_data):
            if di < 2:
                f.write(str("EE PS" + str(di + 1)).ljust(col_width))
            else :
                f.write(str("EXP" + str(di + -1)).ljust(col_width))
                pass
            f.write(exp_info[di].ljust(3 * col_width))
            f.write(str("Run " + str(best_run[di])).ljust(col_width))
            # f.write(str("Run " + str(data[0][0])).ljust(col_width))
            writeStat(data, f, minimizing)
            f.write("\n")
            pass
        pass
    pass


def box_plot(bp, num_splits: int, split_info: []):
    for whisker in bp['whiskers']:
        whisker.set(color='#8B008B', linewidth=1)
        pass

    for cap in bp['caps']:
        cap.set(color='#8B008B', linewidth=1)
        pass

    for median in bp['medians']:
        median.set(color='#FF7F0E', linewidth=1)
        pass

    for flier in bp['fliers']:
        flier.set(marker='.', color='#e7298a', alpha=0.5, markersize=5)
        pass

    for info in split_info:
        for idx in info[0]:
            bp['boxes'][idx].set(facecolor=info[1][idx % len(info[1])])
            pass
        pass
    pass


def calc(data):
    mean = float(np.mean(data))
    std = float(np.std(data, ddof=0))
    diff = 1.96 * std / math.sqrt(30)  # 95% CI
    print(str(mean) + "+-" + str(diff))
    pass


def high_low_deg(el: [], verts: int):
    deg = [int(0) for _ in range(verts)]
    low_deg = []
    high_deg = []
    for ed in el:
        deg[ed[0]] += 1
        deg[ed[1]] += 1
        pass
    max = 0
    for idx, deg in enumerate(deg):
        if deg > max:
            max = deg
            pass
        if deg > 20:
            high_deg.append(idx)
            pass
        elif deg > 10:
            low_deg.append(idx)
            pass
        pass
    print("Max: " + str(max))
    return low_deg, high_deg


def edge_list(linked_list, verts: int):
    adjM = [[0 for _ in range(verts)] for _ in range(verts)]
    edges = 0
    weight = 0
    for from_node, line in enumerate(linked_list):
        line = line.rstrip()
        line = line.split(" ")
        for to_node in line:
            if to_node != '':
                adjM[from_node][int(to_node)] += 1
                # adjM[int(to_node)][from_node] += 1
                pass
            pass
        pass

    for row in range(verts):
        for col in range(row + 1, verts):
            val1 = adjM[row][col]
            val2 = adjM[col][row]
            if max(val1, val2) > 0:
                edges += 1
                weight += max(val1, val2)
                pass

            if val1 != val2:
                adjM[row][col] = max(val1, val2)
                adjM[col][row] = max(val1, val2)
                pass
            pass
        pass

    edge_lists = []
    for row in range(verts):
        for col in range(row + 1, verts):
            if adjM[row][col] > 0:
                edge_lists.append([row, col, adjM[row][col]])
                pass
            pass
        pass
    return edge_lists, adjM, edges, weight


def make_graph(el: [], low_deg: [], high_deg: [], out_file: str, verts: int, p0: int):
    g = Graph(engine='sfdp')
    e_cout = 0

    g.graph_attr.update(dpi='1000', size="6,6", outputorder='edgesfirst', overlap='false', splines='true')
    # g.node_attr.update(color='black', shape='point', width='0.02', height='0.02')
    g.node_attr.update(color='black', shape='circle', fixedsize='true', width='0.25', fontsize='8')
    g.edge_attr.update(color='black', penwidth='2')

    for i in range(verts):
        g.node(str(i))
        pass

    for n in range(verts):
        if n == p0:
            if n in low_deg:
                g.node(str(n), label=str(n), color='red', width='0.03', height='0.03')
                pass
            elif n in high_deg:
                g.node(str(n), label=str(n), color='red', width='0.04', height='0.04')
                pass
            else:
                g.node(str(n), label=str(n), color='red')
                pass
        elif n in low_deg:
            g.node(str(n), label=str(n), width='0.03', height='0.03')
        elif n in high_deg:
            g.node(str(n), label=str(n), width='0.04', height='0.04')
        else:
            g.node(str(n), label=str(n))
        pass

    for idx, d in enumerate(el[0]):
        if d[0] < d[1]:
            if d[2] == 1:
                g.edge(str(d[0]), str(d[1]), color='black')
                pass
            elif d[2] == 2:
                g.edge(str(d[0]), str(d[1]), color='#ff0000')
                pass
            elif d[2] == 3:
                g.edge(str(d[0]), str(d[1]), color='#00ff00')
                pass
            elif d[2] == 4:
                g.edge(str(d[0]), str(d[1]), color='#0000ff')
                pass
            else:
                g.edge(str(d[0]), str(d[1]), color='#87cefa')
                pass
            e_cout += 1
            pass
        pass
    g.render(filename=out_file, directory=outp, cleanup=True, format='png')
    # g.save(filename=out_file, directory=outp)
    # g.clear()
    print("Made network: " + out_file + " with " + str(e_cout) + " Edges")
    pass


def check_vals(dat: [], size: int):
    el, am, edges, weight = edge_list(dat[3], size)
    if edges != dat[4]:
        return False
    if weight != dat[5]:
        return False
    return True


def check_fitness(network: [], fitness: int, net_path: str, penalty: int):
    with open(net_path) as f:
        dub_lines = f.readlines()
        size = int(dub_lines[0].split(" ")[0])
        dub_lines = dub_lines[2:]
        pass

    net_fitness, only_net, only_dub, both = hammy_distance(network, dub_lines, size, penalty)
    # print("Fitness found was " + str(net_fitness) + " and it should be " + str(fitness))
    if net_fitness == fitness:
        return True, only_net, only_dub, both
    else:
        return False, only_net, only_dub, both


def hammy_distance(g1_lines: [], g2_lines: [], size: int, penalty: int):
    cost = 0
    only_net = 0
    only_dub = 0
    both = 0
    edge_list1, adj_matrix1, _, _ = edge_list(g1_lines, size)
    edge_list2, adj_matrix2, _, _ = edge_list(g2_lines, size)

    # print(len(adj_matrix1))
    # print(len(adj_matrix1[0]))
    # print(len(adj_matrix2))
    # print(len(adj_matrix2[0]))

    for row in range(size):
        for col in range(row + 1, size):
            count1 = adj_matrix1[row][col]
            count2 = adj_matrix2[row][col]
            if count1 != count2:
                if count1 == 0:
                    cost += penalty
                    only_dub += 1
                elif count2 == 0:
                    cost += penalty
                    only_net += 1
                else:
                    both += 1
                    cost += min(abs(count1 - count2), penalty)
                    pass
            elif count1 > 0 and count1 == count2:
                both += 1
            pass
        pass
    return cost, only_net, only_dub, both


def combine(dir: str, start: str, end: str):
    if os.path.exists(dir + start + end):
        os.remove(dir + start + end)
        pass
    out_file = open(dir + start + end, "a")
    for idx in range(1, 31):
        with open(dir + start + str(idx).zfill(2) + end, "r") as f:
            out_file.write(f.read())
            pass
        pass
    out_file.close()
    pass


def get_data(dir_path: str, net_path: str, size: int, penalty: int):
    fits = []
    networks = []
    SDAs = []
    edges = []
    weights = []
    hists = []
    fi_str = finame
    with open(dir_path + fi_str) as f:
        lines = f.readlines()
        next_graph = False
        cont_graph = False
        next_SDA = False
        cont_SDA = False
        network = []
        SDA = []
        for line in lines:
            if line.__contains__(str("-fitness")):
                d = line.split(" ")
                fits.append(float(d[0]))
                if cont_graph:
                    networks.append(network)
                    network = []
                cont_graph = False
                pass
            elif line.__contains__(str("Self-Driving Automata")):
                next_SDA = True
                pass
            elif line.__contains__("Graph"):
                cont_SDA = False
                SDAs.append(SDA)
                SDA = []
                pass
            elif line.__contains__("Edges: "):
                edges.append(int(line.rstrip("\n").split(" ")[1]))
                pass
            elif line.__contains__("Tot Weight: "):
                weights.append(int(line.rstrip("\n").split(" ")[2]))
                pass
            elif line.__contains__("W Hist: "):
                hists.append([int(v) for v in line.rstrip("\n").split(" ")[3:8]])
                next_graph = True
                pass
            elif next_SDA:
                next_SDA = False
                SDA.append(line)
                cont_SDA = True
                pass
            elif cont_SDA:
                SDA.append(line)
                pass
            elif next_graph:
                next_graph = False
                network.append(line)
                cont_graph = True
                pass
            elif cont_graph:
                network.append(line)
                pass
            pass
        networks.append(network)
        pass

    # run number, fitness, profileS, dnaS, edge count
    if len(fits) != samps:
        print("ERROR in fits: " + dir_path)
        pass
    if len(networks) != samps:
        print("ERROR in networks: " + dir_path)
        pass
    if len(SDAs) != samps:
        print("ERROR in SDAs: " + dir_path)
        pass
    if len(edges) != samps:
        print("ERROR in edges: " + dir_path)
        pass
    if len(weights) != samps:
        print("ERROR in weights: " + dir_path)
        pass
    if len(hists) != samps:
        print("ERROR in hists: " + dir_path)
        pass
    edge_cmpr = []
    for idx, fit in enumerate(fits):
        okay, only_net, only_dub, both = check_fitness(networks[idx], int(fit), net_path, penalty)
        edge_cmpr.append([only_net, only_dub, both])
        if not okay:
            print("ERROR in fitness: " + dir_path + " on run " + str(idx + 1))
            pass
        pass

    data = [[i + 1, fits[i], SDAs[i], networks[i], edges[i], weights[i], hists[i], edge_cmpr[i]] for i in range(samps)]
    data.sort(key=itemgetter(1))  # Ascending

    for dat in data:
        if not check_vals(dat, size):
            print("ERROR unknown.")
            pass
        pass
    return data


def print_best_info(path: str, info: [], bests: [], size: int):
    with open(path, "w") as f:
        for didx, dat in enumerate(bests):
            f.write(info[didx] + "\n")
            pass
        f.write("\n")

        edge_cmp_info = ["Test Only: ", "True Only: ", "Both: "]
        for didx, data in enumerate(bests):
            dat = data[1][0]
            f.write("Experiment " + str(data[0]) + " Info:\n")
            f.write("Run Number: " + str(dat[0]) + "\n")
            f.write("With Fitness: " + str(dat[1]) + "\n")
            f.write("Edge Comparison: ")
            for idx, val in enumerate(dat[7]):
                f.write(edge_cmp_info[idx])
                f.write(str(val) + "\t")
                pass
            f.write("\nGraph:\n")
            f.write("Nodes Edges Weight: " + str(size) + " " + str(dat[4]) + " " + str(dat[5]) + "\n")
            f.write("Weight Histogram (1-5): ")
            for val in dat[6]:
                f.write(str(val) + " ")
                pass
            f.write("\n")
            f.writelines(dat[3])
            f.write("SDA:\n")
            f.writelines(dat[2])
            f.write("\n\n")
            pass
        pass
    pass


def main():
    print("START")
    folder_names = os.listdir(inp)
    mode_itms = [["200", "NetMatch"], ["80", "NetMatch"]]
    # mode_itms = [["80", "NetMatch"]]
    mode_info = ["NM200", "NM80"]
    # mode_info = ["NM80"]
    dub_graph_files = ["dublin_graph.dat", "dublin_graph80.dat"]
    # dub_graph_files = ["dublin_graph80.dat"]
    # dub_graph_files = ["dublin_graph.dat"]
    dub_graph_size = [200, 80]
    # dub_graph_size = [80]
    mode_dirs = [[] for _ in range(len(mode_itms))]
    muts = [2, 8]
    states = [8, 16, 24]
    # penalty = [1, 5, 50]
    penalty = [1, 5, 50]

    for fold in folder_names:
        combine(inp + fold + "/", "best", ".dat")
        pass

    make_james_nets = False
    if make_james_nets:
        make_nets_in_fold(james_inp)
        pass

    for idx, graph in enumerate(dub_graph_files):
        with open(graph_in + graph) as f:
            dub_lines = f.readlines()
            print(dub_lines[0:2])
            dub_lines = dub_lines[2:]
            size = dub_graph_size[idx]
            # make_graph(edge_list(dub_lines, size), [], [], "DublinNetwork" + str(size), size, 0)
            pass
        pass

    james_exp_lbls = ["EdgeEdit 1", "EdgeEdit 2"]
    exp_lbls = [[] for _ in range(len(penalty))]
    exp_dat = [[] for _ in range(len(penalty))]
    james_exp_descriptions = ["EdgeEdit PS1", "EdgeEdit PS2"]
    exp_descriptions = [[] for _ in range(len(penalty))]
    for pidx, p in enumerate(penalty):
        for idx, dat in enumerate(james_exp_descriptions):
            exp_lbls[pidx].append(james_exp_lbls[idx])
            exp_descriptions[pidx].append(dat)
            pass
        exp_idx = 1
        for s in states:
            for m in muts:
                exp_dat[pidx].append([str(s) + "S", str(m) + "M"])
                exp_lbls[pidx].append(str(exp_idx) + "(" + str(s) + ", " + str(m) + ")")
                exp_descriptions[pidx].append(str(s) + " States, " + str(m) + " Muts, " + str(p) + " Penalty")
                exp_idx += 1
                pass
            pass
        pass

    for midx, itms in enumerate(mode_itms):
        for pidx, p in enumerate(penalty):
            pen_dirs = []
            for dat in exp_dat[pidx]:
                for fld in folder_names:
                    if all(itm in fld for itm in itms) and all(d in fld for d in dat) and str(str(p) + "P") in fld:
                        pen_dirs.append(fld)
                        pass
                    pass
                pass
            mode_dirs[midx].append(pen_dirs)
            pass
        pass

    # mode_data[mode][penalty][exp][run][0] = run num
    # mode_data[mode][penalty][exp][run][1] = run's fit (sorted based on this)
    # mode_data[mode][penalty][exp][run][2][:] = run's SDA
    # mode_data[mode][penalty][exp][run][3][:] = run's network
    # mode_data[mode][penalty][exp][run][4] = run's network's edges
    # mode_data[mode][penalty][exp][run][5] = run's network's weight
    # mode_data[mode][penalty][exp][run][6][:] = run's network's weight hist
    # mode_data[mode][penalty][exp][run][7][:] = run's network's edge comparison (in net, in dub, both)
    mode_data = [[[] for _ in range(len(penalty))] for _ in range(len(mode_itms))]
    for midx, itms in enumerate(mode_itms):
        for pidx, p in enumerate(penalty):
            for fld in mode_dirs[midx][pidx]:
                mode_data[midx][pidx].append(get_data(inp + fld + "/", graph_in + dub_graph_files[midx], dub_graph_size[
                    midx], p))
                pass
            pass
        pass

    # james_data[mode][penalty][exp] = [run's fitness vals]
    james_data = get_james_data()

    # mode_stats[mode][exp] = [run's fitness vals]
    mode_stats = [[[] for _ in range(len(penalty))] for _ in range(len(mode_itms))]
    make_all = False
    make_any = False
    for midx, itms in enumerate(mode_itms):
        for pidx, pen in enumerate(penalty):
            mode_path = mode_info[midx] + "_Pen" + str(pen) + "_Networks/"
            if make_any and not os.path.exists(outp + mode_path):
                os.makedirs(outp + mode_path)
                pass
            mode_stats[midx][pidx].append(james_data[midx][pidx][0][2])
            mode_stats[midx][pidx].append(james_data[midx][pidx][1][2])
            for expidx, exp in enumerate(mode_data[midx][pidx]):
                exp_fits = []
                for runidx, run in enumerate(exp):
                    exp_fits.append(run[1])
                    if make_any:
                        if runidx == 0 or make_all:
                            make_graph(edge_list(run[3], int(itms[0])), [], [], mode_path + "EXP" + str(expidx + 1) +
                                       "Run" + str(run[0]), int(itms[0]), 0)
                        pass
                    pass
                mode_stats[midx][pidx].append(exp_fits)
                pass
            pass
        pass

    titles = ["Network Matching 200", "Network Matching 80"]
    # titles = ["Network Matching 80"]
    names = ["NM200_boxplot", "NM80_boxplot"]
    # names = ["NM80_boxplot"]
    # xsp = [[i for i in range(len(all_data[0]))], [i for i in range(len(all_data[1]))]]
    # xpos = [xsp[0], xsp[1], xsp[0], xsp[1], xsp[0], xsp[1], xsp[0], xsp[1]]
    ylb = ["Fitness", "Fitness"]
    xlb = ["Experiment (SDA States, Max Mutations)",
           "Experiment (SDA States, Max Mutations)"]

    # lxpos = []
    # for i in range(2, len(all_data[0]) - 3, 3):
    #     lxpos.append(i + 0.5)
    #     pass
    base_colors = ["#00AAFF", "#00FFAA"]
    colors = ['#3333FF', '#33FF33']
    for idx in range(len(titles)):
        for pidx, pen in enumerate(penalty):
            plt.rc('xtick', labelsize=10)
            plt.rc('ytick', labelsize=10)

            f = plt.figure()
            f.set_figheight(4.5)
            f.set_figwidth(8)
            plot = f.add_subplot(111)

            bp = plot.boxplot(mode_stats[idx][pidx], patch_artist=True)
            box_plot(bp, 2, [[[i for i in range(0, 2)], base_colors],
                             [[i for i in range(2, 8)], colors]])

            # plot.set_xticks(xpos[idx])
            plot.set_xticklabels(exp_lbls[pidx], rotation=90)

            # f.suptitle(titles[idx] + " with Penalty " + str(pen), fontsize=14)
            plot.set_xlabel(xlb[idx], fontsize=12)
            plot.set_ylabel(ylb[idx], fontsize=12)
            # for x in lxpos:
            #     plot.axvline(x=x, color='black', linestyle='--', linewidth=0.75)
            #     pass
            plot.grid(visible="True", axis="y", which='major', color="darkgray", linewidth=0.75)
            f.tight_layout()
            f.savefig(outp + "Pen" + str(pen) + names[idx] + ".png", dpi=300)
            plt.close()
            pass
        pass

    for midx, mode in enumerate(mode_itms):
        for pidx, pen in enumerate(penalty):
            best_runs = [-1, -1]
            best_of_best_val = sys.maxsize
            bests = []
            for didx, dat in enumerate(mode_data[midx][pidx]):
                best_runs.append(dat[0][0])
                if dat[0][1] < best_of_best_val:
                    best_of_best_val = dat[0][1]
                pass

            for didx, dat in enumerate(mode_data[midx][pidx]):
                if dat[0][1] == best_of_best_val:
                    bests.append([didx + 1, dat])
                    pass
                pass

            make_table(mode_stats[midx][pidx], best_runs, exp_descriptions[pidx], outp +
                       str(dub_graph_size[midx]) + "tablePen" + str(pen) + ".dat", True)

            infos = []
            for dat in bests:
                info = "Best for " + mode_info[midx] + "_Pen" + str(pen) + " setup: " +\
                   "EXP" + str(dat[0]) + ": " + str(exp_descriptions[pidx][dat[0] + 2])
                infos.append(info)
                pass

            print_best_info(outp + mode_info[midx] + "_Pen" + str(pen) + "_Networks/best.dat", infos,
                            bests, dub_graph_size[midx])
            pass
        pass

    # mode_data[mode][exp][run][0] = run num
    # mode_data[mode][exp][run][1] = run's fit (sorted based on this)
    # mode_data[mode][exp][run][2][:] = run's SDA
    # mode_data[mode][exp][run][3][:] = run's network
    # mode_data[mode][exp][run][4] = run's network's edges
    # mode_data[mode][exp][run][5] = run's network's weight
    # mode_data[mode][exp][run][6][:] = run's network's weight hist
    # net_stats = [[[[] for _ in range(2)] for _ in range(len(sizes))] for _ in range(len(mode_itms))]
    # net_bests = [[[[] for _ in range(2)] for _ in range(len(sizes))] for _ in range(len(mode_itms))]
    # for midx, itms in enumerate(mode_itms):
    #     for sidx, size in enumerate(sizes):
    #         for expidx, exp in enumerate(mode_data[midx][sidx]):
    #             exp_edges = []
    #             exp_weights = []
    #             for runidx, run in enumerate(exp):
    #                 if runidx == 0:
    #                     net_bests[midx][sidx][0].append(run[4])
    #                     net_bests[midx][sidx][1].append(run[5])
    #                     pass
    #                 exp_edges.append(run[4])
    #                 exp_weights.append(run[5])
    #                 pass
    #             net_stats[midx][sidx][0].append(exp_edges)
    #             net_stats[midx][sidx][1].append(exp_weights)
    #             pass
    #         pass
    #     pass

    # titles = ["Epidemic Length", "Epidemic Profile Matching P1", "Epidemic Profile Matching P7",
    #           "Epidemic Profile Matching Dublin"]
    # names = ["EL_netstats", "PM1_netstats", "PM7_netstats", "PMDUB_netstats"]
    # ylb = ["Network Edges", "Network Weight"]
    # xlb = ["Experiment (Num. States, Max. Mutations)",
    #        "Experiment (Num. States, Max. Mutations)",
    #        "Experiment (Num. States, Max. Mutations)",
    #        "Experiment (Num. States, Max. Mutations)", ]
    # out_path = outp + "Network Stats Boxplots"
    # if not os.path.exists(out_path):
    #     os.makedirs(out_path)
    #     pass
    # out_path += "/"
    # xsp = [[i for i in range(1, 10)], [i - 0.22 for i in range(1, 10)], [i + 0.22 for i in range(1, 10)]]
    # for idx in range(len(titles)):
    #     for sidx, size in enumerate(sizes):
    #         if idx < 3 or sidx == 0:
    #             if idx >= 3:
    #                 size = 200
    #                 pass
    #             plt.rc('xtick', labelsize=6)
    #             plt.rc('ytick', labelsize=6)
    #
    #             f = plt.figure()
    #             f.set_figheight(5)
    #             f.set_figwidth(8)
    #             plot = f.add_subplot(111)
    #             plot2 = plot.twinx()
    #
    #             # plot2.set_axisbelow(True)
    #             # plot2.grid(visible="True", axis="y", which='major', color="lightgray", linewidth=0.75)
    #
    #             bp = plot.boxplot(net_stats[idx][sidx][0], positions=xsp[1], patch_artist=True, zorder=1, widths=0.4)
    #             box_plot(bp, 1, [[[i for i in range(9)], ["#0000FF"]]])
    #             plot.plot(xsp[1], net_bests[idx][sidx][0], 'x', color="#FF0000")
    #             bp = plot2.boxplot(net_stats[idx][sidx][1], positions=xsp[2], patch_artist=True, zorder=2, widths=0.4)
    #             box_plot(bp, 1, [[[i for i in range(9)], ["#00FF00"]]])
    #             plot2.plot(xsp[2], net_bests[idx][sidx][1], 'x', color="#FF0000")
    #
    #             plot.hlines([size * 1, size * 4], 0.5, 9.5, colors="#0000FF", linestyles="dashed", linewidth=1)
    #             plot2.hlines([size * 4, size * 16], 0.5, 9.5, colors="#00FF00", linestyles="dotted", linewidth=1)
    #
    #             plot.set_xticks(xsp[0])
    #             plot.set_xticklabels(exp_lbls[2:], rotation=90)
    #
    #             if not titles[idx].__contains__("Dublin"):
    #                 f.suptitle(titles[idx] + " w " + str(size) + " Nodes", fontsize=12)
    #             else:
    #                 f.suptitle(titles[idx], fontsize=12)
    #                 pass
    #
    #             plot.set_xlabel(xlb[idx], fontsize=10)
    #             plot.set_ylabel(ylb[0], fontsize=10, color="#0000FF")
    #             plot2.set_ylabel(ylb[1], fontsize=10, rotation=270, labelpad=10, color="#00FF00")
    #
    #             plot.set_axisbelow(True)
    #             plot.grid(visible="True", axis="y", which='major', color="darkgray", linewidth=0.75)
    #             f.tight_layout()
    #
    #             if idx < 3:
    #                 f.savefig(out_path + str(size) + names[idx] + ".png", dpi=450)
    #             else:
    #                 f.savefig(out_path + names[idx] + ".png", dpi=450)
    #                 pass
    #             plt.close()
    #         pass
    #     pass
    print("END")
    pass


main()
