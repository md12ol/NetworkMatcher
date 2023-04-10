import math
import os
from math import cos, pi, sin
from operator import itemgetter

import matplotlib.pyplot as plt
from graphviz import Graph
import numpy as np

inp = "../../Conferences and Papers/2023 CIBCB/WSDA/BitOutDone/"
NP_inp = "../../Conferences and Papers/2023 CIBCB/WSDA/BitOutCheck/"
outp = "./BitProcessed/"
james_inp = "../../Conferences and Papers/2023 CIBCB/WSDA/EditBase/"
james_net_inp = "../../Conferences and Papers/2023 CIBCB/WSDA/bestgraphs/"
finame = "best.lint"
samps = 30
precision = 6
col_width = 6 + precision


def get_data(dir_path: str):
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

    if "NM" not in dir_path:
        p0_list = get_patient_zero(dir_path)
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
    if "NM" not in dir_path and len(p0_list) != samps:
        print("ERROR in p0_list: " + dir_path)
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
    if "NM" in dir_path:
        for idx, fit in enumerate(fits):
            if not check_fitness(network[idx], fit):
                print("ERROR in fitness: " + dir_path + " on run " + str(idx + 1))
                pass
            pass
        pass
    if "NM" in dir_path:
        data = [[i + 1, fits[i], SDAs[i], networks[i], edges[i], weights[i], hists[i]] for i in range(samps)]
    else:
        data = [[i + 1, fits[i], SDAs[i], networks[i], edges[i], weights[i], hists[i], p0_list[i]] for i in range(samps)]
        pass
    data.sort(key=itemgetter(1))  # Ascending
    if "ED" in dir_path:
        data.reverse()
        pass
    return data


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

    to_return = [[i+1, mins[i], data[i]] for i in range(len(data))]
    to_return.sort(key=itemgetter(1)) # Ascending
    return data, to_return


def get_james_data():
    data = []
    # ED
    mode_data = []
    exp_dat, all_data = get_base_data(james_inp + "ED128results.dat")
    mode_data.append(exp_dat)
    exp_dat, all_data = get_base_data(james_inp + "ED256results.dat")
    mode_data.append(exp_dat)
    data.append(mode_data)

    # PM 1
    mode_data = []
    exp_dat, all_data = get_base_data(james_inp + "PM128Prof1results.dat")
    mode_data.append(exp_dat)
    exp_dat, all_data = get_base_data(james_inp + "PM256Prof1results.dat")
    mode_data.append(exp_dat)
    data.append(mode_data)

    # PM 7
    mode_data = []
    exp_dat, all_data = get_base_data(james_inp + "PM128Prof7results.dat")
    mode_data.append(exp_dat)
    exp_dat, all_data = get_base_data(james_inp + "PM256Prof7results.dat")
    mode_data.append(exp_dat)
    data.append(mode_data)

    # DUB
    mode_data = []
    exp_dat, all_data = get_base_data(james_inp + "PMdublinResults.dat")
    mode_data.append(exp_dat)
    data.append(mode_data)

    # NM DUB
    mode_data = []
    exp_dat, all_data = get_base_data(james_inp + "dublinresults.txt")
    exp_dat = [all_data[0][2], all_data[1][2]]
    mode_data.append(exp_dat)
    data.append(mode_data)
    return data


def get_patient_zero(dir_path: str):
    p0_list = []
    with open(dir_path + "patient0.dat") as f:
        lines = f.readlines()
        prev = None
        for line in lines:
            if line.__contains__("mating event"):
                li = line.split(" ")
                prev = int(li[4])
                pass
            if line == '\n':
                p0_list.append(prev)
                pass
            pass
        pass
    return p0_list


def writeStat(data: [], out):
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
    out.write(str(mean).ljust(col_width))
    out.write(str(std).ljust(col_width))
    out.write(str(diff).ljust(col_width))
    out.write(str(maxima).ljust(col_width))
    return mean, maxima


def make_table(many_data: [], exp_info: [], fname: str, minimizing: bool):
    with open(fname, "w") as f:
        f.write("EXP".ljust(col_width))
        f.write("Parameters".ljust(2 * col_width))
        f.write("Best Run".ljust(col_width))
        f.write("Mean".ljust(col_width))
        f.write("SD".ljust(col_width))
        f.write("95% CI".ljust(col_width))
        f.write("Best".ljust(col_width))
        f.write("\n")
        for di, data in enumerate(many_data):
            f.write(str("EXP" + str(di)).ljust(col_width))
            f.write(exp_info[di].ljust(2 * col_width))
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
    for from_node, line in enumerate(linked_list):
        line = line.rstrip()
        line = line.split(" ")
        for to_node in line:
            if to_node != '':
                if from_node < int(to_node):
                    adjM[from_node][int(to_node)] += 1
                    adjM[int(to_node)][from_node] += 1
                    pass
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
    return edge_lists, adjM


def make_ring(verts: int):
    el = []
    adjM = [[0 for _ in range(verts)] for _ in range(verts)]
    for v in range(verts):
        el.append([v, (v+1) % verts, 1])
        el.append([v, (v+2) % verts, 1])
        adjM[v][(v+1) % verts] = 1
        adjM[v][(v + 2) % verts] = 1
        pass

    xs = []
    ys = []

    for v in range(verts):
        xs.append(400* float(cos(2*pi*(v/128))) + 500)
        ys.append(400* float(sin(2*pi*(v/128))) + 500)
        pass

    g = Graph(engine='neato')
    e_cout = 0

    # g.graph_attr.update(dpi='1000', size="6,6", outputorder='edgesfirst', overlap='false', splines='true')
    g.graph_attr.update(dpi='1000', size="6,6", overlap='false', splines='true')
    g.node_attr.update(color='black', pin='true', shape='point', width='0.02', height='0.02')
    g.edge_attr.update(color='black', penwidth='0.25')

    for i in range(verts):
        g.node(str(i), pos=str(str(float(xs[i])) + "," + str(float(ys[i]))) + "!")
        pass

    for idx, d in enumerate(el):
        if d[0] < d[1]:
            g.edge(str(d[0]), str(d[1]), color='black')
            e_cout += 1
            pass
        pass

    g.render(filename="128 Ring Network", directory=outp, cleanup=True, format='png')
    print("Made network: " + "128 Ring Network" + " with " + str(e_cout) + " Edges")
    pass


def make_graph(el: [], low_deg: [], high_deg: [], out_file: str, verts: int, p0: int):
    g = Graph(engine='sfdp')
    e_cout = 0

    g.graph_attr.update(dpi='1000', size="6,6", outputorder='edgesfirst', overlap='false', splines='true')
    g.node_attr.update(color='black', shape='point', width='0.02', height='0.02')
    # g.node_attr.update(color='black', shape='circle', fixedsize='true', width='0.25', fontsize='8')
    g.edge_attr.update(color='black', penwidth='0.25')

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


def make_from_dirs(dirs: [], splits: [], titles: []):
    exp_lbls = []
    data = []

    split_dirs = []
    for spl in splits:
        dir_list = []
        for dir in dirs:
            if spl in dir:
                dir_list.append(dir)
                pass
            pass
        split_dirs.append(dir_list)
        pass

    for split in split_dirs:
        lbl_list = []
        split_data = []
        for dir in split:
            lbl_list.append(dir[:10])
            split_data.append(get_data(inp + dir + "/"))
            pass
        exp_lbls.append(lbl_list)
        data.append(split_data)
        pass

    # mode_stats[mode][size][exp] = [run's fitness vals]
    just_fits = []
    for split in data:
        split_fits = []
        for exp in split:
            exp_fits = []
            for run in exp:
                exp_fits.append(run[1])
                pass
            split_fits.append(exp_fits)
            pass
        just_fits.append(split_fits)
        pass

    for split_idx in range(len(splits)):
        plt.rc('xtick', labelsize=8)
        plt.rc('ytick', labelsize=8)

        f = plt.figure()
        f.set_dpi(900)
        f.set_figheight(8)
        f.set_figwidth(10)

        plot = f.add_subplot(111)

        bp = plot.boxplot(just_fits[split_idx], patch_artist=True)
        box_plot(bp, "#FFFFFF")

        # plot.set_xticks(xpos[idx])
        plot.set_xticklabels(exp_lbls[split_idx], rotation=90)

        f.suptitle(titles[split_idx], fontsize=12)
        plot.set_xlabel("Experiment", fontsize=10)
        plot.set_ylabel("Fitness", fontsize=10)
        # for x in lxpos:
        #     plot.axvline(x=x, color='black', linestyle='--', linewidth=0.75)
        #     pass
        plot.grid(visible="True", axis="y", which='major', color="darkgray", linewidth=0.75)
        f.tight_layout()
        f.savefig(outp + titles[split_idx] + ".png", dpi=900)
        plt.close()
        pass
    pass


def james_patient_zeros(filename: str):
    p0s = []
    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            if line.__contains__("EE"):
                line = line.split(" ")
                p0s.append(int(line[1]))
                pass
            pass
        pass
    return p0s


def make_nets_in_fold(parent_dir: str):
    file_names = os.listdir(parent_dir)
    nets = []
    names = []
    patient_zeros = []

    for fname in file_names:
        if fname.__contains__("patientzeroes.txt"):
            patient_zeros = james_patient_zeros(parent_dir + fname)
        else:
            with open(parent_dir + fname) as f:
                lines = f.readlines()
                nets.append(lines)
                names.append(fname)
                pass
            pass

    for idx, dat in enumerate(nets):
        graph_path = "JamesGraphs/"
        if not os.path.exists(graph_path):
            os.makedirs(graph_path)
            pass
        if names[idx].__contains__(str(128)):
            make_graph(edge_list(dat, 128), [], [], graph_path + names[idx], 128, patient_zeros[idx])
        elif names[idx].__contains__(str(256)):
            make_graph(edge_list(dat, 256), [], [], graph_path + names[idx], 256, patient_zeros[idx])
        else:
            make_graph(edge_list(dat, 200), [], [], graph_path + names[idx], 200, patient_zeros[idx])
            pass
        pass
    pass


def get_data_2(dir_path: str):
    fits = []
    networks = []
    SDAs = []
    edges = []
    weights = []
    hists = []
    fi_str = finame
    cnt = 0
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
                cnt += 1
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

    samps = cnt
    p0_list = get_patient_zero(dir_path)
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
    if len(p0_list) != samps:
        print("ERROR in p0_list: " + dir_path)
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

    data = [[i + 1, fits[i], SDAs[i], networks[i], edges[i], weights[i], hists[i]] for i in range(samps)]
    # data.sort(key=itemgetter(1))  # Ascending
    if "ED" in dir_path:
        data.reverse()
        pass
    return data


def network_matching_work(dirs: []):
    data = []
    exp_lbl = []
    for dir in dirs:
        data.append(get_data_2(NP_inp + dir + "/"))
        exp_lbl.append(str(dir))
        pass

    # dub_lines = []
    with open("./dublin_graph.dat") as f:
        dub_lines = f.readlines()
        dub_lines = dub_lines[2:]
        pass

    for exp_idx, exp_dat in enumerate(data):
        for run_dat in exp_dat:
            fit = hammy_distance(run_dat[3], dub_lines)
            if fit != run_dat[1]:
                print("ERROR!! for " + exp_lbl[exp_idx] + " Run " + str(run_dat[0]))
            pass
        pass

    pass


def check_fitness(network: [], fitness: int):
    with open("./dublin_graph.dat") as f:
        dub_lines = f.readlines()
        dub_lines = dub_lines[2:]
        pass

    net_fitness = hammy_distance(network, dub_lines)
    print("Fitness found was " + str(net_fitness) + " and it should be " + str(fitness))
    if net_fitness == fitness:
        return True
    else:
        return False


def hammy_distance(g1_lines: [], g2_lines: []):
    cost = 0
    edge_list1, adj_matrix1 = edge_list(g1_lines, 200)
    edge_list2, adj_matrix2 = edge_list(g2_lines, 200)

    # print(len(adj_matrix1))
    # print(len(adj_matrix1[0]))
    # print(len(adj_matrix2))
    # print(len(adj_matrix2[0]))

    for row in range(200):
        for col in range(row + 1, 200):
            count1 = adj_matrix1[row][col]
            count2 = adj_matrix2[row][col]
            if count1 != count2:
                if count1 == 0:
                    cost += 5
                elif count2 == 0:
                    cost += 5
                else:
                    cost += abs(count1 - count2)
                    pass
                pass
            pass
        pass
    return cost


def main():
    print("START")
    folder_names = os.listdir(inp)
    mode_itms = [["ED"], ["PM", "P1"], ["PM", "P7"], ["PM", "DUB"]]
    mode_info = ["ED", "PM1", "PM7", "DUBPM"]
    sizes = [128, 256]
    mode_dirs = [[[] for _ in range(len(sizes))] for _ in range(len(mode_itms))]
    muts = [1, 2, 3]
    states = [12, 16, 20]

    # make_from_dirs(folder_names, ["ED", "PM"], ["Epidemic Length Tests", "Profile 1 Matching Tests"])
    # make_nets_in_fold(james_net_inp)

    with open("./dublin_graph.dat") as f:
        dub_lines = f.readlines()
        dub_lines = dub_lines[2:]
        pass
    # make_graph(edge_list(dub_lines, 200), [], [], "Dublin Network", 200, 0)
    # make_graph(make_ring(128), [], [], "128 Ring Network", 128, 0)
    make_ring(128)

    exp_lbls = ["Base1", "Base2"]
    exp_dat = []
    exp_descriptions = ["Edit PS1", "Edit PS2"]
    exp_idx = 1
    for s in states:
        for m in muts:
            exp_dat.append([str(s) + "S", str(m) + "M"])
            exp_lbls.append(str(exp_idx) + "(" + str(s) + ", " + str(m) + ")")
            exp_descriptions.append(str(s) + "States, " + str(m) + "Muts")
            exp_idx += 1
            pass
        pass

    for midx, itms in enumerate(mode_itms):
        for sidx, size in enumerate(sizes):
            for dat in exp_dat:
                for fld in folder_names:
                    if all(itm in fld for itm in itms) and str(size) in fld and all(d in fld for d in dat):
                        mode_dirs[midx][sidx].append(fld)
                        pass
                    elif all(itm in fld for itm in itms) and all(str(si) not in fld for si in sizes) and \
                            all(d in fld for d in dat):
                        mode_dirs[midx][sidx].append(fld)
                        pass
                    pass
                pass
            pass
        pass

    # mode_data[mode][size][exp][run][0] = run num
    # mode_data[mode][size][exp][run][1] = run's fit (sorted based on this)
    # mode_data[mode][size][exp][run][2][:] = run's SDA
    # mode_data[mode][size][exp][run][3][:] = run's network
    # mode_data[mode][size][exp][run][4] = run's network's edges
    # mode_data[mode][size][exp][run][5] = run's network's weight
    # mode_data[mode][size][exp][run][6][:] = run's network's weight hist
    # mode_data[mode][size][exp][run][7] = run's network's best p0
    mode_data = [[[] for _ in range(len(sizes))] for _ in range(len(mode_itms))]
    for midx, itms in enumerate(mode_itms):
        for sidx, size in enumerate(sizes):
            for fld in mode_dirs[midx][sidx]:
                mode_data[midx][sidx].append(get_data(inp + fld + "/"))
                pass
            pass
        pass

    # base_stats[mode][size][exp] = [run's fitness vals]
    base_stats = get_james_data()

    # mode_stats[mode][size][exp] = [run's fitness vals]
    mode_stats = [[[] for _ in range(len(sizes))] for _ in range(len(mode_itms))]
    make_all = False
    make_any = False
    for midx, itms in enumerate(mode_itms):
        for sidx, size in enumerate(sizes):
            if midx < 3:
                mode_path = str(size) + mode_info[midx] + "_Networks/"
                if make_any and not os.path.exists(outp + mode_path):
                    os.makedirs(outp + mode_path)
                    pass
                mode_stats[midx][sidx].append(base_stats[midx][sidx][0])
                mode_stats[midx][sidx].append(base_stats[midx][sidx][1])
            elif midx == 3 and sidx == 0:
                mode_path = mode_info[midx] + "_Networks/"
                if make_any and not os.path.exists(outp + mode_path):
                    os.makedirs(outp + mode_path)
                    pass
                mode_stats[midx][sidx].append(base_stats[midx][sidx][0])
                mode_stats[midx][sidx].append(base_stats[midx][sidx][1])
                pass
            for expidx, exp in enumerate(mode_data[midx][sidx]):
                exp_fits = []
                for runidx, run in enumerate(exp):
                    exp_fits.append(run[1])
                    if make_any:
                        if runidx == 0 or make_all:
                            if midx < 3:
                                make_graph(edge_list(run[3], size), [], [], mode_path +
                                           "EXP" + str(expidx + 1) + "Run" + str(run[0]), size, exp[0][7])
                            elif midx == 3:  # Dublin PM
                                size = 200
                                make_graph(edge_list(run[3], size), [], [], mode_path +
                                           "EXP" + str(expidx + 1) + "Run" + str(run[0]), size, exp[0][7])
                            pass
                        pass
                    pass
                mode_stats[midx][sidx].append(exp_fits)
                pass
            pass
        pass

    titles = ["Epidemic Length", "Epidemic Profile Matching P1", "Epidemic Profile Matching P7",
              "Epidemic Profile Matching Dublin"]
    names = ["EL_boxplot", "PM1_boxplot", "PM7_boxplot", "PMDUB_boxplot"]
    # xsp = [[i for i in range(len(all_data[0]))], [i for i in range(len(all_data[1]))]]
    # xpos = [xsp[0], xsp[1], xsp[0], xsp[1], xsp[0], xsp[1], xsp[0], xsp[1]]
    ylb = ["Fitness", "Fitness", "Fitness", "Fitness", "Fitness"]
    xlb = ["Experiment (Num. States, Max. Mutations)",
           "Experiment (Num. States, Max. Mutations)",
           "Experiment (Num. States, Max. Mutations)",
           "Experiment (Num. States, Max. Mutations)",
           "Experiment (Num. States, Max. Mutations)"]

    # lxpos = []
    # for i in range(2, len(all_data[0]) - 3, 3):
    #     lxpos.append(i + 0.5)
    #     pass
    base_colors = ["#808080", "#696969"]
    colors = ['#0000FF', '#00FF00', '#FFFF00']
    for idx in range(len(titles)):
        for sidx, size in enumerate(sizes):
            if idx < 3 or sidx == 0:
                plt.rc('xtick', labelsize=6)
                plt.rc('ytick', labelsize=6)

                f = plt.figure()
                f.set_figheight(5)
                f.set_figwidth(8)
                plot = f.add_subplot(111)

                bp = plot.boxplot(mode_stats[idx][sidx], patch_artist=True)
                box_plot(bp, 2, [[[0, 1], base_colors], [[i for i in range(2, 11)], colors]])

                # plot.set_xticks(xpos[idx])
                plot.set_xticklabels(exp_lbls, rotation=90)

                if not titles[idx].__contains__("Dublin"):
                    f.suptitle(titles[idx] + " w " + str(size) + " Nodes", fontsize=12)
                else:
                    f.suptitle(titles[idx], fontsize=12)
                    pass
                plot.set_xlabel(xlb[idx], fontsize=10)
                plot.set_ylabel(ylb[idx], fontsize=10)
                # for x in lxpos:
                #     plot.axvline(x=x, color='black', linestyle='--', linewidth=0.75)
                #     pass
                plot.grid(visible="True", axis="y", which='major', color="darkgray", linewidth=0.75)
                f.tight_layout()
                if not titles[idx].__contains__("Dublin"):
                    f.savefig(outp + str(size) + names[idx] + ".png", dpi=300)
                else:
                    f.savefig(outp + names[idx] + ".png", dpi=300)
                    pass
                plt.close()
                pass
            pass
        pass

    # for midx, mode in enumerate(mode_itms):
    #     for sidx, size in enumerate(sizes):
    #         if midx != 3:
    #             make_table(mode_stats[midx][sidx], exp_descriptions, outp + str(size) + names[midx] + ".dat")
    #         elif sidx == 0:
    #             make_table(mode_stats[midx][sidx], exp_descriptions, outp + names[midx] + ".dat")
    #         pass
    #     pass

    # mode_data[mode][size][exp][run][0] = run num
    # mode_data[mode][size][exp][run][1] = run's fit (sorted based on this)
    # mode_data[mode][size][exp][run][2][:] = run's SDA
    # mode_data[mode][size][exp][run][3][:] = run's network
    # mode_data[mode][size][exp][run][4] = run's network's edges
    # mode_data[mode][size][exp][run][5] = run's network's weight
    # mode_data[mode][size][exp][run][6][:] = run's network's weight hist
    # mode_data[mode][size][exp][run][7] = run's network's best p0

    net_stats = [[[[] for _ in range(2)] for _ in range(len(sizes))] for _ in range(len(mode_itms))]
    net_bests = [[[[] for _ in range(2)] for _ in range(len(sizes))] for _ in range(len(mode_itms))]
    for midx, itms in enumerate(mode_itms):
        for sidx, size in enumerate(sizes):
            for expidx, exp in enumerate(mode_data[midx][sidx]):
                exp_edges = []
                exp_weights = []
                for runidx, run in enumerate(exp):
                    if runidx == 0:
                        net_bests[midx][sidx][0].append(run[4])
                        net_bests[midx][sidx][1].append(run[5])
                        pass
                    exp_edges.append(run[4])
                    exp_weights.append(run[5])
                    pass
                net_stats[midx][sidx][0].append(exp_edges)
                net_stats[midx][sidx][1].append(exp_weights)
                pass
            pass
        pass

    titles = ["Epidemic Length", "Epidemic Profile Matching P1", "Epidemic Profile Matching P7",
              "Epidemic Profile Matching Dublin"]
    names = ["EL_netstats", "PM1_netstats", "PM7_netstats", "PMDUB_netstats"]
    ylb = ["Network Edges", "Network Weight"]
    xlb = ["Experiment (Num. States, Max. Mutations)",
           "Experiment (Num. States, Max. Mutations)",
           "Experiment (Num. States, Max. Mutations)",
           "Experiment (Num. States, Max. Mutations)",]
    out_path = outp + "Network Stats Boxplots"
    if not os.path.exists(out_path):
        os.makedirs(out_path)
        pass
    out_path += "/"
    xsp = [[i for i in range(1, 10)], [i - 0.22 for i in range(1, 10)], [i + 0.22 for i in range(1, 10)]]
    for idx in range(len(titles)):
        for sidx, size in enumerate(sizes):
            if idx < 3 or sidx == 0:
                if idx >= 3:
                    size = 200
                    pass
                plt.rc('xtick', labelsize=6)
                plt.rc('ytick', labelsize=6)

                f = plt.figure()
                f.set_figheight(5)
                f.set_figwidth(8)
                plot = f.add_subplot(111)
                plot2 = plot.twinx()

                # plot2.set_axisbelow(True)
                # plot2.grid(visible="True", axis="y", which='major', color="lightgray", linewidth=0.75)

                bp = plot.boxplot(net_stats[idx][sidx][0], positions=xsp[1], patch_artist=True, zorder=1, widths=0.4)
                box_plot(bp, 1, [[[i for i in range(9)], ["#0000FF"]]])
                plot.plot(xsp[1], net_bests[idx][sidx][0], 'x', color="#FF0000")
                bp = plot2.boxplot(net_stats[idx][sidx][1], positions=xsp[2], patch_artist=True, zorder=2, widths=0.4)
                box_plot(bp, 1, [[[i for i in range(9)], ["#00FF00"]]])
                plot2.plot(xsp[2], net_bests[idx][sidx][1], 'x', color="#FF0000")

                plot.hlines([size * 1, size * 4], 0.5, 9.5, colors="#0000FF", linestyles="dashed", linewidth=1)
                plot2.hlines([size * 4, size * 16], 0.5, 9.5, colors="#00FF00", linestyles="dotted", linewidth=1)

                plot.set_xticks(xsp[0])
                plot.set_xticklabels(exp_lbls[2:], rotation=90)

                if not titles[idx].__contains__("Dublin"):
                    f.suptitle(titles[idx] + " w " + str(size) + " Nodes", fontsize=12)
                else:
                    f.suptitle(titles[idx], fontsize=12)
                    pass

                plot.set_xlabel(xlb[idx], fontsize=10)
                plot.set_ylabel(ylb[0], fontsize=10, color="#0000FF")
                plot2.set_ylabel(ylb[1], fontsize=10, rotation=270, labelpad=10, color="#00FF00")

                plot.set_axisbelow(True)
                plot.grid(visible="True", axis="y", which='major', color="darkgray", linewidth=0.75)
                f.tight_layout()

                if idx < 3:
                    f.savefig(out_path + str(size) + names[idx] + ".png", dpi=450)
                else:
                    f.savefig(out_path + names[idx] + ".png", dpi=450)
                    pass
                plt.close()
            pass
        pass
    print("END")
    pass


main()
