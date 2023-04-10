def main():
    fname = "./JamesBestED.dat"
    numNodes = 128
    weights = 0
    edges = 0
    adjM = [[0 for _ in range(numNodes)] for _ in range(numNodes)]
    with open(fname) as f:
        lines = f.readlines()
        for from_node, line in enumerate(lines):
            line = line.rstrip()
            line = line.split(" ")
            for to_node in line:
                if to_node != '':
                    adjM[from_node][int(to_node)] += 1
                    if int(to_node) > from_node:
                        weights += 1
                # adjM[int(to_node)][from_node] += 1
                pass
            pass
        pass
    print(weights)
    w_cnt = [0 for _ in range(6)]
    with open("michaelGraph_forme2.dat", "w") as f:
        for rowIdx in range(numNodes):
            for colIdx in range(rowIdx + 1, numNodes):
                f.write(str(adjM[rowIdx][colIdx]) + '\n')
                if adjM[rowIdx][colIdx]>0:
                    edges += 1
                    w_cnt[adjM[rowIdx][colIdx]] += 1
                    pass
                pass
            pass
        pass
    print(edges)
    print(w_cnt)
    pass


main()
