def main():
    fname = "./michaelGraph.dat"
    numNodes = 128
    adjM = [[0 for _ in range(numNodes)] for _ in range(numNodes)]
    with open(fname) as f:
        lines = f.readlines()
        lines = lines[5:]
        for from_node, line in enumerate(lines):
            line = line.rstrip()
            line = line.split(" ")
            if len(line) > 0 and line[0] != '':
                for to_node in line:
                    weight = to_node.split("[")[1].strip("]")
                    to_node = to_node.split("[")[0]
                    adjM[from_node][int(to_node)] = weight
                    adjM[int(to_node)][from_node] = weight
                    pass
                pass
            pass
        pass

    with open("michaelGraph_forme.dat", "w") as f:
        for rowIdx in range(numNodes):
            for colIdx in range(rowIdx + 1, numNodes):
                    f.write(str(adjM[rowIdx][colIdx]) + "\n")
                    pass
            pass
        pass
    pass


main()