from random import randint


def main():
    path = "./cmake-build-release---remote/NetMatch"
    network = [1,2]
    popsize = 250
    generations = 1000000
    penalty = [1,5,50]
    states = [8, 16, 24]
    muts = [2, 8]
    runs = 1
    totRuns = 30

    for net in network:
        with open("./ComputeCanadaScripts/table" + str(net) + ".dat", "w") as f:
            for pen in penalty:
                for sta in states:
                    for mut in muts:
                        for idx in range(totRuns):
                            f.write(path + " " + str(net) + " " + str(popsize) + " " + str(generations) + " " + str(pen)
                                    + " " + str(sta) + " " + str(mut) + " " + str(runs) + " " + str(idx) + " " +
                                    str(randint(1, 1000)) + "\n")
                            pass
                        pass
                    pass
                pass
            pass
        pass
    pass

main()