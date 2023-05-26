#include "main.h"

int main(int argc, char *argv[]) {
    fstream endStats, runStats, readme;
    char fn[70];
    char outDir[50];

    getArgs(argv);
    initAlg();
    cmdLineIntro(cout);
    dublinWeights = dublin_graph.fill(netRoot);
    verts = Graph::quadForm(1, -1, (-1) * (int) dublinWeights.size() * 2);

    sprintf(outDir, "%sNetMatch%03d w %02dP, %02dS, %dM/", outRoot, verts, penalty, sdaStates, maxMuts);
    filesystem::create_directory(outDir);

    for (int runIdx = runNum; runIdx < runNum + runs; ++runIdx) {
        sprintf(fn, "%sbest%02d.dat", outDir, runIdx + 1);
        endStats.open(fn, ios::out);

        sprintf(fn, "%sreadme.dat", outDir);
        readme.open(fn, ios::out);
        createReadMe(readme);
        readme.close();

        sprintf(fn, "%srun%02d.dat", outDir, runIdx + 1);
        runStats.open(fn, ios::out);
        initpop();
        if (verbose) cmdLineRun(runIdx, cout);
        runStats << left << setw(4) << 0;
        report(runStats);

        int mev = 0;
        double lastBest = MAXFLOAT;
        double curBest;
        bool extraLife = true;
        int lastBestRpt = 0;
        do {
            if (mev != 0) {
                if (mev > generations / 2){
                    extraLife = false;
                }
                lastBestRpt = 0;
            }
            while (mev < generations && lastBestRpt < 15) {
                matingevent();
                if ((mev + 1) % RE == 0) {
                    if (verbose) {
                        cout << left << setw(5) << runIdx + 1;
                        cout << left << setw(4) << (mev + 1) / RE;
                    }
                    runStats << left << setw(4) << (mev + 1) / RE;
                    curBest = report(runStats);
                    if (lastBest == curBest) {
                        lastBestRpt++;
                    } else {
                        lastBestRpt = 0;
                        lastBest = curBest;
                    }
                }
                mev++;
            }
            if (mev < generations && extraLife) {
                doCulling(0.99, false, false);
            }
        } while (extraLife);


//        for (int mev = 0; mev < generations; mev++) {
//            matingevent();
//            if ((mev + 1) % RE == 0) {
//                if (verbose) {
//                    cout << left << setw(5) << runIdx + 1;
//                    cout << left << setw(4) << (mev + 1) / RE;
//                }
//                runStats << left << setw(4) << (mev + 1) / RE;
//                report(runStats);
//            }
//        }

        reportbest(endStats);
        endStats.close();
        runStats.close();
    }

    delete[] bPop;
    delete[] fit;
    delete[] dead;
    return (0); // keep the system happy
}

int doCulling(double percentage, bool rndPick, bool biggerBetter) {
    int numKillings = (int) (popsize * percentage);
    vector<int> contestants;
    vector<int> winners;
    winners.reserve(numKillings);

    if (rndPick) {
        contestants = tournSelect(numKillings, !biggerBetter);
    } else {
        contestants = tournSelect(popsize, !biggerBetter);
    }

    for (int idx = 0; idx < numKillings; idx++) {
        winners.push_back(contestants[idx]);
    }

    double sum1 = 0, sum2 = 0;
    for (int idx: winners) {
        bPop[idx].randomize();
        sum1 += fit[idx];
        fit[idx] = fitness(bPop[idx]);
        sum2 += fit[idx];
    }

    cout<< "CULLING: " << sum1/numKillings << " -> " << sum2/numKillings << endl;

    return 0;
}

bool necroticFilter(Bitsprayer &A) { // True means dead.
    return false;
    int size = verts * (verts - 1) / 2;
    vector<int> vals;
    vals.reserve(size);

    vals = A.getBitsVec(size);

    int totWeight = 0;
    int totEdges = 0;
    for (int v: vals) {
        totWeight += v;
        if (v > 0) {
            totEdges++;
        }
    }
    return false;

    if (totEdges < edgBnds[0] * verts || totEdges > edgBnds[1] * verts) {
        return true;
    }
    if (totWeight < wghtBnds[0] * verts || totWeight > wghtBnds[1] * verts) {
        return true;
    }
    return false;
}

double fitness(Bitsprayer &A) {
    int size = verts * (verts - 1) / 2;
    vector<int> vals;
    vals.reserve(size);

    vals = A.getBitsVec(size);
//    Graph G(verts);
//    G.fill(vals, diagFill);
//
//    double fi = G.hammy_distance(dublin_graph, penalty);
    double fi = quickHammy(vals, dublinWeights, penalty);
    return fi;
}

bool compareFitness(int popIdx1, int popIdx2) {
    return fit[popIdx1] < fit[popIdx2];
}

vector<int> tournSelect(int size, bool decreasing) {
    vector<int> tournIdxs;
    int idxToAdd;

    tournIdxs.reserve(size);
    if (size == popsize) {
        for (int idx = 0; idx < size; idx++) {
            tournIdxs.push_back(idx);
        }
    } else {
        do {
            idxToAdd = (int) lrand48() % popsize;
            if (count(tournIdxs.begin(), tournIdxs.end(), idxToAdd) == 0) {
                tournIdxs.push_back(idxToAdd);
            }
        } while (tournIdxs.size() < tsize);
    }

    sort(tournIdxs.begin(), tournIdxs.end(), compareFitness);
    if (decreasing) {
        reverse(tournIdxs.begin(), tournIdxs.end());
    }
    return tournIdxs;
}

void matingevent() { // runNum a mating event
    int rp;
//    Bitsprayer child1, child2;
//    child1 = Bitsprayer(sdaStates, zeroProb);
//    child2 = Bitsprayer(sdaStates, zeroProb);
    vector<int> tournIdxs;

//    bool firstDead, secondDead;
//    int deaths;
//    do {
    tournIdxs = tournSelect(tsize, true);

    bPop[tournIdxs[0]].copy(bPop[tournIdxs[tsize - 2]]);
    bPop[tournIdxs[1]].copy(bPop[tournIdxs[tsize - 1]]);
    bPop[tournIdxs[0]].twoPtCrossover(bPop[tournIdxs[1]]);
    rp = (int) lrand48() % maxMuts + 1;
    bPop[tournIdxs[0]].mutate(rp);
    rp = (int) lrand48() % maxMuts + 1;
    bPop[tournIdxs[1]].mutate(rp);

//        firstDead = necroticFilter(child1);
//        secondDead = necroticFilter(child2);
//        deaths = 0;
//        if (firstDead || secondDead) {
//            for (int i = 0; i < popsize; ++i) {
//                if (i != tournIdxs[tsize - 1] && i != tournIdxs[tsize - 2] && dead[i]) {
//                    deaths++;
//                }
//            }
//        }
//        if (firstDead) {
//            deaths++;
//        }
//        if (secondDead) {
//            deaths++;
//        }
//    } while (false);
//    } while (deaths > 0.8 * popsize);

//    bPop[tournIdxs[0]].copy(child1);
//    bPop[tournIdxs[1]].copy(child2);
//    dead[tournIdxs[0]] = false;
//    dead[tournIdxs[1]] = false;
//    double fitVal = 0.0;

    fit[tournIdxs[0]] = fitness(bPop[tournIdxs[0]]);
    fit[tournIdxs[1]] = fitness(bPop[tournIdxs[1]]);

//    if (firstDead) {
//        dead[tournIdxs[0]] = true;
//        fit[tournIdxs[0]] = worst_fit;
//    } else {
//        fitVal = fitness(bPop[tournIdxs[0]]);
//        if (fitVal > 1.1 * worst_fit) {
//            cout << fitVal << "\t" << 1.1 * fitVal << endl;
//            worst_fit = 1.1 * fitVal;
//            for (int i = 0; i < popsize; i++) {
//                if (dead[i]) {
//                    fit[i] = worst_fit;
//                }
//            }
//        } else {
//            fit[tournIdxs[0]] = fitVal;
//        }
//    }
//
//    if (secondDead) {
//        dead[tournIdxs[1]] = true;
//        fit[tournIdxs[1]] = worst_fit;
//    } else {
//        fitVal = fitness(bPop[tournIdxs[1]]);
//        if (fitVal > 1.1 * worst_fit) {
//            cout << fitVal << "\t" << 1.1 * fitVal << endl;
//            worst_fit = 1.1 * fitVal;
//            for (int i = 0; i < popsize; i++) {
//                if (dead[i]) {
//                    fit[i] = worst_fit;
//                }
//            }
//        } else {
//            fit[tournIdxs[1]] = fitVal;
//        }
//    }
}

vector<double> calcStats(const vector<int> &goodIdxs, bool biggerBetter) {
    double sum = 0.0;
    double bestVal = (biggerBetter ? 0.0 : MAXFLOAT);
    double worstVal = (biggerBetter ? MAXFLOAT : 0.0);

    for (int idx: goodIdxs) {
        sum += fit[idx];
        if ((biggerBetter && fit[idx] > bestVal) || (!biggerBetter && fit[idx] < bestVal)) {
            bestVal = fit[idx];
        }
        if ((biggerBetter && fit[idx] < worstVal) || (!biggerBetter && fit[idx] > worstVal)) {
            worstVal = fit[idx];
        }
    }

    double mean = sum / (double) goodIdxs.size();
    double stdDevSum = 0.0;
    for (int idx: goodIdxs) {
        stdDevSum += pow(fit[idx] - mean, 2);
    }
    double stdDev = sqrt(stdDevSum / ((double) goodIdxs.size() - 1.0));
    double CI95 = 1.96 * (stdDev / sqrt((double) goodIdxs.size()));

    return {mean, stdDev, CI95, bestVal, worstVal}; // {mean, stdDev, 95CI, best, worst}
}

double report(ostream &aus) { // make a statistical report
    vector<int> goodIdxs;
    int deaths = 0;
//    int cnt = 0;

    for (int i = 0; i < popsize; i++) {
        if (!dead[i]) {
            goodIdxs.push_back(i);
        } else {
            deaths++;
        }
    }

    vector<double> stats = calcStats(goodIdxs, false);
    double mean = stats[0];
    double stdDev = stats[1];
    double CI95 = stats[2];
    double bestVal = stats[3];
    double worstVal = stats[4];

//    worst_fit = 0;
//    if (worstVal > worst_fit) {
//        worst_fit = worstVal * 1.1;
//    }
//    for (int i = 0; i < popsize; i++) {
//        if (dead[i]) {
//            fit[i] = worst_fit;
//        }
//    }

    aus << left << setw(10) << mean;
    aus << left << setw(12) << CI95;
    aus << left << setw(10) << stdDev;
    aus << left << setw(8) << bestVal << endl;
//    aus << "Dead: " << deaths << "\t" << worst_fit << endl;
    if (verbose) {
        cout << left << setw(10) << mean;
        cout << left << setw(12) << CI95;
        cout << left << setw(10) << stdDev;
        cout << left << setw(8) << bestVal << endl;
//        cout << "Dead: " << deaths << "\t" << worst_fit << endl;
    }

//    int b = -1;
//    double best_fit = MAXFLOAT;
//    for (int idx: goodIdxs) {
//        if (fit[idx] < best_fit && !dead[idx]) {
//            best_fit = fit[idx];
//            b = idx;
//        }
//    }

//    if (best_fit != bestVal || bestVal != fitness(bPop[b])) {
//        cout << endl << "ERROR!! The bests do not align!" << endl;
//        aus << endl << "ERROR!! The bests do not align!" << endl;
//        cout << "Best Saved: " << best_fit << "; Calculated: " << fitness(bPop[b]) << "; From Stats: " << bestVal
//             << "\t";
//        aus << "Best Saved: " << best_fit << "; Calculated: " << fitness(bPop[b]) << "; From Stats: " << bestVal
//            << "\t";
//    }
    return bestVal;
}

void reportbest(ostream &aus) {
    int b;
    Graph G(verts);

    b = -1;
    double best_fit = MAXFLOAT;
    for (int i = 0; i < popsize; i++) {
        if (fit[i] < best_fit && !dead[i]) {
            b = i;
            best_fit = fit[i];
        }
    }

    double fitCheck = fitness(bPop[b]);
    if (fit[b] != fitCheck) {
        cout << "ERROR!!!  Fitness not what is expected!" << endl;
        aus << "ERROR!!!  Fitness not what is expected!" << endl;
    }
    if (dead[b]) {
        cout << "ERROR!!!  Best is dead!" << endl;
        aus << "ERROR!!!  Best is dead!" << endl;
    }
    if (necroticFilter(bPop[b])) {
        cout << "ERROR!!!  Best fails necrotic filter!!" << endl;
        aus << "ERROR!!!  Best fails necrotic filter!!" << endl;
    }
    if (fit[b] != fitCheck) {
        cout << "ERROR!!!  Fitness values do not align!" << endl;
        aus << "ERROR!!!  Fitness values do not align!" << endl;
        cout << "Stored Fitness: " << fit[b] << "; Calculated Fitness: " << fitCheck << endl;
        aus << "Stored Fitness: " << fit[b] << "; Calculated Fitness: " << fitCheck << endl;
    }

    aus << fit[b] << " -fitness" << endl;
    // Write the SDA
    aus << "Self-Driving Automata" << endl;
    bPop[b].print(aus);
    int size = verts * (verts - 1) / 2;
    bPop[b].printBitsVec(size, aus);
    vector<int> vals = bPop[b].getBitsVec(size);
    G.fill(vals, diagFill);

    aus << "Graph" << endl;
    G.print(aus);
    aus << endl;
}
