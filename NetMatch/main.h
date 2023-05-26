#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iomanip>
#include <string>

// This block of code includes the filesystem package (which depends on OS)
#if defined(__cplusplus) && __cplusplus >= 201703L && defined(__has_include)
#if __has_include(<filesystem>)
#define GHC_USE_STD_FS

#include <filesystem>
#include <thread>

namespace filesystem = std::filesystem;
#endif
#endif
#ifndef GHC_USE_STD_FS
#include "filesystem.hpp"
namespace filesystem = ghc::filesystem;
#endif

using namespace std;

#include "Bitsprayer.h"
#include "Graph.h"

/*************************algorithm controls******************************/
#define verbose true
#define tsize 5
int runNum;
int runs;
int seed;
int verts;
int sdaStates;
int maxMuts;
int generations;
int popsize;
int penalty;
int network;

#define RIs 100
#define RE (int)5000

bool *dead;
double worst_fit;

static vector<double> edgBnds = {1, 10};
static vector<int> wghtBnds = {4, 20};
static double zeroProb = 0.90;
static bool diagFill = true;

/**************************Variable dictionary************************/
Bitsprayer *bPop;
double *fit;
Graph dublin_graph;
vector<int> dublinWeights;

/****************************Procedures********************************/
void initpop();
void matingevent();
double report(ostream &aus);
void reportbest(ostream &aus);
void createReadMe(ostream &aus);
void cmdLineIntro(ostream &aus);
void cmdLineRun(int run, ostream &aus);
double fitness(Bitsprayer &A);
int getArgs(char *args[]);
void initAlg();
bool necroticFilter(Bitsprayer &A);
vector<int> tournSelect(int size, bool decreasing);
vector<double> calcStats(const vector<int> &goodIdxs, bool biggerBetter);
int doCulling(double percentage, bool rndPick, bool biggerBetter);
int quickHammy(const vector<int>& seqOne, const vector<int>& seqTwo, int pen);

/****************************Main Routine*******************************/

char outRoot[20] = "./NetMatchOut/";
char netRoot[40];

/**
 * System Parameters:
 * 1. Network (1 -> 200, 2 -> 80)
 * 2. Population size
 * 3. Generations
 * 4. Penalty
 * 5. Number of states
 * 6. Maximum mutations
 * 7. Runs to run
 * 8. First run index
 * 9. Seed
 */

void initpop() {
    int count = 0;
    for (int i = 0; i < popsize; i++) {
//        dead[i] = false;
//        do {
            bPop[i].randomize();
            count++;
//        } while (necroticFilter(bPop[i]));
        fit[i] = fitness(bPop[i]);
    }
    cout << "\n" << "Attempts: " << count << endl;
}

void createReadMe(ostream &aus) {
    aus << "This file contains the info about the files in this folder and the parameter settings for this experiment." << endl;
    aus << "Graph Evolution Tool." << endl;
    aus << "Fitness Function: Network Matching" << endl;
    aus << "Hamming Distance Penalty: " << penalty << endl;
    aus << "Network: " << netRoot << endl;
    aus << "Representation: Weighted SDAs" << endl;
    aus << endl;
    aus << "The parameter settings are as follows: " << endl;
    aus << "Mating events: " << generations << endl;
    aus << "Population size: " << popsize << endl;
    aus << "Number of vertices: " << verts << endl;
    aus << "Maximum number of mutations: " << maxMuts << endl;
    aus << "Tournament size: " << tsize << endl;
    aus << "Number of States: " << sdaStates << endl;
    aus << "Zero probability for initial SDAs: " << zeroProb << endl;
    aus << "Adjacency matrix was filled: " << (diagFill ? "diagonally" : "row by row") << endl;
    aus << endl;
    aus << "The file descriptions are as follows: " << endl;
    aus << "best##.lint -> the best fitness and it's associated data for each runNum";
    aus << endl;
    aus << "runNum##.dat -> population statistics for each runNum" << endl;
}

void cmdLineIntro(ostream &aus) {
    aus << "Graph Evolution Tool." << endl;
    aus << "Fitness Function: Network Matching" << endl;
    aus << "Hamming Distance Penalty: " << penalty << endl;
    aus << "Network: " << netRoot << endl;
    aus << "Representation: Weighted SDAs" << endl;
    aus << "Check readme.dat for more information about parameters/output.";
    aus << endl;
}

void cmdLineRun(int run, ostream &aus) {
    aus << endl << "Beginning Run " << run + 1 << " of " << runs << endl;
    aus << left << setw(5) << "Run";
    aus << left << setw(4) << "RI";
    aus << left << setw(10) << "Mean";
    aus << left << setw(12) << "95% CI";
    aus << left << setw(10) << "SD";
    aus << left << setw(8) << "Best";
    aus << endl;
    aus << left << setw(5) << run + 1;
    aus << left << setw(4) << "0";
}

void initAlg() { // initialize the algorithm
    srand48(seed); // read the random number seed
    bPop = new Bitsprayer[popsize];
    for (int i = 0; i < popsize; ++i) {
        bPop[i] = Bitsprayer(sdaStates, zeroProb);
    }
    fit = new double[popsize];
    dead = new bool[popsize];
}

int getArgs(char *args[]) {
    string arg = args[1]; // network
    try {
        size_t pos;
        network = std::stoi(arg, &pos);
        if (pos < arg.size()) {
            std::cerr << "Trailing characters after number: " << arg << endl;
        }
    } catch (std::invalid_argument const &ex) {
        std::cerr << "Invalid number: " << arg << '\n';
    } catch (std::out_of_range const &ex) {
        std::cerr << "Number out of range: " << arg << '\n';
    }

    if (network == 0) {
        sprintf(netRoot, "./NetMatchIn/test_graph.dat");
    } else if (network == 1) {
        sprintf(netRoot, "./NetMatchIn/dublin_graph.dat");
    } else if (network == 2) {
        sprintf(netRoot, "./NetMatchIn/dublin_graph80.dat");
    } else {
        cout << "ERROR!! Could not fill variable netRoot with path to network." << endl;
    }

    arg = args[2]; // popsize
    try {
        size_t pos;
        popsize = std::stoi(arg, &pos);
        if (pos < arg.size()) {
            std::cerr << "Trailing characters after number: " << arg << endl;
        }
    } catch (std::invalid_argument const &ex) {
        std::cerr << "Invalid number: " << arg << '\n';
    } catch (std::out_of_range const &ex) {
        std::cerr << "Number out of range: " << arg << '\n';
    }

    arg = args[3]; // generations
    try {
        size_t pos;
        generations = std::stoi(arg, &pos);
        if (pos < arg.size()) {
            std::cerr << "Trailing characters after number: " << arg << endl;
        }
    } catch (std::invalid_argument const &ex) {
        std::cerr << "Invalid number: " << arg << '\n';
    } catch (std::out_of_range const &ex) {
        std::cerr << "Number out of range: " << arg << '\n';
    }

    arg = args[4]; // penalty
    try {
        size_t pos;
        penalty = std::stoi(arg, &pos);
        if (pos < arg.size()) {
            std::cerr << "Trailing characters after number: " << arg << endl;
        }
    } catch (std::invalid_argument const &ex) {
        std::cerr << "Invalid number: " << arg << '\n';
    } catch (std::out_of_range const &ex) {
        std::cerr << "Number out of range: " << arg << '\n';
    }

    arg = args[5]; // sdaStates
    try {
        size_t pos;
        sdaStates = std::stoi(arg, &pos);
        if (pos < arg.size()) {
            std::cerr << "Trailing characters after number: " << arg << endl;
        }
    } catch (std::invalid_argument const &ex) {
        std::cerr << "Invalid number: " << arg << '\n';
    } catch (std::out_of_range const &ex) {
        std::cerr << "Number out of range: " << arg << '\n';
    }

    arg = args[6]; // maxMuts
    try {
        size_t pos;
        maxMuts = std::stoi(arg, &pos);
        if (pos < arg.size()) {
            std::cerr << "Trailing characters after number: " << arg << endl;
        }
    } catch (std::invalid_argument const &ex) {
        std::cerr << "Invalid number: " << arg << '\n';
    } catch (std::out_of_range const &ex) {
        std::cerr << "Number out of range: " << arg << '\n';
    }

    arg = args[7]; // runs
    try {
        size_t pos;
        runs = std::stoi(arg, &pos);
        if (pos < arg.size()) {
            std::cerr << "Trailing characters after number: " << arg << endl;
        }
    } catch (std::invalid_argument const &ex) {
        std::cerr << "Invalid number: " << arg << '\n';
    } catch (std::out_of_range const &ex) {
        std::cerr << "Number out of range: " << arg << '\n';
    }

    arg = args[8]; // runNum
    try {
        size_t pos;
        runNum = std::stoi(arg, &pos);
        if (pos < arg.size()) {
            std::cerr << "Trailing characters after number: " << arg << endl;
        }
    } catch (std::invalid_argument const &ex) {
        std::cerr << "Invalid number: " << arg << '\n';
    } catch (std::out_of_range const &ex) {
        std::cerr << "Number out of range: " << arg << '\n';
    }

    arg = args[9]; // seed
    try {
        size_t pos;
        seed = std::stoi(arg, &pos);
        if (pos < arg.size()) {
            std::cerr << "Trailing characters after number: " << arg << endl;
        }
    } catch (std::invalid_argument const &ex) {
        std::cerr << "Invalid number: " << arg << '\n';
    } catch (std::out_of_range const &ex) {
        std::cerr << "Number out of range: " << arg << '\n';
    }

    cout << "Arguments Captured!" << endl;
    return 0;
}

int quickHammy(const vector<int>& seqOne, const vector<int>& seqTwo, int pen) {
    int cost = 0;
    int diff;
    for (int idx = 0; idx < seqOne.size(); idx++) {
        diff = abs(seqOne[idx] - seqTwo[idx]);
        if (diff > 0) {
            if (seqOne[idx] == 0 || seqTwo[idx] == 0) {
                cost += pen;
            } else {
                cost += min(diff, pen);
            }
        }
    }
    return cost;
}