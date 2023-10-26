#ifndef CIRCUIT_H
#define CIRCUIT_H

#include <fstream>
#include <climits>
#include <cmath>
#include <string>
#include <vector>
#include <unordered_map>
#include <set>
using namespace std;

struct Cell
{
    string name;
    vector<int> associatedNets;
    int group;
    int sim_group;
    Cell* prev;
    Cell* next;

    Cell() : prev(nullptr), next(nullptr) {};
    Cell(string nm) : name(nm), prev(nullptr), next(nullptr) {};
};

struct Net
{
    string name;
    vector<int> associatedCells;
    int numCellsInGroup[2];
    int sim_numCellsInGroup[2];

    Net() : numCellsInGroup{0, 0} {};
    Net(string nm) : name(nm), numCellsInGroup{0, 0} {};
};

class Circuit
{
public:
    Circuit() {};
    ~Circuit() { freePointers(); };
    void readInput(string ifName);
    void partition();
    void writeOutput(string ofName);

private:
    double balanceFactor;
    vector<Cell*> cellList;
    vector<Net*> netList;
    unordered_map<string, int> name2Idx;

    size_t pMax;
    int sizeOfGroup0;
    vector<Cell> bucketOfGroup[2];  // bucket[0/1][idx] is the list of cells whose gain is (idx - pMax)

    int initialPartition();
    void initializeDataStructures();
    void calculateInitialGain(vector<int>& g);
    void copyCircuitState();
    bool passBalanceCriterion(int size0);
    void moveCellAndUpdateNeighboringGain(int i, vector<bool>& locked, vector<int>& cellGain);
    void updateCellGain(int i, int delta, vector<int>& cellGain);
    void updateCircuitState(const vector<int>& movedCell, int k);
    int cutSize();
    void freePointers();
};

#endif