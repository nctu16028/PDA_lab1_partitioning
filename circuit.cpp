#include "circuit.h"


void Circuit::readInput(string ifName)
{
    ifstream fin(ifName);
    fin >> balanceFactor;

    string token;
    int netIdx;
    while (fin >> token)    // read keyword "NET"
    {
        fin >> token;       // read the net name
        netIdx = netList.size();
        name2Idx.insert({token, netIdx});
        Net* currNet = new Net(token);
        netList.push_back(currNet);

        int cellIdx;
        set<int> cellSet;   // to remove duplicated cells
        while (fin >> token)    // read the cell name or ";"
        {
            if (token == ";")
                break;

            if (name2Idx.find(token) == name2Idx.end())
            {
                cellIdx = cellList.size();
                name2Idx.insert({token, cellIdx});
                Cell* newCell = new Cell(token);
                cellList.push_back(newCell);
            }
            else
                cellIdx = name2Idx[token];
            cellSet.insert(cellIdx);
        }
        netList[netIdx]->associatedCells.assign(cellSet.begin(), cellSet.end());
        for (auto& cellIdx : cellSet)
            cellList[cellIdx]->associatedNets.push_back(netIdx);
    }
    fin.close();
}


void Circuit::partition()
{
    sizeOfGroup0 = initialPartition();

    initializeDataStructures();
    while (true)
    {
        vector<bool> locked(cellList.size(), false);
        vector<int> potentialGainOfCell(cellList.size(), 0);
        calculateInitialGain(potentialGainOfCell);
        int currGain = 0;
        int bestGain = INT_MIN;
        int minUnbalancy = cellList.size();
        int stepToBestGain = 0;
        vector<int> movingOrder;

        copyCircuitState();     // for the following simulation pass
        int sim_sizeOfGroup0 = sizeOfGroup0;
        for (size_t i = 0; i < cellList.size(); i++)
        {
            for (int g = pMax; g + pMax >= 0; g--)
            {
                if (bucketOfGroup[0][g + pMax].next != nullptr && passBalanceCriterion(sim_sizeOfGroup0 - 1))
                {
                    // Remove the cell from the bucket
                    Cell* toMove = bucketOfGroup[0][g + pMax].next;
                    if (toMove->next != nullptr)
                        toMove->next->prev = &(bucketOfGroup[0][g + pMax]);
                    bucketOfGroup[0][g + pMax].next = toMove->next;
                    toMove->next = nullptr;
                    toMove->prev = nullptr;

                    // Move the chosen cell from group0 to group1
                    int cellIdx = name2Idx[toMove->name];
                    locked[cellIdx] = true;
                    moveCellAndUpdateNeighboringGain(cellIdx, locked, potentialGainOfCell);
                    movingOrder.push_back(cellIdx);
                    sim_sizeOfGroup0--;

                    // Record the best gain
                    currGain += g;
                    int unbalancy = abs(sim_sizeOfGroup0 - (int)cellList.size() / 2);
                    if (currGain > bestGain)
                    {
                        bestGain = currGain;
                        minUnbalancy = unbalancy;
                        stepToBestGain = i + 1;
                    }
                    else if (currGain == bestGain && unbalancy < minUnbalancy)
                    {
                        minUnbalancy = unbalancy;
                        stepToBestGain = i + 1;
                    }

                    break;
                }
                else if (bucketOfGroup[1][g + pMax].next != nullptr && passBalanceCriterion(sim_sizeOfGroup0 + 1))
                {
                    // Remove the cell from the bucket
                    Cell* toMove = bucketOfGroup[1][g + pMax].next;
                    if (toMove->next != nullptr)
                        toMove->next->prev = &(bucketOfGroup[1][g + pMax]);
                    bucketOfGroup[1][g + pMax].next = toMove->next;
                    toMove->next = nullptr;
                    toMove->prev = nullptr;

                    // Move the chosen cell from group1 to group0
                    int cellIdx = name2Idx[toMove->name];
                    locked[cellIdx] = true;
                    moveCellAndUpdateNeighboringGain(cellIdx, locked, potentialGainOfCell);
                    movingOrder.push_back(cellIdx);
                    sim_sizeOfGroup0++;

                    // Record the best gain
                    currGain += g;
                    int unbalancy = abs(sim_sizeOfGroup0 - (int)cellList.size() / 2);
                    if (currGain > bestGain)
                    {
                        bestGain = currGain;
                        minUnbalancy = unbalancy;
                        stepToBestGain = i + 1;
                    }
                    else if (currGain == bestGain && unbalancy < minUnbalancy)
                    {
                        minUnbalancy = unbalancy;
                        stepToBestGain = i + 1;
                    }

                    break;
                }
            }
        }

        if (bestGain > 0)
            updateCircuitState(movingOrder, stepToBestGain);
        else
            break;
    }
}


int Circuit::initialPartition()
{
    int size0 = 0;
    for (size_t i = 0; i < cellList.size(); i++)
    {
        size0 += (i % 2 == 0);
        cellList[i]->group = i % 2;
        for (auto& netIdx : cellList[i]->associatedNets)
            netList[netIdx]->numCellsInGroup[i % 2]++;
    }

    return size0;
}


void Circuit::initializeDataStructures()
{
    pMax = 0;
    for (size_t i = 0; i < cellList.size(); i++)
        pMax = max(pMax, cellList[i]->associatedNets.size());

    // Allocate bucket slots
    bucketOfGroup[0].resize(2 * pMax + 1, Cell());
    bucketOfGroup[1].resize(2 * pMax + 1, Cell());
}


void Circuit::calculateInitialGain(vector<int>& g)
{
    for (size_t i = 0; i < cellList.size(); i++)
    {
        // Calculate initial potential gain of Cell i
        int from = cellList[i]->group;
        int to = 1 - from;
        for (auto& netIdx : cellList[i]->associatedNets)
        {
            if (netList[netIdx]->numCellsInGroup[from] == 1)
                g[i]++;
            if (netList[netIdx]->numCellsInGroup[to] == 0)
                g[i]--;
        }

        // Insert Cell i into one of the bucket
        Cell* pTemp = &(bucketOfGroup[cellList[i]->group][g[i] + pMax]);
        if (pTemp->next != nullptr)
            pTemp->next->prev = cellList[i];
        cellList[i]->next = pTemp->next;
        cellList[i]->prev = pTemp;
        pTemp->next = cellList[i];
    }
}


void Circuit::copyCircuitState()
{
    for (size_t i = 0; i < cellList.size(); i++)
        cellList[i]->sim_group = cellList[i]->group;
    for (size_t i = 0; i < netList.size(); i++)
    {
        netList[i]->sim_numCellsInGroup[0] = netList[i]->numCellsInGroup[0];
        netList[i]->sim_numCellsInGroup[1] = netList[i]->numCellsInGroup[1];
    }
}


bool Circuit::passBalanceCriterion(int size0)
{
    double lowerbound = (1 - balanceFactor) / 2;
    double upperbound = (1 + balanceFactor) / 2;
    double ratio = (double)size0 / (double)cellList.size();
    return (lowerbound <= ratio && ratio <= upperbound);
}


void Circuit::moveCellAndUpdateNeighboringGain(int i, vector<bool>& locked, vector<int>& cellGain)
{
    int from = cellList[i]->sim_group;
    int to = 1 - from;
    cellList[i]->sim_group = to;
    for (auto& netIdx : cellList[i]->associatedNets)
    {
        // Check for critical nets before the move
        if (netList[netIdx]->sim_numCellsInGroup[to] == 0)
        {
            // increment gains of all free cells on n
            for (auto& cellIdx : netList[netIdx]->associatedCells)
            {
                if (!locked[cellIdx])
                    updateCellGain(cellIdx, 1, cellGain);
            }
        }
        else if (netList[netIdx]->sim_numCellsInGroup[to] == 1)
        {
            // decrement gain of the only T cell on n, if it is free
            for (auto& cellIdx : netList[netIdx]->associatedCells)
            {
                if (cellList[cellIdx]->sim_group == to && !locked[cellIdx])
                    updateCellGain(cellIdx, -1, cellGain);
            }
        }

        // Change F(n) and T(n) to reflect the move
        netList[netIdx]->sim_numCellsInGroup[from]--;
        netList[netIdx]->sim_numCellsInGroup[to]++;

        // Check for critical nets after the move
        if (netList[netIdx]->sim_numCellsInGroup[from] == 0)
        {
            // decrement gains of all free cells on n
            for (auto& cellIdx : netList[netIdx]->associatedCells)
            {
                if (!locked[cellIdx])
                    updateCellGain(cellIdx, -1, cellGain);
            }
        }
        else if (netList[netIdx]->sim_numCellsInGroup[from] == 1)
        {
            // increment gain of the only F cell on n, if it is free
            for (auto& cellIdx : netList[netIdx]->associatedCells)
            {
                if (cellList[cellIdx]->sim_group == from && !locked[cellIdx])
                    updateCellGain(cellIdx, 1, cellGain);
            }
        }
    }
}


void Circuit::updateCellGain(int i, int delta, vector<int>& cellGain)
{
    // Remove Cell i from the bucket list of old gain
    if (cellList[i]->next != nullptr)
        cellList[i]->next->prev = cellList[i]->prev;
    cellList[i]->prev->next = cellList[i]->next;
    cellList[i]->prev = nullptr;
    cellList[i]->next = nullptr;

    // Update the gain of Cell i
    cellGain[i] += delta;

    // Insert Cell i into the bucket list of new gain
    Cell* pTemp = &(bucketOfGroup[cellList[i]->sim_group][cellGain[i] + pMax]);
    if (pTemp->next != nullptr)
        pTemp->next->prev = cellList[i];
    cellList[i]->next = pTemp->next;
    cellList[i]->prev = pTemp;
    pTemp->next = cellList[i];
}


void Circuit::updateCircuitState(const vector<int>& movedCell, int k)
{
    // Perform the first k steps of cell movement
    for (int j = 0; j < k; j++)
    {
        int cellIdx = movedCell[j];
        int from = cellList[cellIdx]->group;
        int to = 1 - from;
        cellList[cellIdx]->group = to;
        sizeOfGroup0 = (from == 0) ? (sizeOfGroup0 - 1) : (sizeOfGroup0 + 1);
        for (auto& netIdx : cellList[cellIdx]->associatedNets)
        {
            netList[netIdx]->numCellsInGroup[from]--;
            netList[netIdx]->numCellsInGroup[to]++;
        }
    }
}


void Circuit::writeOutput(string ofName)
{
    ofstream fout(ofName);
    fout << "Cutsize = " << cutSize() << endl;

    vector<int> group0, group1;
    for (size_t i = 0; i < cellList.size(); i++)
    {
        if (cellList[i]->group == 0)
            group0.push_back(i);
        else
            group1.push_back(i);
    }

    fout << "G1 " << group0.size() << endl;
    for (auto& cellIdx : group0)
        fout << cellList[cellIdx]->name << " ";
    fout << ";" << endl;

    fout << "G2 " << group1.size() << endl;
    for (auto& cellIdx : group1)
        fout << cellList[cellIdx]->name << " ";
    fout << ";" << endl;

    fout.close();
}


int Circuit::cutSize()
{
    int size = 0;
    for (size_t i = 0; i < netList.size(); i++)
    {
        if (netList[i]->numCellsInGroup[0] != 0 && netList[i]->numCellsInGroup[1] != 0)
            size++;
    }

    return size;
}


void Circuit::freePointers()
{
    for (size_t i = 0; i < cellList.size(); i++)
        delete cellList[i];
    for (size_t i = 0; i < netList.size(); i++)
        delete netList[i];
}