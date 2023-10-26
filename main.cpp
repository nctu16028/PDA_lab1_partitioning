#include "circuit.h"
#include <iostream>
using namespace std;

int main(int argc, char* argv[])
{
    if (argc != 3)
    {
        cout << "Usage: ./Lab1 [input file name] [output file name]" << endl;
        return -1;
    }

    Circuit circuit;
    string inputFileName(argv[1]);
    circuit.readInput(inputFileName);

    circuit.partition();

    string outputFileName(argv[2]);
    circuit.writeOutput(outputFileName);

    return 0;
}