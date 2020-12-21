#ifndef RED_QUEEN_H
#define RED_QUEEN_H
#include <vector>

using namespace std; 

struct Generation
{
    int N; //Number of signal A frames
    int M; //Number of signal B frames
    vector<int> p1; //Parent 1's chromosomes
    vector<int> p2; //Parent 2's chromosomes
    double fitness_p1; //respective fitness scores.
    double fitness_p2;
};

void Red_Queen(vector<vector<double> > A, vector<vector<double> > B, int MAX_ITERATIONS);

#endif // RED_QUEEN
