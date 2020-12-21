/*
Now, here, you see, it takes all the running you can do, to keep in the same place. -Lewis Carrol. 

THE RED QUEEN ALGORITHM. 
AN APPROACH TO DYNAMIC TIME WARPING USING GENETIC PROGRAMMING
KEVIN R. LAWRENCE
FEBRUARY 2020

*/
#include <vector>
#include <iostream>
#include <ctime>
#include <stdlib.h>
#include <algorithm>
#include <numeric>
#include <math.h>
#include <iomanip>
#include <fstream>

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



bool parthenogenesis(vector<int> &parent, int N, int M){ //"Virgin Birth"
    int tracker_n = 0; //Trackers to keep track of position of matrix (N,M)
    int tracker_m = 0;

    if(N <= 0 || M <= 0){
        return false;  //exception handling ;)
    }
    for(;;){
        int g1 = rand() % 2;
        int g2 = rand() % 2; //Random 1 or 0.
        if(tracker_m >= M && tracker_n >= N){break;} //break (N,M) is reached

        tracker_m += g1; //Add more genes.
        tracker_n += g2;
        parent.push_back(g1);
        parent.push_back(g2);
    }
    return true;
}

double euclidian_distance(vector<double> a, vector<double> b){
    vector<double> euc_sum;
    for(size_t i = 0; i < a.size(); i++){
        euc_sum.push_back(pow((b[i] - a[i]), 2)); //Distance formula, given by pythagoras
    }
    return sqrt(accumulate(euc_sum.begin(), euc_sum.end(), 0.0));
}

bool Sa(vector<vector<double> > &m, int degree){
    
}

void proofread(Generation gen, vector<int> &c){
    int m_tracker = 0; //Trackers to keep track of position in matrix.
    int n_tracker = 0;
    if(c.size() % 2 != 0){c.push_back(rand() % 2);} //Adjust for uneven crossovers.

    for(size_t i = 0; i < c.size(); i = i+2){ //Iterate through child, 2 by 2.
        m_tracker += c[i];//track projected path.
        n_tracker += c[i+1];
    }
    /*
     * NOTE: it doesn't matter at this point if we are larger than N or M, the fitness
     * function will take care of that. The following while loop will terminate when the
     * tracker_n and tracker_m contain enough chromosomes to reach c(M,N).
     */
    while(m_tracker < (gen.M - 1) || n_tracker < (gen.N - 1)){
        int g1 = rand() % 2;
        int g2 = rand() % 2;
        c.push_back(g1);
        c.push_back(g2);
        n_tracker += g2;
        m_tracker += g1;
    }
    return;
}

double fitness(vector<vector<double> > euclidian_matrix, vector<int> c){
    int M = euclidian_matrix.size()-1;
    int N = euclidian_matrix[0].size()-1;
    int m_tracker = 0;
    int n_tracker = 0;
    vector<double> f;
    f.push_back(euclidian_matrix[0][0]); //Both signals start at [0],[0] so add back that position

    for(size_t i = 0; i<c.size()-2; i=i+2){ //iterate through chromosome, 2 by 2.
        if(c[i] == 0 && c[i+1] == 0){continue;} //if both are 0 c(0,0) continue.

        bool flag = true; //flag to shut down once c(M,N) is reached.
        if(n_tracker < N){ //Prevent from exceeding N range
            n_tracker += c[i+1];
            flag = false;
        }
        if(m_tracker < M){ //Prevent from exceeding M range
            m_tracker += c[i];
            flag = false;
        }
        if(flag == true){
            break; //Break if flag is thrown.
        }
        f.push_back(euclidian_matrix[m_tracker][n_tracker]);//Add euclidian distance to vector
    }
    return(accumulate(f.begin(), f.end(), 0.0)/f.size()); //Return the mean distance.
}

void perform_mutation(vector<int> &c){
    for(int i = 0; i < c.size()-2; i = i + 2){
        int r = rand() % 300 + 1;//random integer between 1 and 300.
        if(r == 150){ //0.03 means that among 300 chromosomes, 1 will be mutated.
            int mut_opt = rand() % 3; //random int between 0 and 2.
            if(mut_opt == 0){ //insertion
                c.push_back(rand() % 2);
                c.push_back(rand() % 2);
            }
            else if(mut_opt == 1){ //deletion
                c.erase(c.begin()+i, c.begin()+(i+2));
            }
            else{ //substitiution
                c[i] = c[i] == 0 ? 1 : 0;
                c[i+1] = c[i+1] == 0 ? 1 : 0;
            }
        }
    }
}

void crossover(Generation gen, vector<int> &child1, vector<int> &child2){
    double pivot1 = rand() % 50 + 20; //Random double between 10 - 90
    double pivot2;
    do{
        pivot2 = rand() % 60 + 20; //Random double between 10 - 90
    }while(pivot1 == pivot2 || fabs(pivot1 - pivot2) < 4); //Exception handling for identical pivot points.
    pivot1 = pivot1/100.0; // transforming into percentage decimals.
    pivot2 = pivot2/100.0;
     
    int temp_m = pivot1 < pivot2 ? pivot2 : pivot1; //Arranging so that pivot 1 < Pivot 2
    if(pivot1 > pivot2){pivot1 = pivot2; pivot2 = temp_m;}
    //Adding chromosomes from (0 - pivot 1) from parent 1 to child1.
    for(size_t i = 0; i < floor(gen.p1.size()*pivot1); i++){
        child1.push_back(gen.p1[i]);
    }
    //Adding chromosomes from (0 - pivot 1) from parent 2 to child2.
    for(size_t i = 0; i < floor(gen.p2.size()*pivot1); i++){
        child2.push_back(gen.p2[i]);
    }
    //Adding chromosomes from (pivot 1 - pivot 2) from parent 2 to child 1.
    for(size_t i = floor(gen.p1.size()*pivot1); i < floor(gen.p1.size()*pivot2); i++){
        child2.push_back(gen.p1[i]);
    }
    //Adding chromosomes from (pivot 1 - pivot 2) from parent 1 to child 2
    for(size_t i = floor(gen.p2.size()*pivot1); i < floor(gen.p2.size()*pivot2); i++){
        child1.push_back(gen.p2[i]);
    }
    //Adding chromosomes from (pivot 2 - M) from parent 1 to child 1
    for(size_t i = floor(gen.p1.size()*pivot2); i < gen.p1.size(); i++){
        child1.push_back(gen.p1[i]);
    }
    //Adding chromosomes from (pivot 2 - N) from parent 2 to child 2.
    for(size_t i = floor(gen.p2.size()*pivot2); i < gen.p2.size(); i++){
        child2.push_back(gen.p2[i]);
    }
    perform_mutation(child1); //Mutate the children
    perform_mutation(child2);
    proofread(gen, child1); //Check for errors.
    proofread(gen, child2);
   
    return; 
}

Generation capacocha(Generation g, vector<int> c, double f, int &gc){// "Sacrifice of Children"
    Generation next_gen; //Create return value. 
    next_gen.M = g.M;
    next_gen.N = g.N;

    if(f < g.fitness_p2){ //if new fitness out-ranks parent 2
        if(f<g.fitness_p1){ //if new fitness out-ranks parent 1 as well. 
            next_gen.fitness_p1 = f;
            next_gen.p1 = c;
            next_gen.fitness_p2 = g.fitness_p1;
            next_gen.p2 = g.p1;
            gc++;
        }
        else{
            next_gen.fitness_p1 = g.fitness_p1;
            next_gen.p1 = g.p1;
            next_gen.fitness_p2 = f;
            next_gen.p2 = c;
            gc++;
        }
        return next_gen;
    }
    else{
        return g;
    }

}

vector<vector<double> > Form_Matrix(vector<vector<double> > A, vector<vector<double> > B){
    int M = A.size(); //Pull M and N from input vectors.
    int N = B.size();
    vector<vector<double> > euclidian_matrix; //Matrix for return.

    for(int i = 0; i < M; i++){
        vector<double> matrix_row; //Temporary vector

        for(int j = 0; j < N; j++){
            matrix_row.push_back(euclidian_distance(A[i], B[j])); //Add euclidian distances
        }
        euclidian_matrix.push_back(matrix_row); //Add temporary vector to return matrix.
    }
    return euclidian_matrix;
}

void Red_Queen(vector<vector<double> > A, vector<vector<double> > B, int MAX_ITERATIONS){
    int GEN_COUNT = 0;

    Generation fittest; //Generation to hold the fittest through iterative generations.
    fittest.M = A.size()-1; //Set size to adjusted vector sizes.
    fittest.N = B.size()-1;
    vector<int> Parent1; //Create vectors for the chromosomes of parents.
    vector<int> Parent2;
    double current_fit = 0; 
    
 
    parthenogenesis(Parent1, fittest.N, fittest.M); //Virgin Birth of Parent 1
    parthenogenesis(Parent2, fittest.N, fittest.M); //Virgin Birth of Parent 2
    vector<vector<double> > euclidian_matrix = Form_Matrix(A, B); //Euclidian_matrix formation
    double p1_fitness = fitness(euclidian_matrix, Parent1); //Calculate fitness of Parent 1
    double p2_fitness = fitness(euclidian_matrix, Parent2); //Calculate fitness of Parent 2
    /*
        The following lines arange the components of the fittest generation so that the fittest
        parent is parent 1 and the lesser fit parent is parent 2. Keeping them in this order makes
        it easier later on. 
    */
    double max_temp = p1_fitness > p2_fitness ? p1_fitness : p2_fitness;
    if(max_temp == p2_fitness){
        p1_fitness = p2_fitness;
        p2_fitness = max_temp;
        vector<int> temp_v = Parent2;
        Parent2 = Parent1;
        Parent1 = temp_v;
    }
    fittest.p1 = Parent1;
    fittest.p2 = Parent2;
    fittest.fitness_p1 = p1_fitness;
    fittest.fitness_p2 = p2_fitness;
    //Begin to iterate through generations. 
    for(int i = 0; i < MAX_ITERATIONS; i++){
        vector<int> Child1; //each generation has 2 children. 
        vector<int> Child2; 
        crossover(fittest, Child1, Child2); //crossover occurs 
        double c1_fitness = fitness(euclidian_matrix, Child1); //calculate fitness of children
        double c2_fitness = fitness(euclidian_matrix, Child2);
        fittest = capacocha(fittest, Child1, c1_fitness, GEN_COUNT); //slaughter the weakest
        fittest = capacocha(fittest, Child2, c2_fitness, GEN_COUNT);
        if(current_fit != fittest.fitness_p1){ //if new fittest, print value. 
            cout<<"CURRENT GENERATION FITNESS: "<<fittest.fitness_p1<<endl; 
            current_fit = fittest.fitness_p1; 
        }
    }

    cout<<GEN_COUNT<<" GENERATIONS NEEDED"<<endl; //print the total number of parents killed. 


}


int main(){
    srand(time(NULL)); //seed rand
    ifstream in_f("A_0_MFCC.txt"); //sample MFCC files
    ifstream in_fi("A_00_MFCC.txt");
    double temp;
    vector<double> A_row;
    vector<vector<double> > A; 
    int max_gen; 
    cout<<"Enter value max number of generations"<<endl; 
    cin>>max_gen; 
    //Fill in signal A. 
    while(in_f >> temp){
        A_row.push_back(temp);
        if(A_row.size() == 26){
            A.push_back(A_row);
            A_row.clear();
        }
    }

    vector<double> B_row;
    vector<vector<double > > B;
    //Fill in signal B
    while(in_fi >> temp){
        B_row.push_back(temp);
        if(B_row.size() == 26){
            B.push_back(B_row);
            B_row.clear();
        }
    }
    //Begin Algorithm. 
    Red_Queen(A, B, max_gen);
}
