#include "simFid.hpp"
#include <omp.h>
#include <cstdlib>

/* This example file computes the first 10 combinations of loss and double photon creation. 
After compilation you can call it with "./test #n" where #n is the number of core you'd like to use.
The result will be written into the folder "resultfiles". */ 

int main(int argc, char* argv[]){
    int size = atoi(argv[1]);
    if (size==0){
        return 0;
    }
    std::vector<float> angErrs(15, 0.0);
    #pragma omp parallel for
    for (int i=0; i<size; i++){
        schedulerGHZshuffled({0.95, 0.99}, angErrs, "resultfiles/", 0, 10, 0, i, size, "../shuffle.txt");
    }
}