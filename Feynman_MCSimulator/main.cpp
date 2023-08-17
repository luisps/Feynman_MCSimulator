//
//  main.cpp
//  Feynman_MCSimulator
//
//  Created by Luis Paulo Santos on 30/07/2023.
//

#include <iostream>
#include <cstdlib>
#include <stdlib.h>

using namespace std;

// https://www.geeksforgeeks.org/measure-execution-time-function-cpp/
#include <chrono>
using namespace std::chrono;

#include "circuit.h"
#include "all_paths.hpp"

static bool check_command_line(int, const char * []);



static const char *fileName;
static unsigned long long init_state, final_state;
static bool loop_init_states, loop_final_states;
static int n_threads=1;

enum TAlgorithms {
    ALL_PATHS=1
} ;

static TAlgorithms algorithm;

int main(int argc, const char * argv[]) {
    TCircuit *circuit;
    
    if (!check_command_line (argc, argv)) return 1;
    
    printf("Filename = %s\n", fileName);
    
    circuit = read_circuit(fileName);
    if (circuit==NULL) {
        fprintf (stderr, "read_circuit() failed!\n");
        return 1;
    }
    else fprintf (stdout, "read_circuit() OK!\n");
    
    print_circuit_stats (circuit);
    fprintf(stdout, "\n");
    fprintf(stderr, "\n");
    fflush(stderr);
    fflush(stdout);
    
 
    static unsigned long long init_state_range_start, init_state_range_end;
    static unsigned long long final_state_range_start, final_state_range_end;
    const unsigned long long NBR_STATES = 1 << circuit->size->num_qubits;
    if (loop_init_states)  {  // loop over all
        init_state_range_start = 0ull;
        init_state_range_end = NBR_STATES;
    }
    else {
        init_state_range_start = init_state;
        init_state_range_end = init_state+1;
    }
    if (loop_final_states)  {  // loop over all
        final_state_range_start = 0ull;
        final_state_range_end = NBR_STATES;
    }
    else {
        final_state_range_start = final_state;
        final_state_range_end = final_state+1;
    }
    
    for (init_state = init_state_range_start ; init_state < init_state_range_end ; init_state++) {
        
        for (final_state = final_state_range_start ; final_state < final_state_range_end ; final_state++) {
            
            bool ret;
            float estimateR, estimateI;
            // https://www.geeksforgeeks.org/measure-execution-time-function-cpp/
            auto start = high_resolution_clock::now();
            
            // simulate
            switch (algorithm) {
                case ALL_PATHS:
                    ret = all_paths (circuit, init_state, final_state, estimateR, estimateI, n_threads);
                    break;
                    
                default:
                    break;
            }
            
            // https://www.geeksforgeeks.org/measure-execution-time-function-cpp/
            auto stop = high_resolution_clock::now();
            
            fprintf (stdout,"< %llu | U | %llu> = %f + i %f \t", init_state, final_state, estimateR, estimateI);
            
            
            // https://www.geeksforgeeks.org/measure-execution-time-function-cpp/
            // Subtract stop and start timepoints and
            // cast it to required unit. Predefined units
            // are nanoseconds, microseconds, milliseconds,
            // seconds, minutes, hours. Use duration_cast()
            // function.
            auto duration = duration_cast<microseconds>(stop - start);
            
            // To get the value of duration use the count()
            // member function on the duration object
            if (duration.count() <= 1000)
                fprintf (stdout, "T = %lld us\n", duration.count());
            else if (duration.count() <= 1e6)
                fprintf (stdout, "T = %.1f ms\n", ((float)duration.count())/1000.f);
            else
                fprintf (stdout, "T = %.1f s\n", ((float)duration.count())/1000000.f);

        } // iterate over final_states
    }  // iterate over init_states
    
    return 0;
}


static bool check_command_line(int argc, const char * argv[]) {
    
    if (argc<2) {
        fprintf (stderr, "Error; requires filename!");
        return false;
    }
    fileName = argv[1];
    
    if (argc<3) {
        fprintf (stderr, "Error; requires algorithm!");
        return false;
    }
    switch (atoi(argv[2])) {
        case 1:    // ALL_PATHS
            algorithm = ALL_PATHS;
            break;
        default:
            algorithm = ALL_PATHS;
            break;
    }
    
    if (argc<5) {
        fprintf (stderr, "Error; requires initial and final states!");
        return false;
    }
    if (argv[3][0]=='a')  { // loop over all init states
        loop_init_states=true;
    }
    else {
        loop_init_states=false;
        init_state = strtoull(argv[3], NULL, 10);
    }
    if (argv[4][0]=='a')  { // loop over all final states
        loop_final_states = true;
    }
    else {
        loop_final_states = false;
        final_state = strtoull(argv[4], NULL, 10);
    }

    if (argc>5) {
        if (argv[5][0]=='a')  { // autonomously set the nbr of threads
            n_threads=1;
        }
        else {
            n_threads=atoi(argv[5]);
        }
    }
    return true;
}
