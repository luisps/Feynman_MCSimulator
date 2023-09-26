//
//  main.cpp
//  Feynman_MCSimulator
//
//  Created by Luis Paulo Santos on 30/07/2023.
//

#include <iostream>
#include <cstdlib>
#include <stdlib.h>
//#include <string>

using namespace std;

// https://www.geeksforgeeks.org/measure-execution-time-function-cpp/
#include <chrono>
using namespace std::chrono;

#include "circuit.h"
#include "all_paths.hpp"
#include "IS_paths.hpp"
#include "BD_paths.hpp"
#include "complex.h"
#include "csv.hpp"   // from https://github.com/vincentlaucsb/csv-parser

#include "PreProcessorSettings.h"

static bool check_command_line(int, const char * []);
static void print_usage (void);
#ifdef CONVERGENCE_STATS
static void save_stats (std::vector<T_Stats>, bool, float, float);
#endif


static char fileName[1024], csv_amplitude_fileName[1024];
static unsigned long long init_state, final_state;
static bool loop_init_states, loop_final_states;
static int n_threads=1;
static unsigned long long n_samples=1ull<<20;
static bool CSV_amplitude_verification = false;

enum TAlgorithms {
    ALL_PATHS=1,
    IS_FORWARD=2,
    BD_SAMPLE=3
} ;

static TAlgorithms algorithm;

int main(int argc, const char * argv[]) {
    TCircuit *circuit;
#ifdef CONVERGENCE_STATS
    std::vector<T_Stats> stats;
#endif
    
    if (!check_command_line (argc, argv)) return 1;
    
    printf("Filename = %s\n", fileName);
    
    circuit = read_circuit(fileName);
    if (circuit==NULL) {
        fprintf (stderr, "read_circuit() failed!\n");
        return 1;
    }
    else fprintf (stdout, "read_circuit() OK!\n");
    
    #ifdef DEBUG
        print_circuit(circuit);
    #endif
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
            float CSV_Ar=0.f, CSV_Ai=0.f;
            bool CSV_found = false;

            // https://www.geeksforgeeks.org/measure-execution-time-function-cpp/
            auto start = high_resolution_clock::now();
            
            // simulate
            switch (algorithm) {
                case ALL_PATHS:
                    ret = all_paths (circuit, init_state, final_state, estimateR, estimateI, n_threads);
                    break;
                    
                case IS_FORWARD:
#ifdef CONVERGENCE_STATS
                    ret = IS_paths (circuit, init_state, final_state, n_samples, estimateR, estimateI, stats, n_threads);
#else
                    ret = IS_paths (circuit, init_state, final_state, n_samples, estimateR, estimateI, n_threads);
#endif
                    break;
                case BD_SAMPLE:
#ifdef CONVERGENCE_STATS
                    ret = BD_paths (circuit, init_state, final_state, n_samples, estimateR, estimateI, stats, n_threads);
#else
                    ret = BD_paths (circuit, init_state, final_state, n_samples, estimateR, estimateI, n_threads);
#endif
                    break;

                default:
                    break;
            }
            
            // https://www.geeksforgeeks.org/measure-execution-time-function-cpp/
            auto stop = high_resolution_clock::now();
            
            fprintf (stdout,"< %llu | U | %llu > = %f + i %f \t", final_state, init_state, estimateR, estimateI);
                        
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

            if (CSV_amplitude_verification)  {  // compare estimate with true velue from CSV
                using namespace csv;
                CSVReader reader(csv_amplitude_fileName);
                unsigned long long CSV_X;

                for (auto& row: reader) {
                    // Note: Can also use index of column with [] operator
                    CSV_X = row["psi0"].get<unsigned long long>();
                    if (CSV_X==init_state) {
                        CSV_Ar = row[to_string(final_state)+"r"].get<float>();
                        CSV_Ai = row[to_string(final_state)+"i"].get<float>();
                        CSV_found = true;
                        break;
                    }
                }
                if (CSV_found) {
                    fprintf (stdout,"< %llu | U | %llu > = %f + i %f TRUE AMPLITUDE L2 error=%f\n\n", final_state, init_state, CSV_Ar, CSV_Ai, complex_abs(estimateR-CSV_Ar, estimateI-CSV_Ai));
                } else {
                    fprintf (stdout,"< %llu | U | %llu > TRUE AMPLITUDE not found in CSV\n\n", final_state, init_state);
                }

            }  // CSV_amplitude verification
            
            
#ifdef CONVERGENCE_STATS
            // print / save stats
            if ((algorithm==IS_FORWARD || algorithm==BD_SAMPLE) && n_threads > 1) {
                save_stats (stats, CSV_found, CSV_Ar, CSV_Ai);
            }
#endif

        } // iterate over final_states
    }  // iterate over init_states
    
    return 0;
}


static void print_usage (void) {
    
    fprintf (stderr, "\n\n >>> USAGE <<<\n\n");
    fprintf (stderr, "program <circuit_file> <algorithm> <init_state> <final_state> <n_threads> <arg1>\n\n");
    fprintf (stderr, "\t <circuit_file> - name of .data file describing the circuit to simulate (without the extension)\n");
    fprintf (stderr, "\t <algorithm> :\n\t\t1 - ALL_PATHS\n\t\t2 - FORWARD IMPORTANCE SAMPLING\n");
    fprintf (stderr, "\t <init_state> : integer; if 'a' executed for all possible input states!\n");
    fprintf (stderr, "\t <final_state> : integer; if 'a' executed for all possible final states!\n");
    fprintf (stderr, "\t <arg1> :\n\t\tif algorithm = FORWARD IMPORTANCE SAMPLING then log2 number of samples\n");
    fprintf (stderr,"\n\n");
}


static bool check_command_line(int argc, const char * argv[]) {
        
    if (argc<2) {
        fprintf (stderr, "Error; requires filename (without the .data extension)!");
        print_usage();
        return false;
    }
    {
        strcpy(fileName, argv[1]);
        strcat (fileName, ".data");
        FILE *f = fopen(fileName, "rb");
        if (f==NULL) {
            fprintf (stderr, "Error; file %s does not exist!", fileName);
            return false;
        }
        fclose (f);

        // IS there a csv file with the pre computed amplitudes ?
        strcpy(csv_amplitude_fileName, argv[1]);
        strcat (csv_amplitude_fileName, ".csv");
        f = fopen(csv_amplitude_fileName, "rb");
        if (f==NULL) {
            csv_amplitude_fileName[0] = '\0';   // no csv file
            CSV_amplitude_verification = false;
        }
        else {
            fprintf (stderr, "CSV file exists: will check accuracy\n");
            fclose (f);
            CSV_amplitude_verification = true;
        }
    }
    
    if (argc<3) {
        fprintf (stderr, "Error; requires algorithm!");
        print_usage();
        return false;
    }
    switch (atoi(argv[2])) {
        case 1:    // ALL_PATHS
            algorithm = ALL_PATHS;
            break;
        case 2:    // Importance Sampling FORWARD
            algorithm = IS_FORWARD;
            break;
        case 3:    // BiDirectional Sampling
            algorithm = BD_SAMPLE;
            break;
        default:
            fprintf (stderr, "Error; <algorithm> :\n\t\t1 - ALL_PATHS\n\t\t2 - FORWARD IMPORTANCE SAMPLING\n!");
            print_usage();
            return false;
            break;
    }
    
    if (argc<5) {
        fprintf (stderr, "Error; requires initial and final states!");
        print_usage();
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
    fprintf (stderr, "Number of threads = %d\n", n_threads);

    if (argc>6) { // set the nbr of samples
        int exp2 = atoi(argv[6]);
        n_samples = ((unsigned long long)(1ull << exp2));
        fprintf (stderr, "Number of samples = %llu\n", n_samples);
    }
    else if (algorithm==IS_FORWARD) {
        fprintf (stderr, "Number of samples not specified in command line. Defaulting to %llu\n", n_samples);
    }

    return true;
}


#ifdef CONVERGENCE_STATS
static void save_stats (std::vector<T_Stats> stats, bool true_exists, float trueR, float trueI, TAlgorithms algorithm) {
    char csv_stats_fileName[1024], alg_str[16];
    FILE *f;
    float stat_estimateR, stat_estimateI, varR, varI;
    
    switch (algorithm) {
        case IS_FORWARD:
            snprintf(alg_str, 16, "IS");
            break;
        case BD_SAMPLE:
            snprintf(alg_str, 16, "BD");
            break;
        default:
            snprintf(alg_str, 16, "UKNOWN");
            break;
    }
    snprintf (csv_stats_fileName, 1024, "%s_stats_%s_%llu_%llu.csv", fileName, alg_str, init_state, final_state);

    f = fopen(csv_stats_fileName,"wt");

    if (true_exists) {
        fprintf (f, "trueR , trueI\n%f , %f\n", trueR, trueI);
        fprintf (f, "n_samples , estimateR , estimateI, varianceR, varianceI\n");
    } else
        fprintf (f, "n_samples , estimateR , estimateI\n");

    for (auto & stat : stats) {
        // compute running estimate
        stat_estimateR = stat.sumR / ((float)stat.n_samples);
        stat_estimateI = stat.sumI / ((float)stat.n_samples);
        
        if (true_exists) {
            varR = (stat_estimateR-trueR)*(stat_estimateR-trueR);
            varI = (stat_estimateI-trueI)*(stat_estimateI-trueI);
            fprintf (f, "%llu , %f , %f, %f, %f\n", stat.n_samples, stat_estimateR, stat_estimateI, varR, varI);
        }
        else {
            fprintf (f, "%llu , %f , %f\n", stat.n_samples, stat_estimateR, stat_estimateI);
        }
    }
    fclose (f);
}
#endif
