//
//  main.cpp
//  Feynman_MCSimulator
//
//  Created by Luis Paulo Santos on 30/07/2023.
//

#include <iostream>
#include <cstdlib>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

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
#include "myReal.h"
#include "CState.h"

#include "PreProcessorSettings.h"

static bool check_command_line(int, const char * []);
static void print_usage (void);


static char fileName[1024], csv_amplitude_fileName[1024];
static char run_label[32];
static CState init_state, final_state;
static int n_threads=1;
static int samples_exp2=20;
static SampleCounter n_samples=1<<samples_exp2;
static bool CSV_amplitude_verification = false;
static bool true_amplitude_given = false;
static myReal true_a_R, true_a_I;

enum TAlgorithms {
    ALL_PATHS=1,
    IS_FORWARD=2,
    BD_SAMPLE=3,
    BD_MIS=4
} ;

static void output_greetings () {
	fprintf (stdout, "\n** Monte Carlo Feynman simulator LPS@2024\n\n");
#ifdef __FLOAT_MP__
	fprintf (stdout, "** \tfloats AS GNU MP class_mpf\n");
#elif defined(__FLOAT_AS_DOUBLE__)
	fprintf (stdout, "** \tfloats AS double\n");
#else
	fprintf (stdout, "** \tfloats AS float\n");
#endif
#ifdef __CSTATE_MP__
	fprintf (stdout, "** \tstates AS GNU MP class_mpz\n");
#elif defined(__CSTATE_VCHAR__)
	fprintf (stdout, "** \tstates AS CState_VChar\n");
#else
	fprintf (stdout, "** \tstates AS unsigned long long\n");
#endif
	fprintf (stdout, "**\n******** \n\n");
}

#ifdef CONVERGENCE_STATS
static void save_stats (std::vector<T_Stats>, bool, myReal, myReal, TAlgorithms );
#endif

static TAlgorithms algorithm;

int main(int argc, const char * argv[]) {
    TCircuit *circuit;
#ifdef CONVERGENCE_STATS
    std::vector<T_Stats> stats;
#endif
    
    output_greetings ();

    if (!check_command_line (argc, argv)) return 1;
    
    printf("Filename = %s\n", fileName);
    
    circuit = read_circuit(fileName);
    if (circuit==NULL) {
        fprintf (stderr, "read_circuit() failed!\n");
        return 1;
    }
    //else fprintf (stderr, "read_circuit() OK!\n");
    
#ifdef DEBUG
    //    print_circuit(circuit);
#endif
    print_circuit_stats (circuit);
    
    fprintf (stdout, "#different paths from PSI_0 to PSI_f = 2 ^ (%d * %d) = 2 ^ %d\n\n", circuit->size->num_qubits, circuit->size->num_layers-1, circuit->size->num_qubits*(circuit->size->num_layers-1));
        
    fflush (stdout);

    bool ret;
    myReal estimateR, estimateI;
    bool CSV_found = false;
    
    // https://www.geeksforgeeks.org/measure-execution-time-function-cpp/
    auto start = high_resolution_clock::now();
    
    // simulate
    switch (algorithm) {
        case ALL_PATHS:
#if defined(__CSTATE_VCHAR__)
            ret = all_paths (circuit, init_state.get_ull(), final_state.get_ull(), estimateR, estimateI, n_threads);
#elif defined(__CSTATE_MP__)
            ret = all_paths (circuit, (unsigned long long)init_state.get_ui(), (unsigned long long)final_state.get_ui(), estimateR, estimateI, n_threads);
#else
            ret = all_paths (circuit, init_state, final_state, estimateR, estimateI, n_threads);
#endif
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
            ret = BD_paths (circuit, init_state, final_state, n_samples, estimateR, estimateI, stats, n_threads, false);
#else
            ret = BD_paths (circuit, init_state, final_state, n_samples, estimateR, estimateI, n_threads, false);
#endif
            break;
        case BD_MIS:
#ifdef CONVERGENCE_STATS
            ret = BD_paths (circuit, init_state, final_state, n_samples, estimateR, estimateI, stats, n_threads, true);
#else
            ret = BD_paths (circuit, init_state, final_state, n_samples, estimateR, estimateI, n_threads, true);
#endif
            break;
            
        default:
            break;
    }
    
    // https://www.geeksforgeeks.org/measure-execution-time-function-cpp/
    auto stop = high_resolution_clock::now();
    
#if defined(__CSTATE_MP__)
    gmp_fprintf (stdout,"< %Zd | U | %Zd > = ", final_state.get_mpz_t(), init_state.get_mpz_t());
#elif defined(__CSTATE_VCHAR__)
    fprintf (stdout,"< %llu | U | %llu > = ", final_state.get_ull(), init_state.get_ull());
#else
    fprintf (stdout,"< %llu | U | %llu > = ", final_state, init_state);
#endif
    
#if defined(__FLOAT_MP__)
    gmp_fprintf (stdout,"%Ff + i %Ff \t", estimateR.get_mpf_t(), estimateI.get_mpf_t() );
#else
    fprintf (stdout,"%f + i %f \t", estimateR, estimateI);
#endif
    
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
    
    if (CSV_amplitude_verification || true_amplitude_given)  {  // compare estimate with true velue from CSV
        
        if (CSV_amplitude_verification) {
            using namespace csv;
            CSVReader reader(csv_amplitude_fileName);
            unsigned long CSV_X;
            CState CCSV_X;
#if defined(__CSTATE_MP__) || defined(__CSTATE_VCHAR__)
            unsigned long const my_final_state = final_state.get_ui();
#else
            unsigned long long const my_final_state = final_state;
#endif
            
            for (auto& row: reader) {
                // Note: Can also use index of column with [] operator
                CSV_X = (unsigned long)(row["psi0"].get<unsigned long long>());
                CCSV_X = (CState) CSV_X;
                if (CCSV_X==init_state) {
                    true_a_R = row[to_string(my_final_state)+"r"].get<float>();
                    true_a_I = row[to_string(my_final_state)+"i"].get<float>();
                    CSV_found = true;
                    break;
                }
            }
#if defined(__CSTATE_MP__)
            gmp_fprintf (stdout,"< %Zd | U | %Zd > = ", final_state.get_mpz_t(), init_state.get_mpz_t());
#elif defined(__CSTATE_VCHAR__)
            fprintf (stdout,"< %llu | U | %llu > = ", final_state.get_ull(), init_state.get_ull());
#else
            fprintf (stdout,"< %llu | U | %llu > = ", final_state, init_state);
#endif
            if (CSV_found) {
                const myReal varR = (estimateR-true_a_R)*(estimateR-true_a_R);
                const myReal varI = (estimateI-true_a_I)*(estimateI-true_a_I);
#ifdef __FLOAT_MP__
                {
                    myReal aux;
                    aux = complex_abs(estimateR-true_a_R, estimateI-true_a_I);
                    gmp_fprintf (stdout,"%Ff + i %Ff TRUE AMPLITUDE L2 error=%Ff\n", true_a_R.get_mpf_t(), true_a_I.get_mpf_t(), aux.get_mpf_t());
                    aux = varR+varI;
                    gmp_fprintf (stdout,"Estimate variance = %Fe\n\n", aux.get_mpf_t());
                }
#else
                fprintf (stdout,"%f + i %f TRUE AMPLITUDE L2 error=%f\n", true_a_R, true_a_I, complex_abs(estimateR-true_a_R, estimateI-true_a_I));
                fprintf (stdout,"Estimate variance = %e\n\n", varR+varI);
#endif
            } else {
                fprintf (stdout,"TRUE AMPLITUDE not found in CSV\n\n");
            }
        }
        else {  // true amplitude given in the command line
//#if defined(__CSTATE_MP__) || defined(__CSTATE_VCHAR__)
#if defined(__CSTATE_MP__)
            gmp_fprintf (stdout,"< %Zd | U | %Zd > = ", final_state.get_mpz_t(), init_state.get_mpz_t());
#elif defined (__CSTATE_VCHAR__)
            fprintf (stdout,"< %llu | U | %llu > = ", final_state.get_ull(), init_state.get_ull());
#else
            fprintf (stdout,"< %llu | U | %llu > = ", final_state, init_state);
#endif
            const myReal varR = (estimateR-true_a_R)*(estimateR-true_a_R);
            const myReal varI = (estimateI-true_a_I)*(estimateI-true_a_I);
#ifdef __FLOAT_MP__
            {
                myReal aux;
                aux = complex_abs(estimateR-true_a_R, estimateI-true_a_I);
                gmp_fprintf (stdout,"%Ff + i %Ff TRUE AMPLITUDE L2 error=%Ff\n", true_a_R.get_mpf_t(), true_a_I.get_mpf_t(), aux.get_mpf_t());
                aux = varR+varI;
                gmp_fprintf (stdout,"Estimate variance = %Fe\n\n", aux.get_mpf_t());
            }
#else
            fprintf (stdout,"%f + i %f TRUE AMPLITUDE L2 error=%f\n", true_a_R, true_a_I, complex_abs(estimateR-true_a_R, estimateI-true_a_I));
            fprintf (stdout,"Estimate variance = %e\n\n", varR+varI);
#endif
        }
    }  // amplitude verification given true amplitude (CSV or command line)
    
    
#ifdef CONVERGENCE_STATS
    // print / save stats
    if ((algorithm==IS_FORWARD || algorithm==BD_SAMPLE || algorithm==BD_MIS) && n_threads > 1) {
        save_stats (stats, CSV_found || true_amplitude_given, true_a_R, true_a_I, algorithm);
    }
#endif
    
    
    return 0;
}


static void print_usage (void) {
    
    fprintf (stderr, "\n\n >>> USAGE <<<\n\n");
    fprintf (stderr, "program -f <circuit_file> -a <algorithm> -i <init_state> -o <final_output_state> -t <n_threads> -s <log2_samples> -r <amplitude_real> <amplitude_img> -l <label>\n\n");
    fprintf (stderr, "\t <circuit_file> - name of .data file describing the circuit to simulate (without the extension)\n");
    fprintf (stderr, "\t <algorithm> :\n\t\t1 - ALL_PATHS\n\t\t2 - FORWARD IMPORTANCE SAMPLING\n");
    fprintf (stderr, "\t\t3 - BIDIRECTIONAL IMPORTANCE SAMPLING\n\t\t4 - BIDIRECTIONAL IMPORTANCE SAMPLING WITH MIS\n");
    fprintf (stderr, "\t <init_state> : integer; \n");
    fprintf (stderr, "\t <final_state> : integer;\n");
    fprintf (stderr, "\t <n_threads> :\n\t\tNumber of threads (default = 1)\n");
    fprintf (stderr, "\t <log2_samples> :\n\t\tif algorithm != ALL_PATHS then log2 number of samples (default = 20 - 2^20 samples)\n");
    fprintf (stderr, "\t <amplitude_real> <amplitude_img> :\n\t\ttrue value (REAL + IMAG) for the transition amplitude\n");
    fprintf (stderr, "\t <label> :\n\t\tOPTIONAL: textual label to be added to the files generated by this run!\n");
    fprintf (stderr,"\n\n");
}


static bool check_command_line(int argc, const char * argv[]) {
    int i;
    
    if (argc<3) {
        fprintf (stderr, "Error: requires command line parameters!");
        print_usage();
        return false;
    }
    
    run_label[0] = '\0'; // make sure we can detect whether one exists

    for (i=1 ; i< argc ; i+=2) {
        // Parse command line arguments
        int opt;
        
        if ((strlen(argv[i]) != 2) || argv[i][0] != '-') {
            fprintf (stderr, "Error: Unknown option \"%s\"!", argv[i]);
            print_usage();
            return false;
        }
        opt = argv[i][1];

        switch (opt) {
            case 'f':  {  // -f : requires filename
                strcpy(fileName, argv[i+1]);
                strcat (fileName, ".data");
                //fprintf (stdout, "Circuit file: %s\n", fileName);
                FILE *f = fopen(fileName, "rb");
                if (f==NULL) {
                    fprintf (stderr, "Error; file %s does not exist!", fileName);
                    return false;
                }
                fclose (f);
                
                // IS there a csv file with the pre computed amplitudes ?
                strcpy(csv_amplitude_fileName, argv[i+1]);
                strcat (csv_amplitude_fileName, ".csv");
                f = fopen(csv_amplitude_fileName, "rb");
                if (f==NULL) {
                    csv_amplitude_fileName[0] = '\0';   // no csv file
                    CSV_amplitude_verification = false;
                }
                else {
                    fprintf (stdout, "CSV file exists: will check accuracy\n");
                    fclose (f);
                    CSV_amplitude_verification = true;
                }

                break;
            }
            case 'a': {  // -a algorithm
                switch (atoi(argv[i+1])) {
                    case 1:    // ALL_PATHS
                        algorithm = ALL_PATHS;
                        fprintf (stdout, "Algorithm :\tALL_PATHS\n");
                        break;
                    case 2:    // Importance Sampling FORWARD
                        algorithm = IS_FORWARD;
                        fprintf (stdout, "Algorithm :\tFORWARD IMPORTANCE SAMPLING\n");
                        break;
                    case 3:    // BiDirectional Sampling
                        algorithm = BD_SAMPLE;
                        fprintf (stdout, "Algorithm :\tBIDIRECTIONAL IMPORTANCE SAMPLING\n");
                        break;
                    case 4:    // BiDirectional Sampling MIS
                        algorithm = BD_MIS;
                        fprintf (stdout, "Algorithm :\tBIDIRECTIONAL IMPORTANCE SAMPLING WITH MIS\n");
                        break;
                    default:
                        fprintf (stderr, "Error; Unknown algorithm\n!");
                        print_usage();
                        return false;
                        break;
                }
                break;
            }
            case 't':
            {
                n_threads = atoi(argv[i+1]);
                fprintf (stdout, "Number of threads = %d\n", n_threads);
                break;
            }
            case 'i':
            {
#ifdef __CSTATE_MP__
                init_state = argv[i+1];
#else
                init_state = strtoull(argv[i+1], NULL, 10);
#endif
                break;
            }
            case 'o':
            {
#ifdef __CSTATE_MP__
                final_state = argv[i+1];
#else
                final_state = strtoull(argv[i+1], NULL, 10);
#endif
                break;
            }
            case 's':
            {
                samples_exp2 = atoi(argv[i+1]);
                n_samples = ((SampleCounter)(1 << samples_exp2));
                fprintf (stdout, "Number of samples = 2^%d = ", samples_exp2);
                fprintf (stdout, "%llu\n", n_samples);
                break;
            }
            case 'r':
            {
#ifdef __FLOAT_MP__
                true_a_R = argv[i+1];
                true_a_I = argv[i+2];
                gmp_fprintf (stdout, "True amplitude = %Fe + j %Fe\n", true_a_R.get_mpf_t(), true_a_I.get_mpf_t());
#else
                true_a_R = atof(argv[i+1]);
                true_a_I = atof(argv[i+2]);
                fprintf (stdout, "True amplitude = %e + j %e\n", true_a_R, true_a_I);
#endif
                true_amplitude_given = true;
                i = i+1;

                break;
            }
            case 'l':
            {
                if (strlen(argv[i+1])>30) {
                    fprintf (stderr, "Error: Run label cannot be larger than 30 characters: %s!\n", argv[i+1]);
                    print_usage();
                    return false;
                }
                strcpy(run_label, argv[i+1]);
                fprintf (stdout, "Run label: %s\n", run_label);
                break;
            }
            default:
            {
                fprintf (stderr, "Error: Unknown option: %s!\n", argv[i]);
                print_usage();
                return false;
                break;
            }
        }
    }
    return true;
}

#ifdef CONVERGENCE_STATS
static void save_stats (std::vector<T_Stats> stats, bool true_exists, myReal trueR, myReal trueI, TAlgorithms algorithm) {
    char csv_stats_fileName[1024], alg_str[16];
    FILE *f;
    myReal stat_estimateR, stat_estimateI, varR, varI;
    
    switch (algorithm) {
        case IS_FORWARD:
            snprintf(alg_str, 16, "IS");
            break;
        case BD_SAMPLE:
            snprintf(alg_str, 16, "BD");
            break;
        case BD_MIS:
            snprintf(alg_str, 16, "BD_MIS");
            break;
        default:
            snprintf(alg_str, 16, "UKNOWN");
            break;
    }
#if defined(__CSTATE_MP__) 
    if (strlen(run_label)>0) {  // run label exists
        gmp_snprintf (csv_stats_fileName, 1024, "%s_stats_%s_%Zd_%Zd_%d_%d_%s.csv", fileName, alg_str, init_state.get_mpz_t(), final_state.get_mpz_t(), samples_exp2, n_threads, run_label);
    } else {
        gmp_snprintf (csv_stats_fileName, 1024, "%s_stats_%s_%Zd_%Zd_%d_%d.csv", fileName, alg_str, init_state.get_mpz_t(), final_state.get_mpz_t(), samples_exp2, n_threads);
    }
#elif defined(__CSTATE_VCHAR__)
    if (strlen(run_label)>0) {  // run label exists
        snprintf (csv_stats_fileName, 1024, "%s_stats_%s_%llu_%llu_%d_%d_%s.csv", fileName, alg_str, init_state.get_ull(), final_state.get_ull(), samples_exp2, n_threads, run_label);
    } else {
        snprintf (csv_stats_fileName, 1024, "%s_stats_%s_%llu_%llu_%d_%d.csv", fileName, alg_str, init_state.get_ull(), final_state.get_ull(), samples_exp2, n_threads);
    }
#else
    if (strlen(run_label)>0) {  // run label exists
        snprintf (csv_stats_fileName, 1024, "%s_stats_%s_%llu_%llu_%d_%d_%s.csv", fileName, alg_str, init_state, final_state, samples_exp2, n_threads, run_label);
    } else {
        snprintf (csv_stats_fileName, 1024, "%s_stats_%s_%llu_%llu_%d_%d.csv", fileName, alg_str, init_state, final_state, samples_exp2, n_threads);
    }
#endif
    
    f = fopen(csv_stats_fileName,"wt");

    if (true_exists) {
#ifdef __FLOAT_MP__
        gmp_fprintf (f, "trueR , trueI\n%Ff , %Ff\n", trueR.get_mpf_t(), trueI.get_mpf_t());
#else
        fprintf (f, "trueR , trueI\n%f , %f\n", trueR, trueI);
#endif
    }
    fprintf (f, "n_samples , time, n_paths , ");
#ifdef NON_ZERO_PATHS
    fprintf (f, "n_NZero_paths , ");
#endif
    fprintf (f, "estimateR , estimateI");
    if (true_exists) {
        fprintf (f, " , varianceR , varianceI");
    }
    fprintf (f, "\n");

    for (auto & stat : stats) {
        double const dn_Paths = (double)stat.n_Paths;
        stat_estimateR = stat.sumR / dn_Paths;
        stat_estimateI = stat.sumI / dn_Paths;
        // keep in mind that variance is the sum of the real and complex variances
        // we are storing both terms separated here to allow for posterior processing
#ifdef __FLOAT_MP__
        gmp_fprintf (f, "%llu , %lld , %llu , ", stat.n_samples, stat.time_us, stat.n_Paths);
#ifdef NON_ZERO_PATHS
        gmp_fprintf (f, "%llu , ", stat.n_nonZero_paths );
#endif
        gmp_fprintf (f, "%Ff , %Ff", stat_estimateR.get_mpf_t(), stat_estimateI.get_mpf_t());
#else
        fprintf (f, "%llu , %lld , %llu , ", stat.n_samples, stat.time_us, stat.n_Paths);
#ifdef NON_ZERO_PATHS
        fprintf (f, "%llu , ", stat.n_nonZero_paths );
#endif
        fprintf (f, "%f , %f", stat_estimateR, stat_estimateI);
#endif

        if (true_exists) {
            varR = (stat_estimateR-trueR)*(stat_estimateR-trueR);
            varI = (stat_estimateI-trueI)*(stat_estimateI-trueI);
#ifdef __FLOAT_MP__
            gmp_fprintf (f, " , %Ff , %Ff", varR.get_mpf_t(), varI.get_mpf_t());
#else
            fprintf (f, " , %f , %f", varR, varI);
#endif
        }
        fprintf (f, "\n");
    }
    fclose (f);
}
#endif
