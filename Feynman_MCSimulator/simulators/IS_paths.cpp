//
//  IS_paths.cpp
//  Feynman_MCSimulator
//
//  Created by Luis Paulo Santos on 17/08/2023.
//

#include "IS_paths.hpp"
#include <thread>
#include <mutex>
#include <condition_variable>
#include <vector>
#include <random>
#include "pcg_random.hpp"

#include "complex.h"
#include "layer.hpp"
#include "PreProcessorSettings.h"

// For time stats
// https://www.geeksforgeeks.org/measure-execution-time-function-cpp/
#include <chrono>
using namespace std::chrono;


/*
 * MASTER SLAVE
 *
 *  We need:
 *
 *       std::mutex mx_m2s;  // for sending task from the master to slaves
 *       std::condition_variable cv_m2s;
 *       std::mutex mx_mW2s;  // for Signaling the master is waiting for a resut from a slave
 *       std::condition_variable cv_mW2s;
 *  Predicate Variables, such  as
 *       FOR MASTER 2 SLAVE
 *       bool terminate = false;  // set to true for thread to terminate instead of processing
 *       bool taskReady = false  // is there actually a task ?
 *       int samples2proc = 0;   // number of samples for this thread to process (this is the task)
 *       FOR SLAVE2MASTER
 *       bool resAvailable = false;  // is there actually a result?
 *       XXXX res=...;               // the result : problem dependent
 */

static std::mutex mx_m2s, mx_s2m;
static std::condition_variable cv_m2s, cv_s2m;

// Shared predicate variables
// Master to Slave
static bool terminate = false;


#ifdef NON_ZERO_PATHS
static bool IS_paths_NOVEC_NOT (TCircuit *c,
                CState init_state, CState final_state, const SampleCounter n_samples,
                myReal &sumR, myReal &sumI, SampleCounter& non_zero_paths);
static bool IS_paths_NOVEC_T (TCircuit *c,
                              CState init_state, CState final_state, SampleCounter& samples2proc, bool& taskReady, bool& resAvailable, myReal& T_sumR, myReal& T_sumI, SampleCounter& processedSamples, SampleCounter& non_zero_paths);
#else
static bool IS_paths_NOVEC_NOT (TCircuit *c,
                                CState init_state, CState final_state, const SampleCounter n_samples,
                myReal &sumR, myReal &sumI);
static bool IS_paths_NOVEC_T (TCircuit *c,
                              CState init_state, CState final_state, SampleCounter& samples2proc, bool& taskReady, bool& resAvailable, myReal& T_sumR, myReal& T_sumI, SampleCounter& processedSamples);
#endif

#ifdef NON_ZERO_PATHS
static bool IS_paths_NOVEC_kernel (TCircuit *c,
                CState init_state, CState final_state, const SampleCounter n_samples,
                myReal &sumR, myReal &sumI,
                                   pcg32& e,
                                   std::uniform_real_distribution<float>& d,
                                   SampleCounter& non_zero_paths);
#else
    static bool IS_paths_NOVEC_kernel (TCircuit *c,
                                       CState init_state, CState final_state, const SampleCounter n_samples,
                myReal &sumR, myReal &sumI,
                                       pcg32& e,
                                       std::uniform_real_distribution<float>& d);
#endif

#ifdef CONVERGENCE_STATS
bool IS_paths (TCircuit *c, CState init_state,
               CState final_state, const SampleCounter n_samples,
               myReal &estimateR, myReal &estimateI, std::vector<T_Stats>& stats, const int n_threads) {
#else
bool IS_paths (TCircuit *c, CState init_state,
               CState final_state, const SampleCounter n_samples,
                myReal &estimateR, myReal &estimateI, const int n_threads) {
#endif

    bool ret=true;
#ifdef NON_ZERO_PATHS
    SampleCounter non_zero_paths = 0;
#endif
    myReal sumR=0.f, sumI=0.f;
    SampleCounter n_ProcessedSamples = 0;

    
    if (n_threads<=1) {
#ifdef NON_ZERO_PATHS
        ret = IS_paths_NOVEC_NOT (c, init_state, final_state, n_samples, sumR, sumI, non_zero_paths);
#else
        ret = IS_paths_NOVEC_NOT (c, init_state, final_state, n_samples, sumR, sumI);
#endif
        
        double const dn_samples = (double) n_samples;
        estimateR = sumR / dn_samples;
        estimateI = sumI / dn_samples;

        n_ProcessedSamples = n_samples;
    }
    else {
        
        // MULTIPLE THREADS
        // MASTER WORKER
        const SampleCounter TASK_SIZE=n_samples/(8*n_threads);
        SampleCounter remaining_samples = n_samples;
        int busy_threads = 0;   // how many threads have a task

        std::vector<std::thread> threads;
#ifdef NON_ZERO_PATHS
        SampleCounter * l_NzeroP = new SampleCounter [n_threads];
#endif
        bool * taskReady = new bool [n_threads]; //false;
        bool * resAvailable = new bool [n_threads]; //false;
        // returned result from each thread
        myReal * T_sumR = new myReal [n_threads]; // 0.f
        myReal * T_sumI = new myReal [n_threads]; // 0.f

        SampleCounter * samples2proc = new SampleCounter [n_threads]; // = 0ull
        SampleCounter * processedSamples = new SampleCounter [n_threads]; // = 0ull

        bool * idleThreads = new bool[n_threads]; // false
        
#ifdef CONVERGENCE_STATS
        // https://www.geeksforgeeks.org/measure-execution-time-function-cpp/
        auto time_start = high_resolution_clock::now();
#endif

#ifdef DEBUG_THREAD
        fprintf (stderr, "MASTER : Creating %d threads... \n", n_threads); fflush (stderr);
#endif
        // Launch the threads
        // after creation they will be waiting on cv_m2s
        for (int t = 0; t < n_threads; t++) {
            taskReady[t] = false;
            resAvailable[t] = false;
            T_sumR[t] = 0.f;
            T_sumI[t] = 0.f;
            samples2proc[t] = 0;
            processedSamples[t] = 0;
            idleThreads[t] = false;
#ifdef NON_ZERO_PATHS
            l_NzeroP[t] = 0;
            threads.push_back(std::thread(IS_paths_NOVEC_T, c, init_state, final_state, std::ref(samples2proc[t]), std::ref(taskReady[t]), std::ref(resAvailable[t]), std::ref(T_sumR[t]), std::ref(T_sumI[t]), std::ref(processedSamples[t]), std::ref(l_NzeroP[t])));
#else
            threads.push_back(std::thread(IS_paths_NOVEC_T, c, init_state, final_state, std::ref(samples2proc[t]), std::ref(taskReady[t]), std::ref(resAvailable[t]), std::ref(T_sumR[t]), std::ref(T_sumI[t]), std::ref(processedSamples[t])));
#endif
        }
#ifdef DEBUG_THREAD
        fprintf (stderr, "MASTER : ... %d threads created! \n", n_threads); fflush (stderr);
#endif
        
        // Distribute 1st task to each thread
        for (int t = 0; (t < n_threads) && remaining_samples>0 ; t++) {
#ifdef DEBUG_THREAD
            fprintf (stderr, "MASTER : Locking mx_m2s to send task... \n"); fflush (stderr);
#endif
            { // scoped lock
                std::lock_guard<std::mutex> lk_m2s(mx_m2s);
                samples2proc[t] = ((remaining_samples < TASK_SIZE ) ? remaining_samples : TASK_SIZE);
                remaining_samples -= samples2proc[t];
                taskReady[t] = true;
            }  // lk_m2s terminates, thus mx_m2s is released by the lock destructor
            //  notify one thread
            cv_m2s.notify_one();
            busy_threads++;  // one thread got busy

#ifdef DEBUG_THREAD
            fprintf (stderr, "MASTER : mx_m2s unlocked and thread notified (%d busy_threads) \n", busy_threads); fflush (stderr);
#endif
        }
        
        // Threads have work now
        // while there are busy threads wait for them
        // and while there is work send new tasks
        while (busy_threads > 0) {
            // get the result
#ifdef DEBUG_THREAD
            fprintf (stderr, "MASTER : Locking mx_s2m to wait for result... \n"); fflush (stderr);
#endif
            { // scoped lock
                std::unique_lock<std::mutex> lk_s2m(mx_s2m);
                // wait on cv_s2m while resAvailable==false
                cv_s2m.wait(lk_s2m, [&]{bool ret=false; for (int t=0; t<n_threads && !ret; t++) ret= ret || resAvailable[t]; return (ret);});
                // mx_s2m is locked thus
                // loop through all resAvailable and retrieve the results of those that are == true
                for (int t=0 ; t < n_threads ; t++) {
                    if (resAvailable[t]) {
                        sumR += T_sumR[t];  // retrieve results
                        sumI += T_sumI[t];  // retrieve results
                        non_zero_paths += l_NzeroP[t];
                        resAvailable[t] = false;  // results read
                        idleThreads[t] = true;  // this thread is now idle
                        n_ProcessedSamples += processedSamples[t];
                        busy_threads--;
#ifdef CONVERGENCE_STATS
                        T_Stats stat;
                        stat.sumR = sumR;
                        stat.sumI = sumI;
                        stat.n_samples = n_ProcessedSamples;
                        stat.n_Paths = n_ProcessedSamples;
#ifdef NON_ZERO_PATHS
                        stat.n_nonZero_paths = non_zero_paths;
#endif

                        // https://www.geeksforgeeks.org/measure-execution-time-function-cpp/
                        auto time_now = high_resolution_clock::now();
                        auto duration = duration_cast<microseconds>(time_now - time_start);
                        
                        // To get the value of duration use the count()
                        // member function on the duration object
                        stat.time_us = duration.count();
                        stats.push_back (stat);
#endif
                    }
                }
            } // lk_s2m terminates, thus mx_s2m is released by the lock destructor
#ifdef DEBUG_THREAD
            fprintf (stderr, "MASTER : mx_s2m unlocked and results gotten (%d busy_threads) \n", busy_threads); fflush (stderr);
#endif
            
            // Free threads will now be waiting on cv_m2s
            if (remaining_samples > 0) {
#ifdef DEBUG_THREAD
                fprintf (stderr, "MASTER : Locking mx_m2s to send task... \n"); fflush (stderr);
#endif
                { // scoped lock
                    std::lock_guard<std::mutex> lk_m2s(mx_m2s);
                    // loop through all idleThreads and send a task to those that are == true
                    for (int t=0 ; t < n_threads && remaining_samples > 0 ; t++) {
                        if (idleThreads[t]) {
                            samples2proc[t] = ((remaining_samples < TASK_SIZE) ? remaining_samples : TASK_SIZE);
                            remaining_samples -= samples2proc[t];
                            taskReady[t] = true;
                            idleThreads[t] = false;
                            busy_threads++;
                        }
                    }
                }  // lk_m2s terminates, thus mx_m2s is released by the lock destructor
                //  notify all threads
                cv_m2s.notify_all();
#ifdef DEBUG_THREAD
                fprintf (stderr, "MASTER : mx_m2s unlocked and thread notified (%d busy_threads) \n", busy_threads); fflush (stderr);
#endif
            }
        }
            
        // notify all threads to terminate
#ifdef DEBUG_THREAD
        fprintf (stderr, "MASTER : Locking mx_m2s to send terminate... \n"); fflush (stderr);
#endif
        { // scoped lock
            std::lock_guard<std::mutex> lk_m2s(mx_m2s);
            terminate = true;
        }  // lk_m2s terminates, thus mx_m2s is released by the lock destructor
        //  notify all threads
        cv_m2s.notify_all();
#ifdef DEBUG_THREAD
        fprintf (stderr, "MASTER : mx_m2s unlocked and threads notified  \n"); fflush (stderr);
#endif
        
        for (auto &th : threads) {
            th.join();
        }
            
#ifdef NON_ZERO_PATHS
        /*for (int t=0 ; t< n_threads ; t++) {
            non_zero_paths += l_NzeroP[t];
        }*/
        delete[] l_NzeroP;
#endif
        delete[] idleThreads;
        delete[] taskReady;
        delete[] samples2proc;
        delete[] processedSamples;
        delete[] resAvailable;
        delete[] T_sumR;
        delete[] T_sumI;

        double const dn_samples = (double) n_samples;
        estimateR = sumR / dn_samples;
        estimateI = sumI / dn_samples;
    }

    fprintf (stdout, "Total samples: %llu\n", n_ProcessedSamples);
#ifdef NON_ZERO_PATHS
    fprintf (stdout, "Non zero paths: %llu\n", non_zero_paths);
#endif
    fprintf (stdout, "Total evaluated paths: %llu\n", n_samples);

    return ret;
}


#ifdef NON_ZERO_PATHS
static bool IS_paths_NOVEC_NOT (TCircuit *c,
                CState init_state, CState final_state, const SampleCounter n_samples,
                              myReal &sumR, myReal &sumI, SampleCounter& non_zero_paths) {
#else
    static bool IS_paths_NOVEC_NOT (TCircuit *c,
                    CState init_state, CState final_state, const SampleCounter n_samples,
                                  myReal &sumR, myReal &sumI) {
#endif
    // thread local random number generator (seeded by a local random device)
    // see https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2013/n3551.pdf
    std::random_device rdev{};
    thread_local pcg32 e{rdev()};
    std::uniform_real_distribution<float>d{0.0,1.0};  // uniform distribution in[0,1[ (myReal)
        
#ifdef NON_ZERO_PATHS
    return IS_paths_NOVEC_kernel(c, init_state, final_state, n_samples, sumR, sumI, e, d, non_zero_paths);
#else
    return IS_paths_NOVEC_kernel(c, init_state, final_state, n_samples, sumR, sumI, e, d);
#endif
        
}

    
#ifdef NON_ZERO_PATHS
static bool IS_paths_NOVEC_T (TCircuit *c,
                CState init_state, CState final_state, SampleCounter& samples2proc, bool& taskReady, bool& resAvailable, myReal& T_sumR, myReal& T_sumI, SampleCounter& processedSamples, SampleCounter& non_zero_paths) {
#else
static bool IS_paths_NOVEC_T (TCircuit *c,
                CState init_state, CState final_state, SampleCounter& samples2proc, bool& taskReady, bool& resAvailable, myReal& T_sumR, myReal& T_sumI, SampleCounter& processedSamples) {
#endif
    
    // thread local random number generator (seeded by a local random device)
    // see https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2013/n3551.pdf
    std::random_device rdev{};
    thread_local pcg32 e{rdev()};
    std::uniform_real_distribution<float>d{0.0,1.0};  // uniform distribution in[0,1[ (myReal)
        
    SampleCounter n_samples=0;
    myReal sumR, sumI;
    bool T_terminate=false, isTask = false;
        
    // wait for tasks or terminate on cv_m2s
    while (!T_terminate) {
#ifdef DEBUG_THREAD
        fprintf (stderr, "SLAVE : Locking mx_m2s to wait for task or terminate... \n"); fflush (stderr);
#endif
        { // scoped lock
            std::unique_lock<std::mutex> lk_m2s(mx_m2s);
            cv_m2s.wait(lk_m2s, [&]{return (taskReady == true  || terminate == true);});
            // mx_m2s is locked thus
            
            // is it a task ?
            if (taskReady) {
                isTask = true;
                n_samples = samples2proc;
                taskReady = false;
                samples2proc = 0;
#ifdef DEBUG_THREAD
                fprintf (stderr, "SLAVE : ... received task ... \n"); fflush (stderr);
#endif
            }
            else {  // terminate
                T_terminate = true;
#ifdef DEBUG_THREAD
                fprintf (stderr, "SLAVE : ... received terminate ... \n"); fflush (stderr);
#endif
            }
        } // lk_m2s terminates, thus mx_m2s is released by the lock destructor
#ifdef DEBUG_THREAD
        fprintf (stderr, "SLAVE : ... released mx_m2s \n"); fflush (stderr);
#endif
        
        // if this is a task then process it, send result and loop
        if (isTask) {
            isTask = false;
            sumR = sumI = 0.f;
            // process
#ifdef NON_ZERO_PATHS
            IS_paths_NOVEC_kernel(c, init_state, final_state, n_samples, sumR, sumI, e, d, non_zero_paths);
#else
            IS_paths_NOVEC_kernel(c, init_state, final_state, n_samples, sumR, sumI, e, d);
#endif
            
            // we have a result
            // send result to master who will be waiting on cv_s2m
#ifdef DEBUG_THREAD
            fprintf (stderr, "SLAVE : Locking mx_s2m to send result... \n"); fflush (stderr);
#endif
            { // scoped lock
                std::lock_guard<std::mutex> lk_s2m(mx_s2m);
                T_sumR = sumR;  // send results
                T_sumI = sumI;  // send results
                processedSamples = n_samples;
                resAvailable = true;  // there is a result
            }  // lk_s2m terminates, thus mx_s2m is released by the lock destructor
            //  notify one thread
            cv_s2m.notify_one();
#ifdef DEBUG_THREAD
            fprintf (stderr, "SLAVE : ... notified and released mx_s2m \n"); fflush (stderr);
#endif
        }  // if (isTask)
    }  // while (!T_terminate)
    
    return true;
}

    
#ifdef NON_ZERO_PATHS
static bool IS_paths_NOVEC_kernel (TCircuit *c,
                CState init_state, CState final_state, const SampleCounter n_samples,
                myReal &sumR, myReal &sumI,
                pcg32& e,
                std::uniform_real_distribution<float>& d,
                SampleCounter& non_zero_paths) {
#else
    static bool IS_paths_NOVEC_kernel (TCircuit *c,
                CState init_state, CState final_state, const SampleCounter n_samples,
                myReal &sumR, myReal &sumI,
                pcg32& e,
                std::uniform_real_distribution<float>& d) {
#endif
    
    SampleCounter s;   // sample counter
    int l;                  // layer counter
    
    myReal l_sumR=0.f, l_sumI=0.f;  // local summ accumulators for performance reasons
#ifdef NON_ZERO_PATHS
    SampleCounter l_non_zero_paths = 0;      // local counter for performance reasons
#endif
        
    const int L = c->size->num_layers;   // number of layers in the circuit
        
    // iteratively generate samples
    for (s=0; s<n_samples ; s++) {
        myReal path_pdf= 1.f;
        myReal path_wR=1.f, path_wI = 0.f;
        bool zero_power_transition = false;
        
        CState next_state=0, current_state = init_state;  // state before the next layer
        
#if defined(DEBUG) && !defined(__CSTATE_MP__)
        fprintf(stderr, "Sample: %llu out of %llu\n", s, n_samples);
#endif
        
        // generate the path by stochastically sampling each layer from l=0 to l=L-2
        // the last layer (l=L-1) is handled outside the 'for' loop since it is deterministically
        // connected to 'final_state'
        for (l=0 ; l< L-1 && !zero_power_transition ; l++) {
            myReal wR, wI, pdf;
            
#if defined(DEBUG) && !defined(__CSTATE_MP__)
                fprintf(stderr, "\tLayer: %d out of %d\n", l, L);
                fprintf(stderr, "\tCurrent state: %llu\n", current_state);
                fprintf(stderr, "\tLayer transition...\n");
#endif
            // get gate layer l
            TCircuitLayer *layer = &c->layers[l];
            
            // sample this layer for the current state,
            // returning the amplitude, pdf and next state
            pdf = layer_sample (layer, l, current_state, next_state, wR, wI, e, d);

#if defined(DEBUG) && !defined(__CSTATE_MP__)

                fprintf(stderr, "\tLayer w: %f + j %f\n", wR, wI);
                fprintf(stderr, "\tTransition p: %f\n", pdf);
                fprintf(stderr, "\tNext state: %llu\n", next_state);
#endif

            if (complex_abs_square(wR, wI) <= 0.f || pdf <= 0.f) {  // I believe this should never happen
                zero_power_transition = true;
#ifdef DEBUG
                    fprintf(stderr, "\tFinishing sample\n");
#endif
            }
            else {
                path_pdf *= pdf;
                complex_multiply (path_wR, path_wI, path_wR, path_wI, wR, wI);
#if defined(DEBUG) && !defined(__FLOAT_MP__)
                    fprintf(stderr, "\tPath w up to now: %f + j %f\n", path_wR, path_wI);
                    fprintf(stderr, "\tPath Transition p: %f\n", path_pdf);
#endif

            }
            
            current_state = next_state;
#ifdef DEBUG
            fprintf (stderr,"....\n");
#endif
        }  // layers 0 .. L-2 done
        
        // final layer (L-1)
        if (!zero_power_transition) { // path still contributes (I believer it always will)
            myReal wR, wI;
            
#if defined(DEBUG) && !defined(__CSTATE_MP__)
                fprintf(stderr, "\tLast Layer out of %d\n", L);
                fprintf(stderr, "\tCurrent state: %llu\n", current_state);
                fprintf(stderr, "\tFinal state: %llu\n", final_state);
                fprintf(stderr, "\tLayer deterministic transition...\n");
#endif
            
            TCircuitLayer *last_layer;
            /*fprintf (stderr, "VERIFFICATION OF ALL LAYERS\n");
            for (int lll=0 ; lll<L ; lll++){
                // get gate layer L-l
                last_layer = &c->layers[lll];

                fprintf (stderr, "Layer verification: %d out of %d with %d gates", lll, L, last_layer->num_gates);
                for (int gg=0 ; gg<last_layer->num_type_gates[0] ; gg++) {
                    fprintf (stderr, "> gate %d ", (((TGate1P0 *) last_layer->gates[0])[gg]).name);
                }
                fprintf(stderr, "\n");
            }*/
            
            //fprintf (stderr, "CALLING print_circuit...\n");
            //print_circuit(c);
            //fprintf (stderr, "... print_circuit DONE...\n");

            last_layer = &c->layers[L-1];

            // evaluate this layer amplitude thansitioning from the current state to final state
            layer_w (last_layer, L-1, current_state, final_state, wR, wI);

#if defined(DEBUG) && !defined(__FLOAT_MP__)
                fprintf(stderr, "\tLast Layer w: %f + j %f\n", wR, wI);
#endif

            if (complex_abs_square(wR, wI) <= 0.f) {  // This will happen frfequently
                zero_power_transition = true;
#ifdef DEBUG
                    fprintf(stderr, "\tFinishing sample\n");
#endif
            }
            else {
                complex_multiply (path_wR, path_wI, path_wR, path_wI, wR, wI);
#if defined(DEBUG) && !defined(__FLOAT_MP__)
                    fprintf(stderr, "\tPath final w : %f + j %f\n", path_wR, path_wI);
                    fprintf(stderr, "\tPath Transition p: %f\n", path_pdf);
#endif
            }

        }
        
        if (!zero_power_transition) {  // OK, count non_zero paths and accumulate
#ifdef NON_ZERO_PATHS
            l_non_zero_paths++;
#endif
            myReal const pdf_reciprocal = 1.f / path_pdf, auxZERO=0.f;
            myReal path_contR, path_contI;
            complex_multiply (path_contR, path_contI, path_wR, path_wI, pdf_reciprocal, auxZERO);
            
            // accumulate on local sums
            l_sumR += path_contR;
            l_sumI += path_contI;
#if defined(DEBUG) && !defined(__FLOAT_MP__)
                fprintf(stderr, "\tPath contribution : %f + j %f\n", path_contR, path_contI);
                fprintf(stderr, "\tSamples sum: %f + j %f\n", l_sumR, l_sumI);
#endif

        }
        // next sample
#ifdef DEBUG
            fprintf(stderr, "\n");
#endif

    } // end iterating over samples
    
#ifdef NON_ZERO_PATHS
    non_zero_paths = l_non_zero_paths;
#endif
    sumR=l_sumR;
    sumI=l_sumI;

    return true;
}

