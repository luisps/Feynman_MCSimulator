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

#include "complex.h"
#include "layer.hpp"
#include "PreProcessorSettings.h"

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

std::mutex mx_m2s, mx_s2m;
std::condition_variable cv_m2s, cv_s2m;

// Shared predicate variables
// Master to Slave
bool terminate = false;


#ifdef NON_ZERO_PATHS
static bool IS_paths_NOVEC_NOT (TCircuit *c,
                unsigned long long init_state, unsigned long long final_state, const unsigned long long n_samples,
                float &sumR, float &sumI, int& non_zero_paths);
static bool IS_paths_NOVEC_T (TCircuit *c,
                unsigned long long init_state, unsigned long long final_state, unsigned long long& samples2proc, bool& taskReady, bool& resAvailable, float& T_sumR, float& T_sumI, unsigned long long& processedSamples, int& non_zero_paths);
#else
static bool IS_paths_NOVEC_NOT (TCircuit *c,
                unsigned long long init_state, unsigned long long final_state, const unsigned long long n_samples,
                float &sumR, float &sumI);
static bool IS_paths_NOVEC_T (TCircuit *c,
                unsigned long long init_state, unsigned long long final_state, unsigned long long& samples2proc, bool& taskReady, bool& resAvailable, float& T_sumR, float& T_sumI, unsigned long long& processedSamples);
#endif

#ifdef NON_ZERO_PATHS
static bool IS_paths_NOVEC_kernel (TCircuit *c,
                unsigned long long init_state, unsigned long long final_state, const unsigned long long n_samples,
                float &sumR, float &sumI,
                std::default_random_engine& e, std::uniform_real_distribution<float>& d,
                                   int& non_zero_paths);
#else
    static bool IS_paths_NOVEC_kernel (TCircuit *c,
                unsigned long long init_state, unsigned long long final_state, const unsigned long long n_samples,
                float &sumR, float &sumI,
                                       std::default_random_engine& e, std::uniform_real_distribution<float>& d);
#endif

#ifdef CONVERGENCE_STATS
bool IS_paths (TCircuit *c, unsigned long long init_state,
                unsigned long long final_state, const unsigned long long n_samples,
               float &estimateR, float &estimateI, std::vector<T_Stats>& stats, const int n_threads) {
#else
bool IS_paths (TCircuit *c, unsigned long long init_state,
                unsigned long long final_state, const unsigned long long n_samples,
                float &estimateR, float &estimateI, const int n_threads) {
#endif

    bool ret=true;
#ifdef NON_ZERO_PATHS
    int non_zero_paths = 0;
#endif
    float sumR=0.f, sumI=0.f;
    unsigned long long n_ProcessedSamples = 0ull;

    
    if (n_threads<=1) {
#ifdef NON_ZERO_PATHS
        ret = IS_paths_NOVEC_NOT (c, init_state, final_state, n_samples, sumR, sumI, non_zero_paths);
#else
        ret = IS_paths_NOVEC_NOT (c, init_state, final_state, n_samples, sumR, sumI);
#endif
        estimateR = sumR / ((float)n_samples);
        estimateI = sumI / ((float)n_samples);
    }
    else {
        
        // MULTIPLE THREADS
        // MASTER WORKER
        const long long TASK_SIZE=1024ll;
        long long remaining_samples = (long long)n_samples;
        int busy_threads = 0;   // how many threads have a task

        std::vector<std::thread> threads;
#ifdef NON_ZERO_PATHS
        int * l_NzeroP = new int [n_threads];
#endif
        bool * taskReady = new bool [n_threads]; //false;
        bool * resAvailable = new bool [n_threads]; //false;
        // returned result from each thread
        float * T_sumR = new float [n_threads]; // 0.f
        float * T_sumI = new float [n_threads]; // 0.f

        unsigned long long * samples2proc = new unsigned long long [n_threads]; // = 0ull
        unsigned long long * processedSamples = new unsigned long long [n_threads]; // = 0ull

        bool * idleThreads = new bool[n_threads]; // false

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
            samples2proc[t] = 0ull;
            processedSamples[t] = 0ull;
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
        for (int t = 0; (t < n_threads) && remaining_samples>0ll ; t++) {
#ifdef DEBUG_THREAD
            fprintf (stderr, "MASTER : Locking mx_m2s to send task... \n"); fflush (stderr);
#endif
            { // scoped lock
                std::lock_guard<std::mutex> lk_m2s(mx_m2s);
                samples2proc[t] = ((remaining_samples-TASK_SIZE <0ull) ? remaining_samples : TASK_SIZE);
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
                        resAvailable[t] = false;  // results read
                        idleThreads[t] = true;  // this thread is now idle
                        n_ProcessedSamples += processedSamples[t];
                        busy_threads--;
#ifdef CONVERGENCE_STATS
                        T_Stats stat;
                        stat.sumR = sumR;
                        stat.sumI = sumI;
                        stat.n_samples = n_ProcessedSamples;
                        stats.push_back (stat);
#endif
                    }
                }
            } // lk_s2m terminates, thus mx_s2m is released by the lock destructor
#ifdef DEBUG_THREAD
            fprintf (stderr, "MASTER : mx_s2m unlocked and results gotten (%d busy_threads) \n", busy_threads); fflush (stderr);
#endif
            
            // Free threads will now be waiting on cv_m2s
            if (remaining_samples > 0ll) {
#ifdef DEBUG_THREAD
                fprintf (stderr, "MASTER : Locking mx_m2s to send task... \n"); fflush (stderr);
#endif
                { // scoped lock
                    std::lock_guard<std::mutex> lk_m2s(mx_m2s);
                    // loop through all idleThreads and send a task to those that are == true
                    for (int t=0 ; t < n_threads && remaining_samples > 0ll ; t++) {
                        if (idleThreads[t]) {
                            samples2proc[t] = ((remaining_samples-TASK_SIZE <0ull) ? remaining_samples : TASK_SIZE);
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
        for (int t=0 ; t< n_threads ; t++) {
            non_zero_paths += l_NzeroP[t];
        }
        delete[] l_NzeroP;
#endif
        delete[] idleThreads;
        delete[] taskReady;
        delete[] samples2proc;
        delete[] processedSamples;
        delete[] resAvailable;
        delete[] T_sumR;
        delete[] T_sumI;

        estimateR = sumR / n_samples;
        estimateI = sumI / n_samples;
    }

#ifdef NON_ZERO_PATHS
    fprintf (stderr, "Non zero paths: %d\n", non_zero_paths);
#endif
    
    return ret;
}


#ifdef NON_ZERO_PATHS
static bool IS_paths_NOVEC_NOT (TCircuit *c,
                unsigned long long init_state, unsigned long long final_state, const unsigned long long n_samples,
                              float &sumR, float &sumI, int& non_zero_paths) {
#else
    static bool IS_paths_NOVEC_NOT (TCircuit *c,
                    unsigned long long init_state, unsigned long long final_state, const unsigned long long n_samples,
                                  float &sumR, float &sumI) {
#endif
    // thread local random number generator (seeded by a local random device)
    // see https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2013/n3551.pdf
    std::random_device rdev{};
    thread_local std::default_random_engine e{rdev()};
    std::uniform_real_distribution<float>d{0.0,1.0};  // uniform distribution in[0,1[ (float)
        
#ifdef NON_ZERO_PATHS
    return IS_paths_NOVEC_kernel(c, init_state, final_state, n_samples, sumR, sumI, e, d, non_zero_paths);
#else
    return IS_paths_NOVEC_kernel(c, init_state, final_state, n_samples, sumR, sumI, e, d);
#endif
        
}

    
#ifdef NON_ZERO_PATHS
static bool IS_paths_NOVEC_T (TCircuit *c,
                unsigned long long init_state, unsigned long long final_state, unsigned long long& samples2proc, bool& taskReady, bool& resAvailable, float& T_sumR, float& T_sumI, unsigned long long& processedSamples, int& non_zero_paths) {
#else
static bool IS_paths_NOVEC_T (TCircuit *c,
                unsigned long long init_state, unsigned long long final_state, unsigned long long& samples2proc, bool& taskReady, bool& resAvailable, float& T_sumR, float& T_sumI, unsigned long long& processedSamples) {
#endif
    
    // thread local random number generator (seeded by a local random device)
    // see https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2013/n3551.pdf
    std::random_device rdev{};
    thread_local std::default_random_engine e{rdev()};
    std::uniform_real_distribution<float>d{0.0,1.0};  // uniform distribution in[0,1[ (float)
        
    unsigned long long n_samples=0ull;
    float sumR, sumI;
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
                samples2proc = 0ll;
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
                unsigned long long init_state, unsigned long long final_state, const unsigned long long n_samples,
                float &sumR, float &sumI,
                std::default_random_engine& e, std::uniform_real_distribution<float>& d,
                int& non_zero_paths) {
#else
    static bool IS_paths_NOVEC_kernel (TCircuit *c,
                unsigned long long init_state, unsigned long long final_state, const unsigned long long n_samples,
                float &sumR, float &sumI,
                std::default_random_engine& e, std::uniform_real_distribution<float>& d) {
#endif
    
    unsigned long long s;   // sample counter
    int l;                  // layer counter
    
    float l_sumR=0.f, l_sumI=0.f;  // local summ accumulators for performance reasons
#ifdef NON_ZERO_PATHS
    int l_non_zero_paths = 0;      // local counter for performance reasons
#endif
        
    const int L = c->size->num_layers;   // number of layers in the circuit
        
    // iteratively generate samples
    for (s=0ull; s<n_samples ; s++) {
        float path_pdf= 1.f;
        float path_wR=1.f, path_wI = 0.f;
        bool zero_power_transition = false;
        
        unsigned long long next_state=0ull, current_state = init_state;  // state before the next layer
        
#ifdef DEBUG
        fprintf(stderr, "Sample: %llu out of %llu\n", s, n_samples);
#endif
        
        // generate the path by stochastically sampling each layer from l=0 to l=L-2
        // the last layer (l=L-1) is handled outside the 'for' loop since it is deterministically
        // connected to 'final_state'
        for (l=0 ; l< L-1 && !zero_power_transition ; l++) {
            float wR, wI, pdf;
            
#ifdef DEBUG
                fprintf(stderr, "\tLayer: %d out of %d\n", l, L);
                fprintf(stderr, "\tCurrent state: %llu\n", current_state);
                fprintf(stderr, "\tLayer transition...\n");
#endif
            // get gate layer l
            TCircuitLayer *layer = &c->layers[l];
            
            // sample this layer for the current state,
            // returning the amplitude, pdf and next state
            pdf = layer_sample (layer, l, current_state, next_state, wR, wI, e, d);

#ifdef DEBUG

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
#ifdef DEBUG
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
            float wR, wI;
            
#ifdef DEBUG
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

#ifdef DEBUG
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
#ifdef DEBUG
                    fprintf(stderr, "\tPath final w : %f + j %f\n", path_wR, path_wI);
                    fprintf(stderr, "\tPath Transition p: %f\n", path_pdf);
#endif
            }

        }
        
        if (!zero_power_transition) {  // OK, count non_zero paths and accumulate
#ifdef NON_ZERO_PATHS
            l_non_zero_paths++;
#endif
            float pdf_reciprocal = 1.f / path_pdf;
            float path_contR, path_contI;
            complex_multiply (path_contR, path_contI, path_wR, path_wI, pdf_reciprocal, 0.f);
            
            // accumulate on local sums
            l_sumR += path_contR;
            l_sumI += path_contI;
#ifdef DEBUG
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

