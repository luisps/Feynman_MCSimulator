//
//  BD_paths.cpp
//  Feynman_MCSimulator
//
//  Created by Luis Paulo Santos on 22/09/2023.
//

#include "BD_paths.hpp"

#include <thread>
#include <mutex>
#include <condition_variable>
#include <vector>
#include <algorithm>
#include <random>
#include "pcg_random.hpp"

#include "complex.h"
#include "layer.hpp"
#include "PreProcessorSettings.h"
#include "path.h"
#include "PathVertex.h"

// For time stats
// https://www.geeksforgeeks.org/measure-execution-time-function-cpp/
#include <chrono>
using namespace std::chrono;

static std::mutex mx_m2s, mx_s2m;
static std::condition_variable cv_m2s, cv_s2m;

// Shared predicate variables
// Master to Slave
static bool terminate = false;


#ifdef NON_ZERO_PATHS
static bool BD_paths_NOVEC_NOT (TCircuit *c,
                unsigned long long init_state, unsigned long long final_state, const unsigned long long n_samples,
                myReal &sumR, myReal &sumI, unsigned long long& n_Paths, const bool MIS, int& non_zero_paths);
static bool BD_paths_NOVEC_T (TCircuit *c,
                unsigned long long init_state, unsigned long long final_state, unsigned long long& samples2proc, bool& taskReady, bool& resAvailable, myReal& T_sumR, myReal& T_sumI, unsigned long long& n_Paths, unsigned long long& processedSamples, const bool MIS, int& non_zero_paths);
#else
static bool BD_paths_NOVEC_NOT (TCircuit *c,
                unsigned long long init_state, unsigned long long final_state, const unsigned long long n_samples,
                myReal &sumR, myReal &sumI, unsigned long long& n_Paths, const bool MIS);
static bool BD_paths_NOVEC_T (TCircuit *c,
                unsigned long long init_state, unsigned long long final_state, unsigned long long& samples2proc, bool& taskReady, bool& resAvailable, myReal& T_sumR, myReal& T_sumI, unsigned long long& n_Paths, unsigned long long& processedSamples, const bool MIS);
#endif

#ifdef NON_ZERO_PATHS
static bool BD_paths_NOVEC_kernel (TCircuit *c,
                unsigned long long init_state, unsigned long long final_state, const unsigned long long n_samples,
                myReal &sumR, myReal &sumI,
                unsigned long long& n_Paths, const bool MIS,
                           pcg32& e,
                                   std::uniform_real_distribution<myReal>& d,
                                   int& non_zero_paths);
#else
    static bool BD_paths_NOVEC_kernel (TCircuit *c,
                unsigned long long init_state, unsigned long long final_state, const unsigned long long n_samples,
                myReal &sumR, myReal &sumI,
                unsigned long long& n_Paths, const bool MIS,
                           pcg32& e,
                                       std::uniform_real_distribution<myReal>& d);
#endif

static void _ForwardPath (TCircuit *c, unsigned long long init_state,
                          PathVertexVector& Fpath,
                           pcg32& e,
                          std::uniform_real_distribution<myReal>& d);

static void _BackwardPath (TCircuit *c, unsigned long long final_state,
                              PathVertexVector& Bpath,
                           pcg32& e,
                           std::uniform_real_distribution<myReal>& d);

static myReal _ConnectPath (TCircuitLayer* layer, const int l,
                     unsigned long long currentState,
                       PathVertexVector Fpath,
                       unsigned long long nextState,
                        PathVertexVector Bpath,
                       myReal& wR, myReal& wI);

static myReal _ConnectPathMIS (TCircuitLayer* layer, const int l, const int num_layers,
                            PathVertexVector Fpath,
                            PathVertexVector Bpath, 
                            myReal probBuffer[], myReal probFBuffer[],myReal probBBuffer[],
                            myReal& wR, myReal& wI);

#ifdef CONVERGENCE_STATS
bool BD_paths (TCircuit *c, unsigned long long init_state,
                unsigned long long final_state, const unsigned long long n_samples,
               myReal &estimateR, myReal &estimateI, std::vector<T_Stats>& stats, const int n_threads, const bool MIS) {
#else
bool BD_paths (TCircuit *c, unsigned long long init_state,
                unsigned long long final_state, const unsigned long long n_samples,
                myReal &estimateR, myReal &estimateI, const int n_threads, const bool MIS) {
#endif

    bool ret=true;
#ifdef NON_ZERO_PATHS
    int non_zero_paths = 0;
#endif
    myReal sumR=0.f, sumI=0.f;
    unsigned long long n_ProcessedSamples = 0ull;
    unsigned long long n_Paths = 0ull;

    
    if (n_threads<=1) {
#ifdef NON_ZERO_PATHS
        ret = BD_paths_NOVEC_NOT (c, init_state, final_state, n_samples, sumR, sumI, n_Paths, MIS, non_zero_paths);
#else
        ret = BD_paths_NOVEC_NOT (c, init_state, final_state, n_samples, sumR, sumI, n_Paths, MIS);
#endif
        estimateR = sumR / ((myReal)n_Paths);
        estimateI = sumI / ((myReal)n_Paths);
        n_ProcessedSamples = n_samples;
    }
    else {
        
        // MULTIPLE THREADS
        // MASTER WORKER
        const long long TASK_SIZE=n_samples/(8*n_threads);;
        long long remaining_samples = (long long)n_samples;
        int busy_threads = 0;   // how many threads have a task

        std::vector<std::thread> threads;
#ifdef NON_ZERO_PATHS
        int * l_NzeroP = new int [n_threads];
#endif
        bool * taskReady = new bool [n_threads]; //false;
        bool * resAvailable = new bool [n_threads]; //false;
        // returned result from each thread
        myReal * T_sumR = new myReal [n_threads]; // 0.f
        myReal * T_sumI = new myReal [n_threads]; // 0.f

        unsigned long long * samples2proc = new unsigned long long [n_threads]; // = 0ull
        unsigned long long * processedSamples = new unsigned long long [n_threads]; // = 0ull
        unsigned long long * nT_Paths = new unsigned long long [n_threads]; // = 0ull

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
            samples2proc[t] = 0ull;
            processedSamples[t] = 0ull;
            nT_Paths[t] = 0ull;
            idleThreads[t] = false;
#ifdef NON_ZERO_PATHS
            l_NzeroP[t] = 0;
            threads.push_back(std::thread(BD_paths_NOVEC_T, c, init_state, final_state, std::ref(samples2proc[t]), std::ref(taskReady[t]), std::ref(resAvailable[t]), std::ref(T_sumR[t]), std::ref(T_sumI[t]), std::ref(processedSamples[t]), std::ref(nT_Paths[t]), MIS, std::ref(l_NzeroP[t])));
#else
            threads.push_back(std::thread(BD_paths_NOVEC_T, c, init_state, final_state, std::ref(samples2proc[t]), std::ref(taskReady[t]), std::ref(resAvailable[t]), std::ref(T_sumR[t]), std::ref(T_sumI[t]), std::ref(processedSamples[t]), std::ref(nT_Paths[t]), MIS));
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
                samples2proc[t] = ((remaining_samples < TASK_SIZE) ? remaining_samples : TASK_SIZE);
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
                        n_Paths += nT_Paths[t];
                        busy_threads--;
#ifdef CONVERGENCE_STATS
                        T_Stats stat;
                        stat.sumR = sumR;
                        stat.sumI = sumI;
                        stat.n_samples = n_ProcessedSamples;
                        stat.n_Paths = n_Paths;
                        // https://www.geeksforgeeks.org/measure-execution-time-function-cpp/
                        auto time_now = high_resolution_clock::now();
                        auto duration = duration_cast<microseconds>(time_now - time_start);
                        
                        // To get the value of duration use the count()
                        // member function on the duration object
                        stat.time_us = duration.count();
                        stats.push_back (stat);
#endif
#ifdef DEBUG_THREAD
                        fprintf (stderr, "MASTER : mx_s2mlocked: thread %d sent %f + i %f for %llu paths \n", t, T_sumR[t], T_sumI[t], nT_Paths[t]); fflush  (stderr);
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

        estimateR = sumR / n_Paths;
        estimateI = sumI / n_Paths;
    }

    fprintf (stdout, "Total samples: %llu\n", n_ProcessedSamples);
#ifdef NON_ZERO_PATHS
    fprintf (stdout, "Non zero paths: %d\n", non_zero_paths);
#endif
    fprintf (stdout, "Total evaluated paths: %llu\n", n_Paths);

    return ret;
}

#ifdef NON_ZERO_PATHS
static bool BD_paths_NOVEC_NOT (TCircuit *c,
                unsigned long long init_state, unsigned long long final_state, const unsigned long long n_samples,
                              myReal &sumR, myReal &sumI,
                                unsigned long long& n_Paths, const bool MIS,
                                int& non_zero_paths) {
#else
    static bool BD_paths_NOVEC_NOT (TCircuit *c,
                    unsigned long long init_state, unsigned long long final_state, const unsigned long long n_samples,
                                  myReal &sumR, myReal &sumI,
                                    unsigned long long& n_Paths, const bool MIS) {
#endif
    // thread local random number generator (seeded by a local random device)
    // see https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2013/n3551.pdf
    std::random_device rdev{};
    thread_local pcg32 e{rdev()};
    std::uniform_real_distribution<myReal>d{0.0,1.0};  // uniform distribution in[0,1[ (myReal)
        
#ifdef NON_ZERO_PATHS
    return BD_paths_NOVEC_kernel(c, init_state, final_state, n_samples, sumR, sumI, n_Paths, MIS, e, d, non_zero_paths);
#else
    return BD_paths_NOVEC_kernel(c, init_state, final_state, n_samples, sumR, sumI, n_Paths, MIS, e, d);
#endif
        
}

#ifdef NON_ZERO_PATHS
static bool BD_paths_NOVEC_T (TCircuit *c,
                unsigned long long init_state, unsigned long long final_state, unsigned long long& samples2proc, bool& taskReady, bool& resAvailable, myReal& T_sumR, myReal& T_sumI, unsigned long long& processedSamples,
                              unsigned long long& n_Paths, const bool MIS,
                              int& non_zero_paths) {
#else
static bool BD_paths_NOVEC_T (TCircuit *c,
                unsigned long long init_state, unsigned long long final_state, unsigned long long& samples2proc, bool& taskReady, bool& resAvailable, myReal& T_sumR, myReal& T_sumI, unsigned long long& processedSamples,
                              unsigned long long& n_Paths, const bool MIS) {
#endif
    
    // thread local random number generator (seeded by a local random device)
    // see https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2013/n3551.pdf
    std::random_device rdev{};
    thread_local pcg32 e{rdev()};
    std::uniform_real_distribution<myReal>d{0.0,1.0};  // uniform distribution in[0,1[ (myReal)
        
    unsigned long long n_samples=0ull;
    myReal sumR, sumI;
    unsigned long long Tn_Paths;
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
            Tn_Paths = 0ull;
            // process
#ifdef NON_ZERO_PATHS
            BD_paths_NOVEC_kernel(c, init_state, final_state, n_samples, sumR, sumI, Tn_Paths, MIS, e, d, non_zero_paths);
#else
            BD_paths_NOVEC_kernel(c, init_state, final_state, n_samples, sumR, sumI, Tn_Paths, MIS, e, d);
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
                n_Paths = Tn_Paths;
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
static bool BD_paths_NOVEC_kernel (TCircuit *c,
                unsigned long long init_state, unsigned long long final_state, const unsigned long long n_samples,
                myReal &sumR, myReal &sumI,
                unsigned long long& n_Paths, const bool MIS,
                pcg32& e,
                                   std::uniform_real_distribution<myReal>& d,
                                   int& non_zero_paths) {
#else
static bool BD_paths_NOVEC_kernel (TCircuit *c,
                                    unsigned long long init_state, unsigned long long final_state,
                                   const unsigned long long n_samples,
                                    myReal &sumR, myReal &sumI,
                                   unsigned long long& n_Paths, const bool MIS,
                                   pcg32& e,
                                   std::uniform_real_distribution<myReal>& d) {
#endif
    unsigned long long s;   // sample counter
    int l;                  // layer counter
    unsigned long long l_n_Paths = 0ull;   // nbr of paths (local counter)

    myReal l_sumR=0.f, l_sumI=0.f;  // local summ accumulators for performance reasons
#ifdef NON_ZERO_PATHS
    int l_non_zero_paths = 0;      // local counter for performance reasons
#endif
    
    const int L = c->size->num_layers;   // number of layers in the circuit
    
    PathVertexVector Fpath(L), Bpath(L+1);
    myReal *ProbBuffer = new myReal[L];
    myReal *ProbFBuffer = new myReal[L];
    myReal *ProbBBuffer = new myReal[L];

    // iteratively generate samples
    for (s=0ull; s<n_samples ; s++) {
        
#ifdef DEBUG
        fprintf (stderr, "Sample %llu:\n", s);
#endif
        
        Fpath.clear();
        _ForwardPath (c, init_state, Fpath, e, d);

#ifdef DEBUG
/*        //print_path(Fpath);
        fprintf (stderr, "\n\t\tFORWARD\n\n ");
        fprintf (stderr, "%llu ", Fpath[0]);
        for (l=1 ; l< L ; l++) {
            fprintf (stderr, "-> %llu ", Fpath[l]);
        }
        fprintf (stderr, "||\n");
        fprintf (stderr, "\n");
        fprintf (stderr, "%.5f + i %.5f ", Fpath_wR[0], Fpath_wI[0]);
        for (l=1 ; l< L ; l++) {
            fprintf (stderr, "-> %.5f + i %.5f ", Fpath_wR[l], Fpath_wI[l]);
        }
        fprintf (stderr, "||\n");
        fprintf (stderr, "\n");
        fprintf (stderr, "%.5f ", Fpath_prob[0]);
        for (l=1 ; l< L ; l++) {
            fprintf (stderr, "-> %.5f ", Fpath_prob[l]);
        }
        §fprintf (stderr, "||\n"); */
#endif
        Bpath.clear();
        _BackwardPath (c, final_state, Bpath, e, d);

#ifdef DEBUG
        //print_path(Bpath);
       /* fprintf (stderr, "\n\t\tBACKWARD\n\n ");
        fprintf (stderr, "||\t%llu ", Bpath[1]);
        for (l=2 ; l<= L ; l++) {
            fprintf (stderr, "<- %llu ", Bpath[l]);
        }
        fprintf (stderr, "\n");
        fprintf (stderr, "\n");
        fprintf (stderr, "||\t%.5f + i %.5f ", Bpath_wR[1], Bpath_wI[1]);
        for (l=2 ; l<= L ; l++) {
            fprintf (stderr, "<- %.5f + i %.5f ", Bpath_wR[l], Bpath_wI[l]);
        }
        fprintf (stderr, "\n");
        fprintf (stderr, "\n");
        fprintf (stderr, "||\t%.5f ", Bpath_prob[1]);
        for (l=2 ; l<= L ; l++) {
            fprintf (stderr, "<- %.5f ", Bpath_prob[l]);
        }
        fprintf (stderr, "\n");
        fprintf (stderr, "\n");*/
#endif
        
        for (l=0 ; l < L ; l++ ) {  // generate all L different paths
            myReal prob, wR, wI;
            TCircuitLayer *layer = &c->layers[l];
            const unsigned long long currentState = Fpath.data[l].state;
            const unsigned long long nextState = Bpath.data[l+1].state;

#ifdef DEBUG
/*            fprintf (stderr, "Connect state %llu to state %llu through layer %d:\n", currentState, nextState, l);
            fprintf (stderr, "Connecting path is:\n%llu ", Fpath[0]);
            for (int ll=1; ll <=l ; ll++) {
                fprintf (stderr, " -> %llu ", Fpath[ll]);
            }
            fprintf (stderr, " <-> ");
            fprintf (stderr, " %llu ", Bpath[l+1]);
            for (int ll=l+2; ll <=L ; ll++) {
                fprintf (stderr, " <- %llu ", Bpath[ll]);
            }
             fprintf (stderr, "\n");*/
#endif
            if (MIS) {
                // connect deterministically fPath to bPath through layer l using MIS
                prob = _ConnectPathMIS (layer, l, L,
                                     Fpath, Bpath, ProbBuffer, ProbFBuffer, ProbBBuffer,
                                            wR, wI);
                
                // prob is in fact the reciprocal of the MIS weight.
            }
            else {
                // connect deterministically fPath to bPath through layer l
                prob = _ConnectPath (layer, l,
                                 currentState, Fpath,
                                 nextState, Bpath,
                                 wR, wI);
            }
#ifdef DEBUG
            fprintf (stderr, "\tw = %.5f + i %.5f \t prob = %.5f\n\n", wR, wI, prob);
#endif
            if ((complex_abs_square(wR, wI) > 0.f) && prob > 0.f) {  // non_zero_path
                myReal path_throughputR = wR / prob;
                myReal path_throughputI = wI / prob;
                l_sumR += path_throughputR;
                l_sumI += path_throughputI;
#ifdef NON_ZERO_PATHS
                //fprintf (stderr, "\tw = %e + i %e \t prob = %e\n", wR, wI, prob);
                l_non_zero_paths++;
#endif
            }
            l_n_Paths++;
        }  // iterate over layers
    }   // iterate over samples
    
#ifdef NON_ZERO_PATHS
    non_zero_paths = l_non_zero_paths;
#endif
    sumR=l_sumR;
    sumI=l_sumI;
    n_Paths = l_n_Paths;
    
    return true;
}
    
static void _ForwardPath (TCircuit *c, unsigned long long init_state,
                          PathVertexVector& Fpath,
                          pcg32& e,
                          std::uniform_real_distribution<myReal>& d) {
        
        
    unsigned long long current_state, next_state;
    int ndx;
    ndx = Fpath.append();
    // The forward path 1st node is the initial state
    Fpath.data[ndx].state = init_state;
    Fpath.data[ndx].wR = 1.0;
    Fpath.data[ndx].wI = 0.0;
    Fpath.data[ndx].PwR = 1.0;
    Fpath.data[ndx].PwI = 0.0;
    Fpath.data[ndx].prob = 1.0;
    Fpath.data[ndx].Pprob = 1.0;
    current_state = init_state;
        
        // Evolve through all gates' layes (starting at 0) , except the last one
        // The last layer will be later on deterministically connected to the final state
        for (int l=0 ; l<c->size->num_layers-1 ; l++) {
            myReal l_wR, l_wI, l_prob, auxR, auxI;
            
            TCircuitLayer *layer = &(c->layers[l]);
            
            // sample this layer
            l_prob = layer_sample(layer, l, current_state, next_state, l_wR, l_wI, e, d);
            
            // add to Fpath
            ndx = Fpath.append();
            Fpath.data[ndx].state = next_state;
            Fpath.data[ndx].wR = l_wR;
            Fpath.data[ndx].wI = l_wI;
            complex_multiply (auxR, auxI, Fpath.data[ndx-1].PwR, Fpath.data[ndx-1].PwI, l_wR, l_wI);
            Fpath.data[ndx].PwR = auxR;
            Fpath.data[ndx].PwI = auxI;
            Fpath.data[ndx].prob = l_prob;
            Fpath.data[ndx].Pprob = Fpath.data[ndx-1].Pprob * l_prob;

            current_state = next_state;
        }
        return;
    }

    
static void _BackwardPath (TCircuit *c, unsigned long long final_state,
                            PathVertexVector& Bpath,
                           pcg32& e,
                           std::uniform_real_distribution<myReal>& d) {
        
    unsigned long long current_state, next_state;
    int ndx;
    // The backward path 1st node is the final state: to revert later
    ndx = Bpath.prepend();
    Bpath.data[ndx].state = final_state;
    Bpath.data[ndx].wR = 1.0;
    Bpath.data[ndx].wI = 0.0;
    Bpath.data[ndx].PwR = 1.0;
    Bpath.data[ndx].PwI = 0.0;
    Bpath.data[ndx].prob= 1.0;
    Bpath.data[ndx].Pprob = 1.0;
        
    current_state = final_state;
        
    // Evolve through all gates' layes (starting at the last one (num_lkayers-1)) , except the first one
    // The first layer will be later on deterministically connected to the initial state
    for (int l=c->size->num_layers-1 ; l>0 ; l--) {
        myReal l_wR, l_wI, l_prob, auxR, auxI;
            
        TCircuitLayer *layer = &(c->layers[l]);
            
        // sample this layer
        l_prob = layer_sample(layer, l, current_state, next_state, l_wR, l_wI, e, d, false);
            
        // add to Bpath
        ndx = Bpath.prepend();
        Bpath.data[ndx].state = next_state;
        Bpath.data[ndx].wR = l_wR;
        Bpath.data[ndx].wI = l_wI;
        complex_multiply (auxR, auxI, Bpath.data[ndx+1].PwR, Bpath.data[ndx+1].PwI, l_wR, l_wI);
        Bpath.data[ndx].PwR = auxR;
        Bpath.data[ndx].PwI = auxI;
        Bpath.data[ndx].prob= l_prob;
        Bpath.data[ndx].Pprob = Bpath.data[ndx+1].Pprob * l_prob;
            
        current_state = next_state;
    }
    // to make sure that the element index [0] on the backward path means nothing add 0. to the vectors
    ndx = Bpath.prepend();
    Bpath.data[ndx].state = 0ull;
    Bpath.data[ndx].wR = 0.0;
    Bpath.data[ndx].wI = 0.0;
    Bpath.data[ndx].PwR = 0.0;
    Bpath.data[ndx].PwI = 0.0;
    Bpath.data[ndx].prob= 0.0;
    Bpath.data[ndx].Pprob = 0.0;
        
    return;
}
    
static myReal _ConnectPathMIS (TCircuitLayer* layer, const int l, const int L,
                                PathVertexVector Fpath,
                                PathVertexVector Bpath,
                                myReal prob_Buffer[],myReal probFor_Buffer[],myReal probBack_Buffer[],
                                myReal& wR, myReal& wI) {
    myReal lwR, lwI;
    unsigned long long currentState, nextState;
    int i;
            
    currentState = Fpath.data[l].state;
    nextState = Bpath.data[l+1].state;
    const myReal prob = layer_w_prob (layer, l, currentState, nextState, lwR, lwI);
            
    #ifdef DEBUG
        /*fprintf (stderr, "\nCONNECT_MIS\n");
        fprintf (stderr, "Fsegment Pw = %.5f + i %.5f\n", Fpath_PwR, Fpath_PwI);
        fprintf (stderr, "C layer w = %e + i %e\n", lwR, lwI);
        fprintf (stderr, "Bsegment Pw = %.5f + i %.5f\n", Bpath_PwR, Bpath_PwI);
        fprintf (stderr, "C prob = %e \n", prob);*/
    #endif
            
            
    // compute path w
    complex_multiply(lwR, lwI, lwR, lwI, Fpath.data[l].PwR, Fpath.data[l].PwI);
    complex_multiply(lwR, lwI, lwR, lwI, Bpath.data[l+1].PwR, Bpath.data[l+1].PwI);
            
    if (complex_abs_square(lwR, lwI) <=0 || prob <= 0) { // terminate path
        wR = 0.f;
        wI = 0.f;
        return 0.f;
    }

    // store each path segment prob as if it had been generated stochastically
    // i.e. without a det connection at layer l
    for (i=1 ; i<=l ; i++) {
        prob_Buffer[i-1] = Fpath.data[i].prob;
    }
    prob_Buffer[l] = prob;
    for (i=l+1 ; i<L ; i++) {
        prob_Buffer[i] = Bpath.data[i].prob;
    }

    // compute forward and backward products
    probFor_Buffer[0] = prob_Buffer[0];
    probBack_Buffer[L-1] = prob_Buffer[L-1];
    for (i=1 ; i<=L-1 ; i++) {
        probFor_Buffer[i] = probFor_Buffer[i-1] * prob_Buffer[i];
        probBack_Buffer[L-i-1] = probBack_Buffer[L-i] * prob_Buffer[L-i-1];
    }

    #ifdef DEBUG
        fprintf (stderr, "Path w = %.5f + i %.5f\n", lwR, lwI);
        fprintf (stderr, "Connection prob = %e \n", prob);
    #endif

    // compute MIS_weight as the sum of
    // the probabilities of the num_layers alternative deterministic connections
    myReal MIS_weight_reciprocal = 0.f;
    for (int det_connect=0 ; det_connect < L ; det_connect++) {
        myReal path_prob = 1.;
        if (det_connect > 0) path_prob = probFor_Buffer[det_connect-1];
        if (det_connect < (L-1)) path_prob *= probBack_Buffer[det_connect+1];

        MIS_weight_reciprocal += path_prob;

        #ifdef DEBUG
            fprintf (stderr, "(l=%d) Det connection layer %d : NEW Path prob = %e \n", det_connect,path_prob);
            #endif
    }

    // divide the above sum by the number of summands (L)
    // note: we divide because we are computing the reciprocal of the MIS_weight
    MIS_weight_reciprocal /= L;

    #ifdef DEBUG
        fprintf (stderr, "MIS weight reciprocal : NEW = %e\n", MIS_weight_reciprocal);
    #endif

    wR = lwR;
    wI = lwI;
            
    return MIS_weight_reciprocal;

}


static myReal _ConnectPath (TCircuitLayer* layer, const int l,
                        unsigned long long currentState,
                        PathVertexVector Fpath,
                        unsigned long long nextState,
                        PathVertexVector Bpath,
                        myReal& wR, myReal& wI) {
        
    myReal lwR, lwI, prob;
    layer_w(layer, l, currentState, nextState, lwR, lwI);
    #ifdef DEBUG
        /*fprintf (stderr, "\nCONNECT\n");
        fprintf (stderr, "Fsegment Pw = %.5f + i %.5f\n", Fpath_PwR, Fpath_PwI);
        fprintf (stderr, "C layer w = %.5f + i %.5f\n", lwR, lwI);
        fprintf (stderr, "Bsegment Pw = %.5f + i %.5f\n", Bpath_PwR, Bpath_PwI); */
    #endif
    // layer_weight = Fpath_Pw * (layer_weight * Bpath_Pw)
    complex_multiply(lwR, lwI, lwR, lwI, Bpath.data[l+1].PwR, Bpath.data[l+1].PwI);
    complex_multiply(lwR, lwI, lwR, lwI, Fpath.data[l].PwR, Fpath.data[l].PwI);

    // prob = Fpath_Pprob * 1.0 * Bpath_Pprob
    prob = Fpath.data[l].Pprob * Bpath.data[l+1].Pprob;
    wR = lwR;
    wI = lwI;

    return prob;
}
