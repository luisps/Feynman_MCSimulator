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
#include <random>

#include "complex.h"
#include "layer.hpp"
#include "PreProcessorSettings.h"
#include "path.h"

static std::mutex mx_m2s, mx_s2m;
static std::condition_variable cv_m2s, cv_s2m;

// Shared predicate variables
// Master to Slave
static bool terminate = false;


#ifdef NON_ZERO_PATHS
static bool BD_paths_NOVEC_NOT (TCircuit *c,
                unsigned long long init_state, unsigned long long final_state, const unsigned long long n_samples,
                float &sumR, float &sumI, unsigned long long& n_Paths, int& non_zero_paths);
static bool BD_paths_NOVEC_T (TCircuit *c,
                unsigned long long init_state, unsigned long long final_state, unsigned long long& samples2proc, bool& taskReady, bool& resAvailable, float& T_sumR, float& T_sumI, unsigned long long& n_Paths, unsigned long long& processedSamples, int& non_zero_paths);
#else
static bool BD_paths_NOVEC_NOT (TCircuit *c,
                unsigned long long init_state, unsigned long long final_state, const unsigned long long n_samples,
                float &sumR, float &sumI, unsigned long long& n_Paths);
static bool BD_paths_NOVEC_T (TCircuit *c,
                unsigned long long init_state, unsigned long long final_state, unsigned long long& samples2proc, bool& taskReady, bool& resAvailable, float& T_sumR, float& T_sumI, unsigned long long& n_Paths, unsigned long long& processedSamples);
#endif

#ifdef NON_ZERO_PATHS
static bool BD_paths_NOVEC_kernel (TCircuit *c,
                unsigned long long init_state, unsigned long long final_state, const unsigned long long n_samples,
                float &sumR, float &sumI,
                unsigned long long& n_Paths,
               std::default_random_engine& e, std::uniform_real_distribution<float>& d,
                                   int& non_zero_paths);
#else
    static bool BD_paths_NOVEC_kernel (TCircuit *c,
                unsigned long long init_state, unsigned long long final_state, const unsigned long long n_samples,
                float &sumR, float &sumI,
                unsigned long long& n_Paths,
                std::default_random_engine& e, std::uniform_real_distribution<float>& d);
#endif

static void _ForwardPath (TCircuit *c, unsigned long long init_state,
                          std::vector<unsigned long long>& Fpath,
                            std::vector<float>& Fpath_wR, std::vector<float>& Fpath_wI,
                            std::vector<float>& Fpath_PwR, std::vector<float>& Fpath_PwI,
                          std::vector<float>& Fpath_prob, std::vector<float>& Fpath_Pprob,
                          std::default_random_engine& e, std::uniform_real_distribution<float>& d);

static void _BackwardPath (TCircuit *c, unsigned long long final_state,
                              std::vector<unsigned long long>& Bpath,
                                std::vector<float>& Bpath_wR, std::vector<float>& Bpath_wI,
                                std::vector<float>& Bpath_PwR, std::vector<float>& Bpath_PwI,
                              std::vector<float>& Bpath_prob, std::vector<float>& Bpath_Pprob,
                           std::default_random_engine& e, std::uniform_real_distribution<float>& d);

static float _ConnectPath (TCircuitLayer* layer, const int l,
                     unsigned long long currentState,
                       float Fpath_PwR, float Fpath_PwI, float Fpath_Pprob,
                       unsigned long long nextState,
                       float Bpath_PwR, float Bpath_PwI, float Bpath_Pprob,
                       float& wR, float& wI);

#ifdef CONVERGENCE_STATS
bool BD_paths (TCircuit *c, unsigned long long init_state,
                unsigned long long final_state, const unsigned long long n_samples,
               float &estimateR, float &estimateI, std::vector<T_Stats>& stats, const int n_threads) {
#else
bool BD_paths (TCircuit *c, unsigned long long init_state,
                unsigned long long final_state, const unsigned long long n_samples,
                float &estimateR, float &estimateI, const int n_threads) {
#endif

    bool ret=true;
#ifdef NON_ZERO_PATHS
    int non_zero_paths = 0;
#endif
    float sumR=0.f, sumI=0.f;
    unsigned long long n_ProcessedSamples = 0ull;
    unsigned long long n_Paths = 0ull;

    
    if (n_threads<=1) {
#ifdef NON_ZERO_PATHS
        ret = BD_paths_NOVEC_NOT (c, init_state, final_state, n_samples, sumR, sumI, n_Paths, non_zero_paths);
#else
        ret = BD_paths_NOVEC_NOT (c, init_state, final_state, n_samples, sumR, sumI, n_Paths);
#endif
        estimateR = sumR / ((float)n_Paths);
        estimateI = sumI / ((float)n_Paths);
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
        unsigned long long * nT_Paths = new unsigned long long [n_threads]; // = 0ull

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
            nT_Paths[t] = 0ull;
            idleThreads[t] = false;
#ifdef NON_ZERO_PATHS
            l_NzeroP[t] = 0;
            threads.push_back(std::thread(BD_paths_NOVEC_T, c, init_state, final_state, std::ref(samples2proc[t]), std::ref(taskReady[t]), std::ref(resAvailable[t]), std::ref(T_sumR[t]), std::ref(T_sumI[t]), std::ref(nT_Paths[t]), std::ref(processedSamples[t]), std::ref(l_NzeroP[t])));
#else
            threads.push_back(std::thread(BD_paths_NOVEC_T, c, init_state, final_state, std::ref(samples2proc[t]), std::ref(taskReady[t]), std::ref(resAvailable[t]), std::ref(T_sumR[t]), std::ref(T_sumI[t]), std::ref(nT_Paths[t]), std::ref(processedSamples[t])));
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
                        n_Paths += nT_Paths[t];
                        busy_threads--;
#ifdef CONVERGENCE_STATS
                        T_Stats stat;
                        stat.sumR = sumR;
                        stat.sumI = sumI;
                        stat.n_samples = n_ProcessedSamples;
                        stat.n_Paths = n_Paths;
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

        estimateR = sumR / n_Paths;
        estimateI = sumI / n_Paths;
    }

#ifdef NON_ZERO_PATHS
    fprintf (stderr, "Non zero paths: %d\n", non_zero_paths);
#endif
    
    return ret;
}

#ifdef NON_ZERO_PATHS
static bool BD_paths_NOVEC_NOT (TCircuit *c,
                unsigned long long init_state, unsigned long long final_state, const unsigned long long n_samples,
                              float &sumR, float &sumI,
                                unsigned long long& n_Paths,
                                int& non_zero_paths) {
#else
    static bool BD_paths_NOVEC_NOT (TCircuit *c,
                    unsigned long long init_state, unsigned long long final_state, const unsigned long long n_samples,
                                  float &sumR, float &sumI,
                                    unsigned long long& n_Paths) {
#endif
    // thread local random number generator (seeded by a local random device)
    // see https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2013/n3551.pdf
    std::random_device rdev{};
    thread_local std::default_random_engine e{rdev()};
    std::uniform_real_distribution<float>d{0.0,1.0};  // uniform distribution in[0,1[ (float)
        
#ifdef NON_ZERO_PATHS
    return BD_paths_NOVEC_kernel(c, init_state, final_state, n_samples, sumR, sumI, n_Paths, e, d, non_zero_paths);
#else
    return BD_paths_NOVEC_kernel(c, init_state, final_state, n_samples, sumR, sumI, n_Paths, e, d);
#endif
        
}

#ifdef NON_ZERO_PATHS
static bool BD_paths_NOVEC_T (TCircuit *c,
                unsigned long long init_state, unsigned long long final_state, unsigned long long& samples2proc, bool& taskReady, bool& resAvailable, float& T_sumR, float& T_sumI, unsigned long long& processedSamples,
                              unsigned long long& n_Paths,
                              int& non_zero_paths) {
#else
static bool BD_paths_NOVEC_T (TCircuit *c,
                unsigned long long init_state, unsigned long long final_state, unsigned long long& samples2proc, bool& taskReady, bool& resAvailable, float& T_sumR, float& T_sumI, unsigned long long& processedSamples,
                              unsigned long long& n_Paths) {
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
            BD_paths_NOVEC_kernel(c, init_state, final_state, n_samples, sumR, sumI, n_Paths, e, d, non_zero_paths);
#else
            BD_paths_NOVEC_kernel(c, init_state, final_state, n_samples, sumR, sumI, n_Paths, e, d);
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
static bool BD_paths_NOVEC_kernel (TCircuit *c,
                unsigned long long init_state, unsigned long long final_state, const unsigned long long n_samples,
                float &sumR, float &sumI,
                unsigned long long& n_Paths,
                std::default_random_engine& e, std::uniform_real_distribution<float>& d,
                                   int& non_zero_paths) {
#else
static bool BD_paths_NOVEC_kernel (TCircuit *c,
                                    unsigned long long init_state, unsigned long long final_state,
                                   const unsigned long long n_samples,
                                    float &sumR, float &sumI,
                                   unsigned long long& n_Paths,
                                   std::default_random_engine& e, std::uniform_real_distribution<float>& d) {
#endif
    unsigned long long s;   // sample counter
    int l;                  // layer counter
    unsigned long long l_n_Paths = 0ull;   // nbr of paths (local counter)

    float l_sumR=0.f, l_sumI=0.f;  // local summ accumulators for performance reasons
#ifdef NON_ZERO_PATHS
    int l_non_zero_paths = 0;      // local counter for performance reasons
#endif
    
    const int L = c->size->num_layers;   // number of layers in the circuit
    
    // iteratively generate samples
    for (s=0ull; s<n_samples ; s++) {
        
#ifdef DEBUG
        fprintf (stderr, "Sample %llu:\n", s);
#endif
        
        std::vector<unsigned long long> Fpath;
        std::vector<float> Fpath_wR, Fpath_wI, Fpath_PwR, Fpath_PwI, Fpath_prob, Fpath_Pprob;
        _ForwardPath (c, init_state, Fpath, Fpath_wR, Fpath_wI, Fpath_PwR, Fpath_PwI, Fpath_prob, Fpath_Pprob, e, d);
        
#ifdef DEBUG
        //print_path(Fpath);
        fprintf (stderr, "\n\t\tFORWARD\n\n ");
        fprintf (stderr, "%llu ", Fpath[0]);
        for (l=1 ; l< L ; l++) {
            fprintf (stderr, "-> %llu ", Fpath[l]);
        }
        fprintf (stderr, "||\n");
        fprintf (stderr, "\n");
        fprintf (stderr, "%.5f + i %.5f ", Fpath_PwR[0], Fpath_PwI[0]);
        for (l=1 ; l< L ; l++) {
            fprintf (stderr, "-> %.5f + i %.5f ", Fpath_PwR[l], Fpath_PwI[l]);
        }
        fprintf (stderr, "||\n");
        fprintf (stderr, "\n");
        fprintf (stderr, "%.5f ", Fpath_Pprob[0]);
        for (l=1 ; l< L ; l++) {
            fprintf (stderr, "-> %.5f ", Fpath_Pprob[l]);
        }
        fprintf (stderr, "||\n");
#endif
        std::vector<unsigned long long> Bpath;
        std::vector<float> Bpath_wR, Bpath_wI, Bpath_PwR, Bpath_PwI, Bpath_prob, Bpath_Pprob;
        _BackwardPath (c, init_state, Bpath, Bpath_wR, Bpath_wI, Bpath_PwR, Bpath_PwI, Bpath_prob, Bpath_Pprob, e, d);
        
#ifdef DEBUG
        //print_path(Bpath);
        fprintf (stderr, "\n\t\tBACKWARD\n\n ");
        fprintf (stderr, "||\t%llu ", Bpath[1]);
        for (l=2 ; l<= L ; l++) {
            fprintf (stderr, "<- %llu ", Bpath[l]);
        }
        fprintf (stderr, "\n");
        fprintf (stderr, "\n");
        fprintf (stderr, "||\t%.5f + i %.5f ", Bpath_PwR[1], Bpath_PwI[1]);
        for (l=2 ; l<= L ; l++) {
            fprintf (stderr, "<- %.5f + i %.5f ", Bpath_PwR[l], Bpath_PwI[l]);
        }
        fprintf (stderr, "\n");
        fprintf (stderr, "\n");
        fprintf (stderr, "||\t%.5f ", Bpath_Pprob[1]);
        for (l=2 ; l<= L ; l++) {
            fprintf (stderr, "<- %.5f ", Bpath_Pprob[l]);
        }
        fprintf (stderr, "\n");
        fprintf (stderr, "\n");
#endif
        
        for (l=0 ; l < L ; l++ ) {  // generate all L different paths
            float prob, wR, wI;
            TCircuitLayer *layer = &c->layers[l];
            const unsigned long long currentState = Fpath[l];
            const unsigned long long nextState = Bpath[l+1];
            
#ifdef DEBUG
            fprintf (stderr, "Connect state %llu to state %llu through layer %d:\n", currentState, nextState, l);
            fprintf (stderr, "Connecting path is:\n%llu ", Fpath[0]);
            for (int ll=1; ll <=l ; ll++) {
                fprintf (stderr, " -> %llu ", Fpath[ll]);
            }
            fprintf (stderr, " <-> ");
            fprintf (stderr, " %llu ", Bpath[l+1]);
            for (int ll=l+2; ll <=L ; ll++) {
                fprintf (stderr, " <- %llu ", Bpath[ll]);
            }
            fprintf (stderr, "\n");
#endif
            // connect deterministically fPath to bPath through layer l
            prob = _ConnectPath (layer, l,
                             currentState, Fpath_PwR[l], Fpath_PwI[l], Fpath_Pprob[l],
                             nextState, Bpath_PwR[l+1], Bpath_PwI[l+1], Bpath_Pprob[l+1],
                             wR, wI);

#ifdef DEBUG
            fprintf (stderr, "\tw = %.5f + i %.5f \t prob = %.5f\n\n", wR, wI, prob);
#endif
            if ((complex_abs_square(wR, wI) > 0.f) && prob > 0.f) {  // non_zero_path
                float path_throughputR = wR / prob;
                float path_throughputI = wI / prob;
                l_sumR += path_throughputR;
                l_sumI += path_throughputI;
#ifdef NON_ZERO_PATHS
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
                              std::vector<unsigned long long>& Fpath,
                              std::vector<float>& Fpath_wR, std::vector<float>& Fpath_wI,
                              std::vector<float>& Fpath_PwR, std::vector<float>& Fpath_PwI,
                              std::vector<float>& Fpath_prob, std::vector<float>& Fpath_Pprob,
                          std::default_random_engine& e, std::uniform_real_distribution<float>& d) {
    
    
    unsigned long long current_state, next_state;
    // The forward path 1st node is the initial state
    Fpath.push_back(init_state);
    Fpath_wR.push_back(1.0f);     // path amplituide is 1 + i0
    Fpath_wI.push_back(0.0f);     // path amplituide is 1 + i0
    Fpath_PwR.push_back(1.0f);    // product of nodes' amplitudes is 1. + i0
    Fpath_PwI.push_back(0.0f);    // product of nodes' amplitudes is 1. + i0
    Fpath_prob.push_back(1.0f);   // nodes' generation probability is 1.
    Fpath_Pprob.push_back(1.0f);   // product of nodes' generation probabilities is 1.
    
    current_state = init_state;
    
    // Evolve through all gates' layes (starting at 0) , except the last one
    // The last layer will be later on deterministically connected to the final state
    for (int l=0 ; l<c->size->num_layers-1 ; l++) {
        float l_wR, l_wI, l_prob, auxR, auxI;
        
        TCircuitLayer *layer = &(c->layers[l]);
        
        // sample this layer
        l_prob = layer_sample(layer, l, current_state, next_state, l_wR, l_wI, e, d);
        
        // add to Fpath
        Fpath.push_back(next_state);
        Fpath_wR.push_back (l_wR);
        Fpath_wI.push_back (l_wI);
        complex_multiply (auxR, auxI, Fpath_PwR.back(), Fpath_PwI.back(), l_wR, l_wI);
        Fpath_PwR.push_back (auxR);
        Fpath_PwI.push_back (auxI);
        Fpath_prob.push_back (l_prob);
        auxR = Fpath_Pprob.back() * l_prob;
        Fpath_Pprob.push_back (auxR);
        
        current_state = next_state;
    }
    return;
}
    
static void _BackwardPath (TCircuit *c, unsigned long long final_state,
                              std::vector<unsigned long long>& Bpath,
                                std::vector<float>& Bpath_wR, std::vector<float>& Bpath_wI,
                                std::vector<float>& Bpath_PwR, std::vector<float>& Bpath_PwI,
                              std::vector<float>& Bpath_prob, std::vector<float>& Bpath_Pprob,
                           std::default_random_engine& e, std::uniform_real_distribution<float>& d) {
    
    
    unsigned long long current_state, next_state;
    // The backward path 1st node is the final state: to revert later
    Bpath.push_back(final_state);
    Bpath_wR.push_back(1.0f);     // path amplituide is 1 + i0
    Bpath_wI.push_back(0.0f);     // path amplituide is 1 + i0
    Bpath_PwR.push_back(1.0f);    // product of nodes' amplitudes is 1. + i0
    Bpath_PwI.push_back(0.0f);    // product of nodes' amplitudes is 1. + i0
    Bpath_prob.push_back(1.0f);   // nodes' generation probability is 1.
    Bpath_Pprob.push_back(1.0f);   // product of nodes' generation probabilities is 1.
    
    current_state = final_state;
    
    // Evolve through all gates' layes (starting at the last one (num_lkayers-1)) , except the first one
    // The first layer will be later on deterministically connected to the initial state
    for (int l=c->size->num_layers-1 ; l>0 ; l--) {
        float l_wR, l_wI, l_prob, auxR, auxI;
        
        TCircuitLayer *layer = &(c->layers[l]);
        
        // sample this layer
        l_prob = layer_sample(layer, l, current_state, next_state, l_wR, l_wI, e, d, false);
        
        // add to Fpath
        Bpath.push_back(next_state);
        Bpath_wR.push_back (l_wR);
        Bpath_wI.push_back (l_wI);
        complex_multiply (auxR, auxI, Bpath_PwR.back(), Bpath_PwI.back(), l_wR, l_wI);
        Bpath_PwR.push_back (auxR);
        Bpath_PwI.push_back (auxI);
        Bpath_prob.push_back (l_prob);
        auxR = Bpath_Pprob.back() * l_prob;
        Bpath_Pprob.push_back (auxR);
        
        current_state = next_state;
    }
    // reverse the vectors orders: make sure the new[0] means nothing
    
    // to make sure that the element index [0] on the backward path means nothing add 0. to the vectors
    // then reverse
    Bpath.push_back(0ull);
    Bpath_wR.push_back (0.f);
    Bpath_wI.push_back (0.f);
    Bpath_PwR.push_back (0.f);
    Bpath_PwI.push_back (0.f);
    Bpath_prob.push_back (0.f);
    Bpath_Pprob.push_back (0.f);
    
    reverse(Bpath.begin(), Bpath.end());
    reverse(Bpath_wR.begin(), Bpath_wR.end());
    reverse(Bpath_wI.begin(), Bpath_wI.end());
    reverse(Bpath_PwR.begin(), Bpath_PwR.end());
    reverse(Bpath_PwI.begin(), Bpath_PwI.end());
    reverse(Bpath_prob.begin(), Bpath_prob.end());
    reverse(Bpath_Pprob.begin(), Bpath_Pprob.end());
    
    return;
}

static float _ConnectPath (TCircuitLayer* layer, const int l,
                       unsigned long long currentState,
                         float Fpath_PwR, float Fpath_PwI, float Fpath_Pprob,
                         unsigned long long nextState,
                         float Bpath_PwR, float Bpath_PwI, float Bpath_Pprob,
                       float& wR, float& wI) {
    
    float lwR, lwI, prob;
    layer_w(layer, l, currentState, nextState, lwR, lwI);
#ifdef DEBUG
    fprintf (stderr, "\nCONNECT\n");
    fprintf (stderr, "Fsegment Pw = %.5f + i %.5f\n", Fpath_PwR, Fpath_PwI);
    fprintf (stderr, "C layer w = %.5f + i %.5f\n", lwR, lwI);
    fprintf (stderr, "Bsegment Pw = %.5f + i %.5f\n", Bpath_PwR, Bpath_PwI);
#endif
    // layer_weight = Fpath_Pw * (layer_weight * Bpath_Pw)
    complex_multiply(lwR, lwI, lwR, lwI, Bpath_PwR, Bpath_PwI);
    complex_multiply(lwR, lwI, lwR, lwI, Fpath_PwR, Fpath_PwI);

    // prob = Fpath_Pprob * 1.0 * Bpath_Pprob)
    prob = Fpath_Pprob * Bpath_Pprob;
    wR = lwR;
    wI = lwI;

    return prob;
}
