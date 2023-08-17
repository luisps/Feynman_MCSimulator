//
//  main.cpp
//  GEMM
//
//  Created by Luis Paulo Santos on 13/08/2023.
//

#include <iostream>
#include <cstdlib>
#include <vector>
#include <thread>
#include <mutex>

using namespace std;

// https://www.geeksforgeeks.org/measure-execution-time-function-cpp/
#include <chrono>
using namespace std::chrono;

static void GEMM (void);
static void GEMM_T (const int r0, const int rf, const int i, int &, std::mutex &);


const int M_SIZE=2048;
const int N_THREADS=48;

float A[M_SIZE][M_SIZE],B[M_SIZE][M_SIZE],C[M_SIZE][M_SIZE];

int main(int argc, const char * argv[]) {
    // protect access to sum
    std::mutex sumMutex;
    int mySum = 0;


    // initialize
    for (int i=0 ; i < M_SIZE ; i++) {
        for (int j=0 ; j < M_SIZE ; j++) {
            A[i][j] = ((float) rand())/RAND_MAX;
            B[i][j] = ((float) rand())/RAND_MAX;
            C[i][j] = 0.f;
        }
    }
    
    // https://www.geeksforgeeks.org/measure-execution-time-function-cpp/
    auto start = high_resolution_clock::now();

    std::vector<std::thread> threads;
    
    const int N_ROWS = M_SIZE / N_THREADS;
    for (int i = 0; i < N_THREADS; i++) {
        threads.push_back(std::thread(GEMM_T, i*N_ROWS, (i+1)*N_ROWS, i, std::ref(mySum), std::ref(sumMutex)));
    }
     
    for (auto &th : threads) {
        th.join();
    }
    //GEMM ();
    //std::thread myThread (GEMM);
    
   // myThread.join();
    
    auto stop = high_resolution_clock::now();
    
    cout << mySum << std::endl ;
    
    auto duration = duration_cast<microseconds>(stop - start);
    
    // To get the value of duration use the count()
    // member function on the duration object
    cout << duration.count() << " us!" << endl;


    return 0;
}

static void GEMM_T (const int r0, const int rf, const int r, int& sum_ref, std::mutex& sumMutex) {
    int i,j,k;
    int suml = 0;
    
    suml = r+2*r;
    
    // multiply
    for (i=r0 ; i < rf ; i++) {
        for (k=0 ; k < M_SIZE ; k++) {
            const float aik = A[i][k];
            for (j=0 ; j < M_SIZE ; j++) {
                C[i][j] += aik * B[k][j];
            }
        }
    }
    {
        // the access to this function is mutually exclusive
        std::lock_guard<std::mutex> guard(sumMutex);
        sum_ref += suml;
    }
}

static void GEMM (void) {
    int i,j,k;
    
    // multiply
    for (i=0 ; i < M_SIZE ; i++) {
        for (k=0 ; k < M_SIZE ; k++) {
            const float aik = A[i][k];
            for (j=0 ; j < M_SIZE ; j++) {
                C[i][j] += aik * B[k][j];
            }
        }
    }
}


