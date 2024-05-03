//
//  PathVertex.h
//  Feynman_MCSimulator
//
//  Created by Luis Paulo Santos on 25/03/2024.
//

#ifndef PathVertex_h
#define PathVertex_h

#include "myReal.h"
#include "CState.h"

typedef  struct pV{
    CState state;
    myReal wR, wI, PwR, PwI, prob, Pprob;
    pV() {};
} PathVertex, *PathVertexPtr;

class PathVertexVector {
    int nElements,head, tail;
public:
    PathVertexPtr data;
    
    PathVertexVector (int _len) {
        nElements = _len;
        data = new PathVertex[nElements];
        head = 0;
        tail = nElements -1;
    }
    ~PathVertexVector (void) {}
    /*    delete[] data;
        nElements = 0;
        head = 1;
        tail = -1;
    }*/
    int append (void) {  // append from the start
        if (head >= nElements) return -1;
        head++;
        return head-1;
    }
    int prepend (void) {  // prepend from the end
        if (tail < 0) return -1;
        tail--;
        return tail+1;
    }
    void clear (void) {        
        head = 0;
        tail = nElements -1;
    }
};

#endif /* PathVertex_h */
