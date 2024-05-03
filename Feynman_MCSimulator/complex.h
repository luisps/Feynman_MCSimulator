//
//  complex.h
//  Feynman_MCSimulator
//
//  Created by Luis Paulo Santos on 04/08/2023.
//

#ifndef complex_h
#define complex_h

#include <math.h>
#include "myReal.h"

/*
 * Inline complex multiplication
 */

inline void  complex_multiply (myReal &cR, myReal &cI, const myReal aR, const myReal aI, const myReal bR, const myReal bI) {
    cR = (aR*bR - aI*bI);
    cI = (aR*bI + aI*bR);
}

inline void  complex_multiply (myReal &cR, myReal &cI, const myReal aR, const myReal aI, const float bR, const float bI) {
    cR = (aR*bR - aI*bI);
    cI = (aR*bI + aI*bR);
}

/*inline myReal  complex_multiply_real (const myReal aR, const myReal aI, const myReal bR, const myReal bI) {
    return (aR*bR - aI*bI);
}

inline myReal  complex_multiply_imag (const myReal aR, const myReal aI, const myReal bR, const myReal bI) {
    return (aR*bI + aI*bR);
}*/

inline myReal  complex_abs_square (const myReal aR, const myReal aI) {
    return (aR*aR + aI*aI);
}

inline float  complex_abs_square (const float aR, const float aI) {
    return (aR*aR + aI*aI);
}

inline myReal  complex_abs (const myReal aR, const myReal aI) {
    return (sqrt(complex_abs_square(aR, aI)));
}

#endif /* complex_h */
