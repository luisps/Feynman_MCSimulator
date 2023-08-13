//
//  complex.h
//  Feynman_MCSimulator
//
//  Created by Luis Paulo Santos on 04/08/2023.
//

#ifndef complex_h
#define complex_h

#include <math.h>

/*
 * Inline complex multiplication
 */

inline void  complex_multiply (float &cR, float &cI, const float aR, const float aI, const float bR, const float bI) {
    cR = (aR*bR - aI*bI);
    cI = (aR*bI + aI*bR);
}

/*inline float  complex_multiply_real (const float aR, const float aI, const float bR, const float bI) {
    return (aR*bR - aI*bI);
}

inline float  complex_multiply_imag (const float aR, const float aI, const float bR, const float bI) {
    return (aR*bI + aI*bR);
}*/

inline float  complex_abs_square (const float aR, const float aI) {
    return (aR*aR + aI*aI);
}

inline float  complex_abs (const float aR, const float aI) {
    return (sqrtf(complex_abs_square(aR, aI)));
}

#endif /* complex_h */
