//
//  myReal.h
//  Feynman_MCSimulator
//
//  Created by Luis Paulo Santos on 11/03/2024.
//

#ifndef myReal_h
#define myReal_h

//#define __FLOAT_MP__

#ifdef __FLOAT_MP__

#include <gmpxx.h>

#define myReal mpf_class

#else

#define __FLOAT_AS_DOUBLE__

#ifdef __FLOAT_AS_DOUBLE__
#define myReal double
#else
#define myReal float
#endif

#endif

#endif /* myReal_h */
