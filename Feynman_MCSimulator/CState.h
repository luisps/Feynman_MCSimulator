//
//  CState.h
//  old_Feynman
//
//  Created by Luis Paulo Santos on 30/04/2024.
//

#ifndef CState_h
#define CState_h

//
// define macro below to use multiprecision integers


//#define __CSTATE_MP__
#define __CSTATE_VCHAR__

#if defined(__CSTATE_VCHAR__)

#include <gmpxx.h>
#include "CState_VChar.h"
typedef CState_VChar CState;

#elif defined(__CSTATE_MP__)
#include <gmpxx.h>
typedef mpz_class CState;

#else
typedef unsigned long long CState;

#endif

typedef unsigned long long SampleCounter;

/* State qubit handling functions for unsigned long long (all_paths) */

inline int qb_value (int qb, unsigned long long state) {
    const unsigned long long one = 1ull;
    const unsigned long long uvalue = one << qb;
    return (int) ((state & uvalue) >> qb);
}

inline unsigned long long qb_set (int qb, unsigned long long state) {
    const unsigned long long one = 1ull;
    const unsigned long long uvalue = (one << qb);
    return (state | uvalue);
}

inline unsigned long long qb_reset (int qb, unsigned long long state) {
    const unsigned long long one = 1;
    const unsigned long long uvalue = ~(one << qb);
    return state & uvalue;
}

/* State qubit handling functions for CState */

#if defined(__CSTATE_MP__) || defined(__CSTATE_VCHAR__)

inline int qb_value (int qb, CState state) {
#if defined(__CSTATE_VCHAR__)
    return state.qb_value(qb);
#elif defined(__CSTATE_MP__)
    return (int) mpz_tstbit (state.get_mpz_t(), qb);;
#else
    const CState one = 1;
    const CState uvalue = one << qb;
    return (int) ((state & uvalue) >> qb);
#endif
}

inline CState qb_set (int qb, CState state) {
#if defined(__CSTATE_VCHAR__)
    state.qb_set(qb);
    return state;
#else
    const CState one = 1;
    const CState uvalue = (one << qb);
    return (state | uvalue);
#endif
}

inline CState qb_reset (int qb, CState state) {
#if defined(__CSTATE_VCHAR__)
    state.qb_reset(qb);
    return state;
#else
    const CState one = 1;
    const CState uvalue = ~(one << qb);
    return state & uvalue;
#endif
}
#endif
#endif /* CState_h */
