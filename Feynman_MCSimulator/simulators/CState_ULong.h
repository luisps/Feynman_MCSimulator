//
//  CState_ULong.h
//  Feynman_MCSimulator
//
//  Created by Luis Paulo Santos on 06/04/2024.
//

#ifndef CState_ULong_h
#define CState_ULong_h

class CState_ULong {
    int const MAX_QBS=64;
    unsigned long long data;
    int nqb;
public:
    CState_ULong (unsigned long long _val) {nqb = MAX_QBS; data = _val;}
    CState_ULong () {data = 0ull; nqb = MAX_QBS;}
    CState_ULong (int _nqbs) {
        if (_nqbs > MAX_QBS) {
            fprintf (stderr, "ERROR: %d qubits required for state when max is %d\n!", _nqbs, MAX_QBS);
        }
        nqb = MAX_QBS; data = 0ull;}  // for compatibility with CState_VChar
    CState_ULong (int _nqbs, unsigned long long _val) {
        if (_nqbs > MAX_QBS) {
            fprintf (stderr, "ERROR: %d qubits required for state when max is %d\n!", _nqbs, MAX_QBS);
        }
        nqb = MAX_QBS;
        data = _val;
    } // for compatibility with CState_VChar
    CState_ULong&  operator=(unsigned long long n) {
        data = n;
        return *this;
    }
    CState_ULong& operator = (const CState_ULong &t) {
        // Check for self assignment
        if(this != &t) {
            nqb = t.nqb;
            data = t.data;
        }
        return *this;
    }

    bool operator==(CState_ULong o) const {
        return (this->data == o.data);
    }
    /* qubit handling functions for unsigned long long states */
    unsigned char qb_value (int qb) const {
        //const unsigned long long uvalue = 1ull <<qb;
        //return ((state & uvalue) >> qb);
        return (data >> qb) & 1;
    }
    void qb_set_value (int qb, int val) {
        if (val) {
            qb_set(qb);
        }
        else {
            qb_reset(qb);
        }
    }

    void qb_set (int qb) {
        //const unsigned long long uvalue = (1ull <<qb);
        //return (state | uvalue);
        data |=  (1<<qb);
    }

    void qb_reset (int qb) {
        //const unsigned long long uvalue = ~(1ull <<qb);
        data &= (~(1 <<qb));
    }
    unsigned int get_nqbs () const {return nqb;}
    void reset_data () { data = 0ull; }
};

#endif /* CState_ULong_h */
