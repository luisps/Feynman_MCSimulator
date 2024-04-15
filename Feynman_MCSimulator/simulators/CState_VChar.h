//
//  CState_VChar.h
//  Feynman_MCSimulator
//
//  Created by Luis Paulo Santos on 07/04/2024.
//

#ifndef CState_VChar_h
#define CState_VChar_h

class CState_VChar {
private:
    static int const MAX_QBS=512;
    unsigned char data[MAX_QBS];
    unsigned int nqb=0;
    void set_data ( unsigned long long val) {
        for (int i=0 ; i< nqb ; i++) {
            data[i] = val & 0x0000000000000001ull;
            val = val >> 1;
        }
    }
public:
    CState_VChar (unsigned int _nqb, unsigned long long val) {
        if (_nqb>MAX_QBS) {
            fprintf (stderr, "ERROR: required %d qubits. Max is %d.\n", _nqb, MAX_QBS);
            _nqb = MAX_QBS;
        }
        nqb = _nqb;
        set_data(val);
        //info();
    }
    CState_VChar (unsigned int _nqb) {
        if (_nqb>MAX_QBS) {
            fprintf (stderr, "ERROR: required %d qubits. Max is %d.\n", _nqb, MAX_QBS);
            _nqb = MAX_QBS;
        }
        nqb = _nqb;
        set_data(0ull);
        //info();
    }
    CState_VChar () { }
    CState_VChar& operator = (const CState_VChar &t) {
        // Check for self assignment
        if(this != &t) {
            nqb = t.nqb;
            memcpy((void *) data, (void * const) t.data, nqb*sizeof(unsigned char));
        }
        return *this;
    }
    bool operator==(CState_VChar o) const {
        if (this->nqb != o.nqb) return false;
        bool eq=true;
        for (int i=0 ; i<nqb && eq ; i++)
            eq = (eq && (this->data[i]==o.data[i]));
        return eq;
    }
    /* qubit handling functions for unsigned long long states */
    unsigned char qb_value (int const qb) const {
        return (data[qb]);
    }
    void qb_set_value (int const qb, int const val) {
        if (val <0 || val>1) fprintf (stderr, "ERROR qb_set_value: val=%d\n", val);
        data[qb] = val;
    }

    void qb_set (int const qb) {
        data[qb] = 1;
    }

    void qb_reset (int const qb) {
        data[qb] = 0;
    }
    void reset_data (void) {
        memset((void *) data, 0, MAX_QBS*sizeof(unsigned char));
    }
    unsigned int get_nqbs () const {
        return nqb;
    }
    void info () const {
        fprintf (stderr,"%d\t(",nqb);
        for (int qb=0;qb<nqb;qb++) fprintf(stderr, "%d ", data[qb]);
        fprintf (stderr,")\n");
    }
};

#endif /* CState_VChar_h */
