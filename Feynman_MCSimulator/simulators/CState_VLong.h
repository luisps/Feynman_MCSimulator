//
//  CState_VLong.h
//  Feynman_MCSimulator
//
//  Created by Luis Paulo Santos on 14/04/2024.
//

#ifndef CState_VLong_h
#define CState_VLong_h

class CState_VLong {
private:
    static int const MAX_QBS=512;
    static int const N_ULONGS=MAX_QBS >> 6;  // each u long is 64 = 2^6 bits
    unsigned long long data[N_ULONGS];
    unsigned int nqb=0;
    void set_data ( unsigned long long val) {
        data[0] = val;
        memset((void *) &(data[1]), 0, (N_ULONGS-1)*sizeof(unsigned char));
    }
public:
    CState_VLong (unsigned int _nqb, unsigned long long val) {
        if (_nqb>MAX_QBS) {
            fprintf (stderr, "ERROR: required %d qubits. Max is %d.\n", _nqb, MAX_QBS);
            _nqb = MAX_QBS;
        }
        nqb = _nqb;
        set_data(val);
        //info();
    }
    CState_VLong (unsigned int _nqb) {
        if (_nqb>MAX_QBS) {
            fprintf (stderr, "ERROR: required %d qubits. Max is %d.\n", _nqb, MAX_QBS);
            _nqb = MAX_QBS;
        }
        nqb = _nqb;
        set_data(0ull);
        //info();
    }
    CState_VLong () { }
    CState_VLong& operator = (const CState_VLong &t) {
        // Check for self assignment
        if(this != &t) {
            nqb = t.nqb;
            memcpy((void *) data, (void * const) t.data, N_ULONGS*sizeof(unsigned long long));
        }
        return *this;
    }
    bool operator==(CState_VLong o) const {
        if (this->nqb != o.nqb) return false;
        bool eq=true;
        for (int i=0 ; i<N_ULONGS && eq ; i++)
            eq = (eq && (this->data[i]==o.data[i]));
        return eq;
    }
    /* qubit handling functions for unsigned long long states */
    unsigned char qb_value (int const qb) const {
        const int ndx = qb >> 6;
        const int bit = qb & 0x03F;
        return (data[ndx] >> bit) & 1;
    }
    void qb_set_value (int const qb, int const val) {
        if (val) {
            qb_set(qb);
        }
        else {
            qb_reset(qb);
        }
    }

    void qb_set (int const qb) {
        const int ndx = qb >> 6;
        const int bit = qb & 0x03F;
        data[ndx] |=  (1<<bit);
    }

    void qb_reset (int const qb) {
        const int ndx = qb >> 6;
        const int bit = qb & 0x03F;
        data[ndx] &= (~(1 <<bit));
    }
    void reset_data (void) {
        memset((void *) data, 0, N_ULONGS*sizeof(unsigned long long));
    }
    unsigned int get_nqbs () const {
        return nqb;
    }
    void info () const {
        fprintf (stderr,"%d\t(",nqb);
        //for (int qb=0;qb<nqb;qb++) fprintf(stderr, "%d ", data[qb]);
        fprintf (stderr,")\n");
    }
};
#endif /* CState_VLong_h */
