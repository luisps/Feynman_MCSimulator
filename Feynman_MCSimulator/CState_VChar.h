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
    static int const NQB=4096;
    unsigned char mynumber[NQB];
    void set_data ( unsigned long long val) {
        for (int i=0 ; i< NQB ; i++) {
            mynumber[i] = val & 0x0000000000000001ull;
            val = val >> 1;
        }
    }
public:
    CState_VChar (unsigned long long val) {
        set_data(val);
        //info();
    }
    CState_VChar () { 
        reset_data();
    }
    CState_VChar& operator = (const unsigned long long val) {
        set_data(val);
        return *this;
    }
    CState_VChar& operator = (const unsigned long val) {
        set_data(val);
        return *this;
    }
    CState_VChar& operator = (const unsigned val) {
        set_data(val);
        return *this;
    }
    CState_VChar& operator = (const int val) {
        set_data(val);
        return *this;
    }
    CState_VChar& operator = (const CState_VChar &t) {
        // Check for self assignment
        if(this != &t) {
            memcpy((void *) mynumber, (void * const) t.mynumber, NQB*sizeof(unsigned char));
        }
        return *this;
    }
    CState_VChar&  operator=(std::string const n) {
        set_data (strtoull(n.c_str(), NULL, 10));
        return *this;
    }
    bool operator==(CState_VChar const o) const {
        bool eq=true;
        for (int i=0 ; i<NQB && eq ; i++)
            eq = (eq && (this->mynumber[i]==o.mynumber[i]));
        return eq;
    }
    bool operator==(unsigned long long const o) const {
        CState_VChar aux = o;
        return (*this == aux);
    }
    /* qubit handling functions for unsigned long long states */
    unsigned char qb_value (int const qb) const {
        return (mynumber[qb]);
    }
    void qb_set_value (int const qb, int const val) {
        if (val <0 || val>1) fprintf (stderr, "ERROR qb_set_value: val=%d\n", val);
        mynumber[qb] = val;
    }

    void qb_set (int const qb) {
        mynumber[qb] = 1;
    }

    void qb_reset (int const qb) {
        mynumber[qb] = 0;
    }
    void reset_data (void) {
        memset((void *) mynumber, 0, NQB*sizeof(unsigned char));
    }
    void info () const {
        fprintf (stderr,"%d\t(",NQB);
        for (int qb=0;qb<NQB;qb++) fprintf(stderr, "%d ", mynumber[qb]);
        fprintf (stderr,")\n");
    }
    // for support of gmp_printf of the GNU MP library (there is a bug...)
    /*mpz_ptr get_mpz_t() const {
        mpz_class aux = mynumber[NQB-1];
        for (int qb=NQB-2 ; qb >= 0 ; qb--) {
            aux = aux << 1;
            aux += mynumber[qb];
        }
        return aux.get_mpz_t();
    }*/
    // for support of gmp_printf of the GNU MP library
    unsigned long get_ui() const {
        unsigned long aux = mynumber[NQB-1];
        for (int qb=NQB-2 ; qb >= 0 ; qb--) {
            aux = aux << 1;
            aux += mynumber[qb];
        }
        return aux;
    }
    // for compatibility with SampleCounter
    unsigned long long get_ull() const {
        unsigned long long aux = mynumber[NQB-1];
        for (int qb=NQB-2 ; qb >= 0 ; qb--) {
            aux = aux << 1;
            aux += mynumber[qb];
        }
        return aux;
    }
};

#endif /* CState_VChar_h */
