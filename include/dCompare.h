#ifndef STATANALY_D_COMPARE_H_
#define STATANALY_D_COMPARE_H_

#include "fl_comparison.h"
#include "density/probDensFunc.h"

namespace statanaly {    

/* Bit-wise equality. 
 */
inline bool isEqual(const probDensFunc& a, const probDensFunc& b) {
    return a.hash() == b.hash();
}


/* Equality upto a numerical tolerance 
 */
inline bool isEqual_tol(const probDensFunc& a, const probDensFunc& b, const double tol) {
    // a & b must have the same concrete type.
    return a.isEqual_tol(b,tol);
}

/* Equality upto n ULP
 */
inline bool isEqual_ulp(const probDensFunc& a, const probDensFunc& b, const unsigned ulp) {
    return a.isEqual_ulp(b,ulp);
}


}   // namespace statanaly

#endif