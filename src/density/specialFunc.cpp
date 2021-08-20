#include "density/specialFunc.h"


namespace statanaly {


double regLowerGamma(double s, double z) {
    return z <= 1. || z < s ? _regLowerGamma(s,z) : 1. - _regUpperGamma(s,z);
}


double regUpperGamma(double s, double z) {
    return z <= 1. || z < s ? 1. - _regLowerGamma(s,z) : _regUpperGamma(s,z);
}


double upperGamma(double s, double z) {
    return regUpperGamma(s,z) * std::tgamma(s);
}


double lowerGamma(double s, double z) {
    return regLowerGamma(s,z) * std::tgamma(s);
}


}   // namespace