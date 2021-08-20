#include <iostream>
#include "density/probDensFunc.h"

namespace statanaly {

// Human friendly text.
std::ostream& operator << (std::ostream& output, const probDensFunc& distr) {
    distr.print(output);
    return output;
}

}