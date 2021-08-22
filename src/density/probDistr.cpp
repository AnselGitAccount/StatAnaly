#include <iostream>
#include "density/probDistr.h"

namespace statanaly {

// Human friendly text.
std::ostream& operator << (std::ostream& output, const probDistr& distr) {
    distr.print(output);
    return output;
}

}