#include <iostream>
#include "../include/dContainer.h"

namespace statanaly {

// Human friendly text.
std::ostream& operator << (std::ostream& output, const dCtr& distr) {
    distr.print(output);
    return output;
}

}