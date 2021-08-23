#include <iostream>
#include "dContainer.h"

namespace statanaly {

std::ostream& operator << (std::ostream& output, const dCtr& distr) {
    distr.print(output);
    return output;
}

}