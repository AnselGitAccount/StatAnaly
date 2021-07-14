# Make an explicit list of all source files in CMakeDemo_INC. This is important
# because CMake is not a build system: it is a build system generator. Suppose
# you add a file foo.cpp to src/ after running cmake .. . If you set
# CMakeDemo_INC with `file(GLOB ... )`, this is not passed to the makefile; it
# doesn't know that foo.cpp exists and will not re-run cmake. Your
# collaborator's builds will fail and it will be unclear why. Whether you use
# file(GLOB ...) or not, you will need to re-run cmake, but with an explicit
# file list, you know beforehand why your code isn't compiling. 
set(StatAnaly_INC
    density/disChiSq.h
    density/disNormal.h
    density/probDensFunc.h
    density/disUniform.h
    density/disMixture.h
    density/disIrwinHall.h
    density/disCauchy.h
    density/disExponential.h
    density/disGamma.h
    dContainer.h
    dConvolution.h
    const_values.h
    rand_num_gen.h
    statanaly.h
    markov_chain.h
    type_info.h
    double_dispatcher.h
    fl_comparison.h
    )

# Form the full path to the source files...
PREPEND(StatAnaly_INC)
# ... and pass the variable to the parent scope.
set(StatAnaly_INC ${StatAnaly_INC}  PARENT_SCOPE)