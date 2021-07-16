# Make an explicit list of all source files in `CMakeDemo_SRC`. This is important
# because CMake is not a build system: it is a build system generator. Suppose
# you add a file foo.cpp to src/ after running cmake .. . If you set
# `CMakeDemo_SRC` with `file(GLOB ... )`, this is not passed to the makefile;
# the makefile doesn't know that foo.cpp exists and will not re-run cmake. Your
# collaborator's builds will fail and it will be unclear why. Whether you use
# file(GLOB ...) or not, you will need to re-run cmake, but with an explicit
# file list, you know beforehand why your code isn't compiling. 
set(StatAnaly_SRC
    density/disChiSq.cpp
    density/disNormal.cpp
    density/probDensFunc.cpp
    density/specialFunc.cpp
    dContainer.cpp
    dConvolution.cpp
    type_info.cpp
    adjacency.h
    graph.h
    hasher.h
    )

# Form the full path to the source files...
PREPEND(StatAnaly_SRC)
# ... and pass the variable to the parent scope.
set(StatAnaly_SRC ${StatAnaly_SRC}  PARENT_SCOPE)