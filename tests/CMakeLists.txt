cmake_minimum_required(VERSION 3.16)
set(CMAKE_CXX_STANDARD 20)

# Explicitly list the test source code and headers. The Catch header-only unit
# test framework is stored in with the test source.
set(StatAnaly_TEST_SRC
    unit_test/tst_disChiSq.cpp
    unit_test/tst_disNormal.cpp
    unit_test/tst_adjacency_matrix.cpp
    unit_test/tst_graph.cpp
    unit_test/tst_disUniform.cpp
    unit_test/tst_disBase.cpp
    unit_test/tst_disMixture.cpp
    unit_test/tst_dContainer.cpp
    unit_test/tst_disIrwinHall.cpp
    unit_test/tst_dConvolution.cpp
    unit_test/tst_dCompare.cpp
    feature_test/tst_markdov_chain.cpp
    feature_test/tst_rng_unix.cpp
)

set(StatAnaly_TEST_HEADER
)

PREPEND(StatAnaly_TEST_SRC)


## # Enable testing via CTest
## enable_testing()
## # Add our test as runnable via CTest
## add_test(NAME TestCMakeDemo  COMMAND TestCMakeDemo)


# Enable CTest
enable_testing()

# GoogleTest needs thread support
find_package(GTest REQUIRED)

# CMake supports GoogleTest
include(GoogleTest)

# Make an executable target that depends on the test source code we specified
# above.
add_executable(TestStatAnaly ${StatAnaly_TEST_SRC} ${StatAnaly_TEST_HEADER})

# Link our unit tests against the library we compiled
target_link_libraries(TestStatAnaly StatAnaly GTest::GTest GTest::Main)

# Auto-discover tests
gtest_discover_tests(TestStatAnaly)
