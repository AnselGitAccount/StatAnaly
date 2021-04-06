#include "gtest/gtest.h"
#include "../../src/adjacency.h"
#include "../tst_utils_graph.h"


namespace statanaly {


template<class T>
void isEqual_adjmat(const Adjmat<T>& a, const Adjmat<T>& b) {
    EXPECT_EQ(a.len(), b.len());

    const int bsize = sizeof(T) * a.len();
    EXPECT_EQ( 0, memcmp(a[0], b[0], sizeof(bsize)) );
    EXPECT_EQ( 0, memcmp(a[1], b[1], sizeof(bsize)) );
    EXPECT_EQ( 0, memcmp(a[2], b[2], sizeof(bsize)) );
    EXPECT_EQ( 0, memcmp(a[3], b[3], sizeof(bsize)) );
}



TEST(Adjacency_Matrix, construction) {
    // Check for construction

    std::string graphStr {"AB1,BC1,AD1,BD1,DC1"}; 

    auto Adjmat = convert2Adjmat<int>(graphStr);

    int r0[4] {0,1,0,1};
    int r1[4] {1,0,1,1};
    int r2[4] {0,1,0,1};
    int r3[4] {1,1,1,0};

    EXPECT_EQ( 0, memcmp(r0, Adjmat[0], sizeof(r0)) );
    EXPECT_EQ( 0, memcmp(r1, Adjmat[1], sizeof(r1)) );
    EXPECT_EQ( 0, memcmp(r2, Adjmat[2], sizeof(r2)) );
    EXPECT_EQ( 0, memcmp(r3, Adjmat[3], sizeof(r3)) );
    EXPECT_EQ( 4, Adjmat.len() );
}


TEST(Adjacency_Matrix, clone) {
    // Expect the clone to be identical to the original.

    std::string graphStr {"AB1,BC1,AD1,BD1,DC1"}; 

    auto Adjmat = convert2Adjmat<int>(graphStr);    
    decltype(Adjmat) Adjmat_clone;
    Adjmat_clone.clone(Adjmat);

    isEqual_adjmat(Adjmat, Adjmat_clone);
};


TEST(Adjacency_Matrix, move) {
    // Expect to move the ownership of the heap memory.
    // No new memory is consumed.

    std::string graphStr {"AB1,BC1,AD1,BD1,DC1"}; 

    auto Adjmat = convert2Adjmat<int>(graphStr);
    auto Adjmat_move = std::move(Adjmat);

    isEqual_adjmat(Adjmat_move, convert2Adjmat<int>(graphStr));
    EXPECT_EQ( Adjmat_move.len(), 4 );

    EXPECT_EQ( Adjmat.data, nullptr );
    EXPECT_EQ( Adjmat.len(), 0 );
};

}