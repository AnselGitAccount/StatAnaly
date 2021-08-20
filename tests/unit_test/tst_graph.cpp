#include "gtest/gtest.h"
#include "graph.h"
#include "../tst_utils_graph.h"


namespace statanaly {

TEST(Graph, copy) {
    // copy operations clone the adjmat.

    std::string graphStr {"AB1,BC1,AD1,BD1,DC1"}; 
    auto Adjmat = convert2Adjmat<int>(graphStr);    

    /* DWGraph */
    DWGraph<int> graph1(Adjmat);
    DWGraph<int> graph2 = graph1;
    EXPECT_EQ(graph1, graph2);


    /* DUGraph */
    DUGraph<int> graph3(Adjmat);
    DUGraph<int> graph4 = graph3;
    EXPECT_EQ(graph3, graph4);
};


TEST(Graph, move) {
    // Expect to move the ownership of the heap memory.
    // No new memory is consumed.

    std::string graphStr {"AB1,BC1,AD1,BD1,DC1"}; 
    auto Adjmat = convert2Adjmat<int>(graphStr);    

    /* DWGraph */
    DWGraph<int> graph1(Adjmat);
    DWGraph<int> graph2 = std::move(graph1);
    EXPECT_EQ( Adjmat, graph2.adjmat() );
    EXPECT_EQ( 0, graph1.adjmat().len() );
    EXPECT_EQ( 4, graph2.adjmat().len() );


    /* DUGraph */
    DUGraph<int> graph3(Adjmat);
    DUGraph<int> graph4 = std::move(graph3);
    EXPECT_EQ( Adjmat, graph4.adjmat() );
    EXPECT_EQ( 0, graph3.adjmat().len() );
    EXPECT_EQ( 4, graph4.adjmat().len() );

}

}