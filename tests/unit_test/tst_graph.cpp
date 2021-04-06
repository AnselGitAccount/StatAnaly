#include "gtest/gtest.h"
#include "../../src/graph.h"
#include "../tst_utils_graph.h"


namespace statanaly {

void isEqual_DWGraph(const DWGraph<int>& g1, const DWGraph<int>& g2) {
    isEqual_adjmat(g1.adjmat, g2.adjmat);
}

void isEqual_DUGraph(const DUGraph<int>& g1, const DUGraph<int>& g2) {
    isEqual_adjmat(g1.adjmat, g2.adjmat);
}


TEST(Graph, copy) {
    // copy operations clone the adjmat.

    std::string graphStr {"AB1,BC1,AD1,BD1,DC1"}; 
    auto Adjmat = convert2Adjmat<int>(graphStr);    

    /* DWGraph */
    DWGraph<int> graph1(Adjmat);
    DWGraph<int> graph2 = graph1;
    isEqual_DWGraph(graph1, graph2);


    /* DUGraph */
    DUGraph<int> graph3(Adjmat);
    DUGraph<int> graph4 = graph3;
    isEqual_DUGraph(graph3, graph4);
};


TEST(Graph, move) {
    // Expect to move the ownership of the heap memory.
    // No new memory is consumed.

    std::string graphStr {"AB1,BC1,AD1,BD1,DC1"}; 
    auto Adjmat = convert2Adjmat<int>(graphStr);    

    /* DWGraph */
    DWGraph<int> graph1(Adjmat);
    DWGraph<int> graph2 = std::move(graph1);
    isEqual_adjmat(graph2.adjmat, Adjmat);
    EXPECT_EQ( graph1.adjmat.data, nullptr );
    EXPECT_EQ( graph1.adjmat.len(), 0 );
    EXPECT_EQ( graph2.adjmat.len(), 4 );


    /* DUGraph */
    DUGraph<int> graph3(Adjmat);
    DUGraph<int> graph4 = std::move(graph3);
    isEqual_adjmat(graph4.adjmat, Adjmat);
    EXPECT_EQ( graph3.adjmat.data, nullptr );
    EXPECT_EQ( graph3.adjmat.len(), 0 );
    EXPECT_EQ( graph4.adjmat.len(), 4 );

}

}