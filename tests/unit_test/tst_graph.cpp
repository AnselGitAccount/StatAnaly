/*
   Copyright 2022, Ansel Blumers

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/
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