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
#include "adjacency.h"
#include "../tst_utils_graph.h"


namespace statanaly {


TEST(Adjacency_Matrix, construction) {
    // Check for construction

    std::string graphStr {"AB1,BC1,AD1,BD1,DC1"}; 

    auto Adjmat = convert2Adjmat<int>(graphStr);

    std::vector<std::vector<int>> raw(4);
    raw[0] = std::vector<int>{0,1,0,1};
    raw[1] = std::vector<int>{1,0,1,1};
    raw[2] = std::vector<int>{0,1,0,1};
    raw[3] = std::vector<int>{1,1,1,0};

    EXPECT_EQ( raw, Adjmat.data() );
    EXPECT_EQ( 4,   Adjmat.len() );
}


TEST(Adjacency_Matrix, move) {
    // Expect to move the ownership of the heap memory.

    std::string graphStr {"AB1,BC1,AD1,BD1,DC1"}; 

    auto Adjmat = convert2Adjmat<int>(graphStr);
    auto Adjmat_move = std::move(Adjmat);

    EXPECT_EQ( convert2Adjmat<int>(graphStr), Adjmat_move );
    EXPECT_EQ( 4, Adjmat_move.len() );
    EXPECT_EQ( 0, Adjmat.len() );
};

}