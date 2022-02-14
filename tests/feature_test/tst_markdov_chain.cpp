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
#include "markov_chain.h"


namespace statanaly {

TEST(Markov_Chain, Gamblers_Ruin) {
    /* Gambler's Ruin problem */
    // The '1' in diagnal is optional.
    // std::string graphStr {"BA1,BC2,CB1,CD2"};

    std::string graphStr {"AA1,BA1,BC2,CB1,CD2,DD1"};
    Adjmat<int> myAdjmat = convert2DirectedAdjmat<int>(graphStr);

    // Normalize the rows to represent probability.
    auto normAdjmat = normalizeRow(myAdjmat);
    DWGraph myg(normAdjmat);
    Markovchain mc(myg);

    // Find transition probability from State B to all terminal states.
    auto prob_map = mc.find_transit_prob(1);
    EXPECT_EQ( 2, prob_map.size() );
    EXPECT_DOUBLE_EQ( 3./7, prob_map[0] );
    EXPECT_DOUBLE_EQ( 4./7, prob_map[3] );    
}


TEST(Markov_Chain, Google_Foobar) {
    /* Google Foobar question */

    std::string graphStr {"AB1,AF1,BA4,BD3,BE2"}; 
    Adjmat<int> myAdjmat = convert2DirectedAdjmat<int>(graphStr);

    auto normAdjmat = normalizeRow(myAdjmat);
    DWGraph myg(normAdjmat);
    Markovchain mc(myg);

    // Find transition probability from State A to all terminal states.
    auto prob_map = mc.find_transit_prob(0);
    EXPECT_EQ( 4, prob_map.size() );
    EXPECT_DOUBLE_EQ( 0,     prob_map[2] );
    EXPECT_DOUBLE_EQ( 3./14, prob_map[3] );
    EXPECT_DOUBLE_EQ( 1./7,  prob_map[4] );
    EXPECT_DOUBLE_EQ( 9./14, prob_map[5] );
}

}