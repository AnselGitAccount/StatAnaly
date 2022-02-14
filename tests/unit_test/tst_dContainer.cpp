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
#include "dContainer.h"
#include "density/disNormal.h"
#include "density/disUniform.h"

namespace statanaly {

std::unique_ptr<dCtr> constructSmartContainer() {
    std::unique_ptr<disNormal>  d1 = std::make_unique<disNormal>(2,1);
    std::unique_ptr<disNormal>  d2 = std::make_unique<disNormal>(1,2);
    std::unique_ptr<disUniform> d3 = std::make_unique<disUniform>(0,1);

    std::unique_ptr<dCtr> myCtr = std::make_unique<dCtr>();
    myCtr->insert(*d1, 1);
    myCtr->insert(*d2, 1);
    myCtr->insert(*d3, 1);

    return myCtr;
}

dCtr* constructContainer() {
    std::unique_ptr<disNormal>  d1 = std::make_unique<disNormal>(2,1);
    std::unique_ptr<disNormal>  d2 = std::make_unique<disNormal>(1,2);
    std::unique_ptr<disUniform> d3 = std::make_unique<disUniform>(0,1);
    
    dCtr* myRawCtr = new dCtr();
    myRawCtr->insert(*d1.get(), 1);
    myRawCtr->insert(*d2.get(), 1);
    myRawCtr->insert(*d3.get(), 1);

    return myRawCtr;
}


TEST( dContainer_Tests, insert ) {
    std::unique_ptr<dCtr> myCtr = constructSmartContainer();
    dCtr* myRawCtr = constructContainer();    

    EXPECT_TRUE( myCtr->hash() == myRawCtr->hash() );

    int w = 1;
    disNormal* dRaw = new disNormal(2,1);
    myRawCtr->insert( *dRaw, w );
    myRawCtr->insert( disNormal(2,1), w );
    myRawCtr->insert( *(std::move(dRaw)), w );

    std::unique_ptr<disNormal> dSmart = std::make_unique<disNormal>(2,1);
    myCtr->insert( *dSmart, w );
    myCtr->insert( *(std::make_unique<disNormal>(2,1)), w );
    myCtr->insert( *(std::move(dSmart)), w );

    EXPECT_TRUE( myCtr->hash() == myRawCtr->hash() );

    delete myRawCtr;
    delete dRaw;
}

TEST( dContainer_Tests, find ) {
    std::unique_ptr<dCtr> myCtr = constructSmartContainer();

    disNormal* dRaw = new disNormal(2,1);
    auto it1 = myCtr->find( *dRaw );
    auto it2 = myCtr->find( disNormal(2,1) );
    auto it3 = myCtr->find( *(std::make_unique<disNormal>(2,1)) );
    auto it4 = myCtr->find( *(std::move(dRaw)) );

    EXPECT_TRUE( it1 != myCtr->end() );
    EXPECT_TRUE( it2 != myCtr->end() );
    EXPECT_TRUE( it3 != myCtr->end() );
    EXPECT_TRUE( it4 != myCtr->end() );

    auto it5 = myCtr->find( disNormal(2,8) );
    EXPECT_TRUE( it5 == myCtr->end() );

    delete dRaw;
}

TEST( dContainer_Tests, hash ) {
    // Test whether two containers are different.

    std::unique_ptr<dCtr> d1 = constructSmartContainer();
    std::unique_ptr<dCtr> d2 = constructSmartContainer();
    std::unique_ptr<disNormal> d3 = std::make_unique<disNormal>(2,1);

    // Must have the same distribution parameters.
    EXPECT_TRUE( d1->hash() == d2->hash() );
    d2->insert(*d3, 1);
    EXPECT_FALSE( d1->hash() == d2->hash() );

    // Weights must be the same as well.
    d1->insert(*d3, 2);
    EXPECT_FALSE( d1->hash() == d2->hash() );
}


TEST( dContainer_Tests, clone_via_cloneUnique ) {
    std::unique_ptr<dCtr> ctr = constructSmartContainer();
    std::unique_ptr<dCtr> clone = ctr->cloneUnique();
    
    EXPECT_TRUE( ctr->hash() == clone->hash() );
}

TEST( dContainer_Tests, clone_via_make_unique ) {
    std::unique_ptr<dCtr> ctr = constructSmartContainer();
    std::unique_ptr<dCtr> clone = std::make_unique<dCtr>(static_cast<dCtr const&>(*(ctr)));
    
    EXPECT_TRUE( ctr->hash() == clone->hash() );
}

TEST( dContainer_Tests, clone_via_raw_ptr ) {
    dCtr* ctr = constructContainer();
    dCtr* clone = ctr->clone();
    
    EXPECT_TRUE( ctr->hash() == clone->hash()  );

    delete ctr;
    delete clone;
}

TEST( dContainer_Tests, move ) {
    dCtr* distr = constructContainer();
    dCtr* truth = constructContainer();
    dCtr* other (std::move(distr));  // move construction
    
    EXPECT_TRUE( truth->hash() == other->hash() );

    dCtr* other2 = nullptr;
    other2 = std::move(other);   // move assignment

    EXPECT_TRUE( truth->hash() == other2->hash() );
}

TEST( dContainer_Tests, rescale ) {
    // Rescale weight distribution after named distribution insertion and deletion.
    
    std::unique_ptr<disNormal>  d1 = std::make_unique<disNormal>(2,1);
    std::unique_ptr<disNormal>  d2 = std::make_unique<disNormal>(1,2);
    std::unique_ptr<disUniform> d3 = std::make_unique<disUniform>(0,1);

    std::unique_ptr<dCtr> ctr = std::make_unique<dCtr>();
    ctr->insert(*d1, 6.2);
    ctr->insert(*d2, 3.4);
    ctr->insert(*d3, 2.0);

    const auto& ingreds = ctr->get();
    for (const auto& [d,ws] : ingreds) {
        if (d->hash() == d1->hash()) {
            EXPECT_DOUBLE_EQ( 6.2/11.6, ws.second );
        } else if (d->hash() == d2->hash()) {
            EXPECT_DOUBLE_EQ( 3.4/11.6, ws.second );
        } else if (d->hash() == d3->hash()) {
            EXPECT_DOUBLE_EQ( 2.0/11.6, ws.second );
        } else {
            EXPECT_TRUE( false ) << "no match\n";
        }
    }
}

}