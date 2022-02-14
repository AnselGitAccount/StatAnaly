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
#include "density/probDistr.h"
#include "density/disUniform.h"
#include "density/disNormal.h"
#include "density/disChiSq.h"

namespace statanaly {

template<class T, class F>
bool tst_isClone(const T* d1, const F* d2, double x=0) {
    // d1 and d2 MUST occupy different memory space.
    if (d1 == d2) return false;
    
    return d1->hash() == d2->hash();
}

template<class T, class F>
bool tst_isClone(const std::unique_ptr<T>& d1, const std::unique_ptr<F>& d2, double x=0) {
    return tst_isClone(d1.get(), d2.get(), x);
};




TEST(distribution_base_class, clone_via_cloneUnique) {
    // make a clone by calling a class method that calls std::make_unique (deep-copy).

    std::unique_ptr<disStdUniform> distr = std::make_unique<disStdUniform>();
    std::unique_ptr<probDistr> clone = distr->cloneUnique();
    EXPECT_TRUE( tst_isClone(distr, clone) );

    std::unique_ptr<probDistr> distr2 = std::make_unique<disStdUniform>();
    std::unique_ptr<probDistr> clone2 = distr2->cloneUnique();
    EXPECT_TRUE( tst_isClone(distr2, clone2) );

    std::unique_ptr<disUniform> distr3 = std::make_unique<disUniform>(0,1);
    std::unique_ptr<probDistr> clone3 = distr3->cloneUnique();
    EXPECT_TRUE( tst_isClone(distr3, clone3) );

    std::unique_ptr<probDistr> distr4 = std::make_unique<disUniform>(0,1);
    std::unique_ptr<probDistr> clone4 = distr4->cloneUnique();
    EXPECT_TRUE( tst_isClone(distr4, clone4) );

}

TEST(distribution_base_class, clone_via_make_unique) {
    // make a clone by calling std::make_unique (deep-copy).
    // Here std::make_unique uses copy constructor.

    std::unique_ptr<disStdUniform> distr = std::make_unique<disStdUniform>();
    std::unique_ptr<disStdUniform> clone = std::make_unique<disStdUniform>(*(distr.get()));
    EXPECT_TRUE( tst_isClone(distr, clone) );

    std::unique_ptr<disUniform> distr3 = std::make_unique<disUniform>(0,1);
    std::unique_ptr<disUniform> clone3 = std::make_unique<disUniform>(*(distr3.get()));
    EXPECT_TRUE( tst_isClone(distr3, clone3) );


    // const provides extra layer of security.
    std::unique_ptr<disStdUniform> distr2 = std::make_unique<disStdUniform>();
    std::unique_ptr<disStdUniform> clone2 = std::make_unique<disStdUniform>(static_cast<disStdUniform const&>(*(distr2.get())));
    EXPECT_TRUE( tst_isClone(distr2, clone2) );

    std::unique_ptr<disUniform> distr4 = std::make_unique<disUniform>(0,1);
    std::unique_ptr<disUniform> clone4 = std::make_unique<disUniform>(static_cast<disUniform const&>(*(distr4.get())));
    EXPECT_TRUE( tst_isClone(distr4, clone4) );
}


TEST(distribution_base_class, cast_to_derived_via_unique_ptr) {
    // After cloning, cast a Base to a Derived with unique_ptr.

    std::unique_ptr<disStdUniform> distr = std::make_unique<disStdUniform>();
    std::unique_ptr<probDistr>  clone = distr->cloneUnique();
    std::unique_ptr<disStdUniform> cloneDerived (static_cast<disStdUniform*>(clone.release()));
    EXPECT_TRUE( tst_isClone(distr, cloneDerived) );

    std::unique_ptr<disUniform> distr2 = std::make_unique<disUniform>(0,1);
    std::unique_ptr<probDistr>  clone2 = distr2->cloneUnique();
    std::unique_ptr<disUniform> cloneDerived2 (static_cast<disUniform*>(clone2.release()));
    EXPECT_TRUE( tst_isClone(distr2, cloneDerived2) );
}


TEST(distribution_base_class, clone_via_raw_ptr) {
    // make a clone by calling a class method that allocate raw memory.

    disStdUniform* distr = new disStdUniform;
    disStdUniform* clone = distr->clone();
    probDistr* clone2 = distr->clone();    // implicit upcast 
    
    EXPECT_TRUE( tst_isClone(distr, clone) );
    EXPECT_TRUE( tst_isClone(distr, clone2) );

    delete distr;
    delete clone;
    delete clone2;
}


}