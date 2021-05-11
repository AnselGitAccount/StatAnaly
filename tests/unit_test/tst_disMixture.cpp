#include "gtest/gtest.h"
#include "../../include/density/disMixture.h"
#include "../../include/density/disNormal.h"
#include "../../include/density/disUniform.h"


namespace statanaly {


std::unique_ptr<disMixture> constructSmartMixture() {
    std::unique_ptr<disNormal>  d1 = std::make_unique<disNormal>(2,1);
    std::unique_ptr<disNormal>  d2 = std::make_unique<disNormal>(1,2);
    std::unique_ptr<disUniform> d3 = std::make_unique<disUniform>(0,1);

    std::unique_ptr<disMixture> myMixture = std::make_unique<disMixture>();
    myMixture->insert(d1, 1);
    myMixture->insert(d2, 1);
    myMixture->insert(d3, 1);

    return myMixture;
}

disMixture* constructMixture() {
    std::unique_ptr<disNormal>  d1 = std::make_unique<disNormal>(2,1);
    std::unique_ptr<disNormal>  d2 = std::make_unique<disNormal>(1,2);
    std::unique_ptr<disUniform> d3 = std::make_unique<disUniform>(0,1);
    
    disMixture* myRawMixture = new disMixture();
    myRawMixture->insert(d1.get(), 1);
    myRawMixture->insert(d2.get(), 1);
    myRawMixture->insert(d3.get(), 1);

    return myRawMixture;
}


TEST(Mixture_Distribution_Tests, instantiation) {
    std::unique_ptr<disMixture> myMixture = constructSmartMixture();
    disMixture* myRawMixture = constructMixture();    

    // TODO: Compare myMixture and myRawMixture.
    EXPECT_TRUE(1);

    delete myRawMixture;
}

TEST(Mixture_Distribution_Tests, clone_via_cloneUnique) {
    std::unique_ptr<disMixture> distr = constructSmartMixture();
    std::unique_ptr<probDensFunc> clone = distr->cloneUnique();
    
    // TODO: Compare distr and clone.
    EXPECT_TRUE(1);
}

TEST(Mixture_Distribution_Tests, clone_via_make_unique) {
    std::unique_ptr<disMixture> distr = constructSmartMixture();
    std::unique_ptr<probDensFunc> clone = std::make_unique<disMixture>(static_cast<disMixture const&>(*(distr)));
    
    // TODO: Compare distr and clone.
    EXPECT_TRUE(1);
}

TEST(Mixture_Distribution_Tests, clone_via_raw_ptr) {
    disMixture* distr = constructMixture();
    disMixture* clone = distr->clone();
    probDensFunc* clone2 = distr->clone();    // implicit upcast 
    
    // TODO: Compare distr and clone.
    EXPECT_TRUE( 1 );
    EXPECT_TRUE( 1 );

    delete distr;
    delete clone;
    delete clone2;
}

TEST(Mixture_Distribution_Tests, move) {
    disMixture* distr = constructMixture();
    disMixture* truth = constructMixture();

    // move construction
    disMixture* other (std::move(distr));
    
    // TODO: Compare other and truth.
    EXPECT_TRUE( 1 );

    // move assignment
    disMixture* other2 = nullptr;
    other2 = std::move(other);

    // TODO: Compare other and truth.
    EXPECT_TRUE( 1 );
}


}