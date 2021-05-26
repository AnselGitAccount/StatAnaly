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
    myMixture->insert(*d1, 1);
    myMixture->insert(*d2, 1);
    myMixture->insert(*d3, 1);

    return myMixture;
}

disMixture* constructMixture() {
    std::unique_ptr<disNormal>  d1 = std::make_unique<disNormal>(2,1);
    std::unique_ptr<disNormal>  d2 = std::make_unique<disNormal>(1,2);
    std::unique_ptr<disUniform> d3 = std::make_unique<disUniform>(0,1);
    
    disMixture* myRawMixture = new disMixture();
    myRawMixture->insert(*d1, 1);
    myRawMixture->insert(*d2, 1);
    myRawMixture->insert(*d3, 1);

    return myRawMixture;
}


TEST( Mixture_Distribution_Tests, insert ) {
    std::unique_ptr<disMixture> myMixture = constructSmartMixture();
    disMixture* myRawMixture = constructMixture();    

    EXPECT_TRUE( myMixture->hash() == myRawMixture->hash() );

    int w = 1;
    disNormal* dRaw = new disNormal(2,1);
    myRawMixture->insert( *dRaw, w );
    myRawMixture->insert( disNormal(2,1), w );
    myRawMixture->insert( *(std::move(dRaw)), w );

    std::unique_ptr<disNormal> dSmart = std::make_unique<disNormal>(2,1);
    myMixture->insert( *dSmart, w );
    myMixture->insert( *(std::make_unique<disNormal>(2,1)), w );
    myMixture->insert( *(std::move(dSmart)), w );

    EXPECT_TRUE( myMixture->hash() == myRawMixture->hash() );


    delete myRawMixture;
    delete dRaw;
}

TEST( Mixture_Distribution_Tests, find ) {
    std::unique_ptr<disMixture> myMixture = constructSmartMixture();

    disNormal* dRaw = new disNormal(2,1);
    auto it1 = myMixture->find( *dRaw );
    auto it2 = myMixture->find( disNormal(2,1) );
    auto it3 = myMixture->find( *(std::make_unique<disNormal>(2,1)) );
    auto it4 = myMixture->find( *(std::move(dRaw)) );

    EXPECT_TRUE( it1 != myMixture->end() );
    EXPECT_TRUE( it2 != myMixture->end() );
    EXPECT_TRUE( it3 != myMixture->end() );
    EXPECT_TRUE( it4 != myMixture->end() );

    auto it5 = myMixture->find( disNormal(2,8) );
    EXPECT_TRUE( it5 == myMixture->end() );

    delete dRaw;
}

TEST( Mixture_Distribution_Tests, hash ) {
    // Test whether two Mixture distributions are different.

    std::unique_ptr<disMixture> d1 = constructSmartMixture();
    std::unique_ptr<disMixture> d2 = constructSmartMixture();
    EXPECT_TRUE( d1->hash() == d2->hash() );

    std::unique_ptr<disNormal> d3 = std::make_unique<disNormal>(2,1);
    d2->insert(*d3, 1);
    EXPECT_TRUE( d1->hash() != d2->hash() );
}


TEST( Mixture_Distribution_Tests, clone_via_cloneUnique ) {
    std::unique_ptr<disMixture> distr = constructSmartMixture();
    std::unique_ptr<probDensFunc> clone = distr->cloneUnique();
    
    EXPECT_TRUE( distr->hash() == clone->hash() );
}

TEST( Mixture_Distribution_Tests, clone_via_make_unique ) {
    std::unique_ptr<disMixture> distr = constructSmartMixture();
    std::unique_ptr<probDensFunc> clone = std::make_unique<disMixture>(static_cast<disMixture const&>(*(distr)));
    
    EXPECT_TRUE( distr->hash() == clone->hash() );
}

TEST( Mixture_Distribution_Tests, clone_via_raw_ptr ) {
    disMixture* distr = constructMixture();
    disMixture* clone = distr->clone();
    probDensFunc* clone2 = distr->clone();    // implicit upcast 
    
    EXPECT_TRUE( distr->hash() == clone->hash()  );
    EXPECT_TRUE( distr->hash() == clone2->hash() );

    delete distr;
    delete clone;
    delete clone2;
}

TEST( Mixture_Distribution_Tests, move ) {
    disMixture* distr = constructMixture();
    disMixture* truth = constructMixture();
    disMixture* other (std::move(distr));  // move construction
    
    EXPECT_TRUE( truth->hash() == other->hash() );

    disMixture* other2 = nullptr;
    other2 = std::move(other);   // move assignment

    EXPECT_TRUE( truth->hash() == other2->hash() );
}

TEST( Mixture_Distribution_Tests, rescale ) {
    // Rescale weight distribution after named distribution insertion and deletion.
    
    std::unique_ptr<disNormal>  d1 = std::make_unique<disNormal>(2,1);
    std::unique_ptr<disNormal>  d2 = std::make_unique<disNormal>(1,2);
    std::unique_ptr<disUniform> d3 = std::make_unique<disUniform>(0,1);

    std::unique_ptr<disMixture> myMixture = std::make_unique<disMixture>();
    myMixture->insert(*d1, 6.2);
    myMixture->insert(*d2, 3.4);
    myMixture->insert(*d3, 2.0);

    const auto& ingreds = myMixture->get();
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

/* Statistical tests */

std::unique_ptr<disMixture>  construct_M1() {
    std::unique_ptr<disUniform> d1 = std::make_unique<disUniform>(0,1);
    std::unique_ptr<disUniform> d2 = std::make_unique<disUniform>(0.5,3.5);
    std::unique_ptr<disMixture> distr = std::make_unique<disMixture>();
    distr->insert(*d1, 1);
    distr->insert(*d2, 2);
    return distr;
}

TEST( Mixture_Distribution_Tests, pdf_M1 ) {
    auto distr = construct_M1();
    
    double exp_pdf = 1./3+2./3/3;
    EXPECT_DOUBLE_EQ( exp_pdf, distr->pdf(0.7) );
}

TEST( Mixture_Distribution_Tests, cdf_M1 ) {
    auto distr = construct_M1();
    
    double exp_cdf = 0.7/3+2./3*0.2/3;
    EXPECT_DOUBLE_EQ( exp_cdf, distr->cdf(0.7) );
}

TEST( Mixture_Distribution_Tests, mean_M1 ) {
    auto distr = construct_M1();

    double exp_mean = 1./6 + 4./3;
    EXPECT_DOUBLE_EQ( exp_mean, distr->mean() );
}

TEST( Mixture_Distribution_Tests, variance_M1 ) {
    auto distr = construct_M1();

    double exp_variance = 1./36+1./2 - (1./6+4./3)*(1./6+4./3) + (1./12+8./3);
    EXPECT_DOUBLE_EQ( exp_variance, distr->variance() );
}

TEST( Mixture_Distribution_Tests, skewness_M1 ) {
    auto distr = construct_M1();

    double exp_skewness = 1./24 + 3 + 1./24 + 16./3;
    exp_skewness += - 3*distr->mean()*distr->variance()
                    - distr->mean()*distr->mean()*distr->mean();
    exp_skewness /= distr->stddev()*distr->stddev()*distr->stddev();

    // The numerical discrepancy is larger than the tolerance of EXPECT_DOUBLE_EQ.
    EXPECT_FLOAT_EQ( exp_skewness, distr->skewness() );
}


std::unique_ptr<disMixture>  construct_M2() {
    std::unique_ptr<disNormal>  d1 = std::make_unique<disNormal>(1,2);
    std::unique_ptr<disNormal>  d2 = std::make_unique<disNormal>(2,1);
    std::unique_ptr<disMixture> distr = std::make_unique<disMixture>();
    distr->insert(*d1, 1);
    distr->insert(*d2, 2);
    return distr;
}

TEST( Mixture_Distribution_Tests, pdf_M2 ) {
    auto distr = construct_M2();

    double pdf_d1 = exp(-log(sqrt(2)) - SACV_LOG_SQRT_2PI - 0.5*(0.7*0.7/2));
    double pdf_d2 = exp(-log(sqrt(1)) - SACV_LOG_SQRT_2PI - 0.5*(0.3*0.3));
    double exp_pdf = pdf_d1/3 + pdf_d2*2/3;

    EXPECT_FLOAT_EQ( exp_pdf, distr->pdf(1.7) );
}

TEST( Mixture_Distribution_Tests, mean_M2 ) {
    auto distr = construct_M2();

    double exp_mean = 1./3+4./3;
    EXPECT_DOUBLE_EQ( exp_mean, distr->mean() );
}

TEST( Mixture_Distribution_Tests, variance_M2 ) {
    auto distr = construct_M2();

    double exp_variance = 2./3 + 2./3 - 5./3*5/3 + 1./3 + 8./3;
    EXPECT_DOUBLE_EQ( exp_variance, distr->variance() );
}

TEST( Mixture_Distribution_Tests, skewness_M2 ) {
    auto distr = construct_M2();

    double exp_skewness = 2 + 4 + 17./3;
    exp_skewness += - 3*distr->mean()*distr->variance()
                    - distr->mean()*distr->mean()*distr->mean();
    exp_skewness /= distr->stddev()*distr->stddev()*distr->stddev();

    // The numerical discrepancy is larger than the tolerance of EXPECT_DOUBLE_EQ.
    EXPECT_FLOAT_EQ( exp_skewness, distr->skewness() );
}


}