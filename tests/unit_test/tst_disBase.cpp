#include "gtest/gtest.h"
#include "../../include/density/probDensFunc.h"
#include "../../include/density/disUniform.h"
#include "../../include/density/disNormal.h"
#include "../../include/density/disChiSq.h"

namespace statanaly {

template<class T, class F>
bool tst_isClone(const T* d1, const F* d2, double x=0) {
    // d1 and d2 have to occupy different memory.
    if (d1 == d2) return false;

    bool isSame = true;
    
    isSame &= (d1->pdf(x) == d2->pdf(x));
    isSame &= (d1->cdf(x) == d2->cdf(x));
    isSame &= (d1->mean() == d2->mean());
    isSame &= (d1->stddev() == d2->stddev());
    isSame &= (d1->variance() == d2->variance());
    isSame &= (d1->skewness() == d2->skewness());

    return isSame;
}

template<class T, class F>
bool tst_isClone(const std::unique_ptr<T>& d1, const std::unique_ptr<F>& d2, double x=0) {
    return tst_isClone(d1.get(), d2.get(), x);
};




TEST(distribution_base_class, clone_via_cloneUnique) {
    // make a clone by calling a class method that calls std::make_unique (deep-copy).

    std::unique_ptr<disStdUniform> distr = std::make_unique<disStdUniform>();
    std::unique_ptr<probDensFunc> clone = distr->cloneUnique();
    EXPECT_TRUE( tst_isClone(distr, clone) );

    std::unique_ptr<probDensFunc> distr2 = std::make_unique<disStdUniform>();
    std::unique_ptr<probDensFunc> clone2 = distr2->cloneUnique();
    EXPECT_TRUE( tst_isClone(distr2, clone2) );

    std::unique_ptr<disUniform> distr3 = std::make_unique<disUniform>(0,1);
    std::unique_ptr<probDensFunc> clone3 = distr3->cloneUnique();
    EXPECT_TRUE( tst_isClone(distr3, clone3) );

    std::unique_ptr<probDensFunc> distr4 = std::make_unique<disUniform>(0,1);
    std::unique_ptr<probDensFunc> clone4 = distr4->cloneUnique();
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
    std::unique_ptr<probDensFunc>  clone = distr->cloneUnique();
    std::unique_ptr<disStdUniform> cloneDerived (static_cast<disStdUniform*>(clone.release()));
    EXPECT_TRUE( tst_isClone(distr, cloneDerived) );

    std::unique_ptr<disUniform> distr2 = std::make_unique<disUniform>(0,1);
    std::unique_ptr<probDensFunc>  clone2 = distr2->cloneUnique();
    std::unique_ptr<disUniform> cloneDerived2 (static_cast<disUniform*>(clone2.release()));
    EXPECT_TRUE( tst_isClone(distr2, cloneDerived2) );
}


TEST(distribution_base_class, clone_via_raw_ptr) {
    // make a clone by calling a class method that allocate raw memory.

    disStdUniform* distr = new disStdUniform;
    disStdUniform* clone = distr->clone();
    probDensFunc* clone2 = distr->clone();    // implicit upcast 
    
    EXPECT_TRUE( tst_isClone(distr, clone) );
    EXPECT_TRUE( tst_isClone(distr, clone2) );

    delete distr;
    delete clone;
    delete clone2;
}


}