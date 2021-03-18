//
// Created by dam on 3/18/21.
//

#include <gtest/gtest.h>
#include "GWAS.hxx"

TEST(AAA,BBB){ // 12/2/2020 -> 737761
    int a = 3;

    EXPECT_EQ(a,3);
}


int main(int argc, char* argv[])
{
    testing::InitGoogleTest(&argc, argv);
    auto v = RUN_ALL_TESTS();
    return v;
}


