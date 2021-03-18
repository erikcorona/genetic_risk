//
// Created by dam on 3/18/21.
//

#include <gtest/gtest.h>
#include "GWAS.hxx"

TEST(AAA,BBB){ // 12/2/2020 -> 737761
    auto gwas = GWAS("gwas_catalog_v1.0-associations_e100_r2021-02-25.tsv");
    gwas.integrityCheck();
}


int main(int argc, char* argv[])
{
    testing::InitGoogleTest(&argc, argv);
    auto v = RUN_ALL_TESTS();
    return v;
}


