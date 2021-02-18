#include <iostream>
#include <random>
#include <fstream>
#include <sstream>
#include <cassert>
#include <set>
#include <unordered_map>

#include "GWAS.hxx"

//enum Allele {A,G,T,C};
//std::unordered_map<int, Allele> get_allele = {{0, A}, {1, G}, {2, T}, {3, C}};

//using rsNum =   long;
//using risk  = double;





int main() {

    auto gwas = GWAS("gwas_catalog_v1.0-associations_e100_r2021-01-29.tsv");

    gwas.printSummary();
    gwas.printHeader();
    gwas.print(1);
    gwas.integrityCheck();
    gwas.disease_counts();
    auto t2d = gwas.get_disease("Type 2 diabetes");
    t2d.disease_counts();
    t2d.printSummary();
    std::cout << "Unique RSIDs for t2d: " << t2d.uniqueRSIDs().size() << std::endl;

    auto t2d6 = t2d.getChr("6");

    t2d6.printSummary();
    //t2d6.positions();


    return 0;
}
