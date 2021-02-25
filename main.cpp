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


#include <iostream>
#include <fstream>




int main() {

    auto gwas = GWAS("gwas_catalog_v1.0-associations_e100_r2021-01-29.tsv");

    gwas.printSummary();
    gwas.printHeader();
    gwas.print(1);
//    gwas.integrityCheck();

    auto t2d = gwas.get_disease("Type 2 diabetes");

    t2d.printSummary();
    std::cout << "Unique RSIDs for t2d: " << t2d.uniqueRSIDs().size() << std::endl;

    auto t2d6 = t2d.getChr("6");

    t2d6.printSummary();
    auto pos = t2d6.positions();

    for(auto p : pos)
        std::cout << p.first << ",";
    std::cout << std::endl;

    for(auto p : pos)
        std::cout << p.second << ",";
    std::cout << std::endl;

    std::ofstream myfile;
    myfile.open ("adis.csv");

    for(auto& dis_nm : gwas.uniqueDiseases())
    {
        auto dis  = gwas.get_disease(dis_nm);

        std::vector<std::string> chrs = {"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"};
        for(auto& chr : chrs) {
            auto dischr = dis.getChr(chr);
            auto positions = dischr.positions();
            std::cout << dis_nm << ":" << chr << " size is " << positions.size() << std::endl;

            if (positions.size() > 70) {
                for (auto &pos : positions)
                    if (pos.second < 1)
                        myfile << pos.first << "," << pos.second << "," << std::endl;
                return 1;
            }
        }

    }

    return 0;
}
