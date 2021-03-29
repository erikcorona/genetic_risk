#include <iostream>
#include <random>
#include <fstream>
#include <set>

#include "GWAS.hxx"

//enum Allele {A,G,T,C};
//std::unordered_map<int, Allele> get_allele = {{0, A}, {1, G}, {2, T}, {3, C}};


//@todo separate out file reading functionality
//@todo add google unit tests
int main() {

    auto gwas = GWAS("gwas_catalog_v1.0-associations_e100_r2021-02-25.tsv");

    gwas.printSummary();
//    gwas.integrityCheck();

    auto t2d = gwas.subsetter(gwas.file->index_of.at("DISEASE/TRAIT"),"Type 2 diabetes");

    t2d.printSummary();
    std::cout << "Unique RSIDs for t2d: " << t2d.uniqueRSIDs().size() << std::endl;

    auto t2d6 = t2d.subsetter(t2d.file->index_of.at("CHR_ID"),"6");

    t2d6.printSummary();
    auto pos = t2d6.positions_and_effect_size();

    for(auto p : pos)
        std::cout << p.first << ",";
    std::cout << std::endl;

    for(auto p : pos)
        std::cout << p.second << ",";
    std::cout << std::endl;

    std::ofstream myfile;
    myfile.open ("adis.csv");

    //@todo investigate different odds ratios at the exact same position and also see if they have the same risk allele
    //@todo see if some genome regions are have a higher prior to being associated with a disease, more than chance allows. there may be other MHC-type regions
    for(auto& dis_nm : gwas.uniqueDiseases())
    {
        auto dis = gwas.subsetter(gwas.file->index_of.at("DISEASE/TRAIT"), dis_nm);

        std::vector<std::string> chrs = {"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"};
        for(auto& chr : chrs) {
            auto dischr = dis.subsetter(dis.file->index_of.at("CHR_ID"),chr);
            auto pos_ES = dischr.positions_and_effect_size();
            if(!pos_ES.empty())
                std::cout << dis_nm << ":" << chr << " size is " << pos_ES.size() << std::endl;

            if (pos_ES.size() > 1546) {
                for (auto &pes: pos_ES)
                    if (pes.second < 1)
                        myfile << pes.first << "," << pes.second << "," << std::endl;
                return 1;
            }
        }

    }

    return 0;
}
