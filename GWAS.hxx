#include <utility>

//
// Created by dam on 2/13/21.
//

#ifndef GEN_RISK2_GWAS_HXX
#define GEN_RISK2_GWAS_HXX

/**
 * Used to convert a string to a set of tab delimeted tokens.
 * @param line a string reference
 * @return list of tokens
 */
std::vector<std::string> getTokens(std::string& line){
    std::vector<std::string> tokens;

    std::istringstream iss(line);
    std::string token;
    while (std::getline(iss, token, '\t'))
        tokens.push_back(token);
    return tokens;
}

std::vector<std::string> get_lines(const std::string& file){
    std::ifstream infile(file);
    std::vector<std::string> lines;

    std::string line;
    while(getline( infile, line ))
        lines.push_back(line);
    assert(!lines.empty());
    return lines;
}

class GWAS{

    using column_name = std::string;
    using gwas_entry  = std::vector<std::string>;
    using strings = std::vector<std::string>;

    strings header;
    std::unordered_map<column_name, std::size_t> index_of; // maps the column name to its index position
    std::vector<gwas_entry> data;

public:

    GWAS (strings a_header, std::unordered_map<column_name, std::size_t> a_index_of, std::vector<gwas_entry> a_data)
    : header{std::move(a_header)}, index_of{std::move(a_index_of)}, data{std::move(a_data)}
    {

    }
    explicit GWAS(const std::string& file){

        auto lines = get_lines(file);
        header = getTokens(lines[0]);
        for(std::size_t i = 0; i < header.size(); i++)
            index_of[header[i]] = i;

        lines.erase(lines.begin());
        for(auto& line : lines)
            data.push_back(getTokens(line));
    }

    auto size() {return data.size();} // Number of associations

    void printHeader(){
        for(std::string& s : header)
            std::cout << s << std::endl;
    }

    void print(std::size_t i)
    {
        for(std::string& s : data[i])
            std::cout << s << std::endl;
    }

    void integrityCheck()
    {
        for(auto& tokens : data)
            assert(tokens.size() == header.size());
    }

    auto diseases(){

        std::unordered_map<std::string, std::size_t> disease_counts;
        for(auto& d : data)
            disease_counts[d[index_of.at("DISEASE/TRAIT")]]++;
        return disease_counts;
    }

    void disease_counts()
    {
        std::cout << "num associations in each disease" << std::endl;
        auto dis_cnts = this->diseases();
        for(auto& pair : dis_cnts)
            std::cout << pair.first << "\t" << pair.second << std::endl;
    }

    //@todo return a new gwas class but made up of only 1 disease
    void printSummary()
    {
        auto dis_counts = this->diseases();
        std::size_t cnt{0};

        for(auto& pair : dis_counts)
            if(pair.second > 9)
                cnt++;

        std::cout << "associations: " << this->size() <<
                  "\tdiseases > 9 " << cnt << std::endl;
    }

    /**
     * Subsets the data in this object. It will return a version of this object with only the specified disease in the
     * data.
     * @param dis_name the disease by which to subset this object
     * @return A version of this object with only the specified disease or an empty one if the disease specified does
     * not exist within the data of this object.
     */
    GWAS get_disease(const std::string& dis_name){
        auto name_idx = this->index_of.at("DISEASE/TRAIT");

        std::vector<gwas_entry> new_data;
        for(auto& gwas_entry : data)
            if(gwas_entry[name_idx] == dis_name)
                new_data.push_back(gwas_entry);

        auto new_header = header;
        auto new_index_of = index_of;
        return GWAS(new_header, new_index_of, new_data);
    }

    auto uniqueRSIDs()
    {
        auto rsid_i = index_of.at("SNPS");
        std::set<std::string> rsids;

        for(auto& gwas_entry : data) {
            std::string& rsid = gwas_entry[rsid_i];
            if (rsid.starts_with("rs") && rsid.find(' ') == std::string::npos && rsid.find('\t') == std::string::npos && rsid.find(';') == std::string::npos)
                rsids.insert(rsid);
        }

        return rsids;
    }
};


#endif //GEN_RISK2_GWAS_HXX
