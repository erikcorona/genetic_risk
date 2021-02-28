#include <utility>
#include <sstream>
#include <cassert>
//
// Created by dam on 2/13/21.
//

#ifndef GEN_RISK2_GWAS_HXX
#define GEN_RISK2_GWAS_HXX

//@todo implement a masking system so that no data is copied. When a user wants to get an object with only chr 6 data then the same underlying object should be used but with a hidden mask applied. will speed up the processing tremendously as no data will need to be copied
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

template<typename scalar>
auto intersect(std::vector<scalar> a, std::vector<scalar> b)
{
    std::set<scalar> s(a.begin(), a.end());
    std::set<scalar> intersect_mask;

    for(auto i : b)
        if(s.contains(i))
            intersect_mask.insert(i);

    return intersect_mask;
}

/**
 * Get all lines from a file
 * @param file the file to be read into memory
 * @return a vector of strings
 */
std::vector<std::string> get_lines(const std::string& file){
    std::ifstream infile(file);
    std::vector<std::string> lines;

    std::string line;
    while(getline( infile, line ))
        lines.push_back(line);
    assert(!lines.empty());
    return lines;
}

/**
 * The purpose of this class is to provider an interface to all GWAS results in the GWAS catalog.
 */
class GWAS{


    using column_name =             std::string ;
    using gwas_entry  = std::vector<std::string>;
    using strings     = std::vector<std::string>;

    strings header;
    std::unordered_map<column_name, std::size_t> index_of; // maps the column name to its index position
    std::vector<gwas_entry> data;
private:

    const std::string dis_col = "DISEASE/TRAIT";

    void initHeaderIndexMap()
    {
        for(std::size_t i = 0; i < header.size(); i++)
            index_of[header[i]] = i;
    }

    /**
     * Returns index of of valid positions in the GWAS object. Not all positions are valid. Some are missing or do not
     * resolve to valid numbers.
     * @return vector holding a vector of valid positions.
     */
    auto positions_mask()
    {
        auto idx = this->index_of.at("CHR_POS");

        std::vector<std::size_t> mask_pos;
        for(std::size_t i{0}; i < data.size(); i++)
        {
            auto& gwas_entry = data[i];

            unsigned long a_pos;
            try {
                a_pos = std::stoul(gwas_entry[idx]);
            }catch(const std::invalid_argument& ia)
            {
                std::cerr << gwas_entry[idx] << " is not a valid pos ";
                a_pos = -1;
            }

            if(a_pos > 0)
                mask_pos.emplace_back(i);
        }

        return mask_pos;
    }

    auto effect_size_mask()
    {
        auto es_i = index_of.at("OR or BETA");

        std::vector<std::size_t> es;
        for(std::size_t i{0}; i < data.size(); i++)
        {
            auto& gwas_entry = data[i];
            double effect_size;
            try {
                effect_size = std::stod(gwas_entry[es_i]);
            }catch(const std::invalid_argument& ia)
            {
                std::cerr << gwas_entry[es_i] << " is not a valid effect size ";
                effect_size = -1;
            }

            if(effect_size > 0)
                es.emplace_back(i);
        }


        return es;
    }

public:

    /**
     * Initialize a new GWAS object
     * @param a_header the header strings that describe the contents in each column
     * @param a_data the vector of GWAS entries that make up a set of GWAS results
     */
    GWAS (strings a_header, std::vector<gwas_entry> a_data) : header{std::move(a_header)}, data{std::move(a_data)}
    {
        //@todo initialize all column index positions (e.g., dis_col_i = disease column index
        initHeaderIndexMap();
    }

    /**
     * Instantiates a GWAS object from the file location of the GWAS catalog
     * @param file the path to the GWAS catalog TSV file
     */
    explicit GWAS(const std::string& file){

        auto lines = get_lines(file);
        header = getTokens(lines[0]);
        initHeaderIndexMap();

        lines.erase(lines.begin());
        for(auto& line : lines)
            data.push_back(getTokens(line));
    }

    /**
     * The number of GWAS entries in this object.
     * @return number of GWAS associations
     */
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

    //@todo replace with a true unit test
//    void integrityCheck()
//    {
//        for(auto& tokens : data)
//            assert(tokens.size() == header.size());
//    }


    /**
     * Get all diseases in this GWAS object.
     * @return List of all diseases in this GWAS object.
     */
    auto uniqueDiseases()
    {
        std::set<std::string> diseases;
        for(auto& gwas_entry : data)
            diseases.insert(gwas_entry[index_of.at(dis_col)]);

        return diseases;
    }

    void printSummary()
    {
        std::size_t cnt{0};
        for(auto& disease : this->uniqueDiseases())
        {
            auto dis = this->get_disease(disease);
            if(dis.size() > 9)
                cnt++;
        }

        std::cout << "associations: " << this->size() << "\tdiseases > 9 " << cnt << std::endl;
    }

    /**
     * Subsets the data in this object. It will return a version of this object with only the specified disease in the
     * data.
     * @param dis_name the disease by which to subset this object
     * @return A version of this object with only the specified disease or an empty one if the disease specified does
     * not exist within the data of this object.
     */
    GWAS get_disease(const std::string& dis_name){
        return subsetter(dis_col, dis_name);
    }

    GWAS subsetter(const std::string& col_name, const std::string& col_value){
        auto name_idx = this->index_of.at(col_name);

        std::vector<gwas_entry> new_data;
        for(auto& gwas_entry : data)
            if(gwas_entry[name_idx] == col_value)
                new_data.push_back(gwas_entry);

        auto new_header = header;
        return GWAS(new_header, new_data);
    }


    /**
     * Retrieve a subset of this GWAS object containing results only in the specified chromosome.
     * @param chr chromosome by which to subset the data
     * @return the same object but only with results present in the specified chromosome
     */
    GWAS chr_sub(const std::string& chr)
    {
        return subsetter("CHR_ID", chr);
    }

    auto positions()
    {
        auto idx = this->index_of.at("CHR_POS");
        std::vector<unsigned long> pos;

        unsigned long a_pos;
        for(auto i : positions_mask())
        {
            auto& gwas_entry = data[i];
            a_pos = std::stoul(gwas_entry[idx]);
            pos.emplace_back(a_pos);
        }

        return pos;
    }

    auto effect_size()
    {
        auto es_i = index_of.at("OR or BETA");

        std::vector<double> es;
        for(auto i : effect_size_mask())
        {
            auto& gwas_entry = data[i];
            double effect_size;
            effect_size = std::stod(gwas_entry[es_i]);
            es.emplace_back(effect_size);
        }

        return es;
    }

    /**
     * Retrieves the position and effect size of all associations in this object. This function returns all positions
     * and effect size info, even if there are multople diseases and multiple chromosomes mixed into the data of this
     * object.
     * @return position and effect size contents of this object. It only returns cases where both the effect size and
     * position are valid numbers.
     */
    auto positions_and_effect_size() {
        auto idx  = index_of.at("CHR_POS"   );
        auto es_i = index_of.at("OR or BETA");

        std::vector<std::pair<unsigned long, double>> pe;
        for (auto i: intersect(positions_mask(), effect_size_mask())) {
            auto &gwas_entry = data[i];
            auto a_pos        = std::stoul(gwas_entry[idx]);
            auto effect_size  = std::stod(gwas_entry[es_i]);
            pe.emplace_back(a_pos, effect_size);
        }

        return pe;
    }

    /**
     * Get all unique RSIDs in this GWAS object.
     * @return List of all RSIDs in this GWAS object.
     */
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
