#include <utility>
#include <sstream>
#include <cassert>
#include <boost/lexical_cast.hpp>

//
// Created by dam on 2/13/21.
//

#ifndef GEN_RISK2_GWAS_HXX
#define GEN_RISK2_GWAS_HXX

/**
 * Generic parser that takes a string and converts it to some other data type.
 * @tparam T The data type to which the string will be converted.
 * @param v The string that will be parsed
 * @return A new value that was parsed from the string of type T.
 */
template<typename T>
auto parser(const std::string& v) -> T {
    try {
        return boost::lexical_cast<T>(v);
    }
    catch(const boost::wrapexcept<boost::bad_lexical_cast>& bc){
        return std::numeric_limits<T>::quiet_NaN();
    }
}

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

template<typename scalar, typename collection1, typename collection2>
auto intersect(collection1 a, collection2 b)
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


private:

    strings header;
    std::vector<gwas_entry> data;

    inline void initHeaderIndexMap()
    {
        std::unordered_map<column_name, std::size_t> index_of; // maps the column name to its index position
        for(std::size_t i = 0; i < header.size(); i++)
            index_of[header[i]] = i;

        dis_i = index_of.at("DISEASE/TRAIT");
        pos_i = index_of.at("CHR_POS"      );
        es_i  = index_of.at("OR or BETA"   );
        chr_i = index_of.at("CHR_ID"       );
        rsid_i = index_of.at("SNPS"        );
    }

    /**
     * Returns index of of parseable value in the GWAS object. Not all positions are valid. Some are missing or do not
     * resolve to valid numbers.
     * @return vector holding a vector of valid positions.
     */

    /**
     * Returns index of of parseable value in the GWAS object.
     * @tparam SomeFunction A function that returns either a value that is parsed (e.g., int, double, etc.) or nan
     * @param col_name the column in the data that will be parsed
     * @param f Determines if the object is parseable. It will return nan if it is not or the right value if it can be parsed.
     * @return A vector of index positions of all parseable values of interest
     */
    template<typename SomeFunction>
    auto grab_mask(const std::size_t idx, SomeFunction f)
    {
        std::vector<std::size_t> mask_pos;
        for(std::size_t i{0}; i < data.size(); i++)
        {
            auto& gwas_entry = data[i];
            auto a_pos = f(gwas_entry[idx]);
            if(!std::isnan(a_pos))
                mask_pos.emplace_back(i);
        }

        return mask_pos;
    }

public:

    std::size_t dis_i, // disease column index
                pos_i, // chromosomal bp position column index
                es_i , // effect size column index
                chr_i, // chromosome column index
                rsid_i; // SNP ID column index

    /**
     * Initialize a new GWAS object
     * @param a_header the header strings that describe the contents in each column
     * @param a_data the vector of GWAS entries that make up a set of GWAS results
     */
    GWAS (strings a_header, std::vector<gwas_entry> a_data) : header{std::move(a_header)}, data{std::move(a_data)} // NOLINT(cppcoreguidelines-pro-type-member-init)
    {
        initHeaderIndexMap();
    }

    /**
     * Instantiates a GWAS object from the file location of the GWAS catalog
     * @param file the path to the GWAS catalog TSV file
     */
    explicit GWAS(const std::string& file){ // NOLINT(cppcoreguidelines-pro-type-member-init)

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
            diseases.insert(gwas_entry[dis_i]);

        return diseases;
    }

    void printSummary()
    {
        std::size_t cnt{0};
        for(auto& disease : this->uniqueDiseases())
        {
            auto dis = this->subsetter(dis_i, disease);
            if(dis.size() > 9)
                cnt++;
        }

        std::cout << "associations: " << this->size() << "\tdiseases > 9 " << cnt << std::endl;
    }

    GWAS subsetter(const std::size_t name_idx, const std::string& col_value){

        std::vector<gwas_entry> new_data;
        for(auto& gwas_entry : data)
            if(gwas_entry[name_idx] == col_value)
                new_data.push_back(gwas_entry);

        auto new_header = header;
        return GWAS(new_header, new_data);
    }

    /**
     * Retrieves the position and effect size of all associations in this object. This function returns all positions
     * and effect size info, even if there are multople diseases and multiple chromosomes mixed into the data of this
     * object.
     * @return position and effect size contents of this object. It only returns cases where both the effect size and
     * position are valid numbers.
     */
    auto positions_and_effect_size() {

        std::vector<std::pair<unsigned long, double>> pe;
        for (auto i : intersect<unsigned long>(grab_mask(pos_i, parser<unsigned long>), grab_mask(es_i, parser<double>))) {
            auto &gwas_entry = data[i];
            auto a_pos        = boost::lexical_cast<unsigned long>(gwas_entry[pos_i]);
            auto effect_size  = boost::lexical_cast<double       >(gwas_entry[es_i ]);
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
