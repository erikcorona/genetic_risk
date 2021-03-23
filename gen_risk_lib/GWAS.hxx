#include <memory>
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

class FlatFile{

    using col_nm  = std::string             ;
    using strings = std::vector<std::string>;
    using a_row   = std::vector<std::string>;

private:

    strings header;
    std::vector<a_row> data;

public:

    std::size_t dis_i, // disease column index
    pos_i, // chromosomal bp position column index
    es_i , // effect size column index
    chr_i, // chromosome column index
    rsid_i; // SNP ID column index


    a_row& ith_row(std::size_t i)
    {
        return this->data[i];
    }

    void print_header(){
        for(std::string& s : header)
            std::cout << s << std::endl;
    }

    auto num_rows(){ return data.size(); }

    auto cell(std::size_t row, std::size_t col) -> std::string&
    {
        return data[row][col];
    }

    inline void initHeaderIndexMap()
    {
        std::unordered_map<col_nm, std::size_t> index_of; // maps the column name to its index position
        for(std::size_t i = 0; i < header.size(); i++)
            index_of[header[i]] = i;

        dis_i = index_of.at("DISEASE/TRAIT");
        pos_i = index_of.at("CHR_POS"      );
        es_i  = index_of.at("OR or BETA"   );
        chr_i = index_of.at("CHR_ID"       );
        rsid_i = index_of.at("SNPS"        );
    }

    FlatFile(strings a_header, std::vector<a_row> a_data){ // NOLINT(cppcoreguidelines-pro-type-member-init)
        header = std::move(a_header);
        data   = std::move(a_data  );
        initHeaderIndexMap();
    }

//    template<typename AString>
    explicit FlatFile(const std::string& file){ // NOLINT(cppcoreguidelines-pro-type-member-init)

        auto lines = get_lines(file);
        header = getTokens(lines[0]);
        initHeaderIndexMap();

        lines.erase(lines.begin());
        for(auto& line : lines)
            data.push_back(getTokens(line));
    }

    auto unique_col(const std::size_t col_i) -> std::set<std::string>
    {
        std::set<std::string> col_v; // unique col values
        for(auto& gwas_entry : data)
            col_v.insert(gwas_entry[col_i]);

        return col_v;
    }

    std::unique_ptr<FlatFile> subsetter(const std::size_t name_idx, const std::string& col_value){

        std::vector<a_row> new_data;
        for(auto& gwas_entry : this->data)
            if(gwas_entry[name_idx] == col_value)
                new_data.push_back(gwas_entry);

        auto new_header = this->header;

        return std::make_unique<FlatFile>(new_header, new_data);
    }
};

/**
 * The purpose of this class is to provider an interface to all GWAS results in the GWAS catalog.
 */
class GWAS{

    using gwas_entry = std::vector<std::string>;
    using strings    = std::vector<std::string>;

private:

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
        for(std::size_t i{0}; i < file->num_rows(); i++)
            if(!std::isnan(f(file->cell(i,idx)))) // data[i] is a gwas_entry, data[i][idx] is a value in a gwas entry
                mask_pos.emplace_back(i);

        return mask_pos;
    }

    gwas_entry& ith_gwas(std::size_t i){
        return file->ith_row(i);
    }

public:

    std::unique_ptr<FlatFile> file;


    /**
     * Initialize a new GWAS object
     * @param a_header the header strings that describe the contents in each column
     * @param a_data the vector of GWAS entries that make up a set of GWAS results
     */

    GWAS (std::unique_ptr<FlatFile>& f){
        file = std::move(f);
    }

    GWAS (const strings& a_header, const std::vector<gwas_entry>& a_data)
    {
        file = std::make_unique<FlatFile>(a_header, a_data);
    }

    /**
     * Instantiates a GWAS object from the file location of the GWAS catalog
     * @param file the path to the GWAS catalog TSV file
     */
    explicit GWAS(const std::string& file_nm){ // NOLINT(cppcoreguidelines-pro-type-member-init)
        file = std::make_unique<FlatFile>(file_nm);
    }

    /**
     * The number of GWAS entries in this object.
     * @return number of GWAS associations
     */
    [[nodiscard]] auto size() const {return file->num_rows();} // Number of associations


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
        return file->unique_col(file->dis_i);
    }

    void printSummary()
    {
        std::size_t cnt{0};
        for(auto& disease : this->uniqueDiseases())
        {
            auto dis = this->subsetter(file->dis_i, disease);
            if(dis.size() > 9)
                cnt++;
        }

        std::cout << "associations: " << this->size() << "\tdiseases > 9 " << cnt << std::endl;

        file->print_header();
    }

    GWAS subsetter(const std::size_t name_idx, const std::string& col_value){

        std::unique_ptr<FlatFile> new_f = file->subsetter(name_idx, col_value);

        return GWAS(new_f);
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
        for (auto i : intersect<unsigned long>(grab_mask(file->pos_i, parser<unsigned long>), grab_mask(file->es_i, parser<double>))) {
            auto &gwas_entry = this->ith_gwas(i);
            auto a_pos        = boost::lexical_cast<unsigned long>(gwas_entry[file->pos_i]);
            auto effect_size  = boost::lexical_cast<double       >(gwas_entry[file->es_i ]);
            pe.emplace_back(a_pos, effect_size);
        }

        return pe;
    }

    /**
     * Retrieves the position and effect size of all associations in this object. This function returns all positions
     * and effect size info, even if there are multople diseases and multiple chromosomes mixed into the data of this
     * object.
     * @return position and effect size contents of this object. It only returns cases where both the effect size and
     * position are valid numbers.
     */
    template<int a, typename partner>
    auto positions_and_X() {

        if (a == 1)
            std::cout << "3" << std::endl;

        std::vector<std::pair<unsigned long, partner>> pe;
        for (auto i : intersect<unsigned long>(grab_mask(file->pos_i, parser<unsigned long>), grab_mask(file->es_i, parser<double>))) { // for every entry with a valid position and effect size
            auto &gwas_entry = this->ith_gwas(i);//file->data[i];
            auto a_pos        = boost::lexical_cast<unsigned long>(gwas_entry[file->pos_i]);
            auto effect_size  = boost::lexical_cast<double       >(gwas_entry[file->es_i ]);
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
        return file->unique_col(file->rsid_i);
    }

};


#endif //GEN_RISK2_GWAS_HXX
