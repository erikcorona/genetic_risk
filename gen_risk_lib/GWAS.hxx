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

/**
 * Returns the intersection of two interables.
 * @tparam scalar The data type of the iterables
 * @tparam collection1 The data type of the first set
 * @tparam collection2 The data type of the second set
 * @param a The first set
 * @param b The second set
 * @return The intersection betewen a and b
 */
template<typename scalar, typename collection1, typename collection2>
auto intersect(collection1 a, collection2 b) -> std::set<scalar>
{
    std::set<scalar> s(a.begin(), a.end());
    std::set<scalar> intersect_mask;

    for(auto i : b)
        if(s.contains(i)) // log complexity on each lookup
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
 * This class represents a flat file that can be read in. The file should have the same number of columns in each line
 * and should have a header.
 */
class FlatFile{

    using col_nm  = std::string             ;
    using strings = std::vector<std::string>;
    using a_row   = std::vector<std::string>;

private:

    strings          header;
    std::vector<a_row> data;

public:

    std::unordered_map<col_nm, std::size_t> index_of; // maps the column name to its index position


    a_row & ith_row(std::size_t i) { return this->data[i]; }

    void print_header(){
        for(std::string& s : header)
            std::cout << s << std::endl;
    }

    auto num_rows() const{ return data.size(); }

    auto cell(std::size_t row, std::size_t col) -> std::string&
    {
        return data[row][col];
    }

    inline void initHeaderIndexMap()
    {
        for(std::size_t i = 0; i < header.size(); i++)
            index_of[header[i]] = i;
    }

    FlatFile(const FlatFile & f)
    {
        this->header   = f.header;
        this->data     = f.data  ;
        this->index_of = f.index_of;
    }

    FlatFile(FlatFile &&f)  noexcept {
        this->header   = std::move(f.header);
        this->data     = std::move(f.data  );
        this->index_of = std::move(f.index_of);
    }

    FlatFile& operator=(FlatFile&& other) noexcept {
        this->header   = std::move(other.header);
        this->data     = std::move(other.data);
        this->index_of = std::move(other.index_of);
        return *this;
    }

    FlatFile(strings a_header, std::vector<a_row> a_data){ // NOLINT(cppcoreguidelines-pro-type-member-init)
        header = std::move(a_header);
        data   = std::move(a_data  );
        initHeaderIndexMap();
    }

    explicit FlatFile(const std::string& file){ // NOLINT(cppcoreguidelines-pro-type-member-init)

        // Read in header
        auto lines = get_lines(file);
        header = getTokens(lines[0]);
        initHeaderIndexMap();

        // read in the rest of the data
        lines.erase(lines.begin());              // delete header
        for(auto& line : lines)                  // read in every line in flat file
            data.push_back(getTokens(line));
    }

    /**
     * Get all unique values in a column
     * @param col_i the index from which to get all unique values
     * @return a set of all unique calues in the specified column
     */
    auto unique_col(const std::size_t col_i) -> std::set<std::string> {
        std::set<std::string> col_v; // unique col values
        for(auto& gwas_entry : data)
            col_v.insert(gwas_entry[col_i]);

        return col_v;
    }

    /**
     * Creates a smaller version of this object based on some conditions
     * @param name_idx the name of the column that will be matched for a value
     * @param col_value the value that column at name_idx must have for subsetting
     * @return a smaller version of this object where col at name_idx matches a value
     */
    FlatFile subsetter2(const std::size_t name_idx, const std::string& col_value){
        return FlatFile(header, trim(name_idx, col_value));
    }


    std::vector<a_row> trim(const std::size_t name_idx, const std::string& col_value){

        std::vector<a_row> new_data;
        for(auto& gwas_entry : data)
            if(gwas_entry[name_idx] == col_value)
                new_data.push_back(gwas_entry);

        return new_data;
    }
};

/**
 * The purpose of this class is to provider an interface to all GWAS results in the GWAS catalog.
 */
class GWAS{

    using gwas_entry = std::vector<std::string>;

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
        for(std::size_t i{0}; i < file.num_rows(); i++)
            if(!std::isnan(f(file.cell(i,idx)))) // data[i] is a gwas_entry, data[i][idx] is a value in a gwas entry
                mask_pos.emplace_back(i);

        return mask_pos;
    }

    [[nodiscard]] gwas_entry& ith_gwas(std::size_t i) {
        return file.ith_row(i);
    }

public:

    FlatFile file;


//    explicit GWAS(FlatFile& f) : file(std::move(f)){}

    explicit GWAS(FlatFile&& f) : file(std::move(f)){}


    /**
     * Instantiates a GWAS object from the file location of the GWAS catalog
     * @param file the path to the GWAS catalog TSV file
     */
    explicit GWAS(const std::string& file_nm) : file(FlatFile(file_nm)){}

    /**
     * The number of GWAS entries in this object.
     * @return number of GWAS associations
     */
    [[nodiscard]] auto size() const {return file.num_rows();} // Number of associations


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
    [[nodiscard]] auto uniqueDiseases()
    {
        return file.unique_col(file.index_of.at("DISEASE/TRAIT"));
    }

    void printSummary()
    {
        std::size_t cnt{0};
        for(auto& disease : this->uniqueDiseases())
        {
            auto dis = this->subsetter("DISEASE/TRAIT", disease);
            if(dis.size() > 9)
                cnt++;
        }

        std::cout << "associations: " << this->size() << "\tdiseases > 9 " << cnt << std::endl;

        file.print_header();
    }


    GWAS subsetter(const std::string& col_nm, const std::string& col_value) {

        return GWAS(file.subsetter2(file.index_of.at(col_nm), col_value));
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
        for (auto i : intersect<unsigned long>(grab_mask(file.index_of.at("CHR_POS"), parser<unsigned long>), grab_mask(file.index_of.at("OR or BETA"), parser<double>))) {
            auto &gwas_entry = this->ith_gwas(i);
            auto a_pos        = boost::lexical_cast<unsigned long>(gwas_entry[file.index_of.at("CHR_POS")]);
            auto effect_size  = boost::lexical_cast<double       >(gwas_entry[file.index_of.at("OR or BETA")]);
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
//    template<int a, typename partner>
//    auto positions_and_X() {
//
//        if (a == 1)
//            std::cout << "3" << std::endl;
//
//        std::vector<std::pair<unsigned long, partner>> pe;
//        for (auto i : intersect<unsigned long>(grab_mask(file->pos_i, parser<unsigned long>), grab_mask(file->es_i, parser<double>))) { // for every entry with a valid position and effect size
//            auto &gwas_entry = this->ith_gwas(i);//file->data[i];
//            auto a_pos        = boost::lexical_cast<unsigned long>(gwas_entry[file->pos_i]);
//            auto effect_size  = boost::lexical_cast<double       >(gwas_entry[file->es_i ]);
//            pe.emplace_back(a_pos, effect_size);
//        }
//
//        return pe;
//    }

    /**
     * Get all unique RSIDs in this GWAS object.
     * @return List of all RSIDs in this GWAS object.
     */
    [[nodiscard]] auto uniqueRSIDs()
    {
        return file.unique_col(file.index_of.at("SNPS"));
    }

};


#endif //GEN_RISK2_GWAS_HXX
