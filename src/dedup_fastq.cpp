//
// Created by Jon Wright (EI) on 08/09/2017.
//

#include <iostream>
#include <vector>
#include <fstream>
#include <unordered_set>
#include <iomanip>
#include "deps/cxxopts/include/cxxopts.hpp"
#include <chrono>


class dna_qual_pair {
private:
    std::string dna;
    std::string qual;

public:
    const std::string get_dna() const {return dna;}
    const std::string get_qual() const {return qual;}

    // constructor
    dna_qual_pair(std::string &in_dna, std::string &in_qual) { dna = in_dna; qual = in_qual; };

    bool operator==(const dna_qual_pair &another) const {
        return dna == another.dna;
    }

};

namespace std {
    template<>
    struct hash<dna_qual_pair> {
        size_t operator()(const dna_qual_pair& k) const {
            return hash<string>()(k.get_dna());
        }
    };
}

std::string get_rc(const std::string &dna) {
    // reverse complement a DNA string
    std::string rc;
    for (auto c = dna.rbegin(); c != dna.rend(); ++c) {
        switch (*c) {
            case 'A':
            case 'a':
                rc.push_back('T');
                break;
            case 'C':
            case 'c':
                rc.push_back('G');
                break;
            case 'T':
            case 't':
                rc.push_back('A');
                break;
            case 'G':
            case 'g':
                rc.push_back('C');
                break;
            default:
                break;
        }
    }
    return rc;
};

bool is_canonical(const std::string &dna){
    // checks if a DNA string is canonical
    return (dna < get_rc(dna));
}

std::string make_canonical(const std::string &dna){
    // Makes a DNA string canonical if required
    if (!is_canonical(dna)){
        std::string rc = get_rc(dna);
        return rc;
    }
    return dna;
};

int main(int argc, char * argv[]) {

    std::string in_fastq;
    std::string out_prefix;
    bool suppress_rc(false);

    try {
        cxxopts::Options options("dedup_fastq", "FASTQ de-duplicator, part of the w2rap LMP pipeline.");

        options.add_options()
                ("help", "Print help")
                ("i,input", "input FASTQ", cxxopts::value<std::string>(in_fastq))
                ("o,output_prefix", "output prefix", cxxopts::value<std::string>(out_prefix))
                ("s,suppress_rc", "do not output reverse complement FASTQ", cxxopts::value<bool>(suppress_rc));

        options.parse(argc, argv);

        if (options.count("help")) {
            std::cout << options.help({""}) << std::endl;
            exit(0);
        }

        if (options.count("i") != 1 or options.count("o") != 1) {
            std::cout << "Error: please specify input FASTQ and output prefix." << std::endl
            << " Use option --help to check command line arguments." << std::endl;
            exit(1);
        }

    } catch (const cxxopts::OptionException& e) {
        std::cout << "Error parsing options: " << e.what() << std::endl;
        exit(1);
    }

    std::ifstream ifs(in_fastq);
    std::ofstream ofs(out_prefix + ".fastq");
    std::ofstream ofs_rc;

    // check input file exists
    if (!ifs.good()) {
        std::cout << "Cannot open file: " << in_fastq << std::endl;
        exit(1);
    }

    if (!suppress_rc) {
        ofs_rc.open(out_prefix + "_rc.fastq");
    }

    std::chrono::time_point<std::chrono::system_clock> start, end;

    std::cout << "Reading FASTQ..." << std::endl;

    start = std::chrono::system_clock::now();

    int64_t rec_pos;
    std::string id, seq_str, tmp, qual_str;
    std::vector<std::pair<u_int, uint64_t>> file_idx;    // a vector of pairs (length, file_offset)

    while(!ifs.eof()){
        rec_pos = ifs.tellg();
        std::getline(ifs, id);
        std::getline(ifs, seq_str);
        std::getline(ifs, tmp);
        std::getline(ifs, qual_str);

        if (id.length() > 0) {
            file_idx.emplace_back(seq_str.length(), rec_pos);
        }
    }
    ifs.close();

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds;
    elapsed_seconds = end - start;

    std::cout << "completed in " << elapsed_seconds.count() << " secs" << std::endl;

    uint64_t start_read_count = file_idx.size();
    std::cout << start_read_count << " reads to process." << std::endl;
    std::sort(file_idx.begin(), file_idx.end());

    std::cout << "De-duplicating reads..." << std::endl;

    start = std::chrono::system_clock::now();

    // iterate through the index and load sequences of the same length into a vector for processing
    std::size_t prev_length = 0;

    std::unordered_set<dna_qual_pair> seq_set;

    ifs.open(in_fastq, std::ifstream::in);   // reopen the ifstream
    uint64_t rec_count = 0;
    uint64_t rec_id = 1;

    for (auto it_idx: file_idx){
        if (it_idx.first != prev_length) {
            // dump sequences in set to fastq (if any exist)
            if (!seq_set.empty()) {
                rec_count += seq_set.size();

                // dump to fastq
                for (auto it: seq_set){
                    // R1
                    ofs << "@lmp_read_" << rec_id << " 1:N:0" << std::endl
                        << it.get_dna() << std::endl
                        << "+" << std::endl
                        << it.get_qual() << std::endl;



                    // R2, ie. R1 reverse complemented
                    if (!suppress_rc) {
                        ofs_rc << "@lmp_read_" << rec_id << " 2:N:0" << std::endl
                               << get_rc(it.get_dna()) << std::endl
                               << "+" << std::endl;
                        std::string rev_qual_str = it.get_qual();
                        std::reverse(rev_qual_str.begin(), rev_qual_str.end());
                        ofs_rc << rev_qual_str << std::endl;
                    }

                    rec_id++;
                }

                // reinitialise sequence set
                seq_set.clear();
            }
        }

        // load current seq using the filestream position
        ifs.seekg(it_idx.second);
        std::getline(ifs, id);
        std::getline(ifs, seq_str);
        std::getline(ifs, tmp);
        std::getline(ifs, qual_str);

        // make the read canonical if it isn't already
        std::string canon_seq_str = seq_str;
        std::string canon_qual_str = qual_str;
        if(!is_canonical(seq_str)) {
            canon_seq_str = make_canonical(seq_str);
            std::reverse(canon_qual_str.begin(), canon_qual_str.end());
        }

        seq_set.emplace(canon_seq_str, canon_qual_str);

        // set prev_length
        prev_length = it_idx.first;
    }

    // write the records from the last set (if any)
    if (!seq_set.empty()) {
        rec_count += seq_set.size();

        // dump to fastq
        for (auto it: seq_set) {
            // R1
            ofs << "@lmp_read_" << rec_id << " 1:N:0" << std::endl
                << it.get_dna() << std::endl
                << "+" << std::endl
                << it.get_qual() << std::endl;

            // R2, ie. R1 reverse complemented
            if (!suppress_rc) {
                ofs_rc << "@lmp_read_" << rec_id << " 2:N:0" << std::endl
                       << get_rc(it.get_dna()) << std::endl
                       << "+" << std::endl;
                std::string rev_qual_str = it.get_qual();
                std::reverse(rev_qual_str.begin(), rev_qual_str.end());
                ofs_rc << rev_qual_str << std::endl;
            }
            rec_id++;
        }
    }

    ifs.close();
    ofs.close();
    ofs_rc.close();

    end = std::chrono::system_clock::now();
    elapsed_seconds = end - start;
    std::cout << "completed in " << elapsed_seconds.count() << " secs" << std::endl;

    // calc percentage and print
    float percent = (float)rec_count / (float)start_read_count * 100.0;
    std::cout << rec_count << " (" << std::showpoint << std::fixed << std::setprecision(2) << percent << "%) reads remaining." << std::endl;
    std::cout << "DONE." << std::endl;

    return 0;

}
