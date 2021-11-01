#ifndef VNTR_CANDIDATE_H
#define VNTR_CANDIDATE_H

#include <cstdlib>
#include <cstdint>
#include <string>
#include <cmath>
#include <cfloat>
#include <vector>
#include <set>
#include <iostream>
#include <utility>

#include "wphmm.h"

struct candidate_fuzzy_motif {
    Motif_fuzzy_binary* motif;
    int32_t st, ed; // zero based, inclusive. genome position of repeat region
    int32_t n_ru;   // number of ru matched
    double viterbi_score;
    candidate_fuzzy_motif(Motif_fuzzy_binary* _m, int32_t _s, int32_t _e, int32_t _n, double _v) :
        motif(_m), st(_s), ed(_e), n_ru(_n), viterbi_score(_v) {}
    bool operator<(const candidate_fuzzy_motif & rhs) const {
        return (viterbi_score > rhs.viterbi_score);
    }
};

struct candidate_unit {
    std::string ru;
    bool inexact;
    std::vector< std::pair<int32_t, int32_t> > variable_base;
    candidate_unit(std::string _s, bool _i = 0) : ru(_s), inexact(_i) {}
    void check() {
        if (inexact && variable_base.size() != ru.size()) {
            inexact = 0;
            variable_base.clear();
            return;
        }
        if (inexact) {
            int32_t ct = 0;
            for (size_t i = 0; i < ru.size(); ++i) {
                if (variable_base[i].first > variable_base[i].second) {
                    int32_t tmp = variable_base[i].first;
                    variable_base[i].first  = variable_base[i].second;
                    variable_base[i].second = tmp;
                }
                ct += (variable_base[i].first < variable_base[i].second);
            }
            if (ct == 0) {
                inexact = 0;
                variable_base.clear();
            }
        }
    }
    bool operator<(const candidate_unit & rhs) const {
        if (inexact != rhs.inexact) {
            return (inexact < rhs.inexact);
        }
        if (ru.size() != rhs.ru.size()) {
            return (ru.size() < rhs.ru.size());
        }
        if (ru.compare(rhs.ru) == 0) {return false;}
        if (ru.size() > 1) {
            for (size_t i = 1; i < ru.size(); ++i) {
                if (rhs.ru.compare(ru.substr(i)+ru.substr(0,i)) == 0) {
                    return false;
                }
            }
        }
        return (ru.compare(rhs.ru) < 0);
    }
};

/**
 * Class for representing a VNTR.
 *
 * 2 sets of attributes for the exact and fuzzy detections of the repeat region.
 */
class VNTR_CANDIDATE
{
    public:

    //chromosome
    int32_t rid;  //rid, redundant data with Variant. todo: something about this.

    //motif
    std::vector<candidate_unit> candidate_ru;


    std::string motif;         //motif
    std::string ru;            //repeat unit on the reference
    uint32_t mlen;             //length of motif
    float motif_score;         //motif score from motif tree

    //exact repeat tract
    std::string repeat_tract;   //repeat tract
    int32_t rbeg1;              //beginning of repeat tract
    int32_t rend1;              //end of repeat tract
    float rl;                   //number of repeat units on repeat tract
    float ll;                   //number of repeat units on longest allele
    float motif_concordance;    //motif concordance from hmm
    int32_t no_exact_ru;        //number exact repeat units from hmm
    int32_t total_no_ru;        //total no of repeat units from hmm
    std::string lflank;         //left flank
    std::string rflank;         //right flank

    //fuzzy repeat tract
    std::string fuzzy_repeat_tract;   //repeat tract
    int32_t fuzzy_rbeg1;              //beginning of repeat tract
    int32_t fuzzy_rend1;              //end of repeat tract
    float fuzzy_rl;                   //number of repeat units on repeat tract
    float fuzzy_ll;                   //number of repeat units on longest allele
    float fuzzy_motif_concordance;    //motif concordance from hmm
    int32_t fuzzy_no_exact_ru;        //number exact repeat units from hmm
    int32_t fuzzy_total_no_ru;        //total no of repeat units from hmm
    std::string fuzzy_lflank;         //left flank
    std::string fuzzy_rflank;         //right flank

    //large repeat tract
    bool is_large_repeat_tract;

    /**
     * Constructor.
     */
    VNTR_CANDIDATE();

    /**
     * Clear object.
     */
    void clear();

    /**
     * Checks for equality.
     */
    bool equals(VNTR_CANDIDATE& vntr);

    /**
     * Get VNTR representation in string format.
     */
    void get_vntr_allele_string(std::string& var);

    /**
     * Get VNTR fuzzy representation in string format.
     */
    void get_fuzzy_vntr_allele_string(std::string& var);

    /**
     * Print object.
     */
    void print();
};
#endif
