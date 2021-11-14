#ifndef VNTR_CANDIDATE_H
#define VNTR_CANDIDATE_H

#include <cstdlib>
#include <cstdint>
#include <string>
#include <cmath>
#include <cfloat>
#include <vector>
#include <set>
#include <list>
#include <iostream>
#include <utility>

#include "htslib/vcf.h"
#include "htslib/faidx.h"
#include "wphmm.h"
#include "wphmm_ungapped.h"

// inline bool comparePtrToCandyMotify(candidate_fuzzy_motif* a, candidate_fuzzy_motif* b) {
//     return (a->viterbi_score > b->viterbi_score);
// }

struct candidate_fuzzy_motif {
    uint32_t mlen;
    std::string ru;
    std::vector<bool> label;
    bool inexact;
    int32_t st, ed; // zero based, inclusive. genome position of repeat region
    int32_t n_ru;   // number of (noninterrupted) ru matched
    double viterbi_score;
    int32_t l;
    bool insertion_relevant;
    std::string insertion;
    candidate_fuzzy_motif(std::string _u, std::vector<bool> _w, int32_t _s, int32_t _e, int32_t _n, double _v) : ru(_u), label(_w), st(_s), ed(_e), n_ru(_n), viterbi_score(_v) {
            mlen = ru.length();
            inexact = 0;
            for (uint32_t i = 0; i < mlen; ++i) {
                inexact = inexact || label[i];
            }
            l = ed-st+1;
            insertion_relevant=0;
        }
    candidate_fuzzy_motif() {}
    bool operator<(const candidate_fuzzy_motif & rhs) const {
        return (viterbi_score > rhs.viterbi_score);
    }

    int32_t ru_compare(candidate_fuzzy_motif& rhs) {
        if (ru.length() != rhs.ru.length()) {
            return 3; // Repeat unit incompatible
        }
        if (ru == rhs.ru && label == rhs.label) {
            return 0; // Equivalent
        }
        bool equal_seq = (ru == rhs.ru);
        uint32_t phase = 1;
        while (phase < mlen && (!equal_seq)) {
            std::string cru = ru.substr(phase, mlen-phase)+ru.substr(0,phase);
            if (cru == rhs.ru) {
                equal_seq = 1;
                break;
            }
            phase++;
        }
        if (!equal_seq) {
            return 3; // Repeat unit incompatible
        }
        bool equal_cir = 1;
        for (int32_t i = 0; i < mlen; ++i) {
            if (label[(i+phase)%mlen] != rhs.label[i]) {
                equal_cir = 0;
                break;
            }
        }
        if (equal_cir) {
            return 0; // Circular equivalent
        }
        // Compatible, pick the model with higher score
        if (viterbi_score >= rhs.viterbi_score) {
            return 1; // left unit has higher score
        } else {
            return 2; //
        }
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
class VNTR_candidate
{
    public:

    bool debug;
    faidx_t* fai;
    int32_t rid;
    candidate_fuzzy_motif motif;   // store st/ed/len w.r.t. reference
    int32_t st, ed, mlen; // not equal to motif if multiple candidate models are merged
    int32_t rlen_r, n_ru_r; // info w.r.t. reference seq
    int32_t rlen_l, n_ru_l; // info w.r.t. longest mosaic sequence
    double score_r, score_l;
    const char* chrom;
    std::string repeat_ref;     // repeat region in ref
    std::string lflank;         // left flank
    std::string rflank;         // right flank
    std::list<std::pair<int32_t, std::string> > insertions; // maintain sorted
    // std::set<int32_t> inserted_pos;
    std::string merged_longest_rr; // could make this easier by keeping only the longest allele
    bool need_refit;
    // Criteria to decide if merge two candidates
    double critical_ovlp;
    double max_interrupt;
    std::vector<candidate_fuzzy_motif> alternative_models;

    VNTR_candidate(faidx_t* _f, bcf_hdr_t* h, bcf1_t* v, candidate_fuzzy_motif& _m, bool _b);
    VNTR_candidate() {}

    bool operator<(const VNTR_candidate& rhs) const {
        if (motif.st == rhs.motif.st) {
            return (motif.ed <= rhs.motif.ed);
        }
        return (motif.st < rhs.motif.st);
    }

    void set_merging_criteria(double _c, double _i) {
        critical_ovlp = _c;
        max_interrupt = _i;
    }

    void add_insertion(int32_t _p, std::string _s);
    void combine_insertions();
    int32_t get_pos_in_ref(int32_t p_rel);

    /**
     * Merge two VNTR record if possible.
       0 if not merged
       1 if merged and the query does not need to seek further
       2 if merged but the query needs to check with others in the buffer
     */
    int32_t merge(VNTR_candidate& rt);

    /**
     * Number of overlapped bases.
       0 for adjacent
       < 0 for nonoverlapping, with absolute value being the distance
     */
    int32_t intersect(VNTR_candidate& rt);

    /**
     * If two repeat models are compatible.
       0 - Equivalent; 3 - Incompatible;
       1/2 Compatible, 1 Left has higher score.
     */
    int32_t ru_compare(candidate_fuzzy_motif& rhs);

    double fit_alt_model(candidate_fuzzy_motif& rhs);

    /**
     * Fit ungapped local alignment model
     */
    void model_refit();
    bool choose_model();
    void evaluate_models();
    void evaluate_models(uint32_t index);

    /**
     * Prepare this TR for output
     */
    bool finalize(int32_t flk = 10);

    /**
     * Print object.
     */
    void print();
};
#endif
