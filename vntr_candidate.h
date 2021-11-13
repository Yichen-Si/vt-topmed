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
    candidate_fuzzy_motif(std::string _u, std::vector<bool> _w, int32_t _s, int32_t _e, int32_t _n, double _v) : ru(_u), label(_w), st(_s), ed(_e), n_ru(_n), viterbi_score(_v) {
            mlen = ru.length();
            inexact = 0;
            for (uint32_t i = 0; i < mlen; ++i) {
                inexact = inexact || label[i];
            }
            l = ed-st+1;
        }
    candidate_fuzzy_motif() {}
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
class VNTR_candidate
{
    public:

    bool debug;
    int32_t rid;
    candidate_fuzzy_motif motif;   // store st/ed/len w.r.t. reference
    uint32_t mlen, rlen_r, n_ru_r; // duplicated info
    uint32_t rlen_l, n_ru_l;       // info w.r.t. longest mosaic sequence
    double score_r;
    const char* chrom;
    // std::vector<std::pair<std::string, std::vector<bool> > > alternative_ru;
    std::string repeat_ref;     // repeat region in ref
    std::string lflank;         // left flank
    std::string rflank;         // right flank
    std::list<std::pair<int32_t, std::string> > insertions; // maintain sorted
    // std::set<int32_t> inserted_pos;
    std::string merged_longest_rr; // could make this easier by keeping only the longest allele
    bool updated = 0;

    VNTR_candidate(bcf_hdr_t* h, int32_t _r, candidate_fuzzy_motif& _m, bool _b) : debug(_b), rid(_r), motif(_m) {
        mlen   = motif.mlen;
        rlen_l = motif.l;
        rlen_r = motif.ed - motif.st + 1;
        n_ru_r = motif.n_ru;
        n_ru_l = motif.n_ru;
        chrom = h->id[BCF_DT_CTG][rid].key;
    }
    VNTR_candidate() {}

    void add_insertion(int32_t _p, std::string _s);
    void combined_insertions(faidx_t* fai);
    int32_t get_pos_in_ref(int32_t p_rel);

    /**
     * Merge two VNTR record if possible. (1 if merge successfully)
     */
    int32_t merge(VNTR_candidate& rt, int32_t min_overlap = 1);

    /**
     * 0: overlap, -1: <rt, 1: >rt.
     */
    int32_t intersect(VNTR_candidate& rt, int32_t min_overlap = 1);

    /**
     * If two repeat models are compatible.
       0 - Equivalent; 3 - Incompatible;
       1 - Left is better; 2 - Right is better.
     */
    int32_t ru_compair(candidate_fuzzy_motif& rhs);

    bool operator<(const VNTR_candidate& rhs) const {
        if (motif.st == rhs.motif.st) {
            return (motif.ed <= rhs.motif.ed);
        }
        return (motif.st < rhs.motif.st);
    }

    /**
     * Fit an ungapped local alignment if updates have been made
     */
    bool model_refit(faidx_t* fai, int32_t flk=10);

    /**
     * Print object.
     */
    void print();
};
#endif
