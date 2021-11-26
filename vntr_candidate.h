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
#include "ru_candidate.h"

/**
 * Class for representing a VNTR.
 * Hold information of multiple repeat models
 * perhaps generated from different indels
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
    int32_t rlen_r2, n_ru_r2;
    double score_r, score_l;
    double score_r2, concordance_r, concordance_r2;
    const char* chrom;
    std::string repeat_ref;     // repeat region in ref
    std::string lflank;         // left flank
    std::string rflank;         // right flank
    std::list<std::pair<int32_t, std::string> > insertions; // maintain sorted
    // std::set<int32_t> inserted_pos;
    int32_t rel_st, len_mrg, ins_tot;
    std::string query;
    // std::string merged_longest_rr; // could make this easier by keeping only the longest allele
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
    void combine_insertions(int32_t pad = -1);
    int32_t get_pos_in_ref(int32_t p_rel, int32_t rel_st = 0);

    /**
     * Merge two VNTR record if possible.
     * 0 if not merged
     * 1 if merged and the query does not need to seek further
     * 2 if merged but the query needs to check with others in the buffer
     */
    int32_t merge(VNTR_candidate& rt);

    /**
     * Number of overlapped bases.
       0 for adjacent
       < 0 for nonoverlapping, with absolute value being the distance
     */
    int32_t intersect(VNTR_candidate& rt);

    /**
     * Fit an alternative RU to the current sequence
     */
    bool fit_alt_model(candidate_fuzzy_motif& rhs);

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
