#ifndef VNTR_ANNOTATOR_H
#define VNTR_ANNOTATOR_H

#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <complex>
#include <vector>
#include <map>
#include <queue>
#include <list>
#include <algorithm>
#include "hts_utils.h"
#include "htslib/kstring.h"
#include "bcf_ordered_reader.h"
#include "bcf_ordered_writer.h"
#include "program.h"

#include "variant_manip.h"
#include "vntr_candidate.h"
#include "wphmm.h"
#include "seq_ipdft.h"

//definition of STRs
#define LAI_2003_STR       1
#define KELKAR_2008_STR    2
#define FONDON_2012_STR    3
#define ANANDA_2013_STR    4
#define WILLEMS_2014_STR   5
#define TAN_KANG_2015_VNTR 6

#define HOMOPOLYMER_MIN 6

/**
 * Class for determining basic traits of a candidate VNTR
 * motifs, flanks and VNTR type statistics
 * RU,RL,LFLANK,RFLANK,LFLANKPOS,RFLANKPOS,MOTIF_CONCORDANCE,MOTIF_CONCORDANCE
 */
class VNTRAnnotator
{
    public:

    bool debug;
    uint32_t max_mlen; //maximum length for motif
    std::vector<char> alphabet{'A','C','G','T'};
    double min_ecover_indel, min_pcover_indel;
    double min_ecover_extended, min_pcover_extended;
    int32_t extend_bp; // distance in both direction to go when HMM

    // tools
    faidx_t* fai;
    VariantManip *vm;
    // TODO need a separate hmm operator for collapsed motif
    WPHMM* wphmm;

    VNTRAnnotator(std::string& ref_fasta_file, bool _debug=false, uint32_t _m = 16);
    ~VNTRAnnotator();

    void set_ru_detection_threshold(double _e1, double _p1, double _e2, double _p2) {
        min_ecover_indel = _e1; min_pcover_indel = _p1;
        min_ecover_extended = _e2; min_pcover_extended = _p2;
    }

    void annotate(bcf_hdr_t* h, bcf1_t* v, std::string mode);
    /**
     * Try to find VNTR signatures.
     */
    void find_repeat_unit(bcf_hdr_t* h, bcf1_t* v, std::set<candidate_unit>& candidate_ru);
    int32_t rl_find_repeat_unit(std::string& context, std::set<candidate_unit>& candidate_ru, bool flag=0);
    int32_t og_find_repeat_unit(std::string& context, std::set<candidate_unit>& candidate_ru, bool flag=0);
    int32_t if_homopoly(const char* chrom, int32_t left, int32_t right, char&b);
    /**
     * Try to find the boundary of repeat region given a set of candidate RU
     */
     void find_repeat_region(bcf_hdr_t* h, bcf1_t* v, std::set<candidate_unit>& candidate_ru);

    /**
     * Returns true if is to be classified as a VNTR
     */
    bool is_vntr(Variant& variant, int32_t mode, std::string& method);
};

#endif
