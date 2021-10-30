/* The MIT License

   Copyright (c) 2014 Adrian Tan <atks@umich.edu>

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.
*/

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
#include "hts_utils.h"
#include "htslib/kstring.h"
#include "bcf_ordered_reader.h"
#include "bcf_ordered_writer.h"
#include "rfhmm.h"
#include "lfhmm.h"
#include "ahmm.h"
#include "variant_manip.h"
#include "program.h"
#include "motif_tree.h"
#include "vntr.h"
#include "candidate_region_extractor.h"
#include "flank_detector.h"
#include "seq_ipdft.h"

//definition of STRs
#define LAI_2003_STR       1
#define KELKAR_2008_STR    2
#define FONDON_2012_STR    3
#define ANANDA_2013_STR    4
#define WILLEMS_2014_STR   5
#define TAN_KANG_2015_VNTR 6


//forms of alignment
#define REFERENCE                                0
#define EXACT_LEFT_RIGHT_ALIGNMENT               1
#define FUZZY_LEFT_RIGHT_ALIGNMENT               2
#define FUZZY_LEFT_RIGHT_ALIGNMENT_WITH_PENALTY  3

//forms of choosing a motif
#define PICK_BEST_MOTIF             0

#define ALLELE_EXACT  1
#define ALLELE_FUZZY  2

#define CLIP_ENDS 0
#define CLIP_1L2R 1
#define FRAHMM    2

#define HOMOPOLYMER_MIN 6

/**
 * Class for determining basic traits of an indel
 * motifs, flanks and VNTR type statistics.
 * RU,RL,LFLANK,RFLANK,LFLANKPOS,RFLANKPOS,MOTIF_CONCORDANCE,MOTIF_CONCORDANCE
 */
class VNTRAnnotator
{
    public:

    uint32_t max_mlen; //maximum length for motif search in the fast tree.
    std::vector<char> alphabet{'A','C','G','T'};
    double min_ecover_indel, min_pcover_indel;
    double min_ecover_extended, min_pcover_extended;

    //model

    bool debug;
    int32_t max_len;

    ////////
    //raHMMs
    ////////
    std::string qual;
    AHMM* ahmm;

    ///////
    //tools
    ///////
    VariantManip *vm;
    faidx_t* fai;
    CandidateRegionExtractor* cre;
    MotifTree* mt;
    FlankDetector* fd;

    //for retrieving sequences
    int8_t* seq;

    //factors[n][index], for determining what sub repeat units to examine
    int32_t** factors;

    /**
     * Constructor.
     */
    VNTRAnnotator(std::string& ref_fasta_file, bool debug=false);

    /**
     * Destructor.
     */
    ~VNTRAnnotator();

    void set_ru_detection_threshold(double _e1, double _p1, double _e2, double _p2) {
        min_ecover_indel = _e1; min_pcover_indel = _p1;
        min_ecover_extended = _e2; min_pcover_extended = _p2;
    }

    /**
     * Annotates VNTR characteristics.
     * @mode
     *   e - determine by exact alignment
     *   f - determine by fuzzy alignment
     *   p - determine by penalized fuzzy alignment
     *   h - using HMMs
     *   x - integrated models
     */
    void annotate(bcf_hdr_t* h, bcf1_t* v, Variant& variant, std::string mode);
    void find_repeat_unit(bcf_hdr_t* h, bcf1_t* v, std::set<candidate_unit>& candidate_ru);
    int32_t rl_find_repeat_unit(std::string& context, std::set<candidate_unit>& candidate_ru, bool flag=0);
    int32_t og_find_repeat_unit(std::string& context, std::set<candidate_unit>& candidate_ru, bool flag=0);
    int32_t if_homopoly(const char* chrom, int32_t left, int32_t right, char&b);

    /**
     * Returns true if is to be classified as a VNTR
     */
    bool is_vntr(Variant& variant, int32_t mode, std::string& method);
};

#endif
