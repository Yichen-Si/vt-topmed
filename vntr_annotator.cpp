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

#include "vntr_annotator.h"

/**
 * Constructor.
 */
VNTRAnnotator::VNTRAnnotator(std::string& ref_fasta_file, bool debug)
{
    vm = new VariantManip(ref_fasta_file.c_str());

    float delta = 0.0001;
    float epsilon = 0.05;
    float tau = 0.01;
    float eta = 0.01;
    float mismatch_penalty = 3;

    ahmm = new AHMM(false);
    ahmm->set_delta(delta);
    ahmm->set_epsilon(epsilon);
    ahmm->set_tau(tau);
    ahmm->set_eta(eta);
    ahmm->set_mismatch_penalty(mismatch_penalty);
    ahmm->initialize_T();

    fai = fai_load(ref_fasta_file.c_str());
    if (fai==NULL)
    {
        fprintf(stderr, "[%s:%d %s] Cannot load genome index: %s\n", __FILE__, __LINE__, __FUNCTION__, ref_fasta_file.c_str());
        exit(1);
    }

    cre = new CandidateRegionExtractor(ref_fasta_file, debug);
    fd = new FlankDetector(ref_fasta_file, debug);

    max_mlen = 10;
    mt = new MotifTree(max_mlen, debug);

    this->debug = debug;
    qual.assign(2048, 'K');

    min_ecover_indel = 0.5; min_pcover_indel = 0.5;
    min_ecover_extended = 0.3; min_pcover_extended = 0.3;
};

/**
 * Destructor.
 */
VNTRAnnotator::~VNTRAnnotator()
{

if (debug) {
    std::cerr << "VNTRAnnotator Destructor\n";
}
    delete vm;
    fai_destroy(fai);
    if (factors)
    {
        for (size_t i=1; i<=max_len; ++i)
        {
            if (factors[i])
                free(factors[i]);
        }
        free(factors);
    }
}

/**
 * Annotates VNTR characteristics.
 * @mode -
 */
void VNTRAnnotator::annotate(bcf_hdr_t* h, bcf1_t* v, Variant& variant, std::string mode)
{
    if (!(variant.type&VT_INDEL)) {
        return;
    }

    std::set<candidate_unit> candidate_ru;
    // VNTR& vntr = variant.vntr;
    variant.rid = bcf_get_rid(v);
    variant.pos1 = bcf_get_pos1(v);

    //1. detect candidate repeat units
    //2. for each candidate, detect repeat region and evaluate
    //3. choose the best repeat unit and track

    // if (debug) std::cerr << "============================================\n";
    // if (debug) std::cerr << "ANNOTATING INDEL\n";
    //1. detect candidate repeat units
    find_repeat_unit(h, v, candidate_ru);

    // if (debug) std::cerr << "============================================\n";
    return;

}

void VNTRAnnotator::find_repeat_unit(bcf_hdr_t* h, bcf1_t* v, std::set<candidate_unit>& candidate_ru) {

    int32_t pos1 = v->pos + 1;
    int32_t st_rel = 0, ed_rel = 0;
    const char* chrom = bcf_get_chrom(h, v);
    char** alleles = bcf_get_allele(v);
    std::string focal_seq(alleles[0]);
    int32_t indel_end = pos1 + focal_seq.size() - 1; // 1-based, inclusive
    if (strlen(alleles[1]) > focal_seq.size()) {
        focal_seq = alleles[1];
        // indel_end = pos1;
    }
    focal_seq = focal_seq.substr(1); // only keep the inserted or deleted seq
if (debug) {
    std::cerr << "============================================\n";
    std::cerr << "Read indel: " << pos1 << '\t' << alleles[0] << '\t' << alleles[1] << '\t' << focal_seq << std::endl;
}
    // check if it is an indel inside homopolymer
    char b;
    int32_t homopoly = if_homopoly(chrom, pos1-1, indel_end, b);
    if (homopoly > HOMOPOLYMER_MIN) {
        if (focal_seq.size() < 4) {
            bcf_add_filter(h, v, bcf_hdr_id2int(h, BCF_DT_ID, "overlap_homopolymer"));
            return;
        } else {
            // Check if the indel is homopolymer itself
            int32_t same_base = 0;
            for (size_t i = 0; i < focal_seq.size(); ++i) {
                if (focal_seq.at(i) == b) {same_base++;}
            }
            if (same_base == focal_seq.size()) {
                bcf_add_filter(h, v, bcf_hdr_id2int(h, BCF_DT_ID, "overlap_homopolymer"));
                return;
            }
        }
    }
    if (focal_seq.size() == 1) {
        return;
    }

    // extend in both direction (small region)
    int32_t slen = std::min((int32_t)focal_seq.size() * 6, 64);
    int32_t ct_indel = 0;
    ct_indel += og_find_repeat_unit(focal_seq, candidate_ru, 1);
    ct_indel += rl_find_repeat_unit(focal_seq, candidate_ru, 1);
    if (focal_seq.size() < slen) {
        int32_t pad = (int32_t) ((slen - focal_seq.size() + 1) / 2);
        int32_t seq_len;
        std::string extended_seq = faidx_fetch_seq(fai, chrom, std::max(0, pos1 - pad), pos1 - 1, &seq_len);
        st_rel = extended_seq.size();
        extended_seq += focal_seq;
        ed_rel = extended_seq.size() - 1;
        extended_seq += faidx_fetch_seq(fai, chrom, indel_end, indel_end + pad - 1, &seq_len);
if (debug) {
    std::cerr << "Extended sequence: " << slen << '\t' << focal_seq << '\t' << extended_seq << '\t' << st_rel  << ',' << ed_rel << std::endl;
}
        int32_t ct_extnd = 0;
        ct_extnd += rl_find_repeat_unit(extended_seq, candidate_ru, 0);
        ct_extnd += og_find_repeat_unit(extended_seq, candidate_ru, 0);
        if (ct_indel == 0 && focal_seq.size() > 4) {
            if (extended_seq.substr(ed_rel+1,focal_seq.size()).compare(focal_seq) == 0 ||
                extended_seq.find(focal_seq) <= st_rel-focal_seq.size()) {
                    candidate_ru.insert(candidate_unit(focal_seq));
                }
        }
    }


if (debug) {
    std::cerr << "Identified " << candidate_ru.size() << " candidates\n";
    for (const auto& s : candidate_ru) {
        std::cerr << s.ru << ' ' << s.inexact << '\t';
    }
    std::cerr << std::endl;
}
    return;
}

/**
 * Find periodic subsequence
 */
int32_t VNTRAnnotator::og_find_repeat_unit(std::string& context, std::set<candidate_unit>& candidate_ru, bool flag) {
    int32_t added_ct = 0;
    if (context.size() > 2) {
        periodic_seq pseq_obj(context.c_str(), max_mlen, debug);
        if (flag) {
            pseq_obj.get_candidate(min_ecover_indel, min_pcover_indel);
        } else {
            pseq_obj.get_candidate(min_ecover_extended, min_pcover_extended);
        }
        if (pseq_obj.candidate.size() > 0) {
            for (const auto& s: pseq_obj.candidate) {
                if (context.find(s) == std::string::npos) {
                    continue;
                }
                candidate_unit tmp(s,0);
                candidate_ru.insert(tmp);
                added_ct++;
            }
        }
    }
    return added_ct;
}

/**
 * Find periodic subsequence ignoring repeats of a single base
 */
int32_t VNTRAnnotator::rl_find_repeat_unit(std::string& context, std::set<candidate_unit>& candidate_ru, bool flag) {

    int32_t added_ct = 0;
    std::string cseq = "";
    cseq += context.at(0);
    std::vector<int32_t> rl;
    char pre_base = context.at(0);
    std::map<char, int32_t> base_ct;
    for (auto& s : alphabet) {
        base_ct[s] = 0;
    }
    int32_t uninfo = 0;
    int32_t ct = 1;
    for (uint32_t i = 1; i < context.size(); ++i) {
        if (context.at(i) != pre_base) {
            if (base_ct.find(context.at(i))==base_ct.end()) {
                uninfo++;
                continue;
            }
            cseq += context.at(i);
            base_ct[context.at(i)] += 1;
            rl.push_back(ct);
            pre_base = context.at(i);
            ct = 1;
        } else {ct++;}
    }
    rl.push_back(ct);
    int32_t N = cseq.size();

    if (cseq.size() == 1) {
        if (context.size() >= HOMOPOLYMER_MIN) {
            candidate_unit tmp(cseq);
            candidate_ru.insert(tmp);
            added_ct++;
        }
    } else {
        if (cseq.at(0) == cseq.at(cseq.size()-1)) {
            cseq = cseq.substr(1);
        }
        if (cseq.size() > 2) {
            periodic_seq pseq_obj(cseq.c_str(), max_mlen, debug);
            if (flag) {
                pseq_obj.get_candidate(min_ecover_indel, min_pcover_indel);
            } else {
                pseq_obj.get_candidate(min_ecover_extended, min_pcover_extended);
            }
            if (pseq_obj.candidate.size() > 0) {
                for (const auto& s: pseq_obj.candidate) {
                    // Convert back from collapsed sequence signature
                    candidate_unit tmp(s,1);
                    size_t pt = cseq.find(s);
                    std::vector<int32_t> omin(s.size(), cseq.size());
                    std::vector<int32_t> omax(s.size(), 0);
                    if (pt == std::string::npos) {
                        continue;
                    }
                    while (pt != std::string::npos) {
                        for (size_t i = 0; i < s.size(); ++i) {
                            omin[i] = std::min(omin[i], rl[pt+i]);
                            omax[i] = std::max(omax[i], rl[pt+i]);
                        }
                        pt = cseq.find(s, pt+s.size());
                    }
                    std::string ss;
                    for (size_t i = 0; i < s.size(); ++i) {
                        if (omin[i] == omax[i]) {
                            for (int32_t j = 0; j < omin[i]; ++j) {
                                ss += s.at(i);
                                tmp.variable_base.push_back(0);
                            }
                        } else {
                            ss += s.at(i);
                            tmp.variable_base.push_back(omin[i]);
                        }
                    }
                    tmp.ru = ss;
                    tmp.check();
                    candidate_ru.insert(tmp);
                    added_ct++;
                }
            }
        }
    }
    return added_ct;
}

int32_t VNTRAnnotator::if_homopoly(const char* chrom, int32_t left, int32_t right, char&b) {
    int32_t seq_len;
    std::string hseq = faidx_fetch_seq(fai, chrom, std::max(0, left-2*HOMOPOLYMER_MIN), left, &seq_len);
    uint32_t st_rel = hseq.size();
    hseq += faidx_fetch_seq(fai, chrom, right, right+2*HOMOPOLYMER_MIN, &seq_len);
    b = hseq.at(st_rel);
    int32_t ct = 0;
    for (uint32_t i = st_rel; i < hseq.size(); ++i) {
        if (hseq.at(i) == b) {
            ct++;
        } else {
            break;
        }
    }
    for (uint32_t i = st_rel-1; i > 0; --i) {
        if (hseq.at(i) == b) {
            ct++;
        } else {
            break;
        }
    }
    return ct;
}

/**
 * Returns true if is to be classified as an STR
 */
bool VNTRAnnotator::is_vntr(Variant& variant, int32_t mode, std::string& method)
{
    uint32_t mlen = 0;
    uint32_t rlen = 0;
    float motif_concordance = 0;
    uint32_t no_exact_ru = 0;

    if (method == "e")
    {
        mlen = variant.vntr.mlen;
        rlen = variant.vntr.rl;
        motif_concordance = variant.vntr.motif_concordance;
        no_exact_ru = variant.vntr.no_exact_ru;
    }
    else if (method == "f")
    {
        mlen = variant.vntr.mlen;
        rlen = variant.vntr.fuzzy_rl;
        motif_concordance = variant.vntr.fuzzy_motif_concordance;
        no_exact_ru = variant.vntr.fuzzy_no_exact_ru;
    }

    if (mode==TAN_KANG_2015_VNTR)
    {
        if ((rlen - mlen) >= 6 && no_exact_ru>=2)
        {
            if (mlen==1 && motif_concordance>0.9)
            {
                return true;
            }
            else if (mlen>1 || motif_concordance>0.75)
            {
                return true;
            }
        }

        return false;
    }
    else if (mode==WILLEMS_2014_STR)
    {
        return ((mlen==1 && rlen>=6) ||
                (mlen==2 && rlen>=11) ||
                (mlen==3 && rlen>=14) ||
                (mlen==4 && rlen>=14) ||
                (mlen==5 && rlen>=16) ||
                (mlen==6 && rlen>=17) ||
                (mlen>=7 && rlen>=mlen*2));
    }
    else if (mode==ANANDA_2013_STR)
    {
        return ((mlen==1 && rlen>=2) ||
                (mlen==2 && rlen>=4) ||
                (mlen==3 && rlen>=6) ||
                (mlen==4 && rlen>=8) ||
                (mlen==5 && rlen>=10) ||
                (mlen==6 && rlen>=12) ||
                (mlen>=7 && rlen>=mlen*2));
    }
    else if (mode==FONDON_2012_STR)
    {
        return ((mlen==1 && rlen>=6) ||
                (mlen==2 && rlen>=13) ||
                (mlen==3 && rlen>=20) ||
                (mlen==4 && rlen>=23) ||
                (mlen==5 && rlen>=27) ||
                (mlen==6 && rlen>=27));
    }
    else if (mode==KELKAR_2008_STR)
    {
        return ((mlen==1 && rlen>=6) ||
                (mlen==2 && rlen>=10) ||
                (mlen==3 && rlen>=6) ||
                (mlen==4 && rlen>=8) ||
                (mlen==5 && rlen>=10) ||
                (mlen==6 && rlen>=12) ||
                (mlen>=7 && rlen>=mlen*2));
    }
    else if (mode==LAI_2003_STR)
    {
        return ((mlen==1 && rlen>=6) ||
                (mlen==2 && rlen>=8) ||
                (mlen==3 && rlen>=12) ||
                (mlen==4 && rlen>=16) ||
                (mlen==5 && rlen>=20) ||
                (mlen==6 && rlen>=24) ||
                (mlen>=7 && rlen>=mlen*2));
    }
    else
    {
        fprintf(stderr, "[%s:%d %s] STR definition mode does not exist: %d\n", __FILE__,__LINE__,__FUNCTION__, mode);
        exit(1);
    }
}
