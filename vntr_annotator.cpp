#include "vntr_annotator.h"

VNTRAnnotator::VNTRAnnotator(std::string& ref_fasta_file, bool _debug, uint32_t _m) : debug(_debug), max_mlen(_m) {
    vm = new VariantManip(ref_fasta_file.c_str());
    fai = fai_load(ref_fasta_file.c_str());
    if (fai==NULL) {
        fprintf(stderr, "[%s:%d %s] Cannot load genome index: %s\n", __FILE__, __LINE__, __FUNCTION__, ref_fasta_file.c_str());
        exit(1);
    }
    min_ecover_indel = 0.5; min_pcover_indel = 0.5;
    min_ecover_extended = 0.3; min_pcover_extended = 0.3;
    extend_bp = 150;
    // wphmm = NULL;
};

VNTRAnnotator::~VNTRAnnotator() {
    delete vm;
    vm = NULL;
    // delete wphmm;
    // wphmm = NULL;
    fai_destroy(fai);
}

/**
 * Annotates VNTR characteristics.
 */
int32_t VNTRAnnotator::annotate(bcf_hdr_t* h, bcf1_t* v, std::vector<candidate_fuzzy_motif>& candidate_model, int32_t mode)
{
    if (bcf_get_variant_types(v) != VCF_INDEL) {return 0;}
    std::set<candidate_unit> candidate_ru;
    // if (debug) std::cerr << "============================================\n";
    // if (debug) std::cerr << "ANNOTATING INDEL\n";
    //1. detect candidate repeat units
    find_repeat_unit(h, v, candidate_ru);
    //2. for each candidate, detect repeat region (boundary) and evaluate (HMM)
    find_repeat_region(h, v, candidate_ru, candidate_model);
    //3. keep the best repeat model, check if qualified as VNTR
    pick_top_candidates(candidate_model, mode);

    return candidate_model.size();
}

void VNTRAnnotator::pick_top_candidates(std::vector<candidate_fuzzy_motif>& candidate_model, int32_t mode) {
    if (candidate_model.size() == 0) {
        return;
    }
    std::sort(candidate_model.begin(), candidate_model.end());
    auto it = candidate_model.begin();
    bool ex = 0, fz = 0;
    while (it != candidate_model.end()) {
        if (is_vntr(*it, mode)) {
            if (it->inexact && (!fz)) {
                fz = 1;
                it++;
                continue;
            }
            if ((!it->inexact) && (!ex)) {
                ex = 1;
                it++;
                continue;
            }
        }
        it = candidate_model.erase(it);
    }
if (debug) {
    std::cerr << "VNTRAnnotator::pick_top_candidates " << ex << ' ' << fz << ' ' << candidate_model.size() << '\n';
}
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
    }
    focal_seq = focal_seq.substr(1); // only keep the inserted or deleted seq
    // check if it is an indel inside homopolymer
    char b;
    int32_t homopoly = if_homopoly(chrom, v->pos, indel_end, b);
    if (homopoly >= HOMOPOLYMER_MIN) {
        std::string s(1, b);
        candidate_unit ru(s, 0);
        candidate_ru.insert(ru);
    }
    // extend in both direction (small region)
    int32_t slen = std::min((int32_t)focal_seq.size() * 6, 64);
    if (focal_seq.size() == 1) {slen *= 2;}
    int32_t ct_indel = 0;
    ct_indel += og_find_repeat_unit(focal_seq, candidate_ru, 1);
    ct_indel += rl_find_repeat_unit(focal_seq, candidate_ru, 1);
    if (focal_seq.size() < slen) {
        int32_t pad = (int32_t) ((slen - focal_seq.size() + 1) / 2);
        int32_t seq_len;
        int32_t extended_st0 = std::max(0, pos1 - pad);
        int32_t extended_ed0 = indel_end + pad - 1;
        std::string extended_seq = faidx_fetch_seq(fai, chrom, extended_st0, pos1 - 1, &seq_len);
        st_rel = extended_seq.size();
        extended_seq += focal_seq;
        ed_rel = extended_seq.size() - 1;
        extended_seq += faidx_fetch_seq(fai, chrom, indel_end, extended_ed0, &seq_len);
        int32_t ct_extnd = 0;
        ct_extnd += rl_find_repeat_unit(extended_seq, candidate_ru, 0);
        ct_extnd += og_find_repeat_unit(extended_seq, candidate_ru, 0);
        if (ct_indel == 0 && focal_seq.size() > 4 && focal_seq.size() <= max_mlen) {
            if (extended_seq.substr(ed_rel+1,focal_seq.size()).compare(focal_seq) == 0 ||
                extended_seq.find(focal_seq) <= st_rel-focal_seq.size()) {
                candidate_ru.insert(candidate_unit(focal_seq));
            }
        }
if (debug && candidate_ru.size() > 0) {
    std::cerr << "============================================\n";
    std::cerr << "Read indel: " << pos1 << '\t' << alleles[0] << '\t' << alleles[1] << std::endl;
    std::cerr << "Extended sequence: " << slen << '\t' << focal_seq << '\t' << extended_seq << std::endl;
}
    }
if (debug && candidate_ru.size() > 0) {
    std::cerr << "Identified " << candidate_ru.size() << " candidates:\t";
    if (candidate_ru.size() > 0) {
        for (const auto& s : candidate_ru) {
            std::cerr << s.ru << ' ' << s.inexact << '\t';
        }
        std::cerr << std::endl;
    }
}
    return;
}


void VNTRAnnotator::find_repeat_region(bcf_hdr_t* h, bcf1_t* v, std::set<candidate_unit>& candidate_ru, std::vector<candidate_fuzzy_motif>& candidate_model) {
    if (candidate_ru.size() < 1) {return;}
    int32_t pos1 = v->pos + 1;
    const char* chrom = bcf_get_chrom(h, v);
    int32_t start0 = 0, end0 = 0; // genome pos of the whole query
    int32_t spiked_st = 0, spiked_ed = 0; // relative pos of the longest allele
    int32_t seq_len;
    char** alleles = bcf_get_allele(v);
    std::string focal_seq(alleles[0]);
    int32_t indel_end = v->pos + focal_seq.size() - 1; // 0-based, inclusive

    // Take a large enough region
    start0 = (v->pos - extend_bp < 0) ? 0 : v->pos - extend_bp;
    end0   = indel_end + extend_bp;
    if (strlen(alleles[1]) > strlen(alleles[0])) {
        focal_seq = alleles[1];
    }
    focal_seq = focal_seq.substr(1);
    std::string query = faidx_fetch_seq(fai, chrom, start0, pos1-1, &seq_len);
    if (seq_len != extend_bp + 1) {
        start0 = pos1-seq_len;
    }
    spiked_st = query.size(); // 0-based rel. pos of the first base of inserted or deleted sequence
    query += focal_seq;
    spiked_ed = query.size() - 1; // 0-based rel. pos of the last base of the longest allele
    query += faidx_fetch_seq(fai, chrom, indel_end+1, end0, &seq_len);

// if (debug) {
//     std::cerr << "VNTRAnnotator::find_repeat_region - Query sequence:\n" << query << '\n';
// }

    // Process each candidate RU separately
    for (const auto& s : candidate_ru) {
        // Convert candidate repeat unit to pHMM motif, initialize hmm object
        std::string unit = s.ru;
        int32_t m = s.ru.size();
        std::string display_ru = unit;
        std::vector<bool> inexact_label(m, 0);
        std::vector<bool> tmp;
if (debug) {
    std::cerr << "VNTRAnnotator::find_repeat_region - Motif of HMM: " << unit;
}
        if (!s.inexact) {
            if (debug) {
                std::cerr << std::endl;
            }
            tmp = inexact_label;
            // wphmm = new WPHMM(query.c_str(), unit.c_str(), debug);
        } else {
            unit = "";
            for (size_t i = 0; i < s.ru.size(); ++i) {
                int32_t r = s.variable_base[i].first;
                for (int32_t j = 0; j < r; ++j) {
                    unit += s.ru.at(i);
                    tmp.push_back(0);
                }
                if (r < s.variable_base[i].second) {
                    unit += s.ru.at(i);
                    tmp.push_back(1);
                    inexact_label[i] = 1;
                }
            }
            m = unit.length();
            if (debug) {
                std::cerr << ", expanded as " << unit << '\t';
                for (const auto & v : tmp) {std::cerr << v;}
                std::cerr << std::endl;
            }
        }
        bool base_relax[m];
        std::copy(tmp.begin(), tmp.end(), base_relax);
        WPHMM* wphmm = new WPHMM(query.c_str(), unit.c_str(), debug, base_relax);
        // // Run HMM
        wphmm->set_ru(display_ru, inexact_label);
        wphmm->initialize();
        wphmm->viterbi();
        wphmm->detect_range();
        if (wphmm->select_segment(spiked_st, spiked_ed)) {
            seq_segment& rr = wphmm->focal_rr;
            int32_t real_st = start0 + rr.p_st;
            int32_t real_ed = start0 + rr.p_ed;
            if (indel_end - v->pos < focal_seq.size()) {
                real_ed -= focal_seq.size();
            }
            candidate_fuzzy_motif model(display_ru, inexact_label, real_st, real_ed, rr.nr, rr.score);
            candidate_model.push_back(model);
        }
        delete wphmm;
        wphmm = NULL;
    }
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
                candidate_ru.insert(candidate_unit(s,0));
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
            candidate_ru.insert(candidate_unit(cseq));
            added_ct++;
        }
    } else {
        if (cseq.at(0) == cseq.at(cseq.size()-1)) {
            cseq = cseq.substr(0, cseq.size()-1);
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
                                tmp.variable_base.push_back(std::make_pair(1,1));
                            }
                        } else {
                            ss += s.at(i);
                            tmp.variable_base.push_back(std::make_pair(omin[i], omax[i]));
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
    b = hseq.at(st_rel-1);
    if (hseq.at(st_rel) != b) {
        char bl, br;
        int32_t ctl = 0, ctr = 0;
        bl = hseq.at(st_rel-1);
        br = hseq.at(st_rel);
        for (uint32_t i = st_rel; i < hseq.size(); ++i) {
            if (hseq.at(i) == br) {
                ctr++;
            } else {
                break;
            }
        }
        for (uint32_t i = st_rel-1; i > 0; --i) {
            if (hseq.at(i) == bl) {
                ctl++;
            } else {
                break;
            }
        }
        if (ctl >= ctr) {
            b = bl;
            return(ctl);
        }
        b = br;
        return(ctr);
    }
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
 bool VNTRAnnotator::is_vntr(candidate_fuzzy_motif& motif, int32_t mode)
 {
     uint32_t mlen = motif.mlen;
     uint32_t rlen = motif.ed - motif.st + 1;

     if (mode==VITERBI_VNTR) {
         return ((rlen - mlen) >= 6 && motif.n_ru>=2 && motif.viterbi_score >= 3.);
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
