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
int32_t VNTRAnnotator::annotate(bcf_hdr_t* h, bcf1_t* v, std::vector<candidate_fuzzy_motif>& candidate_model, int32_t mode) {
    if (bcf_get_variant_types(v) != VCF_INDEL) {return 0;}
    std::set<candidate_unit> candidate_ru;
    //1. detect candidate repeat units
    find_repeat_unit(h, v, candidate_ru);
    //2. for each candidate, detect repeat region (boundary) and evaluate (HMM)
    find_repeat_region(h, v, candidate_ru, candidate_model);
    //3. keep the best repeat models, check if qualified as VNTR
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
    double max_score = it->viterbi_score;
    while (it != candidate_model.end()) {
        bool dd = is_vntr(*it, mode);
        if (dd && (it->viterbi_score > max_score*0.7) ) {
            if (it->viterbi_score > max_score) {
                max_score = it->viterbi_score;
            }
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
    std::cerr << "VNTRAnnotator::pick_top_candidates " << ex << ' ' << fz << ' ' << candidate_model.size() << '\t';
    if (candidate_model.size() != 0) {
        for (auto & v : candidate_model) {
            std::cerr << v.ru << '\t';
        }
    }
    std::cerr << '\n';
 }
}

void VNTRAnnotator::find_repeat_unit(bcf_hdr_t* h, bcf1_t* v, std::set<candidate_unit>& candidate_ru) {

    int32_t pos1 = v->pos + 1;
    int32_t st_rel = 0, ed_rel = 0;
    const char* chrom = bcf_get_chrom(h, v);
    char** alleles = bcf_get_allele(v);
    bool insertion = 0;
    std::string focal_seq(alleles[0]);
    int32_t indel_end = pos1 + focal_seq.size() - 1; // 1-based, inclusive
    if (strlen(alleles[1]) > focal_seq.size()) {
        insertion = 1;
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
    if (focal_seq.size() >= 6) {
        ct_indel += og_find_repeat_unit(focal_seq, candidate_ru, 1);
        ct_indel += rl_find_repeat_unit(focal_seq, candidate_ru, 1);
    }
    std::string extended_seq = focal_seq;
    if (focal_seq.size() < slen) {
        int32_t pad = (int32_t) ((slen - focal_seq.size() + 1) / 2);
        int32_t seq_len;
        int32_t extended_st0 = std::max(0, pos1 - pad);
        int32_t extended_ed0 = indel_end + pad - 1;
        extended_seq = faidx_fetch_seq(fai, chrom, extended_st0, extended_ed0, &seq_len);
        ct_indel += rl_find_repeat_unit(extended_seq, candidate_ru, 0);
        ct_indel += og_find_repeat_unit(extended_seq, candidate_ru, 0);
        if (insertion && focal_seq.length() > 3) {
            extended_seq = faidx_fetch_seq(fai, chrom, extended_st0, pos1 - 1, &seq_len);
            st_rel = extended_seq.size();
            extended_seq += focal_seq;
            ed_rel = extended_seq.size() - 1;
            extended_seq += faidx_fetch_seq(fai, chrom, indel_end, extended_ed0, &seq_len);
            ct_indel += rl_find_repeat_unit(extended_seq, candidate_ru, 0);
            ct_indel += og_find_repeat_unit(extended_seq, candidate_ru, 0);
            if (focal_seq.size() > 4 && focal_seq.size() <= max_mlen) {
                candidate_unit tmp(focal_seq);
                if (candidate_ru.find(tmp) == candidate_ru.end()) {
                    std::set<char> s1;
                    std::for_each(focal_seq.begin(), focal_seq.end(), [&s1] (char c) -> void { s1.insert(c);});
                    if (s1.size() > 1) {
                        if (extended_seq.substr(ed_rel+1,focal_seq.size()).compare(focal_seq) == 0 ||
                            extended_seq.find(focal_seq) == st_rel-focal_seq.size()) {
                            tmp.check();
                            tmp.reduce();
                            candidate_ru.insert(tmp);
                        }
                    }
                }
            }
        }
    }
    if (debug && candidate_ru.size() > 0) {
        std::cerr << "============================================\n";
        std::cerr << "Read indel: " << pos1 << '\t' << alleles[0] << '\t' << alleles[1] << '\t' << extended_seq << std::endl;
        std::cerr << "Candidate RU:\t";
        for (auto & v : candidate_ru) {
            std::cerr << v.ru << '\t';
        }
        std::cerr << '\n';
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
    bool is_insertion = 0;
    char** alleles = bcf_get_allele(v);
    std::string focal_seq(alleles[0]);
    int32_t indel_end = v->pos + focal_seq.size() - 1; // 0-based, inclusive

    // Take a large enough region
    start0 = v->pos - extend_bp;
    end0   = indel_end + extend_bp;
    if (strlen(alleles[1]) > strlen(alleles[0])) {
        focal_seq = alleles[1];
        is_insertion = 1;
    }
    focal_seq = focal_seq.substr(1);
    std::string query = faidx_fetch_seq(fai, chrom, start0, pos1-1, &seq_len);
    start0 = pos1-seq_len;

    spiked_st = query.size(); // 0-based rel. pos of the first base of inserted or deleted sequence
    query += focal_seq;
    spiked_ed = query.size() - 1; // 0-based rel. pos of the last base of the longest allele
    query += faidx_fetch_seq(fai, chrom, indel_end+1, end0, &seq_len);
    // Process each candidate RU separately
    for (const auto& s : candidate_ru) {
        // Convert candidate repeat unit to pHMM motif, initialize hmm object
        std::string unit = s.ru;
        int32_t m = s.ru.size();
        std::string display_ru = unit;
        std::vector<bool> inexact_label(m, 0);
        std::vector<bool> tmp;
        if (!s.inexact) {
            tmp = inexact_label;
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
        }
        bool base_relax[m];
        std::copy(tmp.begin(), tmp.end(), base_relax);
        WPHMM* wphmm = new WPHMM(query.c_str(), unit.c_str(), debug, base_relax);
        // Run HMM
        wphmm->set_ru(display_ru, inexact_label);
        wphmm->initialize();
        wphmm->viterbi();
        // std::string vpath = wphmm->print_viterbi_path();
        // wphmm->count_ru();
        wphmm->detect_range();
        if (wphmm->select_segment(spiked_st, spiked_ed)) {
            seq_segment& rr = wphmm->focal_rr;
            int32_t real_st = start0 + rr.p_st;
            int32_t real_ed = start0 + rr.p_ed;
            int32_t l = rr.p_ed - rr.p_st + 1;
            if (is_insertion && rr.p_ed > spiked_ed) {
                real_ed -= focal_seq.length();
            }
            if (is_insertion && rr.p_ed < spiked_ed && rr.p_ed >= spiked_st) {
                real_ed = v->pos;
            }
            // if (real_st < v->pos) {
            //     real_st = v->pos;
            // }
            candidate_fuzzy_motif model(display_ru, inexact_label, real_st, real_ed, rr.nr, rr.score, l, rr.match_motif);
            if (is_insertion) { // Decide if the inserted sequence is relevant (need to be explained by final model) to the TR.
                if (wphmm->ru_complete[spiked_ed] - wphmm->ru_complete[spiked_st-1] > 1) {
                    model.insertion_relevant = 1;
                }
                std::set<char> s1;
                int32_t nn = 0;
                std::for_each(model.ru.begin(), model.ru.end(), [&s1] (char c) -> void { s1.insert(c);});
                for (int32_t pt = 0; pt < focal_seq.size(); ++pt) {
                    if (s1.find(focal_seq.at(pt))==s1.end()) {
                        nn++;
                    }
                }
                if (nn>0.2*focal_seq.size()) {
                    model.insertion_relevant = 0;
                } else {
                    for (int32_t pt = 0; pt < model.mlen; ++pt) {
                        std::string cc = model.ru.substr(pt)+model.ru.substr(0,pt);
                        if (focal_seq.find(cc) != std::string::npos || cc.find(focal_seq) != std::string::npos) {
                            model.insertion_relevant = 1;
                            break;
                        }
                    }
                }
                if (model.insertion_relevant) {
                    model.insertion = focal_seq;
                }
            }
            candidate_model.push_back(model);
if (debug) {
    std::cerr << "VNTRAnnotator::find_repeat_region - Motif of HMM: " << display_ru << ',';
    for (auto b : tmp) {std::cerr << b;}
    std::cerr << '\t' << rr.score << '\t' << rr.match_motif << '\t' << rr.p_ed - rr.p_st + 1 << '\t' << rr.nr << '\t' << model.insertion_relevant << '\n';
}
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
            for (auto& w: pseq_obj.candidate) {
                if (context.find(w) == std::string::npos) {
                    continue;
                }
                candidate_unit tmp(w,0);
                tmp.reduce();
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
            candidate_ru.insert(candidate_unit(cseq));
            added_ct++;
        }
    } else {
        // if (cseq.at(0) == cseq.at(cseq.size()-1)) {
        //     cseq = cseq.substr(0, cseq.size()-1);
        // }
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
                    if (s.length() == 1) {
                        continue;
                    }
                    size_t pt = cseq.find(s);
                    if (pt == std::string::npos) {
                        continue;
                    }
                    candidate_unit tmp(s,1);
                    std::vector<std::map<int32_t, int32_t> > ovar(s.size());
                    while (pt != std::string::npos) {
                        for (size_t i = 0; i < s.size(); ++i) {
                            if (ovar[i].find(rl[pt+i]) == ovar[i].end()) {
                                ovar[i][rl[pt+i]] = 1;
                            } else {
                                ovar[i][rl[pt+i]] += 1;
                            }
                        }
                        pt = cseq.find(s, pt+s.size());
                    }
                    std::string ss;
                    bool flag = 1;
                    for (size_t i = 0; i < s.size(); ++i) {
                        auto it = std::prev(ovar[i].end());
                        if (it->second == 1) {
                            ovar[i].erase(it);
                        } else if (ovar[i].begin()->second == 1) {
                            ovar[i].erase(ovar[i].begin());
                        } // At most erase one extreme value
                        if (ovar[i].size() == 0) {
                            flag = 0;
                            break;
                        }
                        int32_t omin = ovar[i].begin()->first;
                        int32_t omax = std::prev(ovar[i].end())->first;
                        if (omin == omax) {
                            for (int32_t j = 0; j < omin; ++j) {
                                ss += s.at(i);
                                tmp.variable_base.push_back(std::make_pair(1,1));
                            }
                        } else {
                            ss += s.at(i);
                            tmp.variable_base.push_back(std::make_pair(omin, omax));
                        }
                    }
                    if (flag) {
                        tmp.ru = ss;
                        tmp.check();
                        candidate_ru.insert(tmp);
                        added_ct++;
                    }
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
 * Returns true if is to be classified as a VNTR
 */
 bool VNTRAnnotator::is_vntr(candidate_fuzzy_motif& motif, int32_t mode)
 {
     uint32_t mlen = motif.mlen;
     uint32_t rlen = motif.ed - motif.st + 1;
     if (motif.insertion_relevant) {
         rlen += motif.insertion.length();
     }
     if (mode==VITERBI_VNTR) {
         if (motif.viterbi_score < 3.) {
             return false;
         }
         if (mlen==1) {
             return motif.concordance >= 0.85;
         }
         return ((rlen - mlen) >= 6 && motif.n_ru>=2 && motif.concordance >= 0.75);
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
