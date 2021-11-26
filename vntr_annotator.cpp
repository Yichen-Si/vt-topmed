#include "vntr_annotator.h"

VNTRAnnotator::VNTRAnnotator(std::string& ref_fasta_file, bool _debug, uint32_t _m) : debug(_debug), max_mlen(_m) {
    fai = fai_load(ref_fasta_file.c_str());
    if (fai==NULL) {
        fprintf(stderr, "[%s:%d %s] Cannot load genome index: %s\n", __FILE__, __LINE__, __FUNCTION__, ref_fasta_file.c_str());
        exit(1);
    }
    min_ecover_indel = 0.5; min_pcover_indel = 0.5;
    min_ecover_extended = 0.3; min_pcover_extended = 0.3;
    extend_bp = 150;
};

VNTRAnnotator::~VNTRAnnotator() {
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
    while (it != candidate_model.end()) {
        bool dd = is_vntr(*it, mode);
        if (dd) {
            if (it->inexact && (!fz)) {
                fz = 1;
            }
            if ((!it->inexact) && (!ex)) {
                ex = 1;
            }
            it++;
        } else {
            it = candidate_model.erase(it);
        }
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
    int32_t seq_len;
    int32_t ct_indel = 0;
    const char* chrom = bcf_get_chrom(h, v);
    char** alleles = bcf_get_allele(v);
    int32_t indel_end = v->pos + strlen(alleles[0]); // 0-based first base after indel
    bool insertion = 0;
    std::string focal_seq(alleles[0]);
    if (strlen(alleles[1]) > focal_seq.size()) {
        insertion = 1;
        focal_seq = alleles[1];
    }
    focal_seq = focal_seq.substr(1); // only keep the inserted or deleted seq

    if (debug) {
        std::cerr << "============================================\n";
        std::cerr << "Read indel: " << pos1 << '\t' << alleles[0] << '\t' << alleles[1] << '\t';
        printf("%s\n", faidx_fetch_seq(fai, chrom, v->pos-20, v->pos+20, &seq_len));
    }

    // check if it is an indel inside homopolymer
    char b;
    int32_t homopoly = if_homopoly(chrom, v->pos, indel_end, b);
    if (homopoly >= HOMOPOLYMER_MIN) {
        std::string s(1, b);
        candidate_unit ru(s, 0);
        candidate_ru.insert(ru);
    }

    if (focal_seq.size() >= 6) {
        // ct_indel += og_find_repeat_unit(focal_seq, candidate_ru, 1);
        ct_indel += rl_find_repeat_unit(focal_seq, candidate_ru, 1);
    }

    // extend in both direction (small region)
    int32_t pad = std::min((int32_t)focal_seq.size() * 6, 32);
    if (focal_seq.size() < 2) {pad = 12;}
    std::string extended_seq = focal_seq;

    ShiftAlign align(fai, h, v, max_mlen, max_mlen*3, debug);
    int32_t st0, ed0, st1, ed1;
    int32_t extended_st0, extended_ed0;
    std::string mseq;
    align.right_shift(st0, st1, extended_ed0, extended_seq);
    if (extended_seq.length() > 3) {
        if (debug) {
            std::cerr << "Query after align " << extended_seq << '\n';
        }
        ct_indel += rl_find_repeat_unit(extended_seq, candidate_ru, 0);
    }
    if (align.delta_exact) {
        candidate_unit tmp(focal_seq);
        tmp.reduce();
        candidate_ru.insert(tmp);
    }
    align.left_shift(ed0, extended_st0, ed1, extended_seq);
    if (extended_seq.length() > 3) {
        if (debug) {
            std::cerr << "Query after align " << extended_seq << '\n';
        }
        ct_indel += rl_find_repeat_unit(extended_seq, candidate_ru, 0);
    }

    int32_t ct0 = ct_indel;
    if (extended_ed0 - indel_end  + 1 < pad) {
        extended_ed0 = indel_end + pad - 1;
        extended_seq = faidx_fetch_seq(fai, chrom, v->pos+1, extended_ed0, &seq_len);
        if (insertion && focal_seq.length() > 1) {
            extended_seq = extended_seq.insert(1, focal_seq);
            if (debug) {
                std::cerr << "Query right " << extended_seq << '\n';
            }
        }
        ct_indel += rl_find_repeat_unit(extended_seq, candidate_ru, 0);

    }
    if (v->pos - extended_st0 < pad) {
        extended_st0 = std::max(0, v->pos - pad);
        extended_seq = faidx_fetch_seq(fai, chrom, extended_st0, v->pos, &seq_len) + focal_seq;
        if (debug) {
            std::cerr << "Query left " << extended_seq << '\n';
        }
        ct_indel += rl_find_repeat_unit(extended_seq, candidate_ru, 0);
    }
    if (ct_indel == ct0) {
        pad = std::min(pad, 12);
        extended_seq = faidx_fetch_seq(fai, chrom, v->pos-pad, v->pos+pad, &seq_len);
        if (debug) {
            std::cerr << "Query joint " << extended_seq << '\n';
        }
        ct_indel += rl_find_repeat_unit(extended_seq, candidate_ru, 0);
    }
    if (debug && candidate_ru.size() > 0) {
        if (candidate_ru.size() > 0) {
            std::cerr << "Candidate RU:\t";
            for (auto & v : candidate_ru) {
                std::cerr << v.ru << '\t';
            }
            std::cerr << '\n';
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

        if (s.ru.length() == 1) {
            int32_t hst, hed, nr;
            char b = s.ru.at(0);
            hst = spiked_ed;
            find_homopoly_region(query, b, hst, hed, nr);
            if (hed < 0) {
                hst = spiked_st-1;
                find_homopoly_region(query, b, hst, hed, nr);
            }
            if (hed < 0) {
                continue;
            }
            int32_t real_st = start0 + hst;
            int32_t real_ed = start0 + hed;
            int32_t l = hed - hst + 1;
            if (is_insertion) {
                if (hed > spiked_ed) {
                    real_ed -= focal_seq.length();
                }
                if (hst > spiked_ed) {
                    real_st -= focal_seq.length();
                }
                if (hed < spiked_ed && hed >= spiked_st) {
                    real_ed = v->pos;
                }
                if (hst < spiked_ed && hst >= spiked_st) {
                    real_st = v->pos;
                }
            }
            candidate_fuzzy_motif model(s.ru, inexact_label, real_st, real_ed, nr, nr*2, l, nr);
            candidate_model.push_back(model);
            continue;
        }

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
        // wphmm->viterbi();
        std::string vpath = wphmm->print_viterbi_path();
        // wphmm->count_ru();
        wphmm->detect_range();
        if (wphmm->select_segment(spiked_st-1, spiked_ed+1)) {
            seq_segment& rr = wphmm->focal_rr;
        // for (auto &rr : wphmm->segments) {
        //     if (rr.match_motif == 0) {
        //         continue;
        //     }
            int32_t real_st = start0 + rr.p_st;
            int32_t real_ed = start0 + rr.p_ed;
            int32_t l = rr.p_ed - rr.p_st + 1;
            if (is_insertion) {
                if (rr.p_ed > spiked_ed) {
                    real_ed -= focal_seq.length();
                }
                if (rr.p_st > spiked_ed) {
                    real_st -= focal_seq.length();
                }
                if (rr.p_ed < spiked_ed && rr.p_ed >= spiked_st) {
                    real_ed = v->pos;
                }
                if (rr.p_st < spiked_ed && rr.p_st >= spiked_st) {
                    real_st = v->pos;
                }
            }
            candidate_fuzzy_motif model(display_ru, inexact_label, real_st, real_ed, rr.nr, rr.score, l, rr.match_motif);
            if (is_insertion && rr.p_ed >= spiked_st-1 && rr.p_st <= spiked_ed) { // Decide if the inserted sequence is relevant (need to be explained by final model) to the TR.
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
                if (nn>0.1*focal_seq.size()) {
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
            if (debug) {
                std::cerr << "VNTRAnnotator::find_repeat_region - Motif of HMM: " << display_ru << ',';
                for (auto b : tmp) {std::cerr << b;}
                std::cerr << '\t' << model.viterbi_score << '\t' << model.concordance << '\t' << rr.match_motif << '\t' << rr.p_ed - rr.p_st + 1 << '\t' << rr.nr << '\t' << model.insertion_relevant << '\t' << real_st << '-' << real_ed << '\n';
                // std::cerr << vpath.substr(0, extend_bp) << '\n';
                // std::cerr << vpath.substr(extend_bp) <<  '\n';
            }
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
        if (cseq.size() > 2) {

            periodic_seq pseq_obj(cseq.c_str(), max_mlen, debug, 1);
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
                                // tmp.mode.push_back(1);
                            }
                        } else {
                            // int32_t otot = 0;
                            // for (auto & v : ovar[i]) {
                            //     otot += v.second;
                            // }
                            // int32_t qt = 0, med1 = 0, med2 = 0;
                            // it = ovar[i].begin();
                            // med1 = it->first;
                            // while (qt < otot/2) {
                            //     qt += it->second;
                            //     if (qt >= otot / 2) {
                            //         break;
                            //     }
                            //     med1 = it->first;
                            //     it++;
                            // }
                            // med2 = it->first;
                            // if (med2 != med1) {
                            //     med1 = (int32_t) ((med1+med2)/2);
                            // }
                            ss += s.at(i);
                            tmp.variable_base.push_back(std::make_pair(omin, omax));
                            // tmp.mode.push_back(med1);
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


void VNTRAnnotator::find_homopoly_region(std::string& query, char b, int32_t& st, int32_t& ed, int32_t& nr) {

    int32_t max_insertion = 3;
    double max_cost = 0.15;
    int32_t n = query.length();
    ed = -1;
    nr = 0;
    if (st >= n || st < 0) {
        return;
    }
    // Find first segment
    int32_t l = st, r = st;
    while (l < n && query.at(l) != b) {
        l++;
    }
    if (l == n) {
        return;
    }
    r = l;
    while (r < n && query.at(r) == b) {
        r++;
    }
    while (l >=0 && query.at(l) == b) {
        l--;
    }
    l++;
    nr = r-l;
    if (nr < HOMOPOLYMER_MIN) {
        return;
    }
    // Extend to the left
    int32_t ct, cst, ced;
    while (l > 0) {
        ct = 0;
        cst = l - 1;
        while (cst >= 0 && query.at(cst) != b) {
            ct++;
            cst--;
        }
        if (ct > max_insertion || cst <= 0) {
            break;
        }
        ced = cst;
        while (ced >= 0 && query.at(ced) == b) {
            ced--;
        }
        if (cst - ced > ct/max_cost) {
            nr += cst-ced;
            l = ced + 1;
        } else {
            break;
        }
    }
    // Extend to the right
    while (r < n) {
        ct = 0;
        cst = r;
        while (cst < n && query.at(cst) != b) {
            ct++;
            cst++;
        }
        if (ct > max_insertion || cst >= n) {
            break;
        }
        ced = cst;
        while (ced < n && query.at(ced) == b) {
            ced++;
        }
        if (ced - cst > ct/max_cost) {
            nr += ced-cst;
            r = ced;
        } else {
            break;
        }
    }
    st = l;
    ed = r-1;
    return;
}

bool VNTRAnnotator::choose_models(VNTR_candidate& vntr, int32_t mode) {
        vntr.combine_insertions();
        if (!vntr.choose_model()) { // Choose model based on the longest allele
            return 0;
        }
        // Keep at most two candidate RUs
        // (if there exists reasonable exact and inexact ones)
        // TODO: this is temporary
        auto it = vntr.alternative_models.begin();
        bool ex = 0, fz = 0;
        while (it != vntr.alternative_models.end()) {
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
            it = vntr.alternative_models.erase(it);
        }
        return 1;
}


/**
 * Returns true if is to be classified as a VNTR
 */
 bool VNTRAnnotator::is_vntr(candidate_fuzzy_motif& motif, int32_t mode)
 {
     uint32_t mlen = motif.ru.length();
     uint32_t rlen = motif.ed - motif.st + 1;
     if (motif.insertion_relevant) {
         rlen += motif.insertion.length();
     }
     if (mode==VITERBI_VNTR) {
         if (motif.viterbi_score < 3.) {
             return false;
         }
         if (mlen==1) {
             return motif.concordance >= 0.85 || (motif.concordance >= 0.8 && motif.n_ru > 10 && motif.viterbi_score > 10);
         }
         if (motif.n_ru < 2 || motif.concordance < 0.75 || rlen - mlen < 6) {
             return false;
         }
         if (motif.inexact) {
             if (motif.n_ru < 3) {
                 return false;
             }
             int32_t ct = 0;
             for (auto v : motif.label) {
                 ct += v;
             }
             if (ct == mlen) {
                 if (motif.n_ru >= 5) {
                     return motif.concordance >= 0.9;
                 }
                 return (motif.concordance >= 0.95);
             }
             return (motif.concordance >= 0.85);
         }
         return true;
         // return ((rlen - mlen) >= 6 && motif.n_ru>=2 && motif.concordance >= 0.75);
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
