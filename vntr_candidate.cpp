#include "vntr_candidate.h"

VNTR_candidate::VNTR_candidate(faidx_t* _f, bcf_hdr_t* h, bcf1_t* v, candidate_fuzzy_motif& _m, bool _b) : debug(_b), fai(_f), rid(v->rid), motif(_m), st(_m.st), ed(_m.ed), mlen(_m.mlen) {
    rlen_l = motif.l;
    n_ru_l = motif.n_ru;
    score_l = motif.viterbi_score;
    if (motif.insertion_relevant || ed - st <= mlen) {
        rlen_r = 0;
        n_ru_r = 0;
        score_r = 0;
    } else {
        rlen_r = rlen_l;
        n_ru_r = n_ru_l;
        score_r = score_l;
    }
    chrom = h->id[BCF_DT_CTG][rid].key;
    need_refit = 0;
    critical_ovlp = 0.5;
    max_interrupt = 0.1;
    if (motif.insertion_relevant) {
        add_insertion(v->pos, motif.insertion);
    }
    alternative_models.push_back(motif);
    combine_insertions();
}

void VNTR_candidate::add_insertion(int32_t _p, std::string _s) {
    auto it = insertions.begin();
    while (it != insertions.end()) {
        if (it->first == _p) {
            if (it->second.length() < _s.length()) {
                it->second = _s;
            }
            return;
        }
        if (it->first > _p) {
            break;
        }
        ++it;
    }
    insertions.insert(it, std::make_pair(_p,_s));
}

void VNTR_candidate::combine_insertions(int32_t pad) {

    int32_t seq_len;
    char* seq = faidx_fetch_seq(fai, chrom, st, ed, &seq_len);
    repeat_ref.assign(seq);
    free(seq);
    // repeat_ref = faidx_fetch_seq(fai, chrom, st, ed, &seq_len);
    std::string merged_longest_rr = repeat_ref;
    // std::sort(insertions.begin(), insertions.end());
    int32_t offset = 1;
    for (auto &v : insertions) {
        if (v.first == st-1) {
            merged_longest_rr.insert(0, v.second);
            offset += v.second.length();
        } else if (v.first <= ed && v.first >= st) {
            merged_longest_rr.insert(v.first-st+offset, v.second);
            offset += v.second.length();
        }
    }
    ins_tot = offset - 1;
    if (pad < 0) {
        pad = 10+repeat_ref.size()/10;
    }
    seq = faidx_fetch_seq(fai, chrom, std::max(st-pad, 0), st-1, &seq_len);
    query.assign(seq);
    free(seq);
    // query = faidx_fetch_seq(fai, chrom, std::max(st-pad, 0), st-1, &seq_len);
    rel_st = query.length();
    query += merged_longest_rr;
    seq = faidx_fetch_seq(fai, chrom, ed+1, ed+pad, &seq_len);
    query.append(seq);
    free(seq);
    // query += faidx_fetch_seq(fai, chrom, ed+1, ed+pad, &seq_len);
    len_mrg = merged_longest_rr.length();
}

int32_t VNTR_candidate::get_pos_in_ref(int32_t p_rel, bool right) {

    if (len_mrg <= repeat_ref.length()) {
        return st + p_rel - rel_st;
    }
    int32_t l_rel = 0, r_rel = 0;
    int32_t offset = 0;
    for (auto &v : insertions) {
        if (v.first <= ed && v.first >= st - 1) {
            l_rel = v.first - st + 1 + rel_st + offset;
            r_rel = l_rel + v.second.length() - 1;
            if (p_rel < l_rel) {
                return st + p_rel - rel_st - offset;
            }
            if (p_rel >= l_rel && p_rel <= r_rel) {
                if (right) {
                    return v.first + 1;
                } else {
                    return v.first;
                }
            }
            offset += v.second.length();
        }
    }
    return st + p_rel - rel_st - offset;
}

int32_t VNTR_candidate::merge(VNTR_candidate& rt) {
    if (need_refit) {
        model_refit();
    }
    if (rt.need_refit) {
        rt.model_refit();
    }
    double cover = 0.85;
    candidate_fuzzy_motif& rhs = rt.motif;
    int32_t cp = motif.ru_compare(rhs);
    int32_t ovlp = intersect(rt);
    int32_t min_l = std::min(ed - st + 1, rt.ed - rt.st + 1);
    double f_ovlp = ((double) ovlp) / min_l;

    // Rule out immediate incompatible
    if (((cp < 0 || cp == 2) && f_ovlp < critical_ovlp) || -f_ovlp > max_interrupt) {
        return 0;
    }
    if (ovlp == 1 && (rt.ed == rt.st || ed == st)) { // Repeat region is purely inserted
        if (cp != 0) {
            return 0;
        }
        return 1;
    }

if (debug) {
    printf( "--- VNTR_candidate::merge - %d, %d, %s_%d with query %d, %d, %s_%d. %d, %.2f --- \n", st,ed,motif.ru.c_str(),motif.inexact, rt.st,rt.ed,rt.motif.ru.c_str(),rt.motif.inexact,cp,f_ovlp );
}
    int32_t u_st = std::min(st, rhs.st);
    int32_t u_ed = std::max(ed, rhs.ed);
    bool novel = 1;
    for (auto &m : alternative_models) {
        int32_t cmp = m.ru_compare(rhs);
        if (cmp == 0 || cmp == 3) {
            novel = 0;
            break;
        }
    }

    // If readily mergeable
    if (cp == 0 || (cp > 0 && f_ovlp >= critical_ovlp) ) {
        need_refit = (st - u_st > motif.mlen || u_ed - ed > motif.mlen);
        st = u_st;
        ed = u_ed;
        for (auto &v : rt.insertions) {
            add_insertion(v.first, v.second);
        }
        if (novel) {
            alternative_models.push_back(rhs);
            candidate_fuzzy_motif& mot = alternative_models.back();
            evaluate_models(mot);
            if (mot.viterbi_score > 3 && mot.l > rlen_l * 0.6 && mot.n_ru >= 2 && mot.concordance > 0.75) {
                std::sort(alternative_models.begin(), alternative_models.end());
            } else {
                alternative_models.pop_back();
            }
            return 1;
        }
        return 1;
    }

    // O.w. (incompatible RU or not enough overlap)
    // Need to check if a model fits both segments well
    // Check if the current model fits the query sequence well
    if (!rt.fit_alt_model(motif)) { // Does not explain the additional region well
        return 0;
    }
    if (st > u_st || ed < u_ed) { // Extend current focal region
        need_refit = (st - u_st > motif.mlen || u_ed - ed > motif.mlen);
        st = u_st;
        ed = u_ed;
    }
    // Check if the query model is worth considering
    if (novel) {
        alternative_models.push_back(rhs);
        candidate_fuzzy_motif& mot = alternative_models.back();
        evaluate_models(mot);
        if (mot.viterbi_score > 3 && mot.l > rlen_l * 0.6 && mot.n_ru >= 2 && mot.concordance > 0.75) {
            std::sort(alternative_models.begin(), alternative_models.end());
        } else {
            alternative_models.pop_back();
        }
    }
    if (f_ovlp >= critical_ovlp) { // Close enough, and can be explained, merge
        return 1;
    }
    return 2;
}

int32_t VNTR_candidate::intersect(VNTR_candidate& rt) {
    // Non-overlapping
    if (ed < rt.st) {
        return ed+1 - rt.st;
    }
    if (st > rt.ed) {
        return rt.ed+1 - st;
    }
    // Overlapping
    int32_t u_st = std::min(st, rt.st);
    int32_t u_ed = std::max(ed, rt.ed);
    return  ed-st+1+rt.ed-rt.st+1- (u_ed-u_st+1);
}

bool VNTR_candidate::fit_alt_model(candidate_fuzzy_motif& rhs) {
    if (need_refit) {
        model_refit();
    }
    if ( std::abs(rhs.st - st) <= mlen && std::abs(rhs.ed - ed) <= mlen ) {
        return (rhs.viterbi_score > 3 && rhs.l > rlen_l * 0.6 && rhs.n_ru >= 2 && rhs.concordance > 0.75);
    }
    std::string unit  = rhs.ru;
    bool bv[rhs.mlen];
    std::copy(rhs.label.begin(), rhs.label.end(), bv);
    WPHMM_UNGAP* wphmm = new WPHMM_UNGAP(query.c_str(), rhs.ru.c_str(), debug, bv, rhs.inexact);
    wphmm->set_ru(rhs.ru, rhs.label);
    wphmm->initialize();
    wphmm->viterbi();
    wphmm->detect_range();
    bool decision = 0;
    if (wphmm->segments.size() == 0) {
        decision = false;
if (debug) {
    std::cerr << "VNTR_candidate::fit_alt_model Rejected\n";
    std::cerr << "Query model " << rhs.ru << "\tOriginal model " << motif.ru << '\n';
}
    } else {
        seq_segment& rr = wphmm->focal_rr;
        decision = (rr.score > 3 && rr.l > rlen_l * 0.6 && rr.nr >= 2 && rr.match_motif > 0.75*rr.l);
if (debug) {
    std::cerr << "VNTR_candidate::fit_alt_model " << decision << '\t';
    std::cerr << "Query model " << rhs.ru << "\tOriginal model " << motif.ru << '\n';
    std::cerr << "New length, nr, log2 OR: " << rr.l << '\t' << rr.nr << '\t' << rr.score << '\t';
    std::cerr << "Original: " << rlen_l << '\t' << n_ru_l << '\t' << score_l << '\n';
}
    }
    delete wphmm;
    wphmm = NULL;
    return decision;
}

void VNTR_candidate::model_refit() {
    if (alternative_models.size() == 0) {
        error("VNTR_candidate::model_refit No candidate model");
    }
    combine_insertions();
    evaluate_models();
    score_l = motif.viterbi_score;
    n_ru_l = motif.n_ru;
    rlen_l = motif.l;
    // Take union of repeat region by alternative models
    int32_t st1 = motif.st;
    int32_t ed1 = motif.ed;
    for (auto & v : alternative_models) {
        if (v.st < st1) {
            st1 = v.st;
        }
        if (v.ed > ed1) {
            ed1 = v.ed;
        }
    }
    if (st1 != st || ed1 != ed) {
        st = st1;
        ed = ed1;
        combine_insertions();
    }
    need_refit = 0;
}

bool VNTR_candidate::finalize(int32_t flk) {

    if (alternative_models.size() == 0) {
        return 0;
    }
    int32_t seq_len;
    char* seq = faidx_fetch_seq(fai, chrom, st, ed, &seq_len);
    repeat_ref.assign(seq);
    free(seq);
    if (len_mrg > repeat_ref.length() && repeat_ref.length() >= mlen) {
        int32_t n_model = 0;
        for (auto & mot : alternative_models) {
            if (mot.mlen == 1) {
                n_model++;
                continue;
            }
            bool bv[mot.mlen];
            std::copy(mot.label.begin(), mot.label.end(), bv);
            WPHMM_UNGAP* wphmm = new WPHMM_UNGAP(repeat_ref.c_str(), mot.ru.c_str(), debug, bv, mot.inexact);
            wphmm->set_ru(mot.ru, mot.label);
            wphmm->initialize();
            wphmm->viterbi();
            wphmm->detect_range();
            if (n_model == 0) {
                if (wphmm->segments.size()>0) {
                    seq_segment& rr = wphmm->focal_rr;
                    score_r       = rr.score;
                    n_ru_r        = rr.nr;
                    rlen_r        = rr.p_ed - rr.p_st + 1;
                    concordance_r = (double) rr.match_motif / rlen_r;
                } else {
                    score_r       = wphmm->viterbi_score;
                    n_ru_r        = 0;
                    rlen_r        = 0;
                    concordance_r = 0.;
                }
            } else {
                if (wphmm->segments.size()>0) {
                    seq_segment& rr = wphmm->focal_rr;
                    score_r2       = rr.score;
                    n_ru_r2        = rr.nr;
                    rlen_r2        = rr.p_ed - rr.p_st + 1;
                    concordance_r2 = (double) rr.match_motif / rlen_r2;
                } else {
                    score_r2       = wphmm->viterbi_score;
                    n_ru_r2        = 0;
                    rlen_r2        = 0;
                    concordance_r2 = 0.;
                }
            }
            n_model++;
            delete wphmm;
            wphmm = NULL;
if (debug) {
    printf("VNTR_candidate::finalize Fit reference only model. length, nr, score: %d, %d, %.3f\n",rlen_r, n_ru_r, score_r);
}
        }
    } else {
        score_r = score_l;
        n_ru_r = n_ru_l;
        rlen_r = rlen_l;
        concordance_r = motif.concordance;
        if (alternative_models.size() > 1) {
            score_r2 = alternative_models[1].viterbi_score;
            n_ru_r2 = alternative_models[1].n_ru;
            rlen_r2 = alternative_models[1].l;
            concordance_r2 = alternative_models[1].concordance;
        }
    }
    seq = faidx_fetch_seq(fai, chrom, st-flk, st-1, &seq_len);
    lflank.assign(seq);
    free(seq);
    seq = faidx_fetch_seq(fai, chrom, ed+1, ed+flk, &seq_len);
    rflank.assign(seq);
    free(seq);

    return 1;
}

bool VNTR_candidate::choose_model() {
    if (alternative_models.size() == 0) {
        return 0;
    }
    if (alternative_models.size() > 1 || need_refit) {
        evaluate_models();
    } else {
        motif = alternative_models[0];
    }
    if (st != motif.st || ed != motif.ed) {
        st = motif.st;
        ed = motif.ed;
    }
    score_l = motif.viterbi_score;
    n_ru_l = motif.n_ru;
    rlen_l = motif.l;
    return 1;
}

void VNTR_candidate::evaluate_models() {
    if (alternative_models.size() == 0) {
        error("VNTR_candidate: No candidate model");
    }
    for (uint32_t i = 0; i < alternative_models.size(); ++i) {
        evaluate_models(alternative_models[i]);
    }
    std::sort(alternative_models.begin(), alternative_models.end());
    motif = alternative_models[0];
    auto it = alternative_models.begin();
    it++;
    while (it != alternative_models.end()) {
        if (it->ru.length() == 1 && it->n_ru < 6) {
            it = alternative_models.erase(it);
            continue;
        }
        if (it->n_ru <= 1 || it->concordance <= 0.75 || it->viterbi_score < 3) {
            it = alternative_models.erase(it);
            continue;
        }
        it++;
    }
}


void VNTR_candidate::evaluate_models(candidate_fuzzy_motif& mot) {
    bool bv[mot.mlen];
    std::copy(mot.label.begin(), mot.label.end(), bv);
    // Run HMM - with the longest allele
    WPHMM_UNGAP* wphmm = new WPHMM_UNGAP(query.c_str(), mot.ru.c_str(), debug, bv, mot.inexact);
    wphmm->set_ru(mot.ru, mot.label);
    wphmm->initialize();
    wphmm->viterbi();
    // std::string v_path = wphmm->print_viterbi_path();
    wphmm->detect_range();
    if (wphmm->segments.size() == 0) { // Rare case
        mot.st = st;
        mot.ed = ed;
        mot.l  = len_mrg;
        mot.n_ru = 0;
        mot.viterbi_score = 0;
        mot.concordance = 0;
    } else {
        seq_segment& rr = wphmm->focal_rr;
        mot.st = get_pos_in_ref(rr.p_st, 1);
        mot.ed = get_pos_in_ref(rr.p_ed, 0);
        mot.l  = rr.p_ed - rr.p_st + 1;
        mot.n_ru = rr.nr; // n_ru in longest mosaic sequence
        mot.viterbi_score = rr.score;
        mot.concordance = (double) rr.match_motif / mot.l;
    }
    if (debug) {
        std::cerr << "VNTR_candidate::evaluate_models " << mot.ru << '\t' << mot.inexact << '\t' << mot.st << ',' << mot.ed << '\t' << st << ',' << ed  << '\t' << mot.l << '\t' << mot.n_ru << '\t' << mot.concordance << '\t' << mot.viterbi_score << '\n';
    }
    delete wphmm;
    wphmm = NULL;
}


/**
 * Print object.
 */
void VNTR_candidate::print()
{
    std::cerr << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n";
    std::cerr << "VNTR Summary\n";
    printf("ru\t%s\n", motif.ru.c_str());
    printf("repeat region (ref)\t%s\n", repeat_ref.c_str());
    printf("repeat region position\t%d,%d\n", st, ed);
    printf("repeat region length\t%d,%d\n", rlen_l, rlen_r);
    printf("number of complete repeat units\t%d,%d\n", n_ru_l, n_ru_r);
    printf("viterbi score\t%.3f,%.3f\n", score_l, score_r);
};
