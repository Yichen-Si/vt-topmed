#include "vntr_candidate.h"

VNTR_candidate::VNTR_candidate(faidx_t* _f, bcf_hdr_t* h, bcf1_t* v, candidate_fuzzy_motif& _m, bool _b) : debug(_b), fai(_f), rid(v->rid), motif(_m), st(_m.st), ed(_m.ed), mlen(_m.mlen) {
    rlen_l = motif.l;
    rlen_r = motif.l;
    n_ru_r = motif.n_ru;
    n_ru_l = motif.n_ru;
    score_l = motif.viterbi_score;
    score_r = score_l;
    chrom = h->id[BCF_DT_CTG][rid].key;
    need_refit = 1;
    critical_ovlp = 0.5;
    max_interrupt = 0.1;
    if (motif.insertion_relevant) {
        add_insertion(v->pos, motif.insertion);
    }
    alternative_models.push_back(motif);
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
            insertions.insert(it, std::make_pair(_p,_s));
            return;
        }
        ++it;
    }
    insertions.push_back(std::make_pair(_p,_s));
}

void VNTR_candidate::combine_insertions() {
    int32_t seq_len;
    repeat_ref = faidx_fetch_seq(fai, chrom, st, ed, &seq_len);
    merged_longest_rr = repeat_ref;
    // std::sort(insertions.begin(), insertions.end());
    int32_t offset = 0;
    for (auto &v : insertions) {
        if (v.first <= motif.ed && v.first >= motif.st) {
            merged_longest_rr.insert(v.first-motif.st+offset, v.second);
            offset += v.second.length();
        }
    }
}

int32_t VNTR_candidate::get_pos_in_ref(int32_t p_rel) {
    if (merged_longest_rr.length() <= repeat_ref.length()) {
        return motif.st + p_rel;
    }
    int32_t l_rel = 0, r_rel = 0;
    int32_t offset = 0;
    int32_t pos = motif.st + p_rel;
    for (auto &v : insertions) {
        if (v.first <= motif.ed && v.first >= motif.st) {
            l_rel = v.first - motif.st + offset;
            r_rel = l_rel + v.second.length();
            if (p_rel < l_rel) {
                return motif.st + p_rel - offset;
            }
            if (p_rel >= l_rel && p_rel <= r_rel) {
                return v.first;
            }
            offset += v.second.length();
        }
    }
    return motif.st + p_rel - offset;
}

int32_t VNTR_candidate::merge(VNTR_candidate& rt) {
    if (need_refit) {
        model_refit();
    }
    if (rt.need_refit) {
        rt.model_refit();
    }
    candidate_fuzzy_motif& rhs = rt.motif;
    int32_t cp = ru_compare(rhs);
    int32_t ovlp = intersect(rt);
    int32_t min_l = std::min(rlen_l, rhs.l);
    double f_ovlp = ((double) ovlp) / min_l;

printf( "--- VNTR_candidate::merge - %d, %d, %s_%d with query %d, %d, %s_%d. %d,%.2f --- \n", st,ed,motif.ru.c_str(),motif.inexact, rt.st,rt.ed,rt.motif.ru.c_str(),rt.motif.inexact,cp,f_ovlp );

    // Rule out immediate incompatible
    if ((cp == 3 && f_ovlp < critical_ovlp) || -f_ovlp > max_interrupt) {
        return 0;
    }
    if (cp == 1 || cp == 2) { // Don't compare scores now, both mean compatible
        cp = 1;
    }
    int32_t u_st = std::min(st, rhs.st);
    int32_t u_ed = std::max(ed, rhs.ed);
    // Rule out immediate mergeable (Type 1 merge)
    if (cp == 0 || (cp == 1 && f_ovlp >= critical_ovlp)) {
        need_refit = (st > u_st || ed < u_ed);
        st = u_st;
        ed = u_ed;
        for (auto &v : rt.insertions) {
            add_insertion(v.first, v.second);
        }
        if (cp == 1) {
            alternative_models.push_back(rhs);
            if (!need_refit) {
                evaluate_models(alternative_models.size()-1);
            }
        }
        return 1;
    }
    // Need to check if a model fits both segments well
    double  cover = 0.85;
    if (cp==1 && rt.st >= st && rt.ed <= ed) {
        // query is nested in existing record
        for (auto &m : alternative_models) {
            int32_t cmp = m.ru_compare(rhs);
            if (cmp != 3) {return 1;} // Already in candidates
        }
    }
    if ((cp == 3 && f_ovlp >= critical_ovlp) || (cp==1 && f_ovlp < critical_ovlp)) {
        double m2r1 = fit_alt_model(rhs);
        double m1r2 = rt.fit_alt_model(motif);
        if (m2r1 >= score_l*cover || m1r2 >= rt.score_l*cover) { // At least one model can explain both
            alternative_models.push_back(rhs);
            if (m1r2 < rt.score_l*cover) { //
                return 0;
            }
            st = u_st;
            ed = u_ed;
            need_refit = 1;
            if (f_ovlp >= critical_ovlp) {
                return 1;
            }
            return 2;
        }
        return 0;
    }
if (debug) {
    error("VNTR_candidate::merge Uncovered situation. Trying to merge st1 %d, ed1 %d, st2 %d, ed2 %d, ru1 %s, ru2 %s\n", motif.st, motif.ed, rhs.st, rhs.ed, motif.ru.c_str(), rhs.ru.c_str());
}
        return 0;
    }

int32_t VNTR_candidate::intersect(VNTR_candidate& rt) {
    candidate_fuzzy_motif& rhs = rt.motif;
    // Non-overlapping
    if (ed < rhs.st) {
        return ed+1 - rhs.st;
    }
    if (st > rhs.ed) {
        return rhs.ed+1 - st;
    }
    // Overlapping
    int32_t u_st = std::min(st, rhs.st);
    int32_t u_ed = std::max(ed, rhs.ed);
    return  rlen_l+rhs.l - (u_ed-u_st+1);
}

int32_t VNTR_candidate::ru_compare(candidate_fuzzy_motif& rhs) {
    if (motif.ru.length() != rhs.ru.length()) {
        return 3; // Repeat unit incompatible
    }
    if (motif.ru == rhs.ru && motif.label == rhs.label) {
        return 0; // Equivalent
    }
    bool equal_seq = (motif.ru == rhs.ru);
    uint32_t phase = 1;
    while (phase < mlen && (!equal_seq)) {
        std::string cru = motif.ru.substr(phase, mlen-phase)+motif.ru.substr(0,phase);
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
        if (motif.label[(i+phase)%mlen] != rhs.label[i]) {
            equal_cir = 0;
            break;
        }
    }
    if (equal_cir) {
        return 0; // Circular equivalent
    }
    // Compatible, pick the model with higher score
    if (motif.viterbi_score >= rhs.viterbi_score) {
        return 1; // left unit has higher score
    } else {
        return 2; //
    }
}


double VNTR_candidate::fit_alt_model(candidate_fuzzy_motif& rhs) {
    if (need_refit) {
        model_refit();
    }
    std::string query = merged_longest_rr;
    std::string unit  = rhs.ru;
    bool bv[rhs.mlen];
    std::copy(rhs.label.begin(), rhs.label.end(), bv);
    WPHMM_UNGAP* wphmm = new WPHMM_UNGAP(query.c_str(), rhs.ru.c_str(), debug, bv, rhs.inexact);
    wphmm->set_ru(rhs.ru, rhs.label);
    wphmm->initialize();
    std::string v_path = wphmm->print_viterbi_path();
    wphmm->detect_range();
    if (debug) {
        std::cerr << "VNTR_candidate::fit_alt_model\n";
        std::cerr << query << std::endl;
        std::cerr << "Query model " << rhs.ru << "\tOriginal model " << motif.ru << '\n';
        // std::cerr << "New viterbi path\n" << v_path << std::endl;
        std::cerr << "WPHMM_UNGAP New log2 OR: " << wphmm->viterbi_score << '\t';
        std::cerr << "Original OR: " << score_l << '\n';
    }
    double new_score = -1;
    if (wphmm->segments.size() > 0) {
        new_score = wphmm->focal_rr.score;
    }
    delete wphmm;
    wphmm = NULL;
    return new_score;
}



void VNTR_candidate::model_refit() {
    if (alternative_models.size() == 0) {
        error("VNTR_candidate::model_refit No candidate model");
    }
    combine_insertions();
    evaluate_models();
    score_l = motif.viterbi_score;
    n_ru_l = motif.n_ru;
    rlen_l = ed - st + 1;
    need_refit = 0;
}

bool VNTR_candidate::finalize(int32_t flk) {
std::cerr << "VNTR_candidate::finalize\n";

    combine_insertions();
    if (!choose_model()) { // Choose model based on the longest allele
        return 0;
    }
std::cerr << "VNTR_candidate::finalize\t" << st << '\t' << ed << '\t' << motif.ru << '\n';
    if (merged_longest_rr.length() > repeat_ref.length() && repeat_ref.length() >= mlen) {
        std::string query = repeat_ref;
        bool bv[motif.mlen];
        std::copy(motif.label.begin(), motif.label.end(), bv);
        WPHMM_UNGAP* wphmm = new WPHMM_UNGAP(query.c_str(), motif.ru.c_str(), debug, bv, motif.inexact);
        wphmm->set_ru(motif.ru, motif.label);
        wphmm->initialize();
        wphmm->viterbi();
        wphmm->detect_range();
        if (wphmm->segments.size()>0) {
            seq_segment& rr = wphmm->focal_rr;
            score_r = rr.score;
            n_ru_r  = rr.nr;
            rlen_r  = rr.p_ed - rr.p_st + 1;
        } else {
            score_r = wphmm->viterbi_score;
            n_ru_r  = 0;
            rlen_r  = 0;
        }
        delete wphmm;
        wphmm = NULL;
if (debug) {
    printf("VNTR_candidate::finalize Fit reference only model. length, nr, score: %d, %d, %.3f\n",rlen_r, n_ru_r, score_r);
}
    }
    int32_t seq_len;
    lflank = faidx_fetch_seq(fai, chrom, st-flk, st-1, &seq_len);
    rflank = faidx_fetch_seq(fai, chrom, ed+1, ed+flk, &seq_len);
    return 1;
}


bool VNTR_candidate::choose_model() {
    if (alternative_models.size() == 0) {
        return 0;
    }
    if (alternative_models.size() > 1 || need_refit) {
        evaluate_models();
    }
    st = motif.st;
    ed = motif.ed;
    score_l = motif.viterbi_score;
    n_ru_l = motif.n_ru;
    rlen_l = ed - st + 1;
    return 1;
}

void VNTR_candidate::evaluate_models() {
    std::string query = merged_longest_rr;
    auto it = alternative_models.begin();
    while (it != alternative_models.end()) {
        motif = *it;
        bool bv[motif.mlen];
        std::copy(motif.label.begin(), motif.label.end(), bv);
        // Run HMM - with the longest allele
        WPHMM_UNGAP* wphmm = new WPHMM_UNGAP(query.c_str(), motif.ru.c_str(), debug, bv, motif.inexact);
        wphmm->set_ru(motif.ru, motif.label);
        wphmm->initialize();
        std::string v_path = wphmm->print_viterbi_path();
        wphmm->detect_range();
        if (wphmm->segments.size() == 0) { // Rare case
            motif.st = st;
            motif.ed = ed;
            motif.n_ru = wphmm->ru_complete.back(); // n_ru in longest mosaic sequence
            motif.viterbi_score = wphmm->viterbi_score;
        } else {
            seq_segment& rr = wphmm->focal_rr;
            motif.st = get_pos_in_ref(rr.p_st);
            motif.ed = get_pos_in_ref(rr.p_ed);
            motif.n_ru = rr.nr; // n_ru in longest mosaic sequence
            motif.viterbi_score = rr.score;
        }
        delete wphmm;
        wphmm = NULL;
        ++it;
    }
    std::sort(alternative_models.begin(), alternative_models.end());
    motif = alternative_models[0];
}


void VNTR_candidate::evaluate_models(uint32_t index) {
    if (index >= alternative_models.size()) {return;}
    std::string query = merged_longest_rr;
    candidate_fuzzy_motif& mot = alternative_models[index];
    bool bv[mot.mlen];
    std::copy(mot.label.begin(), mot.label.end(), bv);
    // Run HMM - with the longest allele
    WPHMM_UNGAP* wphmm = new WPHMM_UNGAP(query.c_str(), mot.ru.c_str(), debug, bv, mot.inexact);
    wphmm->set_ru(mot.ru, mot.label);
    wphmm->initialize();
    std::string v_path = wphmm->print_viterbi_path();
    wphmm->detect_range();
    if (wphmm->segments.size() == 0) { // Rare case
        mot.st = st;
        mot.ed = ed;
        mot.n_ru = wphmm->ru_complete.back(); // n_ru in longest mosaic sequence
        mot.viterbi_score = wphmm->viterbi_score;
    } else {
        seq_segment& rr = wphmm->focal_rr;
        mot.st = get_pos_in_ref(rr.p_st);
        mot.ed = get_pos_in_ref(rr.p_ed);
        mot.n_ru = rr.nr; // n_ru in longest mosaic sequence
        mot.viterbi_score = rr.score;
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
