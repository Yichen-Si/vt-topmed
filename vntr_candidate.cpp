/* The MIT License

   Copyright (c) 2015 Adrian Tan <atks@umich.edu>

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

#include "vntr_candidate.h"

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

void VNTR_candidate::combined_insertions(faidx_t* fai) {
    int32_t seq_len;
    repeat_ref = faidx_fetch_seq(fai, chrom, motif.st, motif.ed, &seq_len);
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


int32_t VNTR_candidate::merge(VNTR_candidate& rt, int32_t min_overlap) {
    candidate_fuzzy_motif& rhs = rt.motif;
    int32_t u_st = std::min(motif.st, rhs.st);
    int32_t u_ed = std::max(motif.ed, rhs.ed);
if (debug) {
    printf("VNTR_candidate::merge st1 %d, ed1 %d, st2 %d, ed2 %d, ru1 %s, ru2 %s\n", motif.st, motif.ed, rhs.st, rhs.ed, motif.ru.c_str(), rhs.ru.c_str());
}
    if ( motif.l+rhs.l - (u_ed-u_st+1) < min_overlap ) {
printf("Non-overlapping\n");

        return 0;
    }
    int32_t cp = ru_compair(rhs);
printf("ru_compair: %d\n",cp);
    if (cp == 3) {
        return 0;
    }
    if (cp == 2) {
        motif.ru = rhs.ru;
        motif.inexact = rhs.inexact;
        motif.label = rhs.label;
        updated = 1;
    }
    if (motif.st > u_st || motif.ed < u_ed) {updated = 1;}
    motif.st = u_st;
    motif.ed = u_ed;
    motif.l  = u_ed-u_st+1;
    for (auto &v : rt.insertions) {
        add_insertion(v.first, v.second);
    }
if (debug) {
    printf("VNTR_candidate::merge successful\n");
}
    return 1;
}

int32_t VNTR_candidate::intersect(VNTR_candidate& rt, int32_t min_overlap) {
    candidate_fuzzy_motif& rhs = rt.motif;
    if (motif.ed < rhs.st + min_overlap) {
        return -1;
    }
    if (motif.st > rhs.ed - min_overlap) {
        return 1;
    }
    int32_t u_st = std::min(motif.st, rhs.st);
    int32_t u_ed = std::max(motif.st, rhs.st);
    if ( motif.l+rhs.l - (u_ed-u_st+1) >= min_overlap ) {
        return 0;
    }
    if (motif.st == rhs.st) {
        return ((motif.ed < rhs.ed) ? -1 : 1);
    }
    return ((motif.st < rhs.st) ? -1 : 1);
}

int32_t VNTR_candidate::ru_compair(candidate_fuzzy_motif& rhs) {
    if (motif.ru.length() != rhs.ru.length()) {
        return 3; // Repeat unit incompatible
    }
    if (motif.ru == rhs.ru && motif.label == rhs.label) {
        return 0; // Equivalent
    }
    // Circular equivalent
    bool equal_seq = (motif.ru == rhs.ru);
    uint32_t phase = 1;
    while (phase < mlen && (!equal_seq)) {
        std::string cru = motif.ru.substr(phase, mlen-phase)+motif.ru.substr(0,phase);
        if (cru == rhs.ru) {
            equal_seq = 1;
        }
        phase++;
    }
    if (equal_seq) {
        // Compatible, pick the model with higher score
        if (motif.viterbi_score >= rhs.viterbi_score) {
            return 1; // Keep left unit
        } else {
            return 2; // Keep right unit
        }
    }
    return 3; // Repeat unit incompatible
}

/**
 * Fit an ungapped local alignment if updates have been made
 */
bool VNTR_candidate::model_refit(faidx_t* fai, int32_t flk) {
    if (!updated) { // Don't do anything
        int32_t seq_len;
        lflank = faidx_fetch_seq(fai, chrom, motif.st-flk, motif.st-1, &seq_len);
        rflank = faidx_fetch_seq(fai, chrom, motif.ed+1, motif.ed+flk, &seq_len);
        return 1;
    }
    if (merged_longest_rr.length() == 0) {
        combined_insertions(fai);
    }
    std::string query = merged_longest_rr;
    std::string unit  = motif.ru;
    std::vector<bool> tmp;
    uint32_t m = mlen;
    std::string display_ru = motif.ru;
    std::vector<bool> inexact_label(m, 0);

    if (!motif.inexact) {
        unit = motif.ru;
        tmp.resize(mlen);
        for (uint32_t i = 0; i < mlen; ++i) {tmp[i] = 0;}
    } else {
        unit = "";
        for (uint32_t i = 0; i < mlen; ++i) {
            unit += motif.ru.at(i);
            tmp.push_back(0);
            if (motif.label[i]) {
                m += 1;
                unit += motif.ru.at(i);
                tmp.push_back(1);
            }
        }
    }
    bool base_relax[m];
    std::copy(tmp.begin(), tmp.end(), base_relax);
    WPHMM_UNGAP* wphmm = new WPHMM_UNGAP(query.c_str(), unit.c_str(), debug, base_relax);

    // Run HMM - with the longest allele
    wphmm->set_ru(display_ru, inexact_label);
    wphmm->initialize();
    std::string v_path = wphmm->print_viterbi_path();
if (debug) {
    std::cerr << v_path << std::endl;
    std::cerr << "WPHMM_UNGAP log2 OR:\t" << wphmm->viterbi_score << '\n';
}
    wphmm->detect_range();
    seq_segment& rr = wphmm->focal_rr;
if (debug) {
    printf("Original st, ed, score: %d, %d, %.3f\n",motif.st, motif.ed, motif.viterbi_score);
}
        motif.st = get_pos_in_ref(rr.p_st);
        motif.ed = get_pos_in_ref(rr.p_ed);
if (debug) {
    printf("Updated st, ed, score: %d, %d, %.3f\n",motif.st, motif.ed, motif.viterbi_score);
}
        motif.viterbi_score = rr.score;
        motif.l = motif.ed - motif.st + 1; // length in ref
        rlen_l  = rr.p_ed - rr.p_st + 1;   // length in longest mosaic sequence
        n_ru_l  = rr.nr;                   // n_ru in longest mosaic sequence

    delete wphmm;
    wphmm = NULL;

    if (merged_longest_rr.length() > repeat_ref.length()) {
        query = repeat_ref;
        wphmm = new WPHMM_UNGAP(query.c_str(), unit.c_str(), debug, base_relax);
        wphmm->set_ru(display_ru, inexact_label);
        wphmm->initialize();
        wphmm->viterbi();
        wphmm->detect_range();
        if (wphmm->segments.size()>0) {
            seq_segment& rr = wphmm->focal_rr;
            score_r = rr.score;
            n_ru_r  = rr.nr;
            rlen_r  = rr.p_ed - rr.p_st + 1;
            if (debug) {
                printf("(Ref) length, nr, score: %d, %d, %.3f\n",rlen_r, n_ru_r, score_r);
            }
        }
        delete wphmm;
        wphmm = NULL;
    }

}










/**
 * Print object.
 */
void VNTR_candidate::print()
{
    // std::cerr << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n";
    // std::cerr << "VNTR Summary\n";
    // std::cerr << "ru           : " << motif.ru << "\n";
    // std::cerr << "\n";
    // std::cerr << "repeat_ref                    : " << repeat_ref << "\n";
    // std::cerr << "position                        : [" << rbeg1 << "," << rend1 << "]\n";
    // std::cerr << "reference repeat unit length    : " << rl << "\n";
    // std::cerr << "longest allele length           : " << ll << "\n";
    // std::cerr << "motif_concordance               : " << motif_concordance << "\n";
    // std::cerr << "repeat units                    : " << rl << "\n";
    // std::cerr << "exact repeat units              : " << no_exact_ru << "\n";
    // std::cerr << "total no. of repeat units       : " << total_no_ru << "\n";

};
