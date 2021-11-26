#include "shift_align.h"

ShiftAlign::ShiftAlign(faidx_t* _f, bcf_hdr_t* h, bcf1_t* v, int32_t _m, int32_t _l, bool _d) : fai(_f), bound(_m), max_l(_l), debug(_d) {
    chrom = bcf_get_chrom(h, v);
    char** alleles = bcf_get_allele(v);
    extension = alleles[0];
    p_st = v->pos;
    p_ed = v->pos + extension.size() - 1; // Last pos of deleted base or same as p_st
    if (strlen(alleles[1]) > extension.size()) {
        extension = alleles[1];
    }
    extension = extension.substr(1);
    s_st = extension.length();

    s0 =  1.;
    s1 = -2.;
    g0 = -2.;
    if (max_l < bound * 3) {
        max_l = bound * 3;
    }
    if (max_l < s_st + bound) {
        max_l = s_st + bound;
    }
    delta_exact = 0;
}

void ShiftAlign::set_parameter(float _m, float _s, float _g) {
    s0 = _m;
    s1 = _s;
    g0 = _g;
}

void ShiftAlign::right_shift(int32_t &st0, int32_t &st, int32_t &ed, std::string& mseq) {
    // Read sequence
    int32_t seq_len;
    int32_t r_most = p_ed + max_l;
    std::string seq = extension + faidx_fetch_seq(fai, chrom, p_ed+1, r_most, &seq_len);
    // [0, s_st-1] is in inserted/deleted sequence, >= s_st is shared
    int32_t min_i, min_j, max_i, max_j;
// if (debug) {
//     printf("ShiftAlign::right_shift - p_st %d, p_ed %d, s_st %d, len %lu, %s\n", p_st, p_ed, s_st, seq.length(), extension.c_str());
//     std::cerr << seq_len << '\t' << seq.substr(s_st) << '\n';
// }
    shift(seq, min_i, min_j, max_i, max_j, mseq);
    ed = p_ed + max_i; // 0-based last matched base position
    st = p_ed + min_i; // 0-based first matched base in the shared region
    if (p_ed > p_st) { // Deletion
        st0 = p_ed - s_st + min_j; // 0-based first matched base in deleted region
    } else {           // Insertion
        st0 = min_j - s_st - 1; // negative if insertion is matched with sequence to the right. -x means the last x bases in insertion is matched
    }
    return;
}

void ShiftAlign::left_shift(int32_t &ed0, int32_t &st, int32_t &ed, std::string& mseq) {
    // Read sequence
    int32_t seq_len;
    int32_t l_most = p_st - max_l;
    // [max_l+1, ] is in inserted/deleted sequence, [0, max_l] is shared
    std::string seq = faidx_fetch_seq(fai, chrom, l_most, p_st, &seq_len) + extension;
    // [0, s_st-1] is in inserted/deleted sequence, >= s_st is shared
    std::reverse(seq.begin(), seq.end());
    int32_t min_i, min_j, max_i, max_j;
// if (debug) {
//     printf("ShiftAlign::left_shift - p_st %d, p_ed %d, s_st %d, len %lu, %s\n", p_st, p_ed, s_st, seq.length(), extension.c_str());
//     std::cerr << seq_len << '\t' << seq.substr(s_st) << '\n';
// }
    shift(seq, min_i, min_j, max_i, max_j, mseq);
    std::reverse(mseq.begin(), mseq.end());
    ed = p_st - min_i + 1; // 0-based last matched base position in shared region
    st = p_st - max_i + 1; // 0-based first matched base in the shared region
    if (p_ed > p_st) { // Deletion
        ed0 = p_st + min_j; // 0-based last matched base in deleted region
    } else {           // Insertion
        ed0 = min_j - s_st - 1; // negative if insertion is matched with sequence to the left. -x means the first x bases in insertion is matched
    }
    return;
}

void ShiftAlign::shift(std::string& seq, int32_t& min_i, int32_t& min_j, int32_t& max_i, int32_t& max_j, std::string& mseq) {
    int32_t m = seq.length();
    int32_t n = m - s_st;
    float lowerbd = -10.*m;
    // DP matrix - Initialize
    float H[n+1][m+1];

// if (debug) {
//     for (int32_t i = 0; i <= n; ++i) {
//         for (int32_t j = 0; j <= m; ++j) {
//             H[i][j] = 0;
//         }
//     }
// }

    for (int32_t i = 0; i <= n; ++i) {
        H[i][0] = 0.;
    }
    for (int32_t j = 0; j <= m; ++j) {
        H[0][j] = 0.;
    }
    for (int32_t i = 1; i <= n; ++i) {
        H[i][i+s_st] = lowerbd;
    }
    for (int32_t i = bound+1; i <= n; ++i) {
        H[i][i-bound] = lowerbd;
    }


    char trace[n+1][m+1];
    max_i = 0;
    max_j = 0;
    float max_score = 0;
    // Bounded DP
    for (int32_t i = 1; i <= n; ++i) {
        int32_t j_st = (i - bound > 0) ? i - bound + 1 : 1;
        for (int32_t j = j_st; j < i + s_st; ++j) {
            if (seq.at(s_st+i-1) == seq.at(j-1)) {
                H[i][j] = H[i-1][j-1] + s0; // M->M
            } else {
                H[i][j] = H[i-1][j-1] + s1;
            }
            trace[i][j] = '1';
            if (H[i][j-1] + g0 > H[i][j]) {
                H[i][j] = H[i][j-1] + g0;
                trace[i][j] = '2';
            }
            if (H[i-1][j] + g0 > H[i][j]) {
                H[i][j] = H[i-1][j] + g0;
                trace[i][j] = '3';
            }
            // if (0 > H[i][j]) {
            //     H[i][j] = 0;
            //     trace[i][j] = '0';
            // }
            if (H[i][j] > max_score) {
                max_i = i;
                max_j = j;
                max_score = H[i][j];
            }
        }
    }

// if (debug) {
//     std::cerr << "   ";
//     for (int32_t j = 1; j <= m; ++j) {
//         std::cerr << ' ' << seq.at(j-1);
//     }
//     std::cerr << '\n';
//     for (int32_t i = 0; i <= n; ++i) {
//         std::cerr << seq.at(s_st+i-1);
//         for (int32_t j = 0; j <= m; ++j) {
//             std::cerr << ' ' << H[i][j];
//         }
//         std::cerr << '\n';
//     }
// }

    if (H[s_st][s_st] >= s0*s_st) {
        delta_exact = 1;
    } else if (H[s_st+1][s_st] >= s0*s_st+g0) {
        int32_t i = s_st+1;
        int32_t j = s_st;
        while (i > 0 && j > 0 && trace[i][j] == '1') {
            i -= 1;
            j -= 1;
        }
        if (j == 0) {
            delta_exact = 1;
        }
    }

    // Backtrack
    min_i = 0;
    min_j = 0;
    mseq = "";
    if (max_i == 0 || max_j == 0) {
        mseq = "";
        return;
    }
    int32_t i = max_i, j = max_j;
    std::string aseq;
    float score = max_score;
    // while (score > 0) {
    while (i > 0 && j > 0) {
        min_j = j;
        min_i = i;
        switch(trace[i][j]) {
            case '1' :
                mseq += seq.at(j-1);
                if (i+s_st > max_j) {
                    aseq += seq.at(s_st+i-1);
                }
                i -= 1;
                j -= 1;
                break;
            case '2' :
                mseq += seq.at(j-1);
                j -= 1;
                break;
            case '3' :
                if (i+s_st > max_j) {
                    aseq += seq.at(s_st+i-1);
                }
                i -= 1;
                break;
        }
        score = H[i][j];
    }
    std::reverse(mseq.begin(), mseq.end());
    std::reverse(aseq.begin(), aseq.end());
// if (debug) {
//     std::cerr << mseq << '\t' << aseq << '\n';
//     std::cerr << min_i << '\t' << min_j << '\t' << max_i << '\t' << max_j << '\n';
//     printf("%s\n", seq.substr(min_j-1, max_j-min_j+1).c_str());
//     printf("%s\n", seq.substr(s_st+min_i-1, max_i-min_i+1).c_str());
// }
    mseq += aseq;
    return;
}
