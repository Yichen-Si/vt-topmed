#include "wphmm_ungapped.h"

WPHMM_UNGAP::WPHMM_UNGAP(const char* _s, const char* _m, bool _debug, bool* _b, bool _c) {
    debug = _debug;
    if (_c && (_b != nullptr)) {
        expanding_motif(_m, _b);
    } else {
        motif = new Motif_fuzzy_binary(_m);
        if (_b != nullptr) {
            motif->set_indicator(_b);
        }
    }
    motif->set_emission(.05); // Mismatch probability
    mlen = motif->mlen;
    L = strlen(_s);
    N = 0;
    B = 1;
    E = 2;
    C = 3;
    n_structure_state = 4;
    std::vector<char> state_list{'N','B','E','C','M','I','D'};
    for (uint32_t i = 0; i < state_list.size(); ++i) {
        state_label[i] = state_list[i];
    }
    int32_t nf = 0;
    for (uint32_t i = 0; i < mlen; ++i) {
        nf += motif->if_soft[i];
    }
    seq = (char*) malloc(sizeof(char) * L);
    memcpy(seq, _s, L);
    delta_B = log2(1./(mlen-nf)); // B -> M, uniform
    delta_E = log2(.5); // M -> E, uniform
    lambda = log2(.05); // insertion openning
    gamma  = log2(.15); // insertion elongation
    eta    = log2(.05); // deletion openning
    zeta   = log2(.15); // deletion elongation
    tmm = log2(.4); // minus probability of indel or end
    tdm = log2(.85);
    tim = log2(.85);
    trr = L / (L + 3.);
    tro = 1. - trr;
    trr = log2(trr);
    tro = log2(tro);
    tec = 0.; // E->C
    W = (double**) malloc(sizeof(double*)*(L+1));
    M = (double**) malloc(sizeof(double*)*(L+1));
    I = (double**) malloc(sizeof(double*)*(L+1));
    D = (double**) malloc(sizeof(double*)*(L+1));
    viterbi_mtx = (uint32_t**) malloc(sizeof(uint32_t*)*(L+1));
    for (int32_t i = 0; i <= L; ++i) {
        W[i] = (double*) malloc(sizeof(double)*(n_structure_state));
        M[i] = (double*) malloc(sizeof(double)*mlen);
        I[i] = (double*) malloc(sizeof(double)*mlen);
        D[i] = (double*) malloc(sizeof(double)*mlen);
        viterbi_mtx[i] = (uint32_t*) malloc(sizeof(uint32_t)*(n_structure_state+3*mlen));
    }
    for (auto &v : alphabet) {
        random_base[v] = log2(1./Alphabet_size); // not used yet
    }
}

WPHMM_UNGAP::~WPHMM_UNGAP() {
    delete motif;
    motif = NULL;
    free(seq);
    for (int32_t i = 0; i <= L; ++i) {
        free(W[i]);
        free(M[i]);
        free(I[i]);
        free(D[i]);
        free(viterbi_mtx[i]);
    }
    free(W); free(M); free(I); free(D); free(viterbi_mtx);
}

void WPHMM_UNGAP::expanding_motif(const char* _m, bool* _b) {
    int32_t m = strlen(_m);
    std::string unit;
    std::vector<bool> tmp;
    for (uint32_t i = 0; i < m; ++i) {
        unit += _m[i];
        tmp.push_back(0);
        if (_b[i]) {
            unit += _m[i];
            tmp.push_back(1);
        }
    }
    m = unit.size();
    bool base_relax[m];
    std::copy(tmp.begin(), tmp.end(), base_relax);
    motif = new Motif_fuzzy_binary(unit.c_str(), base_relax);
}

void WPHMM_UNGAP::initialize() {
    W[0][N] = 0;
    W[0][B] = tro;
    W[0][E] = -DBL_MAX/2;
    W[0][C] = -DBL_MAX/2;
    for (int32_t i = 0; i <= L; ++i) {
        M[i][mlen-1] = -DBL_MAX/2;
        I[i][mlen-1] = -DBL_MAX/2;
        D[i][mlen-1] = -DBL_MAX/2;
    }
    for (int32_t k = 0; k < mlen; ++k) {
        M[0][k] = -DBL_MAX/2;
        I[0][k] = -DBL_MAX/2;
        D[0][k] = -DBL_MAX/2;
    }
    for (int32_t k = 0; k < mlen; ++k) {
        D[1][k] = -DBL_MAX/2;
    }
}

void WPHMM_UNGAP::viterbi() {

    for (uint32_t i = 1; i <= L; ++i) {

        W[i][B] = W[i-1][N] + tro;
        viterbi_mtx[i][B] = N; // N->B

        for (uint32_t k = 0; k < mlen; ++k) {
            uint32_t k_1 = (k==0) ? mlen-1 : k-1; // wraparound index
            // Match
            std::vector<double> tmp(4);
            tmp[0] = W[i][B] + ((motif->if_soft[k]) ? -DBL_MAX/2 : delta_B);
            tmp[2] = I[i-1][k_1] + tim;
            tmp[3] = D[i-1][k_1] + tdm;
            tmp[1] = M[i-1][k_1] + tmm;
            viterbi_mtx[i][n_structure_state+k] = viterbi_mtx[i][B]; // B->M(k)
            M[i][k] = tmp[0];
            for (uint32_t j = 1; j <= 3; ++j) {
                if (tmp[j] > M[i][k]) {
                    viterbi_mtx[i][n_structure_state+k] = n_structure_state+mlen*(j-1)+k_1; // M/I/D(k-1)->M(k)
                    M[i][k] = tmp[j];
                }
            }
            // M[i][k] = *std::max_element(tmp.begin(), tmp.end());
            M[i][k] += (motif->base[k] == seq[i-1]) ? motif->mvec[k] : motif->emtx[k];
            // Insertion
            if (motif->base[k] == seq[i-1] && motif->if_soft[k]) { // Tolerated insertion (treat as match)
                if (I[i-1][k] > M[i-1][k]) {
                    viterbi_mtx[i][n_structure_state+mlen+k] = n_structure_state+mlen+k; // I(k)->I(k)
                    I[i][k] = I[i-1][k] + motif->mvec[k] + tmm;
                } else {
                    viterbi_mtx[i][n_structure_state+mlen+k] = n_structure_state+k; // M(k)->I(k)
                    I[i][k] = M[i-1][k] + motif->mvec[k] + tmm;
                }
                // I[i][k] = (I[i-1][k] > M[i-1][k]) ? I[i-1][k] : M[i-1][k];
            } else { // Intolerated insertion (assume conditional on insertion base is random)
                viterbi_mtx[i][n_structure_state+mlen+k] = n_structure_state+k; // M(k) -> I(k)
                I[i][k] = M[i-1][k] + lambda;
                if (I[i-1][k] + gamma > I[i][k]) {
                    viterbi_mtx[i][n_structure_state+mlen+k] = n_structure_state+mlen+k; // I(k)->I(k)
                    I[i][k] = I[i-1][k] + gamma;
                }
            }
        }

        // Deletion (the inexact case is temporary)
        if (i > 1) {
            for (uint32_t k = 0; k < mlen; ++k) {
                uint32_t k_1 = (k==0) ? mlen-1 : k-1; // wraparound index
                if (k == 0) {
                    viterbi_mtx[i][n_structure_state+mlen*2+k] = n_structure_state+k_1; // M(k-1)->D(k)
                    D[i][k] = M[i][k_1] + eta;
                    if (motif->if_soft[k_1] && D[i][k] < M[i][k_1-1] + eta) {
                        viterbi_mtx[i][n_structure_state+mlen*2+k] = n_structure_state+k_1-1;
                        D[i][k] = M[i][k_1-1] + eta;
                    }
                } else {
                    if (motif->if_soft[k]) { // Tolerated deletion
                        if (D[i][k_1] > M[i][k_1]) {
                            viterbi_mtx[i][n_structure_state+mlen*2+k] = n_structure_state+mlen*2+k_1; // D(k-1)->D(k)
                            D[i][k] = D[i][k_1];
                        } else {
                            viterbi_mtx[i][n_structure_state+mlen*2+k] = n_structure_state+k_1; // M(k-1)->D(k)
                            D[i][k] = M[i][k_1];
                        }
                    } else { // Intolerated deletion
                        if (D[i][k_1] + zeta > M[i][k_1] + eta) {
                            viterbi_mtx[i][n_structure_state+mlen*2+k] = n_structure_state+mlen*2+k_1; // D(k-1)->D(k)
                            D[i][k] = D[i][k_1] + zeta;
                        } else {
                            viterbi_mtx[i][n_structure_state+mlen*2+k] = n_structure_state+k_1; // M(k-1)->D(k)
                            D[i][k] = M[i][k_1] + eta;
                        }
                    }
                }
            }
        }

        W[i][N] = W[i-1][N] + trr;
        viterbi_mtx[i][N] = N; // N->N

        W[i][E] = M[i-1][0];
        viterbi_mtx[i][E] = n_structure_state+0; // M(0)->E
        for (uint32_t k = 1; k < mlen; ++k) {
            double tmp = M[i-1][k] ;
            if (tmp > W[i][E]) {
                viterbi_mtx[i][E] = n_structure_state+k; // M(k)->E
                W[i][E] = tmp;
            }
        }
        W[i][E] += delta_E; // log(P(Mk->E))

        W[i][C] = W[i-1][C] + trr;
        viterbi_mtx[i][C] = C; // C->C
        if (W[i][C] < W[i][E] + tec) {
            viterbi_mtx[i][C] = viterbi_mtx[i][E]; // M->E->C
            W[i][C] = W[i][E] + tec;
        }
    }

    // Backtrack, segment viterbi path into segments
    viterbi_path.clear();
    vpath.clear();// For inexact RU, collapsing placeholder I/D
    vmle.clear(); // vmle[i] = max log(P(x[0..i]|Model)) + const.
    vmle.resize(L);

    // Possible end state: C, M, I
    viterbi_score = W[L][C];
    viterbi_path.push_back(C); // Ends at C
    for (int32_t k = mlen-1; k >= 0; --k) {
        if (I[L][k] >= viterbi_score) {
            viterbi_score = I[L][k];
            viterbi_path[0] = n_structure_state+mlen*1+k;
        }
        if (M[L][k] >= viterbi_score) {
            viterbi_score = M[L][k];
            viterbi_path[0] = n_structure_state+mlen*0+k;
        }
    }
    vmle[L-1] = viterbi_score;
    viterbi_score -= L * log2((double)L/(L+1.)) - log2(1./(L+1.));
    vpath.push_back(viterbi_path.back());
    uint32_t prev_state = viterbi_path.back();
    uint32_t i = L;
    while (i > 1) {
        uint32_t next_state = viterbi_mtx[i][prev_state];
        viterbi_path.push_back(next_state);
        if (prev_state < n_structure_state+mlen*2) {
            i -= 1;
        }
        prev_state = next_state;
        if (next_state >= n_structure_state) {
            uint32_t ptype = (next_state-n_structure_state)/mlen;
            uint32_t k = (next_state-n_structure_state) % mlen;
            if (ptype == 0) {
                vmle[i-1] = M[i][k];
            } else if (ptype == 1) {
                vmle[i-1] = I[i][k];
            }
            if (motif->if_soft[k] && ptype == 2) { // Tolerated deletion
                continue;
            }
            if (motif->if_soft[k] && ptype == 1 && seq[i-1] == motif->base[k]) { // Tolerated insertion, code as M
                vpath.push_back(next_state-mlen);
                continue;
            }
        } else {
            vmle[i-1] = W[i][next_state];
        }
        vpath.push_back(next_state);
    }
    std::reverse(vpath.begin(), vpath.end());
    std::reverse(viterbi_path.begin(), viterbi_path.end());
}

bool WPHMM_UNGAP::count_ru() {

    ru_complete.resize(L);
    std::fill(ru_complete.begin(), ru_complete.end(), -1);
    if (ru.size() == 1) {
        ru_complete[0] = (int32_t) (seq[0] == ru.at(0));
        for (uint32_t i = 1; i < L; ++i) {
            ru_complete[i] = ru_complete[i-1] + (int32_t) (seq[i] == ru.at(0));
        }
        return 0;
    }
    uint32_t i = 0; // points to position in vpath
    uint32_t j = 0; // points to position in seq
    while (i<vpath.size() && vpath[i] < n_structure_state) {
        ru_complete[j] = 0;
        j++;
        i++;
    }
    if (i==vpath.size()) {
        return 1;
    }
    // find anchor
    int32_t acr = (vpath[i]-n_structure_state) % mlen;
    if (motif->if_soft[acr]) {
        acr--;
    }
    int32_t racr = (mlen+acr-1) % mlen;
    if (motif->if_soft[racr]) {
        racr--;
    }
    if (vpath[i] >= n_structure_state+mlen || racr==acr) {
        fprintf(stderr, "[%s:%d %s] Motif is ill defined, cannot find anchor %s, %s; %d, %d; %d, %d.\n", __FILE__, __LINE__, __FUNCTION__,motif->base, ru.c_str(), acr, racr, i, j);
        std::cerr << print_viterbi_path() << '\n';
        exit(1);
    }
    bool pre_state = 0;
    if (j == 0) { // The first base is in M
        ru_complete[j] = 0;
        pre_state = 1;
        i += 1;
        j += 1;
    }
    while (i < vpath.size()) {
        int32_t k = -1;
        int32_t ptype = vpath[i];
        if (vpath[i] >= n_structure_state) {
            k = (vpath[i] - n_structure_state) % mlen;
            ptype = (vpath[i] - n_structure_state) / mlen;
        }
        if (k < 0) {
            pre_state = 0;
        } else if (ptype != 0 || (ptype == 0 && seq[j]!=motif->base[k])) { // I/D/m
            pre_state = 0;
        } else if (pre_state == 1 && k == racr && ptype==0) { // M[m-1], a ru completed
            ru_complete[j] = ru_complete[j-1] + 1;
            pre_state = 0;
        } else if (pre_state == 0 && ptype==0) { // M[0], try to start a ru
            acr = k;
            if (motif->if_soft[acr]) {
                acr--;
            }
            racr = (mlen+acr-1) % mlen;
            if (motif->if_soft[racr]) {
                racr--;
            }
            pre_state = 1;
        }
        if (k>=0 && ptype==2) { // D

        } else {
            if (ru_complete[j]<0) {ru_complete[j] = ru_complete[j-1];}
            j++;
        }
        i++;
    }
    return 0;
}

/**
    Separate the segment consists of only M/I/D, get scores
*/
void WPHMM_UNGAP::detect_range() {

    if (ru_complete.size() < L) {
        if (count_ru()) {
            return;
        }
    }
    int32_t s_ed = -1, s_st = -1;
    int32_t n_ru = 0, n_ins = 0;
    int32_t pre_state = -1;
    uint32_t offset = 10;
    int32_t i = L;
    int32_t j = vpath.size()-1;
    uint32_t next_state = vpath[j];
    uint32_t pre_pos0 = L;
    while (i > 0) {
        uint32_t ptype = next_state;
        int32_t k = -1;
        if (next_state >= n_structure_state) {
            ptype = offset + (next_state-n_structure_state)/mlen;
            k = (next_state-n_structure_state) % mlen;
        }
        uint32_t pos0 = i - 1;
        if (pre_state < 0 && ptype == offset) {
            // Starting a matching segment
            s_ed = pos0;
            pre_state = 1;
            n_ins = 0;
        } else if (pre_state == 1 && ptype == N) {
            // Ending a M segment
            s_st = pre_pos0;
            n_ru = ru_complete[s_ed]-ru_complete[s_st];
            seq_segment seg(s_st, s_ed, n_ru, s_ed-s_st+1-n_ins);
            seg.score = vmle[s_ed] - vmle[s_st-1] - (s_ed - s_st) * tmm;
            segments.push_back(seg);
            pre_state = 0;
            n_ins = 0;
        } else if (pre_state == 1 && ptype > offset) {
            n_ins++;
        }
        pre_pos0 = pos0;
        if (ptype != offset + 2) {
            i -= 1;
        }
        j -= 1;
        next_state = vpath[j];
    }
    if (pre_state == 1) {
        // Ending last M segment
        s_st = 0;
        n_ru = ru_complete[s_ed];
        seq_segment seg(s_st, s_ed, n_ru, s_ed-s_st+1-n_ins);
        seg.score = vmle[s_ed] - (s_ed - s_st) * tmm;
        segments.push_back(seg);
    }
    if (segments.size() == 0) {return;}
    focal_rr = segments[0];
    if (segments.size() > 1) {
        warning("WPHMM_UNGAP::select_segment - multiple segments exist");
        for (uint32_t i = 1; i < segments.size(); ++i) {
            if (segments[i].score > focal_rr.score) {
                focal_rr = segments[i];
            }
        }
    }
}

std::string WPHMM_UNGAP::get_viterbi_path() {
    if (viterbi_path.size() < L) {
        viterbi();
    }
    std::stringstream path;
    uint32_t key = 0;
    for (uint32_t i = 0; i < viterbi_path.size(); ++i) {
        key = (viterbi_path[i] < n_structure_state) ? viterbi_path[i] : n_structure_state + (uint32_t) ((viterbi_path[i]-n_structure_state)/mlen);
        path << state_label[key];
    }
    return path.str();
}

std::string WPHMM_UNGAP::print_viterbi_path() {
    if (vpath.size() < L) {
        viterbi();
    }
    std::stringstream path;
    uint32_t key = 0;
    for (uint32_t i = 0; i < vpath.size(); ++i) {
        key = (vpath[i] < n_structure_state) ? vpath[i] : n_structure_state + (uint32_t) ((vpath[i]-n_structure_state)/mlen);
        path << state_label[key];
    }
    return path.str();
}
