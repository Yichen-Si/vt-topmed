#include "wphmm.h"

WPHMM::WPHMM(const char* _s, const char* _m, bool _debug, bool* _b) {
    debug = _debug;
    motif = new Motif_fuzzy_binary(_m);
    if (_b != nullptr) {
        motif->set_indicator(_b);
    }
    motif->set_emission(log2(.05));
    mlen = motif->mlen;
    L = strlen(_s);
    N = 0;
    B = 1;
    E = 2;
    C = 3;
    J = 4;
    n_structure_state = 5;
    std::vector<char> state_list{'N','B','E','C','J','M','I','D'};
    for (uint32_t i = 0; i < state_list.size(); ++i) {
        state_label[i] = state_list[i];
    }
    seq = (char*) malloc(sizeof(char) * L);
    memcpy(seq, _s, L);
    delta_B = 0.;
    delta_E = 0.;
    lambda = log2(.1);
    gamma  = log2(.2);
    eta    = log2(.1);
    trr = L / (L + 3.);
    tro = 1. - trr;
    trr = log2(trr);
    tro = log2(tro);
    tec = log2(0.5);
    tej = log2(0.5);
    W = (double**) malloc(sizeof(double*)*(L+1));
    M = (double**) malloc(sizeof(double*)*(L+1));
    I = (double**) malloc(sizeof(double*)*(L+1));
    D = (double**) malloc(sizeof(double*)*(L+1));
    for (int32_t i = 0; i <= L; ++i) {
        W[i] = (double*) malloc(sizeof(double)*(n_structure_state));
        M[i] = (double*) malloc(sizeof(double)*mlen);
        I[i] = (double*) malloc(sizeof(double)*mlen);
        D[i] = (double*) malloc(sizeof(double)*mlen);
    }
    viterbi_path = (uint32_t*) malloc(sizeof(uint32_t)*L);
    viterbi_path[0] = mlen*4+n_structure_state;
}

WPHMM::~WPHMM() {
    free(seq);
    for (int32_t i = 0; i <= L; ++i) {
        free(W[i]);
        free(M[i]);
        free(I[i]);
        free(D[i]);
    }
    free(W); free(M); free(I); free(D);
    free(viterbi_path);
}

void WPHMM::initialize() {
    W[0][N] = 0;
    W[0][B] = tro;
    W[0][E] = -DBL_MAX/2;
    W[0][J] = -DBL_MAX/2;
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
}

void WPHMM::viterbi() {
    uint32_t viterbi_mtx[L+1][mlen*3+n_structure_state];
    if (delta_B > 1./mlen) {delta_B = 1./mlen;}
    if (delta_E > 1./mlen) {delta_E = 1./mlen;}
    for (uint32_t i = 1; i <= L; ++i) {

        W[i][B] = W[i-1][N];
        viterbi_mtx[i][B] = N; // N->B
        if (W[i][B] < W[i-1][J]) {
            W[i][B] = W[i-1][J];
            viterbi_mtx[i][B] = J; // J->B
        }
        W[i][B] += tro;

        for (uint32_t k = 0; k < mlen; ++k) {
            uint32_t k_1 = (k==0) ? mlen-1 : k-1; // wraparound index
            // Match
            std::vector<double> tmp(4);
            tmp[0] = W[i][B] + log2(1.-delta_B*k);
            tmp[1] = M[i-1][k_1];
            tmp[2] = I[i-1][k_1];
            tmp[3] = D[i-1][k_1];
            viterbi_mtx[i][n_structure_state+k] = viterbi_mtx[i][B]; // B->M(k)
            M[i][k] = tmp[0];
            for (uint32_t j = 1; j <= 3; ++j) {
                if (tmp[j] > M[i][k]) {
                    viterbi_mtx[i][n_structure_state+k] = n_structure_state+mlen*(j-1)+k_1; // M/I/D(k-1)->M(k)
                    M[i][k] = tmp[j];
                }
            }
            // M[i][k] = *std::max_element(tmp.begin(), tmp.end());
            M[i][k] += (motif->base[k] == seq[i]) ? motif->mvec[k] : motif->emtx[k];
            // Insertion
            if (motif->base[k] == seq[i] && motif->if_soft[k]) {
                if (I[i-1][k] > M[i-1][k]) {
                    viterbi_mtx[i][n_structure_state+mlen+k] = n_structure_state+mlen+k; // I(k)->I(k)
                    I[i][k] = I[i-1][k];
                } else {
                    viterbi_mtx[i][n_structure_state+mlen+k] = n_structure_state+k; // M(k)->I(k)
                    I[i][k] = M[i-1][k];
                }
                // I[i][k] = (I[i-1][k] > M[i-1][k]) ? I[i-1][k] : M[i-1][k];
            } else {
                viterbi_mtx[i][n_structure_state+mlen+k] = n_structure_state+k; // M(k) -> I(k)
                I[i][k] = M[i-1][k] + lambda;
                if (I[i-1][k] + gamma > I[i][k]) {
                    viterbi_mtx[i][n_structure_state+mlen+k] = n_structure_state+mlen+k; // I(k)->I(k)
                    I[i][k] = I[i-1][k] + gamma;
                }
            }
            // Deletion
            if (D[i][k_1] > M[i][k_1]) {
                viterbi_mtx[i][n_structure_state+mlen*2+k] = n_structure_state+mlen*2+k_1; // D(k-1)->D(k)
                D[i][k] = D[i][k_1];
            } else {
                viterbi_mtx[i][n_structure_state+mlen*2+k] = n_structure_state+k_1; // M(k-1)->D(k)
                D[i][k] = M[i][k_1];
            }
            // D[i][k] = (M[i][k_1] > D[i][k_1]) ? M[i][k_1] : D[i][k_1];
            if (!motif->if_soft[k]) {
                D[i][k] += eta;
            }
        }

        W[i][N] = W[i-1][N] + trr;
        viterbi_mtx[i][N] = N; // N->N

        W[i][E] = M[i-1][0] + log2(1.-delta_E*(mlen-1));
        viterbi_mtx[i][E] = n_structure_state+0; // M(0)->E
        for (uint32_t k = 1; k < mlen; ++k) {
            double tmp = M[i-1][k] + log2(1.-delta_E*(mlen-k-1));
            if (tmp > W[i][E]) {
                viterbi_mtx[i][E] = n_structure_state+k; // M(k)->E
                W[i][E] = tmp;
            }
        }

        W[i][C] = W[i-1][C] + trr;
        viterbi_mtx[i][C] = C; // C->C
        if (W[i][C] < W[i][E] + tec) {
            viterbi_mtx[i][C] = viterbi_mtx[i][E]; // M->E->C
            W[i][C] = W[i][E] + tec;
        }

        W[i][J] = W[i-1][J] + trr;
        viterbi_mtx[i][J] = J; // J->J
        if (W[i][J] < W[i][E] + tej) {
            viterbi_mtx[i][J] = viterbi_mtx[i][E]; // M->E->J
            W[i][J] = W[i][E] + tej;
        }


    }
    viterbi_score = W[L][C] + tro - L * log2((double)L/(L+1.)) - log2(1./(L+1.));
    // Backtrack
    viterbi_path[L-1] = viterbi_mtx[L][C]; // Ends at C
    for (uint32_t i = L-1; i > 0; --i) {
        viterbi_path[i-1] = viterbi_mtx[i][viterbi_path[i]];
    }
}

std::string WPHMM::get_viterbi_path() {
    if (viterbi_path[0] > mlen*3+n_structure_state) {
        viterbi();
    }
    std::stringstream path;
    uint32_t key = 0;
    for (uint32_t i = 0; i < L; ++i) {
        key = (viterbi_path[i] < n_structure_state) ? viterbi_path[i] : n_structure_state + (uint32_t) ((viterbi_path[i]-n_structure_state)/mlen);
        path << state_label[key];
    }
    return path.str();
}
