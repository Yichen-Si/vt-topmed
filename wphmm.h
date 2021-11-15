#ifndef WPHMM_H
#define WPHMM_H

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <float.h>
#include <map>
#include <vector>
#include <cstring>
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include "Error.h"

#define Alphabet_size 4 // alphabed size

class Motif {
public:
    char* base;
    uint32_t mlen;

    Motif(const char* _b) {
        mlen = strlen(_b);
        base = strdup(_b);
    }
    ~Motif() {
        free(base);
    }
};

class Motif_exact : public Motif {
public:
    float* mvec; // log score for match
    float* emtx; // log score for mismatch (currently single mismatch prob)

    Motif_exact(const char* _b) : Motif(_b) {
        mvec = new float[mlen];
        emtx = new float[mlen];
    }
    ~Motif_exact() {
        delete mvec;
        delete emtx;
    }

    void set_emission(float* _e) {
        for (int32_t i = 0; i < mlen; ++i) {
            mvec[i] = log2(1.-_e[i]) + log2(Alphabet_size); // log ratio of match vs random
            emtx[i] = log2(_e[i]/(Alphabet_size-1)) + log2(Alphabet_size);
        }
    }
    void set_emission(float _e) {
        for (int32_t i = 0; i < mlen; ++i) {
            mvec[i] = log2(1.-_e) + log2(Alphabet_size);
            emtx[i] = log2(_e/(Alphabet_size-1)) + log2(Alphabet_size);
        }
    }
};

class Motif_fuzzy : public Motif_exact
{
public:
    int8_t** srange; // single base repeat min max

    Motif_fuzzy(const char* _b, int8_t** _s = nullptr) : Motif_exact(_b) {
        srange = (int8_t**) malloc(sizeof(int8_t*) * mlen);
        if (_s == nullptr) {
            for (int32_t i = 0; i < mlen; ++i) {
                srange[i] = (int8_t*) malloc(sizeof(int8_t)*2);
            }
        } else {
            for (int32_t i = 0; i < mlen; ++i) {
                srange[i] = (int8_t*) malloc(sizeof(int8_t)*2);
                memcpy(srange[i], _s[i], 2);
            }
        }
    }
    ~Motif_fuzzy(){
        for (int32_t i = 0; i < mlen; ++i) {
            free(srange[i]);
        }
        free(srange);
    }
    void set_range(int8_t** _s) {
        for (int32_t i = 0; i < mlen; ++i) {
            memcpy(srange[i], _s[i], 2);
        }
    }
};

class Motif_fuzzy_binary : public Motif_exact
{
public:
    bool* if_soft; // If lenient with same base indel at each pos

    Motif_fuzzy_binary(const char* _b) : Motif_exact(_b) {
        if_soft = (bool*) malloc(mlen);
        for (int32_t i = 0; i < mlen; ++i) {
            if_soft[i] = 0;
        }
    }
    Motif_fuzzy_binary(const char* _b, bool* _s) : Motif_exact(_b) {
        if_soft = (bool*) malloc(mlen);
        set_indicator(_s);
    }
    ~Motif_fuzzy_binary() {
        free(if_soft);
    }
    void set_indicator(bool* _s) {
        memcpy(if_soft, _s, mlen);
    }

};

struct seq_segment {
    uint32_t p_st, p_ed; // 0-based, inclusive
    uint32_t nr;
    int32_t match_motif;
    uint32_t l;
    double score;
    seq_segment(uint32_t _s, uint32_t _e, uint32_t _n=0, int32_t _m=0) : p_st(_s), p_ed(_e), nr(_n), match_motif(_m) {
        l = p_ed - p_st + 1;
    }
    seq_segment() {}
};

class WPHMM
{
public:
    bool debug;
    // Repeat model
    uint32_t mlen;
    std::string ru; // collapsed ru
    std::vector<bool> inexact_label;
    Motif_fuzzy_binary* motif;
    float delta_B, delta_E;   // control start/end repeat region
    float lambda, gamma, eta, zeta; // control insertion/deletion open/continue
    float trr, tro, tec, tej, tmm, tim, tdm;
    // Input sequence
    int32_t L; // length
    char* seq;
    // DP matrix
    double** W; // (L+1) x 5, Flanking/Junction
    double** M; // (L+1) x m, Matching
    double** I; // (L+1) x m, Insertion
    double** D; // (L+1) x m, Deletion
    uint32_t** viterbi_mtx; // (L+1) x (5+mlen*3)
    size_t N, B, E, C, J;
    int32_t n_structure_state;
    std::vector<uint32_t> viterbi_path;
    std::vector<uint32_t> vpath; // For display
    std::map<uint32_t, char> state_label;
    // Statistics
    std::vector<int32_t> ru_complete;
    double viterbi_score;
    std::vector<double> vmle;
    uint32_t rlen; // length of repeat region
    uint32_t rbeg, rend; // begin and end of the identified repeat region (relative to seq, 0-based)
    std::vector<seq_segment> segments;
    seq_segment focal_rr;
    // Temporary: iid random null
    std::map<char, float> random_base;
    std::vector<char> alphabet{'A','C','G','T'};

WPHMM(const char* _s, const char* _m, bool _debug = 0, bool* _b = nullptr);
~WPHMM();

void set_ru(std::string& _s, std::vector<bool>& _v) {
    ru = _s;
    inexact_label = _v;
}
void initialize();
void set_emission(float* _e) {motif->set_emission(_e);}
void set_emission(float _e)  {motif->set_emission(_e);}

void viterbi();
bool count_ru();
void detect_range();
int32_t select_segment(int32_t left,int32_t right);

std::string get_viterbi_path();
std::string print_viterbi_path();

};

#endif
