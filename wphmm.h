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
            mvec[i] = log2(1.-_e[i]) + log2(Alphabet_size);
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
        if_soft = (bool*) malloc(sizeof(bool) * mlen);
        for (int32_t i = 0; i < mlen; ++i) {
            if_soft[i] = 0;
        }
    }
    Motif_fuzzy_binary(const char* _b, bool* _s) : Motif_exact(_b) {
        if_soft = (bool*) malloc(sizeof(bool) * mlen);
        set_indicator(_s);
    }
    ~Motif_fuzzy_binary() {
        free(if_soft);
    }
    void set_indicator(bool* _s) {
        memcpy(if_soft, _s, mlen);
    }
};

class WPHMM
{
public:
    bool debug;
    // Repeat model
    uint32_t mlen;
    Motif_fuzzy_binary* motif;
    float delta_B, delta_E;   // control start/end repeat region
    float lambda, gamma, eta; // control insertion open/continue & deletion
    float trr, tro, tec, tej;
    // Input sequence
    int32_t L; // length
    char* seq;
    // DP matrix
    double** W; // (L+1) x 5, Flanking/Junction
    double** M; // (L+1) x m, Matching
    double** I; // (L+1) x m, Insertion
    double** D; // (L+1) x m, Deletion
    size_t N, B, E, C, J;
    int32_t n_structure_state;
    uint32_t* viterbi_path;
    double viterbi_score;
    std::map<uint32_t, char> state_label;

WPHMM(const char* _s, const char* _m, bool _debug = 0, bool* _b = nullptr);
~WPHMM();
void initialize();

void viterbi();
std::string get_viterbi_path();

};

#endif
