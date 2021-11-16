#ifndef RU_PROFILE_H
#define RU_PROFILE_H

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
    int32_t alphabet_size;

    Motif_exact(const char* _b, int32_t _s = 4) : Motif(_b) {
        mvec = new float[mlen];
        emtx = new float[mlen];
        alphabet_size = _s;
    }
    ~Motif_exact() {
        delete mvec;
        delete emtx;
    }

    void set_emission(float* _e) {
        for (int32_t i = 0; i < mlen; ++i) {
            mvec[i] = log2(1.-_e[i]) + log2(alphabet_size); // log ratio of match vs random
            emtx[i] = log2(_e[i]/(alphabet_size-1)) + log2(alphabet_size);
        }
    }
    void set_emission(float _e) {
        for (int32_t i = 0; i < mlen; ++i) {
            mvec[i] = log2(1.-_e) + log2(alphabet_size);
            emtx[i] = log2(_e/(alphabet_size-1)) + log2(alphabet_size);
        }
    }
};

class Motif_fuzzy : public Motif_exact
{
public:
    int8_t** srange; // single base repeat min max

    Motif_fuzzy(const char* _b, int8_t** _s = nullptr, int32_t _a = 4) : Motif_exact(_b, _a) {
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

    Motif_fuzzy_binary(const char* _b, int32_t _a = 4) : Motif_exact(_b, _a) {
        if_soft = (bool*) malloc(mlen);
        for (int32_t i = 0; i < mlen; ++i) {
            if_soft[i] = 0;
        }
    }
    Motif_fuzzy_binary(const char* _b, bool* _s, int32_t _a = 4) : Motif_exact(_b, _a) {
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

#endif
