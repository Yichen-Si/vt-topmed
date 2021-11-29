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
#include "ru_profile.h"

#define Alphabet_size 4 // alphabed size

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

WPHMM(const char* _s, const char* _m, bool _debug = 0, bool* _b = nullptr, bool _c = 0);
~WPHMM();
void expanding_motif(const char* _m, bool* _b);
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
