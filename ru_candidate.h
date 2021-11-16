#ifndef RU_CANDIDATE_H
#define RU_CANDIDATE_H

#include <cstdlib>
#include <cstdint>
#include <string>
#include <cmath>
#include <cfloat>
#include <utility>
#include <vector>
#include <set>
#include <list>
#include <iostream>
#include <utility>
#include <algorithm>

// inline bool comparePtrToCandyMotify(candidate_fuzzy_motif* a, candidate_fuzzy_motif* b) {
//     return (a->viterbi_score > b->viterbi_score);
// }

class candidate_fuzzy_motif
{
public:
    uint32_t mlen;
    std::string ru;
    std::vector<bool> label;
    bool inexact;
    int32_t st, ed; // zero based, inclusive. genome position of repeat region
    int32_t n_ru;   // number of (noninterrupted) ru matched
    double viterbi_score, concordance;
    int32_t l;
    bool insertion_relevant;
    std::string insertion;

    candidate_fuzzy_motif(std::string _u, std::vector<bool> _w, int32_t _s, int32_t _e, int32_t _n, double _v, int32_t _l = -1, int32_t _c = 0) : ru(_u), label(_w), st(_s), ed(_e), n_ru(_n), viterbi_score(_v), l(_l) {
            mlen = ru.length();
            inexact = 0;
            for (uint32_t i = 0; i < mlen; ++i) {
                inexact = inexact || label[i];
            }
            if (l < 0) {
                l = ed-st+1;
            }
            insertion_relevant=0;
            concordance = (double) _c / l;
            canonical();
        }
    candidate_fuzzy_motif() {}

    bool operator<(const candidate_fuzzy_motif & rhs) const {
        return (viterbi_score > rhs.viterbi_score);
    }

    /**
    * Rotate the RU to a canonical phase
    * Caution - currently only this canonical form will be stored
    * not necessarily the first full RU that occurs in the genome
    * because the WPHMM can start and end at any position
    */
    void canonical() {
        std::vector<std::pair<std::string, uint32_t> > w;
        for (uint32_t i = 0; i < mlen; ++i) {
            w.push_back( std::make_pair( ru.substr(i)+ru.substr(0,i), i ) );
        }
        std::sort(w.begin(), w.end());
        if (w[0].second == 0) {
            return;
        }
        uint32_t phase = w[0].second;
        if (inexact) {
            std::vector<bool> tmp(mlen, 0);
            for (uint32_t i = 0; i < mlen; ++i) {
                tmp[i] = label[(phase + i) % mlen];
            }
            label = tmp;
        }
        ru = w[0].first;
    }

    /**
     * If two repeat models are compatible.
     * 0 - Equivalent; 3 - Incompatible;
     * 1 - Compatible, non-inclusive;
     * 2 - Inclusive.
     */
    int32_t ru_compare(candidate_fuzzy_motif& rhs) {
        if (ru == rhs.ru && label == rhs.label) {
            return 0; // Equivalent
        }
        std::set<char> s1, s2;
        std::for_each(ru.begin(), ru.end(), [&s1] (char c) -> void { s1.insert(c);});
        std::for_each(rhs.ru.begin(), rhs.ru.end(), [&s2] (char c) -> void { s2.insert(c);});
        if (s1 != s2) {
            return 3;
        }
        if (ru.length() != rhs.ru.length()) {
            // Check if one includes the other
            if (ru.find(rhs.ru)!=std::string::npos || rhs.ru.find(ru)!=std::string::npos) {
                return 2;
            }
            return 3; // Repeat unit incompatible
        }
        bool equal_seq = (ru == rhs.ru);
        uint32_t phase = 1;
        while (phase < mlen && (!equal_seq)) {
            std::string cru = ru.substr(phase, mlen-phase)+ru.substr(0,phase);
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
            if (label[(i+phase)%mlen] != rhs.label[i]) {
                equal_cir = 0;
                break;
            }
        }
        if (equal_cir) {
            return 0; // Circular equivalent
        }
        // Compatible
        return 1;
    }
};

struct candidate_unit {
    std::string ru;
    bool inexact;
    std::vector< std::pair<int32_t, int32_t> > variable_base;
    candidate_unit(std::string _s, bool _i = 0) : ru(_s), inexact(_i) {}
    void check() {
        if (inexact && variable_base.size() != ru.size()) {
            inexact = 0;
            variable_base.clear();
            return;
        }
        int32_t ct = 0;
        for (size_t i = 0; i < ru.size(); ++i) {
            if (variable_base[i].first > variable_base[i].second) {
                int32_t tmp = variable_base[i].first;
                variable_base[i].first  = variable_base[i].second;
                variable_base[i].second = tmp;
            }
            ct += (variable_base[i].first < variable_base[i].second);
        }
        if (ct == 0) {
            inexact = 0;
            variable_base.clear();
        } else {
            inexact = 1;
        }

    }
    bool operator<(const candidate_unit & rhs) const {
        if (inexact != rhs.inexact) {
            return (inexact < rhs.inexact);
        }
        if (ru.size() != rhs.ru.size()) {
            return (ru.size() < rhs.ru.size());
        }
        if (ru.compare(rhs.ru) == 0) {return false;}
        if (ru.size() > 1) {
            for (size_t i = 1; i < ru.size(); ++i) {
                if (rhs.ru.compare(ru.substr(i)+ru.substr(0,i)) == 0) {
                    return false;
                }
            }
        }
        return (ru.compare(rhs.ru) < 0);
    }
};

#endif
