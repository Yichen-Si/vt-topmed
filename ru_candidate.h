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
    * If this inexact motif covers an exact one
    */
    bool explains(candidate_fuzzy_motif& rhs) {
        if (!inexact || rhs.inexact) {
            return 0;
        }
        if (mlen > rhs.mlen) {
            return 0;
        }
        uint32_t i = 0, j = 0, pre_j = 0;
        bool flag = 1;
        while(i < rhs.mlen) {
            if (j < mlen && rhs.ru[i] == ru[j]) {
                pre_j = j;
                i++;
                j++;
            } else if (label[pre_j] && rhs.ru[i] == ru[pre_j]) {
                i++;
            } else {
                flag = 0;
                break;
            }
        }
        return flag;
    }

    /**
    * If this motif includes rhs as an instance
    * Or this motif consists of k copies of rhs, with k fractional but >= 1
    */
    bool contains(candidate_fuzzy_motif& rhs) {
        if (inexact != rhs.inexact) {
            return explains(rhs);
        }
        if (mlen == rhs.mlen) {
            return ru == rhs.ru;
        }
        if (ru.find(rhs.ru)==std::string::npos) {
            return 0;
        }
        std::string tmp = rhs.ru;
        uint32_t i = 0;
        while (tmp.length() < mlen) {
            tmp += rhs.ru.at((i+rhs.mlen)%rhs.mlen);
            i++;
        }
        return tmp == ru;
    }

    /**
     * If two repeat models are compatible.
     * 0 - Equivalent; 3 - Incompatible;
     * 1 - Compatible, non-inclusive;
     * 2 - Inclusive.
     * 4 - Inclusive and redundant
     */
    int32_t ru_compare(candidate_fuzzy_motif& rhs) {
        if (ru == rhs.ru && label == rhs.label) {
            return 0; // Equivalent
        }
        if (ru == rhs.ru) {
            return 1;
        }
        std::set<char> s1, s2;
        std::for_each(ru.begin(), ru.end(), [&s1] (char c) -> void { s1.insert(c);});
        std::for_each(rhs.ru.begin(), rhs.ru.end(), [&s2] (char c) -> void { s2.insert(c);});
        if (s1 != s2) { // differnet base composition
            return 3;
        }
        uint32_t k1 = ru.length();
        uint32_t k2 = rhs.ru.length();
        if (inexact == rhs.inexact) {
            // Because of canonical form, only need to check inclusion
            if (k1 == k2) {
                return 3;
            }
            if (k1 < k2) {
                return rhs.contains(*this) ? 4 : 3;
            }
            return contains(rhs) ? 4 : 3;
        }
        if (explains(rhs) || rhs.explains(*this)) {
            return 1;
        }
        return 3;
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
        if (variable_base.size() == 0) {
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

    void reduce() {
        if (ru.length() == 1) {
            return;
        }
        if (!inexact) {
            for (uint32_t k = 1; k <= ru.length()/2; ++k) {
                std::string subunit = ru.substr(0,k);
                uint32_t j = k;
                bool flag = 1;
                while(j < ru.length()) {
                    if (j+k > ru.length()) {
                        if (ru.substr(j) + ru.substr(0, j+k-ru.length()) != subunit) {
                            flag = 0;
                            break;
                        }
                    } else if (ru.substr(j, k) != subunit) {
                        flag = 0;
                        break;
                    }
                    j += k;
                }
                if (flag) {
                    ru = subunit;
                    break;
                }
            }
            return;
        }
        if (ru.length() < 4) {
            return;
        }
        // TODO: the inexact part needs rewrite
        std::string vb;
        for (uint32_t i = 0; i < ru.length(); ++i) {
            if (variable_base[i].first == variable_base[i].second) {
                vb += '0';
            } else {
                vb += '1';
            }
        }
        for (uint32_t k = 2; k <= ru.length()/2; ++k) {
            if (ru.length() % k == 0) {
                std::string subunit = ru.substr(0,k);
                std::string subvarb = vb.substr(0,k);
                uint32_t j = k;
                bool flag = 1;
                while(j < ru.length()) {
                    if (ru.substr(j, k) != subunit || vb.substr(j, k) != subvarb) {
                        flag = 0;
                        break;
                    }
                    j += k;
                }
                if (flag) {
                    ru = subunit;
                    variable_base = std::vector< std::pair<int32_t, int32_t> >(variable_base.begin(), variable_base.begin()+k);
                    return;
                }
            }
        }
        return;
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
