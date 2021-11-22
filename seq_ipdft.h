#ifndef BINARY_PERIOD_H
#define BINARY_PERIOD_H

#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <string>
#include <cmath>
#include <cfloat>
#include <iomanip>
#include <complex>
#include <utility>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <stdio.h>
#include <string.h>

/**
 * Recording periodicity eveidence for each base
 */
struct base_evidence {
    char b;
    int32_t n;
    std::vector<double> ev; // energy for each unit length (period)
    std::vector<double> av; // consisitency after shifting by a unit
    std::vector<double> mv; // multiply the above two
    std::vector<size_t> sorted_index;

    void sort_indexes(std::vector<double> & v) {
      std::stable_sort(sorted_index.begin(), sorted_index.end(),
           [&v](size_t i1, size_t i2) {return v[i1] > v[i2];});
    }

    int32_t organize() {
        if (ev.size() < 2) {
            return 0;
        }
        int32_t max_p = ev.size()-1;
        mv.resize(ev.size());
        sorted_index.resize(ev.size());
        double sum_e = 0.;
        for (int32_t i = 2; i <= max_p; ++i) {
            sum_e += ev[i];
        }
        double sum_m = 0.;
        for (int32_t i = 0; i <= max_p; ++i) {
            sorted_index[i] = i;
            ev[i] /= sum_e;
            mv[i] = ev[i] * av[i];
            sum_m += mv[i];
        }
        mv[0] = 0.; mv[1] = 0.;
        for (int32_t i = 0; i <= max_p; ++i) {
            mv[i] /= sum_m;
        }
        sort_indexes(mv);
        return sorted_index[0];
    }

    bool operator < (const base_evidence& rhs) const {
        return (mv[sorted_index[0]] > rhs.mv[rhs.sorted_index[0]]);
    }

};

struct more_frequent
{
    inline bool operator() (const base_evidence& lhs, const base_evidence& rhs)
    {
        return (lhs.n > rhs.n);
    }
};

class periodic_seq {

public:

    bool debug;
    bool clps;
    std::string seq;
    size_t N;
    int32_t max_p;
    std::vector<base_evidence> bases;
    std::map<char, int32_t> base_ct;
    std::vector<char> alphabet{'A','C','G','T'};
    std::vector<std::string> candidate;

    periodic_seq(const char* _s, int32_t _p, bool _d = 0, bool _c = 0) : max_p(_p), debug(_d), clps(_c) {

        for (auto& s : alphabet) {
            base_ct[s] = 0;
        }
        for (size_t i = 0; i < strlen(_s); ++i) {
            auto ptr = base_ct.find(_s[i]);
            if (ptr != base_ct.end()) {
                ptr->second++;
            }
        }
        if (base_ct.size() > 2) {
            for (size_t i = 0; i < 4; ++i) {
                auto ptr = base_ct.find(alphabet[i]);
                if (ptr != base_ct.end()) {
                    if (ptr->second < 3) {
                        base_ct.erase(ptr);
                    }
                }
            }
        }
        seq = "";
        char last_base = ' ';
        for (size_t i = 0; i < strlen(_s); ++i) {
            auto ptr = base_ct.find(_s[i]);
            if (ptr == base_ct.end()) {
                continue;
            }
            if (clps && _s[i] == last_base) {
                continue;
            }
            seq += _s[i];
            last_base = _s[i];
        }
        N = seq.length();
        max_p = std::min((int32_t) (N/2), max_p);
// if (debug) {
//     printf("%lu\t%lu\t%s\t%lu\n", N, base_ct.size(), seq.c_str(), seq.length());
// }
        for (const auto& s : base_ct) {
            if (s.second > 2) {
                base_evidence elmt;
                elmt.b = s.first;
                elmt.n = s.second;
                std::vector<bool> bseq(N, 0);
                for (size_t n = 0; n < N; ++n) {
                    if (seq.at(n) == s.first) {
                        bseq[n] = 1;
                    }
                }
                ipdft(bseq, elmt.ev);
                shiftc(bseq, elmt.av, elmt.b);
                int32_t picked_p = elmt.organize();
                if (picked_p > 0) {
                    bases.push_back(elmt);
                }
// if (debug) {
//     std::cerr << elmt.b << '\t' << elmt.n << '\t' << picked_p << '\n';
//     for (auto v : bseq) {
//         std::cerr << v << ' ';
//     }
//     std::cerr << '\n';
//     for (auto &v : elmt.mv) {
//         std::cerr << v << ' ';
//     }
//     std::cerr << '\n';
// }
            }
        }

    }
    ~periodic_seq() {}

    /**
     * Get candidate repeat units
     * Parameters: threshold to for calling a base as periodic
     */
    void get_candidate(double min_ecover = 0.0, double min_pcover = 0.0) {

        if (bases.size() < 2) {
            return;
        }

// if (debug) {
//     for (auto &v : bases) {
//         std::cerr << v.b << '\t' << v.n << '\t' << v.sorted_index[0] << '\t' << v.ev[v.sorted_index[0]] << '\n';
//     }
// }
        if (bases.size() > 1 && min_ecover > 0. && min_pcover > 0.) {
            auto it = bases.begin();
            while (it != bases.end()) {
                if (it->ev[it->sorted_index[0]] < min_ecover || it->av[it->sorted_index[0]] < min_pcover) {
                    it = bases.erase(it);
                } else {
                    ++it;
                }
            }
        }
        if (bases.size() < 2) {
            return;
        }

        std::sort(bases.begin(), bases.end(), more_frequent());
        // std::sort(bases.begin(), bases.end());
        for ( uint32_t mi = bases.size(); mi > 1; --mi ) {
            int32_t lcm = bases[0].sorted_index[0]; // least common multiple
            std::set<char> kept_base;
            kept_base.insert(bases[0].b);
            for (uint32_t m = 1; m < mi; ++m) {
                int32_t p0 = bases[m].sorted_index[0];
                if (p0 >= lcm) {
                    lcm = lcm*p0/gcd(p0,lcm);
                } else {
                    lcm = lcm*p0/gcd(lcm,p0);
                }
                if (lcm > N/2) {
                    break;
                }
                kept_base.insert(bases[m].b);
            }
            if (kept_base.size() <= 1) {
                return;
            }
            // remove bases that are less periodic - de-noising
            // TODO: may not be the right approach since periodic base identification is not perfect
            std::string subseq = "";
            for (size_t i = 0; i < N; ++i) {
                if (kept_base.find(seq[i]) != kept_base.end()) {
                    subseq+=seq[i];
                }
            }
            if (lcm > subseq.size()/2) {
                return;
            }
            // count # of exact occurence of each kmer and pos of first appearence
            std::map<std::string, std::pair<int32_t, int32_t> > kmer;
            for (size_t i = 0; i < subseq.size()-lcm; ++i) {
                auto ptr = kmer.find(subseq.substr(i,lcm));
                if (ptr != kmer.end()) {
                    ptr->second.first++;
                } else {
                    bool homo = 1;
                    uint32_t j = 0;
                    while (j < lcm) {
                        homo = homo && (subseq.at(i+j) == subseq.at(i));
                        j++;
                    }
                    if (homo) {continue;}
                    kmer[subseq.substr(i,lcm)] = std::make_pair(1,i);
                }
            }
            // pick kmer that occurs most often, break tie by who comes first
            std::string candi = "";
            int32_t max_ct = 0;
            int32_t min_pt = subseq.size();
            for (const auto& ptr: kmer) {
                if (ptr.second.first > max_ct) {
                    candi = ptr.first;
                    max_ct = ptr.second.first;
                    min_pt = ptr.second.second;
                }
            }
            if (max_ct > 1) {
                for (const auto& ptr: kmer) {
                    if (ptr.second.first == max_ct) {
                        candidate.push_back(ptr.first);
                    }
                }
            }
        }


    }

    int32_t gcd(int32_t a, int32_t b) {
        if (b==0) {return a;}
        return gcd(b, a%b);
    }

    /**
     * Integer period discrete Fourier transform (max_p x seq.size())
     */
    void ipdft(std::vector<bool>& seq, std::vector<double>& evec) {

        const double pi = std::acos(-1);
        const std::complex<double> i(0, 1);

        evec.resize(max_p+1);
        std::vector< std::vector<int32_t> > acc_list;
        for (int32_t p = 1; p < max_p+1; ++p) {
            std::vector<int32_t> v(p, 0);
            acc_list.push_back(v);
            evec[p] = 0;
        }
        for (int32_t n = 0; n < N; ++n) {
            if (seq[n]) {
                for (int32_t p = 1; p < max_p+1; ++p) {
                    int32_t l = n % p;
                    acc_list[p-1][l] += 1;
                }
            }
        }
        for (int32_t p = 1; p < max_p+1; ++p) {
            std::complex<double> b(0.,0.);
            for (int32_t l = 0; l < p; ++l) {
                b += ((double) acc_list[p-1][l]) * std::exp( - i * 2. * pi * (double) (l * 1. / p) );
            }
            evec[p] = std::pow(std::abs(b), 2);
        }
    }

    /**
     * Shift invariant (max_p x seq.size())
     */
    void shiftc(std::vector<bool>& seq, std::vector<double>& evec, char& b) {

        evec.resize(max_p+1);
        evec[0] = 0.;
        int32_t tot = base_ct[b];
        for (int32_t p = 1; p < max_p+1; ++p) {
            evec[p] = 0.;
            for (int32_t i = 0; i < N; ++i) {
                evec[p] += seq[i] && seq[(i+p) % N];
            }
            evec[p] /= tot;
        }
    }

};

#endif
