#ifndef SHIFT_ALIGN_H
#define SHIFT_ALIGN_H

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
#include "hts_utils.h"
#include "bcf_ordered_reader.h"
#include "bcf_ordered_writer.h"
#include "program.h"

class ShiftAlign
{
public:

     faidx_t* fai;
     int32_t bound;
     int32_t max_l;
     bool debug;
     const char* chrom;
     int32_t p_st, p_ed, s_st;
     std::string extension;
     float s0, s1, g0;

     ShiftAlign(faidx_t* _f, bcf_hdr_t* _h, bcf1_t* _v, int32_t _m=16, int32_t _l=64, bool _d=0);
     ~ShiftAlign() {}

     void right_shift(int32_t& st0, int32_t& st, int32_t& ed, std::string& mseq);
     void left_shift(int32_t& ed0, int32_t& st, int32_t& ed, std::string& mseq);
     void shift(std::string& seq, int32_t& min_i, int32_t& min_j, int32_t& max_i, int32_t& max_j, std::string& mseq);
};







#endif
