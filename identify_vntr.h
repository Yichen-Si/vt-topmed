#ifndef IDENTIFY_VNTR_H
#define IDENTIFY_VNTR_H

#include <iterator>
#include "htslib/kstring.h"
#include "bcf_ordered_reader.h"
#include "bcf_ordered_writer.h"
#include "variant_manip.h"
#include "program.h"
#include "filter.h"
#include "vntr_annotator.h"

#define BUFFER_BP 2000

void identify_vntr(int argc, char ** argv);

#endif
