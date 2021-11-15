/* The MIT License

   Copyright (c) 2014 Adrian Tan <atks@umich.edu>

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.
*/

#include "identify_vntr.h"

namespace
{

class Igor : Program
{
    public:

    std::string version;

    ///////////
    //options//
    ///////////
    std::string input_vcf_file;
    std::string ref_fasta_file;
    std::string output_vcf_file;
    std::vector<GenomeInterval> intervals;
    std::string interval_list;

    std::string method;           //methods of detection
    std::string annotation_mode;  //modes of annotation
    int32_t vntr_classification;  //vntr classification schemas
    bool override_tag;
    int32_t flk_size;
    int32_t max_mlen;

    //FILTER ids
    // int32_t overlap_homopolymer; // later.

    //INFO tags
    std::string TR;               // chrom:start:end:RU
    std::string RU;               // repeat unit
    std::string RL;               // repeat region length
    std::string CONCORDANCE;   // concordance of the repeat unit
    std::string RU_COUNTS;        // repeat unit counts - in longest allele and in ref
    std::string FLANKS;           // flank positions
    std::string FLANKSEQ;
    std::string VITERBI_SCORE;
    std::string LONGESTMOSAIC;

    //helper variables for populating additional VNTR records
    uint32_t no_samples;
    int32_t* gts;

    //convenience kstring so that we do not have to free memory ALL the time.
    kstring_t s;

    bool debug;

    /////////////
    //vntr buffer
    /////////////
    std::list<VNTR_candidate> vntr_buffer; // sorted by end position, front is the most recent (largest endding point)

    //////////
    //filter//
    //////////
    std::vector<int32_t> req_flt_ids;
    std::vector<std::string> required_filters;
    bool filter_exists;

    ///////
    //i/o//
    ///////
    BCFOrderedReader *odr;
    BCFOrderedWriter *odw;

    /////////
    //stats//
    /////////
    int32_t n_indels, n_vntr, n_vntr_out;

    ////////////////
    //common tools//
    ////////////////
    VariantManip* vm;
    VNTRAnnotator* va;
    faidx_t* fai;

    Igor(int argc, char **argv) {
        version = "0.5";
        //////////////////////////
        //options initialization//
        //////////////////////////
        try {
            std::string desc = "annotates indels with VNTR information - repeat tract length, repeat motif, flank information";

            TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file [-]", false, "-", "str", cmd);
            TCLAP::ValueArg<std::string> arg_ref_fasta_file("r", "r", "reference sequence fasta file []", true, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_annotation_mode("a", "a", "annotation type [v]\n"
                 "              v : a. output VNTR variant (defined by classification).\n"
                 "                     RU    repeat unit on reference sequence (CA)\n"
                 "                     MOTIF canonical representation (AC)\n"
                 "                     RL    repeat tract length in bases (11)\n",
                 false, "v", "str", cmd);
            TCLAP::ValueArg<int32_t> arg_vntr_classification("c", "c", "classification schemas of tandem repeat [6]\n"
                 "              0 : viterbi_score\n"
                 "              1 : lai2003      \n"
                 "              2 : kelkar2008   \n"
                 "              3 : fondon2012   \n"
                 "              4 : ananda2013   \n"
                 "              5 : willems2014  \n"
                 "              6 : tan_kang2015",
                 false, 0, "integer", cmd);
            TCLAP::ValueArg<std::string> arg_method("m", "m", "mode [f]\n"
                 "              e : by exact alignment"
                 "              f : by fuzzy alignment",
                 false, "f", "str", cmd);

            TCLAP::SwitchArg arg_debug("d", "d", "debug [false]", cmd, false);
            TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "","file", cmd);
            TCLAP::MultiArg<std::string> arg_required_filters("f", "apply-filter", "Required filters", false, "str", cmd);
            TCLAP::SwitchArg arg_override_tag("x", "x", "override tags [false]", cmd, false);
            TCLAP::SwitchArg arg_add_vntr_record("v", "v", "add vntr record [false]", cmd, false);

            TCLAP::ValueArg<int32_t> arg_flk_size("k", "padding_size", "bp flanking a TR on each side to annotate", false, 10, "integer", cmd);
            TCLAP::ValueArg<int32_t> arg_max_mlen("n", "max-mlen", "Maximum repeat unit length", false, 16, "integer", cmd);

            cmd.parse(argc, argv);

            input_vcf_file = arg_input_vcf_file.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());
            method = arg_method.getValue();
            annotation_mode = arg_annotation_mode.getValue();
            vntr_classification = arg_vntr_classification.getValue();
            override_tag = arg_override_tag.getValue();
            required_filters = arg_required_filters.getValue();
            debug = arg_debug.getValue();
            ref_fasta_file = arg_ref_fasta_file.getValue();

            flk_size = arg_flk_size.getValue();
            max_mlen = arg_max_mlen.getValue();
        }
        catch (TCLAP::ArgException &e) {
            std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
            abort();
        }
    };

    ~Igor() {
        if (s.m) free(s.s);
    };

    void initialize() {
        ///////////
        //options//
        ///////////
        if (method!="e" && method!="f") {
            fprintf(stderr, "[%s:%d %s] Not a valid mode of VNTR detection: %s\n", __FILE__,__LINE__,__FUNCTION__, method.c_str());
            exit(1);
        }

        if (annotation_mode!="v") {
            fprintf(stderr, "[%s:%d %s] Not a valid mode of annotation: %s\n", __FILE__,__LINE__,__FUNCTION__, annotation_mode.c_str());
            exit(1);
        }

        /////////////////////////
        //filter initialization//
        /////////////////////////
        filter_exists = 0;
        if ( required_filters.size() > 0 ) {
          for(uint32_t i=0; i < required_filters.size(); ++i) {
            req_flt_ids.push_back(bcf_hdr_id2int(odr->hdr, BCF_DT_ID, required_filters[i].c_str()));
          }
          filter_exists = 1;
        }

        //////////////////////
        //i/o initialization//
        //////////////////////
        odr = new BCFOrderedReader(input_vcf_file, intervals);
        odw = new BCFOrderedWriter(output_vcf_file, 10000);
        odw->link_hdr(odr->hdr);

        //////////////////////////////
        //INFO header for new VNTR records in VCF
        //////////////////////////////
        // bcf_hdr_append(odw->hdr, "##FILTER=<ID=overlap_homopolymer,Description=\"Overlaps with homopolymer\">");
        // OVERLAP_HOMOPOLYMER = bcf_hdr_id2int(odw->hdr, BCF_DT_ID, "overlap_homopolymer");
        TR = bcf_hdr_append_info_with_backup_naming(odw->hdr, "TR", "1", "String", "Tandem repeat unique ID.", false);
        RU = bcf_hdr_append_info_with_backup_naming(odw->hdr, "RU", ".", "String", "Repeat unit in a VNTR or homopolymer", false);
        RL = bcf_hdr_append_info_with_backup_naming(odw->hdr, "RL", ".", "Float", "Repeat length (reference and the longest)", false);
        VITERBI_SCORE = bcf_hdr_append_info_with_backup_naming(odw->hdr, "VITERBI_SCORE", ".", "Float", "Viterbi score repeat region.", false);
        CONCORDANCE = bcf_hdr_append_info_with_backup_naming(odw->hdr, "CONCORDANCE", ".", "Float", "Concordance of repeat unit.", false);
        RU_COUNTS = bcf_hdr_append_info_with_backup_naming(odw->hdr, "RU_COUNTS", ".", "Integer", "Number of exact repeat units and total number of repeat units.", false);
        FLANKS = bcf_hdr_append_info_with_backup_naming(odw->hdr, "FLANKS", "2", "Integer", "Left and right most position of the repeat region.", false);
        FLANKSEQ = bcf_hdr_append_info_with_backup_naming(odw->hdr, "FLANKSEQ", "1", "String", "Flanking sequence 10bp on either side of detected repeat region.", false);
        LONGESTMOSAIC = bcf_hdr_append_info_with_backup_naming(odw->hdr, "LONGESTMOSAIC", "1", "String", "Longest allele.", false);

        bcf_hdr_append(odw->hdr, "##INFO=<ID=LARGE_REPEAT_REGION,Number=0,Type=Flag,Description=\"Very large repeat region, vt only detects up to 1000bp long regions.\">");

        //helper variable initialization for adding genotype fields for additional vntr records
        if (annotation_mode=="v")
        {
            no_samples = bcf_hdr_nsamples(odw->hdr);
            gts = (int32_t*) malloc(no_samples*sizeof(int32_t));
            for (uint32_t i=0; i<no_samples; ++i)
            {
                gts[i] = 0;
            }
        }
        else
        {
            no_samples = 0;
            gts = NULL;
        }

        s = {0,0,0};

        ////////////////////////
        //stats initialization//
        ////////////////////////
        n_indels = 0;
        n_vntr = 0;
        n_vntr_out = 0;

        ////////////////////////
        //tools initialization//
        ////////////////////////
        vm = new VariantManip(ref_fasta_file);
        va = new VNTRAnnotator(ref_fasta_file, debug, max_mlen);
        fai = fai_load(ref_fasta_file.c_str());
    }

    void print_options()
    {
        std::clog << "identify_vntr v" << version << "\n";
        std::clog << "\n";
        std::clog << "options:     input VCF file(s)        " << input_vcf_file << "\n";
        std::clog << "         [o] output VCF file          " << output_vcf_file << "\n";
        std::clog << "         [m] method of VNTR detection " << method << "\n";
        std::clog << "         [a] mode of annotation       " << annotation_mode << "\n";
        std::clog << "         [k] max repeat unit length   " << max_mlen << "\n";
        print_boo_op("         [d] debug                    ", debug);
        print_ref_op("         [r] ref FASTA file           ", ref_fasta_file);
        print_boo_op("         [x] override tag             ", override_tag);
        print_int_op("         [i] intervals                ", intervals);
        std::clog << "\n";
    }

    void print_stats()
    {
        std::clog << "\n";
        std::cerr << "Number of indels                     " << n_indels   << "\n";
        std::cerr << "Number of indels identified as vntr  " << n_vntr     << "\n";
        std::cerr << "Number of (merged) vntr output       " << n_vntr_out << "\n";
        std::clog << "\n";
    }

    /**
     * Inserts a VNTR candidate.
     * Returns true if successful.
     */
    void insert_vntr_into_buffer(VNTR_candidate& vntr) {

        if (vntr.need_refit) {
            vntr.model_refit();
        }
        auto it = vntr_buffer.begin();
        // Find overlap
        bool merged = 0;
        while (it != vntr_buffer.end()) {
            VNTR_candidate& cvntr = *it;
            int32_t merge_type = cvntr.merge(vntr);
if (debug && merge_type > 0) {
    std::cerr << "Merge_type " << merge_type << '\n';
}
            if (it != vntr_buffer.begin() && cvntr.ed > std::prev(it)->ed) {
                // End point of existing vntr changes, need to restore order
                auto itt = vntr_buffer.begin();
                auto tmp = std::next(it);
                while (itt->ed > cvntr.ed && itt != vntr_buffer.end()) {
                    itt++;
                }
                vntr_buffer.splice(itt, vntr_buffer, it);
                it = std::prev(tmp);
            }
            if(merge_type == 1) {
                return;
            }
            if (merge_type != 0) {
                merged = 1;
            }
            if (cvntr.motif.ed < vntr.motif.st - vntr.rlen_l*vntr.max_interrupt) {
                break;
            }
            it++;
        }
        if (merged) {
            return;
        }
        // If not merged, insert
        it = vntr_buffer.begin();
        while (it != vntr_buffer.end() && it->ed > vntr.ed) {
            it++;
        }
        vntr_buffer.insert(it, vntr);
if (debug) {
printf("------------- insert_vntr_into_buffer: inserted new: st %d, ed %d, ru %s\n", vntr.motif.st, vntr.motif.ed, vntr.motif.ru.c_str());
}

    }

    /**
     * Try to output old VNTR records
     */
    void flush_vntr_buffer(int32_t rid, int32_t pos) {
        // std::cerr << "Start output vntr\n";
        if (vntr_buffer.empty()) {
            return;
        }
        // int32_t rid = v->rid;
        // int32_t pos = v->pos;
        //search for vntr to start deleting from.
        auto it = vntr_buffer.end();
        if (std::prev(it)->ed > pos-BUFFER_BP && std::prev(it)->rid == rid) {
            return; // The most distance record is indispensible
        }
        while(it != vntr_buffer.begin()) {
            VNTR_candidate& vntr = *std::prev(it);
            if (vntr.rid > rid) {
                //rid < vntr.rid is impossible
                fprintf(stderr, "[%s:%d %s] flush_vntr_buffer - File %s is unordered\n", __FILE__, __LINE__, __FUNCTION__, input_vcf_file.c_str());
                exit(1);
            }
            if (vntr.ed > pos-BUFFER_BP && vntr.rid == rid) {
                break;
            }
            --it;
        }
        std::list<VNTR_candidate> torm(it, vntr_buffer.end());
        torm.sort(); // Sort according to start point
if (debug) {
    std::cerr << "------------- flush_vntr_buffer before " << pos-BUFFER_BP << '\t' << vntr_buffer.size() << '\t' << torm.size() << '\n';
    for (auto &v : torm) {
        printf("%d, %d, %s\n", v.st, v.ed, v.motif.ru.c_str());
    }

}
        // Need to check if they could be further merged
        auto pt = torm.begin();
        while (pt != torm.end()) {
            VNTR_candidate& vntr = *pt;
            auto itt = std::next(pt);
            bool merged = 0;
            while (itt != torm.end()) {
                VNTR_candidate& cvntr = *itt;
                int32_t merge_type = cvntr.merge(vntr);
                if (merge_type != 0) {
                    break;
                }
                itt++;
            }
            if (!merged) { // Output as an independent record
                write_vntr_to_vcf(vntr);
            }
            pt++;
        }
        while (it != vntr_buffer.end()) {
            it = vntr_buffer.erase(it);
        }
    }

    void flush_vntr_buffer() {
        if (vntr_buffer.size() == 0) {
            return;
        }
if (debug) {
    std::cerr << "------------- flush_vntr_buffer output entire buffer " << vntr_buffer.size() << "\n";
}
        int32_t r = vntr_buffer.front().rid;
        int32_t p = vntr_buffer.front().ed + BUFFER_BP*10;
        flush_vntr_buffer(r, p);
    }

    /**
     * Write a VNTR record to be written to vcf.
     */
    void write_vntr_to_vcf(VNTR_candidate& vntr) {

        bcf_hdr_t* h = odw->hdr;
        if (!vntr.finalize()) {
            return;
        }
        bcf1_t* v = bcf_init1();
        v->rid = vntr.rid;
        v->pos = vntr.st;

        std::string tname = vntr.repeat_ref + ",<VNTR>";
        bcf_update_alleles_str(h, v, tname.c_str());

        std::stringstream ss;
        ss << vntr.motif.ru << ':';
        for (uint32_t i = 0; i < vntr.motif.mlen; ++i) {
            ss << vntr.motif.label[i];
        }
        std::string ru_id = ss.str();
        bcf_update_info_string(h, v, RU.c_str(), ss.str().c_str());

        std::string chrom(vntr.chrom);
        ss.clear();
        ss.str("");
        ss << chrom << ':' << vntr.st << ':' << vntr.ed << ':' << ru_id;
        bcf_update_info_string(h, v, TR.c_str(), ss.str().c_str());

        int32_t rlen[2] = {vntr.rlen_r, vntr.rlen_l};
        bcf_update_info_int32(h, v, RL.c_str(), &rlen, 2);

        int32_t ru_count[2] = {vntr.n_ru_r, vntr.n_ru_l};
        bcf_update_info_int32(h, v, RU_COUNTS.c_str(), &ru_count, 2);

        float scores[2] = {(float) vntr.score_r, (float) vntr.score_l};
        bcf_update_info_float(h, v, VITERBI_SCORE.c_str(), &scores, 2);

        int32_t flks[2] = {vntr.st, vntr.ed};
        bcf_update_info_int32(h, v, FLANKS.c_str(), &flks, 2);

        std::string flankseq = vntr.lflank + '[' + vntr.repeat_ref + ']' + vntr.rflank;
        bcf_update_info_string(h, v, FLANKSEQ.c_str(), flankseq.c_str());

        if (vntr.len_mrg > vntr.repeat_ref.length()) {
        bcf_update_info_string(h, v, LONGESTMOSAIC.c_str(), vntr.query.substr(vntr.rel_st, vntr.len_mrg).c_str());
        }

        // if (vntr.is_large_repeat) {
        //     bcf_update_info_flag(h, v, "LARGE_REPEAT_REGION", NULL, 1);
        // }

        //individual fields - just set GT
        bcf_update_genotypes(h, v, gts, no_samples);
        odw->write(v);
        n_vntr_out++;

// std::cerr << "write_vntr_to_vcf: " << ss.str() << '\n';

    }

    bool has_filter(bcf1_t* v) {
        bool flag = !filter_exists;
        if ( !flag ) {
            for(int32_t i=0; i < v->d.n_flt; ++i) {
                for(uint32_t j=0; j < req_flt_ids.size(); ++j) {
                    if ( req_flt_ids[j] == v->d.flt[i] )
                    flag = true;
                }
            }
        }
        return flag;
    }

    void identify_vntr()
    {
        odw->write_hdr();
        bcf1_t* v = bcf_init1();
        // bcf1_t *v = odw->get_bcf1_from_pool();
        bcf_hdr_t *h = odw->hdr;
        const char* chrom = bcf_get_chrom(h, v);
        // Variant variant;
        int32_t no_exact = 0;
        int32_t no_inexact = 0;

        while (odr->read(v)) {
            bcf_unpack(v, BCF_UN_STR);
            if (!has_filter(v)) {continue;}

            // // Normalize
            // if (!vm->is_normalized(v)) {
            //     uint32_t pos1 = bcf_get_pos1(v);
            //     std::vector<std::string> alleles;
            //     for (size_t i=0; i<bcf_get_n_allele(v); ++i) {
            //         char *s = bcf_get_alt(v, i);
            //         while (*s) {
            //             *s = toupper(*s);
            //             ++s;
            //         }
            //         alleles.push_back(std::string(bcf_get_alt(v, i)));
            //     }
            //     uint32_t left_extended = 0, left_trimmed = 0, right_trimmed = 0;
            //     vm->right_trim_or_left_extend(alleles, pos1, bcf_get_chrom(h, v), left_extended, right_trimmed);
            //     vm->left_trim(alleles, pos1, left_trimmed);
            //     if (left_trimmed || left_extended || right_trimmed) {
            //         kstring_t new_alleles = {0,0,0};
            //         bcf_set_pos1(v, pos1);
            //         new_alleles.l = 0;
            //         for (size_t i=0; i<alleles.size(); ++i) {
            //             if (i) kputc(',', &new_alleles);
            //             kputs(alleles[i].c_str(), &new_alleles);
            //         }
            //         bcf_update_alleles_str(odw->hdr, v, new_alleles.s);
            //     }
            // }
            // int32_t vtype = vm->classify_variant(odr->hdr, v, variant);

            if (bcf_get_variant_types(v) == VCF_INDEL) {
                n_indels++;
                std::vector<candidate_fuzzy_motif> candidate_model;
                if (va->annotate(odr->hdr, v, candidate_model) > 0) {
                    n_vntr++;
                    for (auto & vntr_model : candidate_model) {
                        VNTR_candidate vntr(fai, odr->hdr, v, vntr_model, debug);
                        insert_vntr_into_buffer(vntr);
                    }
                }
            }
            flush_vntr_buffer(v->rid, v->pos);
        }
        flush_vntr_buffer();
        odw->close();
        odr->close();
    }
};
}

void identify_vntr(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.identify_vntr();
    igor.print_stats();
};
