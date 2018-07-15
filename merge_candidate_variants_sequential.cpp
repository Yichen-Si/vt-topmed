/* The MIT License

   Copyright (c) 2015 Adrian Tan and Hyun Min Kang <atks@umich.edu>

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

#include "merge_candidate_variants_sequential.h"

namespace
{
  struct bcf1_comp {
    bool operator() (const bcf1_t* lhs, const bcf1_t* rhs) const {
      if ( lhs->rid == rhs->rid ) {
	if ( lhs->pos == rhs->pos ) {
	  if ( lhs->rlen == rhs->rlen ) {
	    if ( lhs->n_allele == rhs->n_allele ) {
	      bcf_unpack((bcf1_t*)lhs, BCF_UN_STR);
	      bcf_unpack((bcf1_t*)rhs, BCF_UN_STR);
	      for(int32_t i=0; i < lhs->n_allele; ++i) {
		int32_t cmp = strcmp(lhs->d.allele[i], rhs->d.allele[i]);
		if ( cmp != 0 ) return ( cmp < 0 );
	      }
	      return false;
	    }
	    else return (lhs->n_allele < rhs->n_allele);
	  }
	  else return (lhs->rlen < rhs->rlen);		  
	}
	else return (lhs->pos < rhs->pos);	
      }
      else return (lhs->rid < rhs->rid);
    }
  };

  static std::vector<double> logfacs;
    
  //binom() { logfacs.push_back(0); }
  static double log_d(int32_t x, int32_t n, double p) {
    if ( x > n ) x = n;
    if ( ( x < 0 ) || ( n < 0 ) || ( x > n ) || ( p <= 0 ) || ( p >= 1 ) ) {
      fprintf(stderr,"[E:%s:%d %s] Cannot calculate binomial density with (%d, %d, %lf)\n", __FILE__, __LINE__, __FUNCTION__, x, n, p);
      exit(1);
    }
    for(int32_t i=(int32_t)logfacs.size(); i <= n; ++i) {
      if ( i == 0 ) logfacs.push_back(0);
      else logfacs.push_back(logfacs[i-1] + log((double)i));
    }
    // n!/(x!(x-n)!)
    //fprintf(stderr, "%d %d %lf %lf\n", x, n, p, logfacs[n] - logfacs[x] - logfacs[n-x] + x * log(p) + (n-x) * log(1.0-p) );
    return ( logfacs[n] - logfacs[x] - logfacs[n-x] + x * log(p) + (n-x) * log(1.0-p) );
  }

  class varMergeInfo {
  public:
    int32_t nsamples;
    //int32_t nhets;
    //int32_t nalts;
    
    //double sumabe;
    //double sqabe;
    //double sumabz;
    //double sqabz;
    float maxqual;

    std::vector<int32_t> ids;
    int32_t ab20[20];
    int32_t dp20[20];    

    //varMergeInfo() : nsamples(0), nhets(0), nalts(0), sumabe(0), sqabe(0), sumabz(0), sqabz(0), maxqual(0) {}
    varMergeInfo() : nsamples(0), maxqual(0) {
      memset(ab20, 0, sizeof(int32_t)*20);
      memset(dp20, 0, sizeof(int32_t)*20);      
    }    
    void add(int32_t id, int32_t e, int32_t n, float qual) {
      if ( n == 0 ) return;

      if ( e > n ) e = n;

      // assuming binomial distribution, compute likelihoods
      double p0 = log_d(e, n, 0.01);
      double p1 = log_d(e, n, 0.50); 
      //double p2 = log_d(e, n, 0.99);

      if ( p0 > p1 ) return;

      /*
      if ( p1 > p0 ) {
	if ( p1 < p2 ) {
	  ++nalts;
	}
	else {
	  ++nhets;
	  sumabe += ((n - e + 0.5)/(n + 1.0));
	  sqabe += ((n - e + 0.5)*(n - e + 0.5)/(n + 1.0)/(n + 1.0));
	  double abz = ((n - e + 1e-9)/(double)(n+2e-9) - 0.5) / sqrt( n*0.5*0.5 + 1e-9);
	  sumabz += abz;
	  sqabz += (abz*abz);
	}
      }
      else if ( p2 > p0 ) {
	  ++nalts;
      }
      else {
	return;
      }
      */

      if ( qual > maxqual ) maxqual = qual;      
      if ( nsamples < 10 ) // store first 10 individuals
	ids.push_back(id);
      ++ab20[(int32_t)(e/(n+1e-6)*20.0)];
      ++dp20[n > 99 ? 19 : n/5];
      ++nsamples;
    }
  };

class Igor : Program
{
public:
  
  std::string version;
  
  ///////////
  //options//
  ///////////
  std::vector<std::string> input_vcf_files;
  std::string input_vcf_file_list;
  std::string output_vcf_file;
  std::vector<GenomeInterval> intervals;
  std::string interval_list;
  float snp_variant_score_cutoff;
  float indel_variant_score_cutoff;
  
  std::map<bcf1_t*, varMergeInfo, struct bcf1_comp> variants;
  ///////
  //i/o//
  ///////
  BCFOrderedReader *odr;
  BCFOrderedWriter *odw;
  bcf1_t *v;

  ///////////////
  //general use//
  ///////////////
  kstring_t variant;
  
  /////////
  //stats//
  /////////
  uint32_t no_samples;
  uint32_t no_candidate_snps;
  uint32_t no_candidate_indels;
  
  /////////
  //tools//
  /////////
  LogTool *lt;
  VariantManip * vm;

  Igor(int argc, char ** argv)
  {
    //////////////////////////
    //options initialization//
    //////////////////////////
    try
      {
	std::string desc =
	  "Merge candidate variants across samples.\n\
Sequentially merge VCF file to reduce the memory consumptions.\n\
Each VCF file is required to have the FORMAT flags E and N and should have exactly one sample.";

	version = "0.1";
	TCLAP::CmdLine cmd(desc, ' ', version);
            VTOutput my; cmd.setOutput(&my);
            TCLAP::ValueArg<std::string> arg_intervals("i", "i", "intervals", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_interval_list("I", "I", "file containing list of intervals []", false, "", "str", cmd);
            TCLAP::ValueArg<std::string> arg_output_vcf_file("o", "o", "output VCF file [-]", false, "-", "", cmd);
            TCLAP::ValueArg<std::string> arg_input_vcf_file_list("L", "L", "file containing list of input VCF files", true, "", "str", cmd);
            TCLAP::ValueArg<float> arg_snp_variant_score_cutoff("c", "c", "SNP variant score cutoff [30]", false, 30, "float", cmd);
            TCLAP::ValueArg<float> arg_indel_variant_score_cutoff("d", "d", "Indel variant score cutoff [30]", false, 30, "float", cmd);

            cmd.parse(argc, argv);

            input_vcf_file_list = arg_input_vcf_file_list.getValue();
            output_vcf_file = arg_output_vcf_file.getValue();
            snp_variant_score_cutoff = arg_snp_variant_score_cutoff.getValue();
            indel_variant_score_cutoff = arg_indel_variant_score_cutoff.getValue();
            parse_intervals(intervals, arg_interval_list.getValue(), arg_intervals.getValue());

            ///////////////////////
            //parse input VCF files
            ///////////////////////
            htsFile *file = hts_open(input_vcf_file_list.c_str(), "r");
            if (file==NULL)
            {
                std::cerr << "cannot open " << input_vcf_file_list.c_str() << "\n";
                exit(1);
            }
            kstring_t *s = &file->line;
            while (hts_getline(file, KS_SEP_LINE, s) >= 0)
            {
                if (s->s[0]!='#')
                {
                    input_vcf_files.push_back(std::string(s->s));
                }
            }
            hts_close(file);
        }
        catch (TCLAP::ArgException &e)
        {
            std::cerr << "error: " << e.error() << " for arg " << e.argId() << "\n";
            abort();
        }
    };

  void initialize()
  {
    /////////
    //tools//
    /////////
    lt = new LogTool();
    vm = new VariantManip();
  }
  
  void merge_candidate_variants_sequential()
  {
    uint32_t nfiles = input_vcf_files.size();
    kstring_t sample_names = {0,0,0};
    float af = 0;
    uint32_t no = 0;
    
    //obtain sample names
    std::vector<std::string> sm_ids;

    odw = new BCFOrderedWriter(output_vcf_file, 0);
    
    int32_t *E = (int32_t*) malloc(1*sizeof(int32_t)); //NULL; // = NULL;
    int32_t *N = (int32_t*) malloc(1*sizeof(int32_t)); // [2] = {0,0}; // *N = NULL;
    int32_t no_E = 1, no_N = 1;    

    bcf1_t* nv = bcf_init();
    bcf1_t* wv = bcf_init();
    std::map<bcf1_t*, varMergeInfo, bcf1_comp>::iterator it;
      
    
    for (uint32_t i=0; i<nfiles; ++i) {
      //////////////////////
      //i/o initialization//
      //////////////////////
      odr = new BCFOrderedReader(input_vcf_files[i], intervals);
      if ( i == 0 ) {
	bcf_hdr_append(odw->hdr, "##fileformat=VCFv4.2");
	bcf_hdr_transfer_contigs(odr->hdr, odw->hdr);
	bcf_hdr_append(odw->hdr, "##QUAL=Maximum variant score of the alternative allele likelihood ratio: -10 * log10 [P(Non variant)/P(Variant)] amongst all individuals.");
	bcf_hdr_append(odw->hdr, "##INFO=<ID=NSAMPLES,Number=1,Type=Integer,Description=\"Number of samples.\">");
	bcf_hdr_append(odw->hdr, "##INFO=<ID=SAMPLES,Number=.,Type=String,Description=\"Samples with evidence. (up to first 10 samples)\">");
	bcf_hdr_append(odw->hdr, "##INFO=<ID=AB20,Number=20,Type=Integer,Description=\"Histogram of allele balance using E and N (0 to 100 using 20 bins)\">");
	bcf_hdr_append(odw->hdr, "##INFO=<ID=DP20,Number=20,Type=Integer,Description=\"Histogram of depth (0-100 with 20 bins)\">");		
	//bcf_hdr_append(odw->hdr, "##INFO=<ID=HET,Number=1,Type=Integer,Description=\"Number of individuals with apparent heterozygous genotypes\">");
	//bcf_hdr_append(odw->hdr, "##INFO=<ID=ALT,Number=1,Type=Integer,Description=\"Number of indivisuals with apparent non-ref-homozygous genotypes\">");
	//bcf_hdr_append(odw->hdr, "##INFO=<ID=ABE_AVG,Number=1,Type=Float,Description=\"Average allele balance statistics\">");
	//bcf_hdr_append(odw->hdr, "##INFO=<ID=ABZ_AVG,Number=1,Type=Float,Description=\"Average allele balance Z-statistics\">");
	//bcf_hdr_append(odw->hdr, "##INFO=<ID=ABE_SD,Number=1,Type=Float,Description=\"Standard deviation of average allele balance statistics\">");
	//bcf_hdr_append(odw->hdr, "##INFO=<ID=ABZ_SD,Number=1,Type=Float,Description=\"Standard deviation of average allele balance Z-statistics\">");
	//bcf_hdr_append(odw->hdr, "##INFO=<ID=EXHET,Number=1,Type=Float,Description=\"Naive HWE Z score\">");		
	odw->write_hdr();
      }

      sm_ids.push_back(bcf_hdr_get_sample_name(odr->hdr, 0));
      
      while( odr->read(nv) ) {
	if ( ( nv->pos + 1 < intervals[0].start1 ) || ( nv->pos + 1 > intervals[0].end1 ) ) continue;
        //fprintf(stderr,"foo %d %d %d\n",nv->pos,intervals[0].start1,intervals[0].end1);
	
	bcf_unpack(nv, BCF_UN_STR);
	float variant_score = bcf_get_qual(nv);
	
	if (bcf_float_is_missing(variant_score)) {
	  variant_score = 0;
	}

	int32_t vtype;
	vtype = bcf_is_snp(nv) ? VT_SNP : VT_INDEL;

	if ((vtype == VT_SNP && variant_score >= snp_variant_score_cutoff) ||
	    (vtype == VT_INDEL && variant_score >= indel_variant_score_cutoff)) {
	  //max_variant_score_gt_cutoff = true;
	  //fprintf(stderr,"%d\n", nv->pos);
	  //fprintf(stderr,"%d\n",bcf_get_format_int32(odr->hdr, nv, "E", &E, &no_E));
	  //fprintf(stderr,"%d\n",E[0]);
	  //fprintf(stderr,"%d\n",bcf_get_format_int32(odr->hdr, nv, "N", &N, &no_N));
	  //fprintf(stderr,"%d\n",N[0]);	  

	  if (bcf_get_format_int32(odr->hdr, nv, "E", &E, &no_E) < 0 ||
	      bcf_get_format_int32(odr->hdr, nv, "N", &N, &no_N) < 0) {
	    fprintf(stderr,"[E:%s:%d %s] cannot get format values E or N from %s\n", __FILE__, __LINE__, __FUNCTION__, input_vcf_files[i].c_str());
	    exit(1);
	  }
	  it = variants.find(nv);
	  if ( it != variants.end() ) {
	    it->second.add(i, E[0], N[0], nv->qual);
	  }
	  else {
	    variants[nv].add(i, E[0], N[0], nv->qual);
	    nv = bcf_init();
	  }
	  //fprintf(stderr,"foo\n");	  
	}
      }
      
      odr->close();
      delete odr;

      //if ( (i + 1) % 100 == 0 )
      fprintf(stderr,"Finished processing %d input BCF/VCF files. Current variant count is %lu\n", i+1, variants.size());
    }

    fprintf(stderr,"Finished processing all %d input BCF/VCF files. Current variant count is %lu\n", nfiles, variants.size());    
    
    // print all variants;
    for(it = variants.begin(); it != variants.end(); ++it) {
      bcf_clear(wv);
      wv->rid = it->first->rid;
      wv->pos = it->first->pos;
      bcf_update_alleles(odw->hdr, wv, (const char**)it->first->d.allele, it->first->n_allele);
      
      wv->qual = it->second.maxqual;
      bcf_update_info_int32(odw->hdr, wv, "NSAMPLES", &it->second.nsamples, 1);
      
      std::string id = it->second.nsamples > 0 ? sm_ids[it->second.ids[0]] : ".";
      for(int32_t i=1; i < (int32_t)it->second.ids.size() ; ++i) {
	id = id + "," + sm_ids[it->second.ids[i]];
      }
      bcf_update_info_string(odw->hdr, wv, "SAMPLES", id.c_str());

      bcf_update_info_int32(odw->hdr, wv, "AB20", &it->second.ab20, 20);
      bcf_update_info_int32(odw->hdr, wv, "DP20", &it->second.dp20, 20);      
      /*
      bcf_update_info_int32(odw->hdr, wv, "HET", &it->second.nhets, 1);
      bcf_update_info_int32(odw->hdr, wv, "ALT", &it->second.nalts, 1);

      float abe_avg = (float)(( it->second.sumabe + 5e-3 ) / ( it->second.nhets + 0.01 ));
      float abz_avg = (float)(( it->second.sumabz ) / ( it->second.nhets + 0.01 ));
      float abe_sd  = (float)(sqrt(( it->second.sqabe ) / ( it->second.nhets + 1e-6 ) - (it->second.sumabe * it->second.sumabe ) / (it->second.nhets + 1e-6) / (it->second.nhets + 1e-6) + 1e-6));
      float abz_sd  = (float)(sqrt(( it->second.sqabz ) / ( it->second.nhets + 1e-6 ) - (it->second.sumabz * it->second.sumabz ) / (it->second.nhets + 1e-6) / (it->second.nhets + 1e-6) + 1e-6));

      bcf_update_info_float(odw->hdr, wv, "ABE_AVG", &abe_avg, 1);
      bcf_update_info_float(odw->hdr, wv, "ABZ_AVG", &abz_avg, 1);
      bcf_update_info_float(odw->hdr, wv, "ABE_SD" , &abe_sd, 1);
      bcf_update_info_float(odw->hdr, wv, "ABZ_SD" , &abz_sd, 1);      

      // chisq
      double n0 = (it->second.nsamples - it->second.nhets - it->second.nalts) + 0.5;
      double n1 = it->second.nhets + 1.0;
      double n2 = it->second.nalts + 0.5;
      double p = (n1 + 2*n2)/(n0+n1+n2)/2;
      double e0 = (n0+n1+n2)*(1-p)*(1-p);
      double e1 = (n0+n1+n2)*(1-p)*(1-p);
      double e2 = (n0+n1+n2)*(1-p)*(1-p);      
      double chisq = (n0-e0)*(n0-e0)/e0 + (n1-e1)*(n1-e1)/e1 + (n2-e2)*(n2-e2)/e2;
      float z = (float)(sqrt(chisq) * (e1 < n1 ? 1 : -1));
      bcf_update_info_float(odw->hdr, wv, "EXHET" , &z, 1);
      */
      odw->write(wv);
    }
    
    odw->close();
  }

  void print_options()
  {
    std::clog << "merge_candidate_variants_sequential v" << version << "\n\n";
    std::clog << "options: [L] input VCF file list         " << input_vcf_file_list << " (" << input_vcf_files.size() << " files)\n";
    std::clog << "         [o] output VCF file             " << output_vcf_file << "\n";
    std::clog << "         [c] SNP variant score cutoff    " << snp_variant_score_cutoff << "\n";
    std::clog << "         [d] Indel variant score cutoff  " << indel_variant_score_cutoff << "\n";
    print_int_op("         [i] intervals                   ", intervals);
    std::clog << "\n";
  }

  void print_stats()
  {
    std::clog << "\n";
    std::clog << "stats: Total Number of Candidate SNPs                 " << no_candidate_snps << "\n";
    std::clog << "       Total Number of Candidate Indels               " << no_candidate_indels << "\n";
    std::clog << "\n";
  };
  
  ~Igor()
  {
  };
  
private:
};

} 

void merge_candidate_variants_sequential(int argc, char ** argv)
{
    Igor igor(argc, argv);
    igor.print_options();
    igor.initialize();
    igor.merge_candidate_variants_sequential();
    igor.print_stats();
}

