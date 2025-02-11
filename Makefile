OPTFLAG = -g -O3 -fopenmp
INCLUDES = -I./lib -I. -I../htslib -I./lib/Rmath -I./lib/pcre2
CFLAGS = -pipe -std=c++0x $(OPTFLAG) $(INCLUDES) -D__STDC_LIMIT_MACROS -no-pie
CXX = g++

SOURCESONLY =

SOURCES = align\
		allele\
		annotate_1000g\
		annotate_dbsnp_rsid\
		annotate_regions\
		annotate_variants\
		augmented_bam_record\
		bam_ordered_reader\
		bcf_genotyping_buffered_reader\
		bcf_ordered_reader\
		bcf_ordered_writer\
		bcf_synced_reader\
		bed\
		candidate_region_extractor\
		cat\
		chmm\
		compute_concordance\
		compute_features\
		compute_features2\
		config\
		consolidate_variants\
		construct_probes\
		decompose\
		decompose2\
		decompose_blocksub\
		discover\
		discover2\
		estimate\
		estimator\
		filter\
		flank_detector\
		gencode\
		genome_interval\
		genotype\
		genotype2\
		genotyping_record\
		hts_utils\
		hfilter\
		index\
		interval_tree\
		interval\
		joint_genotype_sequential\
		joint_genotyping_buffered_reader\
		joint_genotyping_record\
		joint_genotype_block\
		joint_genotype_block_reader\
		joint_genotype_block_record\
		lfhmm\
		lhmm\
		lhmm1\
		log_tool\
		merge\
		merge_candidate_variants\
		merge_candidate_variants2\
		merge_candidate_variants_sequential\
		milk_filter\
		motif_tree\
		motif_map\
		multi_partition\
		needle\
		normalize\
		nuclear_pedigree\
		ordered_bcf_overlap_matcher\
		ordered_region_overlap_matcher\
		partition\
		paste\
		paste_genotypes\
		paste_and_compute_features_sequential\
		pedigree\
		peek\
		pileup\
		pregex\
		profile_afs\
		profile_chm1\
		profile_chrom\
		profile_fic_hwe\
		profile_hwe\
		profile_indels\
		profile_len\
		profile_mendelian\
		profile_na12878\
		profile_snps\
		program\
		reference_sequence\
		remove_overlap\
		repeat_tract_detector\
		rfhmm\
		rminfo\
		rpartition\
		seq\
		sex_ploidy_map\
		sort\
		subset\
		sv_tree\
		test\
		union_variants\
		uniq\
		utils\
		validate\
		variant\
		variant_manip\
		variant_filter\
		view\
		vntr\
		vntr_annotator\
		vntrize\
		tbx_ordered_reader\
		ahmm\
		xcmp\
		Error\
		wphmm\
		wphmm_ungapped\
		shift_align\
		vntr_candidate\
		identify_vntr\

SOURCESONLY = main.cpp

TARGET = vt
TOOLSRC = $(SOURCES:=.cpp) $(SOURCESONLY)
TOOLOBJ = $(TOOLSRC:.cpp=.o)
LIBHTS = ../htslib/libhts.a
LIBRMATH = lib/Rmath/libRmath.a
LIBPCRE2 = lib/pcre2/libpcre2.a

all : $(TARGET)

${LIBHTS} :
	cd ../htslib; $(MAKE) libhts.a || exit 1;

${LIBRMATH} :
	cd lib/Rmath; $(MAKE) libRmath.a || exit 1;

${LIBPCRE2} :
	cd lib/pcre2; $(MAKE) libpcre2.a || exit 1;

$(TARGET) : ${LIBHTS} ${LIBRMATH} ${LIBPCRE2} $(TOOLOBJ)
	$(CXX) $(CFLAGS) -o $@ $(TOOLOBJ) $(LIBHTS) $(LIBRMATH) ${LIBPCRE2} -fopenmp -llzma -lbz2 -lz -lm -lcurl -lcrypto -lpthread

$(TOOLOBJ): $(HEADERSONLY)

.cpp.o :
	$(CXX) $(CFLAGS) -o $@ -c $*.cpp

.PHONY: clean cleanvt test

clean :
	cd ../htslib; $(MAKE) clean
	cd lib/Rmath; $(MAKE) clean
	cd lib/pcre2; $(MAKE) clean
	-rm -rf $(TARGET) $(TOOLOBJ)

cleanvt :
	-rm -rf $(TARGET) $(TOOLOBJ)

test: vt
	for x in ./test/*/_run.sh; do \
		$${x}; \
	done
