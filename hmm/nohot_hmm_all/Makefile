# Mon Apr  9 20:06:15 PDT 2012 Kevin Karplus

#	This Makefile is not intended as a standalone make, but
#	to provide build methods for targets that are requested
#	by prcesses that can be run in parallel, controlled by
#	an external script.

#	Many of these targets will need to have make run with macro definitions
#	overriding the defaults, to customize the target.


#   The overall structure of the process is
#
#    STAGE 1: search for reads that overlap with regions of a homologous sequence
#
#    STAGE 2: for each region (or each pair of adjacent regions), select
#    		a small number of seeds from the region (or that are
#    		in both).
#
#    STAGE 3: for each seed, build a consensus sequence using a few iterations of
#    		search and train HMM with model surgery.
#
#    STAGE 4: Select the HMM and consensus that gets the best scores
#	for each pair of regions.
#
#    STAGE 5: merge the selected sequences in order and
#    		retrain the combined consensus.
#

HOSTNAME:=$(shell hostname)

ifeq (${HOSTNAME},stardance)
   USE_SOE_FILES=1
endif
ifeq (${HOSTNAME},slamdance.soe.ucsc.edu)
   USE_SOE_FILES=1
endif
ifeq (${HOSTNAME},isla)
   USE_SOE_FILES=1
endif
ifeq (${HOSTNAME},cheep.cse.ucsc.edu)
   USE_SOE_FILES=1
endif

CORES := $(shell grep -P 'processor\t:\s\d+' /proc/cpuinfo | wc -l)

ifeq (${CORES},0)
	NUM_CPU = 1
else
ifeq (${CORES},1)
	NUM_CPU = 1
else
	NUM_CPU := $(shell python -c 'print(str(${CORES}-1))')
endif
endif

ifdef USE_SOE_FILES
# Important directories (on isla and cheep)
 PCL = /projects/compbio/lib
 PLUCK_SCRIPTS = /projects/compbio/experiments/pluck/scripts
 GATK = /projects/compbio/bin
 PICARD = /projects/compbio/bin/picard-tools-1.92
 TOP = /projects/compbio/experiments/pluck/Hp
 PACBIO_DIR ?= ${TOP}/PacBio
 READS_DIR ?= ${PACBIO_DIR}/nohotstart/data
 SOLNICK_DIR ?= ${TOP}/solnick-data
else
# Important directories (on campusrocks)
 PCL = /campusdata/BME235/lib/pcl
 PLUCK_SCRIPTS = /campusdata/BME235/bin/pluck-scripts
 READS_DIR ?= /campusdata/karplus/solnick_data
endif

PREFIX ?= c2_nohot
READS_ROOT ?= filtered_subreads

PACBIO_READS := ${READS_DIR}/${READS_ROOT}.fasta

NUM_READS := $(word 1, $(shell grep '>' ${PACBIO_READS} | wc))
$(warning ${PACBIO_READS} has ${NUM_READS} reads.)

# (mostly) verified sequence used for debugging---it should not
#	be part of any real assembly method.
# REAL_SEQ = ${PACBIO_DIR}/forward_gene/draft1.fna
REV_REAL_SEQ = ${PACBIO_DIR}/color-map-hand4/hand5.fasta
REAL_SEQ = probable_real.fasta


# Reference gene used for ordering and orienting reads for first iteration
# It should not be used after the first iteration, nor should it be used
# for anything other than selecting or orienting reads.
J99     ?= ${SOLNICK_DIR}/J99-cagY.fna
J99_SEQ ?= heliPylo_J99_dna

# V42 is the best draft of the SS1 sequence (but without a good cagY gene)
V42 ?= ${TOP}/chrom_v42_mapping1/v42.fa


# What should the combined output be called?
REGION_AB ?= region_1_2

# In stage3, we want the SEED_NUMth read in ${REGION_AB}.ids as a seed,
# but we have to orient it.
# SEED_NUM and SEED_NAME are also used in hmm_iter
SEED_NUM ?= 1
SEED_NAME ?= ${REGION_AB}_seed${SEED_NUM}



# regularizers for SAM's HMM programs
TRANSITION_REGULARIZER ?= ${PACBIO_DIR}/stiff_pacbio2.regularizer
DNA_REGULARIZER ?= ${PCL}/very_pure_dna.4comp

######################################
# Some useful macros for computation
######################################

# subtract one from the numeric argument
define MINUS_ONE
$(shell echo "print $1 -1" | perl)
endef

# $(call LINES, filename)  is the number of lines in the file
define LINES
 $(word 1,$(shell wc $1))
endef

# $(call RANGE, number) is the list of words 1 2 3 ... number
define RANGE
  $(shell python -c 'print(" ".join([str(x) for x in range(1,$1+1)]))')
endef

#lines_test: 
#	echo $(call LINES, Makefile)
#	echo $(call RANGE, $(call LINES, Makefile))



################################################
# Parameters that change on each iteration
################################################

# Instead of selection thresholds by bitscore, 
#	this version uses selection thresholds by approx number of reads.
# TRAINING and ALIGNMENT specify how many reads to include in each iteration
#	for training the HMM and for aligning to make the new consensus

TRAINING:=	30	40	45	50	100	120	150	200	200
ALIGNMENT:=	75	90	120	140	150	170	200	300	3000
#TRAINING:=	150	200	225	250	500	600	750	1000	1000
#ALIGNMENT:=	375	450	600	700	750	850	1000	1500	15000

# Parameters for consensus-from-sam
COVERAGES_MATCH:=	2	3	6	9	12	12	12	12	20
COVERAGES_INDEL:=	2	3	6	9	12	12	12	12	20
#COVERAGES_INDEL:=	0	0	0	0	0	0	0	0	0

# can use longest insert, or do no extension
USE_ALIGN:=	NO	LONG	NO	NO	NO	LONG	NO	NO	NO

MIN_PROB_MATCH = 0.5
MIN_PROB_INDEL = 0.5

#HMM_CUT = 0.5	# what fraction of sequences need to use state to keep it a match state?
		# Increase HMM_CUT if there are too many insertions in the final consensus
		# Decrease HMM_CUT if there are too many deletions in the final consensus

# How much of the end of the maxp consensus sequence should
#	be used for enriching the reads at the ends
END_LENGTH ?= 400

# How many reads should the ends be enriched by for TRAINING
#TRAIN_ENRICH_ENDS ?= 20

# How many reads should the ends be enriched by for ALIGNMENT
#ALIGN_ENRICH_ENDS ?= 40

default:	serial-make

parallel-make:

serial-make:	DATE
	# define the seeds
	${MAKE}  REGION=region1  REGION_BEGIN=1 REGION_END=40 BLASTN_E=1.e-02 STAGE1=ON stage1
	${MAKE}  REGION=region2  REGION_BEGIN=620 REGION_END=1190 BLASTN_E=1.e-02 STAGE1=ON stage1
	${MAKE}  REGION=region3  REGION_BEGIN=3830 REGION_END=5460 BLASTN_E=1.e-02 STAGE1=ON stage1
	${MAKE} REGION_A=region1 REGION_B=region2 REGION_AB=region1_2 REGION_AB_NUM_IDS=5 BLASTN_E=1.e-02 STAGE2=ON stage2
	${MAKE} REGION_A=region2 REGION_B=region3 REGION_AB=region2_3 REGION_AB_NUM_IDS=5 BLASTN_E=1.e-02 STAGE2=ON stage2
	${MAKE} REGION_AB=region1_2 stage3_4_sam
	${MAKE} REGION_AB=region2_3 stage3_4_sam
	${MAKE} 'REGIONS_AND_GUIDES=region1_2-iter3 region2 region2_3'-iter3 MERGE_NAME=all ITER=3 STAGE5=ON stage5
	$(foreach iter, 4 5 6 7 8 9, ${MAKE} SEED_NAME=all REGION_AB=none ITER=${iter} SAM_ITER=ON sam_iter;)

ifneq ($(wildcard ${REGION_AB}_seeds.ids*),)
 REGION_AB_RANGE := $(call RANGE, $(call LINES, ${REGION_AB}_seeds.ids))
 $(warning ${REGION_AB}_seeds.ids has $(call LINES, ${REGION_AB}_seeds.ids) lines.)
 $(warning REGION_AB_RANGE set to ${REGION_AB_RANGE})
endif

stage3_4_sam:	
	$(foreach n,${REGION_AB_RANGE}, ${MAKE} SEED_NUM=${n} STAGE3=ON stage3_seed; ) 
	$(foreach n,${REGION_AB_RANGE}, $(foreach iter,1 2 3, \
		${MAKE} SEED_NUM=${n} ITER=${iter} SAM_ITER=ON sam_iter; ) )
	${MAKE} 'MERGE_SEEDS=${REGION_AB_RANGE}' ITER=3 STAGE4=ON stage4


ifdef STAGE1
#########################################
#    STAGE 1: 
#	Search for reads that overlap with regions of a homologous sequence.
#	The "stage1" target does one search for one region.
########################################


# These three macros should be overwritten to define a particular
#	region name and boundaries when stage1 is made.
REGION ?= all
REGION_BEGIN ?= 1
REGION_END ?= 0


stage1:	${REGION}.fna	\
	Pac_reads_in_${REGION}.blastn.gz	\
	stage1.DONE

${REGION}.fna:	${J99}
	echo "${REGION} ${J99_SEQ} ${REGION_BEGIN} ${REGION_END}" \
	| ${PLUCK_SCRIPTS}/extract-fragments $< \
	> $@


.PRECIOUS:	${REGION}.fna	\
	Pac_reads_in_${REGION}.blastn.gz

endif

ifdef STAGE2
#########################################################################
#    STAGE 2: for each region (or each pair of adjacent regions), select
#    		a small number of seeds from the region (or that are
#    		in both).
#########################################################################

# Which regions are to be combined?
#	(set REGION_A and REGION_B to same value for single region)
REGION_A ?= region1
REGION_B ?= region2

# how many read ids do we want from the combined regions?
REGION_AB_NUM_IDS ?= 4


stage2:	${REGION_AB}_consistent.ids \
	${REGION_AB}.blastn_merged \
	${REGION_AB}_seeds.ids \
	${REGION_AB}_seeds.fna \
	${REGION_AB}_seeds_in_ref.blastn.gz \
	stage2.DONE

${REGION_AB}_consistent.ids: Pac_reads_in_${REGION_A}.blastn.gz Pac_reads_in_${REGION_B}.blastn.gz 
	${PLUCK_SCRIPTS}/direction-from-blast \
	    --only_ids --overlap=0\
	    --split --warning=2m  \
	    --region=${REGION_A} --region=${REGION_B} \
	    $^ \
	    < ${PACBIO_READS} \
	    > $@


# create file whose lines are "read_id region_a_bitscore region_b_bitscore"
# using only those that are in ${REGION_AB}_consistent.ids
${REGION_AB}.blastn_merged: \
		Pac_reads_in_${REGION_A}.blastn.gz \
		Pac_reads_in_${REGION_B}.blastn.gz \
		${REGION_AB}_consistent.ids
	gunzip -c Pac_reads_in_${REGION_A}.blastn.gz \
	| awk '{print $$12, $$1}' \
	| uniq -f 1 \
	| sort -k 2 \
	> ${REGION_A}.tmp
	gunzip -c Pac_reads_in_${REGION_B}.blastn.gz \
	| awk '{print $$12, $$1}' \
	| uniq -f 1 \
	| sort -k 2 \
	> ${REGION_B}.tmp
	join -j 2 -o '1.2 1.1 2.1' ${REGION_A}.tmp ${REGION_B}.tmp \
	| grep -F -f ${REGION_AB}_consistent.ids \
	> $@
	rm ${REGION_A}.tmp ${REGION_B}.tmp

# number of reads to look for in each region
NUM_TOP_READS := $(shell python -c 'print( max(1,int(${NUM_READS}*0.004)))')

# Move the name to the end, and put the two bit scores in front
# Sort on the name and field 1 or 2 (%) (one of the two bitscores, decreasing)
# Remove duplicate names
# Resort on bitscores (decreasing)
# Keep only the top NUM_TOP_READS
# Resort on name
${REGION_AB}_seeds.ids_%.tmp: ${REGION_AB}.blastn_merged 
	cat $< \
	| awk '{print $$2, $$3, $$1}' \
	| sort -k 3,3b  -k $*,$*rg \
	| uniq --skip-fields=2 \
	| sort -k $*,$*rg \
	| head -${NUM_TOP_READS} \
	| sort -k 3,3b \
	> $@
	
# Select those reads in both top sets
#	sort by product of the 2 bitscores,
#	print just the names
#	and select the top ${REGION_AB_NUM_IDS} 	
${REGION_AB}_seeds.ids:	 ${REGION_AB}_seeds.ids_1.tmp ${REGION_AB}_seeds.ids_2.tmp
	join -j 3 -o 1.1,1.2,1.3  $^ \
	| awk '{print $$1*$$2, $$3}' \
	| sort -rg -k 1 \
	| awk '{print $$2}' \
	| head -${REGION_AB_NUM_IDS} \
	> $@
	
.PRECIOUS: ${REGION_AB}.blastn_merged \
	${REGION_AB}_seeds.ids	\
	${REGION_AB}_seeds.fna \
	${REGION_AB}_seeds_in_ref.blastn.gz 

endif

ifdef STAGE3
#############################################################
#    STAGE 3: for each seed, 
#	build a consensus sequence using a few iterations of 
#	search and train .sam alignment with model surgery. 
#############################################################


# stage3_seed selects and orients one seed sequence.

stage3_seed: ${SEED_NAME}_read.ids \
	oriented_${SEED_NAME}.fasta.gz \
	oriented_${SEED_NAME}-length.rdb \
	oriented_${SEED_NAME}-longest.ids \
	${SEED_NAME}-iter0.fna \
	${SEED_NAME}-stage3_seed.DONE

${SEED_NAME}_read.ids: ${REGION_AB}_seeds.ids
	head -${SEED_NUM} < $< \
	| tail -1 \
	> $@

oriented_${SEED_NAME}.fasta.gz: ${SEED_NAME}_read.ids ${REGION_AB}_seeds_in_ref.blastn.gz
	${PLUCK_SCRIPTS}/direction-from-blast \
	--split --overlap=20 \
	--warning=s  \
	--ids=$^   \
	< ${REGION_AB}_seeds.fna \
	| gzip -9 \
	> $@

${SEED_NAME}-iter0.fna: oriented_${SEED_NAME}-longest.ids 	oriented_${SEED_NAME}.fasta.gz
	cat $< 		\
	| awk '{print "$(@:.fna=)", $$1, 1, 0}' \
	| ${PLUCK_SCRIPTS}/extract-fragments \
		oriented_${SEED_NAME}.fasta.gz \
	> $@

.PRECIOUS: ${SEED_NAME}_read.ids \
	oriented_${SEED_NAME}.fasta.gz \
	oriented_${SEED_NAME}-length.rdb \
	oriented_${SEED_NAME}-longest.ids \
	${SEED_NAME}-iter0.fna 

endif

ifdef HMM_ITER
#################################################
#    HMM_ITER: hmm_iter is an iteration using
#	an nhmmer to create a msa alignment to 
#	generate variant calls that are used in 
#	converting ${SEED_NAME}-iter${PREV_ITER}.fna
#	into ${SEED_NAME}-iter${ITER}.fna
#
#	It is used in stages 3, 4, and 5
#################################################

ITER ?= 1
PREV_ITER := $(call MINUS_ONE,${ITER})

# hmm_iter does one iteration of improving the consensus sequence
#	creating ${SEED_NAME}-iter${ITER}.fna

# scan with nhmmer -> Stockholm msa -> hmmbuild -> hmmpress -> iterate with nhmmscan + hmmalign? 

hmm_iter:	${SEED_NAME}-iter${PREV_ITER}.fna \
		${SEED_NAME}-iter${PREV_ITER}.sto \
		${SEED_NAME}-iter${ITER}.fna \
		selection_for_train_${SEED_NAME}-iter${PREV_ITER}.fasta.gz \
		selection_for_align_${SEED_NAME}-iter${PREV_ITER}.fasta.gz \
		${SEED_NAME}-sam_iter.DONE

%.sto: %.fna ${PACBIO_READS}
	mv $< $*.fasta
	nhmmer -A $*.sto -o /dev/null $*.fasta ${PACBIO_READS}
	mv $*.fasta $<

%-iter${ITER}.fna: %-iter${ITER}.sto %-iter${PREV_ITER}.fna ${PACBIO_READS}}
	${PLUCK_SCRIPTS}/training-and-consensus-using-nhmmer \
		--

%-iter${ITER}.merged_sorted.bam: %-iter${PREV_ITER}.sorted.bam %-iter${PREV_ITER}-begin.sorted.bam %-iter${PREV_ITER}-end.sorted.bam
	samtools view \
	    -b \
	    -o $*-iter${ITER}.merged.bam \
	    $*-iter${PREV_ITER}.sorted.bam $*-iter${PREV_ITER}-begin.sorted.bam $*-iter${PREV_ITER}-end.sorted.bam
	samtools sort \
	    $*-iter${ITER}.merged.bam \
	    $*-iter${ITER}.merged_sorted
	samtools index \
	    $@
	rm $*-iter${PREV_ITER}.sorted.bam $*-iter${PREV_ITER}-begin.sorted.bam $*-iter${PREV_ITER}-end.sorted.bam $*-iter${ITER}.merged.bam

ifeq ($(word ${ITER}, ${USE_ALIGN}),LONG)
%-iter${ITER}.fna: %-iter${ITER}.sorted.bam %-iter${PREV_ITER}.fna ${PACBIO_READS}
	${PLUCK_SCRIPTS}/training-and-consensus-from-sam \
	    --alignment=$< \
	    --reads=${PACBIO_READS} \
	    --reference=$*-iter${PREV_ITER}.fna \
	    --parallel=${NUM_CPU} \
	    --limit=10 \
	    --extend \
	    --num_train_reads=$(word ${ITER}, ${TRAINING}) \
	    --num_align_reads=$(word ${ITER}, ${ALIGNMENT}) \
	    --min_coverage=$(word ${ITER}, ${COVERAGES_MATCH}) \
	    --min_prob=${MIN_PROB_MATCH} \
	    --min_coverage_ins=$(word ${ITER}, ${COVERAGES_INDEL}) \
	    --min_prob_ins=${MIN_PROB_INDEL} \
	    --min_coverage_del=$(word ${ITER}, ${COVERAGES_INDEL}) \
	    --min_prob_del=${MIN_PROB_INDEL} \
	    --name=$(@:.fna=) \
	    --rle \
	    --use_evalue \
	> $@
else
%-iter${ITER}.fna: %-iter${ITER}.sorted.bam %-iter${PREV_ITER}.fna ${PACBIO_READS}
	${PLUCK_SCRIPTS}/training-and-consensus-from-sam \
	    --alignment=$< \
	    --reads=${PACBIO_READS} \
	    --reference=$*-iter${PREV_ITER}.fna \
	    --parallel=${NUM_CPU} \
	    --limit=10 \
	    --num_train_reads=$(word ${ITER}, ${TRAINING}) \
	    --num_align_reads=$(word ${ITER}, ${ALIGNMENT}) \
	    --min_coverage=$(word ${ITER}, ${COVERAGES_MATCH}) \
	    --min_prob=${MIN_PROB_MATCH} \
	    --min_coverage_ins=$(word ${ITER}, ${COVERAGES_INDEL}) \
	    --min_prob_ins=${MIN_PROB_INDEL} \
	    --min_coverage_del=$(word ${ITER}, ${COVERAGES_INDEL}) \
	    --min_prob_del=${MIN_PROB_INDEL} \
	    --name=$(@:.fna=) \
	    --rle \
	    --use_evalue \
	> $@
endif

selection_for_train_%-iter${PREV_ITER}.fasta.gz: %-iter${ITER}.fna
	gzip -9 $(@:.gz=)

selection_for_align_%-iter${PREV_ITER}.fasta.gz: %-iter${ITER}.fna
	gzip -9 $(@:.gz=)

.PRECIOUS:	\
		${SEED_NAME}-iter${ITER}.merged_sorted.bam \
		${SEED_NAME}-iter${ITER}.merged_sorted.bai \
		${SEED_NAME}-iter${ITER}.fna \
		selection_for_train_${SEED_NAME}-iter${ITER}.fasta.gz \
		selection_for_align_${SEED_NAME}-iter${ITER}.fasta.gz \

endif

ifdef STAGE4
###################################################################################
# when it works, update description
#    STAGE 4: 
#	for each pair of adjacent regions, select the HMM and
#	consensus sequence that scores reads the best.
###################################################################################

PREV_ITER := $(call MINUS_ONE,${ITER})

# which Pac reads to merge into consensus pool
SELECT_READS ?= $(foreach s,${REGION_AB_RANGE},selection_for_align_${REGION_AB}_seed${s}-iter${PREV_ITER}.fasta.gz)
$(warning SELECT_READS set to ${SELECT_READS})

# select stats to generate from alignments
SELECT_STATS ?= $(foreach s,${REGION_AB_RANGE},${REGION_AB}_seed${s}-iter${ITER}.stats)
$(warning SELECT_STATS set to ${SELECT_STATS})

stage4:	${REGION_AB}-iter${ITER}_merged_pool.fasta \
	${REGION_AB}-iter${ITER}.selected-fna \
	${REGION_AB}-iter${ITER}.fna \
	${REGION_AB}-iter${ITER}.DONE

# Align the merged pool for each consensus sequence using blasr
# Then get .sam error score with Python script
# blasr only likes .fasta file extension, do something about that
%-iter${ITER}.stats: ${REGION_AB}-iter${ITER}_merged_pool.fasta %-iter${ITER}.fna
	mv $*-iter${ITER}.fna $*-iter${ITER}.fasta
	blasr $< $*-iter${ITER}.fasta \
	    -nproc ${NUM_CPU} \
	    -unaligned $*-iter${ITER}_select_seed.unmapped.fna \
	    -sam \
	    -out $*-iter${ITER}_select_seed.sam
	mv $*-iter${ITER}.fasta $*-iter${ITER}.fna
	${PLUCK_SCRIPTS}/errorscore_sam_alignment \
	    -u $*-iter${ITER}_select_seed.unmapped.fna \
	    $*-iter${ITER}_select_seed.sam \
	    > $@

# Combine Pac reads, removing duplicates with fastx_collapser
${REGION_AB}-iter${ITER}_merged_pool.fasta: ${SELECT_READS}
	cat $^ \
	| gunzip -c \
	| sed -e '/>/n;s/ //g' \
	| fasta_formatter \
	| fastx_collapser \
	    -o $@

# Select the top scoring consensus sequence, looking for highest alignment error score (higher is better)
${REGION_AB}-iter${ITER}.selected-fna: ${SELECT_STATS}
	grep "" ${REGION_AB}*.stats \
	| sed -e "s/:/\t/g" \
	| sort -rn -k 2 \
	| head -1 \
	| awk '{print $$1}' \
	| sed -e 's/.stats/.fna/' \
	> $@

# Copy best consensus sequence
${REGION_AB}-iter${ITER}.fna: ${REGION_AB}-iter${ITER}.selected-fna
	cp  `cat $<` $@

.PRECIOUS: ${REGION_AB}-iter${ITER}_merged_pool.fasta \
	%-iter${ITER}.stats \
	${REGION_AB}-iter${ITER}.selected-fna \
	${REGION_AB}-iter${ITER}.fna

endif

ifdef STAGE5
################################################################
#    STAGE 5: merge the superconsensus sequences in order and
#    		retrain the combined consensus.
################################################################

# which regions to merge
REGIONS_AND_GUIDES?= region1_2-iter3 region2 region2_3-iter3

MERGE_NAME ?= foo

stage5: ${MERGE_NAME}-iter${ITER}.fna \
	${MERGE_NAME}-iter${ITER}.DONE

${MERGE_NAME}-iter${ITER}.fna: $(foreach region,${REGIONS_AND_GUIDES},${region}.fna)
	${PLUCK_SCRIPTS}/merge-sequences-with-blat --guides --name=$(@:.fna=) $^  >$@

endif


###################################
# TARGETS FOR MARKING UP LOG FILES
###################################


DATE:	
	date

%.DONE:	
	date
	echo $* finished

##################
# BLASTN searches
##################

BLASTN_E ?= 1.0e-04

BLAST_EXTENSIONS := nhr nin nsd nsi nsq


Pac_reads_in_%.blastn.gz: %.maxp ${PACBIO_READS}
	formatdb -i $< -p F -o T
	blastall -a ${NUM_CPU} -p blastn -m8 -d $< -F F -e ${BLASTN_E} -i ${PACBIO_READS} \
	| gzip -9 \
	> $@
	rm $(foreach ext,${BLAST_EXTENSIONS}, $<.${ext})

Pac_reads_in_%.blastn.gz: %.fna ${PACBIO_READS}
	formatdb -i $< -p F -o T
	blastall -a ${NUM_CPU} -p blastn -m8 -d $< -F F -e ${BLASTN_E} -i ${PACBIO_READS} \
	| gzip -9 \
	> $@
	rm $(foreach ext,${BLAST_EXTENSIONS}, $<.${ext})

%_in_ref.blastn.gz: %.fna ${J99} ${J99}.nsq
	blastall -a ${NUM_CPU} -p blastn -m8 -d ${J99} -F F -e ${BLASTN_E} -i $< \
	| gzip -9 \
	> $@

%_in_real.blastn.gz: %.fna ${REAL_SEQ} ${REAL_SEQ}.nsq
	blastall -a ${NUM_CPU} -p blastn -m8 -d ${REAL_SEQ} -F F -e ${BLASTN_E} -i $< \
	| gzip -9 \
	> $@

%_in_v42.blastn.gz: %.fna ${V42} ${V42}.nsq
	blastall -a ${NUM_CPU} -p blastn -m8 -d ${V42} -F F -e ${BLASTN_E} -i $< \
	| gzip -9 \
	> $@


%.nsq: %
	formatdb -i $* -p F -o T

.PRECIOUS:	\
	Pac_reads_in_%.blastn.gz	\
	%.nsq \
	%.mod %.a2m.gz %.maxp	

####################################
# Consensus Sequence using .saves
####################################

LOGO_WIDTH ?= 200

.PRECIOUS: %.maxp %.saves %-logo.eps %-logo.pdf

%.logo:	%.saves %.maxp %-logo.eps %-logo.pdf
	echo $@ made

%.saves: %.mod 
	makelogo $*.tmp-logo -i $< \
		-logo_rel_entropy 1 \
		-logo_bars_per_line ${LOGO_WIDTH}  -logo_title "$*"  \
		-logo_savings_output $@
	rm $*.tmp-logo.eps

%.maxp: %.saves 
	${PLUCK_SCRIPTS}/consensus-from-saves --name $* < $< > $@

%-begin.fna: %.maxp
	echo "$*-begin $* 1 ${END_LENGTH}" \
	| ${PLUCK_SCRIPTS}/extract-fragments --err=mn $< \
	> $@

%-end.fna: %.maxp
	echo "$*-end $* -${END_LENGTH} -1" \
	| ${PLUCK_SCRIPTS}/extract-fragments --err=mn $< \
	> $@


%-begin.fna: %.fna
	echo "$*-begin $* 1 ${END_LENGTH}" \
	| ${PLUCK_SCRIPTS}/extract-fragments --err=mn $< \
	> $@

%-end.fna: %.fna
	echo "$*-end $* -${END_LENGTH} -1" \
	| ${PLUCK_SCRIPTS}/extract-fragments --err=mn $< \
	> $@


%-logo.eps: %.mod %.maxp
	makelogo $*-logo -i $< \
		-logo_rel_entropy 1 \
		-logo_bars_per_line ${LOGO_WIDTH}  -logo_title "$*"  \
		-logo_under_file $*.maxp

%-logo.pdf:	%-logo.eps
	ps2pdf -dEPSCrop $^ $@

.PRECIOUS: %.maxp \
	%-begin.fna \
	%-end.fna \

#########################################################################
# Split alignment into 3 parts and redo beginning and ending with Muscle
# or kalign
#########################################################################

END_SPLIT ?=50

%-begin.a2m.gz: %.a2m.gz
	gunzip -c $< \
	| ${PLUCK_SCRIPTS}/extract-columns-from-a2m --low=0 --high=${END_SPLIT} --terminal_inserts \
	| gzip -9 > $@

%-mid.a2m.gz: %.a2m.gz
	gunzip -c $< \
	| ${PLUCK_SCRIPTS}/extract-columns-from-a2m --low=${END_SPLIT} --high=-${END_SPLIT} --terminal_inserts \
	| gzip -9 > $@

%-end.a2m.gz: %.a2m.gz
	gunzip -c $< \
	| ${PLUCK_SCRIPTS}/extract-columns-from-a2m --low=-${END_SPLIT} --high=0 --terminal_inserts \
	| gzip -9 > $@
 
%.muscle.a2m.gz: %.a2m.gz
	gunzip -c $< \
	| tr -d '.\-' \
	> $*.tmp
	muscle -in $*.tmp -out $(@:.gz=) -maxiters 2
	gzip -9 $(@:.gz=) 
	rm $*.tmp

%.kalign.a2m.gz: %.a2m.gz
	gunzip -c $< \
	| kalign -format fasta -c input \
	> $(@:.gz=)
	gzip -9 $(@:.gz=) 


.PRECIOUS: %-begin.a2m.gz \
	%-mid.a2m.gz \
	%-end.a2m.gz \
	%.muscle.a2m.gz \
	%.kalign.a2m.gz

###################################
# generic construction techniques
###################################


# cleanup
rm-empty: 
	find . -empty -exec rm -f '{}' \; -print

%-length.rdb: %.fasta.gz
	gunzip -c < $<  \
	| ${PLUCK_SCRIPTS}/make-contig-lengths \
	> $@

%.fna: %.ids ${PACBIO_READS}
	${PLUCK_SCRIPTS}/extract-fragments --whole_seqs ${PACBIO_READS} < $< > $@

%-longest.ids: %-length.rdb
	sorttbl FullLength < $< \
	| /projects/compbio/bin/column ID \
	| tail -1 \
	> $@

# This is a debugging aid, comparing to a (mostly) verified correct sequence.
# It should not be part of any actual assembly effort
%.differences_from_real: %.fna ${REAL_SEQ}
	${PLUCK_SCRIPTS}/find-dna-differences \
		${REAL_SEQ} --short=$@.short \
	< $< > $@

%.differences_from_real: %.maxp ${REAL_SEQ}
	${PLUCK_SCRIPTS}/find-dna-differences \
		${REAL_SEQ} --short=$@.short \
	< $< > $@



#################################
# BELOW HERE IS CURRENTLY UNUSED
#################################


${PREFIX}-cagY.fna: iter$(lastword ${ITERS})-all.a2m.gz
	${PLUCK_SCRIPTS}/consensus-from-a2m \
		--name=$(@:.fna=) \
		--min_coverage=10 \
		--inserts=no \
		--min_prob=0.6 \
		--caps \
		--a2m $^ \
	> $@



###################################
# generic construction techniques
###################################

probable_real.fasta: ${REV_REAL_SEQ}
	${PLUCK_SCRIPTS}/reverse-comp < $< > $@

reversed_%.fna:%.fna
	${PLUCK_SCRIPTS}/reverse-comp --modify_name <$< > $@

%.dotted-a2m: %.a2m
	prettyalign $<  -f > $@
	
%.run_lengths: %.a2m.gz
	gunzip -c $< \
	| ${PLUCK_SCRIPTS}/indel-distribution \
	> $@

%.transitions: %.a2m.gz
	${PLUCK_SCRIPTS}/count-transitions $^ > $@

%.homopolymer-rdb: %.a2m.gz %.fna
	${PLUCK_SCRIPTS}/homopolymer-stats -g $*.fna  $< > $@

%.homopolymer-joint: %.a2m.gz %.fna
	${PLUCK_SCRIPTS}/homopolymer-stats -g $*.fna --joint  $< > $@


%.AGCT : %.fasta.gz
	gunzip -c $< \
	| ${PLUCK_SCRIPTS}/count-AGCT \
	> $@


%.inserts_before: %.a2m.gz
	gunzip -c $< \
	| ${PLUCK_SCRIPTS}/inserts-before \
	> $@
	
%.long_inserts: %.a2m.gz
	gunzip -c $< \
	| ${PLUCK_SCRIPTS}/make-contig-lengths \
	| row FullLength "-" AlignLength '>=' 0.3 '*' AlignLength \
	> $@.ids
	${PLUCK_SCRIPTS}/extract-fragments --whole_seqs  $<  <$@.ids > $@

.PRECIOUS: \
	%-length.rdb \
	%.fna \
	reversed_%.fna \
	%.dotted-a2m \
	%.run_lengths \
	%.transitions \
	%.AGCT \
	%.inserts_before \
	%.long_inserts \
	%.differences_from_real


