#' A small meta data file
#'
#' A small meta data file containing a few samples that
#' exihibit what we want to do in the code and a few intereseting
#' cases. Each row is a sample.
#'
#' @format A data frame with 20 rows and 17 variables:
#' \describe{
#'   \item{HOUSE_ID}{House identifier of sample}
#'   \item{ENROLLID}{personal identifier}
#'   \item{SPECID}{sample identifier}
#'   \item{onset}{Date of symptom onset}
#'   \item{collect}{date of sample collection}
#'   \item{vaccination_status}{vaccine status}
#'   \item{pcr_result}{strain H1N1 H3N2}
#'   \item{LAURING_ID}{Old sample id}
#'   \item{DPI}{Days post infection}
#'   \item{season}{Season of illness}
#'   \item{log_copy_num}{Log copy number of titer in qpcr well}
#'   \item{gc_ul}{genomes per ul in sample isolate}
#'   \item{HIGHSD}{Y or N for high standard deviation in qpcr}
#'   \item{sequenced}{T/F was the sample sequenced}
#'   \item{home_collected}{T/F was the sample taken at home}
#'   \item{snv_qualified}{T/F does the sample qualify for snv calling}
#'   \item{iSNV}{Number of iSNV in sample - not real data here}
#'
#'
#' }
#' @source data-raw/small_meta.R
"small_meta"


#' Sequencing specificity and sensitivity
#'
#' The sensitivity and specificity of the sequencing approach applied here -(frequency cut of 2%)
#' determined emphirically from seqeuncing populations of 20  mutants known a priori.
#'
#' @format A data frame with 6 rows and 7 variables:
#' \describe{
#'   \item{exp.freq}{The expected frequency from dilution series.}
#'   \item{freq}{The observed frequency - set to the expected here.}
#'   \item{TP}{The number of true positives found out of twenty.}
#'   \item{sensitivity}{Sensitivity to detect mutations at this frequency}
#'   \item{FP}{The number of false positives in the sample}
#'   \item{gc_ul}{The genome copies/ul }
#'   \item{PPV}{Positive predictive value}
#' }
#' @source data-raw/accuracy_stringent.R
"accuracy_stringent"

#' iSNV comparison between samples
#'
#' A silly example of what isnv comparisons between samples looks like. Note : this data
#' began as real and then was expanded to inlcude cases that are useful in testing code chunks.
#' It does not represent reality.
#'
#' @format a data frame with 32 rows and 13 variables
#' \describe{
#'    \item{chr}{Chromosome - genomic segment}
#'    \item{pos}{position on genomic segment}
#'    \item{SPECID1}{Specimen id for first sample}
#'    \item{SPECID2}{Specimen id for second sample}
#'    \item{mutation}{mutation as chr_refPosVar}
#'    \item{ref}{nucleotide in the plasmid control}
#'    \item{var}{nucleotide in the sample}
#'    \item{season}{The season of sampling}
#'    \item{pcr_result}{qpcr result for the samples-Strain}
#'    \item{freq1}{The frequency of the variant in the first sample}
#'    \item{freq2}{The frequency of the vairant in the second sample}
#'    \item{pair_id}{Id for the sample pairing}
#'    \item{found}{Is the variant found in the second sample}
#' }
#' @source data-raw/small_community.comp.R
"small_community.comp"

#' iSNV present in samples
#'
#' A small example of the data frames used to hold isnv calls.
#'
#' @format a data frame with 13 rows and 11 variables
#' #' \describe{
#'    \item{HOUSE_ID}{Id of the house of the individual}
#'    \item{ENROLLID}{Id of individual enrolled in cohort}
#'    \item{SPECID}{Specimen Id}
#'    \item{chr}{Chromosome - genomic segment}
#'    \item{pos}{position on genomic segment}
#'    \item{mutation}{mutation as chr_refPosVar}
#'    \item{ref}{nucleotide in the plasmid control}
#'    \item{freq.var}{Frequency of the variant}
#'    \item{var}{nucleotide in the sample}
#'    \item{season}{The season of sampling}
#'    \item{pcr_result}{qpcr result for the samples-Strain}
#' }
#' @source data-raw/small_isnv.R
"small_isnv"

#' Transmission pairs iSNV
#'
#' A small silly example of transmission isnv comparision data.
#'
#' @format a data frame with 12 rows of 35 variables
#' \describe{
#'    \item{chr}{Chromosome - genomic segment}
#'    \item{pos}{position on genomic segment}
#'    \item{SPECID1}{Specimen id for first sample}
#'    \item{SPECID2}{Specimen id for second sample}
#'    \item{ENROLLID1}{Enroll id for first sample}
#'    \item{ENROLLID2}{Enroll id for second sample}
#'    \item{HOUSE_ID}{Id of the house of the individual}
#'    \item{onset1}{The date of symptom onset for the first sample}
#'    \item{onset2}{The date of symptom onset for the second sample}
#'    \item{transmission}{The infered date of transmission}
#'    \item{sequenced1}{Logical - was the first sample sequenced}
#'    \item{sequenced2}{Logical - was the second sample sequenced}
#'    \item{snv_qualified1}{Did the first sample qualify for isnv identification}
#'    \item{snv_qualified1}{Did the second sample qualify for isnv identification}
#'    \item{valid}{Is this a valid transmission pair}
#'    \item{double}{Does this need to be analyzed in both directions}
#'    \item{sequenced_pair}{Did the whole pair get sequenced}
#'    \item{titer_pair}{Are both samples above the titer cut off (10^3)}
#'    \item{snv_qualified_pair}{Did both samples qualify for isnv detection}
#'    \item{pair_id}{Id for the pairing}
#'    \item{quality_distance}{Does the sample meet the distance criteria for a transmission pair}
#'    \item{collect1}{Date of collection for the first sample}
#'    \item{collect2}{Date of collection for the second sample}
#'    \item{mutation}{mutation as chr_refPosVar}
#'    \item{ref}{nucleotide in the plasmid control}
#'    \item{var}{nucleotide in the sample}
#'    \item{pcr_result}{qpcr result for the samples-Strain}
#'    \item{freq1}{The frequency of the variant in the first sample}
#'    \item{freq2}{The frequency of the vairant in the second sample}
#'    \item{season}{The season of sampling}
#'    \item{gc_ul1}{The genomes/ul in the first sample}
#'    \item{gc_ul2}{The genomes/ul in the second sample}
#'    \item{freq}{What benchmarking frequency do we consider this to be at}
#'    \item{found}{Is the variant found in the second sample}
#' }
#' @source data-raw/small_trans.R
"small_trans"

#' Variant data
#'
#' And example of the iSNV data that comes out of the variant_calling_pipeline.
#'
#' @format a data frame with 24 rows and 33 variables.
#' \describe{
#'    \item{Id}{The Id of the indexed sample for sequencing - can be multiple Id/ SPECID}
#'    \item{MapQ}{The average MapQ of reads with the variant}
#'    \item{Phred}{The average Phred of the variant base}
#'    \item{Read_pos}{The average position of the variant base on the read}
#'    \item{chr}{Chromosome - genomic segment}
#'    \item{cov.ctrl.bw}{The coverage of the plasmid control at the position in the reverse read}
#'    \item{cov.ctrl.fw}{The coverage of the plasmid control at the position in the forward read}
#'    \item{cov.tst.bw}{The coverage of the sample at the position in the reverse read}
#'    \item{cov.ctrl.fw}{The coverage of the sample at the position in the forward read}
#'    \item{freq.var}{Frequency of the variant}
#'    \item{mutation}{mutation as chr_refPosVar}
#'    \item{n.ctrl.bw}{The number of variant bases in the plasmid control in the reverse read}
#'    \item{n.ctrl.fw}{The number of variant bases in the plasmid control  in the forwad read}
#'    \item{n.tst.bw}{The number of variant bases in the sample in the reverse read}
#'    \item{n.ctrl.fw}{The number of variant bases in sample in the forward read}
#'    \item{p.val}{DeepSNV p.val}
#'    \item{pos}{position on genomic segment}
#'    \item{raw.p.val}{Uncorrected deepSNV p.val}
#'    \item{ref}{nucleotide in the plasmid control}
#'    \item{sigma2.freq.var}{Estimated variance of frequency measurement}
#'    \item{var}{nucleotide in the sample}
#'    \item{run}{The sequencing run the sample was on}
#'    \item{OR}{A list of the open reading frames overlapping this position}
#'    \item{coding_pos}{A list of the positions relative to the start of the in the ORs}
#'    \item{Ref_AA}{Reference amino acid - either sample or plasmid consensus depending on pipeline variables}
#'    \item{AA_pos}{Codon of the ORs - a list}
#'    \item{Var_AA}{Variant amino acid}
#'    \item{Class}{List of variant classifications in each relative OR}
#'    \item{LAURING_ID}{Sample ID used in processing - akin to SPECID}
#'    \item{dup}{Sequencing sample duplicate - if applicable}
#'    \item{coverage}{Total coverage at this position}
#'    \item{concat.pos}{Position in concatenated genome -
#'    segments in the order used in alignment fasta file}
#'    \item{gc_ul}{genome copy/ul of sample}
#'    }
#'    @source data-raw/variants.R
"variants"


