% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{small_trans}
\alias{small_trans}
\title{Transmission pairs iSNV}
\format{a data frame with 12 rows of 35 variables
\describe{
   \item{chr}{Chromosome - genomic segment}
   \item{pos}{position on genomic segment}
   \item{SPECID1}{Specimen id for first sample}
   \item{SPECID2}{Specimen id for second sample}
   \item{ENROLLID1}{Enroll id for first sample}
   \item{ENROLLID2}{Enroll id for second sample}
   \item{HOUSE_ID}{Id of the house of the individual}
   \item{onset1}{The date of symptom onset for the first sample}
   \item{onset2}{The date of symptom onset for the second sample}
   \item{transmission}{The infered date of transmission}
   \item{sequenced1}{Logical - was the first sample sequenced}
   \item{sequenced2}{Logical - was the second sample sequenced}
   \item{snv_qualified1}{Did the first sample qualify for isnv identification}
   \item{snv_qualified1}{Did the second sample qualify for isnv identification}
   \item{valid}{Is this a valid transmission pair}
   \item{double}{Does this need to be analyzed in both directions}
   \item{sequenced_pair}{Did the whole pair get sequenced}
   \item{titer_pair}{Are both samples above the titer cut off (10^3)}
   \item{snv_qualified_pair}{Did both samples qualify for isnv detection}
   \item{pair_id}{Id for the pairing}
   \item{quality_distance}{Does the sample meet the distance criteria for a transmission pair}
   \item{collect1}{Date of collection for the first sample}
   \item{collect2}{Date of collection for the second sample}
   \item{mutation}{mutation as chr_refPosVar}
   \item{ref}{nucleotide in the plasmid control}
   \item{var}{nucleotide in the sample}
   \item{pcr_result}{qpcr result for the samples-Strain}
   \item{freq1}{The frequency of the variant in the first sample}
   \item{freq2}{The frequency of the vairant in the second sample}
   \item{season}{The season of sampling}
   \item{gc_ul1}{The genomes/ul in the first sample}
   \item{gc_ul2}{The genomes/ul in the second sample}
   \item{freq}{What benchmarking frequency do we consider this to be at}
   \item{found}{Is the variant found in the second sample}
}}
\source{
data-raw/small_trans.R
}
\usage{
small_trans
}
\description{
A small silly example of transmission isnv comparision data.
}
\keyword{datasets}
