\name{immuno}
\alias{immuno}
\alias{groups}
\docType{data}
\title{Gut Microbiome Data for Cancer Immunotherapy Response Prediction}
\description{
 Relative abundance data of gut microbial taxa from 219 cancer patients
 treated with immune checkpoint inhibitor (ICI) therapy. The dataset
 contains both species-level (n = 219 x 190) and genus-level
 features, along with a pre-built group partition \code{groups} that maps
 each genus to its member species. Patients are labelled as responders
 (R) or non-responders (NR) encoded in the row names.
}
\format{
 Loading \code{data(immuno)} places two objects in the workspace:

 \describe{
 \item{\code{immuno}}{A numeric data frame with 219 rows (patients) and
 190 columns (microbial features: species and unclassified taxa).
 Row names follow the pattern
 \code{p<id>_<cohort>_<response>}, where
 \code{<cohort>} is one of \code{Gop}, \code{Mat}, \code{Fra},
 \code{Rou}, or \code{Val}, and \code{<response>} is \code{R}
 (responder) or \code{NR} (non-responder).
 Column values are relative abundances (proportions; non-negative,
 summing to approximately 1 per row).}

 \item{\code{groups}}{A named list with 18 elements. Each element is a
 character vector of species column names belonging to that genus.
 The element names are genus labels (e.g., \code{"Akkermansia"},
 \code{"Bacteroides"}). Together, \code{groups} forms a complete,
 non-overlapping partition of all 190 columns of \code{immuno},
 suitable for direct use as the \code{groups} argument of
 \code{\link{treeowen}}.}
 }
}
\details{
 \subsection{Study design}{
 The dataset pools microbiome samples from five independent ICI cohorts:
 \describe{
 \item{\code{Gop}}{Gopalakrishnan et al. (2018) — melanoma; n = 25.}
 \item{\code{Mat}}{Matson et al. (2018) — melanoma; n = 39.}
 \item{\code{Fra}}{Frankel et al. (2017) — melanoma; n = 14.}
 \item{\code{Rou}}{Routy et al. (2018) — NSCLC / RCC; n = 127.}
 \item{\code{Val}}{Valpione et al. — melanoma; n = 14.}
 }
 Overall: 65 responders (R) and 154 non-responders (NR).
 }

 \subsection{Feature encoding}{
 Species names use dot notation: genus--species separator is \code{.}
 (single dot), and \code{..} encodes \emph{sp.} (unclassified at
 species level). Unclassified taxa at higher ranks follow the
 pattern \code{unknown.<taxon>} or \code{unknown_<phylum>}.
 }

 \subsection{Group partition}{
 The \code{groups} list assigns each of the 190 columns to one of 18
 genera. This partition is used with \code{\link{treeowen}} to compute
 genus-level Owen values, attributing model predictions jointly to all
 species within a genus before distributing within-genus attribution to
 individual species.
 }
}
\source{
 Pooled from publicly available microbiome sequencing data associated with
 the five ICI studies listed in the Details section.
}

\references{
Koh H. Efficient Game-Theoretic Explanations for Tree-Based Ensembles
via Owen Values. (In Review)
}
\seealso{
 \code{\link{treeowen}}, \code{\link{treeowen_importance}},
 \code{\link{treeowen_beeswarm}}, \code{\link{treeowen_hierarchical_beeswarm}}
}
\examples{
data(immuno)

# Dimensions
dim(immuno) # 219 x 190
length(groups) # 18 genera
names(groups)[1:5] # first five genus names

# Response labels from row names
resp_labels <- sub(".*_(R|NR)$", "\\1", rownames(immuno))
table(resp_labels) # NR: 154, R: 65

# First few features
head(colnames(immuno), 8)

# Relative abundance: each row sums to ~1
range(rowSums(immuno))

# Group sizes (species per genus)
sort(lengths(groups), decreasing = TRUE)
}
