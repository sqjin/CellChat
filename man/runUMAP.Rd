% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{runUMAP}
\alias{runUMAP}
\title{Run UMAP}
\usage{
runUMAP(
  data.use,
  n_neighbors = 30L,
  n_components = 2L,
  metric = "correlation",
  n_epochs = NULL,
  learning_rate = 1,
  min_dist = 0.3,
  spread = 1,
  set_op_mix_ratio = 1,
  local_connectivity = 1L,
  repulsion_strength = 1,
  negative_sample_rate = 5,
  a = NULL,
  b = NULL,
  seed.use = 42L,
  metric_kwds = NULL,
  angular_rp_forest = FALSE,
  verbose = FALSE
)
}
\arguments{
\item{data.use}{input data matrix}

\item{n_neighbors}{This determines the number of neighboring points used in
local approximations of manifold structure. Larger values will result in more
global structure being preserved at the loss of detailed local structure. In general this parameter should often be in the range 5 to 50.}

\item{n_components}{The dimension of the space to embed into.}

\item{metric}{This determines the choice of metric used to measure distance in the input space.}

\item{n_epochs}{the number of training epochs to be used in optimizing the low dimensional embedding. Larger values result in more accurate embeddings. If NULL is specified, a value will be selected based on the size of the input dataset (200 for large datasets, 500 for small).}

\item{learning_rate}{The initial learning rate for the embedding optimization.}

\item{min_dist}{This controls how tightly the embedding is allowed compress points together.
Larger values ensure embedded points are moreevenly distributed, while smaller values allow the
algorithm to optimise more accurately with regard to local structure. Sensible values are in the range 0.001 to 0.5.}

\item{spread}{he effective scale of embedded points. In combination with min.dist this determines how clustered/clumped the embedded points are.}

\item{set_op_mix_ratio}{Interpolate between (fuzzy) union and intersection as the set operation used to combine local fuzzy simplicial sets to obtain a global fuzzy simplicial sets.}

\item{local_connectivity}{The local connectivity required - i.e. the number of nearest neighbors
that should be assumed to be connected at a local level. The higher this value the more connected
the manifold becomes locally. In practice this should be not more than the local intrinsic dimension of the manifold.}

\item{repulsion_strength}{Weighting applied to negative samples in low dimensional embedding
optimization. Values higher than one will result in greater weight being given to negative samples.}

\item{negative_sample_rate}{The number of negative samples to select per positive sample in the
optimization process. Increasing this value will result in greater repulsive force being applied, greater optimization cost, but slightly more accuracy.}

\item{a}{More specific parameters controlling the embedding. If NULL, these values are set automatically as determined by min. dist and spread.}

\item{b}{More specific parameters controlling the embedding. If NULL, these values are set automatically as determined by min. dist and spread.}

\item{seed.use}{Set a random seed. By default, sets the seed to 42.}

\item{metric_kwds, angular_rp_forest, verbose}{other parameters used in UMAP}
}
\description{
Run UMAP
}
