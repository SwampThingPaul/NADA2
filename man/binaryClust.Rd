% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/binaryClust.R
\name{binaryClust}
\alias{binaryClust}
\title{Cluster Matrix of Binary Censored Data}
\usage{
binaryClust(
  data,
  method = "ward.D2",
  group = NULL,
  ncluster = NULL,
  plotncluster = TRUE,
  clustIndex = "all"
)
}
\arguments{
\item{data}{A data frame containing only the 0/1 columns, one column per chemical parameter.}

\item{method}{Method of forming clusters.  The default is \code{"ward.D2"}, which is appropriate for a variety of types of data.  Another appropriate option is \verb{“average”} – average distances between cluster centers.  See the vegan package for other possible clustering methods.}

\item{group}{Optional grouping variable. If used, sites being clustered will be represented by their group name, rather than by the row number.}

\item{ncluster}{Optional number of clusters to be differentiated on the graph. Clusters are fenced off with rectangles.}

\item{plotncluster}{default is \code{TRUE} logical flags to add identification of clusters on dendrogram}

\item{clustIndex}{Optional, if not specified, potential number of clusters will be determined based on the mean best number of clusters across all indicies. For a specific index, see details}
}
\value{
Prints a cluster dendrogram based on clustering method and outputs a list of clusters and hierarchical cluster analysis results
}
\description{
Performs clustering of a matrix of 0s and 1s, ie. the censoring indicator columns for multiple variables. If there are multiple detection limits within a column first convert the 0/1 designated to be Below vs Above the highest detection limit in the column.  Detection limits may differ among columns.
}
\details{
If a specific index is desired to determine the best number of clusters see \code{NbClust::NbClust} for index values.
}
\examples{
data(PbHeron)

# without group specified
binaryClust(PbHeron[,4:15])

# With Group argument
binaryClust(PbHeron[,4:15],group=PbHeron$DosageGroup)
}
\references{
Helsel, D.R., 2011. Statistics for Censored Environmental Data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
}
