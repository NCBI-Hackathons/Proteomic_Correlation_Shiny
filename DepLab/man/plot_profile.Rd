% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plot_profile}
\alias{plot_profile}
\title{Plot the intensity profiles across several fractions}
\usage{
plot_profile(long.df, x = "fraction", y = "value", what = c("gene_symbol",
  "expt_id"), color.by = "gene_symbol", split.by = NULL,
  split.by.col = ".", color.palette = NULL, point.shape = 16L,
  point.size = 1L, line.type = 1L, line.size = 1L, y.lab = "value",
  title = NULL, line.smooth = FALSE)
}
\arguments{
\item{long.df}{data.frame fit for ggplot2}

\item{x}{expression indicating the column name to be used for the x-axis
values, e.g. "fraction"}

\item{y}{expression indicating the column name to be used for the x-axis
values, one of:
c("raw.intensity", "MS.MS.count", "peptides.count",
"unique.peptides.only", "razor.and.unique.peptides")}

\item{what}{indicate how the points that will be connected by lines should be
grouped together, e.g. c("gene_symbol", "expt_id") will draw a
separate line for each unique combination of gene_symbol and expt_id}

\item{color.by}{expression indicating the factor that will be used to assign
different colors, e.g. "gene.symbol"}

\item{split.by}{expression indicating the factor that will be used to split
the single plot into several ones, e.g. "date"}

\item{split.by.col}{factor that will be used to split the plot column-wise.
The default '.' is equivalent to 'NA'}

\item{color.palette}{name of RColorBrewer palette used for plotting; 
recommended palettes: "Set1", "Set2", "Set3", "Pastel2", "Paired", "Dark2",
 "Accent"; 
see \url{http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/figure/
unnamed-chunk-14-1.png} for more palettes}

\item{point.shape}{either a \emph{number} defining the shape for the points for all
groups or a \emph{string} which will be interpreted as a factor that the
 point shape should depend on; see
 \url{http://www.cookbook-r.com/Graphs/Shapes_and_line_types/} for shape types
 and the corresponding integers.}

\item{point.size}{number defining the size for the points for all
groups}

\item{line.type}{either a \emph{number} defining the line type for all
groups or a \emph{string} which will be interpreted as a factor that the 
line type should depend on; see \url{http://www.cookbook-r.com/Graphs/Shapes_and_line_types/}
 for line types and the corresponding integers}

\item{line.size}{number defining the width for the line for all
groups}

\item{y.lab}{title for the y axis; default: "value"}

\item{title}{title for the plot}

\item{line.smooth}{Boolean indicating whether a smoothened line should be displayed.}
}
\value{
A line plot visualizing the value of interest for the gene of
interest across all fractions.
}
\description{
Basic plot for visualizing the intensity of individual proteins across 
multiple fractions.
}
\details{
Expects a data.frame in the long format typical of ggplot2.
}
\seealso{
\code{\link{MQ_to_ggplot2}}
}
