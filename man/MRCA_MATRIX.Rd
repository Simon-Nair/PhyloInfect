\name{MRCA_MATRIX}
\alias{MRCA_MATRIX}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MRCAMATRIX
}
\description{
see README.Rmd for how to use package.
}
\usage{
MRCA_MATRIX(table, sample)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{table}{
%%     ~~Describe \code{table} here~~
}
  \item{sample}{
%%     ~~Describe \code{sample} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or standard data sets, see data().

## The function is currently defined as
function (table, sample)
{
    n <- length(sample)
    mrca_matrix <- matrix(NA, nrow = n, ncol = n, dimnames = list(sample,
        sample))
    mrca <- NA
    for (i in seq_along(sample)) {
        for (j in seq_along(sample)) {
            mrca <- MRCA(table, sample[i], sample[j])
            mrca_matrix[i, j] <- mrca
        }
    }
    lower_tri_index <- lower.tri(mrca_matrix)
    mrca_matrix[lower_tri_index] <- 0
    for (i in 1:nrow(mrca_matrix)) {
        for (j in 1:nrow(mrca_matrix)) {
            if (mrca_matrix[i, j] == 0) {
                mrca_matrix[i, j] <- 0
            }
            else {
                mrca_time <- table$time[which(table$affected ==
                  mrca_matrix[i, j])]
                mrca_matrix[i, j] <- mrca_time
            }
        }
    }
    return(mrca_matrix)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
