\name{MRCA}
\alias{MRCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MRCA
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MRCA(table, id1, id2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{table}{
%%     ~~Describe \code{table} here~~
}
  \item{id1}{
%%     ~~Describe \code{id1} here~~
}
  \item{id2}{
%%     ~~Describe \code{id2} here~~
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
function (table, id1, id2)
{
    get_ancestors <- function(id) {
        ancestors <- c(id)
        index <- which(table$affected == id)
        while (length(index) > 0) {
            infector <- table$infector[[index]]
            ancestors <- c(ancestors, infector)
            index <- which(table$affected == infector)
        }
        return(ancestors)
    }
    ancestors1 <- get_ancestors(id1)
    ancestors2 <- get_ancestors(id2)
    intersection <- intersect(ancestors1, ancestors2)
    ancestors1 <- setdiff(ancestors1, intersection)
    ancestors2 <- setdiff(ancestors2, intersection)
    A1 <- ancestors1[length(ancestors1)]
    A2 <- ancestors2[length(ancestors2)]
    ind1 <- 0
    ind2 <- 0
    mrca <- 0
    if (length(A1) == 0 && length(A2) == 0) {
        mrca <- 0
    }
    else if (length(A1) > 0 && length(A2) == 0) {
        mrca <- A1
    }
    else if (length(A1) == 0 && length(A2) > 0) {
        mrca <- A2
    }
    else if (length(A1) > 0 && length(A2) > 0) {
        ind1 <- which(table$affected == ancestors1[length(ancestors1)])
        ind2 <- which(table$affected == ancestors2[length(ancestors2)])
        if (!(length(ind1) == 0 && length(ind2) == 0)) {
            if (ind1 < ind2) {
                mrca <- A1
            }
            else {
                mrca <- A2
            }
        }
    }
    if (length(mrca) == 0) {
        mrca <- 0
    }
    return(mrca)
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
