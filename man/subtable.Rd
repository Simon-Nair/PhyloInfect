\name{subtable}
\alias{subtable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Subtable
}
\description{
see README.Rmd for how to use package.
}
\usage{
subtable(table, mrca_matrix, sample)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{table}{
%%     ~~Describe \code{table} here~~
}
  \item{mrca_matrix}{
%%     ~~Describe \code{mrca_matrix} here~~
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
function (table, mrca_matrix, sample)
{
    p_0 <- as.numeric(table$affected[1])
    if (sample[1] != p_0) {
        sampled_p_0 <- sample[1]
        sample[1] <- p_0
    }
    else {
        sampled_p_0 <- FALSE
    }
    drawn <- c()
    get_drawn_ancestors <- function(id) {
        ancestors <- c(id)
        index <- which(as.numeric(table$affected) == as.numeric(id))
        while (length(index) > 0 && !(table$infector[[index]] \%in\%
            drawn || table$infector[[index]] == p_0)) {
            infector <- table$infector[[index]]
            ancestors <- c(ancestors, infector)
            index <- which(table$affected == infector)
        }
        drawn <<- c(drawn, ancestors[length(ancestors)])
        return(ancestors)
    }
    for (i in 1:length(sample)) {
        get_drawn_ancestors(sample[i])
    }
    if (drawn[1] != p_0) {
        drawn[1] <- p_0
    }
    subtable <- table
    subtable <- subtable[subtable$affected \%in\% drawn | subtable$affected ==
        0, ]
    if (sampled_p_0) {
        sample[1] <- sampled_p_0
    }
    sample_w_0 <- c(0, sample)
    drawn_w_0 <- c(0, drawn)
    id_replacement <- setNames(sample_w_0, drawn_w_0)
    subtable$affected <- id_replacement[as.character(subtable$affected)]
    subtable$infector <- id_replacement[as.character(subtable$infector)]
    for (i in 2:nrow(subtable)) {
        row <- subtable$infector[[i]]
        col <- subtable$affected[[i]]
        row_index <- which(sample == row)
        col_index <- which(sample == col)
        subtable$time[i] <- mrca_matrix[row_index, col_index]
    }
    for (i in 3:nrow(subtable)) {
        time <- subtable$time[i]
        inf <- subtable$infector[i]
        while (inf != sample[1] && time < subtable$time[which(subtable$affected ==
            inf)]) {
            inf <- subtable$infector[which(subtable$affected ==
                inf)]
        }
        subtable$infector[i] <- inf
    }
    subtable <- subtable[order(subtable$time), ]
    for (i in 1:nrow(subtable)) {
        subtable$S[i] <- nrow(subtable) - i
        subtable$I[i] <- i
        subtable$R[i] <- 0
    }
    return(subtable)
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
