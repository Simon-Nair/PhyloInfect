\name{subtree}
\alias{subtree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Subtree
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
subtree(table, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{table}{
%%     ~~Describe \code{table} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
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
function (table, n)
{
    I_list <- end_time_infected(table)
    if (length(I_list) < n) {
        stop("Error: not enough infected individuals remain alive at sample time, \n         try a different simulation or a smaller sample size.")
    }
    filtered_table <- table[2:nrow(table), ]
    filtered_table <- filtered_table[filtered_table$affected !=
        0, ]
    small_table <- rbind(table[1, ], filtered_table)
    sublist <- sublist(small_table, I_list, n)
    table$affected[1] <- 0
    mat <- MRCA_MATRIX(small_table, sublist)
    subtable <- subtable(small_table, mat, sublist)
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
