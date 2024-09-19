\name{phyloInfect}
\alias{phyloInfect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
PhyloInfect
}
\description{
see README.Rmd for how to use package.
}
\usage{
phyloInfect(table, title = "", end_time = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{table}{
%%     ~~Describe \code{table} here~~
}
  \item{title}{
%%     ~~Describe \code{title} here~~
}
  \item{end_time}{
%%     ~~Describe \code{end_time} here~~
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
function (table, title = "", end_time = 10)
{
    height <<- 1
    gap <<- 0.3
    ticks <- seq(0, end_time, by = 0.01)
    label_ticks <- seq(0, end_time, by = 1)
    labels <- as.character(label_ticks)
    plot(1, type = "n", xlim = c(0, end_time + 0.1), ylim = c(0,
        1), main = title, xlab = "Time", ylab = "", axes = FALSE)
    axis(1, at = label_ticks, labels = label_ticks)
    prepped_table <<- add_height(table)
    genPhyloInfect(prepped_table, end_time)
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
