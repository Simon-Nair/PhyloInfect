\name{genPhyloInfect}
\alias{genPhyloInfect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GenPhyloInfect
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
genPhyloInfect(table, end_time = 10, recent = table$affected[1], time_of_infection = 0, counter = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{table}{
%%     ~~Describe \code{table} here~~
}
  \item{end_time}{
%%     ~~Describe \code{end_time} here~~
}
  \item{recent}{
%%     ~~Describe \code{recent} here~~
}
  \item{time_of_infection}{
%%     ~~Describe \code{time_of_infection} here~~
}
  \item{counter}{
%%     ~~Describe \code{counter} here~~
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
function (table, end_time = 10, recent = table$affected[1], time_of_infection = 0,
    counter = 0)
{
    indices <- rev(which(as.numeric(unlist(table$infector)) ==
        as.numeric(unlist(recent))))
    segment_count <- 0
    if (length(indices) > 0 && table$affected[indices[1]] ==
        0) {
        recovery_time <- table$time[indices[1]]
        segments(time_of_infection, height, recovery_time, height,
            lwd = 1)
        segments(time_of_infection, height, time_of_infection,
            height + 0.016 * counter, lwd = 1, )
        text(recovery_time + gap, height + 0.003, labels = recent,
            cex = 0.55)
        indices <- indices[-1]
        segment_count <- segment_count + 1
    }
    else {
        segments(time_of_infection, height, end_time, height,
            lwd = 1)
        segments(time_of_infection, height, time_of_infection,
            height + 0.016 * counter, lwd = 1, )
        text(end_time + gap, height + 0.003, labels = recent,
            cex = 0.55)
        segment_count <- segment_count + 1
    }
    for (i in indices) {
        if (gap == 0.3) {
            gap <<- 0.1
        }
        else {
            gap <<- 0.3
        }
        height <<- height - 0.016
        prepped_table$heights[i] <<- height
        segment_count <- segment_count + genPhyloInfect(table,
            end_time, table$affected[i], table$time[i], segment_count)
    }
    return(segment_count)
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
