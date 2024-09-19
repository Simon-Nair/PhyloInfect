\name{end_time_infected}
\alias{end_time_infected}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
EndTimeInfected
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
end_time_infected(table)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{table}{
%%     ~~Describe \code{table} here~~
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
function (table)
{
    I_list <- c()
    for (i in 1:nrow(table)) {
        if (table$infector[i] \%in\% I_list && table$affected[i] ==
            0) {
            I_list <- I_list[!I_list \%in\% table$infector[i]]
        }
        else {
            I_list <- c(I_list, table$affected[i])
        }
    }
    return(I_list)
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
