\name{sublist}
\alias{sublist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sublist
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sublist(table, I_list, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{table}{
%%     ~~Describe \code{table} here~~
}
  \item{I_list}{
%%     ~~Describe \code{I_list} here~~
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
function (table, I_list, n)
{
    sub_Ilist <- sample(I_list, n)
    sub_Ilist <- unlist(sub_Ilist)
    sampled_indices <- c()
    for (i in sub_Ilist) {
        sampled_indices <- c(sampled_indices, which(table$affected ==
            i))
    }
    sampled_heights <- table$heights[sampled_indices]
    sampled_df <- data.frame(indices = sub_Ilist, heights = sampled_heights)
    sampled_df <- sampled_df[order(-sampled_df$heights), ]
    sorted_sampled_indices <- sampled_df$indices
    return(sorted_sampled_indices)
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
