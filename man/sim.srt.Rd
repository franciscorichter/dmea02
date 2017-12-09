\name{sim.srt}
\alias{sim.srt}
Simulate a set of ...
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sim.srt(tree, pars, parallel = F, n_trees, WT = FALSE, method = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
%%     ~~Describe \code{tree} here~~
}
  \item{pars}{
%%     ~~Describe \code{pars} here~~
}
  \item{parallel}{
%%     ~~Describe \code{parallel} here~~
}
  \item{n_trees}{
%%     ~~Describe \code{n_trees} here~~
}
  \item{WT}{
%%     ~~Describe \code{WT} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
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
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (tree, pars, parallel = F, n_trees, WT = FALSE, method = 1)
{
    if (parallel) {
        no_cores <- detectCores() - 1
        cl <- makeCluster(no_cores)
        registerDoParallel(cl)
        trees <- foreach(i = 1:n_trees, combine = list) \%dopar\%
            dmea02::rec.tree(tree = tree, pars = pars)
        stopCluster(cl)
    }
    else {
        trees = vector("list", length = n_trees)
        for (i in 1:n_trees) {
            if (method == 1)
                rec = rec.tree(tree = tree, pars = pars)
            if (method == 2)
                rec = rec.tree2(tree = tree, pars = pars)
            if (WT) {
                trees[[i]] = rec$wt
            }
            else {
                trees[[i]] = rec
            }
        }
    }
    return(trees)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
