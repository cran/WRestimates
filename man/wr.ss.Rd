\name{wr.ss}
\alias{wr.ss}


\title{Approximate Sample Size of a Win Ratio}
\description{Calculates the approximate required sample size of a win ratio.
\deqn{N ~~ (\sigma^2 * (Z[1-\alpha] + Z[1-\beta])^2)/(ln^2(WR[true]))
}}
\usage{
wr.ss(alpha = 0.025, beta = 0.1, WR.true = 1, k, p.tie, sigma.sqr)
}

\arguments{
  \item{alpha}{Level of significance (Type I error rate); Default: \eqn{\alpha} = 0.025.}
  \item{beta}{Type II error rate; Default: \eqn{\beta} = 0.1.}
  \item{WR.true}{True or assumed win ratio; Default: \code{WR.true} = 1 where \out{H<sub>0</sub>} is assumed true.}
  \item{k}{The proportion of subjects allocated to one group i.e. the proportion of patients allocated to treatment.}
  \item{p.tie}{The proportion of ties.}
  \item{sigma.sqr}{Population variance of the natural log (\eqn{ln}) of the win ratio.}
}
\value{\code{wr.ss} returns an object of \link{class} "\code{list}" containing the following components:
\item{N}{Sample size.}
\item{alpha}{Level of significance (Type I error rate).}
  \item{beta}{Type II error rate.}
  \item{WR.true}{True or assumed win ratio.}
  \item{k}{The proportion of subjects allocated to one group.}
  \item{p.tie}{The proportion of ties.}
  \item{sigma.sqr}{Population variance of the natural log (\eqn{ln}) of the win ratio.}
}
\references{
Yu, R. X. and Ganju, J. (2022). Sample size formula for a win ratio endpoint. \emph{Statistics in medicine}, \emph{41}(6), 950-963. doi: 10.1002/sim.9297.
}
\author{Autumn O'Donnell \email{autumn.research@gmail.com}}
\seealso{\code{\link{wr.sigma.sqr}}}
\examples{
## 1:1 allocation, one-sided alpha = 2.5\%, power = 90\% (beta = 10\%),
## a small proportion of ties p.tie = 0.1, and 50\% more wins on treatment
## than control

### Calculate Sample Size
wr.ss(WR.true = 1.5, k = 0.5, p.tie = 0.1)
}
