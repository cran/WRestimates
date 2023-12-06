\name{wr.power}
\alias{wr.power}
\title{Power of a Win Ratio}
\description{Calculate the power of a win ratio.
\deqn{Power = 1 - \Phi(Z[\alpha] - ln(WR[true])(\sqrt{N}/\sigma))}
}
\usage{
wr.power(N, alpha = 0.025, WR.true = 1, sigma.sqr, k, p.tie)
}

\arguments{
  \item{N}{Sample size.}
  \item{alpha}{Level of significance (Type I error rate); Default: \eqn{\alpha} = 0.025.}
  \item{WR.true}{True or assumed win ratio; Default: \code{WR.true} = 1 where \out{H<sub>0</sub>} is assumed true.}
  \item{sigma.sqr}{Population variance of the natural log (\eqn{ln}) of the win ratio.}
  \item{k}{The proportion of subjects allocated to one group i.e. the proportion of patients allocated to treatment.}
  \item{p.tie}{The proportion of ties.}
}
\value{
\code{wr.power} returns an object of \link{class} "\code{list}" containing the following components:
\item{power}{Power of the win ratio.}
  \item{N}{Sample size.}
  \item{alpha}{Level of significance.}
  \item{WR.true}{True or assumed win ratio.}
  \item{sigma.sqr}{Population variance of the natural log (\eqn{ln}) of the win ratio.}
  \item{k}{The proportion of subjects allocated to one group.}
  \item{p.tie}{The proportion of ties.}
}
\references{
Yu, R. X. and Ganju, J. (2022). Sample size formula for a win ratio endpoint. \emph{Statistics in medicine}, \emph{41}(6), 950-963. doi: 10.1002/sim.9297.
}
\author{Autumn O'Donnell \email{autumn.research@gmail.com}}
\seealso{\code{\link{wr.sigma.sqr}}}
\examples{
## N = 100 patients, 1:1 allocation, one-sided alpha = 2.5\%, small
## proportion of ties p.tie = 0.1, and 50\% more wins on treatment
## than control.

### Calculate the Power
wr.power(N = 100, WR.true = 1.5, k = 0.5, p.tie = 0.1)
}
