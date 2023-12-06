\name{wr.ci}
\alias{wr.ci}
\title{Confidence Interval (CI) for Win Ratio}
\description{
Calculate the confidence interval for a win ratio.
\deqn{CI = exp((ln(WR) +/- Z\sqrt{var})}
Where;

\eqn{ln(WR)} = Natural log of the true or assumed win ratio.

\eqn{Z} = Z-score from normal distribution.

\eqn{\sqrt{var}} = Standard deviation of the natural log of the win ratio.
}
\usage{
wr.ci(WR = 1, Z = 1.96, var.ln.WR, N, sigma.sqr, k, p.tie)
}
\arguments{
  \item{WR}{Win ratio; Default: \code{WR} = 1 for an assumed true win ratio where \out{H<sub>0</sub>} is assumed true.}
  \item{Z}{Z-score from normal distribution; Default: \code{Z} = 1.96 for a 95\% CI.}
  \item{var.ln.WR}{Variance of the natural log (\eqn{ln}) of the win ratio.}
  \item{N}{Sample size.}
  \item{sigma.sqr}{Population variance of the natural log (\eqn{ln}) of the win ratio.}
  \item{k}{The proportion of subjects allocated to one group i.e. the proportion of patients allocated to treatment.}
  \item{p.tie}{The proportion of ties.}
}

\value{
\code{wr.ci} returns an object of \link{class} "\code{list}" containing the following components:
\item{ci}{The confidence interval of a win ratio.}
  \item{WR}{The win ratio.}
  \item{Z}{Z-score from normal distribution.}
  \item{var.ln.WR}{Variance of the natural log (\eqn{ln}) of the win ratio.}
  \item{N}{Sample size.}
  \item{sigma.sqr}{Population variance of the natural log (\eqn{ln}) of the win ratio.}
  \item{k}{The proportion of subjects allocated to one group.}
  \item{p.tie}{The proportion of ties.}
}
\references{
Yu, R. X. and Ganju, J. (2022). Sample size formula for a win ratio endpoint. \emph{Statistics in medicine}, \emph{41}(6), 950-963. doi:10.1002/sim.9297.
}
\author{Autumn O'Donnell \email{autumn.research@gmail.com}}
\seealso{\code{\link{wr.sigma.sqr}}; \code{\link{wr.var}}}
\examples{
## N = 100 patients, 1:1 allocation, one-sided alpha = 2.5\%, power = 90\%
## (beta = 10\%), a small proportion of ties p.tie = 0.1, and 50\% more wins
## on treatment than control.

### Calculation 95\% CI
wr.ci(N = 100, WR = 1.5, k = 0.5, p.tie = 0.1)
}
