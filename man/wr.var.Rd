\name{wr.var}
\alias{wr.var}
\title{Approximate Variance of the Natural Log (\eqn{ln}) of the Win Ratio.}
\description{
Calculating the approximate variance of the natural log (\eqn{ln}) a win ratio.
\deqn{Var(ln(WR)) ~~ \sigma^2/N}
Where;
\deqn{\sigma^2 = (4 * (1 + p[tie]))/(3 * k * (1 - k) * (1 - p[tie])}
}
\usage{
wr.var(N, sigma.sqr, k, p.tie)
}
\arguments{
  \item{N}{Sample size.}
  \item{sigma.sqr}{Population variance of the natural log (\eqn{ln}) of the win ratio.}
  \item{k}{The proportion of subjects allocated to one group i.e. the proportion of patients allocated to treatment.}
  \item{p.tie}{The proportion of ties.}
}
\value{
\code{wr.var} returns an object of \link{class} "\code{list}" containing the following components:
\item{var.ln.WR}{Approximate variance of the natural log (\eqn{ln}) a win ratio.}
  \item{N}{Sample size.}
  \item{sigma.sqr}{Population variance of the natural log (\eqn{ln}) of the win ratio.}
  \item{k}{The proportion of subjects allocated to one group.}
  \item{p.tie}{The proportion of ties.}
}
\references{
Yu, R. X. and Ganju, J. (2022). Sample size formula for a win ratio endpoint. \emph{Statistics in medicine}, \emph{41}(6), 950-963. doi: 10.1002/sim.9297.
}
\author{Autumn O'Donnell \email{autumn.research@gmail.com}}
\seealso{
\code{\link{wr.sigma.sqr}}
}
