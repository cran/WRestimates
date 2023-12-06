\name{wr.sigma.sqr}
\alias{wr.sigma.sqr}
\title{Assumed Population Variance of a Win Ratio}
\description{Calculate the assumed population variance of a win ratio.
\deqn{\sigma^2 = (4 * (1 + p[tie]))/(3 * k * (1 - k) * (1 - p[tie])}
Where;
\deqn{p[tie] = The proportion of ties.}
\deqn{k = The proportion of subjects allocated to one group.}
}
\usage{
wr.sigma.sqr(k, p.tie)
}
\arguments{
    \item{k}{The proportion of subjects allocated to one group i.e. the proportion of patients allocated to treatment.}
  \item{p.tie}{The proportion of ties.}
}
\value{
\code{wr.sigma.sqr} returns an object of \link{class} "\code{list}" containing the following components:
\item{sigma.sqr}{Population variance of the natural log (\eqn{ln}) of the win ratio.}
  \item{k}{The proportion of subjects allocated to one group.}
  \item{p.tie}{The proportion of ties.}
}
\references{
Yu, R. X. and Ganju, J. (2022). Sample size formula for a win ratio endpoint. \emph{Statistics in medicine}, \emph{41}(6), 950-963. doi: 10.1002/sim.9297.
}
\author{Autumn O'Donnell \email{autumn.research@gmail.com}}
\seealso{\code{\link{wr.var}}}
