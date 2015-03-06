\docType{package}
\name{SelvarMix-package}
\alias{SelvarMix-package}
\title{Regularization for variable selection in model-based clustering and discriminant analysis}
\description{
  SelvarMix is a package where a regularization approach 
  of variable selection is considered in model-based 
  clustering and discriminant analysis frameworks. 
  First, this procedure consists of ranking the variables 
  with a lasso-like procedure. Second, the method of
  Maugis et al (2009, 2011) is adapted to define the 
  role of variables in the two frameworks. 
  SelvarMix provides a faster variable selection 
  algorithm than the backward stepwise or forward 
  stepwise algorithms of Maugis et al (2009), allowing us 
  to study high-dimensional datasets.
}
\details{
  \tabular{ll}{ Package: \tab SelvarMix\cr Type: \tab
  Package\cr Version: \tab 1.0\cr Date: \tab
  2014-04-03\cr License: \tab GPL-3 + file LICENSE\cr
  LazyLoad: \tab yes\cr }

  The general purpose of the package is to perform variable 
  selection in model-based clustering and
  discriminant analysis. It focus on model-based 
  clustering, where the clusters are assumed to 
  arise from Gaussian distributions. 
  The most achieved model in model-based clustering
  has been proposed by Maugis et al (2009). This 
  so-called \eqn{SRUW} modeling considers three roles
  of variables: one variable my belong to the relevant 
  clustering set \eqn{S}, the redundant variable set
  \eqn{U} or the independent variable set \eqn{W}. 
  Moreover, the redundant variables may be explained 
  by a subset \eqn{R} of the relevant variables \eqn{S}.
  In order to avoid the slow of this algorithm when data 
  with numerous variables are studied, the SelvarMix procedure
  is proposed. It proceeds in two steps: First, the variables
  are ranked using a lasso-like procedure analogous to the one of 
  Zhou et al (2009); second, the \eqn{SRUW} procedure is run 
  on this  ranked set of  variables. 
  

}
\examples{
\dontrun{
## Simulated data  example as shown in Maugis et al. (2009) (correlated scenario 2) 
## n = 2000 observations, p = 14 variables 
require(Rmixmod)
require(glasso)
data(scenarioCor)
data.cor <- scenarioCor[,1:14]
labels.cor <-scenarioCor[,15]


lambda <- seq(20,  100, by = 10)
rho <- seq(1, 2, length=2)
hybrid.size <- 3
models <- mixmodGaussianModel(family = "spherical", equal.proportions = TRUE)
regModel <- c("LI","LB","LC")
indepModel <- c("LI","LB")

## variable selection in model-based clustering  
nbCluster <-  c(3,4)
criterion <- "BIC"
simulate.cl  <- SelvarClustLasso(data.cor, nbCluster, lambda, rho, hybrid.size, 
                                 criterion, models, regModel, indepModel)



## variables selection in discriminant analysis
## training sample : n = 1900, p = 14 variables 
data.learn <- scenarioCor[1:1900,1:14]
labels.learn <-scenarioCor[1:1900,15]

## testing sample : n = 100, p = 14 variables
data.test <- scenarioCor[1901:2000,1:14]
labels.test <-scenarioCor[1901:2000,15]

lambda <- seq(20,  50, length = 10)

simulate.da <- SelvarLearnLasso(data.learn, labels.learn, lambda, rho, hybrid.size, 
                                models, regModel, indepModel, data.test, labels.test)
}
}
\author{
  Author: Mohammed Sedki, Gilles Celeux and Cathy Maugis-Rabusseau  
  
}
\references{
  Maugis, C., Celeux, G., and Martin-Magniette, M. L., 2009. 
  "Variable selection in model-based clustering: 
  A general variable role modeling". Computational 
  Statistics and Data Analysis, vol. 53/11, pp. 3872-3882.

  Maugis, C., Celeux, G., and Martin-Magniette, M. L., 2011. 
  "Variable selection in model-based discriminant analysis". 
  Journal of Multivariate Analysis, vol. 102, pp. 1374-1387.

  Zhou, H., Pan, W., and Shen, X., 2009. "Penalized model-based 
  clustering with unconstrained covariance matrices". 
  Electronic Journal of Statistics, vol. 3, pp.1473-1496.
  
  Sedki, M., Celeux, G., Maugis-Rabusseau, C., 2014. 
  "SelvarMix: A R package for variable selection in 
  model-based clustering and discriminant analysis with 
  a regularization approach". Inria Research Report 
  available at \url{http://hal.inria.fr/hal-01053784}

}
\keyword{package}

