\name{SelvarClustLasso}
\alias{SelvarClustLasso}
\title{
Regularization for variable selection in model-based clustering
}
\description{
This function implements the variable selection in model-based clustering 
using a lasso ranking on the variables as described in Sedki et al (2014).
The variable ranking step uses the penalized EM algorithm of Zhou et al (2009).
}
\usage{
SelvarClustLasso(data, nbCluster, lambda, rho, hybrid.size, criterion, 
                 models, regModel, indepModel, nbCores)
}

\arguments{

  \item{data}{
  matrix  containing quantitative data. 
  Rows correspond to observations and
  columns correspond to variables
}
  \item{nbCluster}{
 numeric listing of the number of clusters (must be positive integers)
}
  \item{lambda}{
  numeric listing of the tuning parameter for \eqn{\ell_1} mean penalty
}
  \item{rho}{
  numeric listing of the tuning parameter for \eqn{\ell_1} precision matrix penalty
}
\item{hybrid.size}{
optional parameter make less strength  the hybrid forward and backward 
algorithms to select \eqn{S} and \eqn{W} sets
}

\item{criterion}{list of character defining the criterion
  to select the best model. The best model is the one with
  the highest criterion value. Possible values: "BIC",
  "ICL", c("BIC", "ICL"). Default is "BIC"}

\item{models}{a  Rmixmod [\code{\linkS4class{Model}}] object
  defining the list of models to run. The models 
  Gaussian_pk_L_C, Gaussian_pk_Lk_C, Gaussian_pk_L_Ck, 
  and Gaussian_pk_Lk_Ck are called by default 
  (see mixmodGaussianModel() in Rmixmod package to
  specify other models)}

\item{regModel}{list of character defining the covariance matrix form for
  the linear regression of \eqn{U} on the \eqn{R} set of variables. 
  Possible values: "LI" for spherical form,
  "LB" for diagonal form and  "LC" for general form.
  Possible values: "LI", "LB", "LC", c("LI", "LB"), 
  c("LI", "LC"), c("LB", "LC") and c("LI", "LB", "LC"). 
  Default is c("LI", "LB", "LC")}

\item{indepModel}{list of character defining the covariance matrix form for
  independent variables \eqn{W}. Possible values: "LI" for spherical 
  form and "LB" for diagonal form. Possible values: "LI", "LB", c("LI", "LB"). 
  Default is c("LI", LB")}

\item{nbCores}{
  number of CPUs to be used when parallel computing is utilized (default is 2)
  }
}

\value{ for each criterion BIC or ICL
\item{S }{The selected set of relevant clustering variables}
\item{R }{The selected subset of regressors}
\item{U }{The selected set of redundant variables}
\item{W }{The selected set of independent variables}
\item{criterionValue}{The criterion value for the selected model}
\item{nbCluster}{The selected number of clusters}
\item{model}{The selected Gaussian mixture form} %% the selected gaussian mixture form
\item{regModel }{The  selected covariance form for the regression}
\item{indepModel}{The selected covariance form for the independent gaussian distribution}
\item{proba }{Matrix containing the conditional probabilities of belonging to each cluster for all observations}
\item{partition}{Vector of length \emph{n} containing the cluster assignments of the \emph{n} observations according to the Maximum-a-Posteriori rule}
}


\author{
Mohammed Sedki <\url{mohammed.sedki@u-psud.fr}>
}
\references{
  Zhou, H., Pan, W., and Shen, X., 2009. "Penalized model-based 
  clustering with unconstrained covariance matrices". 
  Electronic Journal of Statistics, vol. 3, pp.1473-1496.
  
  Maugis, C., Celeux, G., and Martin-Magniette, M. L., 2009. 
  "Variable selection in model-based clustering: 
  A general variable role modeling". Computational 
  Statistics and Data Analysis, vol. 53/11, pp. 3872-3882.
  

  Sedki, M., Celeux, G., Maugis-Rabusseau, C., 2014. 
  "SelvarMix: A R package for variable selection in 
  model-based clustering and discriminant analysis with 
  a regularization approach". Inria Research Report 
  available at \url{http://hal.inria.fr/hal-01053784}
}

\keyword{model-based clustering, discriminant analysis, variable selection, lasso ranking and graphical lasso}
\seealso{
\link{SelvarLearnLasso}
\link{SortvarClust}
\link{SortvarLearn}
\link{scenarioCor}
}
\examples{
\dontrun{
## Simulated data  example as shown in Maugis et al. (2009) 
## n = 2000 observations, p = 14 variables 
require(Rmixmod)
require(glasso)
data(scenarioCor)
data.cor <- scenarioCor[,1:14]

lambda <- seq(20,  100, by = 10)
rho <- seq(1, 2, length=2)
hybrid.size <- 3
nbCluster <-  c(3,4)
criterion <- "BIC"
models <- mixmodGaussianModel(family = "spherical", equal.proportions = TRUE)
regModel <- c("LI","LB","LC")
indepModel <- c("LI","LB")


simulate.cl <- SelvarClustLasso(data.cor, nbCluster, lambda, rho, hybrid.size, 
                                criterion, models, regModel, indepModel)
 
} 
}
