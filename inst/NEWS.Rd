\name{NEWS}
\title{NEWS for WeMix}

\section{Changes in version 3.1.8}{
  \subsection{BUG FIXEX}{
    \itemize{
      \item \code{mix} now works with factor outcomes in a binomial model.
      \item \code{mix} now works with factor IDs.
      \item the final latent estimates were not converging because the Newton steps were not sufficently controlled. They are now more agressively prevented from diverging.
      \item when the variance estimate of a binomial model converged to the minimum variable (about 0.01) the model never acknowledged convergence, it now does.
      \item thanks to Tan Khai Ee and Prof. Taberez A. Neyazi for identifying and helping develop test code for the above four bugs.
      \item non-linear models summary function was not printing correctly. They are now also correctly showing weighted results, similar to linear models in 3.1.7.
      \item the warning for very low variance models in binomial regression was updated to clearly state the limitation of WeMix (it cannot estimate such models), and possible solutions.
    }
  }
}


\section{Changes in version 3.1.7}{
  \subsection{IMPROVEMENTS}{
    \itemize{
      \item with the extended information in \code{Groups} the warning about conditional weights was extranious and so was removed.
    }
  }
}

\section{Changes in version 3.1.7}{
  \subsection{IMPROVEMENTS}{
    \itemize{
      \item improved infomation about weights by level under \code{Groups} header of an \code{mml} summary.
      \item a conditioning variable cannot contain a dot, which now has a clearer error. Thank you to Felix Schulz for identifying this issue.
    }
  }
}

\section{Changes in version 3.1.6}{
  \subsection{IMPROVEMENTS}{
    \itemize{
      \item it is possible to pass an S3 object that claims to extend data frame but does not (such as a \code{tbl_df}) and so now all data objects are cast as a data frame before procesing beings. The specific issue is that a \code{tbl_df} defaults to \code{drop=FALSE} causing a single vector to maintain its name when it is expected to be unnamed.
    }
  }
}

\section{Changes in version 3.1.5}{
  \subsection{IMPROVEMENTS}{
    \itemize{
      \item Models now use less memory when there is a large number of units and groups after switching to a sparse Matrix for the random effects.
      \item The method of forming a generalized inverse was made more robust. It used to try only a QR decomposition based generalized inverse; this is now lest robust for unknown reasons. Now it also tries an SVD decomposition based generalized inverse.
      \item Updated tests to use EdSurvey 2.6.1.
    }
  }
}

\section{Changes in version 3.1.4}{
  \subsection{BUG FIXES}{
    \itemize{
      \item The vignettes should now appear as pdf files on CRAN.
    }
  }
}

\section{Changes in version 3.1.3}{
  \subsection{BUG FIXES}{
    \itemize{
      \item The method used to determine the rank of a matrix was brought into agreement with the \code{base} package when the matrix was a base matrix. Prior to this it was possible to get an error from a call to \code{mix} that should have returned. Now these calls to \code{mix} should return. Thanks to Christian Kjeldsen of Aarhus University, Danish School of Education, for pointing this issue out.
      \item The Wald test now works for generalized linear models instead of throwing an error. Thanks to Christian Kjeldsen for pointing this issue out as well.
    }
  }
}

\section{Changes in version 3.1.0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{WeMix} can now accept conditional weights. See the \code{cWeights} argument in the \code{mix} function.
      \item the \code{mix} function checks weights and writes a message if they may be conditional and \code{cWeights} is set to \code{FALSE}.
      \item Linear model evaluation is more robust and can handle data with non-invertible Z matrixes within a group.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Linear models now use \code{base::qr} more aggressively because of poor performance of the \code{Matrix::qr.coef} function on a sparse QR when the system is singular. This previously resulted in very large variance estimates. This also fixed an \code{invalid 'times' argument} error.
    }
  }
  \subsection{VIGNETTE}{
    \itemize{
      \item The code in the vignette was not the code used to generate the results and some values were incorrectly entered in the comparison table under the \code{mix} column. These problems were fixed.
    }
  }
}
\section{Changes in version 3.0.0}{
  \subsection{MODEL EVALUATION}{
    \itemize{
      \item Linear models are now solved using an analytical solution based on work by Bates and Pinheiro, (1998). This solution is significantly faster than the previous adaptive quadrature method. 
      \item Non-linear models are still evaluated using adaptive quadrature. 
    }
  }
  \subsection{NEW FEATURES}{
    \itemize{
      \item WeMix can now fit weighted three-level linear models, see the Weighted Linear Mixed-Effects Model vignette for details. Non-linear models are still evaluated using adaptive quadrature and are limited to two-level models. 
      \item Model evaluation is now possible using Wald tests. Wald tests allow users to test both fixed effects and random effects variances. 
    }
  }
}
\section{Changes in version 2.0.0}{
  \subsection{NEW MODEL TYPES}{
    \itemize{
      \item Supports binomial models
    }
  }
  \subsection{MODEL SPECIFICATION}{
    \itemize{
      \item Added ability to perform group and grand mean centering to increase comparability with Hierarchical Linear and Nonlinear Modeling (HLM) software
    }
  }
  \subsection{CODE}{
    \itemize{
      \item Although three-level models are not currently supported, in version 2.0.0, changes were made to the  way groups handled and to the data structures used for integration over random effects so as to be compatible with the future development of three-level models. 
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Corrected the warning message for the fast option (using Rcpp)
    }
  }
  \subsection{OTHER CHANGES}{
    \itemize{
      \item \code{fast} option in \code{mix} defaults to \code{FALSE} now to prioritize accuracy over speed.
    }
  }
}