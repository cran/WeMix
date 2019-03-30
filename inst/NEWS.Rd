\name{NEWS}
\title{NEWS for WeMix}
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