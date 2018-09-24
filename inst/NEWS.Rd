\name{NEWS}
\title{NEWS for WeMix}
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