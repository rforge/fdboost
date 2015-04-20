\name{NEWS}
\title{News for Package 'FDboost'}

\section{Changes in FDboost version 0.0-10 (2015-04-16)}{
  \subsection{User-visible changes}{
    \itemize{
    \item change the default in difference-penalties to first oder difference penalty
    differences=1, as then the effects are better identifiable 
      \item new method cvrisk.FDboost() that uses per default 
      sampling on the levels of curves, which is important for functional response
      \item reorganize documentation of cvrisk() and validateFDboost()
	  \item in bhist(): effect can be standardized 
    }
  }
  \subsection{Miscellaneous}{
    \itemize{
      \item add a CITATION file 
	  \item use mboost 2.4-2 as it exports all important functions 
    }
  }
  \subsection{Bug-fixes}{
    \itemize{
      \item main argument is always passed in plot.FDboost()
      \item bhist() and bconcurrent() work for equal time and s
      \item predict.FDboost() works with tensor-product base-learners bl1 \%X\% bl2
    }
  }
}


