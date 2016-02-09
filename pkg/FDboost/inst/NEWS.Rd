\name{NEWS}
\title{News for Package 'FDboost'}

\section{Changes in FDboost version 0.0-14 (2015-02-05)}{
  \subsection{User-visible changes}{
    \itemize{
	    \item new operator %Xc% for row tensor product where the sum-to-zero constraint is applied to
		  the design matrix resulting from the row-tensor product (experimental!)
	    \item allow newdata to be a list in predict.FDboost() in combination with signal base-learners 
        \item expand coef.FDboost() such that it works for 3-dimensional tensor products
          of with bhistx() the form bhistx() \%X\% bolsc() \%X\% bolsc() (with David Ruegamer)
        \item add a new possibility for scalar-on-function regression: 
          for timeformula=NULL, no Kronecker-product with 1 is used, which 
          changes the penalty as otherwise in the direction of 1 is penalized as well. 
    }
  }
  \subsection{Miscellaneous}{
    \itemize{
	  \item remove dependency on R-package MASS 
      \item use the argument 'prediction' in the internal computation 
         of the base-learners (work in progress) 
      \item throw an error if 'timeLab' of the hmatrix-object in bhistx() is not
         equal to the time-variable in 'timeformula'. 
    }
  }
}


\section{Changes in FDboost version 0.0-13 (2015-11-17)}{
  \subsection{User-visible changes}{
    \itemize{
	\item in function FDboost() the offset is supplied differently, for a scalar offset, use offset = "scalar", the default is still the same offset=NULL
	\item predict.FDboost() has new argument toFDboost (logical)
	\item fitted.FDboost() has argument toFDboost explicitly and not only in ...
    \item new base-learner bhistx() especially suited for effects with %X%, like bhistx %X% bolsc 
    \item coef.FDboost() and plot.FDboost() suited for effects like bhistx %X% bolsc
    \item for predict.FDboost() with effects bhistx() and newdata the latest mboostPatch is necessary
    }
  }
    \subsection{Bug-fixes}{
    \itemize{
      \item check for necessity of smooth offset works for missings in regular response (spotted by Tore Erdmann) 
    }
  }
}

\section{Changes in FDboost version 0.0-12 (2015-09-15)}{
    \itemize{
    \item Internal experimental version.
    } 
}

\section{Changes in FDboost version 0.0-11 (2015-06-01)}{
  \subsection{User-visible changes}{
    \itemize{
    \item integrationWeights() gives equal weights for regular grids
    \item new base-learner bfpc() for a functional covariate where
    functional covariate and the coeffcient are both expanded using fPCA (experimental feature!); 
    only works for regularly observed functional covariate. 
    }
  }
  \subsection{Bug-fixes}{
    \itemize{
      \item the function coef.FDboost() only works for bhist() if the time variable is the same in the timeformula and in bhist()
      \item predict.FDboost() has a check that for newdata only type="link" can be predicted
    }
  }
}

\section{Changes in FDboost version 0.0-10 (2015-04-16)}{
  \subsection{User-visible changes}{
    \itemize{
    \item change the default in difference-penalties to first order difference penalty
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


