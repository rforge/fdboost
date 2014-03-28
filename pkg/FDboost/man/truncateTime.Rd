\name{truncateTime}
\alias{truncateTime}
\title{Function to truncate time in functional data}
\usage{
  truncateTime(funVar, time, newtime, data)
}
\arguments{
  \item{funVar}{names of functional variables that should
  be truncated}

  \item{time}{name of time variable}

  \item{newtime}{new time vector that should be used. Must
  be part of the old time-line.}

  \item{data}{list containing all the data}
}
\value{
  A list with the data containing all variables of the
  original dataset with the variables of \code{funVar}
  truncated according to \code{newtime}.
}
\description{
  Function to truncate time in functional data
}
\note{
  All variables that are not part if \code{funVar}, or
  \code{time} are simply copied into the new data list
}

