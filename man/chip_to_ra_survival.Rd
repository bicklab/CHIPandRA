% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/do_survival.R
\name{chip_to_ra_survival}
\alias{chip_to_ra_survival}
\title{CHIP to RA survival analysis}
\usage{
chip_to_ra_survival(
  bl_data,
  outcomes,
  min_num_events = 100,
  min_num_events_w_chip = 10,
  ra_types = c("ra", "spra", "snra"),
  sensspecs = c("sensitive", "moderate", "specific"),
  chip_types = names(select(bl_data, starts_with("has_chip"))),
  debug = FALSE
)
}
\arguments{
\item{bl_data}{baseline data}

\item{outcomes}{outcome data}

\item{min_num_events}{mnimum number of events (outcomes) needed to be observed to continue with analysis}

\item{min_num_events_w_chip}{minimum number of events (outcomes) needed to be onbserved among people with CHIP to continue analysis}

\item{ra_types}{types of RA to study}

\item{sensspecs}{spectrum of sensitive to specific to study}

\item{chip_types}{types of CHIP to study}

\item{debug}{flag to debug vs run normally}
}
\value{
results tibble
}
\description{
CHIP to RA survival analysis
}
