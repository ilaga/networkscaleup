#' Simulated ARD data set with z and x.
#'
#' A simulated data set to demonstrate and test the NSUM methods. The data was
#' simulated from the basic Killworth Binomial model.
#'
#' @format A named list for an ARD survey from 100 respondents about 5
#'   subpopulations. \describe{ \item{ard}{A `100 x 5` matrix with integer
#'   valued respondents} \item{x}{A `100 x 5` matrix with simulated answers from
#'   a 1-5 Likert scale} \item{z}{A `100 x 4` matrix with answers for each
#'   respondents about 4 demographic questions} \item{N}{An integer specifying
#'   the total population size} \item{subpop_size}{A vector with the 5 true
#'   subpopulation sizes} \item{degrees}{A vector with the 100 true respondent
#'   degrees}}
"example_data"
