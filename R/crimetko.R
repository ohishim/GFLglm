#' The number of recognized crimes committed in the Tokyo area
#'
#' The `crimetko` is a grouped data with 9570 rows and 4 variables.
#' The number of groups is 318 which consists of 6 years and 53 areas.
#'
#' @format
#' \describe{
#'   \item{year}{year}
#'   \item{area}{area id}
#'   \item{crime}{the number of crimes}
#'   \item{pop}{population}
#'   \item{group}{group id}
#' }
#' @source
#' \describe{
#'   \item{`crime`}{the original data was collected by the Metropolitan Police Department,
#'     available at TOKYO OPEN DATA (\url{https://portal.data.metro.tokyo.lg.jp/});
#'     this is arranged and used the following production: Tokyo Metropolitan
#'     Government & Metropolitan Police Department. The number of recognized
#'     cases by region, crime type, and method (yearly total; in Japanese),
#'     \url{https://creativecommons.org/licenses/by/4.0/deed.en.}}
#'   \item{`pop`}{the original data was obtained from the results of the population census,
#'     as provided in e-Stat (\url{https://www.e-stat.go.jp/en}); this is arranged and used}
#' }
"crimetko"
