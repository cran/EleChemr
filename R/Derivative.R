#' Derivative calculation of concentration profile
#'
#' Return a the derivative of the concentration profile simulated
#'
#' @param npoints number of points to be used for the derivative
#' @param h space for the finite difference
#' @param Ox data upon the derivative is calculated
#' @param mode "Forward" or "Backward" the derivative will be calculated for the npoints
#' @param Derivative "First" or "Second" derivative to calculate
#'
#'
#' @return a vector with the derivative requested
#'
#' @examples
#' Derv(npoints = 2, h = 0.13, Ox = matrix(c(1,2), nrow = 1), mode = "Forward", Derivative = "First")
#'
#' @export
#' @importFrom pracma inv
#'

Derv = function(npoints = 2, h, Ox, mode = "Forward", Derivative = "First") {

  if (mode == "Forward") {


    Taylor = c()
    for (i1 in 1:(npoints-1)){
      Taylor = append(Taylor, (1/factorial(i1)))
    }

    A = matrix(data = rep(Taylor, npoints-1), byrow = T, nrow= npoints-1)
    if (npoints > 2){
      for (i1 in 2:npoints){
        for (j1 in 1:ncol(A)){
          A[i1-1,j1] = A[i1-1,j1]*(i1^j1)
        }
      }
    }
    Ainv = inv(A)
    Deriv = matrix(nrow = nrow(Ox), ncol=1)

    for (x2 in 1:nrow(Ox)) {
      b = c()
      for (x1 in 2:npoints){
        b = append(b, (Ox[x2,x1] - Ox[x2,1]))
      }
      if (Derivative == "First") {
        Deriv[x2,1] = (1/h)*Ainv[1,]%*%b
      } else if (Derivative == "Second") {
        Deriv[x2,1] = (1/(h^2))*Ainv[2,]%*%b
      }

    }

    return(Deriv)


  } else if (mode == "Backward") {

    Taylor = c()
    for (i1 in 1:(npoints-1)){
      Taylor = append(Taylor, (-1^i1)*(1/factorial(i1)))
    }

    A = matrix(data = rep(Taylor, npoints-1), byrow = T, nrow= npoints-1)
    if (npoints > 2){
      for (i1 in 2:npoints){
        for (j1 in 1:ncol(A)){
          A[i1-1,j1] = A[i1-1,j1]*(i1^j1)
        }
      }
    }
    Ainv = inv(A)
    Deriv = matrix(nrow = nrow(Ox), ncol=1)

    for (x2 in 1:nrow(Ox)) {
      b = c()
      for (x1 in 2:npoints){
        b = append(b, (Ox[x2,1] - Ox[x2,x1]))
      }
      if (Derivative == "First") {
        Deriv[x2,1] = (1/h)*Ainv[1,]%*%b
      } else if (Derivative == "Second") {
        Deriv[x2,1] = (1/(h^2))*Ainv[2,]%*%b
      }
    }

    return(Deriv)
  }
}
