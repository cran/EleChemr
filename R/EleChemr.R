#' Chrono amperometry digital simulation
#'
#' Return a graph I vs t of the electrochemical process
#'
#' @param Co bulk concentration expressed in Molar
#' @param exptime experimental time to be simulated expressed in seconds
#' @param Dx diffusion coefficient expressed in cm^2/s
#' @param Dm simulation parameter, maximum 0.5 for explicit methods
#' @param Temp temperature in kelvin
#' @param n number of electrons involved in the process
#' @param Area area of the electrode expressed in cm^2
#' @param l number of time steps of the simulation
#' @param DerApprox number of point for the approximation of the first derivative
#' @param errCheck if true the function returns a list with parameters for CottrCheck function
#' @param Method method to be used for the simulation = "Euler" "BI" "RK4" "CN" "BDF"
#'
#'
#' @return if errCheck == F a graph I vs t, if errCheck == T a list
#'
#' @examples
#' ChronAmp(Co = 0.001, exptime = 1, DerApprox = 2, Dm = 0.45, errCheck = FALSE, Method = "Euler")
#'
#' @export
#' @import ggplot2

ChronAmp = function(Co = 0.001, exptime = 1, Dx = 0.00001, Dm = 0.45,
                    Temp = 298.15, n = 1, Area = 1, DerApprox = 2, l = 100,
                    errCheck = FALSE, Method = "Euler") {

  Par = ParCall("ChronAmp", n. = n, Temp. = Temp, Dx1. = Dx, exptime. = exptime, Dm. = Dm, l. = l)
  Ox = OneMat(Par$l, Par$j)
  Jox = ZeroMat(Par$l, 1)

  if (Method == "Euler") {

    for (i1 in 1:(Par$l-1)) {

      Ox[i1,1] = 0

      for (j1 in 2:(Par$j-1)) {
        Ox[i1 + 1,j1] = Ox[i1,j1] + Dm*(Ox[i1, j1 -1] - 2*Ox[i1, j1] + Ox[i1, j1+1])
      }
    }

    Jox = Derv(Ox = Ox, h = Par$h, npoints = DerApprox)

  } else if (Method == "RK4") {

    for (i1 in 1:(Par$l-1)) {
      k1 = ZeroMat(Par$j)
      k2 = ZeroMat(Par$j)
      k3 = ZeroMat(Par$j)
      k4 = ZeroMat(Par$j)
      Ox[i1,1] = 0
      Ox[i1 +1, 1] = 0

      for (j1 in 2:(Par$j-1)) {
        k1[j1] = Dm*(Ox[i1, j1 -1] - 2*Ox[i1, j1] + Ox[i1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + k1[j1]*0.5
      }
      for (j1 in 2:(Par$j-1)) {
        k2[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + k2[j1]*0.5
      }
      for (j1 in 2:(Par$j-1)) {
        k3[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + k3[j1]
      }
      for (j1 in 2:(Par$j-1)) {
        k4[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + (k1[j1] + 2*k2[j1] + 2*k3[j1] + k4[j1])/6
      }
    }

    Jox = Derv(Ox = Ox, h = Par$h, npoints = DerApprox)

  } else if (Method == "BI") {
    al1 = 1/(Par$h^2)
    al2 = -2/(Par$h^2)
    al3 = 1/(Par$h^2)
    a1 = (al2 - 1/Par$dtn)/al1
    a2 = al3/al1

    for (i1 in 1:(Par$l-1)) {
      Y = ZeroMat(Par$j-2,Par$j-2)
      Y[1,1] = a1
      Y[1,2] = a2
      Y[Par$j-2,Par$j-3] = 1
      Y[Par$j-2,Par$j-2] = a1
      for (i in 2:(Par$j-3)) {
        Y[i,i] = a1
        Y[i,i-1] = 1
        Y[i, i +1] = a2
      }

      Ox[i1,1] = 0
      Ox[i1+1,1] = 0
      b = (-Ox[i1,2:(Par$j-1)]/(al1*Par$dtn))
      b[Par$j-2] = b[Par$j-2] - a2*1
      b[1] = b[1] - Ox[i1+1,1]
      Ox[i1+1,2:(Par$j-1)] = solve(Y) %*% b
    }

    Jox = Derv(Ox = Ox, h = Par$h, npoints = DerApprox)

  } else if (Method == "CN") {

    al1 = 1/(Par$h^2)
    al2 = -2/(Par$h^2)
    al3 = 1/(Par$h^2)
    a1 = (al2 - 2/Par$dtn)/al1
    a2 = al3/al1
    a3 = (al2 + 2/Par$dtn)/al1

    for (i1 in 1:(Par$l-1)) {
      Y = ZeroMat(Par$j-2,Par$j-2)
      Y[1,1] = a1
      Y[1,2] = a2
      Y[Par$j-2,Par$j-3] = 1
      Y[Par$j-2,Par$j-2] = a1
      for (i in 2:(Par$j-3)) {
        Y[i,i] = a1
        Y[i,i-1] = 1
        Y[i, i +1] = a2
      }

      Ox[i1,1] = 0
      Ox[i1+1,1] = 0
      b = -a3*Ox[i1,2:(Par$j-1)] - Ox[i1,1:((Par$j-1)-1)] - a2*Ox[i1,3:((Par$j-1)+1)]
      b[Par$j-2] = b[Par$j-2] - a2*1
      b[1] = b[1] - Ox[i1+1,1]
      Ox[i1+1,2:(Par$j-1)] = solve(Y) %*% b
    }

    Jox = Derv(Ox = Ox, h = Par$h, npoints = DerApprox)


  } else if (Method == "BDF") {

    al1 = 1/(Par$h^2)
    al2 = -2/(Par$h^2)
    al3 = 1/(Par$h^2)
    a1 = (al2 - 1.5/Par$dtn)/al1
    a2 = al3/al1

    for (i1 in 1:(Par$l-1)) {
      Y = ZeroMat(Par$j-2,Par$j-2)
      Y[1,1] = a1
      Y[1,2] = a2
      Y[Par$j-2,Par$j-3] = 1
      Y[Par$j-2,Par$j-2] = a1
      for (i in 2:(Par$j-3)) {
        Y[i,i] = a1
        Y[i,i-1] = 1
        Y[i, i +1] = a2
      }

      Ox[i1,1] = 0
      Ox[i1+1,1] = 0
      if (i1 == 1) {
        b = -2*Ox[i1,2:(Par$j-1)]/(Par$dtn*al1) + Ox[1,2:(Par$j-1)]/(2*Par$dtn*al1)
      } else {
        b = -2*Ox[i1,2:(Par$j-1)]/(Par$dtn*al1) + Ox[i1-1,2:(Par$j-1)]/(2*Par$dtn*al1)
      }
      b[Par$j-2] = b[Par$j-2] - a2*1
      b[1] = b[1] - Ox[i1+1,1]
      Ox[i1+1,2:(Par$j-1)] = solve(Y) %*% b
    }

    Jox = Derv(Ox = Ox, h = Par$h, npoints = DerApprox)


  } else if (!(Method %in% c("Euler", "BI", "RK4", "CN", "BDF"))) {
    return("Available methods are Euler, BI, RK4, CN and BDF")
  }

  G = Jox
  i = (n*Par$FA*G*Dx*Area*Co)/(sqrt(Dx*Par$tau))

  graphy = ggplot(data = data.frame(i[1:(length(i)-1)],Par$t[1:(length(i)-1)]),
                  aes(y = i[1:(length(i)-1)], x = Par$t[1:(length(i)-1)])) +
    geom_point() + xlab("t / s") +
    ylab("I / A") + theme_classic()

  if (errCheck == TRUE){
    return(list(G,Dx,Co,Par$dtn,Par$h,Par$l,Par$j,i,n,Area))
  } else {
    return(graphy)
  }
}


#' Chrono amperometry with a finite step digital simulation
#'
#' Return a graph I vs t of the electrochemical process
#'
#' @param Co bulk concentration expressed in Molar
#' @param exptime experimental time to be simulated expressed in seconds
#' @param Dx diffusion coefficient expressed in cm^2/s
#' @param eta overpotential of the step expressed in Volt
#' @param Dm simulation parameter, maximum 0.5 for explicit methods
#' @param Temp temperature in kelvin
#' @param n number of electrons involved in the process
#' @param Area area of the electrode expressed in cm^2
#' @param l number of time steps of the simulation
#' @param DerApprox number of point for the approximation of the first derivative
#' @param errCheck if true the function returns a list with parameters for CottrCheck function
#' @param Method method to be used for the simulation = "Euler" "BI" "RK4" "CN" "BDF"
#'
#'
#' @return if errCheck == F a graph I vs t, if errCheck == T a list
#'
#' @examples
#' PotStep(Co = 0.001, exptime = 1, Dm =0.45, DerApprox = 2, errCheck = FALSE, Method = "Euler")
#'
#' @export
#' @import ggplot2


PotStep = function(Co = 0.001, exptime = 1, Dx = 0.00001, Dm = 0.45,
                   eta = 0, Temp = 298.15, n = 1, Area = 1, l= 100,
                   DerApprox = 2, errCheck = FALSE, Method = "Euler") {

  Par = ParCall("PotStep", n. = n, Temp. = Temp, Dx1. = Dx, exptime. = exptime,
                Dm. = Dm, eta. = eta, l. = l)
  Ox = OneMat(Par$l, Par$j)
  Red = ZeroMat(Par$l, Par$j)
  Jox = ZeroMat(Par$l, 1)

  if (Method == "Euler") {
    for (i1 in 1:(Par$l-1)) {
      B = matrix(data = c(1,-exp(Par$p),Derv(npoints = DerApprox, CoefMat = T)[1],Derv(npoints = DerApprox, CoefMat = T)[1]), byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(0, -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2]
      for (j1 in 2:(Par$j-1)) {
        Ox[i1+1,j1] = Ox[i1,j1] + Dm*(Ox[i1, j1-1] - 2*Ox[i1, j1] + Ox[i1, j1+1])
        Red[i1+1, j1] = Red[i1, j1] + Dm*(Red[i1, j1-1] + Red[i1, j1+1] - 2*Red[i1,j1])
      }
    }

    Jox = Derv(Ox = Ox, h = Par$h, npoints = DerApprox)

  } else if (Method == "RK4") {

    for (i1 in 1:(Par$l-1)) {
      k1 = ZeroMat(Par$j)
      k2 = ZeroMat(Par$j)
      k3 = ZeroMat(Par$j)
      k4 = ZeroMat(Par$j)
      k1red = ZeroMat(Par$j)
      k2red = ZeroMat(Par$j)
      k3red = ZeroMat(Par$j)
      k4red = ZeroMat(Par$j)
      B = matrix(data = c(1,-exp(Par$p),Derv(npoints = DerApprox, CoefMat = T)[1], Derv(npoints = DerApprox, CoefMat = T)[1]), byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(0, -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2]

      for (j1 in 2:(Par$j-1)) {
        k1[j1] = Dm*(Ox[i1, j1 -1] - 2*Ox[i1, j1] + Ox[i1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + k1[j1]*0.5
        k1red[j1] = Dm*(Red[i1, j1 -1] - 2*Red[i1, j1] + Red[i1, j1+1])
        Red[i1 + 1,j1] = Red[i1,j1] + k1red[j1]*0.5
      }

      B = matrix(data = c(1,-exp(Par$p),Derv(npoints = DerApprox, CoefMat = T)[1], Derv(npoints = DerApprox, CoefMat = T)[1]), byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(0, -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1+1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]

      for (j1 in 2:(Par$j-1)) {
        k2[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + k2[j1]*0.5
        k2red[j1] = Dm*(Red[i1 + 1, j1 -1] - 2*Red[i1 + 1, j1] + Red[i1 + 1, j1+1])
        Red[i1 + 1,j1] = Red[i1,j1] + k2red[j1]*0.5
      }

      B = matrix(data = c(1,-exp(Par$p),Derv(npoints = DerApprox, CoefMat = T)[1], Derv(npoints = DerApprox, CoefMat = T)[1]), byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(0, -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1+1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]

      for (j1 in 2:(Par$j-1)) {
        k3[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + k3[j1]
        k3red[j1] = Dm*(Red[i1 + 1, j1 -1] - 2*Red[i1 + 1, j1] + Red[i1 + 1, j1+1])
        Red[i1 + 1,j1] = Red[i1,j1] + k3red[j1]
      }

      B = matrix(data = c(1,-exp(Par$p),Derv(npoints = DerApprox, CoefMat = T)[1], Derv(npoints = DerApprox, CoefMat = T)[1]), byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(0, -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1+1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]

      for (j1 in 2:(Par$j-1)) {
        k4[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + (k1[j1] + 2*k2[j1] + 2*k3[j1] + k4[j1])/6
        k4red[j1] = Dm*(Red[i1 + 1, j1 -1] - 2*Red[i1 + 1, j1] + Red[i1 + 1, j1+1])
        Red[i1 + 1,j1] = Red[i1,j1] + (k1red[j1] + 2*k2red[j1] + 2*k3red[j1] + k4red[j1])/6
      }
    }

    Jox = Derv(Ox = Ox, h = Par$h, npoints = DerApprox)

  }  else if (Method == "BI") {
    al1 = 1/(Par$h^2)
    al2 = -2/(Par$h^2)
    al3 = 1/(Par$h^2)
    a1 = (al2 - 1/Par$dtn)/al1
    a2 = al3/al1

    for (i1 in 1:(Par$l-1)) {

      B = matrix(data = c(1,-exp(Par$p),Derv(npoints = DerApprox, CoefMat = T)[1], Derv(npoints = DerApprox, CoefMat = T)[1]), byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(0, -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2]


      bOx = (-Ox[i1,(2:(Par$j-1))]/(al1*Par$dtn))
      bRed = (-Red[i1,(2:(Par$j-1))]/(al1*Par$dtn))
      A = c(rep(1,Par$j-2))
      A1 = c(rep(a1,Par$j-2))
      A2 = c(rep(a2,Par$j-2))
      bOx[Par$j-2] = bOx[Par$j-2] - A2[Par$j-2]*1
      bRed[Par$j-2] = bRed[Par$j-2] - A2[Par$j-2]*Red[i1,Par$j]
      uox = c(rep(0, DerApprox))
      vox = c(rep(1, DerApprox))
      uRed = c(rep(0, DerApprox))
      vRed = c(rep(1, DerApprox))

      for (j1 in ((Par$j-3):1)) {

        bOx[j1] = bOx[j1] - A2[Par$j-2]*bOx[j1+1]/A1[j1+1]
        bRed[j1] = bRed[j1] - A2[Par$j-2]*bRed[j1+1]/A1[j1+1]
        A1[j1] = A1[j1] - A2[Par$j-2]/A1[j1+1]

      }
      for (m in 2:DerApprox) {
        uox[m] = (bOx[m-1] - uox[m-1])/A1[m-1]
        vox[m] = -vox[m-1]/A1[m-1]
        uRed[m] = (bRed[m-1] - uRed[m-1])/A1[m-1]
        vRed[m] = -vRed[m-1]/A1[m-1]
      }

      B = matrix(data = c(1,-exp(Par$p), Derv(npoints = DerApprox, CoefMat = T)[1] + sum(vox[2:DerApprox]*Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]), Derv(npoints = DerApprox, CoefMat = T)[1] + sum(vRed[2:DerApprox]*Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox])), byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(0, -sum(uox[2:DerApprox]*Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]) - sum(uRed[2:DerApprox]*Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]


      for (j1 in 1:(Par$j-2)) {
        Ox[i1+1,j1+1] = (bOx[j1] -Ox[i1+1,j1])/A1[j1]
        Red[i1+1,j1+1] = (bRed[j1] -Red[i1+1,j1])/A1[j1]
      }
    }

    Jox = Derv(Ox = Ox, h = Par$h, npoints = DerApprox)

  } else if (Method == "CN") {

    al1 = 1/(Par$h^2)
    al2 = -2/(Par$h^2)
    al3 = 1/(Par$h^2)
    a1 = (al2 - 2/Par$dtn)/al1
    a2 = al3/al1
    a3 = (al2 + 2/Par$dtn)/al1

    for (i1 in 1:(Par$l-1)) {

      B = matrix(data = c(1,-exp(Par$p),Derv(npoints = DerApprox, CoefMat = T)[1], Derv(npoints = DerApprox, CoefMat = T)[1]), byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(0, -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2]


      bOx = -a3*Ox[i1,(2:(Par$j-1))] - Ox[i1,(1:(Par$j-2))] - a2*Ox[i1,(3:Par$j)]
      bRed = -a3*Red[i1,(2:(Par$j-1))]- Red[i1,(1:(Par$j-2))] - a2*Red[i1,(3:Par$j)]
      A = c(rep(1,Par$j-2))
      A1 = c(rep(a1,Par$j-2))
      A2 = c(rep(a2,Par$j-2))
      bOx[Par$j-2] = bOx[Par$j-2] - A2[Par$j-2]*1
      bRed[Par$j-2] = bRed[Par$j-2] - A2[Par$j-2]*Red[i1,Par$j]
      uox = c(rep(0, DerApprox))
      vox = c(rep(1, DerApprox))
      uRed = c(rep(0, DerApprox))
      vRed = c(rep(1, DerApprox))

      for (j1 in ((Par$j-3):1)) {

        bOx[j1] = bOx[j1] - A2[Par$j-2]*bOx[j1+1]/A1[j1+1]
        bRed[j1] = bRed[j1] - A2[Par$j-2]*bRed[j1+1]/A1[j1+1]
        A1[j1] = A1[j1] - A2[Par$j-2]/A1[j1+1]

      }
      for (m in 2:DerApprox) {
        uox[m] = (bOx[m-1] - uox[m-1])/A1[m-1]
        vox[m] = -vox[m-1]/A1[m-1]
        uRed[m] = (bRed[m-1] - uRed[m-1])/A1[m-1]
        vRed[m] = -vRed[m-1]/A1[m-1]
      }

      B = matrix(data = c(1,-exp(Par$p), Derv(npoints = DerApprox, CoefMat = T)[1] + sum(vox[2:DerApprox]*Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]), Derv(npoints = DerApprox, CoefMat = T)[1] + sum(vRed[2:DerApprox]*Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox])), byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(0, -sum(uox[2:DerApprox]*Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]) - sum(uRed[2:DerApprox]*Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]


      for (j1 in 1:(Par$j-2)) {
        Ox[i1+1,j1+1] = (bOx[j1] -Ox[i1+1,j1])/A1[j1]
        Red[i1+1,j1+1] = (bRed[j1] -Red[i1+1,j1])/A1[j1]
      }
    }

    Jox = Derv(Ox = Ox, h = Par$h, npoints = DerApprox)


  } else if (Method == "BDF") {

    al1 = 1/(Par$h^2)
    al2 = -2/(Par$h^2)
    al3 = 1/(Par$h^2)
    a1 = (al2 - 1.5/Par$dtn)/al1
    a2 = al3/al1


    for (i1 in 1:(Par$l-1)) {

      B = matrix(data = c(1,-exp(Par$p),Derv(npoints = DerApprox, CoefMat = T)[1], Derv(npoints = DerApprox, CoefMat = T)[1]), byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(0, -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2]

      if (i1 == 1) {
        bOx = -2*Ox[i1,2:(Par$j-1)]/(Par$dtn*al1) + Ox[1,2:(Par$j-1)]/(2*Par$dtn*al1)
        bRed = -2*Red[i1,2:(Par$j-1)]/(Par$dtn*al1) + Red[1,2:(Par$j-1)]/(2*Par$dtn*al1)
      } else {
        bOx = -2*Ox[i1,2:(Par$j-1)]/(Par$dtn*al1) + Ox[i1-1,2:(Par$j-1)]/(2*Par$dtn*al1)
        bRed = -2*Red[i1,2:(Par$j-1)]/(Par$dtn*al1) + Red[i1-1,2:(Par$j-1)]/(2*Par$dtn*al1)
      }
      A = c(rep(1,Par$j-2))
      A1 = c(rep(a1,Par$j-2))
      A2 = c(rep(a2,Par$j-2))
      bOx[Par$j-2] = bOx[Par$j-2] - A2[Par$j-2]*1
      bRed[Par$j-2] = bRed[Par$j-2] - A2[Par$j-2]*Red[i1,Par$j]
      uox = c(rep(0, DerApprox))
      vox = c(rep(1, DerApprox))
      uRed = c(rep(0, DerApprox))
      vRed = c(rep(1, DerApprox))

      for (j1 in ((Par$j-3):1)) {

        bOx[j1] = bOx[j1] - A2[Par$j-2]*bOx[j1+1]/A1[j1+1]
        bRed[j1] = bRed[j1] - A2[Par$j-2]*bRed[j1+1]/A1[j1+1]
        A1[j1] = A1[j1] - A2[Par$j-2]/A1[j1+1]

      }
      for (m in 2:DerApprox) {
        uox[m] = (bOx[m-1] - uox[m-1])/A1[m-1]
        vox[m] = -vox[m-1]/A1[m-1]
        uRed[m] = (bRed[m-1] - uRed[m-1])/A1[m-1]
        vRed[m] = -vRed[m-1]/A1[m-1]
      }
      B = matrix(data = c(1,-exp(Par$p), Derv(npoints = DerApprox, CoefMat = T)[1] + sum(vox[2:DerApprox]*Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]), Derv(npoints = DerApprox, CoefMat = T)[1] + sum(vRed[2:DerApprox]*Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox])), byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(0, -sum(uox[2:DerApprox]*Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]) - sum(uRed[2:DerApprox]*Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]


      for (j1 in 1:(Par$j-2)) {
        Ox[i1+1,j1+1] = (bOx[j1] -Ox[i1+1,j1])/A1[j1]
        Red[i1+1,j1+1] = (bRed[j1] -Red[i1+1,j1])/A1[j1]
      }
    }

    Jox = Derv(Ox = Ox, h = Par$h, npoints = DerApprox)


  } else if (!(Method %in% c("Euler", "BI", "RK4", "CN", "BDF"))) {
    return("Available methods are Euler, BI, RK4, CN and BDF")
  }

  G = Jox
  i = (n*Par$FA*G*Dx*Area*Co)/(sqrt(Dx*Par$tau))

  graphy = ggplot(data = data.frame(i[1:(length(i)-1)],Par$t[1:(length(i)-1)]),
                  aes(y = i[1:(length(i)-1)], x = Par$t[1:(length(i)-1)])) +
    geom_point() + xlab("t / s") +
    ylab("I / A") + theme_classic()

  if (errCheck == TRUE){
    return(list(G,Dx,Co,Par$dtn,Par$h,Par$l,Par$j,i,n,Area))
  } else {
    return(graphy)
  }
}


#' Cottrel current check for the Chronoamperometric simulation
#'
#' Return a graph G/Gcot vs t of the electrochemical process
#'
#' @param Elefun the function to be checked = ChronAmp, PotStep
#'
#'
#' @return A graph G/Gcot vs t for the simulation data selected
#'
#' @examples
#' CottrCheck(ChronAmp(errCheck = TRUE, Method = "BI"))
#'
#' @export
#' @import ggplot2


CottrCheck = function(Elefun) {

  FA = 96485
  R = 8.3145
  Check = Elefun
  if (length(Check) == 9){
    return("ErrCheck inside the called function should be activated")
  } else {
    vt = c(1:Check[[6]])

    if (length(Check) == 10){

      Gcot = 1/sqrt(3.14*Check[[4]]*vt)

    } else if (length(Check) == 12){
      Gcot = (1/sqrt(3.14*Check[[4]]*vt))/(1+ (1/Check[[12]])*exp(Check[[11]]))
    }

    Err = (Check[[1]]/Gcot)
    t = Check[[4]]*vt
    ErrorGraphy = ggplot(data = data.frame(Err[1:(length(Err)-1)],t[1:(length(Err)-1)]), aes(y = Err[1:(length(Err)-1)], x = t[1:(length(Err)-1)])) + geom_point() +xlab("Time(s)") +
      ylab("G/Gcott") + theme_classic()
    return(ErrorGraphy)
  }
}


#' Linear Sweep digitial simulation
#'
#' Return a graph I vs E of the electrochemical process
#'
#' @param Co bulk concentration expressed in Molar
#' @param Dx diffusion coefficient expressed in cm^2/s
#' @param Eo reduction potential of the species expressed in Volt
#' @param Dm simulation parameter, maximum 0.5 for explicit methods
#' @param Vi initial potential of the sweep expressed in Volt
#' @param Vf final potential of the sweep expressed in Volt
#' @param Vs  potential scan rate of the simulation expressed in V/s
#' @param ko heterogeneous electron transfer rate constant expressed in m/s
#' @param alpha charge transfer coefficient
#' @param Temp temperature in kelvin
#' @param n number of electrons involved in the process
#' @param Area area of the electrode expressed in cm^2
#' @param l number of time steps of the simulation
#' @param DerApprox number of point for the approximation of the first derivative
#' @param errCheck if true the function returns a list with parameters for CottrCheck function
#' @param Method method to be used for the simulation = "Euler" "BI" "RK4" "CN" "BDF"
#'
#'
#' @return if errCheck == F a graph I vs E, if errCheck == T a list
#'
#' @examples
#' LinSwp(Co = 0.001, Dm =0.45, DerApprox = 2, errCheck = FALSE, Method = "Euler")
#'
#' @export
#' @import ggplot2
#'


LinSwp = function(Co = 0.001, Dx = 0.00001, Eo = 0, Dm = 0.45,
                  Vi = 0.3, Vf = -0.3, Vs = 0.001, ko = 0.01,
                  alpha = 0.5, Temp = 298.15, n = 1, Area = 1, l = 100,
                  DerApprox = 2, errCheck = FALSE, Method = "Euler"){


  Par = ParCall("LinSwp", n. = n, Temp. = Temp, Dx1. = Dx,
                Eo1. = Eo, Dm. = Dm, Vi. = Vi,
                Vf. = Vf, Vs. = Vs, ko1. = ko, alpha1. = alpha, l. = l)
  Ox = OneMat(Par$l+1, Par$j)
  Red =ZeroMat(Par$l+1, Par$j)
  Jox = ZeroMat(Par$l+1, 1)

  if (Method == "Euler") {
    for (i1 in 1:Par$l) {
      B = matrix(data = c((Par$Kf[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Par$Kb[i1]*Par$h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]), -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2]
      for (j1 in 2:(Par$j-1)) {
        Ox[i1+1,j1] = Ox[i1,j1] + Dm*(Ox[i1, j1-1] - 2*Ox[i1, j1] + Ox[i1, j1+1])
        Red[i1+1, j1] = Red[i1, j1] + Dm*(Red[i1, j1-1] + Red[i1, j1+1] - 2*Red[i1,j1])
      }
    }

    Jox = Derv(Ox = Ox, h = Par$h, npoints = DerApprox)

  } else if (Method == "RK4") {

    for (i1 in 1:(Par$l-1)) {
      k1 = ZeroMat(Par$j)
      k2 = ZeroMat(Par$j)
      k3 = ZeroMat(Par$j)
      k4 = ZeroMat(Par$j)
      k1red = ZeroMat(Par$j)
      k2red = ZeroMat(Par$j)
      k3red = ZeroMat(Par$j)
      k4red = ZeroMat(Par$j)
      B = matrix(data = c((Par$Kf[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Par$Kb[i1]*Par$h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2]

      for (j1 in 2:(Par$j-1)) {
        k1[j1] = Dm*(Ox[i1, j1 -1] - 2*Ox[i1, j1] + Ox[i1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + k1[j1]*0.5
        k1red[j1] = Dm*(Red[i1, j1 -1] - 2*Red[i1, j1] + Red[i1, j1+1])
        Red[i1 + 1,j1] = Red[i1,j1] + k1red[j1]*0.5
      }

      B = matrix(data = c((Par$Kf[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Par$Kb[i1]*Par$h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1+1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]

      for (j1 in 2:(Par$j-1)) {
        k2[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + k2[j1]*0.5
        k2red[j1] = Dm*(Red[i1 + 1, j1 -1] - 2*Red[i1 + 1, j1] + Red[i1 + 1, j1+1])
        Red[i1 + 1,j1] = Red[i1,j1] + k2red[j1]*0.5
      }

      B = matrix(data = c((Par$Kf[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Par$Kb[i1]*Par$h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1+1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]

      for (j1 in 2:(Par$j-1)) {
        k3[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + k3[j1]
        k3red[j1] = Dm*(Red[i1 + 1, j1 -1] - 2*Red[i1 + 1, j1] + Red[i1 + 1, j1+1])
        Red[i1 + 1,j1] = Red[i1,j1] + k3red[j1]
      }

      B = matrix(data = c((Par$Kf[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Par$Kb[i1]*Par$h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1+1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]

      for (j1 in 2:(Par$j-1)) {
        k4[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + (k1[j1] + 2*k2[j1] + 2*k3[j1] + k4[j1])/6
        k4red[j1] = Dm*(Red[i1 + 1, j1 -1] - 2*Red[i1 + 1, j1] + Red[i1 + 1, j1+1])
        Red[i1 + 1,j1] = Red[i1,j1] + (k1red[j1] + 2*k2red[j1] + 2*k3red[j1] + k4red[j1])/6
      }
    }

    Jox = Derv(Ox = Ox, h = Par$h, npoints = DerApprox)

  } else if (Method == "BI") {
    al1 = 1/(Par$h^2)
    al2 = -2/(Par$h^2)
    al3 = 1/(Par$h^2)
    a1 = (al2 - 1/Par$dtn)/al1
    a2 = al3/al1

    for (i1 in 1:(Par$l-1)) {

      B = matrix(data = c((Par$Kf[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Par$Kb[i1]*Par$h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2]

      bOx = (-Ox[i1,(2:(Par$j-1))]/(al1*Par$dtn))
      bRed = (-Red[i1,(2:(Par$j-1))]/(al1*Par$dtn))
      A = c(rep(1,Par$j-2))
      A1 = c(rep(a1,Par$j-2))
      A2 = c(rep(a2,Par$j-2))
      bOx[Par$j-2] = bOx[Par$j-2] - A2[Par$j-2]*1
      bRed[Par$j-2] = bRed[Par$j-2] - A2[Par$j-2]*Red[i1,Par$j]
      uox = c(rep(0, DerApprox))
      vox = c(rep(1, DerApprox))
      uRed = c(rep(0, DerApprox))
      vRed = c(rep(1, DerApprox))

      for (j1 in ((Par$j-3):1)) {

        bOx[j1] = bOx[j1] - A2[Par$j-2]*bOx[j1+1]/A1[j1+1]
        bRed[j1] = bRed[j1] - A2[Par$j-2]*bRed[j1+1]/A1[j1+1]
        A1[j1] = A1[j1] - A2[Par$j-2]/A1[j1+1]

      }

      for (m in 2:DerApprox) {
        uox[m] = (bOx[m-1] - uox[m-1])/A1[m-1]
        vox[m] = -vox[m-1]/A1[m-1]
        uRed[m] = (bRed[m-1] - uRed[m-1])/A1[m-1]
        vRed[m] = -vRed[m-1]/A1[m-1]
      }

      B = matrix(data = c((Par$Kf[i1]*Par$h -  sum(vox*Derv(npoints = DerApprox, CoefMat = T))),
                          -Par$Kb[i1]*Par$h,
                          sum(vox*Derv(npoints = DerApprox, CoefMat = T)),
                          sum(vRed*Derv(npoints = DerApprox, CoefMat = T))),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(uox*Derv(npoints = DerApprox, CoefMat = T)),
                          -sum(uox*Derv(npoints = DerApprox, CoefMat = T)) - sum(uRed*Derv(npoints = DerApprox, CoefMat = T))))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]


      for (j1 in 1:(Par$j-2)) {
        Ox[i1+1,j1+1] = (bOx[j1] -Ox[i1+1,j1])/A1[j1]
        Red[i1+1,j1+1] = (bRed[j1] -Red[i1+1,j1])/A1[j1]
      }
    }

    Jox = Derv(Ox = Ox, h = Par$h, npoints = DerApprox)

  } else if (Method == "CN") {

    al1 = 1/(Par$h^2)
    al2 = -2/(Par$h^2)
    al3 = 1/(Par$h^2)
    a1 = (al2 - 2/Par$dtn)/al1
    a2 = al3/al1
    a3 = (al2 + 2/Par$dtn)/al1

    for (i1 in 1:(Par$l-1)) {

      B = matrix(data = c((Par$Kf[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Par$Kb[i1]*Par$h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2]

      bOx = -a3*Ox[i1,(2:(Par$j-1))] - Ox[i1,(1:(Par$j-2))] - a2*Ox[i1,(3:Par$j)]
      bRed = -a3*Red[i1,(2:(Par$j-1))]- Red[i1,(1:(Par$j-2))] - a2*Red[i1,(3:Par$j)]
      A = c(rep(1,Par$j-2))
      A1 = c(rep(a1,Par$j-2))
      A2 = c(rep(a2,Par$j-2))
      bOx[Par$j-2] = bOx[Par$j-2] - A2[Par$j-2]*1
      bRed[Par$j-2] = bRed[Par$j-2] - A2[Par$j-2]*Red[i1,Par$j]
      uox = c(rep(0, DerApprox))
      vox = c(rep(1, DerApprox))
      uRed = c(rep(0, DerApprox))
      vRed = c(rep(1, DerApprox))

      for (j1 in ((Par$j-3):1)) {

        bOx[j1] = bOx[j1] - A2[Par$j-2]*bOx[j1+1]/A1[j1+1]
        bRed[j1] = bRed[j1] - A2[Par$j-2]*bRed[j1+1]/A1[j1+1]
        A1[j1] = A1[j1] - A2[Par$j-2]/A1[j1+1]

      }

      for (m in 2:DerApprox) {
        uox[m] = (bOx[m-1] - uox[m-1])/A1[m-1]
        vox[m] = -vox[m-1]/A1[m-1]
        uRed[m] = (bRed[m-1] - uRed[m-1])/A1[m-1]
        vRed[m] = -vRed[m-1]/A1[m-1]
      }

      B = matrix(data = c((Par$Kf[i1]*Par$h -  sum(vox*Derv(npoints = DerApprox, CoefMat = T))),
                          -Par$Kb[i1]*Par$h,
                          sum(vox*Derv(npoints = DerApprox, CoefMat = T)),
                          sum(vRed*Derv(npoints = DerApprox, CoefMat = T))),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(uox*Derv(npoints = DerApprox, CoefMat = T)),
                          -sum(uox*Derv(npoints = DerApprox, CoefMat = T)) - sum(uRed*Derv(npoints = DerApprox, CoefMat = T))))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]


      for (j1 in 1:(Par$j-2)) {
        Ox[i1+1,j1+1] = (bOx[j1] -Ox[i1+1,j1])/A1[j1]
        Red[i1+1,j1+1] = (bRed[j1] -Red[i1+1,j1])/A1[j1]
      }
    }

    Jox = Derv(Ox = Ox, h = Par$h, npoints = DerApprox)


  } else if (Method == "BDF") {

    al1 = 1/(Par$h^2)
    al2 = -2/(Par$h^2)
    al3 = 1/(Par$h^2)
    a1 = (al2 - 1.5/Par$dtn)/al1
    a2 = al3/al1


    for (i1 in 1:(Par$l-1)) {

      B = matrix(data = c((Par$Kf[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Par$Kb[i1]*Par$h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2]

      if (i1 == 1) {
        bOx = -2*Ox[i1,2:(Par$j-1)]/(Par$dtn*al1) + Ox[1,2:(Par$j-1)]/(2*Par$dtn*al1)
        bRed = -2*Red[i1,2:(Par$j-1)]/(Par$dtn*al1) + Red[1,2:(Par$j-1)]/(2*Par$dtn*al1)
      } else {
        bOx = -2*Ox[i1,2:(Par$j-1)]/(Par$dtn*al1) + Ox[i1-1,2:(Par$j-1)]/(2*Par$dtn*al1)
        bRed = -2*Red[i1,2:(Par$j-1)]/(Par$dtn*al1) + Red[i1-1,2:(Par$j-1)]/(2*Par$dtn*al1)
      }
      A = c(rep(1,Par$j-2))
      A1 = c(rep(a1,Par$j-2))
      A2 = c(rep(a2,Par$j-2))
      bOx[Par$j-2] = bOx[Par$j-2] - A2[Par$j-2]*1
      bRed[Par$j-2] = bRed[Par$j-2] - A2[Par$j-2]*Red[i1,Par$j]
      uox = c(rep(0, DerApprox))
      vox = c(rep(1, DerApprox))
      uRed = c(rep(0, DerApprox))
      vRed = c(rep(1, DerApprox))

      for (j1 in ((Par$j-3):1)) {

        bOx[j1] = bOx[j1] - A2[Par$j-2]*bOx[j1+1]/A1[j1+1]
        bRed[j1] = bRed[j1] - A2[Par$j-2]*bRed[j1+1]/A1[j1+1]
        A1[j1] = A1[j1] - A2[Par$j-2]/A1[j1+1]

      }

      for (m in 2:DerApprox) {
        uox[m] = (bOx[m-1] - uox[m-1])/A1[m-1]
        vox[m] = -vox[m-1]/A1[m-1]
        uRed[m] = (bRed[m-1] - uRed[m-1])/A1[m-1]
        vRed[m] = -vRed[m-1]/A1[m-1]
      }

      B = matrix(data = c((Par$Kf[i1]*Par$h -  sum(vox*Derv(npoints = DerApprox, CoefMat = T))),
                          -Par$Kb[i1]*Par$h,
                          sum(vox*Derv(npoints = DerApprox, CoefMat = T)),
                          sum(vRed*Derv(npoints = DerApprox, CoefMat = T))),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(uox*Derv(npoints = DerApprox, CoefMat = T)),
                          -sum(uox*Derv(npoints = DerApprox, CoefMat = T)) - sum(uRed*Derv(npoints = DerApprox, CoefMat = T))))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]


      for (j1 in 1:(Par$j-2)) {
        Ox[i1+1,j1+1] = (bOx[j1] -Ox[i1+1,j1])/A1[j1]
        Red[i1+1,j1+1] = (bRed[j1] -Red[i1+1,j1])/A1[j1]
      }
    }

    Jox = Derv(Ox = Ox, h = Par$h, npoints = DerApprox)


  } else if (!(Method %in% c("Euler", "BI", "RK4", "CN", "BDF"))) {
    return("Available methods are Euler, BI, RK4, CN and BDF")
  }

  G = Jox
  i = (n*Par$FA*G*Dx*Area*Co)/(sqrt(Dx*Par$tau))

  graphy = ggplot(data = data.frame(i[1:(length(i)-1)], Par$PotentialScan[1:(length(i)-1)]), aes(y = i[1:(length(i)-1)], x = Par$PotentialScan[1:(length(i)-1)])) +
    geom_point() + scale_x_continuous(trans = "reverse") +
    xlab("E / V") +
    ylab("I / A") +
    theme_classic()

  if (errCheck == TRUE){
    return(list(G,Dx,Co,Par$dtn,Par$h,Par$l,Par$j,i,n,Area,Par$p,Par$Da))
  } else {
    return(graphy)
  }
}


#' Cyclic voltammetry digitial simulation
#'
#' Return a graph I vs E of the electrochemical process
#'
#' @param Co bulk concentration expressed in Molar
#' @param Dx diffusion coefficient expressed in cm^2/s
#' @param Eo reduction potential of the species expressed in Volts
#' @param Dm simulation parameter, maximum 0.5 for explicit methods
#' @param Vi initial potential of the sweep expressed in Volts
#' @param Vf final potential of the sweepexpressed in Volts
#' @param Vs  potential scan rate of the simulation expressed in V/s
#' @param ko heterogeneous electron transfer rate constant expressed in m/s
#' @param alpha charge transfer coefficient
#' @param Temp temperature in kelvin
#' @param n number of electrons involved in the process
#' @param Area area of the electrode expressed in cm^2
#' @param l number of time steps of the simulation
#' @param DerApprox number of point for the approximation of the first derivative
#' @param errCheck if true the function returns a list with parameters for CottrCheck function
#' @param Method method to be used for the simulation = "Euler" "BI" "RK4" "CN" "BDF"
#'
#'
#' @return if errCheck == F a graph I vs E, if errCheck == T a list
#'
#' @examples
#' CV(Co = 0.001, DerApprox = 2, Dm = 0.45, errCheck = FALSE, Method = "Euler")
#'
#' @export
#' @import ggplot2
#'

CV = function(Co = 0.001, Dx = 0.00001, Eo = 0, Dm = 0.45,
              Vi = 0.3, Vf = -0.3, Vs = 0.001, ko = 0.01,
              alpha = 0.5, Temp = 298.15, n = 1, Area = 1, l = 100,
              DerApprox = 2, errCheck = FALSE, Method = "Euler"){


  Par = ParCall("CV", n. = n, Temp. = Temp, Dx1. = Dx,
                Eo1. = Eo, Dm. = Dm, Vi. = Vi,
                Vf. = Vf, Vs. = Vs, ko1. = ko, alpha1. = alpha, l. = l)
  Ox = OneMat(Par$l +1, Par$j)
  Red =ZeroMat(Par$l +1, Par$j)
  Jox = ZeroMat(Par$l+1, 1)

  if (Method == "Euler") {
    for (i1 in 1:Par$l) {
      B = matrix(data = c((Par$Kf[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Par$Kb[i1]*Par$h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]), -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2]
      for (j1 in 2:(Par$j-1)) {
        Ox[i1+1,j1] = Ox[i1,j1] + Dm*(Ox[i1, j1-1] - 2*Ox[i1, j1] + Ox[i1, j1+1])
        Red[i1+1, j1] = Red[i1, j1] + Dm*(Red[i1, j1-1] + Red[i1, j1+1] - 2*Red[i1,j1])
      }
    }

    Jox = Derv(Ox = Ox, h = Par$h, npoints = DerApprox)

  } else if (Method == "RK4") {

    for (i1 in 1:(Par$l-1)) {
      k1 = ZeroMat(Par$j)
      k2 = ZeroMat(Par$j)
      k3 = ZeroMat(Par$j)
      k4 = ZeroMat(Par$j)
      k1red = ZeroMat(Par$j)
      k2red = ZeroMat(Par$j)
      k3red = ZeroMat(Par$j)
      k4red = ZeroMat(Par$j)
      B = matrix(data = c((Par$Kf[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Par$Kb[i1]*Par$h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2]

      for (j1 in 2:(Par$j-1)) {
        k1[j1] = Dm*(Ox[i1, j1 -1] - 2*Ox[i1, j1] + Ox[i1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + k1[j1]*0.5
        k1red[j1] = Dm*(Red[i1, j1 -1] - 2*Red[i1, j1] + Red[i1, j1+1])
        Red[i1 + 1,j1] = Red[i1,j1] + k1red[j1]*0.5
      }

      B = matrix(data = c((Par$Kf[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Par$Kb[i1]*Par$h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1+1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]

      for (j1 in 2:(Par$j-1)) {
        k2[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + k2[j1]*0.5
        k2red[j1] = Dm*(Red[i1 + 1, j1 -1] - 2*Red[i1 + 1, j1] + Red[i1 + 1, j1+1])
        Red[i1 + 1,j1] = Red[i1,j1] + k2red[j1]*0.5
      }

      B = matrix(data = c((Par$Kf[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Par$Kb[i1]*Par$h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1+1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]

      for (j1 in 2:(Par$j-1)) {
        k3[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + k3[j1]
        k3red[j1] = Dm*(Red[i1 + 1, j1 -1] - 2*Red[i1 + 1, j1] + Red[i1 + 1, j1+1])
        Red[i1 + 1,j1] = Red[i1,j1] + k3red[j1]
      }

      B = matrix(data = c((Par$Kf[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Par$Kb[i1]*Par$h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1+1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]

      for (j1 in 2:(Par$j-1)) {
        k4[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + (k1[j1] + 2*k2[j1] + 2*k3[j1] + k4[j1])/6
        k4red[j1] = Dm*(Red[i1 + 1, j1 -1] - 2*Red[i1 + 1, j1] + Red[i1 + 1, j1+1])
        Red[i1 + 1,j1] = Red[i1,j1] + (k1red[j1] + 2*k2red[j1] + 2*k3red[j1] + k4red[j1])/6
      }
    }

    Jox = Derv(Ox = Ox, h = Par$h, npoints = DerApprox)

  } else if (Method == "BI") {
    al1 = 1/(Par$h^2)
    al2 = -2/(Par$h^2)
    al3 = 1/(Par$h^2)
    a1 = (al2 - 1/Par$dtn)/al1
    a2 = al3/al1

    for (i1 in 1:(Par$l-1)) {

      B = matrix(data = c((Par$Kf[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Par$Kb[i1]*Par$h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2]

      bOx = (-Ox[i1,(2:(Par$j-1))]/(al1*Par$dtn))
      bRed = (-Red[i1,(2:(Par$j-1))]/(al1*Par$dtn))
      A = c(rep(1,Par$j-2))
      A1 = c(rep(a1,Par$j-2))
      A2 = c(rep(a2,Par$j-2))
      bOx[Par$j-2] = bOx[Par$j-2] - A2[Par$j-2]*1
      bRed[Par$j-2] = bRed[Par$j-2] - A2[Par$j-2]*Red[i1,Par$j]
      uox = c(rep(0, DerApprox))
      vox = c(rep(1, DerApprox))
      uRed = c(rep(0, DerApprox))
      vRed = c(rep(1, DerApprox))

      for (j1 in ((Par$j-3):1)) {

        bOx[j1] = bOx[j1] - A2[Par$j-2]*bOx[j1+1]/A1[j1+1]
        bRed[j1] = bRed[j1] - A2[Par$j-2]*bRed[j1+1]/A1[j1+1]
        A1[j1] = A1[j1] - A2[Par$j-2]/A1[j1+1]

      }

      for (m in 2:DerApprox) {
        uox[m] = (bOx[m-1] - uox[m-1])/A1[m-1]
        vox[m] = -vox[m-1]/A1[m-1]
        uRed[m] = (bRed[m-1] - uRed[m-1])/A1[m-1]
        vRed[m] = -vRed[m-1]/A1[m-1]
      }

      B = matrix(data = c((Par$Kf[i1]*Par$h - sum(vox*Derv(npoints = DerApprox, CoefMat = T))),
                          -Par$Kb[i1]*Par$h,
                          sum(vox*Derv(npoints = DerApprox, CoefMat = T)),
                          sum(vRed*Derv(npoints = DerApprox, CoefMat = T))),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(uox*Derv(npoints = DerApprox, CoefMat = T)),
                          -sum(uox*Derv(npoints = DerApprox, CoefMat = T)) - sum(uRed*Derv(npoints = DerApprox, CoefMat = T))))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]


      for (j1 in 1:(Par$j-2)) {
        Ox[i1+1,j1+1] = (bOx[j1] -Ox[i1+1,j1])/A1[j1]
        Red[i1+1,j1+1] = (bRed[j1] -Red[i1+1,j1])/A1[j1]
      }
    }

    Jox = Derv(Ox = Ox, h = Par$h, npoints = DerApprox)

  } else if (Method == "CN") {

    al1 = 1/(Par$h^2)
    al2 = -2/(Par$h^2)
    al3 = 1/(Par$h^2)
    a1 = (al2 - 2/Par$dtn)/al1
    a2 = al3/al1
    a3 = (al2 + 2/Par$dtn)/al1

    for (i1 in 1:(Par$l-1)) {

      B = matrix(data = c((Par$Kf[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Par$Kb[i1]*Par$h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2]

      bOx = -a3*Ox[i1,(2:(Par$j-1))] - Ox[i1,(1:(Par$j-2))] - a2*Ox[i1,(3:Par$j)]
      bRed = -a3*Red[i1,(2:(Par$j-1))]- Red[i1,(1:(Par$j-2))] - a2*Red[i1,(3:Par$j)]
      A = c(rep(1,Par$j-2))
      A1 = c(rep(a1,Par$j-2))
      A2 = c(rep(a2,Par$j-2))
      bOx[Par$j-2] = bOx[Par$j-2] - A2[Par$j-2]*1
      bRed[Par$j-2] = bRed[Par$j-2] - A2[Par$j-2]*Red[i1,Par$j]
      uox = c(rep(0, DerApprox))
      vox = c(rep(1, DerApprox))
      uRed = c(rep(0, DerApprox))
      vRed = c(rep(1, DerApprox))

      for (j1 in ((Par$j-3):1)) {

        bOx[j1] = bOx[j1] - A2[Par$j-2]*bOx[j1+1]/A1[j1+1]
        bRed[j1] = bRed[j1] - A2[Par$j-2]*bRed[j1+1]/A1[j1+1]
        A1[j1] = A1[j1] - A2[Par$j-2]/A1[j1+1]

      }

      for (m in 2:DerApprox) {
        uox[m] = (bOx[m-1] - uox[m-1])/A1[m-1]
        vox[m] = -vox[m-1]/A1[m-1]
        uRed[m] = (bRed[m-1] - uRed[m-1])/A1[m-1]
        vRed[m] = -vRed[m-1]/A1[m-1]
      }

      B = matrix(data = c((Par$Kf[i1]*Par$h -  sum(vox*Derv(npoints = DerApprox, CoefMat = T))),
                          -Par$Kb[i1]*Par$h,
                          sum(vox*Derv(npoints = DerApprox, CoefMat = T)),
                          sum(vRed*Derv(npoints = DerApprox, CoefMat = T))),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(uox*Derv(npoints = DerApprox, CoefMat = T)),
                          -sum(uox*Derv(npoints = DerApprox, CoefMat = T)) - sum(uRed*Derv(npoints = DerApprox, CoefMat = T))))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]


      for (j1 in 1:(Par$j-2)) {
        Ox[i1+1,j1+1] = (bOx[j1] -Ox[i1+1,j1])/A1[j1]
        Red[i1+1,j1+1] = (bRed[j1] -Red[i1+1,j1])/A1[j1]
      }
    }

    Jox = Derv(Ox = Ox, h = Par$h, npoints = DerApprox)


  } else if (Method == "BDF") {

    al1 = 1/(Par$h^2)
    al2 = -2/(Par$h^2)
    al3 = 1/(Par$h^2)
    a1 = (al2 - 1.5/Par$dtn)/al1
    a2 = al3/al1

    for (i1 in 1:(Par$l-1)) {

      B = matrix(data = c((Par$Kf[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Par$Kb[i1]*Par$h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2]

      if (i1 == 1) {
        bOx = -2*Ox[i1,2:(Par$j-1)]/(Par$dtn*al1) + Ox[1,2:(Par$j-1)]/(2*Par$dtn*al1)
        bRed = -2*Red[i1,2:(Par$j-1)]/(Par$dtn*al1) + Red[1,2:(Par$j-1)]/(2*Par$dtn*al1)
      } else {
        bOx = -2*Ox[i1,2:(Par$j-1)]/(Par$dtn*al1) + Ox[i1-1,2:(Par$j-1)]/(2*Par$dtn*al1)
        bRed = -2*Red[i1,2:(Par$j-1)]/(Par$dtn*al1) + Red[i1-1,2:(Par$j-1)]/(2*Par$dtn*al1)
      }
      A = c(rep(1,Par$j-2))
      A1 = c(rep(a1,Par$j-2))
      A2 = c(rep(a2,Par$j-2))
      bOx[Par$j-2] = bOx[Par$j-2] - A2[Par$j-2]*1
      bRed[Par$j-2] = bRed[Par$j-2] - A2[Par$j-2]*Red[i1,Par$j]
      uox = c(rep(0, DerApprox))
      vox = c(rep(1, DerApprox))
      uRed = c(rep(0, DerApprox))
      vRed = c(rep(1, DerApprox))

      for (j1 in ((Par$j-3):1)) {

        bOx[j1] = bOx[j1] - A2[Par$j-2]*bOx[j1+1]/A1[j1+1]
        bRed[j1] = bRed[j1] - A2[Par$j-2]*bRed[j1+1]/A1[j1+1]
        A1[j1] = A1[j1] - A2[Par$j-2]/A1[j1+1]

      }

      for (m in 2:DerApprox) {
        uox[m] = (bOx[m-1] - uox[m-1])/A1[m-1]
        vox[m] = -vox[m-1]/A1[m-1]
        uRed[m] = (bRed[m-1] - uRed[m-1])/A1[m-1]
        vRed[m] = -vRed[m-1]/A1[m-1]
      }

      B = matrix(data = c((Par$Kf[i1]*Par$h -  sum(vox*Derv(npoints = DerApprox, CoefMat = T))),
                          -Par$Kb[i1]*Par$h,
                          sum(vox*Derv(npoints = DerApprox, CoefMat = T)),
                          sum(vRed*Derv(npoints = DerApprox, CoefMat = T))),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(uox*Derv(npoints = DerApprox, CoefMat = T)),
                          -sum(uox*Derv(npoints = DerApprox, CoefMat = T)) - sum(uRed*Derv(npoints = DerApprox, CoefMat = T))))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]


      for (j1 in 1:(Par$j-2)) {
        Ox[i1+1,j1+1] = (bOx[j1] -Ox[i1+1,j1])/A1[j1]
        Red[i1+1,j1+1] = (bRed[j1] -Red[i1+1,j1])/A1[j1]
      }
    }

    Jox = Derv(Ox = Ox, h = Par$h, npoints = DerApprox)


  } else if (!(Method %in% c("Euler", "BI", "RK4", "CN", "BDF"))) {
    return("Available methods are Euler, BI, RK4, CN and BDF")
  }

  G = Jox
  i = (n*Par$FA*G*Dx*Area*Co)/(sqrt(Dx*Par$tau))

  graphy = ggplot(data = data.frame(i[1:(length(i)-1)], Par$PotentialScan[1:(length(i)-1)]),
                  aes(y = i[1:(length(i)-1)], x = Par$PotentialScan[1:(length(i)-1)])) +
    geom_point() + scale_x_continuous(trans = "reverse") +
    xlab("E / V") +
    ylab("I / A") +
    theme_classic()

  if (errCheck == TRUE){
    return(list(G,Dx,Co,Par$dtn,Par$h,Par$l,Par$j,i,n,Area,Par$p,Par$Da))
  } else {
    return(graphy)
  }
}

#' EC behaviour cyclic voltammetry simulator
#'
#' Return a graph I vs E of the electrochemical process
#'
#' @param Co bulk concentration expressed in Molar
#' @param Dx diffusion coefficient expressed in cm^2/s
#' @param Eo reduction potential of the species expressed in Volt
#' @param Dm simulation parameter, maximum 0.5 for explicit methods
#' @param Vi initial potential of the sweep expressed in Volt
#' @param Vf final potential of the sweep expressed in Volt
#' @param Vs  potential scan rate of the simulation expressed in V/s
#' @param ko heterogeneous electron transfer rate constant expressed in m/s
#' @param kc rate constant of the reaction Red -> C expressed in s^-1
#' @param alpha charge transfer coefficient
#' @param Temp temperature in kelvin
#' @param n number of electrons involved in the process
#' @param l number of time steps of the simulation
#' @param Area area of the electrode expressed in cm^2
#' @param DerApprox number of point for the approximation of the first derivative
#' @param errCheck if true the function returns a list with parameters for CottrCheck function
#' @param Method method to be used for the simulation = "Euler" "BI" "RK4" "CN "BDF"
#'
#'
#' @return if errCheck == F a graph I vs E, if errCheck == T a list
#'
#' @examples
#' CVEC(Co = 0.001, DerApprox = 2, Dm = 0.45, kc = 0.00001, errCheck = FALSE, Method = "Euler")
#'
#' @export
#' @import ggplot2


CVEC = function(Co = 0.001, Dx = 0.00001, Eo = 0, Dm = 0.45,
                Vi = 0.3, Vf = -0.3, Vs = 0.001, ko = 0.01,
                kc = 0.001, l = 100,
                alpha = 0.5, Temp = 298.15, n = 1, Area = 1,
                DerApprox = 2, errCheck = FALSE, Method = "Euler"){

  Par = ParCall("CVEC", n. = n, Temp. = Temp, Dx1. = Dx,
                Eo1. = Eo, Dm. = Dm, Vi. = Vi, kc. = kc,
                Vf. = Vf, Vs. = Vs, ko1. = ko, alpha1. = alpha, l. = l)
  Ox = OneMat(Par$l +1, Par$j)
  Red =ZeroMat(Par$l +1, Par$j)
  Jox = ZeroMat(Par$l+1, 1)

  if (Method == "Euler") {
    for (i1 in 1:Par$l) {
      B = matrix(data = c((Par$Kf[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Par$Kb[i1]*Par$h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]), -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2] - Par$KC*Red[i1,1]
      for (j1 in 2:(Par$j-1)) {
        Ox[i1+1,j1] = Ox[i1,j1] + Dm*(Ox[i1, j1-1] - 2*Ox[i1, j1] + Ox[i1, j1+1])
        Red[i1+1, j1] = Red[i1, j1] + Dm*(Red[i1, j1-1] + Red[i1, j1+1] - 2*Red[i1,j1]) - Par$KC*Par$dtn*Red[i1,j1]
      }
    }

    Jox = Derv(Ox = Ox, h = Par$h, npoints = DerApprox)

  } else if (Method == "RK4") {

    for (i1 in 1:(Par$l-1)) {
      k1 = ZeroMat(Par$j)
      k2 = ZeroMat(Par$j)
      k3 = ZeroMat(Par$j)
      k4 = ZeroMat(Par$j)
      k1red = ZeroMat(Par$j)
      k2red = ZeroMat(Par$j)
      k3red = ZeroMat(Par$j)
      k4red = ZeroMat(Par$j)
      B = matrix(data = c((Par$Kf[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Par$Kb[i1]*Par$h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2]

      for (j1 in 2:(Par$j-1)) {
        k1[j1] = Dm*(Ox[i1, j1 -1] - 2*Ox[i1, j1] + Ox[i1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + k1[j1]*0.5
        k1red[j1] = Dm*(Red[i1, j1 -1] - 2*Red[i1, j1] + Red[i1, j1+1])
        Red[i1 + 1,j1] = Red[i1,j1] + k1red[j1]*0.5
      }

      B = matrix(data = c((Par$Kf[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Par$Kb[i1]*Par$h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1+1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]

      for (j1 in 2:(Par$j-1)) {
        k2[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + k2[j1]*0.5
        k2red[j1] = Dm*(Red[i1 + 1, j1 -1] - 2*Red[i1 + 1, j1] + Red[i1 + 1, j1+1])
        Red[i1 + 1,j1] = Red[i1,j1] + k2red[j1]*0.5
      }

      B = matrix(data = c((Par$Kf[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Par$Kb[i1]*Par$h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1+1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]

      for (j1 in 2:(Par$j-1)) {
        k3[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + k3[j1]
        k3red[j1] = Dm*(Red[i1 + 1, j1 -1] - 2*Red[i1 + 1, j1] + Red[i1 + 1, j1+1])
        Red[i1 + 1,j1] = Red[i1,j1] + k3red[j1]
      }

      B = matrix(data = c((Par$Kf[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Par$Kb[i1]*Par$h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1+1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2] - Par$KC*Red[i1+1,1]

      for (j1 in 2:(Par$j-1)) {
        k4[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + (k1[j1] + 2*k2[j1] + 2*k3[j1] + k4[j1])/6
        k4red[j1] = Dm*(Red[i1 + 1, j1 -1] - 2*Red[i1 + 1, j1] + Red[i1 + 1, j1+1])
        Red[i1 + 1,j1] = Red[i1,j1] + (k1red[j1] + 2*k2red[j1] + 2*k3red[j1] + k4red[j1])/6 - Par$KC*Par$dtn*Red[i1+1,j1]
      }
    }

    Jox = Derv(Ox = Ox, h = Par$h, npoints = DerApprox)

  } else if (Method == "BI") {
    al1 = 1/(Par$h^2)
    al2 = -2/(Par$h^2)
    al3 = 1/(Par$h^2)
    a1 = (al2 - 1/Par$dtn)/al1
    a2 = al3/al1

    for (i1 in 1:(Par$l-1)) {

      B = matrix(data = c((Par$Kf[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Par$Kb[i1]*Par$h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2] - Par$KC*Red[i1,1]

      bOx = (-Ox[i1,(2:(Par$j-1))]/(al1*Par$dtn))
      bRed = (-Red[i1,(2:(Par$j-1))]/(al1*Par$dtn)) + Par$KC*Red[i1,2:(Par$j-1)]/al1
      A = c(rep(1,Par$j-2))
      A1 = c(rep(a1,Par$j-2))
      A2 = c(rep(a2,Par$j-2))
      bOx[Par$j-2] = bOx[Par$j-2] - A2[Par$j-2]*1
      bRed[Par$j-2] = bRed[Par$j-2] - A2[Par$j-2]*(Red[i1,Par$j] - Par$KC*Red[i1,Par$j])
      uox = c(rep(0, DerApprox))
      vox = c(rep(1, DerApprox))
      uRed = c(rep(0, DerApprox))
      vRed = c(rep(1, DerApprox))

      for (j1 in ((Par$j-3):1)) {

        bOx[j1] = bOx[j1] - A2[Par$j-2]*bOx[j1+1]/A1[j1+1]
        bRed[j1] = bRed[j1] - A2[Par$j-2]*bRed[j1+1]/A1[j1+1]
        A1[j1] = A1[j1] - A2[Par$j-2]/A1[j1+1]

      }

      for (m in 2:DerApprox) {
        uox[m] = (bOx[m-1] - uox[m-1])/A1[m-1]
        vox[m] = -vox[m-1]/A1[m-1]
        uRed[m] = (bRed[m-1] - uRed[m-1])/A1[m-1]
        vRed[m] = -vRed[m-1]/A1[m-1]
      }

      B = matrix(data = c((Par$Kf[i1]*Par$h - sum(vox*Derv(npoints = DerApprox, CoefMat = T))),
                          -Par$Kb[i1]*Par$h,
                          sum(vox*Derv(npoints = DerApprox, CoefMat = T)),
                          sum(vRed*Derv(npoints = DerApprox, CoefMat = T))),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(uox*Derv(npoints = DerApprox, CoefMat = T)),
                          -sum(uox*Derv(npoints = DerApprox, CoefMat = T)) - sum(uRed*Derv(npoints = DerApprox, CoefMat = T))))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]


      for (j1 in 1:(Par$j-2)) {
        Ox[i1+1,j1+1] = (bOx[j1] - Ox[i1+1,j1])/A1[j1]
        Red[i1+1,j1+1] = (bRed[j1] - Red[i1+1,j1])/A1[j1]
      }
    }

    Jox = Derv(Ox = Ox, h = Par$h, npoints = DerApprox)

  } else if (Method == "CN") {

    al1 = 1/(Par$h^2)
    al2 = -2/(Par$h^2)
    al3 = 1/(Par$h^2)
    a1 = (al2 - 2/Par$dtn)/al1
    a2 = al3/al1
    a3 = (al2 + 2/Par$dtn)/al1

    for (i1 in 1:(Par$l-1)) {

      B = matrix(data = c((Par$Kf[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Par$Kb[i1]*Par$h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2] - Par$KC*Red[i1,1]

      bOx = -a3*Ox[i1,(2:(Par$j-1))] - Ox[i1,(1:(Par$j-2))] - a2*Ox[i1,(3:Par$j)]
      bRed = -a3*Red[i1,(2:(Par$j-1))]- Red[i1,(1:(Par$j-2))] - a2*Red[i1,(3:Par$j)] + Par$KC*Red[i1,2:(Par$j-1)]/al1
      A = c(rep(1,Par$j-2))
      A1 = c(rep(a1,Par$j-2))
      A2 = c(rep(a2,Par$j-2))
      bOx[Par$j-2] = bOx[Par$j-2] - A2[Par$j-2]*Ox[i1,Par$j]
      bRed[Par$j-2] = bRed[Par$j-2] - A2[Par$j-2]*(Red[i1,Par$j] - Par$KC*Red[i1,Par$j])
      uox = c(rep(0, DerApprox))
      vox = c(rep(1, DerApprox))
      uRed = c(rep(0, DerApprox))
      vRed = c(rep(1, DerApprox))

      for (j1 in ((Par$j-3):1)) {

        bOx[j1] = bOx[j1] - A2[Par$j-2]*bOx[j1+1]/A1[j1+1]
        bRed[j1] = bRed[j1] - A2[Par$j-2]*bRed[j1+1]/A1[j1+1]
        A1[j1] = A1[j1] - A2[Par$j-2]/A1[j1+1]

      }

      for (m in 2:DerApprox) {
        uox[m] = (bOx[m-1] - uox[m-1])/A1[m-1]
        vox[m] = -vox[m-1]/A1[m-1]
        uRed[m] = (bRed[m-1] - uRed[m-1])/A1[m-1]
        vRed[m] = -vRed[m-1]/A1[m-1]
      }

      B = matrix(data = c((Par$Kf[i1]*Par$h -  sum(vox*Derv(npoints = DerApprox, CoefMat = T))),
                          -Par$Kb[i1]*Par$h,
                          sum(vox*Derv(npoints = DerApprox, CoefMat = T)),
                          sum(vRed*Derv(npoints = DerApprox, CoefMat = T))),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(uox*Derv(npoints = DerApprox, CoefMat = T)),
                          -sum(uox*Derv(npoints = DerApprox, CoefMat = T)) - sum(uRed*Derv(npoints = DerApprox, CoefMat = T))))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]


      for (j1 in 1:(Par$j-2)) {
        Ox[i1+1,j1+1] = (bOx[j1] -Ox[i1+1,j1])/A1[j1]
        Red[i1+1,j1+1] = (bRed[j1] -Red[i1+1,j1])/A1[j1]
      }
    }

    Jox = Derv(Ox = Ox, h = Par$h, npoints = DerApprox)


  } else if (Method == "BDF") {

    al1 = 1/(Par$h^2)
    al2 = -2/(Par$h^2)
    al3 = 1/(Par$h^2)
    a1 = (al2 - 1.5/Par$dtn)/al1
    a2 = al3/al1

    for (i1 in 1:(Par$l-1)) {

      B = matrix(data = c((Par$Kf[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Par$Kb[i1]*Par$h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2] - Par$KC*Red[i1,1]

      if (i1 == 1) {
        bOx = -2*Ox[i1,2:(Par$j-1)]/(Par$dtn*al1) + Ox[1,2:(Par$j-1)]/(2*Par$dtn*al1)
        bRed = -2*Red[i1,2:(Par$j-1)]/(Par$dtn*al1) + Red[1,2:(Par$j-1)]/(2*Par$dtn*al1) + Par$KC*Red[i1,2:(Par$j-1)]/al1
      } else {
        bOx = -2*Ox[i1,2:(Par$j-1)]/(Par$dtn*al1) + Ox[i1-1,2:(Par$j-1)]/(2*Par$dtn*al1)
        bRed = -2*Red[i1,2:(Par$j-1)]/(Par$dtn*al1) + Red[i1-1,2:(Par$j-1)]/(2*Par$dtn*al1) + Par$KC*Red[i1,2:(Par$j-1)]/al1
      }

      A = c(rep(1,Par$j-2))
      A1 = c(rep(a1,Par$j-2))
      A2 = c(rep(a2,Par$j-2))
      bOx[Par$j-2] = bOx[Par$j-2] - A2[Par$j-2]*Ox[i1,Par$j]
      bRed[Par$j-2] = bRed[Par$j-2] - A2[Par$j-2]*(Red[i1,Par$j] - Par$KC*Red[i1,Par$j])
      uox = c(rep(0, DerApprox))
      vox = c(rep(1, DerApprox))
      uRed = c(rep(0, DerApprox))
      vRed = c(rep(1, DerApprox))

      for (j1 in ((Par$j-3):1)) {

        bOx[j1] = bOx[j1] - A2[Par$j-2]*bOx[j1+1]/A1[j1+1]
        bRed[j1] = bRed[j1] - A2[Par$j-2]*bRed[j1+1]/A1[j1+1]
        A1[j1] = A1[j1] - A2[Par$j-2]/A1[j1+1]

      }

      for (m in 2:DerApprox) {
        uox[m] = (bOx[m-1] - uox[m-1])/A1[m-1]
        vox[m] = -vox[m-1]/A1[m-1]
        uRed[m] = (bRed[m-1] - uRed[m-1])/A1[m-1]
        vRed[m] = -vRed[m-1]/A1[m-1]
      }

      B = matrix(data = c((Par$Kf[i1]*Par$h -  sum(vox*Derv(npoints = DerApprox, CoefMat = T))),
                          -Par$Kb[i1]*Par$h,
                          sum(vox*Derv(npoints = DerApprox, CoefMat = T)),
                          sum(vRed*Derv(npoints = DerApprox, CoefMat = T))),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(uox*Derv(npoints = DerApprox, CoefMat = T)),
                          -sum(uox*Derv(npoints = DerApprox, CoefMat = T)) - sum(uRed*Derv(npoints = DerApprox, CoefMat = T))))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]


      for (j1 in 1:(Par$j-2)) {
        Ox[i1+1,j1+1] = (bOx[j1] -Ox[i1+1,j1])/A1[j1]
        Red[i1+1,j1+1] = (bRed[j1] -Red[i1+1,j1])/A1[j1]
      }
    }

    Jox = Derv(Ox = Ox, h = Par$h, npoints = DerApprox)


  } else if (!(Method %in% c("Euler", "BI", "RK4", "CN", "BDF"))) {
    return("Available methods are Euler, BI, RK4, CN and BDF")
  }

  G = Jox
  i = (n*Par$FA*G*Dx*Area*Co)/(sqrt(Dx*Par$tau))

  graphy = ggplot(data = data.frame(i[1:(length(i)-1)], Par$PotentialScan[1:(length(i)-1)]),
                  aes(y = i[1:(length(i)-1)], x = Par$PotentialScan[1:(length(i)-1)])) +
    geom_point() + scale_x_continuous(trans = "reverse") +
    xlab("E / V") +
    ylab("I / A") +
    theme_classic()

  if (errCheck == TRUE){
    return(list(G,Dx,Co,Par$dtn,Par$h,Par$l,Par$j,i,n,Area,Par$p,Par$Da))
  } else {
    return(graphy)
  }
}

#' EE behaviour cyclic voltammetry simulator
#'
#' Return a graph I vs E of the electrochemical process
#'
#' @param Co bulk concentration expressed in Molar
#' @param Dx1 diffusion coefficient of the oxidized species expressed in cm^2/s
#' @param Eo1 reduction potential of the first electrochemical reaction expressed in Volt
#' @param Vi initial potential of the sweep expressed in Volt
#' @param Vf final potential of the sweep expressed in Volt
#' @param Vs  potential scan rate of the simulation expressed in V/s
#' @param ko1 heterogeneous electron transfer rate constant of the first electrochemical reaction expressed in m/s
#' @param alpha1 charge transfer coefficient of the first electrochemical reaction
#' @param Dred diffusion coefficient of the first reduced species expressed in cm^2/s
#' @param Dred2 diffusion coefficient of the second reduced species expressed in cm^2/s
#' @param Eo2 reduction potential of the second electrochemical reaction expressed in Volt
#' @param ko2 heterogeneous electron transfer rate constant of the second electrochemical reaction expressed in m/s
#' @param alpha2 charge transfer coefficient of the second electrochemical reaction
#' @param Dm simulation parameter, maximum 0.5 for explicit methods
#' @param Temp temperature in kelvin
#' @param n number of electrons involved in the process
#' @param Area area of the electrode expressed in cm^2
#' @param l number of time steps of the simulation
#' @param DerApprox number of point for the approximation of the first derivative
#' @param errCheck if true the function returns a list with parameters for CottrCheck function
#' @param Method method to be used for the simulation = "Euler" "BI" "RK4" "CN "BDF"
#'
#'
#' @return if errCheck == F a graph I vs E, if errCheck == T a list
#'
#' @examples
#' CVEE(Co = 0.001, DerApprox = 2, Dm = 0.45, errCheck = FALSE, Method = "Euler")
#' CVEE(Co = 0.001, Eo2 = -0.15, Dm = 0.45)
#'
#' @export
#' @import ggplot2
#'

CVEE = function(Co = 0.001, Dx1 = 0.00001, Eo1 = 0,
                Vi = 0.3, Vf = -0.3, Vs = 0.001, ko1 = 0.01,
                alpha1 = 0.5, Dred = 0.00001, Dred2 = 0.00001, Eo2 = 0,
                ko2 = 0.01, alpha2 = 0.5, Dm = 0.45, l = 100,
                Temp = 298.15, n = 1, Area = 1,
                DerApprox = 2, errCheck = FALSE, Method = "Euler") {

  Par = ParCall("CVEE", n. = n, Temp. = Temp, Dx1. = Dx1, Dred1. = Dred,
                Eo1. = Eo1, Eo2. = Eo2, Dm. = Dm, Vi. = Vi,
                Vf. = Vf, Vs. = Vs, ko1. = ko1, ko2. = ko2,
                alpha1. = alpha1, alpha2. = alpha2, Dred2. = Dred2, l. = l)
  Ox = OneMat(Par$l +1, Par$j)
  Red1 = ZeroMat(Par$l +1, Par$j)
  Jox = ZeroMat(Par$l+1, 1)
  Red2 = ZeroMat(Par$l +1, Par$j)
  Jred1 = ZeroMat(Par$l+1, 1)

  if (Method == "Euler") {

    for (i1 in 1:Par$l) {
      B = matrix(data = c((Par$Kf1[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Par$Kb1[i1]*Par$h, 0, -Par$Kf1[i1]*Par$h, Par$Kb1[i1]*Par$h + Par$Kf2[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1], -Par$Kb2[i1]*Par$h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 3, ncol = 3)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]),
                          Par$DRED*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red1[i1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - Par$DRED*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red1[i1,2:DerApprox]) - Par$DRED2*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red2[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red1[i1,1] = C[2]
      Red2[i1,1] = C[3]

      for (j1 in 2:(Par$j-1)) {
        Ox[i1+1,j1] = Ox[i1,j1] + Dm*(Ox[i1, j1-1] - 2*Ox[i1, j1] + Ox[i1, j1+1])
        Red1[i1+1, j1] = Red1[i1, j1] + Par$DRED*Dm*(Red1[i1, j1-1] + Red1[i1, j1+1] - 2*Red1[i1,j1])
        Red2[i1+1, j1] = Red2[i1, j1] + Par$DRED2*Dm*(Red2[i1, j1-1] + Red2[i1, j1+1] - 2*Red2[i1,j1])
      }
    }

    Jox = Derv(Ox = Ox, h = Par$h, npoints = DerApprox)
    Jred1 = Derv(Ox = Red1, h = Par$h, npoints = DerApprox)

  } else if (Method == "RK4") {

    for (i1 in 1:(Par$l-1)) {
      k1 = ZeroMat(Par$j)
      k2 = ZeroMat(Par$j)
      k3 = ZeroMat(Par$j)
      k4 = ZeroMat(Par$j)
      k1red = ZeroMat(Par$j)
      k2red = ZeroMat(Par$j)
      k3red = ZeroMat(Par$j)
      k4red = ZeroMat(Par$j)
      k1red2 = ZeroMat(Par$j)
      k2red2 = ZeroMat(Par$j)
      k3red2 = ZeroMat(Par$j)
      k4red2 = ZeroMat(Par$j)
      B = matrix(data = c((Par$Kf1[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Par$Kb1[i1]*Par$h, 0, -Par$Kf1[i1]*Par$h, Par$Kb1[i1]*Par$h + Par$Kf2[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1], -Par$Kb2[i1]*Par$h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 3, ncol = 3)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]),
                          Par$DRED*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red1[i1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - Par$DRED*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red1[i1,2:DerApprox]) - Par$DRED2*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red2[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red1[i1,1] = C[2]
      Red2[i1,1] = C[3]

      for (j1 in 2:(Par$j-1)) {
        k1[j1] = Dm*(Ox[i1, j1 -1] - 2*Ox[i1, j1] + Ox[i1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + k1[j1]*0.5
        k1red[j1] = Par$DRED*Dm*(Red1[i1, j1 -1] - 2*Red1[i1, j1] + Red1[i1, j1+1])
        Red1[i1 + 1,j1] = Red1[i1,j1] + k1red[j1]*0.5
        k1red2[j1] = Par$DRED*Dm*(Red2[i1, j1 -1] - 2*Red2[i1, j1] + Red2[i1, j1+1])
        Red2[i1 + 1,j1] = Red2[i1,j1] + k1red2[j1]*0.5
      }

      B = matrix(data = c((Par$Kf1[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Par$Kb1[i1]*Par$h, 0, -Par$Kf1[i1]*Par$h, Par$Kb1[i1]*Par$h + Par$Kf2[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1], -Par$Kb2[i1]*Par$h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 3, ncol = 3)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]),
                          Par$DRED*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red1[i1+1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]) - Par$DRED*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red1[i1+1,2:DerApprox]) - Par$DRED2*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red2[i1+1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red1[i1+1,1] = C[2]
      Red2[i1+1,1] = C[3]

      for (j1 in 2:(Par$j-1)) {
        k2[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + k2[j1]*0.5
        k2red[j1] = Par$DRED*Dm*(Red1[i1+1, j1 -1] - 2*Red1[i1+1, j1] + Red1[i1+1, j1+1])
        Red1[i1 + 1,j1] = Red1[i1,j1] + k2red[j1]*0.5
        k2red2[j1] = Par$DRED*Dm*(Red2[i1+1, j1 -1] - 2*Red2[i1+1, j1] + Red2[i1+1, j1+1])
        Red2[i1 + 1,j1] = Red2[i1,j1] + k2red2[j1]*0.5
      }

      B = matrix(data = c((Par$Kf1[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Par$Kb1[i1]*Par$h, 0, -Par$Kf1[i1]*Par$h, Par$Kb1[i1]*Par$h + Par$Kf2[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1], -Par$Kb2[i1]*Par$h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 3, ncol = 3)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]),
                          Par$DRED*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red1[i1+1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]) - Par$DRED*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red1[i1+1,2:DerApprox]) - Par$DRED2*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red2[i1+1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red1[i1+1,1] = C[2]
      Red2[i1+1,1] = C[3]

      for (j1 in 2:(Par$j-1)) {
        k3[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + k3[j1]
        k3red[j1] = Par$DRED*Dm*(Red1[i1 + 1, j1 -1] - 2*Red1[i1 + 1, j1] + Red1[i1 + 1, j1+1])
        Red1[i1 + 1,j1] = Red1[i1,j1] + k3red[j1]
        k3red2[j1] = Par$DRED*Dm*(Red2[i1 + 1, j1 -1] - 2*Red2[i1 + 1, j1] + Red2[i1 + 1, j1+1])
        Red2[i1 + 1,j1] = Red2[i1,j1] + k3red2[j1]
      }

      B = matrix(data = c((Par$Kf1[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Par$Kb1[i1]*Par$h, 0, -Par$Kf1[i1]*Par$h, Par$Kb1[i1]*Par$h + Par$Kf2[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1], -Par$Kb2[i1]*Par$h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 3, ncol = 3)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]),
                          Par$DRED*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red1[i1+1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]) - Par$DRED*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red1[i1+1,2:DerApprox]) - Par$DRED2*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red2[i1+1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red1[i1+1,1] = C[2]
      Red2[i1+1,1] = C[3]

      for (j1 in 2:(Par$j-1)) {
        k4[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + (k1[j1] + 2*k2[j1] + 2*k3[j1] + k4[j1])/6
        k4red[j1] = Par$DRED*Dm*(Red1[i1 + 1, j1 -1] - 2*Red1[i1 + 1, j1] + Red1[i1 + 1, j1+1])
        Red1[i1 + 1,j1] = Red1[i1,j1] + (k1red[j1] + 2*k2red[j1] + 2*k3red[j1] + k4red[j1])/6
        k4red2[j1] = Par$DRED*Dm*(Red2[i1 + 1, j1 -1] - 2*Red2[i1 + 1, j1] + Red2[i1 + 1, j1+1])
        Red2[i1 + 1,j1] = Red2[i1,j1] + (k1red2[j1] + 2*k2red2[j1] + 2*k3red2[j1] + k4red2[j1])/6
      }
    }

    Jox = Derv(Ox = Ox, h = Par$h, npoints = DerApprox)
    Jred1 = Derv(Ox = Red1, h = Par$h, npoints = DerApprox)

  } else if (Method == "BI") {
    al1 = 1/(Par$h^2)
    al2 = -2/(Par$h^2)
    al3 = 1/(Par$h^2)
    a1 = (al2 - 1/Par$dtn)/al1
    a2 = al3/al1
    al1red = Par$DRED/(Par$h^2)
    al2red = -(2*Par$DRED)/(Par$h^2)
    al3red = Par$DRED/(Par$h^2)
    a1red = (al2red - 1/Par$dtn)/al1red
    a2red = al3red/al1red
    al1red2 = Par$DRED2/(Par$h^2)
    al2red2 = -(2*Par$DRED2)/(Par$h^2)
    al3red2 = Par$DRED2/(Par$h^2)
    a1red2 = (al2red2 - 1/Par$dtn)/al1red2
    a2red2 = al3red2/al1red2

    for (i1 in 1:(Par$l-1)) {

      B = matrix(data = c((Par$Kf1[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Par$Kb1[i1]*Par$h, 0, -Par$Kf1[i1]*Par$h, Par$Kb1[i1]*Par$h + Par$Kf2[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1], -Par$Kb2[i1]*Par$h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 3, ncol = 3)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]),
                          Par$DRED*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red1[i1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - Par$DRED*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red1[i1,2:DerApprox]) - Par$DRED2*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red2[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red1[i1,1] = C[2]
      Red2[i1,1] = C[3]

      bOx = (-Ox[i1,(2:(Par$j-1))]/(al1*Par$dtn))
      bRed = (-Red1[i1,(2:(Par$j-1))]/(al1red*Par$dtn))
      bRed2 = (-Red2[i1,(2:(Par$j-1))]/(al1red2*Par$dtn))
      A = c(rep(1,Par$j-2))
      A1 = c(rep(a1,Par$j-2))
      A2 = c(rep(a2,Par$j-2))
      A1red = c(rep(a1red,Par$j-2))
      A2red = c(rep(a2red,Par$j-2))
      A1red2 = c(rep(a1red2,Par$j-2))
      A2red2 = c(rep(a2red2,Par$j-2))

      bOx[Par$j-2] = bOx[Par$j-2] - A2[Par$j-2]*1
      bRed[Par$j-2] = bRed[Par$j-2] - A2red[Par$j-2]*Red1[i1,Par$j]
      bRed2[Par$j-2] = bRed2[Par$j-2] - A2red2[Par$j-2]*Red2[i1,Par$j]
      uox = c(rep(0, DerApprox))
      vox = c(rep(1, DerApprox))
      uRed = c(rep(0, DerApprox))
      vRed = c(rep(1, DerApprox))
      uRed2 = c(rep(0, DerApprox))
      vRed2 = c(rep(1, DerApprox))

      for (j1 in ((Par$j-3):1)) {

        bOx[j1] = bOx[j1] - A2[Par$j-2]*bOx[j1+1]/A1[j1+1]
        bRed[j1] = bRed[j1] - A2red[Par$j-2]*bRed[j1+1]/A1red[j1+1]
        bRed2[j1] = bRed2[j1] - A2red2[Par$j-2]*bRed2[j1+1]/A1red2[j1+1]
        A1[j1] = A1[j1] - A2[Par$j-2]/A1[j1+1]
        A1red[j1] = A1red[j1] - A2red[Par$j-2]/A1red[j1+1]
        A1red2[j1] = A1red2[j1] - A2red2[Par$j-2]/A1red2[j1+1]

      }

      for (m in 2:DerApprox) {

        uox[m] = (bOx[m-1] - uox[m-1])/A1[m-1]
        vox[m] = -vox[m-1]/A1[m-1]
        uRed[m] = (bRed[m-1] - uRed[m-1])/A1red[m-1]
        vRed[m] = -vRed[m-1]/A1red[m-1]
        uRed2[m] = (bRed2[m-1] - uRed2[m-1])/A1red2[m-1]
        vRed2[m] = -vRed2[m-1]/A1red2[m-1]

      }

      B = matrix(data = c((Par$Kf1[i1]*Par$h - sum(vox*Derv(npoints = DerApprox, CoefMat = T))),
                          -Par$Kb1[i1]*Par$h, 0, -Par$Kf1[i1]*Par$h, Par$Kb1[i1]*Par$h + Par$Kf2[i1]*Par$h - sum(vRed*Derv(npoints = DerApprox, CoefMat = T)), -Par$Kb2[i1]*Par$h,
                          sum(vox*Derv(npoints = DerApprox, CoefMat = T)),
                          sum(vRed*Derv(npoints = DerApprox, CoefMat = T)),
                          sum(vRed2*Derv(npoints = DerApprox, CoefMat = T))),
                 byrow = T, nrow = 3, ncol = 3)
      Y = matrix(data = c(sum(uox*Derv(npoints = DerApprox, CoefMat = T)), sum(uRed*Derv(npoints = DerApprox, CoefMat = T)),
                          -sum(uox*Derv(npoints = DerApprox, CoefMat = T)) - sum(uRed*Derv(npoints = DerApprox, CoefMat = T)) - sum(uRed2*Derv(npoints = DerApprox, CoefMat = T))))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red1[i1+1,1] = C[2]
      Red2[i1+1,1] = C[3]


      for (j1 in 1:(Par$j-2)) {
        Ox[i1+1,j1+1] = (bOx[j1] - Ox[i1+1,j1])/A1[j1]
        Red1[i1+1,j1+1] = (bRed[j1] - Red1[i1+1,j1])/A1red[j1]
        Red2[i1+1,j1+1] = (bRed2[j1] - Red2[i1+1,j1])/A1red2[j1]
      }
    }

    Jox = Derv(Ox = Ox, h = Par$h, npoints = DerApprox)
    Jred1 = Derv(Ox = Red1, h = Par$h, npoints = DerApprox)

  } else if (Method == "CN") {

    al1 = 1/(Par$h^2)
    al2 = -2/(Par$h^2)
    al3 = 1/(Par$h^2)
    a1 = (al2 - 2/Par$dtn)/al1
    a2 = al3/al1
    a3 = (al2 + 2/Par$dtn)/al1
    al1red = Par$DRED/(Par$h^2)
    al2red = -(2*Par$DRED)/(Par$h^2)
    al3red = Par$DRED/(Par$h^2)
    a1red = (al2red - 2/Par$dtn)/al1red
    a2red = al3red/al1red
    a3red = (al2red + 2/Par$dtn)/al1red
    al1red2 = Par$DRED2/(Par$h^2)
    al2red2 = -(2*Par$DRED2)/(Par$h^2)
    al3red2 = Par$DRED2/(Par$h^2)
    a1red2 = (al2red2 - 2/Par$dtn)/al1red2
    a2red2 = al3red2/al1red2
    a3red2 = (al2red2 + 2/Par$dtn)/al1red2

    for (i1 in 1:(Par$l-1)) {

      B = matrix(data = c((Par$Kf1[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Par$Kb1[i1]*Par$h, 0, -Par$Kf1[i1]*Par$h, Par$Kb1[i1]*Par$h + Par$Kf2[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1], -Par$Kb2[i1]*Par$h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 3, ncol = 3)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]),
                          Par$DRED*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red1[i1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - Par$DRED*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red1[i1,2:DerApprox]) - Par$DRED2*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red2[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red1[i1,1] = C[2]
      Red2[i1,1] = C[3]

      bOx = -a3*Ox[i1,(2:(Par$j-1))] - Ox[i1,(1:(Par$j-2))] - a2*Ox[i1,(3:Par$j)]
      bRed = -a3red*Red1[i1,(2:(Par$j-1))]- Red1[i1,(1:(Par$j-2))] - a2red*Red1[i1,(3:Par$j)]
      bRed2 = -a3red2*Red2[i1,(2:(Par$j-1))]- Red2[i1,(1:(Par$j-2))] - a2red2*Red2[i1,(3:Par$j)]
      A = c(rep(1,Par$j-2))
      A1 = c(rep(a1,Par$j-2))
      A2 = c(rep(a2,Par$j-2))
      A1red = c(rep(a1red,Par$j-2))
      A2red = c(rep(a2red,Par$j-2))
      A1red2 = c(rep(a1red2,Par$j-2))
      A2red2 = c(rep(a2red2,Par$j-2))

      bOx[Par$j-2] = bOx[Par$j-2] - A2[Par$j-2]*1
      bRed[Par$j-2] = bRed[Par$j-2] - A2red[Par$j-2]*Red1[i1,Par$j]
      bRed2[Par$j-2] = bRed2[Par$j-2] - A2red2[Par$j-2]*Red2[i1,Par$j]
      uox = c(rep(0, DerApprox))
      vox = c(rep(1, DerApprox))
      uRed = c(rep(0, DerApprox))
      vRed = c(rep(1, DerApprox))
      uRed2 = c(rep(0, DerApprox))
      vRed2 = c(rep(1, DerApprox))

      for (j1 in ((Par$j-3):1)) {

        bOx[j1] = bOx[j1] - A2[Par$j-2]*bOx[j1+1]/A1[j1+1]
        bRed[j1] = bRed[j1] - A2red[Par$j-2]*bRed[j1+1]/A1red[j1+1]
        bRed2[j1] = bRed2[j1] - A2red2[Par$j-2]*bRed2[j1+1]/A1red2[j1+1]
        A1[j1] = A1[j1] - A2[Par$j-2]/A1[j1+1]
        A1red[j1] = A1red[j1] - A2red[Par$j-2]/A1red[j1+1]
        A1red2[j1] = A1red2[j1] - A2red2[Par$j-2]/A1red2[j1+1]

      }

      for (m in 2:DerApprox) {

        uox[m] = (bOx[m-1] - uox[m-1])/A1[m-1]
        vox[m] = -vox[m-1]/A1[m-1]
        uRed[m] = (bRed[m-1] - uRed[m-1])/A1red[m-1]
        vRed[m] = -vRed[m-1]/A1red[m-1]
        uRed2[m] = (bRed2[m-1] - uRed2[m-1])/A1red2[m-1]
        vRed2[m] = -vRed2[m-1]/A1red2[m-1]

      }

      B = matrix(data = c((Par$Kf1[i1]*Par$h - sum(vox*Derv(npoints = DerApprox, CoefMat = T))),
                          -Par$Kb1[i1]*Par$h, 0, -Par$Kf1[i1]*Par$h, Par$Kb1[i1]*Par$h + Par$Kf2[i1]*Par$h - sum(vRed*Derv(npoints = DerApprox, CoefMat = T)), -Par$Kb2[i1]*Par$h,
                          sum(vox*Derv(npoints = DerApprox, CoefMat = T)),
                          sum(vRed*Derv(npoints = DerApprox, CoefMat = T)),
                          sum(vRed2*Derv(npoints = DerApprox, CoefMat = T))),
                 byrow = T, nrow = 3, ncol = 3)
      Y = matrix(data = c(sum(uox*Derv(npoints = DerApprox, CoefMat = T)), sum(uRed*Derv(npoints = DerApprox, CoefMat = T)),
                          -sum(uox*Derv(npoints = DerApprox, CoefMat = T)) - sum(uRed*Derv(npoints = DerApprox, CoefMat = T)) - sum(uRed2*Derv(npoints = DerApprox, CoefMat = T))))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red1[i1+1,1] = C[2]
      Red2[i1+1,1] = C[3]


      for (j1 in 1:(Par$j-2)) {
        Ox[i1+1,j1+1] = (bOx[j1] - Ox[i1+1,j1])/A1[j1]
        Red1[i1+1,j1+1] = (bRed[j1] - Red1[i1+1,j1])/A1red[j1]
        Red2[i1+1,j1+1] = (bRed2[j1] - Red2[i1+1,j1])/A1red2[j1]
      }
    }

    Jox = Derv(Ox = Ox, h = Par$h, npoints = DerApprox)
    Jred1 = Derv(Ox = Red1, h = Par$h, npoints = DerApprox)


  } else if (Method == "BDF") {

    al1 = 1/(Par$h^2)
    al2 = -2/(Par$h^2)
    al3 = 1/(Par$h^2)
    a1 = (al2 - 1.5/Par$dtn)/al1
    a2 = al3/al1
    al1red = Par$DRED/(Par$h^2)
    al2red = -(2*Par$DRED)/(Par$h^2)
    al3red = Par$DRED/(Par$h^2)
    a1red = (al2red - 1.5/Par$dtn)/al1red
    a2red = al3red/al1red
    al1red2 = Par$DRED2/(Par$h^2)
    al2red2 = -(2*Par$DRED2)/(Par$h^2)
    al3red2 = Par$DRED2/(Par$h^2)
    a1red2 = (al2red2 - 1.5/Par$dtn)/al1red2
    a2red2 = al3red2/al1red2

    for (i1 in 1:(Par$l-1)) {

      B = matrix(data = c((Par$Kf1[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Par$Kb1[i1]*Par$h, 0, -Par$Kf1[i1]*Par$h, Par$Kb1[i1]*Par$h + Par$Kf2[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1], -Par$Kb2[i1]*Par$h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 3, ncol = 3)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]),
                          Par$DRED*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red1[i1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - Par$DRED*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red1[i1,2:DerApprox]) - Par$DRED2*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red2[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red1[i1,1] = C[2]
      Red2[i1,1] = C[3]

      if (i1 == 1) {
        bOx = -2*Ox[i1,2:(Par$j-1)]/(Par$dtn*al1) + Ox[1,2:(Par$j-1)]/(2*Par$dtn*al1)
        bRed = -2*Red1[i1,2:(Par$j-1)]/(Par$dtn*al1red) + Red1[1,2:(Par$j-1)]/(2*Par$dtn*al1red)
        bRed2 = -2*Red2[i1,2:(Par$j-1)]/(Par$dtn*al1red2) + Red2[1,2:(Par$j-1)]/(2*Par$dtn*al1red2)
      } else {
        bOx = -2*Ox[i1,2:(Par$j-1)]/(Par$dtn*al1) + Ox[i1-1,2:(Par$j-1)]/(2*Par$dtn*al1)
        bRed = -2*Red1[i1,2:(Par$j-1)]/(Par$dtn*al1red) + Red1[i1-1,2:(Par$j-1)]/(2*Par$dtn*al1red)
        bRed2 = -2*Red2[i1,2:(Par$j-1)]/(Par$dtn*al1red2) + Red2[i1-1,2:(Par$j-1)]/(2*Par$dtn*al1red2)
      }

      A = c(rep(1,Par$j-2))
      A1 = c(rep(a1,Par$j-2))
      A2 = c(rep(a2,Par$j-2))
      A1red = c(rep(a1red,Par$j-2))
      A2red = c(rep(a2red,Par$j-2))
      A1red2 = c(rep(a1red2,Par$j-2))
      A2red2 = c(rep(a2red2,Par$j-2))

      bOx[Par$j-2] = bOx[Par$j-2] - A2[Par$j-2]*1
      bRed[Par$j-2] = bRed[Par$j-2] - A2red[Par$j-2]*Red1[i1,Par$j]
      bRed2[Par$j-2] = bRed2[Par$j-2] - A2red2[Par$j-2]*Red2[i1,Par$j]
      uox = c(rep(0, DerApprox))
      vox = c(rep(1, DerApprox))
      uRed = c(rep(0, DerApprox))
      vRed = c(rep(1, DerApprox))
      uRed2 = c(rep(0, DerApprox))
      vRed2 = c(rep(1, DerApprox))

      for (j1 in ((Par$j-3):1)) {

        bOx[j1] = bOx[j1] - A2[Par$j-2]*bOx[j1+1]/A1[j1+1]
        bRed[j1] = bRed[j1] - A2red[Par$j-2]*bRed[j1+1]/A1red[j1+1]
        bRed2[j1] = bRed2[j1] - A2red2[Par$j-2]*bRed2[j1+1]/A1red2[j1+1]
        A1[j1] = A1[j1] - A2[Par$j-2]/A1[j1+1]
        A1red[j1] = A1red[j1] - A2red[Par$j-2]/A1red[j1+1]
        A1red2[j1] = A1red2[j1] - A2red2[Par$j-2]/A1red2[j1+1]

      }

      for (m in 2:DerApprox) {

        uox[m] = (bOx[m-1] - uox[m-1])/A1[m-1]
        vox[m] = -vox[m-1]/A1[m-1]
        uRed[m] = (bRed[m-1] - uRed[m-1])/A1red[m-1]
        vRed[m] = -vRed[m-1]/A1red[m-1]
        uRed2[m] = (bRed2[m-1] - uRed2[m-1])/A1red2[m-1]
        vRed2[m] = -vRed2[m-1]/A1red2[m-1]

      }

      B = matrix(data = c((Par$Kf1[i1]*Par$h - sum(vox*Derv(npoints = DerApprox, CoefMat = T))),
                          -Par$Kb1[i1]*Par$h, 0, -Par$Kf1[i1]*Par$h, Par$Kb1[i1]*Par$h + Par$Kf2[i1]*Par$h - sum(vRed*Derv(npoints = DerApprox, CoefMat = T)), -Par$Kb2[i1]*Par$h,
                          sum(vox*Derv(npoints = DerApprox, CoefMat = T)),
                          sum(vRed*Derv(npoints = DerApprox, CoefMat = T)),
                          sum(vRed2*Derv(npoints = DerApprox, CoefMat = T))),
                 byrow = T, nrow = 3, ncol = 3)
      Y = matrix(data = c(sum(uox*Derv(npoints = DerApprox, CoefMat = T)), sum(uRed*Derv(npoints = DerApprox, CoefMat = T)),
                          -sum(uox*Derv(npoints = DerApprox, CoefMat = T)) - sum(uRed*Derv(npoints = DerApprox, CoefMat = T)) - sum(uRed2*Derv(npoints = DerApprox, CoefMat = T))))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red1[i1+1,1] = C[2]
      Red2[i1+1,1] = C[3]


      for (j1 in 1:(Par$j-2)) {
        Ox[i1+1,j1+1] = (bOx[j1] - Ox[i1+1,j1])/A1[j1]
        Red1[i1+1,j1+1] = (bRed[j1] - Red1[i1+1,j1])/A1red[j1]
        Red2[i1+1,j1+1] = (bRed2[j1] - Red2[i1+1,j1])/A1red2[j1]
      }
    }

    Jox = Derv(Ox = Ox, h = Par$h, npoints = DerApprox)
    Jred1 = Derv(Ox = Red1, h = Par$h, npoints = DerApprox)

  } else if (!(Method %in% c("Euler", "BI", "RK4", "CN", "BDF"))) {
    return("Available methods are Euler, BI, RK4, CN and BDF")
  }


  G1 = Jox
  G2 = Jox + Jred1
  i = (n*Par$FA*(G1+G2)*Dx1*Area*Co)/(sqrt(Dx1*Par$tau))
  graphy = ggplot(data = data.frame(i[1:(length(i)-1)],Par$PotentialScan[1:(length(i)-1)]),
                  aes(y = i[1:(length(i)-1)], x = Par$PotentialScan[1:(length(i)-1)])) +
    geom_point() + scale_x_continuous(trans = "reverse") +
    xlab("E / V") +
    ylab("I / A") +
    theme_classic()

  if (errCheck == TRUE){
    return(list((G1+G2),Dx1,Dred,Dred2,Co,Par$dtn,Par$h,i,Par$l,Par$j,n,Area,Par$DOx,Par$DRED,Par$DRED2,Par$p1,Par$p2))
  } else {
    return(graphy)
  }
}


#' General Purpose CV simulation
#'
#' Return a graph I vs E of the electrochemical process, up to 4 EE mechanisms and CE mechanisms can be simulated
#'
#' @param Co bulk concentration oxidated speciesexpressed in Molar
#' @param Cred bulk concentration of reduced species expressed in Molar
#' @param kco Chemical rate constant for Ox Species expressed in s^-1
#' @param Dx1 diffusion coefficient of the oxidized species expressed in cm^2/s
#' @param Eo1 reduction potential of the first electrochemical reaction expressed in Volt
#' @param kc1 Chemical rate constant for Red Species expressed in s^-1
#' @param Vi initial potential of the sweep expressed in Volt
#' @param Vf final potential of the sweep expressed in Volt
#' @param Vs  potential scan rate of the simulation expressed in V/s
#' @param ko1 heterogeneous electron transfer rate constant of the first electrochemical reaction expressed in m/s
#' @param alpha1 charge transfer coefficient of the first electrochemical reaction
#' @param Dred diffusion coefficient of the first reduced species expressed in cm^2/S
#' @param Dred2 diffusion coefficient of the second reduced species expressed in cm^2/s
#' @param Eo2 reduction potential of the second electrochemical reaction expressed in Volt
#' @param kc2 Chemical rate constant for second Red Species expressed in s^-1
#' @param ko2 heterogeneous electron transfer rate constant of the second electrochemical reaction expressed in m/s
#' @param alpha2 charge transfer coefficient of the second electrochemical reaction
#' @param Dred3 diffusion coefficient of the third reduced species expressed in cm^2/s
#' @param Dred4 diffusion coefficient of the fourth reduced species cm^2/s
#' @param alpha3 charge transfer coefficient of the third electrochemical reaction
#' @param alpha4 charge transfer coefficient of the fourth electrochemical reaction
#' @param kc3 Chemical rate constant for third Red Species expressed in s^-1
#' @param ko3 heterogeneous electron transfer rate constant of the third electrochemical reaction expressed in m/s
#' @param kc4 Chemical rate constant for fourth Red Species expressed in s^-1
#' @param Eo3 reduction potential of the third electrochemical reaction expressed in Volt
#' @param Eo4 reduction potential of the fourth electrochemical reaction expressed in Volt
#' @param ko4 heterogeneous electron transfer rate constant of the fourth electrochemical reaction expressed in m/s
#' @param Dm simulation parameter, maximum 0.5 for explicit methods
#' @param Temp temperature in kelvin
#' @param n number of electrons involved in the process
#' @param Area area of the electrode expressed in cm^2
#' @param l number of time steps of the simulation
#' @param DerApprox number of point for the approximation of the first derivative
#' @param errCheck if true the function returns a list with parameters for CottrCheck function
#' @param Method method to be used for the simulation = "Euler" "BI" "RK4" "CN "BDF"
#'
#'
#' @return if errCheck == F a graph I vs E, if errCheck == T a list
#'
#' @examples
#' Gen_CV(Co = 0.001, DerApprox = 2, Dm = 0.45, errCheck = FALSE, Method = "Euler")
#' Gen_CV(Co = 0.001, Eo2 = -0.15, Dm = 0.45, kc1 = 0.0001)
#'
#' @export
#' @import ggplot2
#'
#'
#'


Gen_CV = function(Co = 0.001, Cred= 0, kco = 0, Dx1 = 0.00001, Eo1 = 0, kc1 = 0,
                  Vi = 0.3, Vf = -0.3, Vs = 0.001, ko1 = 0.01,
                  alpha1 = 0.5, Dred = 0.00001, Dred2 = 0.00001, Eo2 = 0, kc2 = 0,
                  ko2 = 0, alpha2 = 0.5, Dm = 0.45, Dred3 = 0.00001,
                  Eo3 = 0, kc3 = 0, ko3 = 0, alpha3 = 0.5, Dred4 = 0.00001,
                  Eo4 = 0, kc4 = 0, ko4  = 0, alpha4 = 0.5,
                  Temp = 298.15, n = 1, Area = 1, l = 100,
                  DerApprox = 2, errCheck = FALSE, Method = "Euler") {
  if (kco > 0.001 | kc1 > 0.001 | kc2 > 0.001 | kc3 > 0.001 | kc4 > 0.001 ) {
    warning("Chemical rate costant is too high, this will result in unstable simulation")
  }
  Par = ParCall("Gen_CV", n. = n, Temp. = Temp, Dx1. = Dx1, Dred1. = Dred,
                Dred2. = Dred2, Dred3. = Dred3, Dred4. = Dred4,
                Eo1. = Eo1, Eo2. = Eo2, Eo3. = Eo3, Eo4. = Eo4, Dm. = Dm,
                Vi. = Vi, kco. = kco, kc1. = kc1, kc2. = kc2, kc3. = kc3, kc4. = kc4,
                Vf. = Vf, Vs. = Vs, ko1. = ko1, ko2. = ko2, ko3. = ko3, ko4. = ko4,
                alpha1. = alpha1, alpha2. = alpha2, alpha3. = alpha3, alpha4. = alpha4, l. = l)
  Ox = OneMat(Par$l +1, Par$j)
  if (Co == 0) {
    Co = 0.0000001
  }
  Red1 = (Cred/Co)*OneMat(Par$l +1, Par$j)
  Jox = ZeroMat(Par$l+1, 1)
  Red2 = ZeroMat(Par$l +1, Par$j)
  Jred1 = ZeroMat(Par$l+1, 1)
  Red3 = ZeroMat(Par$l +1, Par$j)
  Jred2 = ZeroMat(Par$l+1, 1)
  Red4 = ZeroMat(Par$l +1, Par$j)
  Jred3 = ZeroMat(Par$l+1, 1)

  if (Method == "Euler") {

    for (i1 in 1:Par$l) {
      B = matrix(data = c((Par$Kf1[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1]), -Par$Kb1[i1]*Par$h, 0, 0, 0,
                          -Par$Kf1[i1]*Par$h, Par$Kb1[i1]*Par$h + Par$Kf2[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1], -Par$Kb2[i1]*Par$h, 0, 0,
                          0, -Par$Kf2[i1]*Par$h, Par$Kb2[i1]*Par$h + Par$Kf3[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1], -Par$Kb3[i1]*Par$h, 0,
                          0, 0, -Par$Kf3[i1]*Par$h, Par$Kb3[i1]*Par$h + Par$Kf4[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1], -Par$Kb4[i1]*Par$h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 5, ncol = 5)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]),
                          Par$DRED*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red1[i1,2:DerApprox]),
                          Par$DRED2*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red2[i1,2:DerApprox]),
                          Par$DRED3*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red3[i1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - Par$DRED*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red1[i1,2:DerApprox]) - Par$DRED2*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red2[i1,2:DerApprox]) - Par$DRED3*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red3[i1,2:DerApprox]) - Par$DRED4*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red4[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1] - Par$KCo*Ox[i1,1]
      Red1[i1,1] = C[2] - Par$KC1*Red1[i1,1]
      Red2[i1,1] = C[3] - Par$KC2*Red2[i1,1]
      Red3[i1,1] = C[4] - Par$KC3*Red3[i1,1]
      Red4[i1,1] = C[5] - Par$KC4*Red4[i1,1]

      for (j1 in 2:(Par$j-1)) {
        Ox[i1+1,j1] = Ox[i1,j1] + Dm*(Ox[i1, j1-1] - 2*Ox[i1, j1] + Ox[i1, j1+1]) - Par$KCo*Par$dtn*Ox[i1,j1]
        Red1[i1+1, j1] = Red1[i1, j1] + Par$DRED*Dm*(Red1[i1, j1-1] + Red1[i1, j1+1] - 2*Red1[i1,j1]) - Par$KC1*Par$dtn*Red1[i1,j1]
        Red2[i1+1, j1] = Red2[i1, j1] + Par$DRED2*Dm*(Red2[i1, j1-1] + Red2[i1, j1+1] - 2*Red2[i1,j1]) - Par$KC2*Par$dtn*Red2[i1,j1]
        Red3[i1+1, j1] = Red3[i1, j1] + Par$DRED3*Dm*(Red3[i1, j1-1] + Red3[i1, j1+1] - 2*Red3[i1,j1]) - Par$KC3*Par$dtn*Red3[i1,j1]
        Red4[i1+1, j1] = Red4[i1, j1] + Par$DRED4*Dm*(Red4[i1, j1-1] + Red4[i1, j1+1] - 2*Red4[i1,j1]) - Par$KC4*Par$dtn*Red4[i1,j1]
      }
    }

    Jox = Derv(Ox = Ox, h = Par$h, npoints = DerApprox)
    Jred1 = Derv(Ox = Red1, h = Par$h, npoints = DerApprox)
    Jred2 = Derv(Ox = Red2, h = Par$h, npoints = DerApprox)
    Jred3 = Derv(Ox = Red3, h = Par$h, npoints = DerApprox)
    Jred4 = Derv(Ox = Red4, h = Par$h, npoints = DerApprox)

  } else if (Method == "RK4") {

    for (i1 in 1:(Par$l-1)) {
      k1 = ZeroMat(Par$j)
      k2 = ZeroMat(Par$j)
      k3 = ZeroMat(Par$j)
      k4 = ZeroMat(Par$j)
      k1red1 = ZeroMat(Par$j)
      k2red1 = ZeroMat(Par$j)
      k3red1 = ZeroMat(Par$j)
      k4red1 = ZeroMat(Par$j)
      k1red2 = ZeroMat(Par$j)
      k2red2 = ZeroMat(Par$j)
      k3red2 = ZeroMat(Par$j)
      k4red2 = ZeroMat(Par$j)
      k1red3 = ZeroMat(Par$j)
      k2red3 = ZeroMat(Par$j)
      k3red3 = ZeroMat(Par$j)
      k4red3 = ZeroMat(Par$j)
      k1red4 = ZeroMat(Par$j)
      k2red4 = ZeroMat(Par$j)
      k3red4 = ZeroMat(Par$j)
      k4red4 = ZeroMat(Par$j)

      B = matrix(data = c((Par$Kf1[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1]), -Par$Kb1[i1]*Par$h, 0, 0, 0,
                          -Par$Kf1[i1]*Par$h, Par$Kb1[i1]*Par$h + Par$Kf2[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1], -Par$Kb2[i1]*Par$h, 0, 0,
                          0, -Par$Kf2[i1]*Par$h, Par$Kb2[i1]*Par$h + Par$Kf3[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1], -Par$Kb3[i1]*Par$h, 0,
                          0, 0, -Par$Kf3[i1]*Par$h, Par$Kb3[i1]*Par$h + Par$Kf4[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1], -Par$Kb4[i1]*Par$h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 5, ncol = 5)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]),
                          Par$DRED*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red1[i1,2:DerApprox]),
                          Par$DRED2*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red2[i1,2:DerApprox]),
                          Par$DRED3*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red3[i1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - Par$DRED*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red1[i1,2:DerApprox]) - Par$DRED2*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red2[i1,2:DerApprox]) - Par$DRED3*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red3[i1,2:DerApprox]) - Par$DRED4*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red4[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red1[i1,1] = C[2]
      Red2[i1,1] = C[3]
      Red3[i1,1] = C[4]
      Red4[i1,1] = C[5]

      for (j1 in 2:(Par$j-1)) {
        k1[j1] = Dm*(Ox[i1, j1 -1] - 2*Ox[i1, j1] + Ox[i1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + k1[j1]*0.5
        k1red1[j1] = Par$DRED*Dm*(Red1[i1, j1 -1] - 2*Red1[i1, j1] + Red1[i1, j1+1])
        Red1[i1 + 1,j1] = Red1[i1,j1] + k1red1[j1]*0.5
        k1red2[j1] = Par$DRED2*Dm*(Red2[i1, j1 -1] - 2*Red2[i1, j1] + Red2[i1, j1+1])
        Red2[i1 + 1,j1] = Red2[i1,j1] + k1red2[j1]*0.5
        k1red3[j1] = Par$DRED3*Dm*(Red3[i1, j1 -1] - 2*Red3[i1, j1] + Red3[i1, j1+1])
        Red3[i1 + 1,j1] = Red3[i1,j1] + k1red3[j1]*0.5
        k1red4[j1] = Par$DRED4*Dm*(Red4[i1, j1 -1] - 2*Red4[i1, j1] + Red4[i1, j1+1])
        Red4[i1 + 1,j1] = Red4[i1,j1] + k1red4[j1]*0.5

      }

      B = matrix(data = c((Par$Kf1[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1]), -Par$Kb1[i1]*Par$h, 0, 0, 0,
                          -Par$Kf1[i1]*Par$h, Par$Kb1[i1]*Par$h + Par$Kf2[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1], -Par$Kb2[i1]*Par$h, 0, 0,
                          0, -Par$Kf2[i1]*Par$h, Par$Kb2[i1]*Par$h + Par$Kf3[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1], -Par$Kb3[i1]*Par$h, 0,
                          0, 0, -Par$Kf3[i1]*Par$h, Par$Kb3[i1]*Par$h + Par$Kf4[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1], -Par$Kb4[i1]*Par$h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 5, ncol = 5)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]),
                          Par$DRED*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red1[i1+1,2:DerApprox]),
                          Par$DRED2*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red2[i1+1,2:DerApprox]),
                          Par$DRED3*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red3[i1+1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]) - Par$DRED*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red1[i1+1,2:DerApprox]) - Par$DRED2*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red2[i1+1,2:DerApprox]) - Par$DRED3*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red3[i1+1,2:DerApprox]) - Par$DRED4*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red4[i1+1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red1[i1+1,1] = C[2]
      Red2[i1+1,1] = C[3]
      Red3[i1+1,1] = C[4]
      Red4[i1+1,1] = C[5]

      for (j1 in 2:(Par$j-1)) {
        k2[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + k2[j1]*0.5
        k2red1[j1] = Par$DRED*Dm*(Red1[i1+1, j1 -1] - 2*Red1[i1+1, j1] + Red1[i1+1, j1+1])
        Red1[i1 + 1,j1] = Red1[i1,j1] + k2red1[j1]*0.5
        k2red2[j1] = Par$DRED2*Dm*(Red2[i1+1, j1 -1] - 2*Red2[i1+1, j1] + Red2[i1+1, j1+1])
        Red2[i1 + 1,j1] = Red2[i1,j1] + k2red2[j1]*0.5
        k2red3[j1] = Par$DRED3*Dm*(Red3[i1+1, j1 -1] - 2*Red3[i1+1, j1] + Red3[i1+1, j1+1])
        Red3[i1 + 1,j1] = Red3[i1,j1] + k2red3[j1]*0.5
        k2red4[j1] = Par$DRED4*Dm*(Red4[i1+1, j1 -1] - 2*Red4[i1+1, j1] + Red4[i1+1, j1+1])
        Red4[i1 + 1,j1] = Red4[i1,j1] + k2red4[j1]*0.5

      }


      B = matrix(data = c((Par$Kf1[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1]), -Par$Kb1[i1]*Par$h, 0, 0, 0,
                          -Par$Kf1[i1]*Par$h, Par$Kb1[i1]*Par$h + Par$Kf2[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1], -Par$Kb2[i1]*Par$h, 0, 0,
                          0, -Par$Kf2[i1]*Par$h, Par$Kb2[i1]*Par$h + Par$Kf3[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1], -Par$Kb3[i1]*Par$h, 0,
                          0, 0, -Par$Kf3[i1]*Par$h, Par$Kb3[i1]*Par$h + Par$Kf4[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1], -Par$Kb4[i1]*Par$h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 5, ncol = 5)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]),
                          Par$DRED*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red1[i1+1,2:DerApprox]),
                          Par$DRED2*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red2[i1+1,2:DerApprox]),
                          Par$DRED3*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red3[i1+1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]) - Par$DRED*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red1[i1+1,2:DerApprox]) - Par$DRED2*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red2[i1+1,2:DerApprox]) - Par$DRED3*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red3[i1+1,2:DerApprox]) - Par$DRED4*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red4[i1+1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red1[i1+1,1] = C[2]
      Red2[i1+1,1] = C[3]
      Red3[i1+1,1] = C[4]
      Red4[i1+1,1] = C[5]

      for (j1 in 2:(Par$j-1)) {
        k3[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + k3[j1]
        k3red1[j1] = Par$DRED*Dm*(Red1[i1 + 1, j1 -1] - 2*Red1[i1 + 1, j1] + Red1[i1 + 1, j1+1])
        Red1[i1 + 1,j1] = Red1[i1,j1] + k3red1[j1]
        k3red2[j1] = Par$DRED2*Dm*(Red2[i1 + 1, j1 -1] - 2*Red2[i1 + 1, j1] + Red2[i1 + 1, j1+1])
        Red2[i1 + 1,j1] = Red2[i1,j1] + k3red2[j1]
        k3red3[j1] = Par$DRED3*Dm*(Red3[i1+1, j1 -1] - 2*Red3[i1+1, j1] + Red3[i1+1, j1+1])
        Red3[i1 + 1,j1] = Red3[i1,j1] + k3red3[j1]
        k3red4[j1] = Par$DRED4*Dm*(Red4[i1+1, j1 -1] - 2*Red4[i1+1, j1] + Red4[i1+1, j1+1])
        Red4[i1 + 1,j1] = Red4[i1,j1] + k3red4[j1]

      }

      B = matrix(data = c((Par$Kf1[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1]), -Par$Kb1[i1]*Par$h, 0, 0, 0,
                          -Par$Kf1[i1]*Par$h, Par$Kb1[i1]*Par$h + Par$Kf2[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1], -Par$Kb2[i1]*Par$h, 0, 0,
                          0, -Par$Kf2[i1]*Par$h, Par$Kb2[i1]*Par$h + Par$Kf3[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1], -Par$Kb3[i1]*Par$h, 0,
                          0, 0, -Par$Kf3[i1]*Par$h, Par$Kb3[i1]*Par$h + Par$Kf4[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1], -Par$Kb4[i1]*Par$h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 5, ncol = 5)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]),
                          Par$DRED*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red1[i1+1,2:DerApprox]),
                          Par$DRED2*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red2[i1+1,2:DerApprox]),
                          Par$DRED3*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red3[i1+1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]) - Par$DRED*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red1[i1+1,2:DerApprox]) - Par$DRED2*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red2[i1+1,2:DerApprox]) - Par$DRED3*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red3[i1+1,2:DerApprox]) - Par$DRED4*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red4[i1+1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1] - Par$KCo*Ox[i1+1,1]
      Red1[i1+1,1] = C[2] - Par$KC1*Red1[i1+1,1]
      Red2[i1+1,1] = C[3] - Par$KC2*Red2[i1+1,1]
      Red3[i1+1,1] = C[4] - Par$KC3*Red3[i1+1,1]
      Red4[i1+1,1] = C[5] - Par$KC4*Red4[i1+1,1]

      for (j1 in 2:(Par$j-1)) {
        k4[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + (k1[j1] + 2*k2[j1] + 2*k3[j1] + k4[j1])/6 - Par$KCo*Par$dtn*Ox[i1+1,j1]
        k4red1[j1] = Par$DRED*Dm*(Red1[i1 + 1, j1 -1] - 2*Red1[i1 + 1, j1] + Red1[i1 + 1, j1+1])
        Red1[i1 + 1,j1] = Red1[i1,j1] + (k1red1[j1] + 2*k2red1[j1] + 2*k3red1[j1] + k4red1[j1])/6 - Par$KC1*Par$dtn*Red1[i1+1,j1]
        k4red2[j1] = Par$DRED2*Dm*(Red2[i1 + 1, j1 -1] - 2*Red2[i1 + 1, j1] + Red2[i1 + 1, j1+1])
        Red2[i1 + 1,j1] = Red2[i1,j1] + (k1red2[j1] + 2*k2red2[j1] + 2*k3red2[j1] + k4red2[j1])/6 - Par$KC2*Par$dtn*Red2[i1+1,j1]
        k4red3[j1] = Par$DRED3*Dm*(Red3[i1 + 1, j1 -1] - 2*Red3[i1 + 1, j1] + Red3[i1 + 1, j1+1])
        Red3[i1 + 1,j1] = Red3[i1,j1] + (k1red3[j1] + 2*k2red3[j1] + 2*k3red3[j1] + k4red3[j1])/6 - Par$KC3*Par$dtn*Red3[i1+1,j1]
        k4red4[j1] = Par$DRED4*Dm*(Red4[i1 + 1, j1 -1] - 2*Red4[i1 + 1, j1] + Red4[i1 + 1, j1+1])
        Red4[i1 + 1,j1] = Red4[i1,j1] + (k1red4[j1] + 2*k2red4[j1] + 2*k3red4[j1] + k4red4[j1])/6 - Par$KC4*Par$dtn*Red4[i1+1,j1]
      }
    }

    Jox = Derv(Ox = Ox, h = Par$h, npoints = DerApprox)
    Jred1 = Derv(Ox = Red1, h = Par$h, npoints = DerApprox)
    Jred2 = Derv(Ox = Red2, h = Par$h, npoints = DerApprox)
    Jred3 = Derv(Ox = Red3, h = Par$h, npoints = DerApprox)
    Jred4 = Derv(Ox = Red4, h = Par$h, npoints = DerApprox)


  } else if (Method == "BI") {
    al1 = 1/(Par$h^2)
    al2 = -2/(Par$h^2)
    al3 = 1/(Par$h^2)
    a1 = (al2 - 1/Par$dtn)/al1
    a2 = al3/al1
    al1red = Par$DRED/(Par$h^2)
    al2red = -(2*Par$DRED)/(Par$h^2)
    al3red = Par$DRED/(Par$h^2)
    a1red = (al2red - 1/Par$dtn)/al1red
    a2red = al3red/al1red
    al1red2 = Par$DRED2/(Par$h^2)
    al2red2 = -(2*Par$DRED2)/(Par$h^2)
    al3red2 = Par$DRED2/(Par$h^2)
    a1red2 = (al2red2 - 1/Par$dtn)/al1red2
    a2red2 = al3red2/al1red2
    al1red3 = Par$DRED3/(Par$h^2)
    al2red3 = -(2*Par$DRED3)/(Par$h^2)
    al3red3 = Par$DRED3/(Par$h^2)
    a1red3 = (al2red3 - 1/Par$dtn)/al1red3
    a2red3 = al3red3/al1red3
    al1red4 = Par$DRED4/(Par$h^2)
    al2red4 = -(2*Par$DRED4)/(Par$h^2)
    al3red4 = Par$DRED4/(Par$h^2)
    a1red4 = (al2red4 - 1/Par$dtn)/al1red4
    a2red4 = al3red4/al1red4

    for (i1 in 1:(Par$l-1)) {

      B = matrix(data = c((Par$Kf1[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1]), -Par$Kb1[i1]*Par$h, 0, 0, 0,
                          -Par$Kf1[i1]*Par$h, Par$Kb1[i1]*Par$h + Par$Kf2[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1], -Par$Kb2[i1]*Par$h, 0, 0,
                          0, -Par$Kf2[i1]*Par$h, Par$Kb2[i1]*Par$h + Par$Kf3[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1], -Par$Kb3[i1]*Par$h, 0,
                          0, 0, -Par$Kf3[i1]*Par$h, Par$Kb3[i1]*Par$h + Par$Kf4[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1], -Par$Kb4[i1]*Par$h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 5, ncol = 5)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]),
                          Par$DRED*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red1[i1,2:DerApprox]),
                          Par$DRED2*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red2[i1,2:DerApprox]),
                          Par$DRED3*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red3[i1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - Par$DRED*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red1[i1,2:DerApprox]) - Par$DRED2*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red2[i1,2:DerApprox]) - Par$DRED3*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red3[i1,2:DerApprox]) - Par$DRED4*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red4[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1] - Par$KCo*Ox[i1,1]
      Red1[i1,1] = C[2] - Par$KC1*Red1[i1,1]
      Red2[i1,1] = C[3] - Par$KC2*Red2[i1,1]
      Red3[i1,1] = C[4] - Par$KC3*Red3[i1,1]
      Red4[i1,1] = C[5] - Par$KC4*Red4[i1,1]

      bOx = (-Ox[i1,(2:(Par$j-1))]/(al1*Par$dtn)) + Par$KCo*Ox[i1,2:(Par$j-1)]/al1
      bRed = (-Red1[i1,(2:(Par$j-1))]/(al1red*Par$dtn)) + Par$KC1*Red1[i1,2:(Par$j-1)]/al1red
      bRed2 = (-Red2[i1,(2:(Par$j-1))]/(al1red2*Par$dtn)) + Par$KC2*Red2[i1,2:(Par$j-1)]/al1red2
      bRed3 = (-Red3[i1,(2:(Par$j-1))]/(al1red3*Par$dtn)) + Par$KC3*Red3[i1,2:(Par$j-1)]/al1red3
      bRed4 = (-Red4[i1,(2:(Par$j-1))]/(al1red4*Par$dtn)) + Par$KC4*Red4[i1,2:(Par$j-1)]/al1red4

      A = c(rep(1,Par$j-2))
      A1 = c(rep(a1,Par$j-2))
      A2 = c(rep(a2,Par$j-2))
      A1red = c(rep(a1red,Par$j-2))
      A2red = c(rep(a2red,Par$j-2))
      A1red2 = c(rep(a1red2,Par$j-2))
      A2red2 = c(rep(a2red2,Par$j-2))
      A1red3 = c(rep(a1red3,Par$j-2))
      A2red3 = c(rep(a2red3,Par$j-2))
      A1red4 = c(rep(a1red4,Par$j-2))
      A2red4 = c(rep(a2red4,Par$j-2))

      bOx[Par$j-2] = bOx[Par$j-2] - A2[Par$j-2]*(1 - Par$KCo*Ox[i1,Par$j])
      bRed[Par$j-2] = bRed[Par$j-2] - A2red[Par$j-2]*(Red1[i1,Par$j] - Par$KC1*Red1[i1,Par$j])
      bRed2[Par$j-2] = bRed2[Par$j-2] - A2red2[Par$j-2]*(Red2[i1,Par$j] - Par$KC2*Red2[i1,Par$j])
      bRed3[Par$j-2] = bRed3[Par$j-2] - A2red3[Par$j-2]*(Red3[i1,Par$j] - Par$KC3*Red3[i1,Par$j])
      bRed4[Par$j-2] = bRed4[Par$j-2] - A2red4[Par$j-2]*(Red4[i1,Par$j] - Par$KC4*Red4[i1,Par$j])

      uox = c(rep(0, DerApprox))
      vox = c(rep(1, DerApprox))
      uRed = c(rep(0, DerApprox))
      vRed = c(rep(1, DerApprox))
      uRed2 = c(rep(0, DerApprox))
      vRed2 = c(rep(1, DerApprox))
      uRed3 = c(rep(0, DerApprox))
      vRed3 = c(rep(1, DerApprox))
      uRed4 = c(rep(0, DerApprox))
      vRed4 = c(rep(1, DerApprox))

      for (j1 in ((Par$j-3):1)) {

        bOx[j1] = bOx[j1] - A2[Par$j-2]*bOx[j1+1]/A1[j1+1]
        bRed[j1] = bRed[j1] - A2red[Par$j-2]*bRed[j1+1]/A1red[j1+1]
        bRed2[j1] = bRed2[j1] - A2red2[Par$j-2]*bRed2[j1+1]/A1red2[j1+1]
        bRed3[j1] = bRed3[j1] - A2red3[Par$j-2]*bRed3[j1+1]/A1red3[j1+1]
        bRed4[j1] = bRed4[j1] - A2red4[Par$j-2]*bRed4[j1+1]/A1red4[j1+1]

        A1[j1] = A1[j1] - A2[Par$j-2]/A1[j1+1]
        A1red[j1] = A1red[j1] - A2red[Par$j-2]/A1red[j1+1]
        A1red2[j1] = A1red2[j1] - A2red2[Par$j-2]/A1red2[j1+1]
        A1red3[j1] = A1red3[j1] - A2red3[Par$j-2]/A1red3[j1+1]
        A1red4[j1] = A1red4[j1] - A2red4[Par$j-2]/A1red4[j1+1]

      }

      for (m in 2:DerApprox) {

        uox[m] = (bOx[m-1] - uox[m-1])/A1[m-1]
        vox[m] = -vox[m-1]/A1[m-1]
        uRed[m] = (bRed[m-1] - uRed[m-1])/A1red[m-1]
        vRed[m] = -vRed[m-1]/A1red[m-1]
        uRed2[m] = (bRed2[m-1] - uRed2[m-1])/A1red2[m-1]
        vRed2[m] = -vRed2[m-1]/A1red2[m-1]
        uRed3[m] = (bRed3[m-1] - uRed3[m-1])/A1red3[m-1]
        vRed3[m] = -vRed3[m-1]/A1red3[m-1]
        uRed4[m] = (bRed4[m-1] - uRed4[m-1])/A1red4[m-1]
        vRed4[m] = -vRed4[m-1]/A1red4[m-1]
      }

      B = matrix(data = c((Par$Kf1[i1]*Par$h - sum(vox*Derv(npoints = DerApprox, CoefMat = T))), -Par$Kb1[i1]*Par$h, 0, 0, 0,
                          -Par$Kf1[i1]*Par$h, Par$Kb1[i1]*Par$h + Par$Kf2[i1]*Par$h - sum(vRed*Derv(npoints = DerApprox, CoefMat = T)), -Par$Kb2[i1]*Par$h, 0, 0,
                          0, -Par$Kf2[i1]*Par$h, Par$Kb2[i1]*Par$h + Par$Kf3[i1]*Par$h - sum(vRed2*Derv(npoints = DerApprox, CoefMat = T)), -Par$Kb3[i1]*Par$h, 0,
                          0, 0, -Par$Kf3[i1]*Par$h, Par$Kb3[i1]*Par$h + Par$Kf4[i1]*Par$h - sum(vRed3*Derv(npoints = DerApprox, CoefMat = T)), -Par$Kb4[i1]*Par$h,
                          sum(vox*Derv(npoints = DerApprox, CoefMat = T)),
                          sum(vRed*Derv(npoints = DerApprox, CoefMat = T)),
                          sum(vRed2*Derv(npoints = DerApprox, CoefMat = T)),
                          sum(vRed3*Derv(npoints = DerApprox, CoefMat = T)),
                          sum(vRed4*Derv(npoints = DerApprox, CoefMat = T))),
                 byrow = T, nrow = 5, ncol = 5)
      Y = matrix(data = c(sum(uox*Derv(npoints = DerApprox, CoefMat = T)), sum(uRed*Derv(npoints = DerApprox, CoefMat = T)), sum(uRed2*Derv(npoints = DerApprox, CoefMat = T)), sum(uRed3*Derv(npoints = DerApprox, CoefMat = T)),
                          -sum(uox*Derv(npoints = DerApprox, CoefMat = T)) - sum(uRed*Derv(npoints = DerApprox, CoefMat = T)) - sum(uRed2*Derv(npoints = DerApprox, CoefMat = T)) - sum(uRed3*Derv(npoints = DerApprox, CoefMat = T)) - sum(uRed4*Derv(npoints = DerApprox, CoefMat = T))))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red1[i1+1,1] = C[2]
      Red2[i1+1,1] = C[3]
      Red3[i1+1,1] = C[4]
      Red4[i1+1,1] = C[5]

      for (j1 in 1:(Par$j-2)) {
        Ox[i1+1,j1+1] = (bOx[j1] - Ox[i1+1,j1])/A1[j1]
        Red1[i1+1,j1+1] = (bRed[j1] - Red1[i1+1,j1])/A1red[j1]
        Red2[i1+1,j1+1] = (bRed2[j1] - Red2[i1+1,j1])/A1red2[j1]
        Red3[i1+1,j1+1] = (bRed3[j1] - Red3[i1+1,j1])/A1red3[j1]
        Red4[i1+1,j1+1] = (bRed4[j1] - Red4[i1+1,j1])/A1red4[j1]
      }
    }

    Jox = Derv(Ox = Ox, h = Par$h, npoints = DerApprox)
    Jred1 = Derv(Ox = Red1, h = Par$h, npoints = DerApprox)
    Jred2 = Derv(Ox = Red2, h = Par$h, npoints = DerApprox)
    Jred3 = Derv(Ox = Red3, h = Par$h, npoints = DerApprox)
    Jred4 = Derv(Ox = Red4, h = Par$h, npoints = DerApprox)

  } else if (Method == "CN") {

    al1 = 1/(Par$h^2)
    al2 = -2/(Par$h^2)
    al3 = 1/(Par$h^2)
    a1 = (al2 - 2/Par$dtn)/al1
    a2 = al3/al1
    a3 = (al2 + 2/Par$dtn)/al1
    al1red = Par$DRED/(Par$h^2)
    al2red = -(2*Par$DRED)/(Par$h^2)
    al3red = Par$DRED/(Par$h^2)
    a1red = (al2red - 2/Par$dtn)/al1red
    a2red = al3red/al1red
    a3red = (al2red + 2/Par$dtn)/al1red
    al1red2 = Par$DRED2/(Par$h^2)
    al2red2 = -(2*Par$DRED2)/(Par$h^2)
    al3red2 = Par$DRED2/(Par$h^2)
    a1red2 = (al2red2 - 2/Par$dtn)/al1red2
    a2red2 = al3red2/al1red2
    a3red2 = (al2red2 + 2/Par$dtn)/al1red2
    al1red3 = Par$DRED3/(Par$h^2)
    al2red3 = -(2*Par$DRED3)/(Par$h^2)
    al3red3 = Par$DRED3/(Par$h^2)
    a1red3 = (al2red3 - 2/Par$dtn)/al1red3
    a2red3 = al3red3/al1red3
    a3red3 = (al2red3 + 2/Par$dtn)/al1red3
    al1red4 = Par$DRED4/(Par$h^2)
    al2red4 = -(2*Par$DRED4)/(Par$h^2)
    al3red4 = Par$DRED4/(Par$h^2)
    a1red4 = (al2red4 - 2/Par$dtn)/al1red4
    a2red4 = al3red4/al1red4
    a3red4 = (al2red4 + 2/Par$dtn)/al1red4

    for (i1 in 1:(Par$l-1)) {

      B = matrix(data = c((Par$Kf1[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1]), -Par$Kb1[i1]*Par$h, 0, 0, 0,
                          -Par$Kf1[i1]*Par$h, Par$Kb1[i1]*Par$h + Par$Kf2[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1], -Par$Kb2[i1]*Par$h, 0, 0,
                          0, -Par$Kf2[i1]*Par$h, Par$Kb2[i1]*Par$h + Par$Kf3[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1], -Par$Kb3[i1]*Par$h, 0,
                          0, 0, -Par$Kf3[i1]*Par$h, Par$Kb3[i1]*Par$h + Par$Kf4[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1], -Par$Kb4[i1]*Par$h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 5, ncol = 5)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]),
                          Par$DRED*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red1[i1,2:DerApprox]),
                          Par$DRED2*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red2[i1,2:DerApprox]),
                          Par$DRED3*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red3[i1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - Par$DRED*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red1[i1,2:DerApprox]) - Par$DRED2*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red2[i1,2:DerApprox]) - Par$DRED3*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red3[i1,2:DerApprox]) - Par$DRED4*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red4[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1] - Par$KCo*Ox[i1,1]
      Red1[i1,1] = C[2] - Par$KC1*Red1[i1,1]
      Red2[i1,1] = C[3] - Par$KC2*Red2[i1,1]
      Red3[i1,1] = C[4] - Par$KC3*Red3[i1,1]
      Red4[i1,1] = C[5] - Par$KC4*Red4[i1,1]

      bOx = -a3*Ox[i1,(2:(Par$j-1))] - Ox[i1,(1:(Par$j-2))] - a2*Ox[i1,(3:Par$j)] + Par$KCo*Ox[i1,2:(Par$j-1)]/al1
      bRed = -a3red*Red1[i1,(2:(Par$j-1))]- Red1[i1,(1:(Par$j-2))] - a2red*Red1[i1,(3:Par$j)] + Par$KC1*Red1[i1,2:(Par$j-1)]/al1red
      bRed2 = -a3red2*Red2[i1,(2:(Par$j-1))]- Red2[i1,(1:(Par$j-2))] - a2red2*Red2[i1,(3:Par$j)] + Par$KC2*Red2[i1,2:(Par$j-1)]/al1red2
      bRed3 = -a3red3*Red3[i1,(2:(Par$j-1))]- Red3[i1,(1:(Par$j-2))] - a2red3*Red3[i1,(3:Par$j)] + Par$KC3*Red3[i1,2:(Par$j-1)]/al1red3
      bRed4 = -a3red4*Red4[i1,(2:(Par$j-1))]- Red4[i1,(1:(Par$j-2))] - a2red4*Red4[i1,(3:Par$j)] + Par$KC4*Red4[i1,2:(Par$j-1)]/al1red4
      A = c(rep(1,Par$j-2))
      A1 = c(rep(a1,Par$j-2))
      A2 = c(rep(a2,Par$j-2))
      A1red = c(rep(a1red,Par$j-2))
      A2red = c(rep(a2red,Par$j-2))
      A1red2 = c(rep(a1red2,Par$j-2))
      A2red2 = c(rep(a2red2,Par$j-2))
      A1red3 = c(rep(a1red3,Par$j-2))
      A2red3 = c(rep(a2red3,Par$j-2))
      A1red4 = c(rep(a1red4,Par$j-2))
      A2red4 = c(rep(a2red4,Par$j-2))

      bOx[Par$j-2] = bOx[Par$j-2] - A2[Par$j-2]*(1 - Par$KCo*Ox[i1,Par$j])
      bRed[Par$j-2] = bRed[Par$j-2] - A2red[Par$j-2]*(Red1[i1,Par$j] - Par$KC1*Red1[i1,Par$j])
      bRed2[Par$j-2] = bRed2[Par$j-2] - A2red2[Par$j-2]*(Red2[i1,Par$j] - Par$KC2*Red2[i1,Par$j])
      bRed3[Par$j-2] = bRed3[Par$j-2] - A2red3[Par$j-2]*(Red3[i1,Par$j] - Par$KC3*Red3[i1,Par$j])
      bRed4[Par$j-2] = bRed4[Par$j-2] - A2red4[Par$j-2]*(Red4[i1,Par$j] - Par$KC4*Red4[i1,Par$j])
      uox = c(rep(0, DerApprox))
      vox = c(rep(1, DerApprox))
      uRed = c(rep(0, DerApprox))
      vRed = c(rep(1, DerApprox))
      uRed2 = c(rep(0, DerApprox))
      vRed2 = c(rep(1, DerApprox))
      uRed3 = c(rep(0, DerApprox))
      vRed3 = c(rep(1, DerApprox))
      uRed4 = c(rep(0, DerApprox))
      vRed4 = c(rep(1, DerApprox))

      for (j1 in ((Par$j-3):1)) {

        bOx[j1] = bOx[j1] - A2[Par$j-2]*bOx[j1+1]/A1[j1+1]
        bRed[j1] = bRed[j1] - A2red[Par$j-2]*bRed[j1+1]/A1red[j1+1]
        bRed2[j1] = bRed2[j1] - A2red2[Par$j-2]*bRed2[j1+1]/A1red2[j1+1]
        bRed3[j1] = bRed3[j1] - A2red3[Par$j-2]*bRed3[j1+1]/A1red3[j1+1]
        bRed4[j1] = bRed4[j1] - A2red4[Par$j-2]*bRed4[j1+1]/A1red4[j1+1]

        A1[j1] = A1[j1] - A2[Par$j-2]/A1[j1+1]
        A1red[j1] = A1red[j1] - A2red[Par$j-2]/A1red[j1+1]
        A1red2[j1] = A1red2[j1] - A2red2[Par$j-2]/A1red2[j1+1]
        A1red3[j1] = A1red3[j1] - A2red3[Par$j-2]/A1red3[j1+1]
        A1red4[j1] = A1red4[j1] - A2red4[Par$j-2]/A1red4[j1+1]

      }

      for (m in 2:DerApprox) {

        uox[m] = (bOx[m-1] - uox[m-1])/A1[m-1]
        vox[m] = -vox[m-1]/A1[m-1]
        uRed[m] = (bRed[m-1] - uRed[m-1])/A1red[m-1]
        vRed[m] = -vRed[m-1]/A1red[m-1]
        uRed2[m] = (bRed2[m-1] - uRed2[m-1])/A1red2[m-1]
        vRed2[m] = -vRed2[m-1]/A1red2[m-1]
        uRed3[m] = (bRed3[m-1] - uRed3[m-1])/A1red3[m-1]
        vRed3[m] = -vRed3[m-1]/A1red3[m-1]
        uRed4[m] = (bRed4[m-1] - uRed4[m-1])/A1red4[m-1]
        vRed4[m] = -vRed4[m-1]/A1red4[m-1]

      }

      B = matrix(data = c((Par$Kf1[i1]*Par$h - sum(vox*Derv(npoints = DerApprox, CoefMat = T))), -Par$Kb1[i1]*Par$h, 0, 0, 0,
                          -Par$Kf1[i1]*Par$h, Par$Kb1[i1]*Par$h + Par$Kf2[i1]*Par$h - sum(vRed*Derv(npoints = DerApprox, CoefMat = T)), -Par$Kb2[i1]*Par$h, 0, 0,
                          0, -Par$Kf2[i1]*Par$h, Par$Kb2[i1]*Par$h + Par$Kf3[i1]*Par$h - sum(vRed2*Derv(npoints = DerApprox, CoefMat = T)), -Par$Kb3[i1]*Par$h, 0,
                          0, 0, -Par$Kf3[i1]*Par$h, Par$Kb3[i1]*Par$h + Par$Kf4[i1]*Par$h - sum(vRed3*Derv(npoints = DerApprox, CoefMat = T)), -Par$Kb4[i1]*Par$h,
                          sum(vox*Derv(npoints = DerApprox, CoefMat = T)),
                          sum(vRed*Derv(npoints = DerApprox, CoefMat = T)),
                          sum(vRed2*Derv(npoints = DerApprox, CoefMat = T)),
                          sum(vRed3*Derv(npoints = DerApprox, CoefMat = T)),
                          sum(vRed4*Derv(npoints = DerApprox, CoefMat = T))),
                 byrow = T, nrow = 5, ncol = 5)
      Y = matrix(data = c(sum(uox*Derv(npoints = DerApprox, CoefMat = T)), sum(uRed*Derv(npoints = DerApprox, CoefMat = T)), sum(uRed2*Derv(npoints = DerApprox, CoefMat = T)), sum(uRed3*Derv(npoints = DerApprox, CoefMat = T)),
                          -sum(uox*Derv(npoints = DerApprox, CoefMat = T)) - sum(uRed*Derv(npoints = DerApprox, CoefMat = T)) - sum(uRed2*Derv(npoints = DerApprox, CoefMat = T)) - sum(uRed3*Derv(npoints = DerApprox, CoefMat = T)) - sum(uRed4*Derv(npoints = DerApprox, CoefMat = T))))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red1[i1+1,1] = C[2]
      Red2[i1+1,1] = C[3]
      Red3[i1+1,1] = C[4]
      Red4[i1+1,1] = C[5]


      for (j1 in 1:(Par$j-2)) {
        Ox[i1+1,j1+1] = (bOx[j1] - Ox[i1+1,j1])/A1[j1]
        Red1[i1+1,j1+1] = (bRed[j1] - Red1[i1+1,j1])/A1red[j1]
        Red2[i1+1,j1+1] = (bRed2[j1] - Red2[i1+1,j1])/A1red2[j1]
        Red3[i1+1,j1+1] = (bRed3[j1] - Red3[i1+1,j1])/A1red3[j1]
        Red4[i1+1,j1+1] = (bRed4[j1] - Red4[i1+1,j1])/A1red4[j1]
      }
    }

    Jox = Derv(Ox = Ox, h = Par$h, npoints = DerApprox)
    Jred1 = Derv(Ox = Red1, h = Par$h, npoints = DerApprox)
    Jred2 = Derv(Ox = Red2, h = Par$h, npoints = DerApprox)
    Jred3 = Derv(Ox = Red3, h = Par$h, npoints = DerApprox)
    Jred4 = Derv(Ox = Red4, h = Par$h, npoints = DerApprox)

  } else if (Method == "BDF") {

    al1 = 1/(Par$h^2)
    al2 = -2/(Par$h^2)
    al3 = 1/(Par$h^2)
    a1 = (al2 - 1.5/Par$dtn)/al1
    a2 = al3/al1
    al1red = Par$DRED/(Par$h^2)
    al2red = -(2*Par$DRED)/(Par$h^2)
    al3red = Par$DRED/(Par$h^2)
    a1red = (al2red - 1.5/Par$dtn)/al1red
    a2red = al3red/al1red
    al1red2 = Par$DRED2/(Par$h^2)
    al2red2 = -(2*Par$DRED2)/(Par$h^2)
    al3red2 = Par$DRED2/(Par$h^2)
    a1red2 = (al2red2 - 1.5/Par$dtn)/al1red2
    a2red2 = al3red2/al1red2
    al1red3 = Par$DRED3/(Par$h^2)
    al2red3 = -(2*Par$DRED3)/(Par$h^2)
    al3red3 = Par$DRED3/(Par$h^2)
    a1red3 = (al2red3 - 1.5/Par$dtn)/al1red3
    a2red3 = al3red3/al1red3
    al1red4 = Par$DRED4/(Par$h^2)
    al2red4 = -(2*Par$DRED4)/(Par$h^2)
    al3red4 = Par$DRED4/(Par$h^2)
    a1red4 = (al2red4 - 1.5/Par$dtn)/al1red4
    a2red4 = al3red4/al1red4


    for (i1 in 1:(Par$l-1)) {

      B = matrix(data = c((Par$Kf1[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1]), -Par$Kb1[i1]*Par$h, 0, 0, 0,
                          -Par$Kf1[i1]*Par$h, Par$Kb1[i1]*Par$h + Par$Kf2[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1], -Par$Kb2[i1]*Par$h, 0, 0,
                          0, -Par$Kf2[i1]*Par$h, Par$Kb2[i1]*Par$h + Par$Kf3[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1], -Par$Kb3[i1]*Par$h, 0,
                          0, 0, -Par$Kf3[i1]*Par$h, Par$Kb3[i1]*Par$h + Par$Kf4[i1]*Par$h - Derv(npoints = DerApprox, CoefMat = T)[1], -Par$Kb4[i1]*Par$h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 5, ncol = 5)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]),
                          Par$DRED*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red1[i1,2:DerApprox]),
                          Par$DRED2*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red2[i1,2:DerApprox]),
                          Par$DRED3*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red3[i1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - Par$DRED*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red1[i1,2:DerApprox]) - Par$DRED2*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red2[i1,2:DerApprox]) - Par$DRED3*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red3[i1,2:DerApprox]) - Par$DRED4*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red4[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1] - Par$KCo*Ox[i1,1]
      Red1[i1,1] = C[2] - Par$KC1*Red1[i1,1]
      Red2[i1,1] = C[3] - Par$KC2*Red2[i1,1]
      Red3[i1,1] = C[4] - Par$KC3*Red3[i1,1]
      Red4[i1,1] = C[5] - Par$KC4*Red4[i1,1]

      if (i1 == 1) {
        bOx = -2*Ox[i1,2:(Par$j-1)]/(Par$dtn*al1) + Ox[1,2:(Par$j-1)]/(2*Par$dtn*al1) + Par$KCo*Ox[i1,2:(Par$j-1)]/al1
        bRed = -2*Red1[i1,2:(Par$j-1)]/(Par$dtn*al1red) + Red1[1,2:(Par$j-1)]/(2*Par$dtn*al1red) + Par$KC1*Red1[i1,2:(Par$j-1)]/al1red
        bRed2 = -2*Red2[i1,2:(Par$j-1)]/(Par$dtn*al1red2) + Red2[1,2:(Par$j-1)]/(2*Par$dtn*al1red2) + Par$KC2*Red2[i1,2:(Par$j-1)]/al1red2
        bRed3 = -2*Red3[i1,2:(Par$j-1)]/(Par$dtn*al1red3) + Red3[1,2:(Par$j-1)]/(2*Par$dtn*al1red3) + Par$KC3*Red3[i1,2:(Par$j-1)]/al1red3
        bRed4 = -2*Red4[i1,2:(Par$j-1)]/(Par$dtn*al1red4) + Red4[1,2:(Par$j-1)]/(2*Par$dtn*al1red4) + Par$KC4*Red4[i1,2:(Par$j-1)]/al1red4

      } else {
        bOx = -2*Ox[i1,2:(Par$j-1)]/(Par$dtn*al1) + Ox[i1-1,2:(Par$j-1)]/(2*Par$dtn*al1) + Par$KCo*Ox[i1,2:(Par$j-1)]/al1
        bRed = -2*Red1[i1,2:(Par$j-1)]/(Par$dtn*al1red) + Red1[i1-1,2:(Par$j-1)]/(2*Par$dtn*al1red) + Par$KC1*Red1[i1,2:(Par$j-1)]/al1red
        bRed2 = -2*Red2[i1,2:(Par$j-1)]/(Par$dtn*al1red2) + Red2[i1-1,2:(Par$j-1)]/(2*Par$dtn*al1red2) + Par$KC2*Red2[i1,2:(Par$j-1)]/al1red2
        bRed3 = -2*Red3[i1,2:(Par$j-1)]/(Par$dtn*al1red3) + Red3[i1-1,2:(Par$j-1)]/(2*Par$dtn*al1red3) + Par$KC3*Red3[i1,2:(Par$j-1)]/al1red3
        bRed4 = -2*Red4[i1,2:(Par$j-1)]/(Par$dtn*al1red4) + Red4[i1-1,2:(Par$j-1)]/(2*Par$dtn*al1red4) + Par$KC4*Red4[i1,2:(Par$j-1)]/al1red4
      }

      A = c(rep(1,Par$j-2))
      A1 = c(rep(a1,Par$j-2))
      A2 = c(rep(a2,Par$j-2))
      A1red = c(rep(a1red,Par$j-2))
      A2red = c(rep(a2red,Par$j-2))
      A1red2 = c(rep(a1red2,Par$j-2))
      A2red2 = c(rep(a2red2,Par$j-2))
      A1red3 = c(rep(a1red3,Par$j-2))
      A2red3 = c(rep(a2red3,Par$j-2))
      A1red4 = c(rep(a1red4,Par$j-2))
      A2red4 = c(rep(a2red4,Par$j-2))

      bOx[Par$j-2] = bOx[Par$j-2] - A2[Par$j-2]*(1 - Par$KCo*Ox[i1,Par$j])
      bRed[Par$j-2] = bRed[Par$j-2] - A2red[Par$j-2]*(Red1[i1,Par$j] - Par$KC1*Red1[i1,Par$j])
      bRed2[Par$j-2] = bRed2[Par$j-2] - A2red2[Par$j-2]*(Red2[i1,Par$j] - Par$KC1*Red1[i1,Par$j])
      bRed3[Par$j-2] = bRed3[Par$j-2] - A2red3[Par$j-2]*(Red3[i1,Par$j] - Par$KC2*Red2[i1,Par$j])
      bRed4[Par$j-2] = bRed4[Par$j-2] - A2red4[Par$j-2]*(Red4[i1,Par$j] - Par$KC3*Red3[i1,Par$j])

      uox = c(rep(0, DerApprox))
      vox = c(rep(1, DerApprox))
      uRed = c(rep(0, DerApprox))
      vRed = c(rep(1, DerApprox))
      uRed2 = c(rep(0, DerApprox))
      vRed2 = c(rep(1, DerApprox))
      uRed3 = c(rep(0, DerApprox))
      vRed3 = c(rep(1, DerApprox))
      uRed4 = c(rep(0, DerApprox))
      vRed4 = c(rep(1, DerApprox))

      for (j1 in ((Par$j-3):1)) {

        bOx[j1] = bOx[j1] - A2[Par$j-2]*bOx[j1+1]/A1[j1+1]
        bRed[j1] = bRed[j1] - A2red[Par$j-2]*bRed[j1+1]/A1red[j1+1]
        bRed2[j1] = bRed2[j1] - A2red2[Par$j-2]*bRed2[j1+1]/A1red2[j1+1]
        bRed3[j1] = bRed3[j1] - A2red3[Par$j-2]*bRed3[j1+1]/A1red3[j1+1]
        bRed4[j1] = bRed4[j1] - A2red4[Par$j-2]*bRed4[j1+1]/A1red4[j1+1]
        A1[j1] = A1[j1] - A2[Par$j-2]/A1[j1+1]
        A1red[j1] = A1red[j1] - A2red[Par$j-2]/A1red[j1+1]
        A1red2[j1] = A1red2[j1] - A2red2[Par$j-2]/A1red2[j1+1]
        A1red3[j1] = A1red3[j1] - A2red3[Par$j-2]/A1red3[j1+1]
        A1red4[j1] = A1red4[j1] - A2red4[Par$j-2]/A1red4[j1+1]

      }

      for (m in 2:DerApprox) {

        uox[m] = (bOx[m-1] - uox[m-1])/A1[m-1]
        vox[m] = -vox[m-1]/A1[m-1]
        uRed[m] = (bRed[m-1] - uRed[m-1])/A1red[m-1]
        vRed[m] = -vRed[m-1]/A1red[m-1]
        uRed2[m] = (bRed2[m-1] - uRed2[m-1])/A1red2[m-1]
        vRed2[m] = -vRed2[m-1]/A1red2[m-1]
        uRed3[m] = (bRed3[m-1] - uRed3[m-1])/A1red3[m-1]
        vRed3[m] = -vRed3[m-1]/A1red3[m-1]
        uRed4[m] = (bRed4[m-1] - uRed4[m-1])/A1red4[m-1]
        vRed4[m] = -vRed4[m-1]/A1red4[m-1]

      }

      B = matrix(data = c((Par$Kf1[i1]*Par$h - sum(vox*Derv(npoints = DerApprox, CoefMat = T))), -Par$Kb1[i1]*Par$h, 0, 0, 0,
                          -Par$Kf1[i1]*Par$h, Par$Kb1[i1]*Par$h + Par$Kf2[i1]*Par$h - sum(vRed*Derv(npoints = DerApprox, CoefMat = T)), -Par$Kb2[i1]*Par$h, 0, 0,
                          0, -Par$Kf2[i1]*Par$h, Par$Kb2[i1]*Par$h + Par$Kf3[i1]*Par$h - sum(vRed2*Derv(npoints = DerApprox, CoefMat = T)), -Par$Kb3[i1]*Par$h, 0,
                          0, 0, -Par$Kf3[i1]*Par$h, Par$Kb3[i1]*Par$h + Par$Kf4[i1]*Par$h - sum(vRed3*Derv(npoints = DerApprox, CoefMat = T)), -Par$Kb4[i1]*Par$h,
                          sum(vox*Derv(npoints = DerApprox, CoefMat = T)),
                          sum(vRed*Derv(npoints = DerApprox, CoefMat = T)),
                          sum(vRed2*Derv(npoints = DerApprox, CoefMat = T)),
                          sum(vRed3*Derv(npoints = DerApprox, CoefMat = T)),
                          sum(vRed4*Derv(npoints = DerApprox, CoefMat = T))),
                 byrow = T, nrow = 5, ncol = 5)
      Y = matrix(data = c(sum(uox*Derv(npoints = DerApprox, CoefMat = T)), sum(uRed*Derv(npoints = DerApprox, CoefMat = T)), sum(uRed2*Derv(npoints = DerApprox, CoefMat = T)), sum(uRed3*Derv(npoints = DerApprox, CoefMat = T)),
                          -sum(uox*Derv(npoints = DerApprox, CoefMat = T)) - sum(uRed*Derv(npoints = DerApprox, CoefMat = T)) - sum(uRed2*Derv(npoints = DerApprox, CoefMat = T)) - sum(uRed3*Derv(npoints = DerApprox, CoefMat = T)) - sum(uRed4*Derv(npoints = DerApprox, CoefMat = T))))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red1[i1+1,1] = C[2]
      Red2[i1+1,1] = C[3]
      Red3[i1+1,1] = C[4]
      Red4[i1+1,1] = C[5]


      for (j1 in 1:(Par$j-2)) {
        Ox[i1+1,j1+1] = (bOx[j1] - Ox[i1+1,j1])/A1[j1]
        Red1[i1+1,j1+1] = (bRed[j1] - Red1[i1+1,j1])/A1red[j1]
        Red2[i1+1,j1+1] = (bRed2[j1] - Red2[i1+1,j1])/A1red2[j1]
        Red3[i1+1,j1+1] = (bRed3[j1] - Red3[i1+1,j1])/A1red3[j1]
        Red4[i1+1,j1+1] = (bRed4[j1] - Red4[i1+1,j1])/A1red4[j1]
      }
    }

    Jox = Derv(Ox = Ox, h = Par$h, npoints = DerApprox)
    Jred1 = Derv(Ox = Red1, h = Par$h, npoints = DerApprox)
    Jred2 = Derv(Ox = Red2, h = Par$h, npoints = DerApprox)
    Jred3 = Derv(Ox = Red3, h = Par$h, npoints = DerApprox)
    Jred4 = Derv(Ox = Red4, h = Par$h, npoints = DerApprox)

  } else if (!(Method %in% c("Euler", "BI", "RK4", "CN", "BDF"))) {
    return("Available methods are Euler, BI, RK4, CN and BDF")
  }

  G1 = Jox
  G2 = Jox + Jred1
  G3 = Jox + Jred1 + Jred2
  G4 = Jox + Jred1 + Jred2 + Jred3
  G5 = Jox + Jred1 + Jred2 + Jred3 + Jred4

  i = (n*Par$FA*(G1+G2+G3+G4+G5)*Dx1*Area*Co)/(sqrt(Dx1*Par$tau))

  graphy = ggplot(data = data.frame(i[1:(length(i)-1)],Par$PotentialScan[1:(length(i)-1)]),
                  aes(y = i[1:(length(i)-1)], x = Par$PotentialScan[1:(length(i)-1)])) +
    geom_point() + scale_x_continuous(trans = "reverse") +
    xlab("E / V") +
    ylab("I / A") +
    theme_classic()

  if (errCheck == TRUE){
    return(list((G1+G2),Dx1,Dred,Dred2,Co,
                Par$dtn,Par$h,i,Par$l,Par$j,n,
                Area,Par$DOx,Par$DRED,Par$DRED2,
                Par$p1,Par$p2,Par$p3,Par$p4))
  } else {
    return(graphy)
  }

}
