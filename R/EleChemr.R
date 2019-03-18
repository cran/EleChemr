#' Chrono amperometry digital simulation
#'
#' Return a graph I vs t of the electrochemical process
#'
#' @param Co bulk concentration
#' @param exptime experimental time to be simulated
#' @param Dx diffusion coefficient
#' @param Temp temperature in kelvin
#' @param n number of electrons involved in the process
#' @param Area area of the electrode
#' @param DerApprox number of point for the approximation of the first derivative
#' @param errCheck if true the function returns a list with parameters for CottrCheck function
#' @param Method method to be used for the simulation = "Euler" "BI" "RK4"
#'
#'
#' @return if errCheck == F a graph I vs t, if errCheck == T a list
#'
#' @examples
#' ChronAmp(Co = 0.001, exptime = 1, DerApprox = 2, errCheck = FALSE, Method = "Euler")
#'
#' @export
#' @import ggplot2
#' @importFrom matlab zeros ones
#' @importFrom pracma inv

ChronAmp = function(Co = 0.001, exptime = 1, Dx = 0.00001,
                    Temp = 298.15, n = 1, Area = 1, DerApprox = 2,
                    errCheck = FALSE, Method = "Euler") {

  FA = 96485
  R = 8.3145
  f = ((FA*n)/(R*Temp))
  dt = 0.01
  Da = Dx/Dx
  Dm = 0.45
  l = exptime/dt
  h = sqrt((Da*dt)/Dm)
  j = ceiling(6*(l)^0.5)
  vt = c(1:l)
  t = dt*vt
  Ox = ones(l, j)
  Jox = zeros(l, 1)

  DerivMatrix = matrix(data = c(-1,1,0,0,0,0,
                                -3,4,-1,0,0,0,
                                -11,18,-9,2,0,0,
                                -25,48,-36,16,-3,0,
                                -137,300,-300,200,-75,12),
                       byrow = T, nrow = 5, ncol = 6)

  NormalFact = matrix(data = c(1,2,6,12,60), nrow = 5)

  if (Method == "Euler") {

    for (i1 in 1:(l-1)) {

      Ox[i1,1] = 0

      for (j1 in 2:(j-1)) {
        Ox[i1 + 1,j1] = Ox[i1,j1] + Dm*(Ox[i1, j1 -1] - 2*Ox[i1, j1] + Ox[i1, j1+1])
      }
      Jox[i1,1] = sum(DerivMatrix[DerApprox -1,1:DerApprox]*Ox[i1,1:DerApprox])/(h*NormalFact[DerApprox-1])
    }
  } else if (Method == "RK4") {

    for (i1 in 1:(l-1)) {
      k1 = zeros(j)
      k2 = zeros(j)
      k3 = zeros(j)
      k4 = zeros(j)
      Ox[i1,1] = 0
      Ox[i1 +1, 1] = 0

      for (j1 in 2:(j-1)) {
        k1[j1] = Dm*(Ox[i1, j1 -1] - 2*Ox[i1, j1] + Ox[i1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + k1[j1]*0.5
      }
      for (j1 in 2:(j-1)) {
        k2[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + k2[j1]*0.5
      }
      for (j1 in 2:(j-1)) {
        k3[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + k3[j1]
      }
      for (j1 in 2:(j-1)) {
        k4[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + (k1[j1] + 2*k2[j1] + 2*k3[j1] + k4[j1])/6
      }
      Jox[i1,1] = sum(DerivMatrix[DerApprox -1,1:DerApprox]*Ox[i1,1:DerApprox])/(h*NormalFact[DerApprox-1])
    }
  } else if (Method == "BI") {
    al1 = 1/(h^2)
    al2 = -2/(h^2)
    al3 = 1/(h^2)
    a1 = (al2 - 1/dt)/al1
    a2 = al3/al1

    for (i1 in 1:(l-1)) {

      Ox[i1,1] = 0
      Ox[i1+1,1] = 0
      b = (-Ox[i1,]/(al1*dt))
      A = c(rep(1,j))
      A1 = c(rep(a1,j))
      A2 = c(rep(a2,j))
      b[j] = b[j] - A2[j]*1

      for (j1 in ((j-1):2)) {

        b[j1] = b[j1] - A2[j]*b[j1+1]/A1[j1+1]
        A1[j1] = A1[j1] - A2[j]/A1[j1+1]

      }

      for (j1 in 2:j) {
        Ox[i1+1,j1] = (b[j1] -Ox[i1+1,j1-1])/A1[j1]
      }
      Jox[i1,1] = sum(DerivMatrix[DerApprox -1,1:DerApprox]*Ox[i1,1:DerApprox])/(h*NormalFact[DerApprox-1])
    }
  } else if (!(Method %in% c("Euler", "BI", "RK4"))) {
      return("Available methods are Euler, BI and RK4")
  }

  G = Jox
  i = (n*FA*G*Dx*Area*Co)/(6.4*(h*(Dx^0.5)))

  graphy = ggplot(data = data.frame(i[1:(length(i)-1)],t[1:(length(i)-1)]),
                  aes(y = i[1:(length(i)-1)], x = t[1:(length(i)-1)])) +
    geom_point() + xlab("Time(s)") +
    ylab("Current (A)") + theme_classic()

  if (errCheck == TRUE){
    return(list(G,Dx,Co,dt,h,l,n,Area))
  } else {
    return(graphy)
  }
}


#' Chrono amperometry with a finite step digital simulation
#'
#' Return a graph I vs t of the electrochemical process
#'
#' @param Co bulk concentration
#' @param exptime experimental time to be simulated
#' @param Dx diffusion coefficient
#' @param eta overpotential of the step
#' @param Temp temperature in kelvin
#' @param n number of electrons involved in the process
#' @param Area area of the electrode
#' @param DerApprox number of point for the approximation of the first derivative
#' @param errCheck if true the function returns a list with parameters for CottrCheck function
#' @param Method method to be used for the simulation = "Euler" "BI" "RK4"
#'
#'
#' @return if errCheck == F a graph I vs t, if errCheck == T a list
#'
#' @examples
#' PotStep(Co = 0.001, exptime = 1, DerApprox = 2, errCheck = FALSE, Method = "Euler")
#'
#' @export
#' @import ggplot2
#' @importFrom matlab zeros ones
#' @importFrom pracma inv
#'

PotStep = function(Co = 0.001, exptime = 1, Dx = 0.00001,
                   eta = 0.1, Temp = 298.15, n = 1, Area = 1,
                   DerApprox = 2, errCheck = FALSE, Method = "Euler") {

  FA = 96485
  R = 8.3145
  f = ((FA*n)/(R*Temp))
  dt = 0.01
  p = eta*f
  Da = Dx/Dx
  Dm = 0.45
  l = exptime/dt
  h = sqrt((Da*dt)/Dm)
  j = ceiling(6*(Da*l)^0.5)
  vt = c(1:l)
  t = dt*vt
  Ox = ones(l, j)
  Red = zeros(l, j)
  Jox = zeros(l, 1)

  DerivMatrix = matrix(data = c(-1,1,0,0,0,0,
                                -3,4,-1,0,0,0,
                                -11,18,-9,2,0,0,
                                -25,48,-36,16,-3,0,
                                -137,300,-300,200,-75,12),
                       byrow = T, nrow = 5, ncol = 6)

  NormalFact = matrix(data = c(1,2,6,12,60), nrow = 5)
  if (Method == "Euler") {
    for (i1 in 1:(l-1)) {
      B = matrix(data = c(1,-exp(p),DerivMatrix[DerApprox-1,1],DerivMatrix[DerApprox-1,1]), byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(0, -sum(DerivMatrix[DerApprox -1,2:DerApprox]*Ox[i1,2:DerApprox]) - sum(DerivMatrix[DerApprox -1,2:DerApprox]*Red[i1,2:DerApprox])))
      C = inv(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2]
      for (j1 in 2:(j-1)) {
        Ox[i1+1,j1] = Ox[i1,j1] + Dm*(Ox[i1, j1-1] - 2*Ox[i1, j1] + Ox[i1, j1+1])
        Red[i1+1, j1] = Red[i1, j1] + Dm*(Red[i1, j1-1] + Red[i1, j1+1] - 2*Red[i1,j1])
      }
      Jox[i1,1] = sum(DerivMatrix[DerApprox -1,1:DerApprox]*Ox[i1,1:DerApprox])/(h*NormalFact[DerApprox-1])
    }
  } else if (Method == "RK4") {

    for (i1 in 1:(l-1)) {
      k1 = zeros(j)
      k2 = zeros(j)
      k3 = zeros(j)
      k4 = zeros(j)
      k1red = zeros(j)
      k2red = zeros(j)
      k3red = zeros(j)
      k4red = zeros(j)
      B = matrix(data = c(1,-exp(p),DerivMatrix[DerApprox-1,1],DerivMatrix[DerApprox-1,1]), byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(0, -sum(DerivMatrix[DerApprox -1,2:DerApprox]*Ox[i1,2:DerApprox]) - sum(DerivMatrix[DerApprox -1,2:DerApprox]*Red[i1,2:DerApprox])))
      C = inv(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2]


      B = matrix(data = c(1,-exp(p),DerivMatrix[DerApprox-1,1],DerivMatrix[DerApprox-1,1]), byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(0, -sum(DerivMatrix[DerApprox -1,2:DerApprox]*Ox[i1+1,2:DerApprox]) - sum(DerivMatrix[DerApprox -1,2:DerApprox]*Red[i1+1,2:DerApprox])))
      C = inv(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]

      for (j1 in 2:(j-1)) {
        k1[j1] = Dm*(Ox[i1, j1 -1] - 2*Ox[i1, j1] + Ox[i1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + k1[j1]*0.5
        k1red[j1] = Dm*(Red[i1, j1 -1] - 2*Red[i1, j1] + Red[i1, j1+1])
        Red[i1 + 1,j1] = Red[i1,j1] + k1red[j1]*0.5
      }
      for (j1 in 2:(j-1)) {
        k2[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + k2[j1]*0.5
        k2red[j1] = Dm*(Red[i1 + 1, j1 -1] - 2*Red[i1 + 1, j1] + Red[i1 + 1, j1+1])
        Red[i1 + 1,j1] = Red[i1,j1] + k2red[j1]*0.5
      }
      for (j1 in 2:(j-1)) {
        k3[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + k3[j1]
        k3red[j1] = Dm*(Red[i1 + 1, j1 -1] - 2*Red[i1 + 1, j1] + Red[i1 + 1, j1+1])
        Red[i1 + 1,j1] = Red[i1,j1] + k3red[j1]
      }
      for (j1 in 2:(j-1)) {
        k4[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + (k1[j1] + 2*k2[j1] + 2*k3[j1] + k4[j1])/6
        k4red[j1] = Dm*(Red[i1 + 1, j1 -1] - 2*Red[i1 + 1, j1] + Red[i1 + 1, j1+1])
        Red[i1 + 1,j1] = Red[i1,j1] + (k1red[j1] + 2*k2red[j1] + 2*k3red[j1] + k4red[j1])/6
      }
      Jox[i1,1] = sum(DerivMatrix[DerApprox -1,1:DerApprox]*Ox[i1,1:DerApprox])/(h*NormalFact[DerApprox-1])
    }
  }  else if (Method == "BI") {
    al1 = 1/(h^2)
    al2 = -2/(h^2)
    al3 = 1/(h^2)
    a1 = (al2 - 1/dt)/al1
    a2 = al3/al1

    for (i1 in 1:(l-1)) {

      B = matrix(data = c(1,-exp(p),DerivMatrix[DerApprox-1,1],DerivMatrix[DerApprox-1,1]), byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(0, -sum(DerivMatrix[DerApprox -1,2:DerApprox]*Ox[i1,2:DerApprox]) - sum(DerivMatrix[DerApprox -1,2:DerApprox]*Red[i1,2:DerApprox])))
      C = inv(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2]


      B = matrix(data = c(1,-exp(p),DerivMatrix[DerApprox-1,1],DerivMatrix[DerApprox-1,1]), byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(0, -sum(DerivMatrix[DerApprox -1,2:DerApprox]*Ox[i1+1,2:DerApprox]) - sum(DerivMatrix[DerApprox -1,2:DerApprox]*Red[i1+1,2:DerApprox])))
      C = inv(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]

      bOx = (-Ox[i1,]/(al1*dt))
      bRed = (-Red[i1,]/(al1*dt))
      A = c(rep(1,j))
      A1 = c(rep(a1,j))
      A2 = c(rep(a2,j))
      bOx[j] = bOx[j] - A2[j]*1
      bRed[j] = bRed[j] - A2[j]*1

      for (j1 in ((j-1):2)) {

        bOx[j1] = bOx[j1] - A2[j]*bOx[j1+1]/A1[j1+1]
        bRed[j1] = bRed[j1] - A2[j]*bRed[j1+1]/A1[j1+1]
        A1[j1] = A1[j1] - A2[j]/A1[j1+1]

      }

      for (j1 in 2:j) {
        Ox[i1+1,j1] = (bOx[j1] -Ox[i1+1,j1-1])/A1[j1]
        Red[i1+1,j1] = (bRed[j1] -Red[i1+1,j1-1])/A1[j1]
      }
      Jox[i1,1] = sum(DerivMatrix[DerApprox -1,1:DerApprox]*Ox[i1,1:DerApprox])/(h*NormalFact[DerApprox-1])
    }

  } else if (!(Method %in% c("Euler", "BI", "RK4"))) {
    return("Available methods are Euler, BI and RK4")
  }

  G = Jox
  i = (n*FA*G*Dx*Area*Co)/(6.4*h*sqrt(Dx))

  graphy = ggplot(data = data.frame(i[1:(length(i)-1)],t[1:(length(i)-1)]),
                  aes(y = i[1:(length(i)-1)], x = t[1:(length(i)-1)])) +
    geom_point() + xlab("Time(s)") +
    ylab("Current (A)") + theme_classic()

  if (errCheck == TRUE){
    return(list(G,Dx,Co,dt,h,l,n,Area,p,Da))
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
#' @importFrom matlab zeros ones
#' @importFrom pracma inv
#'



CottrCheck = function(Elefun) {

  FA = 96485
  R = 8.3145
  Check = Elefun
  if (length(Check) == 9){
    return("ErrCheck inside the called function should be activated")
  } else {
    vt = c(1:Check[[6]])

    if (length(Check) == 8){

      Gcot = 1/sqrt(3.14*Check[[4]]*vt)

    } else if (length(Check) == 10){
      Gcot = (1/sqrt(3.14*Check[[4]]*vt))/(1+ (1/Check[[10]])*exp(Check[[9]]))
    }

    Err = (Check[[1]]/Gcot)
    t = Check[[4]]*vt
    ErrorGraphy = ggplot(data = data.frame(Err[1:(length(Err)-1)],t[1:(length(Err)-1)]),
                         aes(y = Err[1:(length(Err)-1)], x = t[1:(length(Err)-1)])) + geom_point() +xlab("Time(s)") +
      ylab("G/Gcott") + theme_classic()
    return(ErrorGraphy)
  }
}


#' Linear Sweep digitial simulation
#'
#' Return a graph I vs E of the electrochemical process
#'
#' @param Co bulk concentration
#' @param Dx diffusion coefficient
#' @param Eo reduction potential of the species
#' @param Vi initial potential of the sweep
#' @param Vf final potential of the sweep
#' @param Vs  potential scan rate of the simulation
#' @param ko heterogeneous electron transfer rate constant
#' @param alpha charge transfer coefficient
#' @param Temp temperature in kelvin
#' @param n number of electrons involved in the process
#' @param Area area of the electrode
#' @param DerApprox number of point for the approximation of the first derivative
#' @param errCheck if true the function returns a list with parameters for CottrCheck function
#' @param Method method to be used for the simulation = "Euler" "BI" "RK4"
#'
#'
#' @return if errCheck == F a graph I vs E, if errCheck == T a list
#'
#' @examples
#' LinSwp(Co = 0.001, DerApprox = 2, errCheck = FALSE, Method = "Euler")
#'
#' @export
#' @import ggplot2
#' @importFrom matlab zeros ones
#' @importFrom pracma inv
#'

LinSwp = function(Co = 0.001, Dx = 0.00001, Eo = 0.1,
                  Vi = 0.3, Vf = -0.3, Vs = 0.001, ko = 0.01,
                  alpha = 0.5, Temp = 298.15, n = 1, Area = 1,
                  DerApprox = 2, errCheck = FALSE, Method = "Euler"){

  FA = 96485
  R = 8.3145
  f = ((FA*n)/(R*Temp))
  Da = Dx/Dx
  Dm = 0.45
  exptime = abs(Vf-Vi)/Vs
  Tmax = ceiling(abs(f*(Vf-Vi)))
  dt = (1/(f*Vs))
  l = ceiling(exptime/dt)
  j = ceiling(6*(Tmax)^0.5)
  h = sqrt((Da*dt)/Dm)
  vt = c(1:(l+1))
  t = dt*vt
  PotentialScan = Vi-Vs*t
  p = f*(Vi - Vs*t - Eo)
  KO = ko*sqrt(dt/Dx)
  Kf = KO*exp(-alpha*p)
  Kb = KO*exp((1-alpha)*p)
  Ox = ones(l +1, j)
  Red =zeros(l +1, j)
  Jox = zeros(l+1, 1)

  DerivMatrix = matrix(data = c(-1,1,0,0,0,0,
                                -3,4,-1,0,0,0,
                                -11,18,-9,2,0,0,
                                -25,48,-36,16,-3,0,
                                -137,300,-300,200,-75,12),
                       byrow = T, nrow = 5, ncol = 6)

  NormalFact = matrix(data = c(1,2,6,12,60), nrow = 5)

  if (Method == "Euler") {
    for (i1 in 1:l) {
      B = matrix(data = c((Kf[i1]*h - DerivMatrix[DerApprox-1,1]),
                          -Kb[i1]*h,
                          DerivMatrix[DerApprox-1,1],
                          DerivMatrix[DerApprox-1,1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(DerivMatrix[DerApprox -1,2:DerApprox]*Ox[i1,2:DerApprox]), -sum(DerivMatrix[DerApprox -1,2:DerApprox]*Ox[i1,2:DerApprox]) - sum(DerivMatrix[DerApprox -1,2:DerApprox]*Red[i1,2:DerApprox])))
      C = inv(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2]
      for (j1 in 2:(j-1)) {
        Ox[i1+1,j1] = Ox[i1,j1] + Dm*(Ox[i1, j1-1] - 2*Ox[i1, j1] + Ox[i1, j1+1])
        Red[i1+1, j1] = Red[i1, j1] + Dm*(Red[i1, j1-1] + Red[i1, j1+1] - 2*Red[i1,j1])
      }
      Jox[i1,1] = sum(DerivMatrix[DerApprox -1,1:DerApprox]*Ox[i1,1:DerApprox])/(h*NormalFact[DerApprox-1])
    }
  } else if (Method == "RK4") {

    for (i1 in 1:(l-1)) {
      k1 = zeros(j)
      k2 = zeros(j)
      k3 = zeros(j)
      k4 = zeros(j)
      k1red = zeros(j)
      k2red = zeros(j)
      k3red = zeros(j)
      k4red = zeros(j)
      B = matrix(data = c((Kf[i1]*h - DerivMatrix[DerApprox-1,1]),
                          -Kb[i1]*h,
                          DerivMatrix[DerApprox-1,1],
                          DerivMatrix[DerApprox-1,1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(DerivMatrix[DerApprox -1,2:DerApprox]*Ox[i1,2:DerApprox]),
                          -sum(DerivMatrix[DerApprox -1,2:DerApprox]*Ox[i1,2:DerApprox]) - sum(DerivMatrix[DerApprox -1,2:DerApprox]*Red[i1,2:DerApprox])))
      C = inv(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2]

      B = matrix(data = c((Kf[i1+1]*h - DerivMatrix[DerApprox-1,1]),
                          -Kb[i1+1]*h,
                          DerivMatrix[DerApprox-1,1],
                          DerivMatrix[DerApprox-1,1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(DerivMatrix[DerApprox -1,2:DerApprox]*Ox[i1+1,2:DerApprox]),
                          -sum(DerivMatrix[DerApprox -1,2:DerApprox]*Ox[i1+1,2:DerApprox]) - sum(DerivMatrix[DerApprox -1,2:DerApprox]*Red[i1+1,2:DerApprox])))

      C = inv(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]

      for (j1 in 2:(j-1)) {
        k1[j1] = Dm*(Ox[i1, j1 -1] - 2*Ox[i1, j1] + Ox[i1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + k1[j1]*0.5
        k1red[j1] = Dm*(Red[i1, j1 -1] - 2*Red[i1, j1] + Red[i1, j1+1])
        Red[i1 + 1,j1] = Red[i1,j1] + k1red[j1]*0.5
      }
      for (j1 in 2:(j-1)) {
        k2[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + k2[j1]*0.5
        k2red[j1] = Dm*(Red[i1 + 1, j1 -1] - 2*Red[i1 + 1, j1] + Red[i1 + 1, j1+1])
        Red[i1 + 1,j1] = Red[i1,j1] + k2red[j1]*0.5
      }
      for (j1 in 2:(j-1)) {
        k3[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + k3[j1]
        k3red[j1] = Dm*(Red[i1 + 1, j1 -1] - 2*Red[i1 + 1, j1] + Red[i1 + 1, j1+1])
        Red[i1 + 1,j1] = Red[i1,j1] + k3red[j1]
      }
      for (j1 in 2:(j-1)) {
        k4[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + (k1[j1] + 2*k2[j1] + 2*k3[j1] + k4[j1])/6
        k4red[j1] = Dm*(Red[i1 + 1, j1 -1] - 2*Red[i1 + 1, j1] + Red[i1 + 1, j1+1])
        Red[i1 + 1,j1] = Red[i1,j1] + (k1red[j1] + 2*k2red[j1] + 2*k3red[j1] + k4red[j1])/6
      }
      Jox[i1,1] = sum(DerivMatrix[DerApprox -1,1:DerApprox]*Ox[i1,1:DerApprox])/(h*NormalFact[DerApprox-1])
    }
  } else if (Method == "BI") {
    al1 = 1/(h^2)
    al2 = -2/(h^2)
    al3 = 1/(h^2)
    a1 = (al2 - 1/dt)/al1
    a2 = al3/al1

    for (i1 in 1:(l-1)) {

      B = matrix(data = c((Kf[i1]*h - DerivMatrix[DerApprox-1,1]),
                          -Kb[i1]*h,
                          DerivMatrix[DerApprox-1,1],
                          DerivMatrix[DerApprox-1,1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(DerivMatrix[DerApprox -1,2:DerApprox]*Ox[i1,2:DerApprox]),
                          -sum(DerivMatrix[DerApprox -1,2:DerApprox]*Ox[i1,2:DerApprox]) - sum(DerivMatrix[DerApprox -1,2:DerApprox]*Red[i1,2:DerApprox])))
      C = inv(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2]

      B = matrix(data = c((Kf[i1+1]*h - DerivMatrix[DerApprox-1,1]),
                          -Kb[i1+1]*h,
                          DerivMatrix[DerApprox-1,1],
                          DerivMatrix[DerApprox-1,1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(DerivMatrix[DerApprox -1,2:DerApprox]*Ox[i1+1,2:DerApprox]),
                          -sum(DerivMatrix[DerApprox -1,2:DerApprox]*Ox[i1+1,2:DerApprox]) - sum(DerivMatrix[DerApprox -1,2:DerApprox]*Red[i1+1,2:DerApprox])))

      C = inv(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]
      bOx = (-Ox[i1,]/(al1*dt))
      bRed = (-Red[i1,]/(al1*dt))
      A = c(rep(1,j))
      A1 = c(rep(a1,j))
      A2 = c(rep(a2,j))
      bOx[j] = bOx[j] - A2[j]*1
      bRed[j] = bRed[j] - A2[j]*1

      for (j1 in ((j-1):2)) {

        bOx[j1] = bOx[j1] - A2[j]*bOx[j1+1]/A1[j1+1]
        bRed[j1] = bRed[j1] - A2[j]*bRed[j1+1]/A1[j1+1]
        A1[j1] = A1[j1] - A2[j]/A1[j1+1]

      }

      for (j1 in 2:j) {
        Ox[i1+1,j1] = (bOx[j1] -Ox[i1+1,j1-1])/A1[j1]
        Red[i1+1,j1] = (bRed[j1] -Red[i1+1,j1-1])/A1[j1]
      }
      Jox[i1,1] = sum(DerivMatrix[DerApprox -1,1:DerApprox]*Ox[i1,1:DerApprox])/(h*NormalFact[DerApprox-1])
    }
  } else if (!(Method %in% c("Euler", "BI", "RK4"))) {
    return("Available methods are Euler, BI and RK4")
  }

  G = Jox

  i = (n*FA*G*Dx*Area*Co)/(6.4*h*sqrt(Dx))


  graphy = ggplot(data = data.frame(i[1:(length(i)-1)],PotentialScan[1:(length(i)-1)]), aes(y = i[1:(length(i)-1)], x = PotentialScan[1:(length(i)-1)])) +
    geom_point() + scale_x_continuous(trans = "reverse") +
    xlab("Potential (V)") +
    ylab("Current (A)") +
    theme_classic()

  if (errCheck == TRUE){
    return(list(G,Dx,Co,dt,h,l,n,Area,p,Da))
  } else {
    return(graphy)
  }

}

#' Cyclic voltammetry digitial simulation
#'
#' Return a graph I vs E of the electrochemical process
#'
#' @param Co bulk concentration
#' @param Dx diffusion coefficient
#' @param Eo reduction potential of the species
#' @param Vi initial potential of the sweep
#' @param Vf final potential of the sweep
#' @param Vs  potential scan rate of the simulation
#' @param ko heterogeneous electron transfer rate constant
#' @param alpha charge transfer coefficient
#' @param Temp temperature in kelvin
#' @param n number of electrons involved in the process
#' @param Area area of the electrode
#' @param DerApprox number of point for the approximation of the first derivative
#' @param errCheck if true the function returns a list with parameters for CottrCheck function
#' @param Method method to be used for the simulation = "Euler" "BI" "RK4"
#'
#'
#' @return if errCheck == F a graph I vs E, if errCheck == T a list
#'
#' @examples
#' CV(Co = 0.001, DerApprox = 2, errCheck = FALSE, Method = "Euler")
#'
#' @export
#' @import ggplot2
#' @importFrom matlab zeros ones
#' @importFrom pracma inv
#'

CV = function(Co = 0.001, Dx = 0.00001, Eo = 0.1,
              Vi = 0.3, Vf = -0.3, Vs = 0.001, ko = 0.01,
              alpha = 0.5, Temp = 298.15, n = 1, Area = 1,
              DerApprox = 2, errCheck = FALSE, Method = "Euler"){

  FA = 96485
  R = 8.3145
  f = ((FA*n)/(R*Temp))
  Da = Dx/Dx
  Dm = 0.45
  exptime = 2*abs(Vf-Vi)/Vs
  Tmax = ceiling(abs(f*(Vf-Vi)))
  dt = (1/(f*Vs))
  l = ceiling(exptime/dt)
  j = ceiling(6*(Tmax)^0.5)
  h = sqrt((Da*dt)/Dm)
  vt = c(1:(l+1))
  t = dt*vt
  forwardScan = Vi-Vs*t[1:((l/2) +1)]
  backwardScan = Vf + Vs*t[1:((l/2) +1)]
  PotentialScan = c(forwardScan, backwardScan)
  p = f*(PotentialScan - Eo)
  KO = ko*sqrt(dt/Dx)
  Kf = KO*exp(-alpha*p)
  Kb = KO*exp((1-alpha)*p)
  Ox = ones(l +1, j)
  Red =zeros(l +1, j)
  Jox = zeros(l+1, 1)

  DerivMatrix = matrix(data = c(-1,1,0,0,0,0,
                                -3,4,-1,0,0,0,
                                -11,18,-9,2,0,0,
                                -25,48,-36,16,-3,0,
                                -137,300,-300,200,-75,12),
                       byrow = T, nrow = 5, ncol = 6)

  NormalFact = matrix(data = c(1,2,6,12,60), nrow = 5)

  if (Method == "Euler") {
    for (i1 in 1:l) {
      B = matrix(data = c((Kf[i1]*h - DerivMatrix[DerApprox-1,1]),
                          -Kb[i1]*h,
                          DerivMatrix[DerApprox-1,1],
                          DerivMatrix[DerApprox-1,1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(DerivMatrix[DerApprox -1,2:DerApprox]*Ox[i1,2:DerApprox]), -sum(DerivMatrix[DerApprox -1,2:DerApprox]*Ox[i1,2:DerApprox]) - sum(DerivMatrix[DerApprox -1,2:DerApprox]*Red[i1,2:DerApprox])))
      C = inv(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2]
      for (j1 in 2:(j-1)) {
        Ox[i1+1,j1] = Ox[i1,j1] + Dm*(Ox[i1, j1-1] - 2*Ox[i1, j1] + Ox[i1, j1+1])
        Red[i1+1, j1] = Red[i1, j1] + Dm*(Red[i1, j1-1] + Red[i1, j1+1] - 2*Red[i1,j1])
      }
      Jox[i1,1] = sum(DerivMatrix[DerApprox -1,1:DerApprox]*Ox[i1,1:DerApprox])/(h*NormalFact[DerApprox-1])
    }
  } else if (Method == "RK4") {

    for (i1 in 1:(l-1)) {
      k1 = zeros(j)
      k2 = zeros(j)
      k3 = zeros(j)
      k4 = zeros(j)
      k1red = zeros(j)
      k2red = zeros(j)
      k3red = zeros(j)
      k4red = zeros(j)
      B = matrix(data = c((Kf[i1]*h - DerivMatrix[DerApprox-1,1]),
                          -Kb[i1]*h,
                          DerivMatrix[DerApprox-1,1],
                          DerivMatrix[DerApprox-1,1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(DerivMatrix[DerApprox -1,2:DerApprox]*Ox[i1,2:DerApprox]),
                          -sum(DerivMatrix[DerApprox -1,2:DerApprox]*Ox[i1,2:DerApprox]) - sum(DerivMatrix[DerApprox -1,2:DerApprox]*Red[i1,2:DerApprox])))
      C = inv(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2]

      B = matrix(data = c((Kf[i1+1]*h - DerivMatrix[DerApprox-1,1]),
                          -Kb[i1+1]*h,
                          DerivMatrix[DerApprox-1,1],
                          DerivMatrix[DerApprox-1,1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(DerivMatrix[DerApprox -1,2:DerApprox]*Ox[i1+1,2:DerApprox]),
                          -sum(DerivMatrix[DerApprox -1,2:DerApprox]*Ox[i1+1,2:DerApprox]) - sum(DerivMatrix[DerApprox -1,2:DerApprox]*Red[i1+1,2:DerApprox])))
      C = inv(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]

      for (j1 in 2:(j-1)) {
        k1[j1] = Dm*(Ox[i1, j1 -1] - 2*Ox[i1, j1] + Ox[i1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + k1[j1]*0.5
        k1red[j1] = Dm*(Red[i1, j1 -1] - 2*Red[i1, j1] + Red[i1, j1+1])
        Red[i1 + 1,j1] = Red[i1,j1] + k1red[j1]*0.5
      }
      for (j1 in 2:(j-1)) {
        k2[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + k2[j1]*0.5
        k2red[j1] = Dm*(Red[i1 + 1, j1 -1] - 2*Red[i1 + 1, j1] + Red[i1 + 1, j1+1])
        Red[i1 + 1,j1] = Red[i1,j1] + k2red[j1]*0.5
      }
      for (j1 in 2:(j-1)) {
        k3[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + k3[j1]
        k3red[j1] = Dm*(Red[i1 + 1, j1 -1] - 2*Red[i1 + 1, j1] + Red[i1 + 1, j1+1])
        Red[i1 + 1,j1] = Red[i1,j1] + k3red[j1]
      }
      for (j1 in 2:(j-1)) {
        k4[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + (k1[j1] + 2*k2[j1] + 2*k3[j1] + k4[j1])/6
        k4red[j1] = Dm*(Red[i1 + 1, j1 -1] - 2*Red[i1 + 1, j1] + Red[i1 + 1, j1+1])
        Red[i1 + 1,j1] = Red[i1,j1] + (k1red[j1] + 2*k2red[j1] + 2*k3red[j1] + k4red[j1])/6
      }
      Jox[i1,1] = sum(DerivMatrix[DerApprox -1,1:DerApprox]*Ox[i1,1:DerApprox])/(h*NormalFact[DerApprox-1])
    }
  } else if (Method == "BI") {
    al1 = 1/(h^2)
    al2 = -2/(h^2)
    al3 = 1/(h^2)
    a1 = (al2 - 1/dt)/al1
    a2 = al3/al1

    for (i1 in 1:(l-1)) {

      B = matrix(data = c((Kf[i1]*h - DerivMatrix[DerApprox-1,1]),
                          -Kb[i1]*h,
                          DerivMatrix[DerApprox-1,1],
                          DerivMatrix[DerApprox-1,1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(DerivMatrix[DerApprox -1,2:DerApprox]*Ox[i1,2:DerApprox]),
                          -sum(DerivMatrix[DerApprox -1,2:DerApprox]*Ox[i1,2:DerApprox]) - sum(DerivMatrix[DerApprox -1,2:DerApprox]*Red[i1,2:DerApprox])))
      C = inv(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2]

      B = matrix(data = c((Kf[i1+1]*h - DerivMatrix[DerApprox-1,1]),
                          -Kb[i1+1]*h,
                          DerivMatrix[DerApprox-1,1],
                          DerivMatrix[DerApprox-1,1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(DerivMatrix[DerApprox -1,2:DerApprox]*Ox[i1+1,2:DerApprox]),
                          -sum(DerivMatrix[DerApprox -1,2:DerApprox]*Ox[i1+1,2:DerApprox]) - sum(DerivMatrix[DerApprox -1,2:DerApprox]*Red[i1+1,2:DerApprox])))
      C = inv(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]
      bOx = (-Ox[i1,]/(al1*dt))
      bRed = (-Red[i1,]/(al1*dt))
      A = c(rep(1,j))
      A1 = c(rep(a1,j))
      A2 = c(rep(a2,j))
      bOx[j] = bOx[j] - A2[j]*1
      bRed[j] = bRed[j] - A2[j]*1

      for (j1 in ((j-1):2)) {

        bOx[j1] = bOx[j1] - A2[j]*bOx[j1+1]/A1[j1+1]
        bRed[j1] = bRed[j1] - A2[j]*bRed[j1+1]/A1[j1+1]
        A1[j1] = A1[j1] - A2[j]/A1[j1+1]

      }

      for (j1 in 2:j) {
        Ox[i1+1,j1] = (bOx[j1] -Ox[i1+1,j1-1])/A1[j1]
        Red[i1+1,j1] = (bRed[j1] -Red[i1+1,j1-1])/A1[j1]
      }
      Jox[i1,1] = sum(DerivMatrix[DerApprox -1,1:DerApprox]*Ox[i1,1:DerApprox])/(h*NormalFact[DerApprox-1])
    }
  } else if (!(Method %in% c("Euler", "BI", "RK4"))) {
    return("Available methods are Euler, BI and RK4")
  }

  G = Jox

  i = (n*FA*G*Dx*Area*Co)/(6.4*h*sqrt(Dx))

  graphy = ggplot(data = data.frame(i,PotentialScan[1:length(i)]), aes(y = i, x = PotentialScan[1:length(i)])) +
    geom_point() + scale_x_continuous(trans = "reverse") +
    xlab("Potential (V)") +
    ylab("Current (A)") +
    theme_classic()

  if (errCheck == TRUE){
    return(list(i,Dx,Co,dt,h,l,n,Area,p,Da))
  } else {
    return(graphy)
  }

}


#' EC behaviour cyclic voltammetry simulator
#'
#' Return a graph I vs E of the electrochemical process
#'
#' @param Co bulk concentration
#' @param Dx diffusion coefficient
#' @param Eo reduction potential of the species
#' @param Vi initial potential of the sweep
#' @param Vf final potential of the sweep
#' @param Vs  potential scan rate of the simulation
#' @param ko heterogeneous electron transfer rate constant
#' @param kc rate constant of the reaction Red -> C
#' @param alpha charge transfer coefficient
#' @param Temp temperature in kelvin
#' @param n number of electrons involved in the process
#' @param Area area of the electrode
#' @param DerApprox number of point for the approximation of the first derivative
#' @param errCheck if true the function returns a list with parameters for CottrCheck function
#' @param Method method to be used for the simulation = "Euler" "BI" "RK4"
#'
#'
#' @return if errCheck == F a graph I vs E, if errCheck == T a list
#'
#' @examples
#' CVEC(Co = 0.001, DerApprox = 2, kc = 0.00001, errCheck = FALSE, Method = "Euler")
#'
#' @export
#' @import ggplot2
#' @importFrom matlab zeros ones
#' @importFrom pracma inv
#'

CVEC = function(Co = 0.001, Dx = 0.00001, Eo = 0.1,
                Vi = 0.3, Vf = -0.3, Vs = 0.001, ko = 0.01,
                kc = 0.0001,
                alpha = 0.5, Temp = 298.15, n = 1, Area = 1,
                DerApprox = 2, errCheck = FALSE, Method = "Euler"){

  FA = 96485
  R = 8.3145
  f = ((FA*n)/(R*Temp))
  Da = Dx/Dx #normalized coefficient of diffusion
  Dm = 0.45 #Da*dt/h^2
  exptime = 2*abs(Vf-Vi)/Vs #one second as simulation time
  Tmax = ceiling(abs(f*(Vf-Vi)))
  dt = (1/(f*Vs)) #se non cambiasse con Vs la corrente sarebbe dipendente dal tempo osservato, quello però che possiamo fare per aumentare il numero di punti è dividere per una quantità K e poi rimoltiplicarci la corrente
  l = ceiling(exptime/dt) #time interval
  j = ceiling(6*(Tmax)^0.5) #diffusion layer reached in maximum time space (maximum expantion of space layer, by obtaining the space interval h from a fixed Dm we can obtain the number of space sections)
  h = sqrt((Da*dt)/Dm) #1 is one unit of time, by setting l = 100 we divided time in 100 spaces
  vt = c(1:(l+1)) #time index vector
  t = dt*vt
  forwardScan = Vi-Vs*t[1:((l/2) +1)]
  backwardScan = Vf + Vs*t[1:((l/2) +1)]
  PotentialScan = c(forwardScan, backwardScan)
  p = f*(PotentialScan - Eo)
  KC = kc*dt
  KO = ko*sqrt(dt/Dx)
  Kf = KO*exp(-alpha*p)
  Kb = KO*exp((1-alpha)*p)
  Ox = ones(l +1, j) #matrix of the C specie normalized by bulk concentration
  Red =zeros(l +1, j)
  Jox = zeros(l+1, 1) #flux matrix

  DerivMatrix = matrix(data = c(-1,1,0,0,0,0,
                                -3,4,-1,0,0,0,
                                -11,18,-9,2,0,0,
                                -25,48,-36,16,-3,0,
                                -137,300,-300,200,-75,12),
                       byrow = T, nrow = 5, ncol = 6)

  NormalFact = matrix(data = c(1,2,6,12,60), nrow = 5)

  if (Method == "Euler") {
    for (i1 in 1:l) {
      B = matrix(data = c((Kf[i1]*h - DerivMatrix[DerApprox-1,1]), -Kb[i1]*h ,DerivMatrix[DerApprox-1,1],DerivMatrix[DerApprox-1,1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(DerivMatrix[DerApprox -1,2:DerApprox]*Ox[i1,2:DerApprox]), -sum(DerivMatrix[DerApprox -1,2:DerApprox]*Ox[i1,2:DerApprox]) - sum(DerivMatrix[DerApprox -1,2:DerApprox]*Red[i1,2:DerApprox])))
      C = inv(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2] - KC*Red[i1,1]
      for (j1 in 2:(j-1)) {
        Ox[i1+1,j1] = Ox[i1,j1] + Dm*(Ox[i1, j1-1] - 2*Ox[i1, j1] + Ox[i1, j1+1])
        Red[i1+1, j1] = Red[i1, j1] + Dm*(Red[i1, j1-1] + Red[i1, j1+1] - 2*Red[i1,j1]) - KC*dt*Red[i1,j1]
      }

      Jox[i1,1] = sum(DerivMatrix[DerApprox -1,1:DerApprox]*Ox[i1,1:DerApprox])/(h*NormalFact[DerApprox-1])

    }
  } else if (Method == "BI") {
    al1 = 1/(h^2)
    al2 = -2/(h^2)
    al3 = 1/(h^2)
    a1 = (al2 - 1/dt)/al1
    a2 = al3/al1

    for (i1 in 1:(l-1)) {

      B = matrix(data = c((Kf[i1]*h - DerivMatrix[DerApprox-1,1]),
                          -Kb[i1]*h,
                          DerivMatrix[DerApprox-1,1],
                          DerivMatrix[DerApprox-1,1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(DerivMatrix[DerApprox -1,2:DerApprox]*Ox[i1,2:DerApprox]),
                          -sum(DerivMatrix[DerApprox -1,2:DerApprox]*Ox[i1,2:DerApprox]) - sum(DerivMatrix[DerApprox -1,2:DerApprox]*Red[i1,2:DerApprox])))
      C = inv(B) %*% Y

      Ox[i1,1] = C[1]
      Red[i1,1] = C[2] - KC*Red[i1,1]


      B = matrix(data = c((Kf[i1+1]*h - DerivMatrix[DerApprox-1,1]),
                          -Kb[i1+1]*h,
                          DerivMatrix[DerApprox-1,1],
                          DerivMatrix[DerApprox-1,1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(DerivMatrix[DerApprox -1,2:DerApprox]*Ox[i1+1,2:DerApprox]),
                          -sum(DerivMatrix[DerApprox -1,2:DerApprox]*Ox[i1+1,2:DerApprox]) - sum(DerivMatrix[DerApprox -1,2:DerApprox]*Red[i1+1,2:DerApprox])))
      C = inv(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2] - KC*Red[i1+1,1]

      bOx = (-Ox[i1,]/(al1*dt))
      bRed = (-Red[i1,]/(al1*dt))
      A = c(rep(1,j))
      A1 = c(rep(a1,j))
      A1red = c(rep(a1,j))
      A2 = c(rep(a2,j))
      bOx[j] = bOx[j] - A2[j]*1
      bRed[j] = bRed[j] - A2[j]*(Red[i1+1,j] - KC*Red[i1+1,j]/al1)
      A1red[1] = A1red[1] - KC/al1


      for (j1 in ((j-1):2)) {

        bOx[j1] = bOx[j1] - A2[j]*bOx[j1+1]/A1[j1+1]
        bRed[j1] = bRed[j1] - A2[j]*bRed[j1+1]/A1red[j1+1]
        A1[j1] = A1[j1] - A2[j]/A1[j1+1]
        A1red[j1] = A1red[j1] - A2[j]/A1red[j1+1] -KC/al1
      }

      for (j1 in 2:j) {
        Ox[i1+1,j1] = (bOx[j1] -Ox[i1+1,j1-1])/A1[j1]
        Red[i1+1,j1] = (bRed[j1] -Red[i1+1,j1-1])/A1red[j1]
      }
      Jox[i1,1] = sum(DerivMatrix[DerApprox -1,1:DerApprox]*Ox[i1,1:DerApprox])/(h*NormalFact[DerApprox-1])
    }
  } else if (!(Method %in% c("Euler", "BI", "RK4"))) {
    return("Available methods are Euler, BI and RK4")
  }

  G = Jox

  i = (n*FA*G*Dx*Area*Co)/(10*h*sqrt(Dx))

  graphy = ggplot(data = data.frame(i,PotentialScan[1:length(i)]), aes(y = i, x = PotentialScan[1:length(i)])) +
    geom_point() + scale_x_continuous(trans = "reverse") +
    xlab("Overpotential (V)") +
    ylab("Current (A)") +
    theme_classic()

  if (errCheck == TRUE){
    return(list(G,Dx,Co,dt,h,l,n,Area,p,Da))
  } else {
    return(graphy)
  }

}
