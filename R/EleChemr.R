#' Chrono amperometry digital simulation
#'
#' Return a graph I vs t of the electrochemical process
#'
#' @param Co bulk concentration
#' @param exptime experimental time to be simulated
#' @param Dx diffusion coefficient
#' @param Dm simulation parameter, maximum 0.5 for explicit methods
#' @param Temp temperature in kelvin
#' @param n number of electrons involved in the process
#' @param Area area of the electrode
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
                    Temp = 298.15, n = 1, Area = 1, DerApprox = 2,
                    errCheck = FALSE, Method = "Euler") {

  FA = 96485
  R = 8.3145
  f = ((FA*n)/(R*Temp))
  dt = 0.01
  Da = Dx/Dx
  l = exptime/dt
  h = sqrt((Da*dt)/Dm)
  j = ceiling(6*(l)^0.5)
  vt = c(1:l)
  t = dt*vt
  Ox = OneMat(l, j)
  Jox = ZeroMat(l, 1)

  if (Method == "Euler") {

    for (i1 in 1:(l-1)) {

      Ox[i1,1] = 0

      for (j1 in 2:(j-1)) {
        Ox[i1 + 1,j1] = Ox[i1,j1] + Dm*(Ox[i1, j1 -1] - 2*Ox[i1, j1] + Ox[i1, j1+1])
      }
    }

    Jox = Derv(Ox = Ox, h = h, npoints = DerApprox)

  } else if (Method == "RK4") {

    for (i1 in 1:(l-1)) {
      k1 = ZeroMat(j)
      k2 = ZeroMat(j)
      k3 = ZeroMat(j)
      k4 = ZeroMat(j)
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
    }

    Jox = Derv(Ox = Ox, h = h, npoints = DerApprox)

  } else if (Method == "BI") {
    al1 = 1/(h^2)
    al2 = -2/(h^2)
    al3 = 1/(h^2)
    a1 = (al2 - 1/dt)/al1
    a2 = al3/al1

    for (i1 in 1:(l-1)) {
      Y = ZeroMat(j-2,j-2)
      Y[1,1] = a1
      Y[1,2] = a2
      Y[j-2,j-3] = 1
      Y[j-2,j-2] = a1
      for (i in 2:(j-3)) {
        Y[i,i] = a1
        Y[i,i-1] = 1
        Y[i, i +1] = a2
      }

      Ox[i1,1] = 0
      Ox[i1+1,1] = 0
      b = (-Ox[i1,2:(j-1)]/(al1*dt))
      b[j-2] = b[j-2] - a2*1
      b[1] = b[1] - Ox[i1+1,1]
      Ox[i1+1,2:(j-1)] = solve(Y) %*% b
    }

    Jox = Derv(Ox = Ox, h = h, npoints = DerApprox)

  } else if (Method == "CN") {

    al1 = 1/(h^2)
    al2 = -2/(h^2)
    al3 = 1/(h^2)
    a1 = (al2 - 2/dt)/al1
    a2 = al3/al1
    a3 = (al2 + 2/dt)/al1

    for (i1 in 1:(l-1)) {
      Y = ZeroMat(j-2,j-2)
      Y[1,1] = a1
      Y[1,2] = a2
      Y[j-2,j-3] = 1
      Y[j-2,j-2] = a1
      for (i in 2:(j-3)) {
        Y[i,i] = a1
        Y[i,i-1] = 1
        Y[i, i +1] = a2
      }

      Ox[i1,1] = 0
      Ox[i1+1,1] = 0
      b = -a3*Ox[i1,2:(j-1)] - Ox[i1,1:((j-1)-1)] - a2*Ox[i1,3:((j-1)+1)]
      b[j-2] = b[j-2] - a2*1
      b[1] = b[1] - Ox[i1+1,1]
      Ox[i1+1,2:(j-1)] = solve(Y) %*% b
    }

    Jox = Derv(Ox = Ox, h = h, npoints = DerApprox)


  } else if (Method == "BDF") {

    al1 = 1/(h^2)
    al2 = -2/(h^2)
    al3 = 1/(h^2)
    a1 = (al2 - 1.5/dt)/al1
    a2 = al3/al1

    for (i1 in 1:(l-1)) {
      Y = ZeroMat(j-2,j-2)
      Y[1,1] = a1
      Y[1,2] = a2
      Y[j-2,j-3] = 1
      Y[j-2,j-2] = a1
      for (i in 2:(j-3)) {
        Y[i,i] = a1
        Y[i,i-1] = 1
        Y[i, i +1] = a2
      }

      Ox[i1,1] = 0
      Ox[i1+1,1] = 0
      if (i1 == 1) {
        b = -2*Ox[i1,2:(j-1)]/(dt*al1) + Ox[1,2:(j-1)]/(2*dt*al1)
      } else {
        b = -2*Ox[i1,2:(j-1)]/(dt*al1) + Ox[i1-1,2:(j-1)]/(2*dt*al1)
      }
      b[j-2] = b[j-2] - a2*1
      b[1] = b[1] - Ox[i1+1,1]
      Ox[i1+1,2:(j-1)] = solve(Y) %*% b
    }

    Jox = Derv(Ox = Ox, h = h, npoints = DerApprox)


  } else if (!(Method %in% c("Euler", "BI", "RK4", "CN", "BDF"))) {
    return("Available methods are Euler, BI, RK4, CN and BDF")
  }

  G = Jox
  i = (n*FA*G*Dx*Area*Co)/(6.4*(h*(Dx^0.5)))

  graphy = ggplot(data = data.frame(i[1:(length(i)-1)],t[1:(length(i)-1)]),
                  aes(y = i[1:(length(i)-1)], x = t[1:(length(i)-1)])) +
    geom_point() + xlab("Time(s)") +
    ylab("Current (A)") + theme_classic()

  if (errCheck == TRUE){
    return(list(G,Dx,Co,dt,h,l,j,i,n,Area))
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
#' @param Dm simulation parameter, maximum 0.5 for explicit methods
#' @param Temp temperature in kelvin
#' @param n number of electrons involved in the process
#' @param Area area of the electrode
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
                   eta = 0.1, Temp = 298.15, n = 1, Area = 1,
                   DerApprox = 2, errCheck = FALSE, Method = "Euler") {

  FA = 96485
  R = 8.3145
  f = ((FA*n)/(R*Temp))
  dt = 0.01
  p = eta*f
  Da = Dx/Dx
  l = exptime/dt
  h = sqrt((Da*dt)/Dm)
  j = ceiling(6*(Da*l)^0.5)
  vt = c(1:l)
  t = dt*vt
  Ox = OneMat(l, j)
  Red = ZeroMat(l, j)
  Jox = ZeroMat(l, 1)

  if (Method == "Euler") {
    for (i1 in 1:(l-1)) {
      B = matrix(data = c(1,-exp(p),
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(0, -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2]
      for (j1 in 2:(j-1)) {
        Ox[i1+1,j1] = Ox[i1,j1] + Dm*(Ox[i1, j1-1] - 2*Ox[i1, j1] + Ox[i1, j1+1])
        Red[i1+1, j1] = Red[i1, j1] + Dm*(Red[i1, j1-1] + Red[i1, j1+1] - 2*Red[i1,j1])
      }
    }

    Jox = Derv(Ox = Ox, h = h, npoints = DerApprox)

  } else if (Method == "RK4") {

    for (i1 in 1:(l-1)) {
      k1 = ZeroMat(j)
      k2 = ZeroMat(j)
      k3 = ZeroMat(j)
      k4 = ZeroMat(j)
      k1red = ZeroMat(j)
      k2red = ZeroMat(j)
      k3red = ZeroMat(j)
      k4red = ZeroMat(j)
      B = matrix(data = c(1,-exp(p),
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(0, -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2]

      for (j1 in 2:(j-1)) {
        k1[j1] = Dm*(Ox[i1, j1 -1] - 2*Ox[i1, j1] + Ox[i1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + k1[j1]*0.5
        k1red[j1] = Dm*(Red[i1, j1 -1] - 2*Red[i1, j1] + Red[i1, j1+1])
        Red[i1 + 1,j1] = Red[i1,j1] + k1red[j1]*0.5
      }

      B = matrix(data = c(1,-exp(p),
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(0, -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1+1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]

      for (j1 in 2:(j-1)) {
        k2[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + k2[j1]*0.5
        k2red[j1] = Dm*(Red[i1 + 1, j1 -1] - 2*Red[i1 + 1, j1] + Red[i1 + 1, j1+1])
        Red[i1 + 1,j1] = Red[i1,j1] + k2red[j1]*0.5
      }

      B = matrix(data = c(1,-exp(p),
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(0, -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1+1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]

      for (j1 in 2:(j-1)) {
        k3[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + k3[j1]
        k3red[j1] = Dm*(Red[i1 + 1, j1 -1] - 2*Red[i1 + 1, j1] + Red[i1 + 1, j1+1])
        Red[i1 + 1,j1] = Red[i1,j1] + k3red[j1]
      }

      B = matrix(data = c(1,-exp(p),Derv(npoints = DerApprox, CoefMat = T)[1], Derv(npoints = DerApprox, CoefMat = T)[1]), byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(0, -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1+1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]

      for (j1 in 2:(j-1)) {
        k4[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + (k1[j1] + 2*k2[j1] + 2*k3[j1] + k4[j1])/6
        k4red[j1] = Dm*(Red[i1 + 1, j1 -1] - 2*Red[i1 + 1, j1] + Red[i1 + 1, j1+1])
        Red[i1 + 1,j1] = Red[i1,j1] + (k1red[j1] + 2*k2red[j1] + 2*k3red[j1] + k4red[j1])/6
      }
    }

    Jox = Derv(Ox = Ox, h = h, npoints = DerApprox)

  }  else if (Method == "BI") {
    al1 = 1/(h^2)
    al2 = -2/(h^2)
    al3 = 1/(h^2)
    a1 = (al2 - 1/dt)/al1
    a2 = al3/al1

    for (i1 in 1:(l-1)) {

      B = matrix(data = c(1,-exp(p),
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(0, -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2]


      bOx = (-Ox[i1,(2:(j-1))]/(al1*dt))
      bRed = (-Red[i1,(2:(j-1))]/(al1*dt))
      A = c(rep(1,j-2))
      A1 = c(rep(a1,j-2))
      A2 = c(rep(a2,j-2))
      bOx[j-2] = bOx[j-2] - A2[j-2]*1
      bRed[j-2] = bRed[j-2] - A2[j-2]*Red[i1,j]
      uox = c(rep(0, DerApprox))
      vox = c(rep(1, DerApprox))
      uRed = c(rep(0, DerApprox))
      vRed = c(rep(1, DerApprox))

      for (j1 in ((j-3):1)) {

        bOx[j1] = bOx[j1] - A2[j-2]*bOx[j1+1]/A1[j1+1]
        bRed[j1] = bRed[j1] - A2[j-2]*bRed[j1+1]/A1[j1+1]
        A1[j1] = A1[j1] - A2[j-2]/A1[j1+1]

      }
      for (m in 2:DerApprox) {
        uox[m] = (bOx[m-1] - uox[m-1])/A1[m-1]
        vox[m] = -vox[m-1]/A1[m-1]
        uRed[m] = (bRed[m-1] - uRed[m-1])/A1[m-1]
        vRed[m] = -vRed[m-1]/A1[m-1]
      }

      B = matrix(data = c(1,-exp(p),
                          Derv(npoints = DerApprox, CoefMat = T)[1] + sum(vox[2:DerApprox]*Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]),
                          Derv(npoints = DerApprox, CoefMat = T)[1] + sum(vRed[2:DerApprox]*Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox])),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(0, -sum(uox[2:DerApprox]*Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]) - sum(uRed[2:DerApprox]*Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]


      for (j1 in 1:(j-2)) {
        Ox[i1+1,j1+1] = (bOx[j1] -Ox[i1+1,j1])/A1[j1]
        Red[i1+1,j1+1] = (bRed[j1] -Red[i1+1,j1])/A1[j1]
      }
    }

    Jox = Derv(Ox = Ox, h = h, npoints = DerApprox)

  } else if (Method == "CN") {

    al1 = 1/(h^2)
    al2 = -2/(h^2)
    al3 = 1/(h^2)
    a1 = (al2 - 2/dt)/al1
    a2 = al3/al1
    a3 = (al2 + 2/dt)/al1

    for (i1 in 1:(l-1)) {

      B = matrix(data = c(1,-exp(p),
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(0, -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2]


      bOx = -a3*Ox[i1,(2:(j-1))] - Ox[i1,(1:(j-2))] - a2*Ox[i1,(3:j)]
      bRed = -a3*Red[i1,(2:(j-1))]- Red[i1,(1:(j-2))] - a2*Red[i1,(3:j)]
      A = c(rep(1,j-2))
      A1 = c(rep(a1,j-2))
      A2 = c(rep(a2,j-2))
      bOx[j-2] = bOx[j-2] - A2[j-2]*1
      bRed[j-2] = bRed[j-2] - A2[j-2]*Red[i1,j]
      uox = c(rep(0, DerApprox))
      vox = c(rep(1, DerApprox))
      uRed = c(rep(0, DerApprox))
      vRed = c(rep(1, DerApprox))

      for (j1 in ((j-3):1)) {

        bOx[j1] = bOx[j1] - A2[j-2]*bOx[j1+1]/A1[j1+1]
        bRed[j1] = bRed[j1] - A2[j-2]*bRed[j1+1]/A1[j1+1]
        A1[j1] = A1[j1] - A2[j-2]/A1[j1+1]

      }
      for (m in 2:DerApprox) {
        uox[m] = (bOx[m-1] - uox[m-1])/A1[m-1]
        vox[m] = -vox[m-1]/A1[m-1]
        uRed[m] = (bRed[m-1] - uRed[m-1])/A1[m-1]
        vRed[m] = -vRed[m-1]/A1[m-1]
      }

      B = matrix(data = c(1,-exp(p),
                          Derv(npoints = DerApprox, CoefMat = T)[1] + sum(vox[2:DerApprox]*Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]),
                          Derv(npoints = DerApprox, CoefMat = T)[1] + sum(vRed[2:DerApprox]*Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox])),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(0, -sum(uox[2:DerApprox]*Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]) - sum(uRed[2:DerApprox]*Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]


      for (j1 in 1:(j-2)) {
        Ox[i1+1,j1+1] = (bOx[j1] -Ox[i1+1,j1])/A1[j1]
        Red[i1+1,j1+1] = (bRed[j1] -Red[i1+1,j1])/A1[j1]
      }
    }

    Jox = Derv(Ox = Ox, h = h, npoints = DerApprox)


  } else if (Method == "BDF") {

    al1 = 1/(h^2)
    al2 = -2/(h^2)
    al3 = 1/(h^2)
    a1 = (al2 - 1.5/dt)/al1
    a2 = al3/al1


    for (i1 in 1:(l-1)) {

      B = matrix(data = c(1,-exp(p),
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(0, -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2]

      if (i1 == 1) {
        bOx = -2*Ox[i1,2:(j-1)]/(dt*al1) + Ox[1,2:(j-1)]/(2*dt*al1)
        bRed = -2*Red[i1,2:(j-1)]/(dt*al1) + Red[1,2:(j-1)]/(2*dt*al1)
      } else {
        bOx = -2*Ox[i1,2:(j-1)]/(dt*al1) + Ox[i1-1,2:(j-1)]/(2*dt*al1)
        bRed = -2*Red[i1,2:(j-1)]/(dt*al1) + Red[i1-1,2:(j-1)]/(2*dt*al1)
      }
      A = c(rep(1,j-2))
      A1 = c(rep(a1,j-2))
      A2 = c(rep(a2,j-2))
      bOx[j-2] = bOx[j-2] - A2[j-2]*1
      bRed[j-2] = bRed[j-2] - A2[j-2]*Red[i1,j]
      uox = c(rep(0, DerApprox))
      vox = c(rep(1, DerApprox))
      uRed = c(rep(0, DerApprox))
      vRed = c(rep(1, DerApprox))

      for (j1 in ((j-3):1)) {

        bOx[j1] = bOx[j1] - A2[j-2]*bOx[j1+1]/A1[j1+1]
        bRed[j1] = bRed[j1] - A2[j-2]*bRed[j1+1]/A1[j1+1]
        A1[j1] = A1[j1] - A2[j-2]/A1[j1+1]

      }
      for (m in 2:DerApprox) {
        uox[m] = (bOx[m-1] - uox[m-1])/A1[m-1]
        vox[m] = -vox[m-1]/A1[m-1]
        uRed[m] = (bRed[m-1] - uRed[m-1])/A1[m-1]
        vRed[m] = -vRed[m-1]/A1[m-1]
      }

      B = matrix(data = c(1,-exp(p),
                          Derv(npoints = DerApprox, CoefMat = T)[1] + sum(vox[2:DerApprox]*Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]),
                          Derv(npoints = DerApprox, CoefMat = T)[1] + sum(vRed[2:DerApprox]*Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox])),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(0, -sum(uox[2:DerApprox]*Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]) - sum(uRed[2:DerApprox]*Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]


      for (j1 in 1:(j-2)) {
        Ox[i1+1,j1+1] = (bOx[j1] -Ox[i1+1,j1])/A1[j1]
        Red[i1+1,j1+1] = (bRed[j1] -Red[i1+1,j1])/A1[j1]
      }
    }

    Jox = Derv(Ox = Ox, h = h, npoints = DerApprox)


  } else if (!(Method %in% c("Euler", "BI", "RK4", "CN", "BDF"))) {
    return("Available methods are Euler, BI, RK4, CN and BDF")
  }

  G = Jox
  i = (n*FA*G*Dx*Area*Co)/(6.4*h*sqrt(Dx))

  graphy = ggplot(data = data.frame(i[1:(length(i)-1)],t[1:(length(i)-1)]),
                  aes(y = i[1:(length(i)-1)], x = t[1:(length(i)-1)])) +
    geom_point() + xlab("Time(s)") +
    ylab("Current (A)") + theme_classic()

  if (errCheck == TRUE){
    return(list(G,Dx,Co,dt,h,l,j,i,n,Area,p,Da))
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
#' @param Co bulk concentration
#' @param Dx diffusion coefficient
#' @param Eo reduction potential of the species
#' @param Dm simulation parameter, maximum 0.5 for explicit methods
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


LinSwp = function(Co = 0.001, Dx = 0.00001, Eo = 0.1, Dm = 0.45,
                  Vi = 0.3, Vf = -0.3, Vs = 0.001, ko = 0.01,
                  alpha = 0.5, Temp = 298.15, n = 1, Area = 1,
                  DerApprox = 2, errCheck = FALSE, Method = "Euler"){


  FA = 96485
  R = 8.3145
  f = ((FA*n)/(R*Temp))
  Da = Dx/Dx
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
  Ox = OneMat(l +1, j)
  Red =ZeroMat(l +1, j)
  Jox = ZeroMat(l+1, 1)

  if (Method == "Euler") {
    for (i1 in 1:l) {
      B = matrix(data = c((Kf[i1]*h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Kb[i1]*h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]), -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2]
      for (j1 in 2:(j-1)) {
        Ox[i1+1,j1] = Ox[i1,j1] + Dm*(Ox[i1, j1-1] - 2*Ox[i1, j1] + Ox[i1, j1+1])
        Red[i1+1, j1] = Red[i1, j1] + Dm*(Red[i1, j1-1] + Red[i1, j1+1] - 2*Red[i1,j1])
      }
    }

    Jox = Derv(Ox = Ox, h = h, npoints = DerApprox)

  } else if (Method == "RK4") {

    for (i1 in 1:(l-1)) {
      k1 = ZeroMat(j)
      k2 = ZeroMat(j)
      k3 = ZeroMat(j)
      k4 = ZeroMat(j)
      k1red = ZeroMat(j)
      k2red = ZeroMat(j)
      k3red = ZeroMat(j)
      k4red = ZeroMat(j)
      B = matrix(data = c((Kf[i1]*h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Kb[i1]*h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2]

      for (j1 in 2:(j-1)) {
        k1[j1] = Dm*(Ox[i1, j1 -1] - 2*Ox[i1, j1] + Ox[i1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + k1[j1]*0.5
        k1red[j1] = Dm*(Red[i1, j1 -1] - 2*Red[i1, j1] + Red[i1, j1+1])
        Red[i1 + 1,j1] = Red[i1,j1] + k1red[j1]*0.5
      }

      B = matrix(data = c((Kf[i1]*h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Kb[i1]*h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1+1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]

      for (j1 in 2:(j-1)) {
        k2[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + k2[j1]*0.5
        k2red[j1] = Dm*(Red[i1 + 1, j1 -1] - 2*Red[i1 + 1, j1] + Red[i1 + 1, j1+1])
        Red[i1 + 1,j1] = Red[i1,j1] + k2red[j1]*0.5
      }

      B = matrix(data = c((Kf[i1]*h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Kb[i1]*h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1+1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]

      for (j1 in 2:(j-1)) {
        k3[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + k3[j1]
        k3red[j1] = Dm*(Red[i1 + 1, j1 -1] - 2*Red[i1 + 1, j1] + Red[i1 + 1, j1+1])
        Red[i1 + 1,j1] = Red[i1,j1] + k3red[j1]
      }

      B = matrix(data = c((Kf[i1]*h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Kb[i1]*h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1+1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]

      for (j1 in 2:(j-1)) {
        k4[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + (k1[j1] + 2*k2[j1] + 2*k3[j1] + k4[j1])/6
        k4red[j1] = Dm*(Red[i1 + 1, j1 -1] - 2*Red[i1 + 1, j1] + Red[i1 + 1, j1+1])
        Red[i1 + 1,j1] = Red[i1,j1] + (k1red[j1] + 2*k2red[j1] + 2*k3red[j1] + k4red[j1])/6
      }
    }

    Jox = Derv(Ox = Ox, h = h, npoints = DerApprox)

  } else if (Method == "BI") {
    al1 = 1/(h^2)
    al2 = -2/(h^2)
    al3 = 1/(h^2)
    a1 = (al2 - 1/dt)/al1
    a2 = al3/al1

    for (i1 in 1:(l-1)) {

      B = matrix(data = c((Kf[i1]*h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Kb[i1]*h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2]

      bOx = (-Ox[i1,(2:(j-1))]/(al1*dt))
      bRed = (-Red[i1,(2:(j-1))]/(al1*dt))
      A = c(rep(1,j-2))
      A1 = c(rep(a1,j-2))
      A2 = c(rep(a2,j-2))
      bOx[j-2] = bOx[j-2] - A2[j-2]*1
      bRed[j-2] = bRed[j-2] - A2[j-2]*Red[i1,j]
      uox = c(rep(0, DerApprox))
      vox = c(rep(1, DerApprox))
      uRed = c(rep(0, DerApprox))
      vRed = c(rep(1, DerApprox))

      for (j1 in ((j-3):1)) {

        bOx[j1] = bOx[j1] - A2[j-2]*bOx[j1+1]/A1[j1+1]
        bRed[j1] = bRed[j1] - A2[j-2]*bRed[j1+1]/A1[j1+1]
        A1[j1] = A1[j1] - A2[j-2]/A1[j1+1]

      }

      for (m in 2:DerApprox) {
        uox[m] = (bOx[m-1] - uox[m-1])/A1[m-1]
        vox[m] = -vox[m-1]/A1[m-1]
        uRed[m] = (bRed[m-1] - uRed[m-1])/A1[m-1]
        vRed[m] = -vRed[m-1]/A1[m-1]
      }

      B = matrix(data = c((Kf[i1]*h -  sum(vox*Derv(npoints = DerApprox, CoefMat = T))),
                          -Kb[i1]*h,
                          sum(vox*Derv(npoints = DerApprox, CoefMat = T)),
                          sum(vRed*Derv(npoints = DerApprox, CoefMat = T))),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(uox*Derv(npoints = DerApprox, CoefMat = T)),
                          -sum(uox*Derv(npoints = DerApprox, CoefMat = T)) - sum(uRed*Derv(npoints = DerApprox, CoefMat = T))))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]


      for (j1 in 1:(j-2)) {
        Ox[i1+1,j1+1] = (bOx[j1] -Ox[i1+1,j1])/A1[j1]
        Red[i1+1,j1+1] = (bRed[j1] -Red[i1+1,j1])/A1[j1]
      }
    }

    Jox = Derv(Ox = Ox, h = h, npoints = DerApprox)

  } else if (Method == "CN") {

    al1 = 1/(h^2)
    al2 = -2/(h^2)
    al3 = 1/(h^2)
    a1 = (al2 - 2/dt)/al1
    a2 = al3/al1
    a3 = (al2 + 2/dt)/al1

    for (i1 in 1:(l-1)) {

      B = matrix(data = c((Kf[i1]*h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Kb[i1]*h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2]

      bOx = -a3*Ox[i1,(2:(j-1))] - Ox[i1,(1:(j-2))] - a2*Ox[i1,(3:j)]
      bRed = -a3*Red[i1,(2:(j-1))]- Red[i1,(1:(j-2))] - a2*Red[i1,(3:j)]
      A = c(rep(1,j-2))
      A1 = c(rep(a1,j-2))
      A2 = c(rep(a2,j-2))
      bOx[j-2] = bOx[j-2] - A2[j-2]*1
      bRed[j-2] = bRed[j-2] - A2[j-2]*Red[i1,j]
      uox = c(rep(0, DerApprox))
      vox = c(rep(1, DerApprox))
      uRed = c(rep(0, DerApprox))
      vRed = c(rep(1, DerApprox))

      for (j1 in ((j-3):1)) {

        bOx[j1] = bOx[j1] - A2[j-2]*bOx[j1+1]/A1[j1+1]
        bRed[j1] = bRed[j1] - A2[j-2]*bRed[j1+1]/A1[j1+1]
        A1[j1] = A1[j1] - A2[j-2]/A1[j1+1]

      }

      for (m in 2:DerApprox) {
        uox[m] = (bOx[m-1] - uox[m-1])/A1[m-1]
        vox[m] = -vox[m-1]/A1[m-1]
        uRed[m] = (bRed[m-1] - uRed[m-1])/A1[m-1]
        vRed[m] = -vRed[m-1]/A1[m-1]
      }

      B = matrix(data = c((Kf[i1]*h -  sum(vox*Derv(npoints = DerApprox, CoefMat = T))),
                          -Kb[i1]*h,
                          sum(vox*Derv(npoints = DerApprox, CoefMat = T)),
                          sum(vRed*Derv(npoints = DerApprox, CoefMat = T))),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(uox*Derv(npoints = DerApprox, CoefMat = T)),
                          -sum(uox*Derv(npoints = DerApprox, CoefMat = T)) - sum(uRed*Derv(npoints = DerApprox, CoefMat = T))))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]


      for (j1 in 1:(j-2)) {
        Ox[i1+1,j1+1] = (bOx[j1] -Ox[i1+1,j1])/A1[j1]
        Red[i1+1,j1+1] = (bRed[j1] -Red[i1+1,j1])/A1[j1]
      }
    }

    Jox = Derv(Ox = Ox, h = h, npoints = DerApprox)


  } else if (Method == "BDF") {

    al1 = 1/(h^2)
    al2 = -2/(h^2)
    al3 = 1/(h^2)
    a1 = (al2 - 1.5/dt)/al1
    a2 = al3/al1


    for (i1 in 1:(l-1)) {

      B = matrix(data = c((Kf[i1]*h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Kb[i1]*h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2]

      if (i1 == 1) {
        bOx = -2*Ox[i1,2:(j-1)]/(dt*al1) + Ox[1,2:(j-1)]/(2*dt*al1)
        bRed = -2*Red[i1,2:(j-1)]/(dt*al1) + Red[1,2:(j-1)]/(2*dt*al1)
      } else {
        bOx = -2*Ox[i1,2:(j-1)]/(dt*al1) + Ox[i1-1,2:(j-1)]/(2*dt*al1)
        bRed = -2*Red[i1,2:(j-1)]/(dt*al1) + Red[i1-1,2:(j-1)]/(2*dt*al1)
      }
      A = c(rep(1,j-2))
      A1 = c(rep(a1,j-2))
      A2 = c(rep(a2,j-2))
      bOx[j-2] = bOx[j-2] - A2[j-2]*1
      bRed[j-2] = bRed[j-2] - A2[j-2]*Red[i1,j]
      uox = c(rep(0, DerApprox))
      vox = c(rep(1, DerApprox))
      uRed = c(rep(0, DerApprox))
      vRed = c(rep(1, DerApprox))

      for (j1 in ((j-3):1)) {

        bOx[j1] = bOx[j1] - A2[j-2]*bOx[j1+1]/A1[j1+1]
        bRed[j1] = bRed[j1] - A2[j-2]*bRed[j1+1]/A1[j1+1]
        A1[j1] = A1[j1] - A2[j-2]/A1[j1+1]

      }

      for (m in 2:DerApprox) {
        uox[m] = (bOx[m-1] - uox[m-1])/A1[m-1]
        vox[m] = -vox[m-1]/A1[m-1]
        uRed[m] = (bRed[m-1] - uRed[m-1])/A1[m-1]
        vRed[m] = -vRed[m-1]/A1[m-1]
      }

      B = matrix(data = c((Kf[i1]*h -  sum(vox*Derv(npoints = DerApprox, CoefMat = T))),
                          -Kb[i1]*h,
                          sum(vox*Derv(npoints = DerApprox, CoefMat = T)),
                          sum(vRed*Derv(npoints = DerApprox, CoefMat = T))),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(uox*Derv(npoints = DerApprox, CoefMat = T)),
                          -sum(uox*Derv(npoints = DerApprox, CoefMat = T)) - sum(uRed*Derv(npoints = DerApprox, CoefMat = T))))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]


      for (j1 in 1:(j-2)) {
        Ox[i1+1,j1+1] = (bOx[j1] -Ox[i1+1,j1])/A1[j1]
        Red[i1+1,j1+1] = (bRed[j1] -Red[i1+1,j1])/A1[j1]
      }
    }

    Jox = Derv(Ox = Ox, h = h, npoints = DerApprox)


  } else if (!(Method %in% c("Euler", "BI", "RK4", "CN", "BDF"))) {
    return("Available methods are Euler, BI, RK4, CN and BDF")
  }

  G = Jox
  i = (n*FA*G*Dx*Area*Co)/(6.4*h*sqrt(Dx))

  graphy = ggplot(data = data.frame(i[1:(length(i)-1)],PotentialScan[1:(length(i)-1)]), aes(y = i[1:(length(i)-1)], x = PotentialScan[1:(length(i)-1)])) +
    geom_point() + scale_x_continuous(trans = "reverse") +
    xlab("Potential (V)") +
    ylab("Current (A)") +
    theme_classic()

  if (errCheck == TRUE){
    return(list(G,Dx,Co,dt,h,i,l,j,n,Area,p,Da))
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
#' @param Dm simulation parameter, maximum 0.5 for explicit methods
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

CV = function(Co = 0.001, Dx = 0.00001, Eo = 0.1, Dm = 0.45,
              Vi = 0.3, Vf = -0.3, Vs = 0.001, ko = 0.01,
              alpha = 0.5, Temp = 298.15, n = 1, Area = 1,
              DerApprox = 2, errCheck = FALSE, Method = "Euler"){


  FA = 96485
  R = 8.3145
  f = ((FA*n)/(R*Temp))
  Da = Dx/Dx
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
  Ox = OneMat(l +1, j)
  Red =ZeroMat(l +1, j)
  Jox = ZeroMat(l+1, 1)

  if (Method == "Euler") {
    for (i1 in 1:l) {
      B = matrix(data = c((Kf[i1]*h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Kb[i1]*h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]), -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2]
      for (j1 in 2:(j-1)) {
        Ox[i1+1,j1] = Ox[i1,j1] + Dm*(Ox[i1, j1-1] - 2*Ox[i1, j1] + Ox[i1, j1+1])
        Red[i1+1, j1] = Red[i1, j1] + Dm*(Red[i1, j1-1] + Red[i1, j1+1] - 2*Red[i1,j1])
      }
    }

    Jox = Derv(Ox = Ox, h = h, npoints = DerApprox)

  } else if (Method == "RK4") {

    for (i1 in 1:(l-1)) {
      k1 = ZeroMat(j)
      k2 = ZeroMat(j)
      k3 = ZeroMat(j)
      k4 = ZeroMat(j)
      k1red = ZeroMat(j)
      k2red = ZeroMat(j)
      k3red = ZeroMat(j)
      k4red = ZeroMat(j)
      B = matrix(data = c((Kf[i1]*h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Kb[i1]*h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2]

      for (j1 in 2:(j-1)) {
        k1[j1] = Dm*(Ox[i1, j1 -1] - 2*Ox[i1, j1] + Ox[i1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + k1[j1]*0.5
        k1red[j1] = Dm*(Red[i1, j1 -1] - 2*Red[i1, j1] + Red[i1, j1+1])
        Red[i1 + 1,j1] = Red[i1,j1] + k1red[j1]*0.5
      }

      B = matrix(data = c((Kf[i1]*h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Kb[i1]*h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1+1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]

      for (j1 in 2:(j-1)) {
        k2[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + k2[j1]*0.5
        k2red[j1] = Dm*(Red[i1 + 1, j1 -1] - 2*Red[i1 + 1, j1] + Red[i1 + 1, j1+1])
        Red[i1 + 1,j1] = Red[i1,j1] + k2red[j1]*0.5
      }

      B = matrix(data = c((Kf[i1]*h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Kb[i1]*h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1+1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]

      for (j1 in 2:(j-1)) {
        k3[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + k3[j1]
        k3red[j1] = Dm*(Red[i1 + 1, j1 -1] - 2*Red[i1 + 1, j1] + Red[i1 + 1, j1+1])
        Red[i1 + 1,j1] = Red[i1,j1] + k3red[j1]
      }

      B = matrix(data = c((Kf[i1]*h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Kb[i1]*h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1+1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]

      for (j1 in 2:(j-1)) {
        k4[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + (k1[j1] + 2*k2[j1] + 2*k3[j1] + k4[j1])/6
        k4red[j1] = Dm*(Red[i1 + 1, j1 -1] - 2*Red[i1 + 1, j1] + Red[i1 + 1, j1+1])
        Red[i1 + 1,j1] = Red[i1,j1] + (k1red[j1] + 2*k2red[j1] + 2*k3red[j1] + k4red[j1])/6
      }
    }

    Jox = Derv(Ox = Ox, h = h, npoints = DerApprox)

  } else if (Method == "BI") {
    al1 = 1/(h^2)
    al2 = -2/(h^2)
    al3 = 1/(h^2)
    a1 = (al2 - 1/dt)/al1
    a2 = al3/al1

    for (i1 in 1:(l-1)) {

      B = matrix(data = c((Kf[i1]*h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Kb[i1]*h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2]

      bOx = (-Ox[i1,(2:(j-1))]/(al1*dt))
      bRed = (-Red[i1,(2:(j-1))]/(al1*dt))
      A = c(rep(1,j-2))
      A1 = c(rep(a1,j-2))
      A2 = c(rep(a2,j-2))
      bOx[j-2] = bOx[j-2] - A2[j-2]*1
      bRed[j-2] = bRed[j-2] - A2[j-2]*Red[i1,j]
      uox = c(rep(0, DerApprox))
      vox = c(rep(1, DerApprox))
      uRed = c(rep(0, DerApprox))
      vRed = c(rep(1, DerApprox))

      for (j1 in ((j-3):1)) {

        bOx[j1] = bOx[j1] - A2[j-2]*bOx[j1+1]/A1[j1+1]
        bRed[j1] = bRed[j1] - A2[j-2]*bRed[j1+1]/A1[j1+1]
        A1[j1] = A1[j1] - A2[j-2]/A1[j1+1]

      }

      for (m in 2:DerApprox) {
        uox[m] = (bOx[m-1] - uox[m-1])/A1[m-1]
        vox[m] = -vox[m-1]/A1[m-1]
        uRed[m] = (bRed[m-1] - uRed[m-1])/A1[m-1]
        vRed[m] = -vRed[m-1]/A1[m-1]
      }

      B = matrix(data = c((Kf[i1]*h - sum(vox*Derv(npoints = DerApprox, CoefMat = T))),
                          -Kb[i1]*h,
                          sum(vox*Derv(npoints = DerApprox, CoefMat = T)),
                          sum(vRed*Derv(npoints = DerApprox, CoefMat = T))),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(uox*Derv(npoints = DerApprox, CoefMat = T)),
                          -sum(uox*Derv(npoints = DerApprox, CoefMat = T)) - sum(uRed*Derv(npoints = DerApprox, CoefMat = T))))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]


      for (j1 in 1:(j-2)) {
        Ox[i1+1,j1+1] = (bOx[j1] -Ox[i1+1,j1])/A1[j1]
        Red[i1+1,j1+1] = (bRed[j1] -Red[i1+1,j1])/A1[j1]
      }
    }

    Jox = Derv(Ox = Ox, h = h, npoints = DerApprox)

  } else if (Method == "CN") {

    al1 = 1/(h^2)
    al2 = -2/(h^2)
    al3 = 1/(h^2)
    a1 = (al2 - 2/dt)/al1
    a2 = al3/al1
    a3 = (al2 + 2/dt)/al1

    for (i1 in 1:(l-1)) {

      B = matrix(data = c((Kf[i1]*h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Kb[i1]*h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2]

      bOx = -a3*Ox[i1,(2:(j-1))] - Ox[i1,(1:(j-2))] - a2*Ox[i1,(3:j)]
      bRed = -a3*Red[i1,(2:(j-1))]- Red[i1,(1:(j-2))] - a2*Red[i1,(3:j)]
      A = c(rep(1,j-2))
      A1 = c(rep(a1,j-2))
      A2 = c(rep(a2,j-2))
      bOx[j-2] = bOx[j-2] - A2[j-2]*1
      bRed[j-2] = bRed[j-2] - A2[j-2]*Red[i1,j]
      uox = c(rep(0, DerApprox))
      vox = c(rep(1, DerApprox))
      uRed = c(rep(0, DerApprox))
      vRed = c(rep(1, DerApprox))

      for (j1 in ((j-3):1)) {

        bOx[j1] = bOx[j1] - A2[j-2]*bOx[j1+1]/A1[j1+1]
        bRed[j1] = bRed[j1] - A2[j-2]*bRed[j1+1]/A1[j1+1]
        A1[j1] = A1[j1] - A2[j-2]/A1[j1+1]

      }

      for (m in 2:DerApprox) {
        uox[m] = (bOx[m-1] - uox[m-1])/A1[m-1]
        vox[m] = -vox[m-1]/A1[m-1]
        uRed[m] = (bRed[m-1] - uRed[m-1])/A1[m-1]
        vRed[m] = -vRed[m-1]/A1[m-1]
      }

      B = matrix(data = c((Kf[i1]*h -  sum(vox*Derv(npoints = DerApprox, CoefMat = T))),
                          -Kb[i1]*h,
                          sum(vox*Derv(npoints = DerApprox, CoefMat = T)),
                          sum(vRed*Derv(npoints = DerApprox, CoefMat = T))),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(uox*Derv(npoints = DerApprox, CoefMat = T)),
                          -sum(uox*Derv(npoints = DerApprox, CoefMat = T)) - sum(uRed*Derv(npoints = DerApprox, CoefMat = T))))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]


      for (j1 in 1:(j-2)) {
        Ox[i1+1,j1+1] = (bOx[j1] -Ox[i1+1,j1])/A1[j1]
        Red[i1+1,j1+1] = (bRed[j1] -Red[i1+1,j1])/A1[j1]
      }
    }

    Jox = Derv(Ox = Ox, h = h, npoints = DerApprox)


  } else if (Method == "BDF") {

    al1 = 1/(h^2)
    al2 = -2/(h^2)
    al3 = 1/(h^2)
    a1 = (al2 - 1.5/dt)/al1
    a2 = al3/al1

    for (i1 in 1:(l-1)) {

      B = matrix(data = c((Kf[i1]*h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Kb[i1]*h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2]

      if (i1 == 1) {
        bOx = -2*Ox[i1,2:(j-1)]/(dt*al1) + Ox[1,2:(j-1)]/(2*dt*al1)
        bRed = -2*Red[i1,2:(j-1)]/(dt*al1) + Red[1,2:(j-1)]/(2*dt*al1)
      } else {
        bOx = -2*Ox[i1,2:(j-1)]/(dt*al1) + Ox[i1-1,2:(j-1)]/(2*dt*al1)
        bRed = -2*Red[i1,2:(j-1)]/(dt*al1) + Red[i1-1,2:(j-1)]/(2*dt*al1)
      }
      A = c(rep(1,j-2))
      A1 = c(rep(a1,j-2))
      A2 = c(rep(a2,j-2))
      bOx[j-2] = bOx[j-2] - A2[j-2]*1
      bRed[j-2] = bRed[j-2] - A2[j-2]*Red[i1,j]
      uox = c(rep(0, DerApprox))
      vox = c(rep(1, DerApprox))
      uRed = c(rep(0, DerApprox))
      vRed = c(rep(1, DerApprox))

      for (j1 in ((j-3):1)) {

        bOx[j1] = bOx[j1] - A2[j-2]*bOx[j1+1]/A1[j1+1]
        bRed[j1] = bRed[j1] - A2[j-2]*bRed[j1+1]/A1[j1+1]
        A1[j1] = A1[j1] - A2[j-2]/A1[j1+1]

      }

      for (m in 2:DerApprox) {
        uox[m] = (bOx[m-1] - uox[m-1])/A1[m-1]
        vox[m] = -vox[m-1]/A1[m-1]
        uRed[m] = (bRed[m-1] - uRed[m-1])/A1[m-1]
        vRed[m] = -vRed[m-1]/A1[m-1]
      }

      B = matrix(data = c((Kf[i1]*h -  sum(vox*Derv(npoints = DerApprox, CoefMat = T))),
                          -Kb[i1]*h,
                          sum(vox*Derv(npoints = DerApprox, CoefMat = T)),
                          sum(vRed*Derv(npoints = DerApprox, CoefMat = T))),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(uox*Derv(npoints = DerApprox, CoefMat = T)),
                          -sum(uox*Derv(npoints = DerApprox, CoefMat = T)) - sum(uRed*Derv(npoints = DerApprox, CoefMat = T))))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]


      for (j1 in 1:(j-2)) {
        Ox[i1+1,j1+1] = (bOx[j1] -Ox[i1+1,j1])/A1[j1]
        Red[i1+1,j1+1] = (bRed[j1] -Red[i1+1,j1])/A1[j1]
      }
    }

    Jox = Derv(Ox = Ox, h = h, npoints = DerApprox)


  } else if (!(Method %in% c("Euler", "BI", "RK4", "CN", "BDF"))) {
    return("Available methods are Euler, BI, RK4, CN and BDF")
  }

  G = Jox
  i = (n*FA*G*Dx*Area*Co)/(6.4*h*sqrt(Dx))

  graphy = ggplot(data = data.frame(i[1:(length(i)-1)],PotentialScan[1:(length(i)-1)]), aes(y = i[1:(length(i)-1)], x = PotentialScan[1:(length(i)-1)])) +
    geom_point() + scale_x_continuous(trans = "reverse") +
    xlab("Potential (V)") +
    ylab("Current (A)") +
    theme_classic()

  if (errCheck == TRUE){
    return(list(G,Dx,Co,dt,h,i,l,j,n,Area,p,Da))
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
#' @param Dm simulation parameter, maximum 0.5 for explicit methods
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


CVEC = function(Co = 0.001, Dx = 0.00001, Eo = 0.1, Dm = 0.45,
                Vi = 0.3, Vf = -0.3, Vs = 0.001, ko = 0.01,
                kc = 0.0001,
                alpha = 0.5, Temp = 298.15, n = 1, Area = 1,
                DerApprox = 2, errCheck = FALSE, Method = "Euler"){

  FA = 96485
  R = 8.3145
  f = ((FA*n)/(R*Temp))
  Da = Dx/Dx
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
  KC = kc*dt
  KO = ko*sqrt(dt/Dx)
  Kf = KO*exp(-alpha*p)
  Kb = KO*exp((1-alpha)*p)
  Ox = OneMat(l +1, j)
  Red =ZeroMat(l +1, j)
  Jox = ZeroMat(l+1, 1)

  if (Method == "Euler") {
    for (i1 in 1:l) {
      B = matrix(data = c((Kf[i1]*h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Kb[i1]*h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]), -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2] - KC*Red[i1,1]
      for (j1 in 2:(j-1)) {
        Ox[i1+1,j1] = Ox[i1,j1] + Dm*(Ox[i1, j1-1] - 2*Ox[i1, j1] + Ox[i1, j1+1])
        Red[i1+1, j1] = Red[i1, j1] + Dm*(Red[i1, j1-1] + Red[i1, j1+1] - 2*Red[i1,j1]) - KC*dt*Red[i1,j1]
      }
    }

    Jox = Derv(Ox = Ox, h = h, npoints = DerApprox)

  } else if (Method == "RK4") {

    for (i1 in 1:(l-1)) {
      k1 = ZeroMat(j)
      k2 = ZeroMat(j)
      k3 = ZeroMat(j)
      k4 = ZeroMat(j)
      k1red = ZeroMat(j)
      k2red = ZeroMat(j)
      k3red = ZeroMat(j)
      k4red = ZeroMat(j)
      B = matrix(data = c((Kf[i1]*h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Kb[i1]*h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2]

      for (j1 in 2:(j-1)) {
        k1[j1] = Dm*(Ox[i1, j1 -1] - 2*Ox[i1, j1] + Ox[i1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + k1[j1]*0.5
        k1red[j1] = Dm*(Red[i1, j1 -1] - 2*Red[i1, j1] + Red[i1, j1+1])
        Red[i1 + 1,j1] = Red[i1,j1] + k1red[j1]*0.5
      }

      B = matrix(data = c((Kf[i1]*h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Kb[i1]*h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1+1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]

      for (j1 in 2:(j-1)) {
        k2[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + k2[j1]*0.5
        k2red[j1] = Dm*(Red[i1 + 1, j1 -1] - 2*Red[i1 + 1, j1] + Red[i1 + 1, j1+1])
        Red[i1 + 1,j1] = Red[i1,j1] + k2red[j1]*0.5
      }

      B = matrix(data = c((Kf[i1]*h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Kb[i1]*h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1+1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]

      for (j1 in 2:(j-1)) {
        k3[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + k3[j1]
        k3red[j1] = Dm*(Red[i1 + 1, j1 -1] - 2*Red[i1 + 1, j1] + Red[i1 + 1, j1+1])
        Red[i1 + 1,j1] = Red[i1,j1] + k3red[j1]
      }

      B = matrix(data = c((Kf[i1]*h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Kb[i1]*h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1+1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2] - KC*Red[i1+1,1]

      for (j1 in 2:(j-1)) {
        k4[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + (k1[j1] + 2*k2[j1] + 2*k3[j1] + k4[j1])/6
        k4red[j1] = Dm*(Red[i1 + 1, j1 -1] - 2*Red[i1 + 1, j1] + Red[i1 + 1, j1+1])
        Red[i1 + 1,j1] = Red[i1,j1] + (k1red[j1] + 2*k2red[j1] + 2*k3red[j1] + k4red[j1])/6 - KC*dt*Red[i1+1,j1]
      }
    }

    Jox = Derv(Ox = Ox, h = h, npoints = DerApprox)

  } else if (Method == "BI") {
    al1 = 1/(h^2)
    al2 = -2/(h^2)
    al3 = 1/(h^2)
    a1 = (al2 - 1/dt)/al1
    a2 = al3/al1

    for (i1 in 1:(l-1)) {

      B = matrix(data = c((Kf[i1]*h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Kb[i1]*h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2] - KC*Red[i1,1]

      bOx = (-Ox[i1,(2:(j-1))]/(al1*dt))
      bRed = (-Red[i1,(2:(j-1))]/(al1*dt)) + KC*Red[i1,2:(j-1)]/al1
      A = c(rep(1,j-2))
      A1 = c(rep(a1,j-2))
      A2 = c(rep(a2,j-2))
      bOx[j-2] = bOx[j-2] - A2[j-2]*1
      bRed[j-2] = bRed[j-2] - A2[j-2]*(Red[i1,j] - KC*Red[i1,j])
      uox = c(rep(0, DerApprox))
      vox = c(rep(1, DerApprox))
      uRed = c(rep(0, DerApprox))
      vRed = c(rep(1, DerApprox))

      for (j1 in ((j-3):1)) {

        bOx[j1] = bOx[j1] - A2[j-2]*bOx[j1+1]/A1[j1+1]
        bRed[j1] = bRed[j1] - A2[j-2]*bRed[j1+1]/A1[j1+1]
        A1[j1] = A1[j1] - A2[j-2]/A1[j1+1]

      }

      for (m in 2:DerApprox) {
        uox[m] = (bOx[m-1] - uox[m-1])/A1[m-1]
        vox[m] = -vox[m-1]/A1[m-1]
        uRed[m] = (bRed[m-1] - uRed[m-1])/A1[m-1]
        vRed[m] = -vRed[m-1]/A1[m-1]
      }

      B = matrix(data = c((Kf[i1]*h - sum(vox*Derv(npoints = DerApprox, CoefMat = T))),
                          -Kb[i1]*h,
                          sum(vox*Derv(npoints = DerApprox, CoefMat = T)),
                          sum(vRed*Derv(npoints = DerApprox, CoefMat = T))),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(uox*Derv(npoints = DerApprox, CoefMat = T)),
                          -sum(uox*Derv(npoints = DerApprox, CoefMat = T)) - sum(uRed*Derv(npoints = DerApprox, CoefMat = T))))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]


      for (j1 in 1:(j-2)) {
        Ox[i1+1,j1+1] = (bOx[j1] - Ox[i1+1,j1])/A1[j1]
        Red[i1+1,j1+1] = (bRed[j1] - Red[i1+1,j1])/A1[j1]
      }
    }

    Jox = Derv(Ox = Ox, h = h, npoints = DerApprox)

  } else if (Method == "CN") {

    al1 = 1/(h^2)
    al2 = -2/(h^2)
    al3 = 1/(h^2)
    a1 = (al2 - 2/dt)/al1
    a2 = al3/al1
    a3 = (al2 + 2/dt)/al1

    for (i1 in 1:(l-1)) {

      B = matrix(data = c((Kf[i1]*h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Kb[i1]*h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2] - KC*Red[i1,1]

      bOx = -a3*Ox[i1,(2:(j-1))] - Ox[i1,(1:(j-2))] - a2*Ox[i1,(3:j)]
      bRed = -a3*Red[i1,(2:(j-1))]- Red[i1,(1:(j-2))] - a2*Red[i1,(3:j)] + KC*Red[i1,2:(j-1)]/al1
      A = c(rep(1,j-2))
      A1 = c(rep(a1,j-2))
      A2 = c(rep(a2,j-2))
      bOx[j-2] = bOx[j-2] - A2[j-2]*Ox[i1,j]
      bRed[j-2] = bRed[j-2] - A2[j-2]*(Red[i1,j] - KC*Red[i1,j])
      uox = c(rep(0, DerApprox))
      vox = c(rep(1, DerApprox))
      uRed = c(rep(0, DerApprox))
      vRed = c(rep(1, DerApprox))

      for (j1 in ((j-3):1)) {

        bOx[j1] = bOx[j1] - A2[j-2]*bOx[j1+1]/A1[j1+1]
        bRed[j1] = bRed[j1] - A2[j-2]*bRed[j1+1]/A1[j1+1]
        A1[j1] = A1[j1] - A2[j-2]/A1[j1+1]

      }

      for (m in 2:DerApprox) {
        uox[m] = (bOx[m-1] - uox[m-1])/A1[m-1]
        vox[m] = -vox[m-1]/A1[m-1]
        uRed[m] = (bRed[m-1] - uRed[m-1])/A1[m-1]
        vRed[m] = -vRed[m-1]/A1[m-1]
      }

      B = matrix(data = c((Kf[i1]*h -  sum(vox*Derv(npoints = DerApprox, CoefMat = T))),
                          -Kb[i1]*h,
                          sum(vox*Derv(npoints = DerApprox, CoefMat = T)),
                          sum(vRed*Derv(npoints = DerApprox, CoefMat = T))),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(uox*Derv(npoints = DerApprox, CoefMat = T)),
                          -sum(uox*Derv(npoints = DerApprox, CoefMat = T)) - sum(uRed*Derv(npoints = DerApprox, CoefMat = T))))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]


      for (j1 in 1:(j-2)) {
        Ox[i1+1,j1+1] = (bOx[j1] -Ox[i1+1,j1])/A1[j1]
        Red[i1+1,j1+1] = (bRed[j1] -Red[i1+1,j1])/A1[j1]
      }
    }

    Jox = Derv(Ox = Ox, h = h, npoints = DerApprox)


  } else if (Method == "BDF") {

    al1 = 1/(h^2)
    al2 = -2/(h^2)
    al3 = 1/(h^2)
    a1 = (al2 - 1.5/dt)/al1
    a2 = al3/al1

    for (i1 in 1:(l-1)) {

      B = matrix(data = c((Kf[i1]*h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Kb[i1]*h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red[i1,1] = C[2] - KC*Red[i1,1]

      if (i1 == 1) {
        bOx = -2*Ox[i1,2:(j-1)]/(dt*al1) + Ox[1,2:(j-1)]/(2*dt*al1)
        bRed = -2*Red[i1,2:(j-1)]/(dt*al1) + Red[1,2:(j-1)]/(2*dt*al1) + KC*Red[i1,2:(j-1)]/al1
      } else {
        bOx = -2*Ox[i1,2:(j-1)]/(dt*al1) + Ox[i1-1,2:(j-1)]/(2*dt*al1)
        bRed = -2*Red[i1,2:(j-1)]/(dt*al1) + Red[i1-1,2:(j-1)]/(2*dt*al1) + KC*Red[i1,2:(j-1)]/al1
      }

      A = c(rep(1,j-2))
      A1 = c(rep(a1,j-2))
      A2 = c(rep(a2,j-2))
      bOx[j-2] = bOx[j-2] - A2[j-2]*Ox[i1,j]
      bRed[j-2] = bRed[j-2] - A2[j-2]*(Red[i1,j] - KC*Red[i1,j])
      uox = c(rep(0, DerApprox))
      vox = c(rep(1, DerApprox))
      uRed = c(rep(0, DerApprox))
      vRed = c(rep(1, DerApprox))

      for (j1 in ((j-3):1)) {

        bOx[j1] = bOx[j1] - A2[j-2]*bOx[j1+1]/A1[j1+1]
        bRed[j1] = bRed[j1] - A2[j-2]*bRed[j1+1]/A1[j1+1]
        A1[j1] = A1[j1] - A2[j-2]/A1[j1+1]

      }

      for (m in 2:DerApprox) {
        uox[m] = (bOx[m-1] - uox[m-1])/A1[m-1]
        vox[m] = -vox[m-1]/A1[m-1]
        uRed[m] = (bRed[m-1] - uRed[m-1])/A1[m-1]
        vRed[m] = -vRed[m-1]/A1[m-1]
      }

      B = matrix(data = c((Kf[i1]*h -  sum(vox*Derv(npoints = DerApprox, CoefMat = T))),
                          -Kb[i1]*h,
                          sum(vox*Derv(npoints = DerApprox, CoefMat = T)),
                          sum(vRed*Derv(npoints = DerApprox, CoefMat = T))),
                 byrow = T, nrow = 2, ncol = 2)
      Y = matrix(data = c(sum(uox*Derv(npoints = DerApprox, CoefMat = T)),
                          -sum(uox*Derv(npoints = DerApprox, CoefMat = T)) - sum(uRed*Derv(npoints = DerApprox, CoefMat = T))))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red[i1+1,1] = C[2]


      for (j1 in 1:(j-2)) {
        Ox[i1+1,j1+1] = (bOx[j1] -Ox[i1+1,j1])/A1[j1]
        Red[i1+1,j1+1] = (bRed[j1] -Red[i1+1,j1])/A1[j1]
      }
    }

    Jox = Derv(Ox = Ox, h = h, npoints = DerApprox)


  } else if (!(Method %in% c("Euler", "BI", "RK4", "CN", "BDF"))) {
    return("Available methods are Euler, BI, RK4, CN and BDF")
  }

  G = Jox
  i = (n*FA*G*Dx*Area*Co)/(6.4*h*sqrt(Dx))

  graphy = ggplot(data = data.frame(i[1:(length(i)-1)],PotentialScan[1:(length(i)-1)]), aes(y = i[1:(length(i)-1)], x = PotentialScan[1:(length(i)-1)])) +
    geom_point() + scale_x_continuous(trans = "reverse") +
    xlab("Potential (V)") +
    ylab("Current (A)") +
    theme_classic()

  if (errCheck == TRUE){
    return(list(G,Dx,Co,dt,h,i,l,j,n,Area,p,Da))
  } else {
    return(graphy)
  }
}


#' EE behaviour cyclic voltammetry simulator
#'
#' Return a graph I vs E of the electrochemical process
#'
#' @param Co bulk concentration
#' @param Dx1 diffusion coefficient of the oxidized species
#' @param Eo1 reduction potential of the first electrochemical reaction
#' @param Vi initial potential of the sweep
#' @param Vf final potential of the sweep
#' @param Vs  potential scan rate of the simulation
#' @param ko1 heterogeneous electron transfer rate constant of the first electrochemical reaction
#' @param alpha1 charge transfer coefficient of the first electrochemical reaction
#' @param Dred diffusion coefficient of the first reduced species
#' @param Dred2 diffusion coefficient of the second reduced species
#' @param Eo2 reduction potential of the second electrochemical reaction
#' @param ko2 heterogeneous electron transfer rate constant of the second electrochemical reaction
#' @param alpha2 charge transfer coefficient of the second electrochemical reaction
#' @param Dm simulation parameter, maximum 0.5 for explicit methods
#' @param Temp temperature in kelvin
#' @param n number of electrons involved in the process
#' @param Area area of the electrode
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

CVEE = function(Co = 0.001, Dx1 = 0.00001, Eo1 = 0.1,
                Vi = 0.3, Vf = -0.3, Vs = 0.001, ko1 = 0.01,
                alpha1 = 0.5, Dred = 0.00001, Dred2 = 0.00001, Eo2 = 0.05,
                ko2 = 0.01, alpha2 = 0.5, Dm = 0.45,
                Temp = 298.15, n = 1, Area = 1,
                DerApprox = 2, errCheck = FALSE, Method = "Euler") {

  FA = 96485
  R = 8.3145
  f = ((FA*n)/(R*Temp))
  DOx = Dx1/Dx1
  DRED = Dred/Dx1
  DRED2 = Dred2/Dx1
  exptime = 2*abs(Vf-Vi)/Vs
  Tmax = ceiling(abs(f*(Vf-Vi)))
  dt = (1/(f*Vs))
  l = ceiling(exptime/dt)
  j = ceiling(6*(Tmax)^0.5)
  h = sqrt(dt/Dm)
  vt = c(1:(l+1))
  t = dt*vt
  forwardScan = Vi-Vs*t[1:((l/2) +1)]
  backwardScan = Vf + Vs*t[1:((l/2) +1)]
  PotentialScan = c(forwardScan, backwardScan)
  p1 = f*(PotentialScan - Eo1)
  p2 = f*(PotentialScan - Eo2)
  KO1 = ko1*sqrt(dt/Dx1)
  KO2 = ko2*sqrt(dt/Dx1)
  Kf1 = KO1*exp(-alpha1*p1)
  Kf2 = KO2*exp(-alpha2*p2)
  Kb1 = KO1*exp((1-alpha1)*p1)
  Kb2 = KO2*exp((1-alpha2)*p2)
  Ox = OneMat(l +1, j)
  Red1 =ZeroMat(l +1, j)
  Jox = ZeroMat(l+1, 1)
  Red2 =ZeroMat(l +1, j)
  Jred1 = ZeroMat(l+1, 1)

  if (Method == "Euler") {

    for (i1 in 1:l) {
      B = matrix(data = c((Kf1[i1]*h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Kb1[i1]*h, 0, -Kf1[i1]*h, Kb1[i1]*h + Kf2[i1]*h - Derv(npoints = DerApprox, CoefMat = T)[1], -Kb2[i1]*h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 3, ncol = 3)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]),
                          DRED*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red1[i1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - DRED*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red1[i1,2:DerApprox]) - DRED2*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red2[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red1[i1,1] = C[2]
      Red2[i1,1] = C[3]

      for (j1 in 2:(j-1)) {
        Ox[i1+1,j1] = Ox[i1,j1] + Dm*(Ox[i1, j1-1] - 2*Ox[i1, j1] + Ox[i1, j1+1])
        Red1[i1+1, j1] = Red1[i1, j1] + DRED*Dm*(Red1[i1, j1-1] + Red1[i1, j1+1] - 2*Red1[i1,j1])
        Red2[i1+1, j1] = Red2[i1, j1] + DRED2*Dm*(Red2[i1, j1-1] + Red2[i1, j1+1] - 2*Red2[i1,j1])
      }
    }

    Jox = Derv(Ox = Ox, h = h, npoints = DerApprox)
    Jred1 = Derv(Ox = Red1, h = h, npoints = DerApprox)

  } else if (Method == "RK4") {

    for (i1 in 1:(l-1)) {
      k1 = ZeroMat(j)
      k2 = ZeroMat(j)
      k3 = ZeroMat(j)
      k4 = ZeroMat(j)
      k1red = ZeroMat(j)
      k2red = ZeroMat(j)
      k3red = ZeroMat(j)
      k4red = ZeroMat(j)
      k1red2 = ZeroMat(j)
      k2red2 = ZeroMat(j)
      k3red2 = ZeroMat(j)
      k4red2 = ZeroMat(j)
      B = matrix(data = c((Kf1[i1]*h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Kb1[i1]*h, 0, -Kf1[i1]*h, Kb1[i1]*h + Kf2[i1]*h - Derv(npoints = DerApprox, CoefMat = T)[1], -Kb2[i1]*h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 3, ncol = 3)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]),
                          DRED*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red1[i1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - DRED*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red1[i1,2:DerApprox]) - DRED2*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red2[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red1[i1,1] = C[2]
      Red2[i1,1] = C[3]

      for (j1 in 2:(j-1)) {
        k1[j1] = Dm*(Ox[i1, j1 -1] - 2*Ox[i1, j1] + Ox[i1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + k1[j1]*0.5
        k1red[j1] = DRED*Dm*(Red1[i1, j1 -1] - 2*Red1[i1, j1] + Red1[i1, j1+1])
        Red1[i1 + 1,j1] = Red1[i1,j1] + k1red[j1]*0.5
        k1red2[j1] = DRED2*Dm*(Red2[i1, j1 -1] - 2*Red2[i1, j1] + Red2[i1, j1+1])
        Red2[i1 + 1,j1] = Red2[i1,j1] + k1red2[j1]*0.5
      }

      B = matrix(data = c((Kf1[i1]*h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Kb1[i1]*h, 0, -Kf1[i1]*h, Kb1[i1]*h + Kf2[i1]*h - Derv(npoints = DerApprox, CoefMat = T)[1], -Kb2[i1]*h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 3, ncol = 3)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]),
                          DRED*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red1[i1+1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]) - DRED*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red1[i1+1,2:DerApprox]) - DRED2*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red2[i1+1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red1[i1+1,1] = C[2]
      Red2[i1+1,1] = C[3]

      for (j1 in 2:(j-1)) {
        k2[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + k2[j1]*0.5
        k2red[j1] = DRED*Dm*(Red1[i1, j1 -1] - 2*Red1[i1, j1] + Red1[i1, j1+1])
        Red1[i1 + 1,j1] = Red1[i1,j1] + k2red[j1]*0.5
        k2red2[j1] = DRED2*Dm*(Red2[i1, j1 -1] - 2*Red2[i1, j1] + Red2[i1, j1+1])
        Red2[i1 + 1,j1] = Red2[i1,j1] + k2red2[j1]*0.5
      }

      B = matrix(data = c((Kf1[i1]*h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Kb1[i1]*h, 0, -Kf1[i1]*h, Kb1[i1]*h + Kf2[i1]*h - Derv(npoints = DerApprox, CoefMat = T)[1], -Kb2[i1]*h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 3, ncol = 3)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]),
                          DRED*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red1[i1+1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]) - DRED*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red1[i1+1,2:DerApprox]) - DRED2*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red2[i1+1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red1[i1+1,1] = C[2]
      Red2[i1+1,1] = C[3]

      for (j1 in 2:(j-1)) {
        k3[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + k3[j1]
        k3red[j1] = DRED*Dm*(Red1[i1 + 1, j1 -1] - 2*Red1[i1 + 1, j1] + Red1[i1 + 1, j1+1])
        Red1[i1 + 1,j1] = Red1[i1,j1] + k3red[j1]
        k3red2[j1] = DRED2*Dm*(Red2[i1 + 1, j1 -1] - 2*Red2[i1 + 1, j1] + Red2[i1 + 1, j1+1])
        Red2[i1 + 1,j1] = Red2[i1,j1] + k3red2[j1]
      }

      B = matrix(data = c((Kf1[i1]*h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Kb1[i1]*h, 0, -Kf1[i1]*h, Kb1[i1]*h + Kf2[i1]*h - Derv(npoints = DerApprox, CoefMat = T)[1], -Kb2[i1]*h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 3, ncol = 3)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]),
                          DRED*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red1[i1+1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1+1,2:DerApprox]) - DRED*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red1[i1+1,2:DerApprox]) - DRED2*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red2[i1+1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red1[i1+1,1] = C[2]
      Red2[i1+1,1] = C[3]

      for (j1 in 2:(j-1)) {
        k4[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
        Ox[i1 + 1,j1] = Ox[i1,j1] + (k1[j1] + 2*k2[j1] + 2*k3[j1] + k4[j1])/6
        k4red[j1] = DRED*Dm*(Red1[i1 + 1, j1 -1] - 2*Red1[i1 + 1, j1] + Red1[i1 + 1, j1+1])
        Red1[i1 + 1,j1] = Red1[i1,j1] + (k1red[j1] + 2*k2red[j1] + 2*k3red[j1] + k4red[j1])/6
        k4red2[j1] = DRED2*Dm*(Red2[i1 + 1, j1 -1] - 2*Red2[i1 + 1, j1] + Red2[i1 + 1, j1+1])
        Red2[i1 + 1,j1] = Red2[i1,j1] + (k1red2[j1] + 2*k2red2[j1] + 2*k3red2[j1] + k4red2[j1])/6
      }
    }

    Jox = Derv(Ox = Ox, h = h, npoints = DerApprox)
    Jred1 = Derv(Ox = Red1, h = h, npoints = DerApprox)

  } else if (Method == "BI") {
    al1 = 1/(h^2)
    al2 = -2/(h^2)
    al3 = 1/(h^2)
    a1 = (al2 - 1/dt)/al1
    a2 = al3/al1
    al1red = DRED/(h^2)
    al2red = -(2*DRED)/(h^2)
    al3red = DRED/(h^2)
    a1red = (al2red - 1/dt)/al1red
    a2red = al3red/al1red
    al1red2 = DRED2/(h^2)
    al2red2 = -(2*DRED2)/(h^2)
    al3red2 = DRED2/(h^2)
    a1red2 = (al2red2 - 1/dt)/al1red2
    a2red2 = al3red2/al1red2

    for (i1 in 1:(l-1)) {

      B = matrix(data = c((Kf1[i1]*h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Kb1[i1]*h, 0, -Kf1[i1]*h, Kb1[i1]*h + Kf2[i1]*h - Derv(npoints = DerApprox, CoefMat = T)[1], -Kb2[i1]*h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 3, ncol = 3)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]),
                          DRED*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red1[i1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - DRED*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red1[i1,2:DerApprox]) - DRED2*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red2[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red1[i1,1] = C[2]
      Red2[i1,1] = C[3]

      bOx = (-Ox[i1,(2:(j-1))]/(al1*dt))
      bRed = (-Red1[i1,(2:(j-1))]/(al1red*dt))
      bRed2 = (-Red2[i1,(2:(j-1))]/(al1red2*dt))
      A = c(rep(1,j-2))
      A1 = c(rep(a1,j-2))
      A2 = c(rep(a2,j-2))
      A1red = c(rep(a1red,j-2))
      A2red = c(rep(a2red,j-2))
      A1red2 = c(rep(a1red2,j-2))
      A2red2 = c(rep(a2red2,j-2))

      bOx[j-2] = bOx[j-2] - A2[j-2]*1
      bRed[j-2] = bRed[j-2] - A2red[j-2]*Red1[i1,j]
      bRed2[j-2] = bRed2[j-2] - A2red2[j-2]*Red2[i1,j]
      uox = c(rep(0, DerApprox))
      vox = c(rep(1, DerApprox))
      uRed = c(rep(0, DerApprox))
      vRed = c(rep(1, DerApprox))
      uRed2 = c(rep(0, DerApprox))
      vRed2 = c(rep(1, DerApprox))

      for (j1 in ((j-3):1)) {

        bOx[j1] = bOx[j1] - A2[j-2]*bOx[j1+1]/A1[j1+1]
        bRed[j1] = bRed[j1] - A2red[j-2]*bRed[j1+1]/A1red[j1+1]
        bRed2[j1] = bRed2[j1] - A2red2[j-2]*bRed2[j1+1]/A1red2[j1+1]
        A1[j1] = A1[j1] - A2[j-2]/A1[j1+1]
        A1red[j1] = A1red[j1] - A2red[j-2]/A1red[j1+1]
        A1red2[j1] = A1red2[j1] - A2red2[j-2]/A1red2[j1+1]

      }

      for (m in 2:DerApprox) {

        uox[m] = (bOx[m-1] - uox[m-1])/A1[m-1]
        vox[m] = -vox[m-1]/A1[m-1]
        uRed[m] = (bRed[m-1] - uRed[m-1])/A1red[m-1]
        vRed[m] = -vRed[m-1]/A1red[m-1]
        uRed2[m] = (bRed2[m-1] - uRed2[m-1])/A1red2[m-1]
        vRed2[m] = -vRed2[m-1]/A1red2[m-1]

      }

      B = matrix(data = c((Kf1[i1]*h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Kb1[i1]*h, 0, -Kf1[i1]*h, Kb1[i1]*h + Kf2[i1]*h - sum(vRed*Derv(npoints = DerApprox, CoefMat = T)), -Kb2[i1]*h,
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


      for (j1 in 1:(j-2)) {
        Ox[i1+1,j1+1] = (bOx[j1] - Ox[i1+1,j1])/A1[j1]
        Red1[i1+1,j1+1] = (bRed[j1] - Red1[i1+1,j1])/A1red[j1]
        Red2[i1+1,j1+1] = (bRed2[j1] - Red2[i1+1,j1])/A1red2[j1]
      }
    }

    Jox = Derv(Ox = Ox, h = h, npoints = DerApprox)
    Jred1 = Derv(Ox = Red1, h = h, npoints = DerApprox)

  } else if (Method == "CN") {

    al1 = 1/(h^2)
    al2 = -2/(h^2)
    al3 = 1/(h^2)
    a1 = (al2 - 2/dt)/al1
    a2 = al3/al1
    a3 = (al2 + 2/dt)/al1
    al1red = DRED/(h^2)
    al2red = -(2*DRED)/(h^2)
    al3red = DRED/(h^2)
    a1red = (al2red - 2/dt)/al1red
    a2red = al3red/al1red
    a3red = (al2red + 2/dt)/al1red
    al1red2 = DRED2/(h^2)
    al2red2 = -(2*DRED2)/(h^2)
    al3red2 = DRED2/(h^2)
    a1red2 = (al2red2 - 2/dt)/al1red2
    a2red2 = al3red2/al1red2
    a3red2 = (al2red2 + 2/dt)/al1red2

    for (i1 in 1:(l-1)) {

      B = matrix(data = c((Kf1[i1]*h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Kb1[i1]*h, 0, -Kf1[i1]*h, Kb1[i1]*h + Kf2[i1]*h - Derv(npoints = DerApprox, CoefMat = T)[1], -Kb2[i1]*h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 3, ncol = 3)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]),
                          DRED*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red1[i1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - DRED*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red1[i1,2:DerApprox]) - DRED2*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red2[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red1[i1,1] = C[2]
      Red2[i1,1] = C[3]

      bOx = -a3*Ox[i1,(2:(j-1))] - Ox[i1,(1:(j-2))] - a2*Ox[i1,(3:j)]
      bRed = -a3red*Red1[i1,(2:(j-1))]- Red1[i1,(1:(j-2))] - a2red*Red1[i1,(3:j)]
      bRed2 = -a3red2*Red2[i1,(2:(j-1))]- Red2[i1,(1:(j-2))] - a2red2*Red2[i1,(3:j)]
      A = c(rep(1,j-2))
      A1 = c(rep(a1,j-2))
      A2 = c(rep(a2,j-2))
      A1red = c(rep(a1red,j-2))
      A2red = c(rep(a2red,j-2))
      A1red2 = c(rep(a1red2,j-2))
      A2red2 = c(rep(a2red2,j-2))

      bOx[j-2] = bOx[j-2] - A2[j-2]*1
      bRed[j-2] = bRed[j-2] - A2red[j-2]*Red1[i1,j]
      bRed2[j-2] = bRed2[j-2] - A2red2[j-2]*Red2[i1,j]
      uox = c(rep(0, DerApprox))
      vox = c(rep(1, DerApprox))
      uRed = c(rep(0, DerApprox))
      vRed = c(rep(1, DerApprox))
      uRed2 = c(rep(0, DerApprox))
      vRed2 = c(rep(1, DerApprox))

      for (j1 in ((j-3):1)) {

        bOx[j1] = bOx[j1] - A2[j-2]*bOx[j1+1]/A1[j1+1]
        bRed[j1] = bRed[j1] - A2red[j-2]*bRed[j1+1]/A1red[j1+1]
        bRed2[j1] = bRed2[j1] - A2red2[j-2]*bRed2[j1+1]/A1red2[j1+1]
        A1[j1] = A1[j1] - A2[j-2]/A1[j1+1]
        A1red[j1] = A1red[j1] - A2red[j-2]/A1red[j1+1]
        A1red2[j1] = A1red2[j1] - A2red2[j-2]/A1red2[j1+1]

      }

      for (m in 2:DerApprox) {

        uox[m] = (bOx[m-1] - uox[m-1])/A1[m-1]
        vox[m] = -vox[m-1]/A1[m-1]
        uRed[m] = (bRed[m-1] - uRed[m-1])/A1red[m-1]
        vRed[m] = -vRed[m-1]/A1red[m-1]
        uRed2[m] = (bRed2[m-1] - uRed2[m-1])/A1red2[m-1]
        vRed2[m] = -vRed2[m-1]/A1red2[m-1]

      }

      B = matrix(data = c((Kf1[i1]*h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Kb1[i1]*h, 0, -Kf1[i1]*h, Kb1[i1]*h + Kf2[i1]*h - sum(vRed*Derv(npoints = DerApprox, CoefMat = T)), -Kb2[i1]*h,
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


      for (j1 in 1:(j-2)) {
        Ox[i1+1,j1+1] = (bOx[j1] - Ox[i1+1,j1])/A1[j1]
        Red1[i1+1,j1+1] = (bRed[j1] - Red1[i1+1,j1])/A1red[j1]
        Red2[i1+1,j1+1] = (bRed2[j1] - Red2[i1+1,j1])/A1red2[j1]
      }
    }

    Jox = Derv(Ox = Ox, h = h, npoints = DerApprox)
    Jred1 = Derv(Ox = Red1, h = h, npoints = DerApprox)


  } else if (Method == "BDF") {

    al1 = 1/(h^2)
    al2 = -2/(h^2)
    al3 = 1/(h^2)
    a1 = (al2 - 1.5/dt)/al1
    a2 = al3/al1
    al1red = DRED/(h^2)
    al2red = -(2*DRED)/(h^2)
    al3red = DRED/(h^2)
    a1red = (al2red - 1.5/dt)/al1red
    a2red = al3red/al1red
    al1red2 = DRED2/(h^2)
    al2red2 = -(2*DRED2)/(h^2)
    al3red2 = DRED2/(h^2)
    a1red2 = (al2red2 - 1.5/dt)/al1red2
    a2red2 = al3red2/al1red2

    for (i1 in 1:(l-1)) {

      B = matrix(data = c((Kf1[i1]*h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Kb1[i1]*h, 0, -Kf1[i1]*h, Kb1[i1]*h + Kf2[i1]*h - Derv(npoints = DerApprox, CoefMat = T)[1], -Kb2[i1]*h,
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1],
                          Derv(npoints = DerApprox, CoefMat = T)[1]),
                 byrow = T, nrow = 3, ncol = 3)
      Y = matrix(data = c(sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]),
                          DRED*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red1[i1,2:DerApprox]),
                          -sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Ox[i1,2:DerApprox]) - DRED*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red1[i1,2:DerApprox]) - DRED2*sum(Derv(npoints = DerApprox, CoefMat = T)[2:DerApprox]*Red2[i1,2:DerApprox])))
      C = invMat(B) %*% Y
      Ox[i1,1] = C[1]
      Red1[i1,1] = C[2]
      Red2[i1,1] = C[3]

      if (i1 == 1) {
        bOx = -2*Ox[i1,2:(j-1)]/(dt*al1) + Ox[1,2:(j-1)]/(2*dt*al1)
        bRed = -2*Red1[i1,2:(j-1)]/(dt*al1red) + Red1[1,2:(j-1)]/(2*dt*al1red)
        bRed2 = -2*Red2[i1,2:(j-1)]/(dt*al1red2) + Red2[1,2:(j-1)]/(2*dt*al1red2)
      } else {
        bOx = -2*Ox[i1,2:(j-1)]/(dt*al1) + Ox[i1-1,2:(j-1)]/(2*dt*al1)
        bRed = -2*Red1[i1,2:(j-1)]/(dt*al1red) + Red1[i1-1,2:(j-1)]/(2*dt*al1red)
        bRed2 = -2*Red2[i1,2:(j-1)]/(dt*al1red2) + Red2[i1-1,2:(j-1)]/(2*dt*al1red2)
      }

      A = c(rep(1,j-2))
      A1 = c(rep(a1,j-2))
      A2 = c(rep(a2,j-2))
      A1red = c(rep(a1red,j-2))
      A2red = c(rep(a2red,j-2))
      A1red2 = c(rep(a1red2,j-2))
      A2red2 = c(rep(a2red2,j-2))

      bOx[j-2] = bOx[j-2] - A2[j-2]*1
      bRed[j-2] = bRed[j-2] - A2red[j-2]*Red1[i1,j]
      bRed2[j-2] = bRed2[j-2] - A2red2[j-2]*Red2[i1,j]
      uox = c(rep(0, DerApprox))
      vox = c(rep(1, DerApprox))
      uRed = c(rep(0, DerApprox))
      vRed = c(rep(1, DerApprox))
      uRed2 = c(rep(0, DerApprox))
      vRed2 = c(rep(1, DerApprox))

      for (j1 in ((j-3):1)) {

        bOx[j1] = bOx[j1] - A2[j-2]*bOx[j1+1]/A1[j1+1]
        bRed[j1] = bRed[j1] - A2red[j-2]*bRed[j1+1]/A1red[j1+1]
        bRed2[j1] = bRed2[j1] - A2red2[j-2]*bRed2[j1+1]/A1red2[j1+1]
        A1[j1] = A1[j1] - A2[j-2]/A1[j1+1]
        A1red[j1] = A1red[j1] - A2red[j-2]/A1red[j1+1]
        A1red2[j1] = A1red2[j1] - A2red2[j-2]/A1red2[j1+1]

      }

      for (m in 2:DerApprox) {

        uox[m] = (bOx[m-1] - uox[m-1])/A1[m-1]
        vox[m] = -vox[m-1]/A1[m-1]
        uRed[m] = (bRed[m-1] - uRed[m-1])/A1red[m-1]
        vRed[m] = -vRed[m-1]/A1red[m-1]
        uRed2[m] = (bRed2[m-1] - uRed2[m-1])/A1red2[m-1]
        vRed2[m] = -vRed2[m-1]/A1red2[m-1]

      }

      B = matrix(data = c((Kf1[i1]*h - Derv(npoints = DerApprox, CoefMat = T)[1]),
                          -Kb1[i1]*h, 0, -Kf1[i1]*h, Kb1[i1]*h + Kf2[i1]*h - sum(vRed*Derv(npoints = DerApprox, CoefMat = T)), -Kb2[i1]*h,
                          sum(vox*Derv(npoints = DerApprox, CoefMat = T)),
                          sum(vRed*Derv(npoints = DerApprox, CoefMat = T)),
                          sum(vRed2*Derv(npoints = DerApprox, CoefMat = T))),
                 byrow = T, nrow = 3, ncol = 3)
      Y = matrix(data = c(sum(uox*Derv(npoints = DerApprox, CoefMat = T)),
                          sum(uRed*Derv(npoints = DerApprox, CoefMat = T)),
                          -sum(uox*Derv(npoints = DerApprox, CoefMat = T)) - sum(uRed*Derv(npoints = DerApprox, CoefMat = T)) - sum(uRed2*Derv(npoints = DerApprox, CoefMat = T))))
      C = invMat(B) %*% Y
      Ox[i1+1,1] = C[1]
      Red1[i1+1,1] = C[2]
      Red2[i1+1,1] = C[3]


      for (j1 in 1:(j-2)) {
        Ox[i1+1,j1+1] = (bOx[j1] - Ox[i1+1,j1])/A1[j1]
        Red1[i1+1,j1+1] = (bRed[j1] - Red1[i1+1,j1])/A1red[j1]
        Red2[i1+1,j1+1] = (bRed2[j1] - Red2[i1+1,j1])/A1red2[j1]
      }
    }

    Jox = Derv(Ox = Ox, h = h, npoints = DerApprox)
    Jred1 = Derv(Ox = Red1, h = h, npoints = DerApprox)

  } else if (!(Method %in% c("Euler", "BI", "RK4", "CN", "BDF"))) {
    return("Available methods are Euler, BI, RK4, CN and BDF")
  }


  G1 = Jox
  G2 = Jox + Jred1
  i = (n*FA*(G1+G2)*Dx1*Area*Co)/(6.4*h*sqrt(Dx1))

  graphy = ggplot(data = data.frame(i[1:(length(i)-1)],PotentialScan[1:(length(i)-1)]), aes(y = i[1:(length(i)-1)], x = PotentialScan[1:(length(i)-1)])) +
    geom_point() + scale_x_continuous(trans = "reverse") +
    xlab("Potential (V)") +
    ylab("Current (A)") +
    theme_classic()

  if (errCheck == TRUE){
    return(list((G1+G2),Dx1,Dred,Dred2,Co,dt,h,i,l,j,n,Area,p1,p2,DOx,DRED,DRED2))
  } else {
    return(graphy)
  }
}

