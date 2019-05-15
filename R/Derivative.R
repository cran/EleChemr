#' Derivative calculation of concentration profile
#'
#' Return a the derivative of the concentration profile simulated
#'
#' @param npoints number of points to be used for the derivative
#' @param h space for the finite difference
#' @param Ox data upon the derivative is calculated
#' @param mode "Forward" or "Backward" the derivative will be calculated for the npoints
#' @param Derivative "First" or "Second" derivative to calculate
#' @param CoefMat if T return the derivative coefficient matrix for selected derivative
#'
#'
#' @return a vector with the derivative requested or the coefficient of such derivative
#'
#' @examples
#' Derv(npoints = 2, h = 0.13, Ox = matrix(c(1,2), nrow = 1), mode = "Forward", Derivative = "First")
#'
#' @export


Derv = function(npoints = 2, h, Ox, mode = "Forward", Derivative = "First", CoefMat = FALSE) {

  if (mode == "Forward") {

    counter = 1
    Taylor = ZeroMat(npoints-1, npoints-1)
    for (i1 in 1:(npoints-1)){
      for (j1 in 1:(npoints-1)) {
        Taylor[i1,j1] = (1/factorial(j1))*counter^(j1)
      }
      counter = counter + 1
    }

    Ainv = invMat(Taylor)

    if (CoefMat == T) {
      if (mode == "Forward" & Derivative == "First") {
        Coeff = c(-sum(Ainv[1,]), Ainv[1,1:npoints-1])
        return(Coeff)

      } else if (mode == "Forward" & Derivative == "Second") {
        Coeff = c(-sum(Ainv[2,]), Ainv[2,1:npoints-1])
        return(Coeff)
      }
    }

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
    counter = 1
    Taylor = ZeroMat(npoints-1, npoints-1)
    for (i1 in 1:(npoints-1)){
      for (j1 in 1:(npoints-1)) {
        Taylor[i1,j1] = ((-1)^(j1))*(1/factorial(j1))*counter^(j1)
      }
      counter = counter + 1
    }

    Ainv = invMat(Taylor)

    if (mode == "Backward" & Derivative == "First") {
      Coeff = c(-1*Ainv[1,(npoints-1):1], sum(Ainv[1,]))
      return(Coeff)

    } else if (mode == "Backward" & Derivative == "Second") {
      Coeff = c(-1*Ainv[2,(npoints-1):1], sum(Ainv[2,]))
      return(Coeff)
    }

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
  }
}

#' Starting Matrix of oxidazed species
#'
#' Return a matrix ixj filled with 1 value
#'
#' @param i number of rows
#' @param j number of columns
#'
#'
#' @return a matrix of dimention ixj filled with 1 value
#'
#' @examples
#' OneMat(2,2)
#'
#' @export

OneMat = function(i,j=i) {
  Y = matrix(data = c(rep(1, i*j)), nrow = i, ncol = j)
  return(Y)
}


#' Starting Matrix of reduces species and fluxes
#'
#' Return a matrix ixj filled with 0 value
#'
#' @param i number of rows
#' @param j number of columns
#'
#'
#' @return a matrix of dimention ixj filled with 1 value
#'
#' @examples
#' ZeroMat(2,2)
#'
#' @export

ZeroMat = function(i,j=i){
  Y = matrix(data = c(rep(0, i*j)), nrow = i, ncol = j)
  return(Y)
}


#' Inverse matrix
#'
#' Returns the inverse matrix of the selected one
#'
#' @param A matrix to be inverted
#'
#'
#' @return inverse matrix of the selected
#'
#' @examples
#' invMat(A = matrix(c(1,2,6,14), nrow = 2))
#'
#' @export

invMat = function(A) {
  if (is.matrix(A) == FALSE) {
    return("Pointed object is not a matrix")
  } else if (det(A) == 0) {
    return("Inversion can't be performed on a singular matrix")
  }
  INV = ZeroMat(nrow(A), ncol(A))

  for (i1 in 1:nrow(A)) {
    for (j1 in 1:ncol(A)) {
      B = A[-i1,-j1]
      if (is.matrix(B) == FALSE) {
        INV[j1,i1] = (-1)^(i1+j1)*B
      } else {
        INV[j1,i1] = (-1)^(i1+j1)*det(B)
      }
    }
  }
  INV = (1/det(A))*INV
  return(INV)
}

#' Parameters call
#'
#' Returns a list with the parameters necessary for the simulation
#'
#' @param Fun Name of the function this function is called to. Must be a string.
#' @param n. Number of electrons
#' @param Temp. Temperature for the simulation
#' @param Dx1. Diffusion coefficient of species One
#' @param eta. OverPotential for potential step
#' @param exptime. experimental time for the simulation
#' @param Eo1. reduction potential of the first electrochemical reaction
#' @param Eo2. reduction potential of the second electrochemical reaction
#' @param Eo3. reduction potential of the third electrochemical reaction
#' @param Eo4. reduction potential of the fourth electrochemical reaction
#' @param Dred1. diffusion coefficient of the first reduced species
#' @param Dred2. diffusion coefficient of the second reduced species
#' @param Dred3. diffusion coefficient of the third reduced species
#' @param Dred4. diffusion coefficient of the fourth reduced species
#' @param alpha1. charge transfer coefficient of the first electrochemical reaction
#' @param alpha2. charge transfer coefficient of the second electrochemical reaction
#' @param alpha3. charge transfer coefficient of the third electrochemical reaction
#' @param alpha4. charge transfer coefficient of the fourth electrochemical reaction
#' @param kc. Chemical rate constant for first Ox Species, used in simulation with just one species
#' @param kco. Chemical rate constant for first Ox Species
#' @param kc1. Chemical rate constant for first Red Species
#' @param ko1. heterogeneous electron transfer rate constant of the first electrochemical reaction
#' @param ko2. heterogeneous electron transfer rate constant of the second electrochemical reaction
#' @param ko3. heterogeneous electron transfer rate constant of the third electrochemical reaction
#' @param ko4. heterogeneous electron transfer rate constant of the fourth electrochemical reaction
#' @param kc2. Chemical rate constant for second Red Species
#' @param kc3. Chemical rate constant for third Red Species
#' @param kc4. Chemical rate constant for fourth Red Species
#' @param Dm. Simulation parameter, maximum 0.5 for explicit methods
#' @param Vi. Initial potential of the sweep
#' @param Vf. Final potential of the sweep
#' @param Vs. Scan rate of the simulation
#'
#'
#' @return inverse matrix of the selected
#'
#' @examples
#' ParCall("ChronAmp", n. = 1, Temp. = 298, Dx1. = 0.0001, exptime. = 1, Dm. = 0.45)
#'
#' @export


ParCall = function(Fun, n., Temp., Dx1.,
                   eta., exptime., Eo1., ko1., ko2., kc.,
                   Dm., Vf., Vi., Vs., alpha1., Eo2., Dred1., Dred2.,
                   alpha2., Dred3., Dred4., ko3., ko4., kco., kc1., kc2.,
                   kc3., kc4., alpha3., alpha4., Eo3., Eo4.){
  if (!(Fun %in% c("ChronAmp", "PotStep", "LinSwp", "CV", "CVEC", "CVEE", "Gen_CV" )) ) {
    return("Not suitable function was called for parameter calculation")
  }
  if (Fun == "ChronAmp") {
    FA = 96485
    R = 8.3145
    f = ((FA*n.)/(R*Temp.))
    Da = Dx1./Dx1.
    l = 100
    tau = exptime.
    dt = exptime./l
    dtn = dt/tau
    h = sqrt((Da*dtn)/Dm.)
    j = ceiling(6*(l)^0.5)
    vt = c(1:l)
    t = dt*vt
    Par = list(FA,R,f,dtn,Da,l,h,j,t,tau)
    names(Par) = c("FA", "R", "f", "dtn", "Da", "l", "h", "j", "t", "tau")
    return(Par)

  } else if (Fun == "PotStep") {

    FA = 96485
    R = 8.3145
    f = ((FA*n.)/(R*Temp.))
    Da = Dx1./Dx1.
    l = 100
    tau = exptime.
    dt = exptime./l
    dtn = dt/tau
    p = eta.*f
    h = sqrt((Da*dtn)/Dm.)
    j = ceiling(6*(l)^0.5)
    vt = c(1:l)
    t = dt*vt
    Par = list(FA,R,f,dtn,p,Da,l,h,j,t,tau)
    names(Par) = c("FA", "R", "f", "dtn", "p", "Da", "l", "h", "j", "t","tau")
    return(Par)

  } else if (Fun == "LinSwp") {

    FA = 96485
    R = 8.3145
    f = ((FA*n.)/(R*Temp.))
    Da = Dx1./Dx1.
    exptime = abs(Vf.-Vi.)/Vs.
    l = 100
    dt = exptime/l
    tau = (1/(f*Vs.))
    dtn = dt/tau
    j = ceiling(6*(l)^0.5)
    h = sqrt((Da*dtn)/Dm.)
    vt = c(1:(l+1))
    t = dt*vt
    PotentialScan = Vi.-Vs.*t
    p = f*(Vi.- Eo1.) - t/tau
    KO = ko1.*sqrt(tau/Dx1.)
    Kf = KO*exp(-alpha1.*p)
    Kb = KO*exp((1-alpha1.)*p)
    Par = list(FA,R,f,dtn,Da,l,h,j,t,PotentialScan,KO,Kf,Kb,tau)
    names(Par) = c("FA", "R", "f", "dtn", "Da", "l", "h",
                   "j", "t", "PotentialScan", "KO", "Kf", "Kb", "tau")
    return(Par)

  } else if (Fun == "CV" | Fun == "CVEC") {

    FA = 96485
    R = 8.3145
    f = ((FA*n.)/(R*Temp.))
    Da = Dx1./Dx1.
    exptime = 2*abs(Vf.-Vi.)/Vs.
    l = 100
    dt = exptime/l
    tau = (1/(f*Vs.))
    dtn = dt/tau
    j = ceiling(6*(l)^0.5)
    h = sqrt((Da*dtn)/Dm.)
    vt = c(1:(l+1))
    t = dt*vt
    forwardScan = Vi.-Vs.*t[1:((l/2) +1)]
    backwardScan = Vf. + Vs.*t[1:((l/2) +1)]
    PotentialScan = c(forwardScan, backwardScan)
    pf = f*(Vi.- Eo1.) - t[1:((l/2) +1)]/tau
    pb = f*(Vf.- Eo1.) + t[1:((l/2) +1)]/tau
    p = c(pf,pb)
    KO = ko1.*sqrt(tau/Dx1.)
    Kf = KO*exp(-alpha1.*p)
    Kb = KO*exp((1-alpha1.)*p)

    if (Fun == "CVEC") {

      KC = kc.*tau
      Par = list(FA,R,f,dtn,Da,l,h,j,t,PotentialScan,KO,Kf,Kb,KC,tau)
      names(Par) = c("FA", "R", "f", "dtn", "Da", "l", "h",
                     "j", "t", "PotentialScan", "KO", "Kf", "Kb", "KC","tau")
      return(Par)

    }

    Par = list(FA,R,f,dtn,Da,l,h,j,t,PotentialScan,KO,Kf,Kb,tau)
    names(Par) = c("FA", "R", "f", "dtn", "Da", "l", "h",
                   "j", "t", "PotentialScan", "KO", "Kf", "Kb","tau")
    return(Par)
  } else if (Fun == "CVEE"){

    FA = 96485
    R = 8.3145
    f = ((FA*n.)/(R*Temp.))
    DOx = Dx1./Dx1.
    DRED = Dred1./Dx1.
    DRED2 = Dred2./Dx1.
    exptime = 2*abs(Vf.-Vi.)/Vs.
    l = 100
    dt = exptime/l
    tau = (1/(f*Vs.))
    dtn = dt/tau
    j = ceiling(6*(l)^0.5)
    h = sqrt((DOx*dtn)/Dm.)
    vt = c(1:(l+1))
    t = dt*vt
    forwardScan = Vi.-Vs.*t[1:((l/2) +1)]
    backwardScan = Vf. + Vs.*t[1:((l/2) +1)]
    PotentialScan = c(forwardScan, backwardScan)
    p1f = f*(Vi.- Eo1.) - t[1:((l/2) +1)]/tau
    p1b = f*(Vf.- Eo1.) + t[1:((l/2) +1)]/tau
    p2f = f*(Vi.- Eo2.) - t[1:((l/2) +1)]/tau
    p2b = f*(Vf.- Eo2.) + t[1:((l/2) +1)]/tau
    p1 = c(p1f,p1b)
    p2 = c(p2f,p2b)
    KO1 = ko1.*sqrt(tau/Dx1.)
    KO2 = ko2.*sqrt(tau/Dx1.)
    Kf1 = KO1*exp(-alpha1.*p1)
    Kf2 = KO2*exp(-alpha2.*p2)
    Kb1 = KO1*exp((1-alpha1.)*p1)
    Kb2 = KO2*exp((1-alpha2.)*p2)

    Par = list(FA,R,f,dtn,DOx,DRED,DRED2,l,h,j,t,PotentialScan,
               KO1,KO2,Kf1,Kf2,Kb1,Kb2,tau)
    names(Par) = c("FA", "R", "f", "dtn", "DOx", "DRED", "DRED2", "l", "h",
                   "j", "t", "PotentialScan", "KO1",
                   "KO2", "Kf1", "Kf2", "Kb1", "Kb2","tau")
    return(Par)

  } else if (Fun == "Gen_CV"){

    FA = 96485
    R = 8.3145
    f = ((FA*n.)/(R*Temp.))
    DOx = Dx1./Dx1.
    DRED = Dred1./Dx1.
    DRED2 = Dred2./Dx1.
    DRED3 = Dred3./Dx1.
    DRED4 = Dred4./Dx1.
    exptime = 2*abs(Vf.-Vi.)/Vs.
    l = 100
    dt = exptime/l
    tau = (1/(f*Vs.))
    dtn = dt/tau
    j = ceiling(6*(l)^0.5)
    h = sqrt((DOx*dtn)/Dm.)
    vt = c(1:(l+1))
    t = dt*vt
    forwardScan = Vi.-Vs.*t[1:((l/2) +1)]
    backwardScan = Vf. + Vs.*t[1:((l/2) +1)]
    PotentialScan = c(forwardScan, backwardScan)
    p1f = f*(Vi.- Eo1.) - t[1:((l/2) +1)]/tau
    p1b = f*(Vf.- Eo1.) + t[1:((l/2) +1)]/tau
    p2f = f*(Vi.- Eo2.) - t[1:((l/2) +1)]/tau
    p2b = f*(Vf.- Eo2.) + t[1:((l/2) +1)]/tau
    p3f = f*(Vi.- Eo3.) - t[1:((l/2) +1)]/tau
    p3b = f*(Vf.- Eo3.) + t[1:((l/2) +1)]/tau
    p4f = f*(Vi.- Eo4.) - t[1:((l/2) +1)]/tau
    p4b = f*(Vf.- Eo4.) + t[1:((l/2) +1)]/tau
    p1 = c(p1f,p1b)
    p2 = c(p2f,p2b)
    p3 = c(p3f,p3b)
    p4 = c(p4f,p4b)
    KO1 = ko1.*sqrt(tau/Dx1.)
    KO2 = ko2.*sqrt(tau/Dx1.)
    KO3 = ko3.*sqrt(tau/Dx1.)
    KO4 = ko4.*sqrt(tau/Dx1.)
    Kf1 = KO1*exp(-alpha1.*p1)
    Kf2 = KO2*exp(-alpha2.*p2)
    Kf3 = KO3*exp(-alpha3.*p3)
    Kf4 = KO4*exp(-alpha4.*p4)
    Kb1 = KO1*exp((1-alpha1.)*p1)
    Kb2 = KO2*exp((1-alpha2.)*p2)
    Kb3 = KO3*exp((1-alpha3.)*p3)
    Kb4 = KO4*exp((1-alpha4.)*p4)
    KCo = kco.*tau
    KC1 = kc1.*tau
    KC2 = kc2.*tau
    KC3 = kc3.*tau
    KC4 = kc4.*tau

    Par = list(FA,R,f,dtn,DOx,DRED,DRED2,DRED3,DRED4,l,h,j,t,PotentialScan,
               KO1,KO2,KO3,KO4,Kf1,Kf2,Kf3,Kf4,
               Kb1,Kb2,Kb3,Kb4,KCo,KC1,KC2,KC3,KC4,tau)
    names(Par) = c("FA", "R", "f", "dtn", "DOx", "DRED", "DRED2", "DRED3",
                   "DRED4", "l", "h",
                   "j", "t", "PotentialScan", "KO1",
                   "KO2", "KO3", "KO4", "Kf1", "Kf2",
                   "Kf3", "Kf4", "Kb1", "Kb2", "Kb3", "Kb4",
                   "KCo", "KC1", "KC2", "KC3", "KC4", "tau")
    return(Par)
  }

}

