# -*- coding: utf-8 -*-
"""
  This code was supplied in Fortran 77 format by
    Prof. Laura Toran
    Dept. Earth and Environmental Science
    Temple University, College of Science and Technology
    Buery Hall
    1901 North 13th Street
    Philadelphia, PA 19122 USA

  The original code was included as a part of the program CRAFIT
  published in:
    Toran, L. (2000), 'CRAFIT: A computer program for calibrating
    breakthrough curves of CRAFLUSH, a one-dimensional fracture
    flow and transport model'.
    Groundwater, 38: 430-434.
    https://doi.org/10.1111/j.1745-6584.2000.tb00229.x

  CRAFLUSH was converted to Fortran 90 and Python
  by Ignasi Puigdomenech in April 2024

  From the comments written by Prof. Laura Toran:
  The CRAFLUSH version is CRAVS3 dated 4/22/92.
  CRAFLUSH was provided by Ed Sudicky of the Univ of Waterloo.
  Updates be obtained from the University of Waterloo.

  ------------------------------------------------------------

         **************************************
         *                                    *
         *              CRAFLUSH              *
         *                                    *
         **************************************

               Solved and programmed by:
                     E.A. Sudicky
        Waterloo Centre for Goundwater Research
                University of Waterloo
              Waterloo, Ontario, N2L 3G1

              Copyright 1988 E.A. Sudicky

                        WARNING

         Duplication of this program or any part
         thereof without the express permission
         of E.A. Sudicky of the institute for
         Groundwater Research is strictly forbidden

         DISCLAIMER: While care has been exercised to ensure
                     the accuracy of the coding, the author is
                     not responsible for inadvertent errors.


         Transport in a system of parallel fractures
         with matrix diffusion  D > 0 solution
         This version allows the simulation of the flushing
         of the initially contaminated matrix (and fractures)
         with a solution of lower (or higher) concentration
         entering the fractures.

         This version uses the CRUMP algorithm to numerically
         invert the Laplace transformed solution.
"""
import cmath
import math
from typing import List

def CraFlush(C0: float, Cim: float, v: float, alpha: float, b: float, sep: float,
             theta: float, tau: float, D0: float, R: float, Rp: float,
             THalf: float, X: float, Z: float, T: float,
             Ns: int, Ts: List[float], Css: List[float],
             AlphR: float, relErr: float, N: int):
    """
     Definition of variables to enter in datafile.

     NOTE: Use consistent units for input data

         C0= Source concentration at fracture origin
             NOTE:  When using variable source let C0=0
         Cim=Initial concentration in matrix (and in fractures)
         v= Velocity in fracture
         alpha=Longitudinal dispersivity in fracture
         b= Fracture aperture
         sep= Fracture spacing
         theta= Matrix porosity
             NOTE: theta=0 is alllowed, in which case a solution
             with no matrix diffusion will be obtained
         tau= Matrix tortuosity
         D0= Diffusion coeff in water
         R= Fracture retardation factor
         Rp= Matrix retardation factor
         THalf= Half-life. If you don't want decay, then set THalf=0
         X= Distance along fracture
         Z= Distance into matrix (B/2 < Z <= SEP/2)
         T= Time
         Ns= Number of points used to define variable source
         Ts[]= Times used to define variable source
         Css[]= Concentrations used in variable source at Ts
     Numerical inversion parameters used in 'Crump':
         AlphR= alpha
         relErr= relative error
         N= decides the number of terms in the EPAL (EPsilon ALgorithm) summation
         ...Recommend: AlphR=0, relErr=1.E-6, N >= 11
         ...where 2N+1 is the number of terms in the Fourier series for numerical
         ...inversion of the Laplace transform.
         ...NOTE: If the solution oscillates about the steady state concentration,
         ...try relErr=1.E-5 or 1.E-4 etc., or increase N.
    """
    if(C0 < 0):
        print('Error: C0 = '+str(C0)+' must be >=0')
        raise SystemExit
    if(Cim < 0):
        print('Error: Cim = '+str(Cim)+' must be >=0')
        raise SystemExit
    if(v <= 0):
        print('Error: v = '+str(v)+' must be >0')
        raise SystemExit
    if(alpha < 0):
        print('Error: alpha = '+str(alpha)+' must be >=0')
        raise SystemExit
    if(b <= 0):
        print('Error: b = '+str(b)+' must be >0')
        raise SystemExit
    if(sep <= 0):
        print('Error: sep = '+str(sep)+' must be >0')
        raise SystemExit
    if(theta < 0):
        print('Error: sep = '+str(theta)+' must be >=0')
        raise SystemExit
    if(tau <= 0):
        print('Error: tau = '+str(tau)+' must be >0')
        raise SystemExit
    if(D0 < 0):
        print('Error: D0 = '+str(D0)+' must be >0')
        raise SystemExit
    if(THalf < 0):
        print('Error: THalf = '+str(THalf)+' must be >=0')
        raise SystemExit
    if(X < 0):
        print('Error: X = '+str(X)+' must be >=0')
        raise SystemExit
    if(T < 0):
        print('Error: T = '+str(T)+' must be >=0')
        raise SystemExit
    if(Ns <= 0):
        print('Error: Ns = '+str(Ns)+' must be >0')
        raise SystemExit
    for i in range(0,Ns):
        if(Ts[i] < 0):
            print('Error: Ts['+str(i)+'] = '+str(Ts[i])+' must be >=0')
            raise SystemExit
        if(Css[i] < 0):
            print('Error: Css['+str(i)+'] = '+str(Css[i])+' must be >=0')
            raise SystemExit
    if(relErr <= 0):
        print('Error: relErr = '+str(relErr)+' must be >0')
        raise SystemExit
    if(N <= 10):
        print('Error: N = '+str(N)+' must be >10')
        raise SystemExit
    b = b/2.
    sep = sep/2.
    if(Z < b or Z > sep):
        print('Error: Z = '+str(Z)+' must be >='+str(b)+' (1/2 fracture aperture) and <='+str(sep)+' (1/2 fract. spacing)')
        raise SystemExit

    D_p = D0 * tau
    D = alpha * v + D0
    #...Calculate parameters used in the Laplace transformed solution L[C]
    #...They should be stored as complex values
    ZA = complex(theta*math.sqrt(Rp*D_p)/(b*R), 0.)
    Zkappa = complex(4.*R*D/v/v, 0.)
    RpDp = math.sqrt(Rp/D_p)
    Zsig = complex(RpDp*(sep-b), 0.)
    lamda = 0.
    if(THalf >= 1.e-10): lamda = math.log(2.) / THalf
    Zdecay = complex(lamda, 0.)
    vv = v/2./D
    # ...Initialization for CRUMP
    Tfact, Aterm = initC(T,relErr,AlphR)
    vvX = vv * X
    Znu = complex(vvX, 0.)
    if(X < 1.e-10):
        if(Z <= b):
            C = C0
            Cs = C0
    else:
        RpDpZ = complex(RpDp*(sep-Z), 0.)
        # Steady-state solution given by:
        # Limit as P --> 0 of P*L[C]
        #   or C0 if no decay
        if(THalf >= 1.e-10):
            P = complex(0., 0.)
            Cs = FS(P,Znu,Zkappa,ZA,Zsig,Zdecay,RpDpZ,C0,Cim,Cs,Ts,Ns).real
        else:
            Cs = C0
        # Transient solution
        FT = Crump(T,Aterm,Tfact,Znu,Zkappa,ZA,Zsig,Zdecay,RpDpZ,N,C0,Cim,Ns,Css,Ts)
        C = FT
    return C, Cs  # return a tuple

def Crump(T: float, Aterm: float, Tfact: float, Znu: complex,
          Zkappa: complex, ZA: complex, Zsig: complex, Zdecay: complex,
          RPDPZ: complex,
          N: int, C0: float, Cim: float,
          Ns: int, CSS: List[float], Ts: List[float]):
    """
 Routine for numerical inversion of Laplace
 transformations using the CRUMP algorithm.
 This uses the 'EPAL' (EPsilon ALgorithm) method
 of series summation as suggested in the CRUMP paper.

 Reference:  Kenny S. Crump, 'Numerical
   inversion of Laplace transforms using a
   Fourier series approximation', Journal of
   the Association for Computing Machinery
   Vol 23, No. 1, Jan 1976, pp. 89-96.
   https://doi.org/10.1145/321921.321931

 Programmed by:  Frank Letniowski
   Institute for Groundwater Research
   University of Waterloo
   May, 1987

 The CRUMP subroutine numerically inverts a Laplace
 transformed function FS at a particular time T
 returning the value in FT.  FS must be a complex
 function of a complex argument.  Here, X is the Laplace
 variable 'P' which is complex.  For maximum efficiency,
 you may want to store the FS values for each X in an
 array (i.e. precompute them, along with the needed
 X-values) outside your time loop and then substitute
 the FS array into the CRUMP Fourier Series that involves
 time.  That is, it is not necessary to compute FS(X)
 for each value of the time 'T' since X depends only
 on TFACT through TMAX.

 The parameters A and TFACT should be as suggested
 above or in the paper.  N is related to the number
 of terms used in the summation i.e. 2N+1 terms are
 calculated in the 'EPAL' method of summation.
    
 Dimensioning hints:
            F(MAXTRM),Q(MAXQ)
 where: MAXTRM = 2N+1
        MAXQ   = 4N+4
    """
    M: int = 2*N + 1

    X = complex(Aterm, 0.)
    F0: float = FS(X,Znu,Zkappa,ZA,Zsig,Zdecay,RPDPZ,C0,Cim,CSS,Ts,Ns).real

    # Calculate the function values needed
    F = [complex(0.,0.)] * M  # create an 'empty' list of size 'M'; M = 2*N+1
    for k in range(0, M):
        X = complex(Aterm, float(k+1)*math.pi/Tfact)
        F[k] = FS(X,Znu,Zkappa,ZA,Zsig,Zdecay,RPDPZ,C0,Cim,CSS,Ts,Ns)
    # Initialise the queue needed for calculation of the SUM.
    # Mx is the size of the queue needed.
    # The queue stores EPSILON values needed in the
    # EPAL method to calculate the series.
    Mx: int = 2*M + 2
    cur, fin = initQ()    # set  cur=0  and  fin=0
    Q = [0.] * (4*N+4)    # create an 'empty' list of size MaxQ = (4*N+4)
    # Calculate the series by use of the EPAL method.
    # SER(i,F[i],TEMP) calculates the Ith term in the series.
    temp = math.pi * T / Tfact
    Sum = SER(F[0], temp)
    # Insert into the queue
    Q, fin = insert(Sum,Q,cur,fin,Mx-1)
    for k in range(1, M):
        Sum = Sum + SER(F[k], float(k+1)*temp)
        #  SUM is the Kth partial sum of the series
        #  Calculate EPS sub K sup (1) which is TOT
        #  upon leaving the 'j' loop
        Tot = Sum
        Q, fin = insert(Tot,Q,cur,fin,Mx-1)
        P1 = 0.
        for j in range(1,k+1):
            P0 = P1
            P1, cur = remove(Q,cur,fin,Mx-1)
            if(abs(Tot-P1) > 0.):  Tot = P0 + 1./(Tot-P1)
            Q, fin = insert(Tot,Q,cur,fin,Mx-1)
    # 'Tot' is now EPS sub M sup (1), the required series approximation.
    #
    # Calculate the final function value.
    FT = math.exp(Aterm*T) / Tfact * (0.5*F0 + Tot)
    return FT

def initC(Tmax: float, relErr: float, Alpha: float):
    """
 There are two parameters needed for the CRUMP
 routine.  They are related to the maximum time
 used in the inversion, the relative error required,
 and the largest real value of the poles of the 
 transformed function.  'INITC' claculates these
 values given the above information.
    """
    #   TMAX = The maximum time needed in the inversion
    #          (ie. use your largest time value)
    #   RELERR = The maximum relative error required
    #           (1E-6 is a reasonable value for most applications)
    #   ALPHA = The approximate exponential order of the
    #           untransformed function.  This is approx.
    #           the largest real value of the poles of the
    #           transformed function.
    #           (ALPHA=0.0 is a reasonable value for most applications)
    #   A and TFACT are calculated here as suggested
    #   by the Crump (1976) paper.
    Tfact = 0.7931 * Tmax
    Aterm = Alpha - math.log(relErr)/(2.*Tfact)
    return Tfact, Aterm   # return a tuple

def FS(P: complex, Znu: complex, Zkappa: complex, ZA: complex,
       Zsig: complex, Zdecay: complex, RPDPZ: complex,
       C0: float, Cim: float, 
       CSS: List[float], Ts: List[float], Ns: int) -> complex:
    # FS = L[C/C0]
    con1 = complex(1.,0.)
    PL = P + Zdecay
    PLS = cmath.sqrt(PL)
    SP = Zsig * PLS
    PRPDPZ = RPDPZ * PLS
    PTH = PLS * ctanh(SP)
    ATC = ZA*PTH + PL
    argE = Znu * (con1 - cmath.sqrt(con1 + Zkappa * ATC))
    argEr = argE.real
    argEi = argE.imag
    Zexp = complex(0.,0.)
    if(argEr > -75. and argEi > -75.): Zexp = cmath.exp(argE)
    if(argEr > -75. and argEi <= -75.):
        argE = complex(argEr,-75.)
        Zexp = cmath.exp(argE)
    if(argEr <= -75. and argEi > -75.):
        argE = complex(-75.,argEi)
        Zexp = cmath.exp(argE)
    if(abs(P) > 1.e-20):
        FS = Cim/PL + (C0bc(P,CSS,Ts,Ns) - Cim/PL) * Zexp
    else:
        FS = C0 * Zexp
    # Check for very small FS; set to zero if real and imaginary parts less then 1.E-20
    checkR = FS.real
    checkI = FS.imag
    AchR = abs(checkR)
    AchI = abs(checkI)
    if(AchR < 1.e-20  and AchI < 1.e-20):  FS = complex(0.,0.)
    if(AchR >= 1.e-20 and AchI < 1.e-20):  FS = complex(checkR,0.)
    if(AchR < 1.e-20  and AchI >= 1.e-20): FS = complex(0.,checkI)
    if(abs(P) > 1.e-10):
        FS = (FS - Cim/PL) * coshXY(PRPDPZ, SP) + Cim/PL
    else:
        FS = FS * coshXY(PRPDPZ, SP)
    return FS

def SER(Fcalc: complex, fact: float) -> float:
    # Function to calculate the Kth term in the series.
    # 'Fcalc' is the FS(A+K*PI*i/TFACT)
    # 'fact' is K*PI*T/TFACT
    SER = Fcalc.real * math.cos(fact) - Fcalc.imag * math.sin(fact)
    return SER

def C0bc(P: complex, CSS: List[float], TS: List[float], Ns: int) -> complex:
    # Function to evaluate time variable boundary condition for parent
    Sum = complex(0.,0.)
    for i in range(1,Ns):
        Sum = Sum + complex((CSS[i]-CSS[i-1]),0.) * cmath.exp(-P*complex(TS[i],0.))
    C0bc = (complex(CSS[0],0.)+Sum)/P
    return C0bc

def coshXY(X: complex, Y: complex) -> complex:
    # Function to evaluate COSH(X)/COSH(Y) with X, Y complex
    A = X.real
    B = X.imag
    C = Y.real
    D = Y.imag
    csB = math.cos(B)
    ssB = math.sin(B)
    csD = math.cos(D)
    ssD = math.sin(D)
    part1 = complex(csB,ssB)
    part2 = complex(csB,-ssB)
    part3 = complex(csD,ssD)
    part4 = complex(csD,-ssD)
    E1 = math.exp(A-C)
    AC = A + C
    if(AC < 30.):
        E2 = math.exp(-A-C)
    else:
        E2 = 0.
    E3 = 1.
    C2 = 2. * C
    if(C2 < 30.):
        E4 = math.exp(-C2)
    else:
        E4 = 0.
    coshXY = (E1*part1 + E2*part2) / (E3*part3 + E4*part4)
    return coshXY

def ctanh(X: complex) -> complex:
    # Function to evaluate hyperbolic tangent
    U = abs(X)
    XR = X.real
    XI = X.imag
    if((U > 1.e-6) and (U < 15.)):
        A1 = cmath.exp(X)
        A2 = cmath.exp(-X)
        ctanh = (A1-A2) / (A1+A2)
    # For large X, CTANH(X) -> 1.
    if(U >= 15.):
        ctanh = complex(1.,0.)
    # For small X, CTANH(X) -> X-X*X*X/3.
    if(U <= 1.e-6):
        ctanh = complex(0.,0.)
        AXR = abs(XR)
        AXI = abs(XI)
        if(AXR > 1.e-10  and AXI > 1e-10): ctanh = X - X*X*X / 3.
        if(AXR <= 1.e-10 and AXI > 1e-10):
            X_ = complex(0.,XI)
            ctanh = X_ - X_*X_*X_ / 3.
        if(AXR > 1.e-10  and AXI <= 1e-10):
            X_ = complex(XR,0.)
            ctanh = X_ - X_*X_*X_ / 3.
    return ctanh

# -------------------------------------------------
# ---- PROCEDURES TO HANDLE THE QUEUE PROCESSING
#
def initQ():
    # Initialise to an empty queue
    cur = 0
    fin = 0
    return cur,fin # return a tuple

def insert(val: float, Q: List[float], cur: int, fin: int, Mx:int):
    # Insert a value into the queue
    # Values are inserted at the end of the list
    nxt = fin+1
    if(nxt > Mx):
        nxt = 0
    #print("insert initially fin = "+str(fin)+"  Mx="+str(Mx)+"  cur="+str(cur)+" nxt="+str(nxt))
    if (nxt != cur):
        Q[fin] = val
        fin = nxt
    #print("       finally fin = "+str(fin))
    return Q,fin # return a tuple

def remove(Q: List[float], cur: int, fin: int, Mx:int):
    # Subroutine to get a value from the queue
    # Values are removed from the front of the list
    #print("remove initially fin = "+str(fin)+"  Mx="+str(Mx)+"  cur="+str(cur))
    if (cur != fin):
        val = Q[cur]
        cur = cur+1
        if(cur > Mx):
            cur = 0
    #print("       finally cur = "+str(cur))
    return val,cur # return a tuple

# C0 = 1.     # Source concentration at infinite time
# Cim = 0.    # Initial concentration in matrix and fractures
# v = 8.68056e-6  # Velocity in fracture (m/s)
# alpha = 0.  # Fracture dispersivity (m)
# b = 0.00012 # Fracture aperture (m)
# sep = 0.1   # Fracture spacing (m)
# theta = 0.35  # Matrix porosity
# tau = 1.      # Matrix tortuosity
# D0 = 1.e-11   # Diffusion coefficient in water (m2/s)
# R = 1.      # Fracture retardation factor
# Rp = 1.     # Matrix retardation factor
# THalf = 0.  # non-dcaying solute
# # Numerical inversion parameters used in CRUMP:
# AlphR = 0.   # alpha
# relErr = 1.e-6  # relatve error
# N = 15       # decides the number of terms in the EPAL summation
# #Variable source concentration
# Ns = 1
# Ts = [0.]
# Css = [1.]
# X = 0.6 # distance along fracture (m)
# Z = b/2.  # diatance from the centre of the fracture into matrix (m)
# T = float(4*24*60*60) # time = 4 days
# C, Cs = CraFlush(C0,Cim,v,alpha,b,sep,theta,tau,D0,R,Rp,THalf,X,Z,T,
#                  Ns,Ts,Css,AlphR,relErr,N)

# print("C = "+str(C)+"  Cs="+str(Cs)+"\nshould be C = +0.08641")