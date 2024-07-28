PROGRAM CRAFLUSH
Implicit NONE
Integer, PARAMETER :: dp = kind(0.d0) ! double precision
Integer, PARAMETER :: cdp = kind(DCMPLX(0.d0,0.d0)) ! complex double precision

!  This code was supplied in Fortran 77 format by:
!    Prof. Laura Toran
!    Dept. Earth and Environmental Science
!    Temple University, College of Science and Technology
!    Buery Hall
!    1901 North 13th Street
!    Philadelphia, PA 19122 USA
!
!  The original code was contained within the program
!    CRAFIT published in:
!    Toran, L. (2000), 'CRAFIT: A computer program for calibrating
!    breakthrough curves of CRAFLUSH, a one-dimensional fracture
!    flow and transport model'
!    Groundwater, 38: 430-434.
!    https://doi.org/10.1111/j.1745-6584.2000.tb00229.x
!
!  CRAFLUSH converted to Fortran 90 by Ignasi Puigdomenech
!  in April 2024
!
!  From the comments written by Laura Toran:
!  Hard-wired file names:
!  craf.in = CRAFLUSH input file
!  conc.out = output file of selected concentrations created by CRAFLUSH.
!  The CRAFLUSH version is CRAVS3 dated 4/22/92. CRAFLUSH was
!  provided by Ed Sudicky of the Univ of Waterloo.  Updates
!  be obtained from the University of Waterloo.
!  ------------------------------------------------------------
!
!         **************************************
!         *                                    *
!         *              CRAFLUSH              *
!         *                                    *
!         **************************************
!
!               Solved and programmed by:
!                     E.A. Sudicky
!        Waterloo Centre for Goundwater Research
!                University of Waterloo
!              Waterloo, Ontario, N2L 3G1
!
!              Copyright 1988 E.A. Sudicky
!
!                        WARNING
!
!         Duplication of this program or any part
!         thereof without the express permission
!         of E.A. Sudicky of the institute for
!         Groundwater Research is strictly forbidden
!
!         DISCLAIMER: While care has been exercised to ensure
!                     the accuracy of the coding, the author is
!                     not responsible for inadvertent errors.
!
!
!         Transport in a system of parallel fractures
!         with matrix diffusion.
!         This version allows the simulation of the flushing
!         of the initially contaminated matrix (and fractures)
!         with a solution of lower (or higher) concentration
!         entering the fractures.
!
!         This version uses the CRUMP algorithm to numerically
!         invert the Laplace transformed solution.
!
!     Definition of variables to enter in datafile.
!
!     NOTE: Use consistent units for input data
!
!         C0= Source concentration at fracture origin
!             NOTE:  When using variable source let C0=0
!         CIM=Initial concentration in matrix (and in fractures)
!         V= Velocity in fracture
!         ALPHA=Longitudinal dispersivity in fracture
!         B= Fracture aperture
!         SEP= Fracture spacing
!         THETA= Matrix porosity
!             NOTE: THETA=0.0 is alllowed, in which case a solution
!             with no matrix diffusion will be obtained
!         TAU= Matrix tortuosity
!         D0= Diffusion coeff in water
!         R= Fracture retardation factor
!         RP= Matrix retardation factor
!         THALF= Half-life
!
!         NX= Number of distance points along fracture
!         X= Distance along fracture
!         NZ= Number of distance points in matrix
!         ZM= Distance in matrix (B < ZM <= SEP)
!         NT= Number of time points
!         T= Time
!         NS= Number of points used to define variable source
!         TS= Times used to define variable source
!         CSS= Concentrations used in variable source at TS
!         2N+1=Number of terms used in CRUMP inversion
!             NOTE: The code reads "N"
!
COMPLEX(cdp) :: ZNU,ZKAPPA,ZA,ZSIG,ZDECAY,RPDPZ,F(103),P
REAL(dp) :: Q(208),FT,C,CS, DLAMDA
REAL(dp) :: ALPHA,ALPHR,APER,ATERM,B,C0,CIM,D,D0,R,RELERR,RP,RPDP
REAL(dp) :: SEP,SPCE,TAU,TFACT,THALF,THETA,TT,V,VV,VVX
REAL(dp) :: D_p  ! Note: "dp" is used as double precision KIND parameter
REAL(dp) X(51),ZM(32),T(150),CSS(51),TS(51)  ! dimensions changed by Prof. Toran
REAL(dp), PARAMETER :: Z = 0.d0
Integer :: MAXTRM, MAXQ, I,II,J,K,N,NS,NT,NX,NZ, IERR,IOUT,INP
! Dimensioning hints for increasing array sizes for larger problems
!       X(NX),Z(NZ+1),T(NT)
!       F(MAXTRM),Q(MAXQ)
! where: MAXTRM = 2N+1 = Number of terms in Fourier series
!                        for numerical inversion of the
!                        Laplace-transformed solution
!        MAXQ   = 4N+4
!
! For passing of arrays F and Q...
! the values below permit a maximum CRUMP "N" equal to 51
MAXTRM = 103
MAXQ = 208
!
WRITE(*,*) ' '
WRITE(*,*) ' ************************************************'
WRITE(*,*) ' *                Model: CRAFLUSH               *'
WRITE(*,*) ' *           Parallel crack model D > 0         *'
WRITE(*,*) ' *            (C) E.A. Sudicky, 1988            *'
WRITE(*,*) ' *   Waterloo Centre for Groundwater Research   *'
WRITE(*,*) ' *            University of Waterloo            *'
WRITE(*,*) ' *       Waterloo, Ontario, CANADA N2L 3G1      *'
WRITE(*,*) ' ************************************************'
WRITE(*,*) ' '
! WRITE(*,*) ' Enter input datafile name...'
! READ(*,'(A32)') FDATA
! WRITE(*,*) ' Give filename for listing...'
! READ(*,'(A32)') FLIST
!
! Read data from file 'craf.in'
INP = 17
OPEN(INP,FILE='craf.in',STATUS='OLD',IOSTAT=IERR)
IF(IERR /= 0) THEN
    WRITE(*,*) 'Error opening file "craf.in" for input.'
    IF (IERR == 2) WRITE(*,*) 'File not found'
    STOP
END IF
901 FORMAT(I5)
902 FORMAT(5F10.3)
903 FORMAT(F10.4)
904 FORMAT(E10.4)
READ(INP,903) C0     ! Source concentration at infinite time
READ(INP,903) CIM    ! Initial concentration in matrix and fractures
READ(INP,903) V      ! Velocity in fracture
READ(INP,903) ALPHA  ! Fracture dispersivity
READ(INP,904) B      ! Fracture aperture
READ(INP,903) SEP    ! Fracture spacing
READ(INP,903) THETA  ! Matrix porosity
READ(INP,903) TAU    ! Matrix tortuosity
READ(INP,904) D0     ! Diff. coeff. in water
READ(INP,903) R      ! Fracture retardation factor
READ(INP,903) RP     ! Matrix retardation factor
!...If you don't want decay, then read 0.0 for THALF
READ(INP,903) THALF  ! Half-life
! To help finding input errors, the input data is written to the screen
! as soon as it is read...
!...Write values of variables to screen.
!   NOTE: This format does not work in unix for screen write 
WRITE(*,8001) C0
WRITE(*,8009) CIM
WRITE(*,8000) V
WRITE(*,8002) ALPHA
WRITE(*,8005) B
WRITE(*,8010) SEP
WRITE(*,8015) THETA
WRITE(*,8025) TAU
WRITE(*,8030) D0
WRITE(*,8035) R
WRITE(*,8040) RP 
IF(THALF < 1.D-10) WRITE(*,8044)
IF(THALF >=1.D-10) WRITE(*,8045) THALF

!...Read the number of distance points along fracture.
READ(INP,901) NX     ! Number of distance points along fracture
WRITE(*,'("  NX=",I0,", fracture points:")') NX
IF(NX<=0 .OR. NX>SIZE(X)) THEN ! check for input errors
    WRITE(*,'(" Error: NX must be >0 and <",I0)') SIZE(X)
    STOP
ENDIF
!...Read the NX values of distance points along fracture.
READ(INP,902) (X(I),I=1,NX)       ! distance points along fracture
DO II=1,NX
    WRITE(*,'(1X,I5,1PG11.2)') II,X(II)
END DO

!...Read the number of distance points into matrix
!...Note: If you only want to evaluate the fracture solution, then
!         read NZ=0 and then don't read ZM(I)
READ(INP,901) NZ     ! number of distance points into matrix
WRITE(*,'("  NZ=",I0)') NZ
IF(NZ<0 .OR. NX>SIZE(ZM)) THEN ! check for input errors
    WRITE(*,'(" Error: NZ must be >=0 and <",I0)') SIZE(ZM)
    STOP
ENDIF
! Note: ZM(1) corresponds to the interface between
!       the fracture surface and the rock matrix.
NZ = NZ+1
!...Read the NZ values of distance points in the matrix.
!...Note: The values should be in the range B < ZM(I) <= SEP
IF(NZ > 1) THEN
    READ(INP,902) (ZM(I),I=2,NZ)  ! distance points into matrix
    WRITE(*,'("  Matrix points:")')
    DO II=2,NZ
        WRITE(*,'(1X,I5,1PG11.2)') (II-1),ZM(II)
        IF(ZM(II) < (B/2.d0) .OR. ZM(II) > (SEP/2.d0)) THEN
            WRITE(*,'(" Error: ZM(",I0,") must be >",1PG10.4," (1/2 fracture aperture) and <", &
                      1PG10.4," (1/2 fracture spacing)")') (II-1),(B/2.d0),(SEP/2.d0)
            STOP
        ENDIF
    END DO
ENDIF

!...Read the NT values of time points.
READ(INP,901) NT     ! number of time points
WRITE(*,'("  NT=",I0,", time points:")') NT
IF(NT<=0 .OR. NT>SIZE(T)) THEN  ! check for input errors
    WRITE(*,'(" Error: NT must be >0 and <",I0)') SIZE(T)
    STOP
ENDIF
READ(INP,902) (T(I),I=1,NT)       ! time points
DO II=1,NT
    WRITE(*,'(1X,I5,1PG11.2)') II,T(II)
    IF(T(II) <= 0.d0) THEN
        WRITE(*,'(" Error: time must be >0")')
        STOP
    ENDIF
END DO

!...Read parameter values for CRUMP inversion
!...Recommend: ALPHR=0.0, RELERR=1.E-6, N >= 11
!...where 2N+1 is the number of terms in
!...the Fourier series for numerical inversion of the Laplace transform
!...NOTE: If the solution oscillates about the steady state
!         concentration, try RELERR=1.E-5 or 1.E-4 etc., or increase N.
READ(INP,903) ALPHR  ! alpha
READ(INP,904) RELERR ! relative error
READ(INP,901) N      ! decides the number of terms in the EPAL summation.
WRITE(*,8046) ALPHR,RELERR,N
IF(ALPHR<0.d0 .OR. RELERR<0.d0) THEN  ! check for input errors
    WRITE(*,'(" Error: both  ALPHR and RELERR must be >=0")')
    STOP
ENDIF
IF(N<=0) THEN  ! check for input errors
    WRITE(*,'(" Error: N must be >0")')
    STOP
ENDIF

!...Variable source: Read number of points, time and concentration
READ(INP,901) NS     ! number of points in the variable source concentration
IF(NS<=0 .OR. NS>SIZE(TS)) THEN  ! check for input errors
    WRITE(*,'("  NS=",I0)') NS
    WRITE(*,'(" Error: NS must be >0 and <",I0)') SIZE(TS)
    STOP
ENDIF
READ(INP,902) (TS(I),I=1,NS)      ! the time for source variation
READ(INP,902) (CSS(I),I=1,NS)     ! the new source concentration
WRITE(*,8950)
DO II=1,NS
    WRITE(*,8951) TS(II),CSS(II)
8951 FORMAT(1X,1PG10.4,5X,0PF10.5)
    IF(TS(II) < 0.d0) THEN
        WRITE(*,'(" Error: time must be >=0")')
        STOP
    ENDIF
END DO
CLOSE(INP)
!
!...Decay constant
DLAMDA = 0.d0
IF(THALF >=1.D-10) DLAMDA = DLOG(2.d0)/THALF
APER = B
SPCE = SEP
8001 FORMAT('  Source concentration at infinite time= ',F10.4)
8009 FORMAT('  Initial concentration in matrix and fractures= ',F10.4)
8000 FORMAT('  Velocity in fracture= ',1PG12.5)
8002 FORMAT('  Fracture dispersivity= ',F10.3)
8005 FORMAT('  Fracture aperture= ',1PG11.4)
8010 FORMAT('  Fracture spacing= ',F10.3)
8015 FORMAT('  Matrix porosity= ',F10.3)
8025 FORMAT('  Matrix tortuosity= ',F10.3)
8030 FORMAT('  Diff. coeff. in water= ',1PG11.4)
8035 FORMAT('  Fracture retardation factor= ',F10.3)
8040 FORMAT('  Matrix retardation factor= ',F10.3)
8044 FORMAT('  Non-decaying solute')
8045 FORMAT('  Half-life= ',1PG11.4)
8046 FORMAT(//,'  Numerical inversion parameters used in CRUMP',//, &
     '  ALPHR= ',F6.3,/,'  RELERR= ',1PG11.4,/,'  N= ',I3)
8950 FORMAT(//,'  Variable source concentration',//,'     Time  ',5X,'Concentration')

WRITE(*,8200)
WRITE(*,8210)
8200 FORMAT(//'    Time      Distance: along Fracture,  in Matrix     C(t)     Steady C')
8210 FORMAT('    ----      --------- --------------   ---------     ----     --------')

IOUT=16
OPEN(UNIT=IOUT,FILE='conc.out',STATUS='unknown')
WRITE(IOUT,8041)
8041 FORMAT(/20X,'************************************************', &
            /20X,'*                Model: CRAFLUSH               *', &
            /20X,'*           Parallel Crack Model D > 0         *', &
            /20X,'*            (C) E.A. Sudicky, 1988            *', &
            /20X,'*   Waterloo Centre for Groundwater Research   *', &
            /20X,'*            University of Waterloo            *', &
            /20X,'*       Waterloo, Ontario, CANADA N2L 3G1      *', &
            /20X,'************************************************',//)
WRITE(IOUT,8001) C0
WRITE(IOUT,8009) CIM
WRITE(IOUT,8000) V
WRITE(IOUT,8002) ALPHA
WRITE(IOUT,8005) APER
WRITE(IOUT,8010) SPCE
WRITE(IOUT,8015) THETA
WRITE(IOUT,8025) TAU
WRITE(IOUT,8030) D0
WRITE(IOUT,8035) R
WRITE(IOUT,8040) RP 
IF(THALF < 1.D-10) WRITE(IOUT,8044)
IF(THALF >=1.D-10) WRITE(IOUT,8045) THALF
WRITE(IOUT,8046) ALPHR,RELERR,N
WRITE(IOUT,8950)
DO II=1,NS
    WRITE(IOUT,8951) TS(II),CSS(II)
END DO
WRITE(IOUT,8200)
WRITE(IOUT,8210)
!
B = B/2.d0
! Note that even if NZ=0, at least one ZM is used in the calculations
ZM(1) = B
SEP = SEP/2.d0
D_p = D0*TAU
D = ALPHA*V + D0
!...Calculate parameters used in the Laplace transformed solution L[C]
!...They should be stored as complex values
ZA = DCMPLX(THETA*DSQRT(RP*D_p)/(B*R), Z)
ZKAPPA = DCMPLX(4.d0*R*D/V/V, Z)
RPDP = DSQRT(RP/D_p)
ZSIG = DCMPLX(RPDP*(SEP-B), Z)
ZDECAY = DCMPLX(DLAMDA, Z)
VV = V/2.d0/D

DO I=1,NT
!  ...Initialization for CRUMP
    TT=T(I)
    CALL INITC(TT,RELERR,ALPHR,ATERM,TFACT)
    DO J=1,NX
        VVX = VV*X(J)
        ZNU = DCMPLX(VVX, Z)
        DO K=1,NZ
            IF(X(J) < 1.D-10) THEN
                IF(ZM(K) <=B) THEN ! goto 8043
                    C = C0
                    CS = C0
                ENDIF
            ELSE  ! not 8043
                RPDPZ = DCMPLX(RPDP*(SEP-ZM(K)), Z)
                ! Steady-state solution given by:
                ! Limit as P --> 0 of P*L[C]
                !   or C0 if no decay
                IF(THALF >=1.D-10) THEN
                    P = DCMPLX(Z, Z)
                    CS = DBLE(FS(P,ZNU,ZKAPPA,ZA,ZSIG,ZDECAY,RPDPZ,C0,CIM,CSS,TS,NS))
                ELSE
                    CS = C0
                END IF
                ! Transient solution
                CALL CRUMP(FT,TT,ATERM,TFACT,ZNU,ZKAPPA,ZA,ZSIG,ZDECAY,RPDPZ,F,Q,N,MAXTRM,MAXQ,C0,CIM,CSS,TS,NS)
                C = FT
            END IF ! 8043 / 8042

        ! write to file 'conc.out'
        !8042 WRITE(8,*) T(I),C
        WRITE(*,8060) T(I),X(J),ZM(K),C,CS
        WRITE(IOUT,8060) T(I),X(J),ZM(K),C,CS  ! 8042
        8060 FORMAT(1X,1PG10.4,0PF27.3,F12.3,F9.5,F13.5)
        END DO
    END DO
END DO

CLOSE(IOUT)

CONTAINS

!---------------------------------------------------------------------
! Routine for numerical inversion of Laplace
! transformations using the CRUMP algorithm.
! This uses the 'EPAL' (EPsilon ALgorithm) method
! of series summation as suggested in the CRUMP paper.
!
! Reference:  Kenny S. Crump, 'Numerical
!   inversion of Laplace transforms using a
!   Fourier series approximation', Journal of
!   the Association for Computing Machinery
!   Vol 23, No. 1, Jan 1976, pp. 89-96.
!   https://doi.org/10.1145/321921.321931
!
! Programmed by:  Frank Letniowski
!   Institute for Groundwater Research
!   University of Waterloo
!   May, 1987
!
!----------------------------------------------------------------
! The CRUMP subroutine numerically inverts a Laplace
! transformed function FS at a particular time T
! returning the value in FT.  FS must be a complex 
! function of a complex argument.  Here, X is the Laplace
! variable 'P' which is complex.  For maximum efficiency,
! you may want to store the FS values for each X in an
! array (i.e. precompute them, along with the needed
! X-values) outside your time loop and then substitute
! the FS array into the CRUMP Fourier Series that involves
! time.  That is, it is not necessary to compute FS(X)
! for each value of the time 'T' since X depends only
! on TFACT through TMAX.
!
! The parameters A and TFACT should be as suggested
! above or in the paper.  N is related to the number
! of terms used in the summation i.e. 2N+1 terms are
! calculated in the 'EPAL' method of summation.
!
! Dimensioning hints:
!                    F(MAXTRM),Q(MAXQ)
! where: MAXTRM = 2N+1
!        MAXQ   = 4N+4
SUBROUTINE CRUMP (FT,T,ATERM,TFACT,ZNU,ZKAPPA,ZA,ZSIG,ZDECAY,RPDPZ,F,Q,N,MAXTRM,MAXQ,C0,CIM,CSS,TS,NS)
Implicit NONE
REAL(dp), INTENT(OUT) :: FT
COMPLEX(cdp), INTENT(OUT) :: F(MAXTRM) ! a working array
REAL(dp), INTENT(IN) :: T, ATERM, TFACT, C0, CIM
COMPLEX(cdp), INTENT (IN) :: ZNU,ZKAPPA,ZA,ZSIG,ZDECAY,RPDPZ
Integer, INTENT (IN) :: N, MAXTRM, MAXQ, NS
REAL(dp), INTENT (OUT) :: Q(MAXQ) ! a working array
REAL(dp) :: SUM,P0,P1,TEMP,TOT,F0
Integer :: M,J,K,Mx,CUR,FIN
COMPLEX(cdp) :: X
REAL(dp), INTENT (IN) :: CSS(51),TS(51)

REAL(dp), PARAMETER :: PI = 3.1415926535897932385d0

M = 2*N+1
X = DCMPLX(ATERM,0.d0)

F0 = DBLE(FS(X,ZNU,ZKAPPA,ZA,ZSIG,ZDECAY,RPDPZ,C0,CIM,CSS,TS,NS))
! Calculate the function values needed
DO K=1,M
    X = DCMPLX(ATERM, DBLE(K)*PI/TFACT)
    F(K) = FS(X,ZNU,ZKAPPA,ZA,ZSIG,ZDECAY,RPDPZ,C0,CIM,CSS,TS,NS)
END DO

! Initialise the queue needed for calculation of the SUM.
! Mx is the size of the queue needed.
! The queue stores EPSILON values needed in the
! EPAL method to calculate the series.
Mx = 2*M + 2
CALL INITQ(CUR,FIN)  ! set CUR=1 and FIN=1

! Calculate the series by use of the EPAL method.
! SER(I,F(I),TEMP) calculates the Ith term in the series.
TEMP = DBLE(PI*T/TFACT)
SUM = SER(F(1),TEMP)

! Insert into the queue
CALL INSERT(SUM,Q,CUR,FIN,Mx,MAXQ)

DO K=2,M
    SUM = SUM + SER(F(K),DBLE(K)*TEMP)
!   SUM is the Kth partial sum of the series
!   Calculate EPS sub K sup (1) which is TOT
!   upon leaving the 'J' loop
    TOT = SUM
    CALL INSERT(TOT,Q,CUR,FIN,Mx,MAXQ)
    P1 = DBLE(0.d0)
    DO J=2,K
        P0=P1
        CALL REMOVE(P1,Q,CUR,FIN,Mx,MAXQ)
        IF(DABS(TOT-P1) > 0.D0) THEN
            TOT = P0 + 1.d0/(TOT-P1)
        END IF
        CALL INSERT(TOT,Q,CUR,FIN,Mx,MAXQ)
    END DO
END DO
!
! TOT is now EPS sub M sup (1), the required series approximation.
!
! Calculate the final function value.
FT = DEXP(DBLE(ATERM*T)) / TFACT * (0.5d0*F0 + TOT)
RETURN
END SUBROUTINE CRUMP

! There are two parameters needed for the CRUMP
! routine.  They are related to the maximum time
! used in the inversion, the relative error required,
! and the largest real value of the poles of the 
! transformed function.  'INITC' claculates these
! values given the above information.
!
SUBROUTINE INITC(TMAX,RELERR,ALPHA,ATERM,TFACT)
!          -----
!     TMAX = The maximum time needed in the inversion
!            (ie. use your largest time value)
!     RELERR = The maximum relative error required
!            (1E-6 is a reasonable value for most applications)
!     ALPHA = The approximate exponential order of the
!             untransformed function.  This is approx.
!             the largest real value of the poles of the
!             transformed function.
!             (ALPHA=0.0 is a reasonable value for most applications)
!     A and TFACT are calculated here as suggested
!     by the Crump (1976) paper.
!
Implicit NONE
REAL(dp), INTENT(IN) :: TMAX,RELERR,ALPHA
REAL(dp), INTENT(OUT) :: ATERM,TFACT
         TFACT = 0.7931d0*TMAX
         ATERM = ALPHA - DLOG(RELERR)/(2.d0*TFACT)
RETURN
END SUBROUTINE INITC

!---------------------------------------------------------------------
COMPLEX(cdp) FUNCTION FS (P,ZNU,ZKAPPA,ZA,ZSIG,ZDECAY,RPDPZ,C0,CIM,CSS,TS,NS)
Implicit NONE
COMPLEX(cdp), INTENT (IN) :: P, ZNU, ZKAPPA, ZA, ZSIG, ZDECAY, RPDPZ
REAL(dp), INTENT (IN) :: CSS(51),TS(51)
REAL(dp), INTENT (IN) :: C0, CIM
Integer, INTENT (IN) :: NS
REAL(dp) :: ARGER,ARGEI, CHECKR,CHECKI, ACHR,ACHI
REAL(dp), PARAMETER :: Z = 0.d0
COMPLEX(cdp) :: ARGE,SP,PL,CON1,PLS,PRPDPZ,PTH,ATC,ZEXP
! COMPLEX(cdp) :: BTC,ABTC
! FS = L[C/C0]
CON1 = DCMPLX(1.d0,0.d0)
PL = P + ZDECAY
PLS = CDSQRT(PL)
SP = ZSIG*PLS
PRPDPZ = RPDPZ*PLS
PTH = PLS*CTANH(SP)
ATC = ZA*PTH + PL
!  BTC=(ZA+PTH)/PTH
!  ABTC=BTC/ATC
ARGE = ZNU * (CON1-CDSQRT(CON1 + ZKAPPA*ATC))
ARGER = DBLE(ARGE)
ARGEI = DIMAG(ARGE)
ZEXP = DCMPLX(Z,Z)
IF (ARGER > -75.d0 .AND. ARGEI > -75.d0) THEN
    ZEXP = CDEXP(ARGE)
END IF
IF(ARGER > -75.d0 .AND. ARGEI <=-75.d0) THEN
    ARGE = DCMPLX(ARGER,-75.d0)
    ZEXP = CDEXP(ARGE)
END IF
IF(ARGER <=-75.d0 .AND. ARGEI > -75.d0) THEN
    ARGE=DCMPLX(-75.d0,ARGEI)
    ZEXP=CDEXP(ARGE)
END IF
IF(CDABS(P) > 1.D-20) THEN
    FS = CIM/PL + (C0BC(P,CSS,TS,NS) - CIM/PL) * ZEXP
ELSE
    FS = C0*ZEXP
END IF
! Check for very small FS; set to zero if real and imaginary parts less then 1.E-20
CHECKR = DBLE(FS)
CHECKI = DIMAG(FS)
ACHR = DABS(CHECKR)
ACHI = DABS(CHECKI)
IF(ACHR < 1.D-20 .AND. ACHI < 1.D-20) FS = DCMPLX(Z,Z)
IF(ACHR >=1.D-20 .AND. ACHI < 1.D-20) FS = DCMPLX(CHECKR,Z)
IF(ACHR < 1.D-20 .AND. ACHI >=1.D-20) FS = DCMPLX(Z,CHECKI)
IF(CDABS(P) > 1.D-20) THEN
    FS = (FS-CIM/PL)*COSHXY(PRPDPZ,SP) + CIM/PL
ELSE
    FS = FS * COSHXY(PRPDPZ,SP)
ENDIF
RETURN
END FUNCTION FS

!---------------------------------------------------------------------
REAL(dp) FUNCTION SER (FCALC,FACT)
!  Function to calculate the Kth term in the series.
!  FCALC is the FS(A+K*PI*i/TFACT)
!  FACT is K*PI*T/TFACT
Implicit NONE
REAL(dp), INTENT (IN) :: FACT
COMPLEX(cdp), INTENT(IN) :: FCALC
SER = DBLE(FCALC)*DCOS(FACT) - DBLE(DIMAG(FCALC))*DSIN(FACT)
RETURN
END FUNCTION SER

!---------------------------------------------------------------------
COMPLEX(cdp) FUNCTION C0BC (P,CSS,TS,NS)
! Function to evaluate time variable boundary condition for parent
Implicit NONE
COMPLEX(cdp), INTENT (IN) :: P
REAL(dp), INTENT (IN) :: CSS(51),TS(51)
Integer, INTENT (IN) :: NS
COMPLEX(cdp) :: SUM
Integer :: I
REAL(dp), PARAMETER :: Z = 0.d0
SUM = DCMPLX(0.d0,0.d0)
DO I=2,NS
    SUM = SUM + DCMPLX((CSS(I)-CSS(I-1)),Z)*CDEXP(-P*DCMPLX(TS(I),Z))
END DO
C0BC = (DCMPLX(CSS(1),Z)+SUM)/P
RETURN
END FUNCTION C0BC

!---------------------------------------------------------------------
! Function to evaluate COSH(X)/COSH(Y) with X, Y complex
COMPLEX(cdp) FUNCTION COSHXY(X,Y)
Implicit NONE
COMPLEX(cdp), INTENT(IN) :: X,Y
COMPLEX(cdp) :: PART1,PART2,PART3,PART4
REAL(dp) :: A,B,C,D,E1,E2,E3,E4,AC,C2,CSB,SSB,CSD,SSD
A = DBLE(X)
B = DBLE(DIMAG(X))
C = DBLE(Y)
D = DBLE(DIMAG(Y))
CSB = DCOS(B)
SSB = DSIN(B)
CSD = DCOS(D)
SSD = DSIN(D)
PART1 = DCMPLX(CSB,SSB)
PART2 = DCMPLX(CSB,-SSB)
PART3 = DCMPLX(CSD,SSD)
PART4 = DCMPLX(CSD,-SSD)
E1 = DEXP(A-C)
AC = A + C
IF(AC < 30.D0) THEN
    E2 = DEXP(-A-C)
ELSE
    E2 = 0.D0
ENDIF
E3 = 1.D0
C2 = 2.D0*C
IF(C2 < 30.D0) THEN
    E4 = DEXP(-C2)
ELSE
    E4=0.D0
ENDIF
!ES1 = SNGL(E1)
!ES2 = SNGL(E2)
!ES3 = SNGL(E3)
!ES4 = SNGL(E4)
COSHXY = (E1*PART1 + E2*PART2) / (E3*PART3 + E4*PART4)
RETURN
END FUNCTION COSHXY

!---------------------------------------------------------------------
! Function to evaluate hyperbolic tangent
COMPLEX(cdp) FUNCTION CTANH(X)
Implicit NONE
COMPLEX(cdp), INTENT(INOUT) :: X
COMPLEX(cdp) :: A1,A2
REAL(dp) :: U, XR,XI,AXR,AXI
U = CDABS(X)
XR = DBLE(X)
XI = DIMAG(X)
IF(U > 1.D-06 .AND. U < 15.D0) THEN
    A1 = CDEXP(X)
    A2 = CDEXP(-X)
    CTANH = (A1-A2)/(A1+A2)
END IF
! For large X, CTANH(X) -> 1.
IF(U >= 15.D0) THEN
    CTANH = DCMPLX(1.D0,0.D0)
END IF
! For small X, CTANH(X) -> X-X*X*X/3.
IF(U <=1.D-06) THEN
    CTANH = DCMPLX(0.D0,0.D0)
    AXR = DABS(XR)
    AXI = DABS(XI)
    IF(AXR > 1.D-10 .AND. AXI > 1.D-10) CTANH = X-X*X*X/3.
    IF(AXR <=1.D-10 .AND. AXI > 1.D-10) THEN
        X = DCMPLX(0.D0,XI)
        CTANH = X-X*X*X/3.
    END IF
    IF(AXR > 1.D-10 .AND. AXI <=1.D-10) THEN
        X = DCMPLX(XR,0.D0)
        CTANH = X-X*X*X/3.
    END IF
END IF
RETURN
END FUNCTION CTANH

!---------------------------------------------------------------------
! ------ PROCEDURES TO HANDLE THE QUEUE PROCESSING
!
!----------------
SUBROUTINE INITQ (CUR,FIN)
!          -----     
! Initialise to an empty queue
Implicit NONE
Integer, INTENT (OUT) :: CUR,FIN
    CUR=1
    FIN=1
RETURN
END SUBROUTINE INITQ

!----------------
SUBROUTINE INSERT (VAL,Q,CUR,FIN,Mx,MAXQ)
!          ------
! Insert a value into the queue
! Values are inserted at the end of the list
Implicit NONE
Integer, INTENT (IN) :: MAXQ
REAL(dp), INTENT (INOUT) :: Q(MAXQ)
REAL(dp), INTENT (IN) :: VAL
Integer, INTENT (IN) :: CUR, Mx
Integer, INTENT (INOUT) :: FIN
Integer :: NEXT
NEXT=FIN+1
IF(NEXT > Mx)THEN
    NEXT=1
END IF
IF(NEXT /= CUR)THEN
    Q(FIN)=VAL
    FIN=NEXT
END IF
RETURN
END SUBROUTINE INSERT

!----------------
SUBROUTINE REMOVE (VAL,Q,CUR,FIN,Mx,MAXQ)
!          ------
! Subroutine to get a value from the queue
! Values are removed from the front of the list
Implicit NONE
Integer, INTENT (IN) :: MAXQ
REAL(dp), INTENT (IN) :: Q(MAXQ)
REAL(dp), INTENT (OUT) :: VAL
Integer, INTENT (INOUT) :: CUR
Integer, INTENT (IN) :: FIN,Mx
IF(CUR /= FIN)THEN
    VAL = Q(CUR)
    CUR = CUR+1
    IF(CUR > Mx)THEN
        CUR=1
    END IF
END IF
RETURN
END SUBROUTINE REMOVE

END PROGRAM CRAFLUSH