c
c  This program uses Latin hypercube sampling (LHS) to
c  generate fracture parameters, 
c  creates craflush files, runs craflush, then calculates
c  an error summary.  The Latin hypercube program was
c  obtained from Max Morris at Oak Ridge National Laboratory
c
c  Variables created by LHS are b,sep,theta and
c  they MUST be listed in that order in lhs.input
c
c   max runs=5000, max obs data=25 (changed max runs in LHS code
c   from 200 to 5000), max parameters in LHS is 50, in craflush
c   the array sizes are generally 51
c
c  Hard-wired file names:
c  lhs.seeds = seeds for lhs input (unit 1) USER INPUT 
c              this file is re-written each run
c  lhs.input = lhs input file (parameter ranges) (unit 1) USER INPUT
c              parameter order as given in CRAFLUSH
c              see "READ LHS.INPUT" for instructions
c  lhs.output = table created by latin (again unit 1)**
c  craf.in = craflush input file created by latin (overwrites 
c            for each run) unit 7
c  conc.out = output file of selected concentrations created
c             by craflush, used for debugging (will overwrite for each run)
c             unit 8 in this program
c  interp.data = data from observed bthru curves (nobs points) 
c             (unit 9) USER INPUT
c             see "READ INTERP.DATA"  for further input instructions
c  start.cra = initial data for craflush input, based on CRAFLUSH, 
               but modified for column description (unit 10)
c              see "READ START.CRA" for further input instructions
c  err.sum = summary of errors by run number (accumulate for
c            each run...unit 11)**
c
c  ** These two files are used to process the output.  Load into excel
c     or some other spreadsheet.  The err.sum entries will match up
c     with the lhs.output parameter values to show which values gave
c     which answers.
c
c
c     See READ LHS.INPUT and .
c
c-----------------SUBROUTINE DESCRIPTIONS/MODS-------------------
c  SUBROUTINE CRAFWRIT uses formatting directions from CRAFLUSH
c  to read default parameters and write new input files.
c  The CRAFLUSH version is cravs3 dated 4/22/92.
c  If you use another version of CRAFLUSH you should
c  substitute read and write statements for that version
c  number.  Look for the comment on changing the count for
c  NZ and make that change in any substitution.
c
c  SUBROUTINE CRAFLUSH runs the craflush program  CRAFLUSH was
c  provided by Ed Sudicky of the Univ of Waterloo.  Updates
c  be obtained from the University of Waterloo.
c
c  THE FOLLOWING MODIFICATIONS SHOULD BE MADE TO CRAFLUSH TO
c  INSERT INTO CRAFIT:
c
c  1. Increase the dimension for time so that more times can
c  be printed.  The new dimension statement should look like this:
c        DIMENSION X(51),ZM(32),T(150),CSS(51),TS(51)
c
c  2.  Comment out all of the existing write statements (files
c  created are too big; screen print out will slow the Crafit
c  program).  
c
c  3.  Substitute the following simple write statment after
c  the ENDIF statement that begins at 8043.  Unit 8 is
c  the file "conc.out"
c 
c
c      8042  WRITE(8,*) T(I),C
c
c  4.  Comment out the read and writes that have the CRAFLUSH 
c   program prompt the user for input and output file names.  All
c   file names are hard-wired in CRAFIT.
c
c    SUBROUTINE FINDMIN calculates the minimum separation
c    between the observed and modeled breakthrough. 
c    This distance is calculated by computing the hypotenuse 
c    of the x and y distance for each interpolation point and 
c    saving the minimum distance between the modeled and observed 
c    breakthrough curves.  The entire breakthrough curve is 
c    searched for a minimum distance.  A more efficient method 
c    would be to limit the number of points searched, but this 
c    would require additional user input that specifies the time 
c    range sampled.  To minimize the amount of user input, the 
c    program does the extra calculations all along the breakthrough 
c    curve.
c
c
c
c------------------Main Program----------------------------------
c
      implicit real*8(a-h,o-z)
      common/rand/ix,iy,iz
      dimension x(5000,50)
      dimension xlo(50),xhi(50),logflag(50)
c
c set number of interpolation points 
c right now the program prompts for this, but
c it could be put into lhs.input
c
c nobs will be passed to subroutine findmin
c
      WRITE(*,*) ' How many interpolation points (max 25)...'
      READ(*,4682) nobs 
4682  FORMAT(i4)
c Read seeds for random number generator, usually 5-digit integers
      open(unit=1,file='lhs.seeds')
      read(1,*) ix,iy,iz
      close(unit=1)
      idummy=12345
c
c   READ LHS.INPUT
c   Read the ranges for Latin hypercube sampling
c   n=number of runs, k=number of variables
c   then k more lines, with low value, high value and whether value should be
c         distributed uniformly on untransformed (0) 
c         or log-transformed (1) scale.
c   The order of the k variables in the same as the order in CRAFLUSH
c
      open(unit=1,file='lhs.input')
      read(1,*) n,k
      do 90 i=1,k
      read(1,*) xlo(i),xhi(i),logflag(i)
      if(logflag(i).eq.2) then
        xlo(i)=dlog10(xlo(i))
        xhi(i)=dlog10(xhi(i))
        endif
90    continue
c Generate integer (1 through n) LHS.
      do 10 ik=1,k
      do 12 in=1,n
      x(in,ik)= random(idummy)-3.d0
12    continue
      do 14 in=1,n
      test=-1.0d0
      jtest=0
      do 16 jn=1,n
      if(x(jn,ik).lt.test) then
        jtest=jn
        test=x(jn,ik)
        endif
16    continue
      x(jtest,ik)=in
14    continue
10    continue
c Convert to desired variables.
      do 30 ik=1,k
      do 30 in=1,n
      x(in,ik)=xlo(ik)+(xhi(ik)-xlo(ik))*(x(in,ik)-1.d0)/(n-1)
      if(logflag(ik).eq.2)
     $   x(in,ik)=10.d0**x(in,ik)
30    continue
c Report result.
      open(unit=1,file='lhs.output')
      write(1,*) ' ap         spac        porm        r           rp'
      do 987 ii=1,n  
987   write(1,876) (x(ii,jj),jj=1,k)
876   format(10e12.5)
      close(unit=1)
c
c ------------------------------------------------------------
c pass the randomly generated variable and write craflush files
c the order of the variables in LHS.INPUT must match the
c order of the variables in craf.in (and shown below)
c
      do 40 ii=1,n
      ap=x(ii,1)
      spac=x(ii,2)
      porm=x(ii,3)
      r=x(ii,4)
      rp=x(ii,5)
c
c  could add other variables here by adding more to the 
c  "do 40" loop and passing more variables through "call crafwrit"
c
      open(8,file='conc.out',status='unknown')
      call crafwrit(ap,spac,porm,r,rp)
c -------------------------------------------------------------
c   run craflush
c
      call craflush
c -------------------------------------------------------------
c   Using the conc.out file written in craflush,
c   here the results must be compared with observed 
c   (or interpolated) data.
      call findmin(nobs,ii)
      close(8)
40    continue
c Replace seeds for random number generator.
      open(unit=1,file='lhs.seeds')
      write(1,*) ix,iy,iz
      close(unit=1)
c
      stop
      end
c
c Latin function
c
      FUNCTION RANDOM(L)
C
C	   ALGORITHM AS 183 APPL. STATIST. (1982) VOL.31, NO.2
C
C	   RETURNS A PSEUDO-RANDOM NUMBER RECTANGULARLY DISTRIBUTED
C	   BETWEEN 0 AND 1.
C
C	   IX, IY AND IZ SHOULD BE SET TO INTEGER VALUES BETWEEN
C	   1 AND 30000 BEFORE FIRST ENTRY
C
C	   INTEGER ARITHMETIC UP TO 30323 IS REQUIRED
C
        implicit real*8(a-h,o-z)
	COMMON/RAND/IX, IY, IZ
	IX = 171 * MOD(IX, 177) - 2  * (IX/177)
	IY = 172 * MOD(IY, 176) - 35 * (IY/176)
	IZ = 170 * MOD(IZ, 178) - 63 * (IZ/178)
C
	IF (IX .LT. 0) IX = IX + 30269
	IF (IY .LT. 0) IY = IY + 30307
	IF (IZ .LT. 0) IZ = IZ + 30323
C
C	   IF INTEGER ARITHMETIC UP TO 5212632 IS AVAILABLE,
C	   THE PRECEDING 6 STATEMENTS MAY BE REPLACED BY
C
C	IX = MOD(171 * IX, 30269)
C	IY = MOD(172 * IY, 30307)
C	IZ = MOD(170 * IZ, 30323)
C
C	   ON SOME MACHINES, THIS MAY SLIGHTLY INCREASE
C	   THE SPEED. THE RESULTS WILL BE IDENTICAL.
C
	RANDOM = AMOD(FLOAT(IX)/30269.0 +
     *                FLOAT(IY)/30307.0 + FLOAT(IZ)/30323.0, 1.0)
	RETURN
	END
c
c--------------------------------------------------------------
c
      SUBROUTINE CRAFWRIT(ap,spac,porm,r,rp)
c
      implicit real*8(a-h,o-z)
      DIMENSION T(150),X(51),CSS(51),TS(51),ZM(32)
C
C
C  READ START.CRA:  Read a file with initial CRAFLUSH data and new header. 
c  This piece of code is modified from the Craflush
c  program and follows craflush formats...except
c  ADDED LINE AT TOP TO READ COLUMN DIMENSIONS & FLOW RATE
C    ri=rate of injection, cl=column length, cv=column volume
C
C
      open(10,file='start.cra',status='old')
      read(10,*) ri,cl,cv
      READ(10,903) C0
      READ(10,903) CIM
      READ(10,903) V
      READ(10,903) ALPHA
903   FORMAT(F10.4)
      READ(10,904) B
      READ(10,903) SEP
      READ(10,903) THETA
      READ(10,903) TAU
      READ(10,904) D0
904   FORMAT(E10.4)
c  ...rename the retardation variables so that the
c  ...latin hypercube variables are not overwritten
      READ(10,903) Rold
      READ(10,903) RPold
C
C  ...IF YOU DON'T WANT DECAY, THEN READ 0.0 FOR THALF
      READ(10,903) THALF
C
C  ...DECAY CONSTANT
      DLAMDA=0.0
      IF(THALF.GE.1.E-10) DLAMDA=ALOG(2.)/THALF
C  ...Read the number of distance points along fracture.
      READ(10,901) NX
901   FORMAT(I5)
C  ...Read the NX values of distance points along fracture.
        READ(10,902) (X(I),I=1,NX)
902     FORMAT(5F10.3)
C  ...Read the number of distance points into matrix
C  ...Note: If you only want to evaluate the fracture solution, then
C           read NZ=0 and then don't read ZM(I)
      READ(10,901) NZ
      NZ=NZ+1
C  ...Read the NZ values of distance points in the matrix.
C  ...Note: The values should be in the range B < ZM(I) <= SEP
      IF(NZ.GT.1) THEN
        READ(10,902) (ZM(I),I=2,NZ)
      ENDIF
c
c  *** I added a statement here modifying the Craflush
c  *** input instructions to return nz to zero.  This will cause
c  *** trouble if multiple nz's (distance along fracture) are read.
c  *** This statement had to be added to write the files
c  *** correctly; the original CRAFLUSH section of code was
c  *** ment for reading and the NZ increase was ok there.
c      
	nz=nz-1
      READ(10,901) NT
C  ...Read the NT values of time points.
        READ(10,902) (T(I),I=1,NT)
C
C  ...READ PARAMETER VALUES FOR CRUMP INVERSION
C  ...RECOMMEND: ALPHR=0.0, RELERR=1.E-6, N >= 11
C  ...WHERE 2N+1 IS THE NUMBER OF TERMS IN
C  ...THE FOURIER SERIES FOR NUMERICAL INVERSION OF THE
C  ...LAPLACE TRANSFORM
C  ...NOTE: IF THE SOLUTION OSCILLATES ABOUT THE STEADY STATE
C           CONCENTRATION, TRY RELERR=1.E-5 OR 1.E-4 ETC. OR
C           INCREASE N.
C
      READ(10,903) ALPHR
      READ(10,904) RELERR
      READ(10,901) N
C
C  ...READ NUMBER OF POINTS, TIME AND CONCENTRATION OF 
C  ...VARIABLE SOURCE
C
      READ(10,901) NS
      READ(10,902) (TS(I),I=1,NS)
      READ(10,902) (CSS(I),I=1,NS)
      CLOSE(10)
C
C  Reassign the passed variables to b,sep,theta
c  This is where the random runs are created!
c
      b=ap
      sep=spac
      theta=porm
C
C
C  CALCULATE VELOCITY AND WRITE A NEW INPUT FILE USING LATIN DATA
C
C
c  ...ri=rate of injection (4 ml/min=5670ml/day)
c  ...fn=fracture porosity
c  ...pv=pore volume
c  ...cl is length of column
c  ...cv IS THE VOLUME OF THE COLUMN
      fn=b/sep
      pv=fn*cv
      v=ri*cl/pv
c  ...specify a filename for writing data...I chose craf.in
      OPEN(7,FILE='craf.in',STATUS='OLD')
C  ...Write data.
c  ...All parameters are fixed in start.cra except for 
c  ...V(calculated above),B,SEP,THETA,R,RP
      WRITE(7,903) C0
      WRITE(7,903) CIM
      WRITE(7,903) V
      WRITE(7,903) ALPHA
      WRITE(7,904) B
      WRITE(7,903) SEP
      WRITE(7,903) THETA
      WRITE(7,903) TAU
      WRITE(7,904) D0
      WRITE(7,903) R
      WRITE(7,903) RP
C
C  ...IF YOU DON'T WANT DECAY, THEN WRIT 0.0 FOR THALF
      WRITE(7,903) THALF
C
C  ...DECAY CONSTANT
      DLAMDA=0.0
      IF(THALF.GE.1.E-10) DLAMDA=ALOG(2.)/THALF
C  ...Writ the number of distance points along fracture.
      WRITE(7,901) NX
C  ...Writ the NX values of distance points along fracture.
        WRITE(7,902) (X(I),I=1,NX)
C  ...Writ the number of distance points into matrix
C  ...Note: If you only want to evaluate the fracture solution, then
C           read NZ=0 and then don't read ZM(I)
      WRITE(7,901) NZ
      NZ=NZ+1
C  ...Writ the NZ values of distance points in the matrix.
C  ...Note: The values should be in the range B < ZM(I) <= SEP
      IF(NZ.GT.1) THEN
        WRITE(7,902) (ZM(I),I=2,NZ)
      ENDIF
      WRITE(7,901) NT
C  ...Write the NT values of time points
        WRITE(7,902) (T(I),I=1,NT)
C
C  ...WRIT PARAMETER VALUES FOR CRUMP INVERSION
C  ...RECOMMEND: ALPHR=0.0, RELERR=1.E-6, N >= 11
C  ...WHERE 2N+1 IS THE NUMBER OF TERMS IN
C  ...THE FOURIER SERIES FOR NUMERICAL INVERSION OF THE
C  ...LAPLACE TRANSFORM
C  ...NOTE: IF THE SOLUTION OSCILLATES ABOUT THE STEADY STATE
C           CONCENTRATION, TRY RELERR=1.E-5 OR 1.E-4 ETC. OR
C           INCREASE N.
C
      WRITE(7,903) ALPHR
      WRITE(7,904) RELERR
      WRITE(7,901) N
C
C  ...WRIT NUMBER OF POINTS, TIME AND CONCENTRATION OF 
C  ...VARIABLE SOURCE
C
      WRITE(7,901) NS
      WRITE(7,902) (TS(I),I=1,NS)
      WRITE(7,902) (CSS(I),I=1,NS)
      CLOSE(7)
      end
c
c  ------------------------------------------------------------
c
       SUBROUTINE CRAFLUSH 
c     Commented out the write statements and modified the
c     write to keep writing to the same file for multiple runs
c     a selected amount of output (see lines around 295).
c  
c
C         **************************************
C         *                                    *
C         *              CRAFLUSH              *
C         *                                    *
C         **************************************
C         
C               SOLVED AND PROGRAMMED BY:
C                     E.A. SUDICKY
C        WATERLOO CENTRE FOR GOUNDWATER RESEARCH
C                UNIVERSITY OF WATERLOO
C              WATERLOO, ONTARIO, N2L 3G1
C         
C              COPYRIGHT 1988 E.A. SUDICKY
C         
C                        WARNING
C        
C         DUPLICATION OF THIS PROGRAM OR ANY PART
C         THEREOF WITHOUT THE EXPRESS PERMISSION
C         OF E.A. SUDICKY OF THE INSTITUTE FOR
C         GROUNDWATER RESEARCH IS STRICTLY FORBIDDEN
C
C         DISCLAIMER: WHILE CARE HAS BEEN EXERCISED TO ENSURE
C                     THE ACCURACY OF THE CODING, THE AUTHOR IS
C                     NOT RESPONSIBLE FOR INADVERTENT ERRORS.
C         
C          
C         TRANSPORT IN A SYSTEM OF PARALLEL FRACTURES
C         WITH MATRIX DIFFUSION
C         D > 0 SOLUTION
C         THIS VERSION ALLOWS THE SIMULATION OF THE FLUSHING
C         OF THE INITIALLY CONTAMINATED MATRIX (AND FRACTURES)
C         WITH A SOLUTION OF LOWER (OR HIGHER) CONCENTRATION
C         ENTERING THE FRACTURES.
C
C         THIS VERSION USES THE CRUMP ALGORITHM TO NUMERICALLY
C         INVERT THE LAPLACE TRANSFORMED SOLUTION.
C     
C     DEFINITION OF VARIABLES TO ENTER IN DATAFILE.
C
C     NOTE: USE CONSISTENT UNITS FOR INPUT DATA
C
C         C0= SOURCE CONCENTRATION AT FRACTURE ORIGIN
C         NOTE:  WHEN USING VARIABLE SOURCE LET C0=0
C         CIM=INITIAL CONCENTRATION IN MATRIX (AND IN FRACTURES)
C         V= VELOCITY IN FRACTURE
C         ALPHA=LONGITUDINAL DISPERSIVITY IN FRACTURE
C         B= FRACTURE APERTURE
C         SEP= FRACTURE SPACING
C         THETA= MATRIX POROSITY
C         NOTE: THETA=0.0 IS ALLLOWED, IN WHICH CASE A SOLUTION
C               WITH NO MATRIX DIFFUSION WILL BE OBTAINED
C         TAU= MATRIX TORTUOSITY
C         D0= DIFFUSION COEFF IN WATER
C         R= FRACTURE RETARDATION FACTOR    
C         RP= MATRIX RETARDATION FACTOR
C         THALF= HALF-LIFE
C
C         NX= NUMBER OF DISTANCE POINTS ALONG FRACTURE
C         X= DISTANCE ALONG FRACTURE
C         NZ= NUMBER OF DISTANCE POINTS IN MATRIX
C         ZM= DISTANCE IN MATRIX (B < ZM <= SEP)
C         NT= NUMBER OF TIME POINTS
C         T= TIME
C         NS= NUMBER OF POINTS USED TO DEFINE VARIABLE SOURCE 
C         TS= TIMES USED TO DEFINE VARIABLE SOURCE
C         CSS= CONCENTRATIONS USED IN VARIABLE SOURCE AT TS
C         2N+1=NUMBER OF TERMS USED IN CRUMP INVERSION
C         NOTE: THE CODE READS "N"
C
      COMPLEX ZNU,ZKAPPA,ZA,ZSIG,ZDECAY,RPDPZ,F(103),P,FS
      COMPLEX C0BC
      DOUBLE PRECISION Q(208),FT,C,CS
      CHARACTER*32 FDATA,FLIST
      DIMENSION X(51),ZM(32),T(150),CSS(51),TS(51)
      DATA Z/0.0/
C
C  ...Dimensioning hints for increasing array sizes for larger problems
C
C         X(NX),Z(NZ+1),T(NT)
C         F(MAXTRM),Q(MAXQ)
C     WHERE: MAXTRM = 2N+1 = NUMBER OF TERMS IN FOURIER SERIES
C                            FOR NUMERICAL INVERSION OF THE
C                            LAPLACE-TRANSFORMED SOLUTION
C            MAXQ   = 4N+4
C
C     FOR PASSING OF ARRAYS F AND Q...
C     THE VALUES BELOW PERMIT A MAXIMUM CRUMP "N" EQUAL TO 51
      MAXTRM=103
      MAXQ=208
C
C      WRITE(*,*) ' '
C      WRITE(*,*) ' ************************************************'
C      WRITE(*,*) ' *                MODEL: CRAFLUSH               *'
C      WRITE(*,*) ' *           PARALLEL CRACK MODEL D > 0         *'
C      WRITE(*,*) ' *            (c) E.A. SUDICKY, 1988            *'
C      WRITE(*,*) ' *   WATERLOO CENTRE FOR GROUNDWATER RESEARCH   *'
C      WRITE(*,*) ' *            UNIVERSITY OF WATERLOO            *'
C      WRITE(*,*) ' *       WATERLOO, ONTARIO, CANADA N2L 3G1      *'
C      WRITE(*,*) ' ************************************************'
C      WRITE(*,*) ' '
C      WRITE(*,*) ' ENTER INPUT DATAFILE NAME...'
C      READ(*,4682) FDATA
C4682  FORMAT(A32)
       OPEN(7,FILE='craf.in',STATUS='OLD')
C      WRITE(*,*) ' GIVE FILENAME FOR LISTING...'
C      READ(*,4682) FLIST
C  ...Read data from craf.in.
      READ(7,903) C0
      READ(7,903) CIM
      READ(7,903) V
      READ(7,903) ALPHA
903   FORMAT(F10.4)
      READ(7,904) B
      READ(7,903) SEP
      READ(7,903) THETA
      READ(7,903) TAU
      READ(7,904) D0
904   FORMAT(E10.4)
      READ(7,903) R
      READ(7,903) RP
C
C  ...IF YOU DON'T WANT DECAY, THEN READ 0.0 FOR THALF
      READ(7,903) THALF
C
C  ...DECAY CONSTANT
      DLAMDA=0.0
      IF(THALF.GE.1.E-10) DLAMDA=ALOG(2.)/THALF
C  ...Read the number of distance points along fracture.
      READ(7,901) NX
901   FORMAT(I5)
C  ...Read the NX values of distance points along fracture.
        READ(7,902) (X(I),I=1,NX)
902     FORMAT(5F10.3)
C  ...Read the number of distance points into matrix
C  ...Note: If you only want to evaluate the fracture solution, then
C           read NZ=0 and then don't read ZM(I)
      READ(7,901) NZ
      NZ=NZ+1
C  ...Read the NZ values of distance points in the matrix.
C  ...Note: The values should be in the range B < ZM(I) <= SEP
      IF(NZ.GT.1) THEN
        READ(7,902) (ZM(I),I=2,NZ)
      ENDIF
      READ(7,901) NT
C  ...Read the NT values of time points.
        READ(7,902) (T(I),I=1,NT)
C
C  ...READ PARAMETER VALUES FOR CRUMP INVERSION
C  ...RECOMMEND: ALPHR=0.0, RELERR=1.E-6, N >= 11
C  ...WHERE 2N+1 IS THE NUMBER OF TERMS IN
C  ...THE FOURIER SERIES FOR NUMERICAL INVERSION OF THE
C  ...LAPLACE TRANSFORM
C  ...NOTE: IF THE SOLUTION OSCILLATES ABOUT THE STEADY STATE
C           CONCENTRATION, TRY RELERR=1.E-5 OR 1.E-4 ETC. OR
C           INCREASE N.
C
      READ(7,903) ALPHR
      READ(7,904) RELERR
      READ(7,901) N
C
C  ...READ NUMBER OF POINTS, TIME AND CONCENTRATION OF 
C  ...VARIABLE SOURCE
C
      READ(7,901) NS
      READ(7,902) (TS(I),I=1,NS)
      READ(7,902) (CSS(I),I=1,NS)
      CLOSE(7)
C
      APER=B
      SPCE=SEP
C  ...Write values of variables to screen.
C
C  NOTE: THIS FORMAT DOES NOT WORK IN UNIX FOR SCREEN WRITE 
C
C      WRITE(*,8001) C0
 8001 FORMAT('  SOURCE CONCENTRATION AT INFINITE TIME= ',F10.4)
C      WRITE(*,8009) CIM
 8009 FORMAT('  INITIAL CONCENTRATION IN MATRIX AND FRACTURES= ',F10.4)
C      WRITE(*,8000) V
 8000 FORMAT('  VELOCITY IN FRACTURE= ',F10.3)
C      WRITE(*,8002) ALPHA
 8002 FORMAT('  FRACTURE DISPERSIVITY= ',F10.3)
C      WRITE(*,8005) APER
 8005 FORMAT('  FRACTURE APERTURE= ',E10.4)
C      WRITE(*,8010) SPCE
 8010 FORMAT('  FRACTURE SPACING= ',F10.3)
C      WRITE(*,8015) THETA
 8015 FORMAT('  MATRIX POROSITY= ',F10.3)
C      WRITE(*,8025) TAU
 8025 FORMAT('  MATRIX TORTUOSITY= ',F10.3)
C      WRITE(*,8030) D0
 8030 FORMAT('  DIFF. COEFF. IN WATER= ',E10.4)
C      WRITE(*,8035) R
 8035 FORMAT('  FRACTURE RETARDATION FACTOR= ',F10.3)
C      WRITE(*,8040) RP 
 8040 FORMAT('  MATRIX RETARDATION FACTOR= ',F10.3)
C      IF(THALF.LT.1.E-10) WRITE(*,8044)
 8044 FORMAT('  NON-DECAYING SOLUTE')
C      IF(THALF.GE.1.E-10) WRITE(*,8045) THALF
 8045 FORMAT('  HALF-LIFE= ',E10.3)
C      WRITE(*,8046) ALPHR,RELERR,N
 8046 FORMAT(//,'  NUMERICAL INVERSION PARAMETERS USED IN CRUMP',//,
     1'  ALPHR= ',F6.3,/,'  RELERR= ',E12.4,/,'  N= ',I3)
C      WRITE(*,8950)
 8950 FORMAT(//,'  VARIABLE SOURCE CONCENTRATION',//,
     1'    TIME  ',5X,'CONCENTRATION')
C      DO 8955 II=1,NS
C      WRITE(*,8951) TS(II),CSS(II)
C 8951 FORMAT(1X,F10.5,5X,F10.5)
C 8955 CONTINUE
C      WRITE(*,8200)
C      WRITE(*,8210)
C      WRITE(6,8041)
8041  FORMAT(/20X,'************************************************',
     &/20X,'*                MODEL: CRAFLUSH               *',/20X,
     &'*           PARALLEL CRACK MODEL D > 0         *',/20X,
     &'*            (c) E.A. SUDICKY, 1988            *',/20X,
     &'*   WATERLOO CENTRE FOR GROUNDWATER RESEARCH   *',/20X,
     &'*            UNIVERSITY OF WATERLOO            *',/20X,
     &'*       WATERLOO, ONTARIO, CANADA N2L 3G1      *',/20X,
     &'************************************************',//)
C      WRITE(6,8001) C0
C      WRITE(6,8009) CIM
C      WRITE(6,8000) V
C      WRITE(6,8002) ALPHA
C      WRITE(6,8005) APER
C      WRITE(6,8010) SPCE
C      WRITE(6,8015) THETA
C      WRITE(6,8025) TAU
C      WRITE(6,8030) D0
C      WRITE(6,8035) R
C      WRITE(6,8040) RP 
C      IF(THALF.LT.1.E-10) WRITE(6,8044)
C      IF(THALF.GE.1.E-10) WRITE(6,8045) THALF
C      WRITE(6,8046) ALPHR,RELERR,N
C      WRITE(6,8950)
      DO 8956 II=1,NS
C      WRITE(6,8951) TS(II),CSS(II)
8956  CONTINUE
C      WRITE(6,8200)
8200  FORMAT(//'    TIME      DISTANCE: ALONG FRACTURE,  IN MATRIX     C
     &(t)     STEADY C')
C      WRITE(6,8210)
8210  FORMAT('    ----      --------- --------------   ---------     ---
     &-     --------')
      B=B/2.
      ZM(1)=B
      SEP=SEP/2.
      DP=D0*TAU
      D=ALPHA*V+D0
C
C  ...CALCULATE PARAMETERS USED IN THE LAPLACE TRANSFORMED
C  ...SOLUTION L[C]. THEY SHOULD BE STORED AS COMPLEX VALUES
      ZA=CMPLX(THETA*SQRT(RP*DP)/(B*R),Z)
      ZKAPPA=CMPLX(4.0*R*D/V/V,Z)
      RPDP=SQRT(RP/DP)
      ZSIG=CMPLX(RPDP*(SEP-B),Z)
      ZDECAY=CMPLX(DLAMDA,Z)
      VV=V/2./D
      DO 44 I=1,NT
C
C  ...INITIALIZATION FOR CRUMP
      TT=T(I)
      CALL INITC(TT,RELERR,ALPHR,ATERM,TFACT)
      DO 55 J=1,NX
      VVX=VV*X(J)
      ZNU=CMPLX(VVX,Z)
      DO 555 K=1,NZ
      IF(X(J).LT.1.E-10) GOTO 8043
      RPDPZ=CMPLX(RPDP*(SEP-ZM(K)),Z)
C
C     STEADY-STATE SOLUTION GIVEN BY:
C     LIMIT AS P --> 0 OF P*L[C]
C     OR C0 IF NO DECAY
      IF(THALF.GE.1.E-10) THEN
      P=CMPLX(Z,Z)
      CS=REAL(FS(P,ZNU,ZKAPPA,ZA,ZSIG,ZDECAY,RPDPZ,C0,CIM,CSS,TS,NS))
      ELSE
      CS=C0
      END IF
C
C     TRANSIENT SOLUTION
      CALL CRUMP(FT,TT,ATERM,TFACT,ZNU,ZKAPPA,ZA,ZSIG,ZDECAY,RPDPZ,F,Q,
     &N,MAXTRM,MAXQ,C0,CIM,CSS,TS,NS)
      C=FT
      GOTO 8042
8043  IF(ZM(K).LE.B) THEN
      C=C0
      CS=C0
      ENDIF
c
c   write to conc.out 
c
8042  WRITE(8,*) T(I),C
c
c
C 8042  WRITE(6,8060) T(I),X(J),ZM(K),C,CS
C      WRITE(*,8060) T(I),X(J),ZM(K),C,CS
C      WRITE(*,8060) T(I),X(J),ZM(K),C,CS
c 8060  FORMAT(1X,I5,F10.2,12X,F8.3)
555   CONTINUE

55    CONTINUE
44    CONTINUE
      CLOSE(6)
      return
      END
      COMPLEX FUNCTION FS(P,ZNU,ZKAPPA,ZA,ZSIG,ZDECAY,RPDPZ,C0,CIM,CSS
     &,TS,NS)
      COMPLEX P,ZNU,ZKAPPA,ZA,ZSIG,ZDECAY,RPDPZ,CTANH,COSHXY,ARGE
      COMPLEX SP,PL,CON1,PLS,PRPDPZ,PTH,ATC,ZEXP
      COMPLEX BTC,ABTC,C0BC
      DIMENSION CSS(51),TS(51)
      DATA Z/0.0/
C
C     FS=L[C/C0]
C
      CON1=CMPLX(1.0,0.0)
      PL=P+ZDECAY
      PLS=CSQRT(PL)
      SP=ZSIG*PLS
      PRPDPZ=RPDPZ*PLS
      PTH=PLS*CTANH(SP)
      ATC=ZA*PTH+PL
C     BTC=(ZA+PTH)/PTH
C     ABTC=BTC/ATC
      ARGE=ZNU*(CON1-CSQRT(CON1+ZKAPPA*ATC))
      ARGER=REAL(ARGE)
      ARGEI=AIMAG(ARGE)
      ZEXP=CMPLX(Z,Z)
      IF(ARGER.GT.-75.0.AND.ARGEI.GT.-75.0) THEN
      ZEXP=CEXP(ARGE)
      END IF
      IF(ARGER.GT.-75.0.AND.ARGEI.LE.-75.0) THEN
      ARGE=CMPLX(ARGER,-75.0)
      ZEXP=CEXP(ARGE)
      END IF
      IF(ARGER.LE.-75.0.AND.ARGEI.GT.-75.0) THEN
      ARGE=CMPLX(-75.0,ARGEI)
      ZEXP=CEXP(ARGE)
      END IF
      IF(CABS(P).GT.1.E-20) THEN
      FS=CIM/PL+(C0BC(P,CSS,TS,NS)-CIM/PL)*ZEXP
      ELSE
      FS=C0*ZEXP
      END IF
C
C     CHECK FOR VERY SMALL FS; SET TO 0.0 IF REAL AND IMAGINARY
C     PARTS LESS THEN 1.E-20
      CHECKR=REAL(FS)
      CHECKI=AIMAG(FS)
      ACHR=ABS(CHECKR)
      ACHI=ABS(CHECKI)
      IF(ACHR.LT.1.E-20.AND.ACHI.LT.1.E-20) FS=CMPLX(Z,Z)
      IF(ACHR.GE.1.E-20.AND.ACHI.LT.1.E-20) FS=CMPLX(CHECKR,Z)
      IF(ACHR.LT.1.E-20.AND.ACHI.GE.1.E-20) FS=CMPLX(Z,CHECKI)
      IF(CABS(P).GT.1.E-20) THEN
      FS=(FS-CIM/PL)*COSHXY(PRPDPZ,SP)+CIM/PL
      ELSE
      FS=FS*COSHXY(PRPDPZ,SP)
      ENDIF
      RETURN
      END
C
C---------------------------------------------------------------------
C     FUNCTION TO EVALUATE HYPERBOLIC TANGENT
      COMPLEX FUNCTION CTANH(X)
      COMPLEX X,A1,A2
      REAL XR,XI,AXR,AXI
      U=CABS(X)
      XR=REAL(X)
      XI=AIMAG(X)
      IF(U.GT.1.E-06.AND.U.LT.15.) THEN
      A1=CEXP(X)
      A2=CEXP(-X)
      CTANH=(A1-A2)/(A1+A2)
      END IF
C
C     FOR LARGE X, CTANH(X) -> 1.
      IF(U.GE.15.) THEN
      CTANH=CMPLX(1.,0.)
      END IF
C
C     FOR SMALL X, CTANH(X) -> X-X*X*X/3.
      IF(U.LE.1.E-06) THEN
      CTANH=CMPLX(0.0,0.0)
      AXR=ABS(XR)
      AXI=ABS(XI)
      IF(AXR.GT.1.E-10.AND.AXI.GT.1.E-10) CTANH=X-X*X*X/3.
      IF(AXR.LE.1.E-10.AND.AXI.GT.1.E-10) THEN
      X=CMPLX(0.0,XI)
      CTANH=X-X*X*X/3.
      END IF
      IF(AXR.GT.1.E-10.AND.AXI.LE.1.E-10) THEN
      X=CMPLX(XR,0.0)
      CTANH=X-X*X*X/3.
      END IF
      END IF
      RETURN
      END
C
C---------------------------------------------------------------------
C     FUNCTION TO EVALUATE COSH(X)/COSH(Y) WITH X, Y COMPLEX
      COMPLEX FUNCTION COSHXY(X,Y)
      COMPLEX X,Y,PART1,PART2,PART3,PART4
      DOUBLE PRECISION A,B,C,D,E1,E2,E3,E4,AC,C2,CSB,SSB,CSD,SSD
      A=DBLE(REAL(X))
      B=DBLE(AIMAG(X))
      C=DBLE(REAL(Y))
      D=DBLE(AIMAG(Y))
      CSB=DCOS(B)
      SSB=DSIN(B)
      CSD=DCOS(D)
      SSD=DSIN(D)
      PART1=CMPLX(SNGL(CSB),SNGL(SSB))
      PART2=CMPLX(SNGL(CSB),-SNGL(SSB))
      PART3=CMPLX(SNGL(CSD),SNGL(SSD))
      PART4=CMPLX(SNGL(CSD),-SNGL(SSD))
      E1=DEXP(A-C)
      AC=A+C
      IF(AC.LT.30.D0) THEN
      E2=DEXP(-A-C)
      ELSE
      E2=0.D0
      ENDIF
      E3=1.D0
      C2=2.D0*C
      IF(C2.LT.30.D0) THEN
      E4=DEXP(-C2)
      ELSE
      E4=0.D0
      ENDIF
      ES1=SNGL(E1)
      ES2=SNGL(E2)
      ES3=SNGL(E3)
      ES4=SNGL(E4)
      COSHXY=(ES1*PART1+ES2*PART2)/(ES3*PART3+ES4*PART4)
      RETURN
      END
C
C---------------------------------------------------------------------
C        ROUTINE FOR NUMERICAL INVERSION OF LAPLACE
C        TRANSFORMATIONS USING THE CRUMP ALGORITHM.
C        THIS USES THE 'EPAL' METHOD OF SERIES
C        SUMMATION AS SUGGESTED IN THE CRUMP PAPER.
C
C        REFERENCE:  KENNY S. CRUMP, 'NUMERICAL
C        INVERSION OF LAPLACE TRANSFORMS USING A
C        FOURIER SERIES APPROXIMATION', JOURNAL OF
C        THE ASSOCIATION FOR COMPUTING MACHINERY
C        VOL 23, NO. 1, JAN 1976, PP. 89-96.
C
C        PROGRAMMED BY:  FRANK LETNIOWSKI
C                        INSTITUTE FOR GROUNDWATER RESEARCH
C                        UNIVERSITY OF WATERLOO
C
C                        MAY, 1987
C
C
C        THERE ARE TWO PARAMETERS NEEDED FOR THE CRUMP
C        ROUTINE.  THEY ARE RELATED TO THE MAXIMUM TIME
C        USED IN THE INVERSION, THE RELATIVE ERROR REQUIRED,
C        AND THE LARGEST REAL VALUE OF THE POLES OF THE 
C        TRANSFORMED FUNCTION.  'INITC' CLACULATES THESE
C        VALUES GIVEN THE ABOVE INFORMATION.
C
      SUBROUTINE INITC(TMAX,RELERR,ALPHA,ATERM,TFACT)
C                -----
C     TMAX = THE MAXIMUM TIME NEEDED IN THE INVERSION
C            (ie. USE YOUR LARGEST TIME VALUE)
C     RELERR = THE MAXIMUM RELATIVE ERROR REQUIRED
C            (1E-6 IS A REASONABLE VALUE FOR MOST APPLICATIONS)
C     ALPHA = THE APPROXIMATE EXPONENTIAL ORDER OF THE
C             UNTRANSFORMED FUNCTION.  THIS IS APPROX.
C             THE LARGEST REAL VALUE OF THE POLES OF THE
C             TRANSFORMED FUNCTION.
C             (ALPHA=0.0 IS A REASONABLE VALUE FOR MOST APPLICATIONS)
C     A AND TFACT ARE CALCULATED HERE AS SUGGESTED
C     BY THE CRUMP PAPER.
C
         TFACT = 0.7931*TMAX
         ATERM = ALPHA - ALOG(RELERR)/(2.0*TFACT)
C
      RETURN
      END
C
C----------------------------------------------------------------
C     THE CRUMP SUBROUTINE NUMERICALLY INVERTS A LAPLACE
C     TRANSFORMED FUNCTION FS AT A PARTICULAR TIME T
C     RETURNING THE VALUE IN FT.  FS MUST BE A COMPLEX 
C     FUNCTION OF A COMPLEX ARGUMENT. HERE, X IS THE LAPLACE
C     VARIABLE 'P' WHICH IS COMPLEX. FOR MAXIMUM EFFICIENCY,
C     YOU MAY WANT TO STORE THE FS VALUES FOR EACH X IN AN
C     ARRAY (ie. PRECOMPUTE THEM (ALONG WITH THE NEEDED 
C     X-VALUES) OUTSIDE YOUR TIME LOOP AND THEN SUBSTITUTE
C     THE FS ARRAY INTO THE CRUMP FOURIER SERIES THAT INVOLVES
C     TIME. THAT IS, IT IS NOT NECESSARY TO COMPUTE FS(X)
C     FOR EACH VALUE OF THE TIME 'T' SINCE X DEPENDS ONLY
C     ON TFACT THROUGH TMAX.
C
C     THE PARAMETERS A AND TFACT SHOULD BE AS SUGGESTED
C     ABOVE OR IN THE PAPER.  N IS RELATED TO THE NUMBER 
C     OF TERMS USED IN THE SUMMATION ie. 2N+1 TERMS ARE
C     CALCULATED IN THE 'EPAL' METHOD OF SUMMATION.
C
C
C     DIMENSIONING HINTS:
C                        F(MAXTRM),Q(MAXQ)
C     WHERE: MAXTRM = 2N+1
C            MAXQ   = 4N+4
      SUBROUTINE CRUMP(FT,T,ATERM,TFACT,ZNU,ZKAPPA,ZA,ZSIG,ZDECAY,
     &RPDPZ,F,Q,N,MAXTRM,MAXQ,C0,CIM,CSS,TS,NS)
C
         DOUBLE PRECISION SUM,SER,Q(MAXQ),P0,P1,TEMP,
     1                    TOT,F0,FT
         INTEGER M,J,K,MAX,CUR,FIN
         COMPLEX FS,F(MAXTRM),X,ZNU,ZKAPPA,ZA,ZSIG,ZDECAY,RPDPZ
         COMPLEX C0BC
         DIMENSION CSS(51),TS(51)
C
         DATA PI/3.14159265/
C
         M = 2*N+1
         X = CMPLX(ATERM,0.0)
         F0 = DBLE(FS(X,ZNU,ZKAPPA,ZA,ZSIG,ZDECAY,RPDPZ,C0,CIM,CSS,TS
     &  ,NS))
C       CALCULATE THE FUNCION VALUES NEEDED
         DO 800 K=1,M
            X=CMPLX(ATERM,REAL(K)*PI/TFACT)
            F(K)=FS(X,ZNU,ZKAPPA,ZA,ZSIG,ZDECAY,RPDPZ,C0,CIM,CSS,TS,NS)
800      CONTINUE
C
C        INITIALISE THE QUEUE NEEDED FOR CALCULATION OF 
C        THE SUM.  MAX IS THE SIZE OF THE QUEUE NEEDED.
C        THE QUEUE STORES EPSILON VALUES NEEDED IN THE
C        EPAL METHOD TO CALCULATE THE SERIES.
         MAX=2*M+2     
         CALL INITQ(Q,CUR,FIN,MAX,MAXQ)
C
C        CALCULATE THE SERIES BY USE OF THE EPAL METHOD.
C          SER(I,F(I),TEMP) CALCULATES THE Ith TERM IN
C          THE SERIES.
         TEMP=DBLE(PI*T/TFACT)
         SUM=SER(1,F(1),TEMP)
C        INSERT INTO THE QUEUE
         CALL INSERT(SUM,Q,CUR,FIN,MAX,MAXQ)
C
         DO 810 K=2,M
C           SUM IS THE Kth PARTIAL SUM OF THE SERIES
            SUM = SUM + SER(K,F(K),DBLE(K)*TEMP)
C           CALCULATE EPS SUB K SUP (1) WHICH IS TOT
C           UPON LEAVING THE J LOOP
            TOT = SUM
            CALL INSERT(TOT,Q,CUR,FIN,MAX,MAXQ)
            P1 = DBLE(0.0)
            DO 820 J=2,K
               P0=P1
               CALL REMOVE(P1,Q,CUR,FIN,MAX,MAXQ)
               IF(DABS(TOT-P1).GT.0.0D0)THEN
                  TOT=P0+DBLE(1.0)/(TOT-P1)
               END IF
               CALL INSERT(TOT,Q,CUR,FIN,MAX,MAXQ)
820         CONTINUE
810      CONTINUE   
C
C        TOT IS NOW EPS SUB M SUP (1), THE REQUIRED
C        SERIES APPROXIMATION.
C
C        CALCULATE THE FINAL FUNCTION VALUE.
C
         FT=DEXP(DBLE(ATERM*T))/TFACT*(DBLE(0.5)*F0+TOT)
C
      RETURN
      END
C
C     FUNCTION TO CALCULATE THE Kth TERM IN THE SERIES.
C        FCALC IS THE FS(A+K*PI*i/TFACT)
C        FACT IS K*PI*T/TFACT
      DOUBLE PRECISION FUNCTION SER(K,FCALC,FACT)
C                               ---
         DOUBLE PRECISION FACT
         COMPLEX FCALC
         INTEGER K
C
         SER=DBLE(FCALC)*DCOS(FACT)-
     1       DBLE(AIMAG(FCALC))*DSIN(FACT)
C
      RETURN
      END
C
C     PROCEDURES TO HANDLE THE QUEUE PROCESSING
C
C     INITIALISE TO AN EMPTY QUEUE
      SUBROUTINE INITQ(Q,CUR,FIN,MAX,MAXQ)
C                -----     
         DOUBLE PRECISION Q(MAXQ)
         INTEGER CUR,FIN,MAX
C
         CUR=1
         FIN=1
C
      RETURN
      END
C
C     INSERT A VALUE INTO THE QUEUE
C        VALUES ARE INSERTED AT THE END OF THE LIST
      SUBROUTINE INSERT(VAL,Q,CUR,FIN,MAX,MAXQ)
C                ------
         DOUBLE PRECISION Q(MAXQ),VAL
         INTEGER CUR,FIN,MAX,NEXT
C
         NEXT=FIN+1
         IF(NEXT.GT.MAX)THEN
            NEXT=1
         END IF
C
         IF(NEXT.EQ.CUR)THEN
            CONTINUE
         ELSE
            Q(FIN)=VAL
            FIN=NEXT
         END IF
C
      RETURN
      END
C
C     SUBROUTINE TO GET A VALUE FROM THE QUEUE.
C        VALUES ARE REMOVED FROM THE FRONT OF THE LIST
      SUBROUTINE REMOVE(VAL,Q,CUR,FIN,MAX,MAXQ)
C                ------
         DOUBLE PRECISION Q(MAXQ),VAL
         INTEGER CUR,FIN,MAX
C
         IF(CUR.EQ.FIN)THEN
            CONTINUE
         ELSE
            VAL=Q(CUR)
            CUR=CUR+1
            IF(CUR.GT.MAX)THEN
               CUR=1
            END IF
         END IF
C
      RETURN
      END
C
C       FUNCTION TO EVALUATE TIME VARIABLE BOUNDARY CONDITION FOR PARENT
C
      COMPLEX FUNCTION C0BC(P,CSS,TS,NS)
      COMPLEX P,SUM
      DIMENSION CSS(51),TS(51)
      DATA Z /0.0/
      SUM=CMPLX(0.0,0.0)
      DO 50 I=2,NS
      SUM=SUM+CMPLX((CSS(I)-CSS(I-1)),Z)*CEXP(-P*CMPLX(TS(I),Z))
  50  CONTINUE
      C0BC=(CMPLX(CSS(1),Z)+SUM)/P
      RETURN
      END
c  ------------------------------------------------------------
c  Subroutine for calculating minimum error
c
      subroutine findmin(nobs,ii)
      dimension t(150),c(150),to(25),co(25)
      errsum=0
c  READ INTERP.DATA (read the interpolation data)
c  This file containes time, concentration data used
c  to represent the measured breakthrough curve.
c  Number of line in the file is given by nobs (screen input)
c
      open (9,file='interp.data',status='old')
      do 10 i=1,nobs
      read (9,*) to(i),co(i)       
10    continue
      rewind(8)
      open (8,file='conc.out',status='old')
      do 20 k=1,150
      read (8,*) t(k),c(k) 
20    continue
c
c find the minimum error in any direction
c

      do 40 i=1,nobs
      error=1.e10
      do 30 k=1,150
c
c Here you could select modeled times within a certain range
c of the observed data.  Not done in present version because
c it requires additional input (ie what is the range of to
c to bound the search?)
c 
c      if (t(k).lt.(to(i)-1) then
c          goto 30
c       endif
c      if (t(k).gt.(to(i)+1) then
c          goto 30
c       endif
c calculate hypotenuse
      s1=t(k)-to(i)
      s2=c(k)-co(i)
      dist=sqrt((s1**2)+(s2**2))
c
c find the smallest dist and put it into error
c
      if (dist.lt.error) then
          error=dist
      endif
30    continue
      open(22,file='error.out',status='unknown')
      write(22,*) i,error
c
c sum the errors
c
      errsum=errsum+error  
40    continue
      open(11,file='err.sum',status='unknown')
c
c  write the run no. and the error sum divided by nobs
c  (number of interpolation points)
c
      write(11,*)ii,errsum/nobs
      rewind(9)
      close(9)
      rewind(8)
      close(8)
      return
      end