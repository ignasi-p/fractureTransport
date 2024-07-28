# -*- coding: utf-8 -*-
"""
Script to reproduce the model described in:
   Muskus and Falta 2018. Semi-analytical method for matrix diffusion in heterogeneous
   and fractured systems with parent-daughter reactions.
   Journal of Contaminant Hydrology, 218, 94-109.
   https://doi.org/10.1016/j.jconhyd.2018.10.002.

This script will create a PHREEQC input file, using the parameters given in the script
itself, with the values given in the text and in Table 1 of Muskus and Falta (2018).
Variable 'modelType' decides the type of discretization used for the rock matrix.

The script will then run the selected model using [IPHREEQC](https://www.usgs.gov/software/phreeqc-version-3)
      Charlton and Parkhurst (2011) http://dx.doi.org/10.1016/j.cageo.2011.02.005
Except that if the input tile (*.pqi) and the corresponding selected output
file (*.tsv) exist and the output file is newer than the input file, then the
calculations are NOT performed, and the esisting tsv-file is read instead.
To force the calculation, just remove either the pqi- or the tsv-files.
The reason for not repeating calculations is that a PHREEQC run
may take quite a long time, about 5.5 hours for 150 cells on a 2.4 GHz processor.

All figures are saved in subfolder 'plots'.

The script will run under Windows.

@author: Ignasi Puigdomenech
"""
import os, sys
import shutil
import math
from win32com.client import Dispatch
import matplotlib.pyplot as plt
import pandas as pd
from time import time, strftime, localtime
sys.path.insert(0, '../Analytical_eqn/CraFlush')
from craflush import CraFlush
# The Phreeqc database
database = '../phreeqc_3.7.3.dat'

# The model is 3 m long. The value of 'n_cells' controls:
#    - Delta_x, the cell x-length
#    - Delta_t, the time step
#    - the number of shifts to reach the end time
# larger values give lower numerical dispersion,
# smaller values give shorter execution times,
# should be between 10 and perhaps 400
n_cells = 125

# The model type selects the type of rock-matrix model:
# - 0 (zero): 5 stagnant cells of equal size
# - 1 (one): 4 stagnant cells of equal size
# - 2 (two): 5 stagnant cells of increasing size
# - 3 (three): 4 stagnant cells of increasing size
# Increasing sizes require a larger value for 'n_cells'
# Four cells (instead of five) require shorter run time, but might be less accurate.
modelType = 1

# Boundary cell factor, see eqn.(127) in Phreeqc user's guide (Parkhurst and Appelo 1999)
fbc = 1 # it can be =1, =2 or =0 which means calculated according to eqn.(127)

# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------

def runPhreeqc(input_string):
    """Load database via linked library object and run input string.
    The results (selected_output) are in a "tuple" of tuples.
    """
    print("--- starting Phreeqc at "+strftime("%H:%M:%S", localtime()))
    t1 = time() # the start time to check how long does the optimisation take
    iphreeqc = Dispatch('IPhreeqcCOM.Object')
    iphreeqc.LoadDatabase(database)
    iphreeqc.OutputFileOn = True
    iphreeqc.SelectedOutputFileOn = True
    iphreeqc.SelectedOutputStringOn = True
    try:
        iphreeqc.RunString(input_string)
    except :
        errTxt = iphreeqc.GetErrorString()
        if len(errTxt) > 0: print("IPhreeqc.GetErrorString():\n"+errTxt)
        raise SystemExit
    # add a header so that Phreeqc Interactive can show the output file sections
    headStr = "   Input file: \n  Output file: \nDatabase file: \n\n"
    headStr = headStr + "------------------\nReading data base.\n------------------\n\n"
    with open('header_pqo','w') as f:
        f.writelines(headStr)
    # copy the two files together into 'pqo', remove the original files
    with open(pqo,'wb') as wfd:
        for f in ['header_pqo',iphreeqc.OutputFileName]:
            with open(f,'rb') as fd:
                shutil.copyfileobj(fd, wfd)
    os.remove('header_pqo')
    print("    (removing file",iphreeqc.OutputFileName,")")
    os.remove(iphreeqc.OutputFileName)

    # Get the selected_output results
    # The results (selected_output) are in a "tuple" of tuples.
    selOutput = iphreeqc.GetSelectedOutputArray()
    # Create a DataFrame
    res = pd.DataFrame(selOutput[1:])
    res.columns = selOutput[0]
    # remove rows containing the inflowing solution (at x=0)
    res = res[(res["dist_x"] > 0)]
    t2 = time() - t1
    minutes, seconds = divmod(t2, 60)
    hours, minutes = divmod(minutes, 60)
    print("    completed after %d:%02d:%02d" % (hours, minutes, seconds))
    return res

def pqi_tsv_ok(pqi,tsv):
    # find out if files 'pqi' and 'tsv' both exist and they are files (not directories)
    # and return ok if 'pqi' is older than 'tsv'
    ok = False
    exists_pqi = os.path.exists(pqi)
    if exists_pqi: exists_pqi = os.path.isfile(pqi)
    exists_tsv = False
    if exists_pqi:
        exists_tsv = os.path.exists(tsv)
        if exists_tsv: exists_tsv = os.path.isfile(tsv)    
    if exists_pqi and exists_tsv:
        ok = os.path.getmtime(pqi) < os.path.getmtime(tsv)
    return ok

def c_c0_S(tyr,X,Z,N,Rp):
    # Implements the analytical equation for matrix diffusion from Sudicky et al.
    # N decides the number of terms in the EPAL summation
    # units: x and z in 'm', tyr in 'years'
    t = tyr * 365.25*86400 # convert years to 's'
    # Data for the PhreeqC 'model':
    # b = fracture aperture
    C0 = 0.     # Source concentration at infinite time
    Cim = 0.    # Initial concentration in matrix and fractures
    alpha = alpha_num + alpha_L  # Fracture dispersivity (m)
    sep = 2.*z_tot   # Fracture spacing (m)
    theta = phi_im  # Matrix porosity
    D0 = D_0    # Diffusion coefficient in water (m2/s)
    if D0 <= 0: D0 = 1.e-9
    tau = D_im/D0    # Matrix tortuosity
    R = 1.      # Fracture retardation factor
    # Rp        # Matrix retardation factor
    THalf = 0.  # non-dcaying solute
    # Numerical inversion parameters used in CRUMP:
    AlphR = 0.   # alpha
    relErr = 1.e-6  # relatve error
    #Variable source concentration
    Ns = 2
    Ts =  [0.,1.585e9]
    Css = [1.,0.]
    # ---- the equations by Sudicky
    if(X <= 0): X=1e-3
    if(Z <= 0): Z=b/2.
    cc0, Cs = CraFlush(C0,Cim,v,alpha,b,sep,theta,tau,D0,R,Rp,THalf,
                             X,Z,t, Ns,Ts,Css,AlphR,relErr,N)
    if(cc0 > 1.): cc0 =1.
    if(cc0 < 0.): cc0 =0.
    return cc0

def annotTxt():
# The model type selects the type of rock-matrix model:
# - 0 (zero): 5 stagnant cells of equal size
# - 1 (one): 4 stagnant cells of equal size
# - 2 (two): 5 stagnant cells of increasing size
# - 3 (three): 4 stagnant cells of increasing size
    modelTxt = 'Rock matrix diffusion?$_{}$'
    if modelType == 0 or modelType == 1:
        modelTxt = str(n_z)+' equal stagnant layers$_{}$'
    elif modelType == 2 or modelType == 3:
        modelTxt = str(n_z)+' stagnant layers of\n    increasing size$_{}$'
    annotateTxt = (
        modelTxt+"\n"+
        "$n\mathsf{_{cells}}$ = "+str(n_cells)+"\n"+
        "$\\alpha_L$ = "+"{:.3f}".format(alpha_L)+" m\n"+
        "$D_{im}$ = "+"{:.2e}".format(D_im)+" m$^2$/s\n"+  # used in PHREEQC
        "$f_{bc}$ = "+"{:.2f}".format(f_bc_ij[1])+"\n"+
        "Markers: PHREEQC$_{ }$\n"+
        "Lines: analytical model\n  Sudicky et al.$^{ }$\n"+
        "using: $\\alpha_{tot}$ = "+"{:.3f}".format(alpha_num + alpha_L)+" m\n"+
        "$D'$ = {:.2e}".format(D_im)+" m$^2$/s\n"+
        "$R'$ = 2 for Na\n   and = 1 for NO$_3$")
    return annotateTxt

def plot(rsl,include):
    # Plot results.  'rsl' is the DataFrame with the data to plot.
    # The figures are saved in folder 'plots'
    # Do not plot all data, plots only symbols for every 'include'
    if include < 1: include = 1
    print("--- making plots along fracture and matrix.")

    # Two figures will be created:
    #    1- concentrations along the fracture at different times
    #    2- concentrations into the rock matrix for different times
    #       and at three distances from the inflow

    # 1- First create a figure with
    # concentrations along the fracture at different times

    fig, ax=plt.subplots(nrows=1,ncols=2,layout='constrained',figsize=(7.087,3.543),dpi=300)

    # colors from
    # https://davidmathlogic.com/colorblind/#%23000000-%23E69F00-%2356B4E9-%23009E73-%23F0E442-%230072B2-%23D55E00-%23CC79A7
    clr1 = "#88CCEE"  # light blue
    clr2 = "#44AA99"  # green
    clr3 = "#117733"  # dark green
    clr4 = "#332288"  # dark purple
    # markers: see the Notes section in
    # https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html
    mrk1 = "v"  # triangle down
    mrk2 = "D"  # diamond
    mrk3 = "^"  # triangle up
    mrk4 = "s"  # square

    # t_1, t_2 ... are times at which to plot results, for example, 50, 100 ... years
    Yr_1 = rsl[(rsl["Years"]==t_1) & (rsl["dist_x"]>0) & (rsl["dist_x"]<L) & (rsl["soln"]<=n_cells)].iloc[::include,:]
    Yr_2 = rsl[(rsl["Years"]==t_2) & (rsl["dist_x"]>0) & (rsl["dist_x"]<L) & (rsl["soln"]<=n_cells)].iloc[::include,:]
    Yr_3 = rsl[(rsl["Years"]==t_3) & (rsl["dist_x"]>0) & (rsl["dist_x"]<L) & (rsl["soln"]<=n_cells)].iloc[::include,:]
    Yr_4 = rsl[(rsl["Years"]==t_4) & (rsl["dist_x"]>0) & (rsl["dist_x"]<L) & (rsl["soln"]<=n_cells)].iloc[::include,:]

    xLabel = 'Distance x along fracture (m)'

    elements = ['Na','NO3']
    for i in range(2):
        # The retardation factor in the rock matrix
        ele = elements[i]
        Rp = 2. if (ele == 'Na') else 1.  # Matrix retardation factor

        #ax.axis([0, L, 0, 1.])
        colmn = ele+'_fr'
        ax[i].plot(Yr_1["dist_x"],Yr_1[colmn],c=clr1,marker=mrk1,markersize=9,linestyle='None',label="{:.0f}".format(t_1)+" yr")
        ax[i].plot(Yr_2["dist_x"],Yr_2[colmn],c=clr2,marker=mrk2,markersize=6,linestyle='None',label="{:.0f}".format(t_2)+" yr")
        ax[i].plot(Yr_3["dist_x"],Yr_3[colmn],c=clr3,marker=mrk3,markersize=9,linestyle='None',label="{:.0f}".format(t_3)+" yr")
        ax[i].plot(Yr_4["dist_x"],Yr_4[colmn],c=clr4,marker=mrk4,markersize=6,linestyle='None',label="{:.0f}".format(t_4)+" yr")

        yLabel = ele+'/'+ele+'$_0$'
        ax[i].set_ylabel(yLabel)
        ax[i].set_xlabel(xLabel)

        # calculate and plot the analytical (theoretical) curves
        distx = [0]*51
        y1 = [1]*51
        y2 = [0]*51
        y3 = [0]*51
        y4 = [0]*51
        for j in range(1,51):
            distx[j] = L*j/50
            y1[j] = c_c0_S(t_1,distx[j],0,12,Rp)
            y2[j] = c_c0_S(t_2,distx[j],0,12,Rp)
            y3[j] = c_c0_S(t_3,distx[j],0,12,Rp)
            y4[j] = c_c0_S(t_4,distx[j],0,12,Rp)
        ax[i].plot(distx,y1,c=clr1,marker="None")
        ax[i].plot(distx,y2,c=clr2,marker="None")
        ax[i].plot(distx,y3,c=clr3,marker="None")
        ax[i].plot(distx,y4,c=clr4,marker="None")

        if i == 1:
            ax[i].legend(fontsize='small', bbox_to_anchor=(1.1,0.), loc='lower left', borderaxespad=0)
            # write an explanation text
            txt = annotTxt()
            ax[i].annotate(txt, xy=(1.07, 1.), xycoords='axes fraction',
                           fontsize=8, verticalalignment='top')

    plotName=fileName+"-1.png"
    plt.savefig('plots/'+plotName)
    plt.show()
    plt.close()

    # 2- Create a figure containing three subfigures in columns.
    # - in each of the two subfigures at the left, three plots
    #   ('axes') showing concentrations versus distance into the matrix
    # - to the right a text box with the parameters for the simulation
    #   and a legend (the same one for all plot)

    fig = plt.figure(layout='constrained',figsize=(7.087,3.543),dpi=300)
    subfigs = fig.subfigures(nrows=1,ncols=3, width_ratios=[2,2,1])

    clr1 = "#DC267F" # mauve
    clr2 = "#FE6100"  # orange
    clr3 = "#FFB000"  # yellow
    mrk1 = "v"  # triangle down
    mrk2 = "D"  # diamond
    mrk3 = "^" # triangle up
    xLabel = 'Distance z in matrix (m)'

    # Left subfigures: in each three plots (axes) in a single column
    for i in range(2):
        # The retardation factor in the rock matrix
        ele = elements[i]
        Rp = 2. if (ele == 'Na') else 1.  # Matrix retardation factor

        axsLeft = subfigs[i].subplots(3, 1, sharex=True)

        # Create a plot with concentrations versus distance into the matrix for
        # three distances along the fracture and three different
        # simulation times (t_1, t_2 and t_3).
        yLabel = ele+'/'+ele+'$_0$'

        # t_1, t_2 ... are times at which to plot results, for example, 50, 100 ... years
        years = [t_1,t_2,t_3]
        for j in range(0,3):        
            tYr = years[j]
            xMn=rsl[(rsl["Years"]==tYr) & (rsl["dist_x"]==xmn)]
            xMd=rsl[(rsl["Years"]==tYr) & (rsl["dist_x"]==xmd)]
            xMx=rsl[(rsl["Years"]==tYr) & (rsl["dist_x"]==xmx)]
            colmn = ele+'_fr'
            axsLeft[j].plot(xMn['dist_z'],xMn[colmn],c=clr1,marker=mrk1,linestyle='None',markersize=9)
            axsLeft[j].plot(xMd['dist_z'],xMd[colmn],c=clr2,marker=mrk2,linestyle='None',markersize=6)
            axsLeft[j].plot(xMx['dist_z'],xMx[colmn],c=clr3,marker=mrk3,linestyle='None',markersize=9)
            axsLeft[j].set_ylim([-0.1, 1.1])
            axsLeft[j].set_ylabel(yLabel)
            axsLeft[j].set_title("{:.0f}".format(tYr)+" yr", x=0.7, y=0.63, fontsize=10)
            # calculate and plot the analytical (theoretical) curves
            distz = [0]*21
            xMnA = [1]*21;  xMdA = [1]*21;  xMxA = [1]*21
            N = 11
            for k in range(0,21):
                distz[k]= z_tot*k/20
                if k==0:
                    distz[k]=1e-3
                    N=13
                xMnA[k] = c_c0_S(tYr,xmn,distz[k],N,Rp)
                xMdA[k] = c_c0_S(tYr,xmd,distz[k],N,Rp)
                xMxA[k] = c_c0_S(tYr,xmx,distz[k],N,Rp)
            axsLeft[j].plot(distz,xMnA,c=clr1,marker="None")
            axsLeft[j].plot(distz,xMdA,c=clr2,marker="None")
            axsLeft[j].plot(distz,xMxA,c=clr3,marker="None")
        
            if(j==2): axsLeft[j].set_xlabel(xLabel)

    # Third subfigure to the right:  Two subplots (axes) in a single column.
    axsRight = subfigs[2].subplots(nrows=2, ncols=1)

    # do not show the axes
    axsRight[0].axis('off')
    axsRight[1].axis('off')

    # First 'axes': a text box listing parameters used.
    axsRight[0].annotate(txt, fontsize=8,
             xy=(0.03,0.97), xycoords='subfigure fraction',  verticalalignment='top',
             #bbox=dict(boxstyle="square,pad=0.5",linewidth=1,facecolor='white')
             )

    # Second 'axes': an empty plot defining the legend texts
    axsRight[1].plot([],[],c=clr1,marker=mrk1,linestyle='None',markersize=9,label="{:.2f}".format(Delta_x/2.)+" m")
    axsRight[1].plot([],[],c=clr2,marker=mrk2,linestyle='None',markersize=6,label="{:.2f}".format(xmd)+" m")
    axsRight[1].plot([],[],c=clr3,marker=mrk3,linestyle='None',markersize=9,label="{:.2f}".format(xmx)+" m")
    axsRight[1].legend(fontsize='small', loc=(0.,0.3), borderaxespad=0) #, bbox_to_anchor=(0.,2.)

    plotName = fileName+"-2.png"
    plt.savefig('plots/'+plotName)
    plt.show()
    
    plt.close()

def make_input_string():
    # Set up the input for Phreeqc
    input_string = """
#DATABASE """+database+"""

TITLE .
!------------------------------------------------------------------------
Muskus and Falta 2018. Semi-analytical method for matrix diffusion in heterogeneous
and fractured systems with parent-daughter reactions. Journal of Contaminant
Hydrology, 218, 94-109. https://doi.org/10.1016/j.jconhyd.2018.10.002.
.
The inflow source is maintained for a period of 50 years, before it is removed,
and the system is flushed with clean water for an additional 150 years.
The advective transport velocity is 100 m/y.
In the paper the model was set up as a 1-D system with 200 gridblocks using
a gridblock spacing of 1 m, and a time-step of 0.1 yr.    The dispersivity in
the Sudicky and Frind (1982) model was 0.5 m, but zero in their own model.
.
*** Here: """+str(n_cells)+""" cells """+"{:.4f}".format(Delta_t)+"""m long. Time step: """+"{:.2f}".format(Delta_t/(60*60*24))+""" days. ***
    numerical dispersivity: """+"{:.3f}".format(alpha_num)+""" m
    dispersivity used in Phreeqc: """+"{:.3f}".format(alpha_L)+""" m
.
!------------------------------------------------------------------------

SOLUTION_MASTER_SPECIES
Cl_	    Cl_O4-    0.0     99.45       35.453
N_      N_O3-     0.0     62.005      14.007

SOLUTION_SPECIES
Cl_O4- = Cl_O4-
    log_k	0
    -gamma  4.57  0.15
N_O3- = N_O3-
    log_k	0
    -gamma  3.0  0
CO3-2 + 10 H+ + 8 e- = CH4 + 3 H2O
    -log_k   -9999   # eliminate the formation of methane
    -delta_h  0

EXCHANGE_SPECIES
# For linear exchange, make ech exch. coeff. equal
    Na+ + X- = NaX
      -log_k   0.
      -gamma   4.08   0.082
    K+ + X- = KX
      -log_k   0.
      -gamma   3.5   0.015
#
USER_PRINT
-start
 10 Vol_soln = SOLN_VOL # (liters)
 20 density = RHO # kg_s/L
 30 Water_mass = TOT("water") # kg_w
 40 Mass_soln = Vol_soln * density
100 PRINT ".*********************************************************************************"
110 PRINT "density = ",TRIM(STR_F$(density,20,5))," kg/L"
120 PRINT "Vol_soln = ",TRIM(STR_F$(Vol_soln,20,5))," L"
130 PRINT "Mass_soln = ",TRIM(STR_F$(Mass_soln,20,5))," kg"
140 PRINT "H2O = ",TRIM(STR_F$(Water_mass,20,5))," kg"
200 yr = CALC_VALUE("Years")
210 IF(yr > 3.155e-7) THEN PRINT "Years =",TRIM(STR_F$(yr,20,5)) 
220 PRINT "STEP_NO =",STEP_NO, "   CELL_NO =",CELL_NO
250 IF(TOT("X") <= 0) THEN GOTO 999
260 PRINT "TOT(X)=",TOT("X"),", mol(X-)=",MOL("X-")
270 PRINT "  mol(KX)=",MOL("KX"),",  mol(NaX)=",MOL("NaX"),",  mol(RbX)=",MOL("RbX")
280 PRINT "  tot(K)= ",TOT("K"), ",  tot(Na)= ",TOT("Na"), ",  mol(Rb)= ",MOL("Rb")
999 PRINT "-*********************************************************************************"
-end

CALCULATE_VALUES
Hours
  -start
    10 SAVE TOTAL_TIME/3600
  -end
Days
  -start
    10 SAVE TOTAL_TIME/86400
  -end
Years
  -start
    10 SAVE TOTAL_TIME/31556952
  -end

SOLUTION 1  Initial mobile (fracture) solution: Pure water with a KCl tracer
    temp      25
    pH        7 charge
    pe        0 O2(g) -28
    redox     pe
    units     mol/kgw
    density   1 calculate
    O(0)      0
    K         1e-6
    Cl        1e-6
    -water    1 # kg
END

SOLUTION """+str(n_cells+2)+"""  Initial stagnant (rock matrix) solution: Pure water with a KClO4 tracer
    temp      25
    pH        7 charge
    pe        0 O2(g) -28
    redox     pe
    units     mol/kgw
    density   1 calculate
    O(0)      0
    K         1e-6
    Cl_       1e-6
    -water    1 # kg
EXCHANGE """+str(n_cells+2)+"""
   -equilibrate  """+str(n_cells+2)+"""
   X          1.e-6
END

TITLE .
Make copies of the exchanger and of the solutions

COPY solution 1 1-"""+str(n_cells+1)+"""
COPY solution """+str(n_cells+2)+" "+str(n_cells+2)+"-"+str((1+n_z)*n_cells+1)+"""
COPY exchange """+str(n_cells+2)+" "+str(n_cells+2)+"-"+str((1+n_z)*n_cells+1)+"""

END

SOLUTION 0  Inflowing water: Pure water with a NaNO3 tracer
    temp      25
    pH        7 charge
    pe        0 O2(g) -28
    redox     pe
    units     mol/kgw
    density   1 calculate
    O(0)      0
    Na        1e-6
    N_        1e-6
    -water    1 # kg
END

TITLE .
!------------------------------------------------------------------------
50 years:
The source is maintained for a period of 50 years, before it is removed,
and the system is flushed with clean water for an additional 150 years.
!------------------------------------------------------------------------

SELECTED_OUTPUT
   -file   """+tsv+"""
   -reset  false
   -distance  true
   -solution  true
   -pH                   true
   -ionic_strength       true
   -totals               Na N_ K Cl Cl_ X
   -molalities           NaX KX X-

USER_PUNCH
   -headings dist_z  Na_fr  NO3_fr  K_fr  Cl_fr  ClO4_fr  Years
-start
  10 IF (CELL_NO >= 0 AND CELL_NO <= """+str(n_cells+1)+""") THEN PUNCH 0"""
    lineNr = 10
    for i in range(n_z):
        lineNr = lineNr + 5
        input_string = input_string + '\n  '+str(lineNr)+' IF (CELL_NO >= '+str((i+1)*n_cells+2)+' AND CELL_NO <= '+str((i+2)*n_cells+1)+') THEN PUNCH '+" {:.4f}".format(z[i+1])
    input_string = input_string +'\n  '+str(lineNr+5)+' PUNCH TOT("Na")/1e-6,TOT("N_")/1e-6,TOT("K")/1e-6,TOT("Cl")/1e-6,TOT("Cl_")/1e-6,CALC_VALUE("Years")\n-end\n'
    input_string = input_string + """

#PRINT
#   -reset false
#   -echo_input true
#   -status false

INCLUDE$ """+pqi_mix+"""

TRANSPORT
   -cells  """+str(n_cells)+"""
   -lengths  """+str(n_cells)+"""*"""+str(Delta_x)+"""  # m
   -dispersivities  """+str(n_cells)+"""*"""+str(alpha_L)+""" # m
   -time_step """+"{:.0f}".format(Delta_t)+""" # s
   -shifts """+str(shifts)+"""
   -flow_direction  forward
   -boundary_conditions   flux  flux
   -correct_disp true
   -diffusion_coefficient """+str(D_0)+""" # m2/s
   -porosity """+str(n_cells)+"""*"""+"{:.5f}".format(phi_m)+"""  """+str(n_z*n_cells)+"""*"""+"{:.5f}".format(phi_im)+"""
   -stagnant   """+str(n_z)+"""
   -punch_frequency   """+str(punch_fr)+"""
   -print_frequency   """+str(print_fr)+"""

END

TITLE .
!------------------------------------------------------------------------
150 years:
The source is maintained for a period of 50 years, before it is removed,
and the system is flushed with clean water for an additional 150 years.
!------------------------------------------------------------------------

SOLUTION 0 Inflowing water: Initial mobile (fracture) solution. Pure water with a KCl tracer
    temp      25
    pH        7 charge
    pe        0 O2(g) -28
    redox     pe
    units     mol/kgw
    density   1 calculate
    O(0)      0
    K         1e-6
    Cl        1e-6
    -water    1 # kg
SAVE solution 0
END

TRANSPORT
   -shifts """+str(shifts2)+"""
   -punch_frequency   """+str(punch_fr2)+"""
   -print_frequency   """+str(print_fr2)+"""
   -initial_time  """+str(shifts*Delta_t)+""" # s
END
PRINT
   -reset true
END
"""
    return input_string

# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------

print("--- Model from Muskus and Falta (2018).")
# --- The model parameters:
b = 0.0001 # m (fracture aperture, decides porosity, etc)
spacing = 0.5 # m the distance between two fractures
phi_im = 0.1 # rock matrix porosity, Table 1 ("im" for immobile)
phi_m = b / (b + spacing) # (equivalent fracture porosity; "m" for mobile)
print('    phi_m = '+" {:.5f}".format(phi_m)+'   phi_im = '+str(phi_im))
D_0 = 1e-9 # m2/s diffusion coefficient—mobile zone, Table 1
tau = 0.1  # tortuosity in the rock matrix
D_im = D_0*tau # m2/s solute diffusion coefficient , Table 1
print('    D_im='+" {:.3e}".format(D_im))
print("    n_cells = "+str(n_cells))

L = 50. # m (length of rock model)
Delta_x = L/n_cells # m
print('    Delta_x = '+"{:.5f}".format(Delta_x)+' m')
v = 3.168874e-6 # = 100 m/yr (groundwater velocity)
print('    groundwater velocity, v = '+"{:.3e}".format(v)+' m/s,  ='+"{:.3e}".format(v*(60*60*24*365.25))+' m/yr')
Delta_t = Delta_x / v
alpha_L0 = 0.5 # m (dispersivity) Table 1
beta = phi_im/phi_m # plume attenuation factor in Lipson et al. (2007)
# estimated numerical spatial and temporal dispersivity
alpha_num = (Delta_x/2.) + (v*Delta_t/(2.*beta))
print('    given dispersivity: '+str(alpha_L0)+' m,  numerical dispersivity: '+"{:.4f}".format(alpha_num)+' m')
alpha_L = 0. # m (dispersivity)  # max(0,alpha_L0 - alpha_num)
print('    dispersivity used in Phreeqc: '+"{:.4f}".format(alpha_L)+' m')
Delta_t_D = (Delta_x**2) / (3*(D_im + alpha_L*v))
print('    Delta_t = '+"{:.0f}".format(Delta_t)+' s;  Delta_t_D = '+"{:.0f}".format(Delta_t_D)+' s')
if Delta_t_D > Delta_t:
    print('    (setting Delta_t_D = Delta_t in the finite-diffs. for the matrix)')
    Delta_t_D = Delta_t
endTime = 50 # years, the initial injection period
shifts = round(endTime*(60*60*24*365.25) / Delta_t)
endTime2 = 150 # years, the elution period
shifts2 = round(endTime2*(60*60*24*365.25) / Delta_t)

# Number of cells in the rock matrix
n_z = 5 if (modelType == 0 or modelType == 2) else 4

# the base-name of the Phreeqc input and output files
fileNameBase = 'M-F'
# modelType:
# - 0 (zero): 5 stagnant cells of equal size
# - 1 (one): 4 stagnant cells of equal size
# - 2 (two): 5 stagnant cells of increasing size
# - 3 (three): 4 stagnant cells of increasing size
if modelType == 0 or modelType == 1:  fileName = fileNameBase + '_'+str(n_z)+'L_'
elif modelType == 2 or modelType == 3:  fileName = fileNameBase + '_'+str(n_z)+'Lv_'
else:  fileName = fileNameBase + '_!!_'
fileName = fileName+str(n_cells)+"_"+str(fbc)
aTxt = "_a={:.3f}".format(alpha_L)
fileName = fileName + aTxt
# the input, output and selected output files:
pqi = fileName+'.pqi'  # the name of the Phreeqc input file
pqo = fileName+'.pqo'  # the name of the Phreeqc output file
tsv = fileName+'.tsv'  # the name of the Phreeqc selected output file
# The name of the Phreeqc input file containing the MIX entries
pqi_mix = fileName+'_MIX.pqi'

mixruns = ((D_im + alpha_L*v)*Delta_t) / (0.3333*Delta_x*Delta_x)
print('    shifts='+str(shifts)+'    end time = '+"{:.1f}".format((Delta_t*shifts)/(60*60*24*365.25))+' yr\n    mixruns='+"{:.3f}".format(mixruns))
# there is no need to punch more than about 100 time-step results
# 25 in the initial 50 years and 75 in the remaining 150 years
n_punch = 25
punch_fr = math.floor(shifts/n_punch) # punch frequency
if punch_fr < 1: punch_fr = 1
n_punch2 = 100 - n_punch
punch_fr2 = math.floor(shifts2/n_punch2) # punch frequency
if punch_fr2 < 1: punch_fr2 = 1
# for print, we nned 20 time-steps, 5 in the initial 50y period, and 15 for the rest
print_fr = math.floor(shifts/5)
if print_fr < 1: print_fr = 1
print_fr2 = math.floor(shifts2/15)
if print_fr2 < 1: print_fr2 = 1

# -------------------------------------------------------------------
# --- Modelling with RMD (rock matrix diffusion)
# The dimensions of the rock matrix model
z_tot = spacing / 2. # m (max diffusion depth)
# The dimensions of the rock matrix model
# modelType:
# - 0 (zero): 5 stagnant cells of equal size
# - 1 (one): 4 stagnant cells of equal size
# - 2 (two): 5 stagnant cells of increasing size
# - 3 (three): 4 stagnant cells of increasing size
n_z1 = n_z + 1
Delta_z = [0]*n_z1  # m (size (depth) of the rock matrix layer)
if modelType == 0 or modelType == 1: # stagnant cells of equal size (either 4 or 5 cells)
    for i in range(1, n_z1):
        Delta_z[i] = z_tot / n_z
elif modelType == 2: # Five (5) rock matrix cells of increasing size
    z_ = z_tot/15.
    for i in range(len(Delta_z)):
        Delta_z[i] = i * z_
elif modelType == 3: # Four (4) rock matrix cells of increasing size
    z_ = z_tot/10.
    for i in range(len(Delta_z)):
        Delta_z[i] = i * z_
# check that the Delta_z add up to z_tot
sum_z = 0
for i in range(len(Delta_z)):
    sum_z = sum_z + Delta_z[i]
if abs(sum_z - z_tot) > 1e-5:
        print("Error: sum of Delta_z...")
        raise SystemExit

# psi is the ratio of advective to diffusive solute flux in the fracture (Lipson et al 2007)
psi = (v/D_im)*(b/phi_im)*(Delta_z[1]/Delta_x)
print('    psi = '+str(psi))

# contact area between cells taking into account porosity, cell nbr zero is the fracture
A = [None]*(n_z1) # m2
for i in range(1, n_z1):
    A[i] = (Delta_x**2) * phi_im
# volumes of cells, cell nbr zero is the fracture
V = [None]*n_z1 # m3
V[0] = 0.5*b * Delta_x**2
for i in range(1, n_z1):
    V[i] = A[i] * Delta_z[i]

# midpoint distance between cells
# Note that the midpoint distance between the "fracture" and the first matrix cell 
# is decisive for the overall degree of matrix diffusion in the model.
h = [None]*n_z # m
h[0] = (0.5*b + Delta_z[1] ) /2.
for i in range(1, n_z):
    h[i] = (Delta_z[i]+Delta_z[i+1]) /2.

# Distace along matrix of cell center, used in USER_PUNCH only, for plotting
z = [None]*n_z1  # m
z[0] = b/4.
for i in range(1, n_z1):
    z[i] = z[i-1] + 0.5*(Delta_z[i-1]+Delta_z[i])

# boundary condition factor
if fbc != 0 and fbc !=1 and fbc !=2:
    print('wrong "fbc"')
    raise SystemExit()
f_bc = 2.*V[0]/(V[0]+V[1])  # eqn.(127) in Phreeqc user's guide
f_bc_ij = [1] * n_z
f_bc_jk = [1] * n_z1
if fbc == 0:
    f_bc_ij[1] = f_bc
else:
    f_bc_ij[1] = fbc
print("    b_fc = {:.3f}".format(f_bc)+"   using f_bc = {:.3f}".format(f_bc_ij[1]))

# Do we have alread a 'mix' file?
ok = pqi_tsv_ok(pqi_mix,tsv)
if ok:
    ok = pqi_tsv_ok(pqi,tsv)
if ok:
    print('    file: "'+pqi_mix+'" already exists...')
else:
    print("--- making input 'pqi' MIX file: '"+pqi_mix+"'")
    # --- Calculate the mixing factors:
    #   mixf_ij = D_im Delta_t A_ij f_bc / (h_ij V_j)
    #   mixf_jk = D_im Delta_t A_jk f_bc / (h_jk V_j)
    # if 'n' the cell, with 1=fracture, 2=first matrix cell, 3=second matrix cell, etc
    # so between fracture and 1st matrix cell 'n' = 2
    # A_ij = A[n-1], A_jk = A[n]
    # h_ij = h[n-2], h_jk = h[n-1]
    # f_bc_ij = f_bc_ij[n-2]
    # f_bc_jk = f_bc_ij[n-1]
    D_Delt = D_im * Delta_t_D
    # calculate 'mixes' and write them into file 'pqi_mix'
    # which is included into the main input file at run time.
    with open(pqi_mix,'w') as f:
        #print("cell  n   i   j   k  mixfij   mixfjj    mixfjk")  # for debugging
        #for cell_i in range(1,2): # for debugging
        for cell_i in range(1, n_cells+1):
            i = -1
            for n in range(1, (n_z+2)):
                k = -1
                if n < n_z1: k = n * n_cells + 1 + cell_i
                j = cell_i
                if n > 1:
                    i = (n-2) * n_cells + 1 + cell_i
                    j = (n-1) * n_cells + 1 + cell_i
                if i <= n_cells+1: i = i-1
                mixf_ij = 0.
                if i > 0:
                    mixf_ij = D_Delt * A[n-1]*f_bc_ij[n-2] / (h[n-2]*V[n-1])
                mixf_jk = 0.
                if k > 0:
                    mixf_jk = D_Delt * A[n]*f_bc_jk[n-1] / (h[n-1]*V[n-1])
                mixf_jj = 1. - mixf_ij - mixf_jk
                #print("   "+str(cell_i)+"  "+str(n)+'{:4}'.format(i)+'{:4}'.format(j)+'{:4}'.format(k)+" "+" {:.6f}".format(mixf_ij)+" {:.6f}".format(mixf_jj)+" {:.6f}".format(mixf_jk))  # for debugging
                if cell_i == 1:  # test that mixing factors are between zero and one
                    ok = True
                    if (mixf_ij < 0.) | (mixf_ij > 1.): ok=False; print('cell=1, n='+str(n)+' i='+str(i)+' j='+str(j)+' k='+str(k)+' mixf_ij='+str(mixf_ij))
                    if (mixf_jk < 0.) | (mixf_jk > 1.): ok=False; print('cell=1, n='+str(n)+' i='+str(i)+' j='+str(j)+' k='+str(k)+' mixf_ij='+str(mixf_jk))
                    if (mixf_jj < 0.) | (mixf_jj > 1.): ok=False; print('cell=1, n='+str(n)+' i='+str(i)+' j='+str(j)+' k='+str(k)+' mixf_ij='+str(mixf_jj))
                    if not ok: raise SystemExit('Increase number of cells...')
                mix_str = "MIX "+'{:4}'.format(j)
                if i > 0: mix_str = mix_str +"; "+'{:5}'.format(i)+" {:.6f}".format(mixf_ij)
                mix_str = mix_str +"; "+'{:5}'.format(j)+" {:.6f}".format(mixf_jj)
                if k > 0: mix_str = mix_str +"; "+'{:5}'.format(k)+" {:.6f}".format(mixf_jk)
                #print(mix_str)   # for debugging
                f.writelines(mix_str+"\n")
    print('    written: "'+pqi_mix+'"')

ok = pqi_tsv_ok(pqi,tsv)
if ok:
    print('    file: "'+tsv+'" already exists\n    reading previous results...')
    res=pd.read_csv(tsv,sep='\t',encoding='ANSI',decimal='.')
    res.rename(columns=lambda x: x.strip(), inplace=True)
    res.drop(res.columns[len(res.columns)-1],axis="columns", inplace=True)
else:
    print('    making input pqi-file: "'+pqi+'"')
    input_string = make_input_string()
    # Save the input file in case you want to run PhreeqC outside Python
    with open(pqi,'w') as f:
        f.writelines(input_string)
    print("    'pqi' file written.")
    
    answ = input(">>> Press 'q' to quit or <Enter> to continue and run Phreeqc...")
    if (answ == 'q') | (answ == 'Q'): raise SystemExit('quit...')
    
    # Make the Phreeqc calculations
    res = runPhreeqc(input_string)

# analyse the results...
print('--- Results:')
time_1 = 50; time_2 = 100; time_3 = 150; time_4 = 200
print("--- finding results close to "+str(time_1)+", "+str(time_2)+", "+str(time_3)+" and "+str(time_4)+" years")
t_1 = -1; t_2 = -1; t_3 = -1; t_4 = -1
dist_1 = 1e20; dist_2 = 1e20; dist_3 = 1e20; dist_4 = 1e20
times = res['Years']
for t in times:
    if(abs(t-time_1) < dist_1):
        t_1 = t
        dist_1 = abs(t-time_1)
    if(abs(t-time_2) < dist_2):
        t_2 = t
        dist_2 = abs(t-time_2)
    if(abs(t-time_3) < dist_3):
        t_3 = t
        dist_3 = abs(t-time_3)
    if(abs(t-time_4) < dist_4):
        t_4 = t
        dist_4 = abs(t-time_4)
print('    '+str(time_1)+' years = '+str(t_1)+'\n    '+str(time_2)+' yr = '+str(t_2)+'\n    '+str(time_3)+' yr = '+str(t_3)+'\n    '+str(time_4)+' yr = '+str(t_4))
if (t_1 <0) | (t_2 <0) | (t_3 < 0) | (t_3 < 0):
    print('Error, one of these values is negative...')
    raise SystemExit

print("--- x-coordinate: finding 0, 7.5 and 15 m")
# x-coordinate: get three distances, the minimum, 15 m and in between
xmn = 1e30
xmx = -1e20
dx_5L = res.loc[(res['dist_x']>0) & (res['Years']==0) & (res['soln']>0) & (res['soln']<(n_cells+1)),'dist_x']
for x in dx_5L:
    if(x < xmn):
        xmn = x
    if((x > xmx) & (x < 15.001)):
        xmx = x
mid = (xmx - xmn)/2.
dist = 1e20
xmd = -1
for x in dx_5L:
    if(abs(x-mid) < dist):
        xmd = x
        dist = abs(x-mid)
print(' X-distance:  min='+str(xmn)+',  mid='+str(xmd)+',  max='+str(xmx))
if (xmn >1e25) | (xmd <0) | (xmx < 0):
    print('Error, one of these three values is wrong...')
    raise SystemExit

print("--- making plots.")
# 'include' is the frequency of markers on the plots.
# there should be around n_cells steps in the output file.
# but it is enough to plot around 10 markers
incl_x = round(n_cells/10)
plot(res,incl_x)
