# -*- coding: utf-8 -*-
"""
Script to reproduce the explicit single fracture model described in:
   Tang, D.H., E.O. Frind, and E.A. Sudicky. 1981. Contaminant transport in
   fractured porous media: Analytical solution for a single fracture,
   Water Resources Research 17, 555-564.
   https://doi.org/10.1029/WR017i003p00555
See Figs 9 and 10 in the paper.
This example is also described in:
  - Weatherill, Graf, Simmons, Cook, Therrien and Reynolds (2008)
    Discretizing the Fracture-Matrix Interface to Simulate Solute Transport,
    Groundwater, 46, 606-615. https://doi.org/10.1111/j.1745-6584.2007.00430.x
  - Watanabe and Kolditz (2015), Numerical stability analysis of two-dimensional
    solute transport along a discrete fracture in a porous rock matrix,
    Water Resour. Res., 51, 5855-5868, https://doi.org/10.1002/2015WR017164
  
This script will create two PHREEQC input files
   - without RMD (rock matrix diffusion)
   - with RMD
using the parameters given in the script itself.  The discretization used for
the rock matrix model may be selected through variable 'modelType'.
The script will then run the selected model using
[IPHREEQC](https://www.usgs.gov/software/phreeqc-version-3):
   Charlton and Parkhurst (2011) http://dx.doi.org/10.1016/j.cageo.2011.02.005
Except that if the input tile (*.pqi) and the corresponding 'selected-output'
file (*.tsv) both exist, and the selected-output file is newer than the input file,
then the calculations are NOT performed, and the esisting tsv-file is instead
used to create figures with the results.  To force the calculation, just remove
either the pqi- or the tsv-files.
The reason for not repeating calculations is that a PHREEQC run
may take quite a long time.

Finally the results from the two simulations are plotted to create two figures,
where the simulation results are compared with those from an analytical equation.
See Figs. 9 and 10 of Tang et al. (1981), Fig. 4 in Weatherill et al. (2008),
Figs. 6-8 in
   Grisak and Pickens (1980). Solute transport through fractured media:
   1. The effect of matrix diffusion. Water Resources Research 16, 719–730.
   https://doi.org/10.1029/WR016i004p00719


The figures are saved in subfolder 'plots'.
An analytical equation is also plotted for comparison.

The script will run under Windows.

@author: Ignasi Puigdomenech
"""
import os, sys
import shutil
import math
from win32com.client import Dispatch
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import pandas as pd
from time import time, strftime, localtime
sys.path.insert(0, '../Analytical_eqn/CraFlush')
from craflush import CraFlush
# The Phreeqc database
database = '../phreeqc_3.7.3.dat'

# The model is 0.8 m long. The value of 'n_cells' controls:
#    - Delta_x, the cell x-length
#    - Delta_t, the time step
#    - the number of shifts to reach the desired number of days
# larger values give lower numerical dispersion,
# smaller values give shorter execution times,
# should be between 20 and perhaps 1000
n_cells = 100

# The model type selects the type of rock-matrix model:
# - 0 (zero): 5 stagnant cells of equal size
# - 1 (one): 5 stagnant cells of increasing size
# - 2 (two): 10 stagnant cells of increasing size
# Increasing sizes require a larger value for 'n_cells'
# Five cells (instead of 10) require shorter run time, but might be less accurate.
modelType = 1

# Analytical model to be used in the plots
modelA = 1 # 0 = Neretnieks (1980);  1 = Sudicky et al.

# Boundary cell factor, see eqn.(127) in Phreeqc user's guide (Parkhurst and Appelo 1999)
fbc = 1 # it can be =1, =2 or =0 which means calculated according to eqn.(127)

# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------

def runPhreeqc(input_string,outFileName):
    """Load database via linked library object and run the input string.
    The results (selected_output) are in a "tuple" of tuples, converted
    into a Pandas DataFrame (https://pandas.pydata.org/docs/getting_started/index.html).
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
    with open(outFileName,'wb') as wfd:
        for f in ['header_pqo',iphreeqc.OutputFileName]:
            with open(f,'rb') as fd:
                shutil.copyfileobj(fd, wfd)
    os.remove('header_pqo')
    print("    (removing file",iphreeqc.OutputFileName,")")
    os.remove(iphreeqc.OutputFileName)

    # Get the selected_output results
    # The results (selected_output) are in a "tuple" of tuples.
    selOutput = iphreeqc.GetSelectedOutputArray()
    # Create a Pandas DataFrame
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

def c_c0_S(tdays,theta,X,Z,N):
    # Implements the analytical equation for matrix diffusion from Sudicky et al.
    # N decides the number of terms in the EPAL summation
    # theta = Matrix porosity
    # units: x and z in 'm', tdays in 'days'
    t = tdays * 86400 # convert days to 's'
    # Data for the PhreeqC 'model':
    # b = fracture aperture
    C0 = 1.     # Source concentration at infinite time
    Cim = 0.    # Initial concentration in matrix and fractures
    alpha = alpha_A # alpha_num + alpha_L  # Fracture dispersivity (m)
    sep = spacing   # Fracture spacing (m)
    # theta
    D0 = D_0    # Diffusion coefficient in water (m2/s)
    if D0 <= 0: D0 = 1.e-9
    tau = D_im/D0    # Matrix tortuosity
    R = 1.      # Fracture retardation factor
    Rp = 1.     # Matrix retardation factor
    THalf = 0.  # non-decaying solute
    # Numerical inversion parameters used in CRUMP:
    AlphR = 0.   # alpha
    relErr = 1.e-6  # relatve error
    #Variable source concentration
    Ns = 1
    Ts =  [0.]
    Css = [1.]
    # ---- the equations by Sudicky
    if(X <= 0): X=1e-6
    if(Z <= 0): Z=b/2.
    cc0, Cs = CraFlush(C0,Cim,v,alpha,b,sep,theta,tau,D0,R,Rp,THalf,
                             X,Z,t, Ns,Ts,Css,AlphR,relErr,N)
    if(cc0 > 1.): cc0 =1.
    if(cc0 < 0.): cc0 =0.
    return cc0

def c_c0_N(tdays,x,z):
    # Implements the analytical equation for matrix diffusion from Neretnieks (1980),
    # which excludes the effect of dispersion, and can not be used without
    # rock matrix diffusion.
    # Units: x and z in 'm', tdays in 'days'
    t = tdays * 86400 # convert days to 's'

    # Data for the PhreeqC 'model':
    v_f = Delta_x / Delta_t # groundwater velocity in the Phreeqc model
    D_a_ = D_im # the diffusion coefficient in the rock matrix
    # Neretnieks states that D_a = D_e / (Kd rho_p)
    # and that for non-sorbing species (Kd rho_p) = matrix porosity
    # D_e_ = D_a_ * phi_im
    # ---- the equation by Neretnieks (1980)
    if x <= 0: x=1e-6
    t_w = x / v_f # water travel time: the time for groundwater to reach position 'x'
    G_x_z = (((D_a_ * phi_im) + (0.5*v_f*b*z)/x)/(b*math.sqrt(D_a_)))*t_w
    if (t-t_w) <= 0:
        cc0 = 0
    else:
       #print('tdays='+str(tdays)+' t='+str(t)+'  t_w='+str(t_w))
        cc0 = math.erfc(G_x_z/(math.sqrt(t-t_w)))
    return cc0

def annotTxt(modelA,RMD):
    # modelA: 0 = Neretnieks (1980);  1 = Sudicky et al.
    modelTxt = 'No rock matrix diffusion.$_{}$'
    if RMD:
        if modelType == 0:
            modelTxt = str(n_z)+' equal stagnant layers$_{}$'
        else:
            modelTxt = str(n_z)+' stagnant layers of\n    increasing size$_{}$'
    annotateTxt = (
        modelTxt+"\n"+
        "$n\mathsf{_{cells}}$ = "+str(n_cells)+"\n"+
        "$\\alpha_L$ = "+"{:.3f}".format(alpha_L)+" m\n"+
        "$D_{im}$ = "+"{:.2e}".format(D_im)+" m$^2$/s\n")  # used in PHREEQC
    if RMD: annotateTxt = annotateTxt + "$f_{bc}$ = "+"{:.2f}".format(f_bc_ij[1])+"\n"
    annotateTxt = annotateTxt +("Markers: PHREEQC$_{ }$\n"+
        "Lines: analytical model")
    if modelA == 0:
        annotateTxt = annotateTxt + "\n  Neretnieks (1980)"
    else:
        annotateTxt = annotateTxt + "\n  Sudicky et al."
    annotateTxt = annotateTxt + "$^{ }$\n"+"using: $\\alpha_{tot}$ = "+"{:.3f}".format(alpha_A)+" m\n"
    if modelA == 0:
        annotateTxt = annotateTxt + "$D_a$ = "+"{:.2e}".format(D_im)+" m$^2$/s"
    else:
        annotateTxt = annotateTxt + "$D'$ = "+"{:.2e}".format(D_im)+" m$^2$/s"
    return annotateTxt

def plot(rsl,RMD,incl_x,incl_t):
    # Plot results into a single figure containing two plots.
    # The figure is saved in folder 'plots'
    # Not all points are plotted, the variables incl_x and incl_t are used to select how many points are plotted.
    if incl_x < 1: incl_x = 1
    print("--- making plots, RMD = "+str(RMD))

    # Create a figure containing two rows.
    fig, ax=plt.subplots(2,layout='constrained',figsize=(5.12,3.84),dpi=100)

    D_1 = rsl[(rsl["Days"]>=t_1*0.999) & (rsl["Days"]<=t_1*1.001) & (rsl["dist_x"]>=0) & (rsl["dist_x"]<L) & (rsl["soln"]<=n_cells)].iloc[::incl_x,:]
    D_2 = rsl[(rsl["Days"]>=t_2*0.999) & (rsl["Days"]<=t_2*1.001) & (rsl["dist_x"]>=0) & (rsl["dist_x"]<L) & (rsl["soln"]<=n_cells)].iloc[::incl_x,:]
    D_3 = rsl[(rsl["Days"]>=t_3*0.999) & (rsl["Days"]<=t_3*1.001) & (rsl["dist_x"]>=0) & (rsl["dist_x"]<L) & (rsl["soln"]<=n_cells)].iloc[::incl_x,:]
    # colors from
    # https://davidmathlogic.com/colorblind/#%23000000-%23E69F00-%2356B4E9-%23009E73-%23F0E442-%230072B2-%23D55E00-%23CC79A7
    clr1 = "#56B4E9"  # light blue
    clr2 = "#117733"  # green
    clr3 = "#332288"  # dark purple
    # markers: see the Notes section in
    # https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html
    mrk1 = "v"  # triangle down
    mrk2 = "D"  # diamond
    mrk3 = "^" # triangle up
    xLabel = 'Distance x along fracture (m)'
    yLabel = 'C/C$_0$'
    # ax[0].axis([0, 3, 0, 1.])
    ax[0].plot(D_1["dist_x"],D_1["Ti_fr"],c=clr1,marker=mrk1,markersize=9,linestyle='None',label="{:.1f}".format(t_1)+" days")
    ax[0].plot(D_2["dist_x"],D_2["Ti_fr"],c=clr2,marker=mrk2,markersize=6,linestyle='None',label="{:.1f}".format(t_2)+" d")
    ax[0].plot(D_3["dist_x"],D_3["Ti_fr"],c=clr3,marker=mrk3,markersize=9,linestyle='None',label="{:.1f}".format(t_3)+" d")
    # calculate and plot the analytical (theoretical) curves
    theta = phi_im if RMD else 0.  # if theta=0 rock matrix diffusion is not included
    distx = [0]*101
    y1 = [1]*101
    y2 = [1]*101
    y3 = [1]*101
    for i in range(1,101):
        distx[i] = L*i/100
        if modelA == 0:
            if RMD:
                y1[i] = c_c0_N(t_1,distx[i],0)
                y2[i] = c_c0_N(t_2,distx[i],0)
                y3[i] = c_c0_N(t_3,distx[i],0)
        else:
            N = 11
            y1[i] = c_c0_S(t_1,theta,distx[i],0,N)
            N = 12
            y2[i] = c_c0_S(t_2,theta,distx[i],0,N)
            y3[i] = c_c0_S(t_3,theta,distx[i],0,N)
    if modelA != 0 or RMD: 
        ax[0].plot(distx,y1,c=clr1,marker="None")
        ax[0].plot(distx,y2,c=clr2,marker="None")
        ax[0].plot(distx,y3,c=clr3,marker="None")
    ax[0].set_ylabel(yLabel)
    ax[0].set_xlabel(xLabel)
    
    ax[0].legend(fontsize='small', bbox_to_anchor=(1.07,0.), loc='lower left', borderaxespad=0)

    #clr1 = "#AA4499" # mauve
    clr2 = "#CC6677"  # orange
    #clr3 = "#FFB000"  # yellow
    #mrk1 = "v"  # triangle down
    mrk2 = "D"  # diamond
    #mrk3 = "^" # triangle up
    xLabel = 'Time (days)'
    if incl_t < 1: incl_t = 1
    d=rsl[(rsl["dist_x"]>=x076*0.999) & (rsl["dist_x"]<=x076*1.001) & (rsl["dist_z"]==0)].iloc[::incl_t,:]
    #ax[1].axis([0, 0.25, 0, 1.])
    ax[1].set_ylim([-0.1, 1.1])
    ax[1].plot(d['Days'],d['Ti_fr'],c=clr2,marker=mrk2,linestyle='None',markersize=4,label="{:.2f}".format(x076)+" m")
    ax[1].set_ylabel(yLabel)
    ax[1].set_xlabel(xLabel)
    ax[1].set_title("$x=$"+"{:.2f}".format(x076)+" m", x=0.7, y=0.7)
    # plot the analytical (theoretical) curve
    timex = [0]*101
    y1 = [0]*101
    X = 0.76 # distance along fracture (m)
    N = 14
    for i in range(1,101):
        timex[i] = 4*i/100
        if modelA == 0:
            if RMD: y1[i] = c_c0_N(timex[i],X,0)
        else:
            y1[i] = c_c0_S(timex[i],theta,X,0,N)
    if modelA != 0 or RMD: ax[1].plot(timex,y1,c=clr2,marker="")

    txt = annotTxt(modelA,RMD)
    ax[1].annotate(txt, xy=(1.07, 1.), xycoords='axes fraction',
            fontsize=8, verticalalignment='top')

    plotName = fileName+"_Fract"+".png"
    if not RMD: plotName = fileName+"_noRMD"+".png"
    plt.savefig('plots/'+plotName)
    plt.show()
    
    plt.close()

def plotMtrx(rsl):
    # Plot results into a single figure containing one plots.
    # The figure is saved in folder 'plots'
    print("--- making plots along the matrix")

    # Create a figure containing one row and one column.
    fig, ax=plt.subplots(2,layout='constrained',figsize=(5.12,3.84),dpi=100)
    # colors from
    # https://davidmathlogic.com/colorblind/#%23000000-%23E69F00-%2356B4E9-%23009E73-%23F0E442-%230072B2-%23D55E00-%23CC79A7
    clr1 = "#DC267F" # mauve
    clr2 = "#FE6100"  # orange
    clr3 = "#FFB000"  # yellow
    # markers: see the Notes section in
    # https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.plot.html
    mrk1 = "v"  # triangle down
    mrk2 = "D"  # diamond
    mrk3 = "^" # triangle up
    xLabel = 'Distance z in matrix (m)'
    yLabel = 'C/C$_0$'

    days = [t_2,t_3]
    for i in range(0,2):        
        d = days[i]
        dMn=res[(res["Days"]==d) & (res["dist_x"]==Delta_x/2.)]
        dMd=res[(res["Days"]==d) & (res["dist_x"]==xmd)]
        dMx=res[(res["Days"]==d) & (res["dist_x"]==xmx)]
        # ax[0].axis([0, 3, 0, 1.])
        ax[i].yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
        ax[i].plot(dMn["dist_z"],dMn["Ti_fr"],c=clr1,marker=mrk1,markersize=9,linestyle='None',label="x= {:.2f}".format(Delta_x/2.)+" m")
        ax[i].plot(dMd["dist_z"],dMd["Ti_fr"],c=clr2,marker=mrk2,markersize=6,linestyle='None',label="x= {:.2f}".format(xmd)+" m")
        ax[i].plot(dMx["dist_z"],dMx["Ti_fr"],c=clr3,marker=mrk3,markersize=9,linestyle='None',label="x= {:.2f}".format(xmx)+" m")
        # calculate and plot the analytical (theoretical) curves
        theta = phi_im  # if theta=0 rock matrix diffusion is not included
        distz = [0]*21
        yXmn = [1]*21;  yXmd = [1]*21;  yXmx = [1]*21
        N = 11
        for j in range(0,21):
            distz[j]= z_tot_max*j/20
            if modelA == 0:
                yXmn[j] = c_c0_N(d,xmn,distz[j])
                yXmd[j] = c_c0_N(d,xmd,distz[j])
                yXmx[j] = c_c0_N(d,xmx,distz[j])
            else:
                yXmn[j] = c_c0_S(d,theta,xmn,distz[j],N)
                yXmd[j] = c_c0_S(d,theta,xmd,distz[j],N)
                yXmx[j] = c_c0_S(d,theta,xmx,distz[j],N)
        ax[i].plot(distz,yXmn,c=clr1,marker="None")
        ax[i].plot(distz,yXmd,c=clr2,marker="None")
        ax[i].plot(distz,yXmx,c=clr3,marker="None")
        if i == 1: ax[i].set_xlabel(xLabel)
        ax[i].set_ylabel(yLabel)
        ax[i].set_title("$t=$"+"{:.2f}".format(d)+" days", x=0.7, y=0.7)

    ax[1].legend(fontsize='small', bbox_to_anchor=(1.07,0.), loc='lower left', borderaxespad=0)

    txt = annotTxt(modelA,True)
    ax[0].annotate(txt, xy=(1.07, 1.), xycoords='axes fraction',
            fontsize=8, verticalalignment='top')

    plotName = fileName+"_Matrix"+".png"
    plt.savefig('plots/'+plotName)
    plt.show()
    
    plt.close()

def make_input_string(pqi_mix_name,tsv_name,RMD):
    # Set up the input for Phreeqc
    # 'pqi_mix_name' is the name of the input file containing the mixing fractions to calculate diffusion in the rock matrix
    # 'tsv_name' is the name of the selected output file
    # 'RMD' is True or False to include or not rock matrix diffusion
    txt = """
#DATABASE D:/OneDrive/Documents/Geokemi/Modelling/PhreeqC/Databases/phreeqc_3.7.3.dat

TITLE .
!------------------------------------------------------------------------
Script to reproduce the explicit single fracture model described in:
   Tang, D.H., E.O. Frind, and E.A. Sudicky. 1981. Contaminant transport
   in fractured porous media: Analytical solution for a single fracture,
   Water Resources Research 17, no. 3: 555-564.
   https://doi.org/10.1029/WR017i003p00555
See Figs 9 and 10 in the reference.  This example is also described in:
  - Weatherill, Graf, Simmons, Cook, Therrien and Reynolds (2008)
    Discretizing the Fracture-Matrix Interface to Simulate Solute Transport,
    Groundwater, 46: 606-615. https://doi.org/10.1111/j.1745-6584.2007.00430.x
  - Watanabe and Kolditz (2015), Numerical stability analysis of two-dimensional
    solute transport along a discrete fracture in a porous rock matrix,
    Water Resour. Res., 51, 5855-5868, https://doi.org/10.1002/2015WR017164
.
A single row of mobile cells (fracture) and the bedrock matrix adjacent to
the fracture is simulated as an array of five layers of stagnant cells.
.
The bedrock model is """+str(L)+""" m long, discretized into """+str(n_cells)+""" cells, Delta-x = """+str(Delta_x)+""" m.
As seen in Fig.2 of Weatherill et al (2008) and in Fig.1 of Watanabe and Kolditz (2015)
the matrix depth is 0.05 m.
.
The groundwater velocity is 0.75 m/day.
.
Three tracers: input, fracture and rock matrix
'Ti_', 'Tf_', and 'Tm_'.
!------------------------------------------------------------------------

SOLUTION_MASTER_SPECIES
Ti_                 Ti_         0.0     Ti_                 1.
Tf_                 Tf_         0.0     Tf_                 1.
Tm_                 Tm_         0.0     Tm_                 1.

SOLUTION_SPECIES
Ti_ = Ti_
        log_k   0.0
Tf_ = Tf_
        log_k   0.0
Tm_ = Tm_
        log_k   0.0
CO3-2 + 10 H+ + 8 e- = CH4 + 3 H2O
	-log_k   -9999   # eliminate the formation of methane
	-delta_h  0

USER_PRINT
-start
 10 Vol_soln = SOLN_VOL # (liters)
 20 density = RHO # kg_s/L
 30 Water_mass = TOT("water") # kg_w
 40 Vol_soln = SOLN_VOL # L
 50 Mass_soln = Vol_soln * density
100 PRINT ".*********************************************************************************"
110 PRINT "density = ",TRIM(STR_F$(density,20,5))," kg/L"
120 PRINT "Vol_soln = ",TRIM(STR_F$(Vol_soln,20,5))," L"
130 PRINT "Mass_soln = ",TRIM(STR_F$(Mass_soln,20,5))," kg"
140 PRINT "H2O = ",TRIM(STR_F$(Water_mass,20,5))," kg"
200 yr = CALC_VALUE("Years")
210 IF(yr > 3.155e-7) THEN PRINT "Years =",TRIM(STR_F$(yr,20,5)) 
220 PRINT "STEP_NO =",STEP_NO, "   CELL_NO =",CELL_NO
999 PRINT "-*********************************************************************************"
-end

CALCULATE_VALUES
Hours
  -start
    10 SAVE SIM_TIME/3600
  -end
Days
  -start
    10 SAVE SIM_TIME/86400
  -end
Years
  -start
    10 SAVE SIM_TIME/31556952
  -end

SOLUTION 0  Inflowing water: Pure water with a tracer
    temp      25
    pH        7 charge
    pe        0 O2(g) -28
    redox     pe
    units     mol/kgw
    density   1 calculate
    O(0)      0
    Ti_       1e-6
    -water    1 # kg
END

SOLUTION 1-"""+str(n_cells+1)+""" Fracture water: Pure water with a tracer
    temp      25
    pH        7 charge
    pe        0 O2(g) -28
    redox     pe
    units     mol/kgw
    density   1 calculate
    O(0)      0
    Tf_       1e-6
    -water    1 # kg
END
"""
    if RMD:
        txt = txt + """
SOLUTION """+str(n_cells+2)+"""-"""+str((1+n_z)*n_cells+1)+""" Matrix water: Pure water with a tracer
    temp      25
    pH        7 charge
    pe        0 O2(g) -28
    redox     pe
    units     mol/kgw
    density   1 calculate
    O(0)      0
    Tm_       1e-6
    -water    1 # kg
END"""
    txt = txt + """

#PRINT
#   -reset false
#   -echo_input true
#   -status false

INCLUDE$ """+pqi_mix_name+"""

TRANSPORT
   -cells  """+str(n_cells)+"""
   -lengths  """+str(n_cells)+"""*"""+str(Delta_x)+"""  # m, total length = """+str(L)+"""
   -dispersivities  """+str(n_cells)+"""*"""+str(alpha_L)+""" # m
   -time_step """+"{:.2f}".format(Delta_t)+""" # s
   -shifts """+str(shifts)+"""
   -flow_direction  forward
   -boundary_conditions   flux  flux
   -correct_disp true
   -diffusion_coefficient """+str(D_0)+""" # m2/s
   -porosity """+str(n_cells)+"""*"""+"{:.5f}".format(phi_m)
    if RMD:
        txt = txt + "  "+str(n_z*n_cells)+"*"+"{:.5f}".format(phi_im)
    txt = txt + """
   -stagnant   """+str(n_z)+"""
   -punch_frequency   """+str(punch_fr)+"""
   -print_frequency   """+str(print_fr)+"""

SELECTED_OUTPUT
   -file   """+tsv_name+"""
   -reset  false
   -distance  true
   -solution

USER_PUNCH
   -headings dist_z Ti_fr Tf_fr Tm_fr Days
-start
  10 IF (CELL_NO >= 0 AND CELL_NO <= """+str(n_cells+1)+""") THEN PUNCH 0"""
    lineNr = 10
    if RMD:
        for i in range(n_z):
            lineNr = lineNr + 5
            txt = txt + '\n  '+str(lineNr)+' IF (CELL_NO >= '+str((i+1)*n_cells+2)+' AND CELL_NO <= '+str((i+2)*n_cells+1)+') THEN PUNCH '+" {:.4f}".format(z[i+1])
    txt = txt +'\n  '+str(lineNr+5)+' PUNCH TOT("Ti_")/1e-6,TOT("Tf_")/1e-6,TOT("Tm_")/1e-6,CALC_VALUE("Days")\n-end\n'
    txt = txt + """
END
PRINT
   -reset true
END
"""
    return txt

# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------

print("--- Model from Tang et al (1981).")
b = 0.00012 # m (fracture aperture, decides porosity, etc)
spacing = 0.1 # m the distance between two fractures from Fig.2 in Weatherill et al (2008) and Watanabe and Kolditz (2015)
phi_im = 0.35 # (rock matrix porosity; "im" for immobile)
phi_m = b / (b + spacing) # (equivalent fracture porosity; "m" for mobile)
print('    phi_m = '+" {:.5f}".format(phi_m)+'   phi_im = '+str(phi_im))
D_0 = 1e-9 # m2/s diffusion coefficient — mobile zone
# We wish to reproduce the lines in Figs.9 and 10 of Tang et al (1981) for D'=1e-8 cm2/s (1e-12 m2/s)
# In this case R' = 1 (non-sorbing solute).
D_im = 1e-12 # m2/s, diffusion coefficient in the rock matrix ("im" for immobile)
print('    D_im='+" {:.3e}".format(D_im)+' m2/s')
print("    n_cells = "+str(n_cells))

L = 0.8 # m (length of rock model)
Delta_x = L/n_cells # m
print('    Delta_x = '+"{:.5f}".format(Delta_x)+' m')
v = 0.75 /(24*60*60) # 0.75 m/day = 8.68056e-6 m/s (groundwater velocity)
print('    groundwater velocity, v = '+"{:.5e}".format(v)+' m/s   ('+"{:.4f}".format(v*(60*60*24))+' m/day)')
Delta_t = Delta_x / v # = 0.01d in the original paper, with n_cells = 200
alpha_L0 = 0.76 # m dispersivity given in Tang et al (1981).
alpha_L = 0.01 # m dispersivity used in PhreeqC

Delta_t_D = (Delta_x**2) / (3*(D_im + alpha_L*v))
print('    Delta_t = '+"{:.0f}".format(Delta_t)+' s;  Delta_t_D = '+"{:.0f}".format(Delta_t_D)+' s')
if Delta_t_D > Delta_t:
    print('    (setting Delta_t_D = Delta_t in the finite-diffs. for the matrix)')
    Delta_t_D = Delta_t

#beta = phi_im/phi_m # plume attenuation factor in Lipson et al. (2007)  ---------------
# estimated numerical spatial and temporal dispersivity
alpha_num = (Delta_x/2.) + (v*Delta_t/(2.)) 
print('    given dispersivity: '+str(alpha_L0)+' m,  numerical dispersivity: '+"{:.4f}".format(alpha_num)+' m')
print('    dispersivity used in Phreeqc: '+"{:.4f}".format(alpha_L)+' m')
alpha_A = alpha_L + alpha_num  # m dispersivity used in the analytical model

endTime = 4 # days
shifts = 1 + math.floor(endTime * (60*60*24) / Delta_t)

# The model type selects the type of rock-matrix model:
# - 0 (zero): 5 stagnant cells of equal size
# - 1 (one): 5 stagnant cells of increasing size
# - 2 (two): 10 stagnant cells of increasing size
n_z = 5 if (modelType == 0 or modelType ==1) else 10

# the base-name of the Phreeqc input and output files
fileNameBase = '1981Tang'
if modelType == 0:  fileName = fileNameBase + '-'+str(n_z)+'L_'
else:  fileName = fileNameBase + '-'+str(n_z)+'Lv_'
fileName = fileName+str(n_cells)+"_"+str(fbc)
aTxt = "_a={:.3f}".format(alpha_L)
fileName = fileName + aTxt
# for the simulations without RMD: the base-name of the Phreeqc input and output files
fileName_noRMD = fileNameBase + '-noRMD_'
fileName_noRMD = fileName_noRMD + str(n_cells)
fileName_noRMD = fileName_noRMD + aTxt

# the input, output and selected output files:
pqi = fileName+'.pqi'  # the name of the Phreeqc input file
pqi_mix = fileName+'_MIX.pqi'  # the name of the Phreeqc input file containing the MIX entries
pqo = fileName+'.pqo'  # the name of the Phreeqc output file
tsv = fileName+'.tsv'  # the name of the Phreeqc selected output file
# for the simulations without RMD:
pqi_noRMD = fileName_noRMD+'.pqi'
pqo_noRMD = fileName_noRMD+'.pqo'
tsv_noRMD = fileName_noRMD+'.tsv'

mixruns = ((D_im + alpha_L*v)*Delta_t) / (0.3333*Delta_x*Delta_x)
print("    shifts = "+str(shifts)+"  mixruns = {:.3f}".format(mixruns))
# there is no need to punch more than about one thousand time-step results
n_punch = 1000
punch_fr = round(shifts/n_punch) # punch frequency
if punch_fr < 1: punch_fr = 1
# for print, we nned 5 time-steps
print_fr = math.floor(shifts/5)
if print_fr < 1: print_fr = 1

# -------------------------------------------------------------------
# --- Modelling without RMD (rock matrix diffusion)
#     to compare results as in Fig. 3-2 of Applegate et al (2020)
RMD = False; n_z = 0

print("--- Modeling without RMD ---")
# check if the results from previous simulations exist and are ok...
ok = pqi_tsv_ok(pqi_noRMD,tsv_noRMD)
if ok:
    print('    file: "'+tsv_noRMD+'" already exists\n    reading previous results...')
    res_noRMD = pd.read_csv(tsv_noRMD,sep='\t',encoding='ANSI',decimal='.')
    # Note that the headings in the PHREEWC selected_output file contain spaces to the left
    # of the names, so they have to be trimmed or new names must be given.
    # Each line ends with a tab, so there is an extra column, the last one filled with NaN.
    res_noRMD.rename(columns=lambda x: x.strip(), inplace=True)
    # remove the last column containing 'NaN'
    res_noRMD.drop(res_noRMD.columns[len(res_noRMD.columns)-1],axis="columns", inplace=True)

else:
    print('    making input pqi-file: "'+pqi_noRMD+'"')
    input_string = make_input_string("",tsv_noRMD,RMD)
    # Save the input file in case you want to run PhreeqC outside Python
    with open(pqi_noRMD,'w') as f:
        f.writelines(input_string)
    print("    'pqi' file written.")
    answ = input(">>> Press <Enter> to continue and run Phreeqc ('q' to quit)...")
    if (answ == 'q') | (answ == 'Q'): raise SystemExit('quit...')
    # Make the Phreeqc calculations
    res_noRMD = runPhreeqc(input_string,pqo_noRMD)

# -------------------------------------------------------------------
# --- Modelling with RMD (rock matrix diffusion)
RMD = True
# The model type selects the type of rock-matrix model:
# - 0 (zero): 5 stagnant cells of equal size
# - 1 (one): 5 stagnant cells of increasing size
# - 2 (two): 10 stagnant cells of increasing size
n_z = 5 if (modelType == 0 or modelType ==1) else 10

print("--- Modelling with RMD ---")
# The dimensions of the rock matrix model
# The maximum diffusion depth is not given in the original references listed above.
# But the distance between two fractures from figuress in Weatherill et al (2008) and Watanabe and Kolditz (2015)
# is 0.1 m.  Fig.3 in Weatherill et al (2008) shows that the penetration during 4 days
# into the matrix is much smaller...
z_tot_max = 0.003 # m (max diffusion depth)
if z_tot_max > spacing /2.:
    raise SystemExit('z_tot_max is too large')
n_z1 = n_z + 1
Delta_z = [0]*n_z1
# The model type selects the type of rock-matrix model:
# - 0 (zero): 5 stagnant cells of equal size
# - 1 (one): 5 stagnant cells of increasing size
# - 2 (two): 10 stagnant cells of increasing size
if modelType == 0:  # stagnant cells of equal size
    for i in range(1, n_z1):
        Delta_z[i] = z_tot_max / n_z
elif modelType == 1: # Five (5) rock matrix cells of increasing size
    z_ = z_tot_max/15.
    for i in range(len(Delta_z)):
        Delta_z[i] = i * z_
elif modelType == 2: # 10 rock matrix cells of increasing size
    m_x = 1.1923738
    for i in range(1,11):
        Delta_z[i] = b * m_x**(i-1)
# check that the Delta_z add up to z_tot
sum_z = 0
for i in range(len(Delta_z)):
    sum_z = sum_z + Delta_z[i]
if abs(sum_z - z_tot_max) > 1e-5:
        print("Error: sum of Delta_z...")
        raise SystemExit

# Distace of cell center along matrix, used in USER_PUNCH only, for plotting
z = [None]*n_z1  # m
z[0] = b/2.
for i in range(1, n_z1):
    z[i] = z[i-1] + 0.5*(Delta_z[i-1]+Delta_z[i])
z[0] = z[0]/2.

psi = (v/D_im)*(b/phi_im)*(Delta_z[1]/Delta_x) # psi is the ratio of advective to diffusive solute flux in the fracture.
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
h = [None]*n_z # m
h[0] = (0.5*b + Delta_z[1] ) /2.
for i in range(1, n_z):
    h[i] = (Delta_z[i]+Delta_z[i+1]) /2.

# boundary condition factor
if fbc < 0 or fbc >2:
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
    print('    file: "'+pqi_mix+'" already exists...')
else:
    print("--- making input 'pqi' MIX file: '"+pqi_mix+"'")
    # Calculate the mixing factors:
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
    input_string = make_input_string(pqi_mix,tsv,RMD)
    # Save the input file in case you want to run PhreeqC outside Python
    with open(pqi,'w') as f:
        f.writelines(input_string)
    print("    'pqi' file written.")

    answ = input(">>> Press 'q' to quit or <Enter> to continue and run Phreeqc...")
    if (answ == 'q') | (answ == 'Q'): raise SystemExit('quit...')
    
    # Make the Phreeqc calculations
    res = runPhreeqc(input_string,pqo)

# analyse the results...
print('--- Results:')
# For the purposes of plotting, that is, to select the results for example midway along the column,
# find out the minimum, maximum and middle distances along the column
print("--- x-coordinate: finding min., mid and max.")
# x-coordinate: get the minimum, maximum and mid distance
xmn = 1e30
xmx = -1e20
distances = res.loc[(res['dist_x']>=0) & (res['Days']==0) & (res['soln']<=(n_cells)),'dist_x']
for x in distances:
    if(x < xmn):
        xmn = x
    if(x > xmx):
        xmx = x
mid = xmx/2.
dist = 1e20
xmd = -1; Dx2 = Delta_x/2.
for x in distances:
    if(abs((x-Dx2)-mid) < dist):
        xmd = x
        dist = abs((x-Dx2)-mid)
print('    X-distance:  min='+str(xmn)+',  mid='+str(xmd)+',  max='+str(xmx))
if (xmn >1e25) | (xmd <0) | (xmx < 0):
    print('Error, one of these three values is wrong...')
    raise SystemExit

# For the purposes of plotting, that is, to select the results at a given simulation time...
time_1 = 0.5; time_2 = 1.5; time_3 = endTime  # days
print("--- finding results close to "+str(time_1)+", "+str(time_2)+" and "+str(time_3)+" days")
t_1 = -1; t_2 = -1; t_3 = -1
dist_1 = 1e20; dist_2 = 1e20; dist_3 = 1e20
times = res.loc[(res['soln']==1),'Days']
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
print('    '+str(time_1)+' days = '+str(t_1)+'\n    '+str(time_2)+' days = '+str(t_2)+'\n    '+str(time_3)+' days = '+str(t_3))
if (t_1 <0) | (t_2 <0) | (t_3 < 0):
    print('Error, one of these three values is negative...')
    raise SystemExit

print("--- x-coordinate: finding 0.76m")
# x-coordinate: get exact distance closest to 0.76 m
dx_5L = res.loc[(res['dist_x']>0) & (res['Days']==0) & (res['soln']>0) & (res['soln']<(n_cells+1)),'dist_x']
target = 0.76 # m to compare with the results in Fig.4b of Weatherill et al (2008)
dist = 1e20
x076 = -1
for x in dx_5L:
    if(abs(x-target) < dist):
        x076 = x
        dist = abs(x-target)
print('    X-distance:  0.76 ='+str(x076)+' m')
if (x076 <0):
    print('Error, did not find 0.76m...')
    raise SystemExit

print("--- making plots.")
# 'include' is the frequency of markers on the plots.
# there should be around n_cells steps in the output file.
# but it is enough to plot around 10 markers
incl_x = round(n_cells/10)
# it is enough to plot around 15 markers
incl_t = round(shifts/15)
plot(res_noRMD,False,incl_x,incl_t)
plot(res,True,incl_x,incl_t)
plotMtrx(res)
