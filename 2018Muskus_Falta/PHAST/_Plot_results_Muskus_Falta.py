# -*- coding: utf-8 -*-
"""
Script to plot the results from PHAST output file in 'subFolder'
Plots are saved in folder: subFolder+"plots"
@author: Ignasi
"""
import sys
import matplotlib.pyplot as plt
import pandas as pd
sys.path.insert(0, 'D:/OneDrive/Documents/Geokemi/Modelling/PhreeqC/1D_transport/Analytical_eqn/CraFlush')
from craflush import CraFlush

# subfolder with the output file from 'PHAST'
tsv = 'Muskus_Falta_1.chem.xyz.tsv'

# Data used in the PHAST calculation
alpha_L = 0.001 # dispersivity used in PHAST
alpha_A = 0.1 # dispersivity used in the analytical model
D_0 = 1e-9 # m2/s
tau = 0.1
D_im = D_0*tau # m2/s, diffusion coefficient in the rock matrix ("im" for immobile)
b = 0.0001 # m (fracture aperture)
spacing = 0.50 # m the distance between two fractures from Fig.2 in Weatherill et al (2008) and Watanabe and Kolditz (2015)
phi_m = 1 # porosity in the fracture
phi_im = 0.1 # (rock matrix porosity; "im" for immobile)
K_x = 0.0158444 # m/s
Delta_head = 0.01 # m
Delta_t = 30*86400 # s = 30 days

L = 50 # m, the length of the rock model
# the interstitial groundwater velocity is calculated in PHAST
# eqn.(D.3) in the user's guide
v = (K_x/phi_m) * (Delta_head/L)
print('    velocity = {:.3e}'.format(v)+' m/s;  {:.4f}'.format(v*86400)+' m/day;  {:.1f}'.format(v*365.25*86400)+' m/yr')
# v = 0.75 /(24*60*60) # 0.75 m/day = 8.68056e-6 m/s (groundwater velocity)


def c_c0_S(tyr,X,Z,N,Rp):
    # Implements the analytical equation for matrix diffusion from Sudicky et al.
    # N decides the number of terms in the EPAL summation
    # units: x and z in 'm', tyr in 'years'
    t = tyr * 365.25*86400 # convert days to 's'
    # Data for the PhreeqC 'model':
    # b = fracture aperture
    C0 = 0.     # Source concentration at infinite time
    Cim = 0.    # Initial concentration in matrix and fractures
    alpha = alpha_A  # alpha_num + alpha_L  # Fracture dispersivity (m)
    sep = spacing    # Fracture spacing (m)
    theta = phi_im  # Matrix porosity
    D0 = 1.e-9    # Diffusion coefficient in water (m2/s)
    tau = D_im/D0    # Matrix tortuosity
    R = 1.      # Fracture retardation factor
    # Rp        # Matrix retardation factor
    THalf = 0.  # non-decaying solute
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
    modelTxt = '5 stagnant layers of\n    increasing size$_{}$'
    annotateTxt = (
        modelTxt+"\n"+
        "$n\mathsf{_{cells}}$ = "+str(n_cells)+"\n"+
        "$\\alpha_L$ = "+"{:.3f}".format(alpha_L)+" m\n"+
        "$D_{im}$ = "+"{:.2e}".format(D_im)+" m$^2$/s\n")  # used in PHREEQC
    annotateTxt = annotateTxt +("Markers: PHAST\n"+
        "Lines: analytical model\n"+
        "  Sudicky et al.\n"+
        "using: $\\alpha$ = "+"{:.3f}".format(alpha_A)+" m\n"+
        "$D'$ = "+"{:.2e}".format(D_im)+" m$^2$/s\n"+
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
    Yr_1 = rsl[(rsl["time"]==t_1) & (rsl["x"]>=0) & (rsl["z"]==0)].iloc[::include,:]
    Yr_2 = rsl[(rsl["time"]==t_2) & (rsl["x"]>=0) & (rsl["z"]==0)].iloc[::include,:]
    Yr_3 = rsl[(rsl["time"]==t_3) & (rsl["x"]>=0) & (rsl["z"]==0)].iloc[::include,:]
    Yr_4 = rsl[(rsl["time"]==t_4) & (rsl["x"]>=0) & (rsl["z"]==0)].iloc[::include,:]

    xLabel = 'Distance x along fracture (m)'

    elements = ['Na','NO3']
    for i in range(2):
        # The retardation factor in the rock matrix
        ele = elements[i]
        Rp = 2. if (ele == 'Na') else 1.  # Matrix retardation factor

        #ax.axis([0, L, 0, 1.])
        ax[i].plot(Yr_1["x"],Yr_1[ele],c=clr1,marker=mrk1,markersize=9,linestyle='None',label="{:.0f}".format(t_1)+" yr")
        ax[i].plot(Yr_2["x"],Yr_2[ele],c=clr2,marker=mrk2,markersize=6,linestyle='None',label="{:.0f}".format(t_2)+" yr")
        ax[i].plot(Yr_3["x"],Yr_3[ele],c=clr3,marker=mrk3,markersize=9,linestyle='None',label="{:.0f}".format(t_3)+" yr")
        ax[i].plot(Yr_4["x"],Yr_4[ele],c=clr4,marker=mrk4,markersize=6,linestyle='None',label="{:.0f}".format(t_4)+" yr")

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
            y1[j] = c_c0_S(t_1,distx[j],0,11,Rp)
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

    plotName="Muskus_Falta_PHAST-1.png"
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
            xMn=rsl[(rsl["time"]==tYr) & (rsl["x"]==xmn)]
            xMd=rsl[(rsl["time"]==tYr) & (rsl["x"]==xmd)]
            xMx=rsl[(rsl["time"]==tYr) & (rsl["x"]==xmx)]
            axsLeft[j].plot(xMn['z'],xMn[ele],c=clr1,marker=mrk1,linestyle='None',markersize=9)
            axsLeft[j].plot(xMd['z'],xMd[ele],c=clr2,marker=mrk2,linestyle='None',markersize=6)
            axsLeft[j].plot(xMx['z'],xMx[ele],c=clr3,marker=mrk3,linestyle='None',markersize=9)
            axsLeft[j].set_ylim([-0.1, 1.1])
            axsLeft[j].set_ylabel(yLabel)
            axsLeft[j].set_title("{:.0f}".format(tYr)+" yr", x=0.7, y=0.63, fontsize=10)
            # calculate and plot the analytical (theoretical) curves
            distz = [0]*21
            xMnA = [1]*21;  xMdA = [1]*21;  xMxA = [1]*21
            N = 11
            for k in range(0,21):
                distz[k]= z_tot_max*k/20
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

    plotName = "Muskus_Falta_PHAST-2.png"
    plt.savefig('plots/'+plotName)
    plt.show()
    
    plt.close()

def plot2(rsl,incl_x):
    # Plot results.  'rsl' is the DataFrame with the data to plot.
    # The figures are saved in folder 'plots'
    # Not all points are plotted, the variable incl_x is used
    # to select how many points are plotted.
    if incl_x < 1: incl_x = 1
    # 'element' is a chemical element name in the column header of the DataFrame

    # The retardation factor in the rock matrix
    R = 2 if element == 'Na' else 1

    print('--- making plots for "'+element+'"')

    # Two figures will be created:
    #    1- concentrations along the fracture at different times
    #    2- concentrations into the rock matrix for different times
    #       and at three distances from the inflow

    # First create a figure with one 'axes' to plot
    # concentrations along the fracture at different times
    fig, ax=plt.subplots(1,layout='constrained',figsize=(5.12,3.84),dpi=100)

    # t_1, t_2 ... are times at which to plot results, for example, 50, 100 ... years
    D_1 = rsl[(rsl["time"]==t_1) & (rsl["x"]>=0) & (rsl["z"]==0)].iloc[::incl_x,:]
    D_2 = rsl[(rsl["time"]==t_2) & (rsl["x"]>=0) & (rsl["z"]==0)].iloc[::incl_x,:]
    D_3 = rsl[(rsl["time"]==t_3) & (rsl["x"]>=0) & (rsl["z"]==0)].iloc[::incl_x,:]
    D_4 = rsl[(rsl["time"]==t_4) & (rsl["x"]>=0) & (rsl["z"]==0)].iloc[::incl_x,:]
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
    xLabel = 'Distance x along fracture (m)'
    yLabel = element+'/'+element+'$_0$'
    ax.set_ylabel(yLabel)
    ax.set_xlabel(xLabel)
    # ax[0].axis([0, 3, 0, 1.])
    ax.plot(D_1["x"],D_1[element],c=clr1,marker=mrk1,markersize=9,linestyle='None',label="{:.2f}".format(t_1)+" yr")
    ax.plot(D_2["x"],D_2[element],c=clr2,marker=mrk2,markersize=6,linestyle='None',label="{:.2f}".format(t_2)+" yr")
    ax.plot(D_3["x"],D_3[element],c=clr3,marker=mrk3,markersize=9,linestyle='None',label="{:.2f}".format(t_3)+" yr")
    ax.plot(D_4["x"],D_4[element],c=clr4,marker=mrk4,markersize=6,linestyle='None',label="{:.2f}".format(t_4)+" yr")
    # calculate and plot the analytical (theoretical) curves
    distx = [0]*51
    y1 = [1]*51; y2 = [0]*51; y3 = [0]*51; y4 = [0]*51
    for i in range(1,51):
        distx[i] = L*i/50
        N = 11
        y1[i] = c_c0_S(t_1,distx[i],0,N,R)
        y2[i] = c_c0_S(t_2,distx[i],0,N,R)
        y3[i] = c_c0_S(t_3,distx[i],0,N,R)
        y4[i] = c_c0_S(t_4,distx[i],0,N,R)
    ax.plot(distx,y1,c=clr1,marker="None")
    ax.plot(distx,y2,c=clr2,marker="None")
    ax.plot(distx,y3,c=clr3,marker="None")
    ax.plot(distx,y4,c=clr4,marker="None")
    
    ax.legend(fontsize='small', bbox_to_anchor=(1.07,0.), loc='lower left', borderaxespad=0)

    txt = annotTxt()
    ax.annotate(txt, xy=(1.07, 1.), xycoords='axes fraction',
            fontsize=8, verticalalignment='top')

    plotName = "_"+element+"_Fract.png"
    plt.savefig('plots/'+plotName)
    plt.show()
    plt.close()

    # 2- Create a figure containing two subfigures in two columns.
    # - in the left subfigure three plots ('axes') showing concentrations
    #   versus distance into the matrix
    # - to the right a text box with the parameters for the simulation
    #   and a legend (the same one for all plot)
    print('--- making plots along the matrix for "'+element+'"')
    fig = plt.figure(layout='constrained',figsize=(5.12,3.84),dpi=100)
    subfigs = fig.subfigures(1,2, wspace = 0.04, width_ratios=[2, 1])

    # Left subfigure:
    # Three plots (axes) in a single column
    axsLeft = subfigs[0].subplots(3, 1, sharex=True)

    # Create a plot with concentrations versus distance into the matrix for
    # three distances along the fracture and three different
    # simulation times (t_1, t_2 and t_3).
    clr1 = "#DC267F" # mauve
    clr2 = "#FE6100"  # orange
    clr3 = "#FFB000"  # yellow
    mrk1 = "v"  # triangle down
    mrk2 = "D"  # diamond
    mrk3 = "^" # triangle up
    xLabel = 'Distance z in matrix (m)'
    yLabel = element+'/'+element+'$_0$'

    # t_1, t_2 ... are times at which to plot results, for example, 50, 100 ... years
    years = [t_1,t_2,t_3]
    for i in range(0,3):        
        tYr = years[i]
        xMn=res[(res["time"]==tYr) & (res["x"]==xmn)]
        xMd=res[(res["time"]==tYr) & (res["x"]==xmd)]
        xMx=res[(res["time"]==tYr) & (res["x"]==xmx)]
        # ax[0].axis([0, 3, 0, 1.])
        axsLeft[i].set_ylim([-0.1, 1.1])
        #ax[i].yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
        axsLeft[i].plot(xMn["z"],xMn[element],c=clr1,marker=mrk1,markersize=9,linestyle='None',label="x= {:.2f}".format(Delta_x/2.)+" m")
        axsLeft[i].plot(xMd["z"],xMd[element],c=clr2,marker=mrk2,markersize=6,linestyle='None',label="x= {:.2f}".format(xmd)+" m")
        axsLeft[i].plot(xMx["z"],xMx[element],c=clr3,marker=mrk3,markersize=9,linestyle='None',label="x= {:.2f}".format(xmx)+" m")
        # calculate and plot the analytical (theoretical) curves
        distz = [0]*21
        yXmn = [1]*21;  yXmd = [1]*21;  yXmx = [1]*21
        N = 11
        for j in range(0,21):
            distz[j]= z_tot_max*j/20
            yXmn[j] = c_c0_S(tYr,xmn,distz[j],N,R)
            yXmd[j] = c_c0_S(tYr,xmd,distz[j],N,R)
            yXmx[j] = c_c0_S(tYr,xmx,distz[j],N,R)
        axsLeft[i].plot(distz,yXmn,c=clr1,marker="None")
        axsLeft[i].plot(distz,yXmd,c=clr2,marker="None")
        axsLeft[i].plot(distz,yXmx,c=clr3,marker="None")
        if i == 1: axsLeft[i].set_xlabel(xLabel)
        axsLeft[i].set_ylabel(yLabel)
        axsLeft[i].set_title("{:.1f}".format(tYr)+" yr", x=0.7, y=0.63, fontsize=10)

    # Subfigure to the right:
    # Two subplots (axes) in a single column.
    axsRight = subfigs[1].subplots(2, 1)

    # do not show the axes
    axsRight[0].axis('off')
    axsRight[1].axis('off')

    # First 'axes': a text box listing parameters used.
    axsRight[0].annotate(txt, fontsize=8,
             xy=(0.06,1.), xycoords='subfigure fraction',  verticalalignment='top',
             #bbox=dict(boxstyle="square,pad=0.5",linewidth=1,facecolor='white')
             )

    # Second 'axes': an empty plot defining the legend texts
    axsRight[1].plot([],[],c=clr1,marker=mrk1,linestyle='None',markersize=9,label="{:.2f}".format(Delta_x/2.)+" m")
    axsRight[1].plot([],[],c=clr2,marker=mrk2,linestyle='None',markersize=6,label="{:.2f}".format(xmd)+" m")
    axsRight[1].plot([],[],c=clr3,marker=mrk3,linestyle='None',markersize=9,label="{:.2f}".format(xmx)+" m")
    axsRight[1].legend(fontsize='small', loc=(0.,0.3), borderaxespad=0) #, bbox_to_anchor=(0.,2.)

    plotName = "_"+element+"-Mtrx.png"
    plt.savefig('plots/'+plotName)
    plt.show()
    plt.close()

# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------

print('--- Reading file "'+tsv+'"...')
res=pd.read_csv(tsv,sep='\t',encoding='ANSI',decimal='.')
# remove blank space from column texts
res.rename(columns=lambda x: x.strip(), inplace=True)
# remove last column with 'NaN'
res.drop(res.columns[len(res.columns)-1],axis="columns", inplace=True)

Delta_x = res.loc[(res['z']==0) & (res['time']==60),'x'][1]  # 60 secs is the first output

# x-coordinate: get three distances, the minimum, 15 m and in between
xmn = 1e30
xmx = -1e20
distances = res.loc[(res['x']>=0) & (res['time']==60) & (res['z']<=0),'x']
for x in distances:
    if(x < xmn):
        xmn = x
    if(x > xmx) & (x < 15.001):
        xmx = x
mid = xmx/2.
dist = 1e20
xmd = -1
for x in distances:
    if(abs(x-mid) < dist):
        xmd = x
        dist = abs(x-mid)
print('    X-distance:  min='+str(xmn)+',  mid='+str(xmd)+',  max='+str(xmx))
if (xmn >1e25) | (xmd <0) | (xmx < 0):
    print('Error, one of these three values is wrong...')
    raise SystemExit

depths = res.loc[(res['z']>=0) & (res['time']==60) & (res['x']==0),'z']
zmx = -1e20
for z in depths:
    if(z > zmx):
        zmx = z
print('    max Z = '+str(zmx))
if (zmx < 0):
    print('Error, value is wrong...')
    raise SystemExit
z_tot_max = zmx # m (max diffusion depth)

# convert seconds to years
res['time']=res['time'].apply(lambda x: x/(365.25*86400))

# time: get the closest time to 50, 100, 150 and 200 years (for plotting)
time_1 = 49.5; time_2 = 100; time_3 = 150; time_4 = 200
print("--- finding results close to "+str(time_1)+", "+str(time_2)+", "+str(time_3)+" and "+str(time_4)+" years")
t_1 = -1; t_2 = -1; t_3 = -1; t_4 = -1
dist_1 = 1e20; dist_2 = 1e20; dist_3 = 1e20; dist_4 = 1e20
times = res.loc[(res['z']==0) & (res['x']==0),'time']
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
if (t_1 <0) | (t_2 <0) | (t_3 < 0) | (t_4 < 0):
    print('Error, one of these values is negative...')
    raise SystemExit

n_times = times.shape[0]-1
print('    n_times='+str(n_times))

res0 = res[(res["time"]==t_1) & (res["x"]>=0) & (res["z"]==0)]
n_cells = res0.shape[0]-1
Delta_x = L/n_cells
print('    n_cells='+str(n_cells)+'   Delta_x={:.3f}'.format(Delta_x))

# estimated numerical spatial and temporal dispersivity
#             eqn.(D.7)        eqn.(D.10)
alpha_num = (Delta_x/2.) + (v*Delta_t/(2.)) 
print('    given dispersivity: '+str(alpha_L)+' m,  numerical dispersivity: '+"{:.4f}".format(alpha_num)+' m')

Pe = Delta_x / alpha_L # Peclet number, eqn.(D.5)
print('    Pe = {:.3f}'.format(Pe)+' (=Delta_x/alpha_L) Peclet nbr must be < 2')
if Pe >= 2: print('        decrease "Delta_x" or increase "alpha"')

Cr = v * Delta_t / Delta_x # Courant number, eqn.(D.8)
print('    Cr = {:.3f}'.format(Cr)+' (=v*Delta_t/Delta_x) Courant nbr should be <= Pe')
if Cr > Pe: print('        decrease either "v" or "Delta_t" or increase "Delta_x"')

# 'include' is the frequency of markers on the plots.
# there should be around n_cells steps in the output file.
# but it is enough to plot around 15 markers
include_x = max(1,round(n_cells/15))
plot(res,include_x)
