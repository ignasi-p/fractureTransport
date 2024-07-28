# -*- coding: utf-8 -*-
"""
Script to plot the results from a PHAST output file
Plots are saved in folder: subFolder+"plots"
@author: Ignasi
"""
import sys
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import pandas as pd
sys.path.insert(0, 'D:/OneDrive/Documents/Geokemi/Modelling/PhreeqC/1D_transport/Analytical_eqn/CraFlush')
from craflush import CraFlush

# Data used in the PHAST calculation
alpha_L = 0.01 # dispersivity
alpha_A = 0.01 # dispersivity in the analytical model
D_im = 1e-12 # m2/s, diffusion coefficient in the rock matrix ("im" for immobile)
b = 0.00012 # m (fracture aperture)
spacing = 0.1 # m the distance between two fractures from Fig.2 in Weatherill et al (2008) and Watanabe and Kolditz (2015)
phi_m = 1 # porosity in the fracture
phi_im = 0.35 # (rock matrix porosity; "im" for immobile)
K_x = 0.00069445 # m/s
Delta_head = 0.01 # m
Delta_t = 3600 # s = 1 hour

L = 0.8 # m, the length of the rock model
# the interstitial groundwater velocity is calculated in PHAST
# eqn.(D.3) in the user's guide
v = (K_x/phi_m) * (Delta_head/L)
print('    velocity = {:.3e}'.format(v)+' m/s;  {:.4f}'.format(v*86400)+' m/day;  {:.1f}'.format(v*365.25*86400)+' m/yr')
# v  = 1.77e-5 m/s (groundwater velocity)

def c_c0_S(tdays,X,Z,N):
    # Implements the analytical equation for matrix diffusion from Sudicky et al.
    # N decides the number of terms in the EPAL summation
    # theta = Matrix porosity
    # units: x and z in 'm', tdays in 'days'
    t = tdays * 86400 # convert days to 's'
    # Data for the PhreeqC 'model':
    # b = fracture aperture
    C0 = 1.     # Source concentration at infinite time
    Cim = 0.    # Initial concentration in matrix and fractures
    alpha = alpha_A  # Fracture dispersivity (m)
    theta = phi_im
    D0 = 1.e-9    # Diffusion coefficient in water (m2/s)
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
    cc0, Cs = CraFlush(C0,Cim,v,alpha,b,spacing,theta,tau,D0,R,Rp,THalf,
                             X,Z,t, Ns,Ts,Css,AlphR,relErr,N)
    if(cc0 > 1.): cc0 =1.
    if(cc0 < 0.): cc0 =0.
    return cc0

def annotTxt():
    # modelA: 0 = Neretnieks (1980);  1 = Sudicky et al.
    modelTxt = 'No rock matrix diffusion.$_{}$'
    modelTxt = '5 stagnant layers of\n    increasing size$_{}$'
    annotateTxt = (
        modelTxt+"\n"+
        "$n\mathsf{_{cells}}$ = "+str(n_cells)+"\n"+
        "$\\alpha_L$ = "+"{:.3f}".format(alpha_L)+" m\n"+
        "$D_{im}$ = "+"{:.2e}".format(D_im)+" m$^2$/s\n")  # used in PHREEQC
    annotateTxt = annotateTxt +("Markers: PHREEQC$_{ }$\n"+
        "Lines: analytical model")
    annotateTxt = annotateTxt + "\n  Sudicky et al."
    annotateTxt = annotateTxt + "$^{ }$\n"+"using: $\\alpha_{tot}$ = "+"{:.3f}".format(alpha_A)+" m\n"
    annotateTxt = annotateTxt + "$D'$ = "+"{:.2e}".format(D_im)+" m$^2$/s"
    return annotateTxt

def plot(rsl,incl_x,incl_t):
    # Plot results into a single figure containing two plots.
    # The figure is saved in folder 'plots'
    # Not all points are plotted, the variables incl_x and incl_t are used to select how many points are plotted.
    if incl_x < 1: incl_x = 1
    print("--- making plots")

    # Create a figure containing two rows.
    fig, ax=plt.subplots(2,layout='constrained',figsize=(5.12,3.84),dpi=100)

    D_1 = rsl[(rsl["time"]==t_1) & (rsl["x"]>=0) & (rsl["x"]<L) & (rsl["z"]==0)].iloc[::incl_x,:]
    D_2 = rsl[(rsl["time"]==t_2) & (rsl["x"]>=0) & (rsl["x"]<L) & (rsl["z"]==0)].iloc[::incl_x,:]
    D_3 = rsl[(rsl["time"]==t_3) & (rsl["x"]>=0) & (rsl["x"]<L) & (rsl["z"]==0)].iloc[::incl_x,:]
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
    ax[0].plot(D_1["x"],D_1["A"],c=clr1,marker=mrk1,markersize=9,linestyle='None',label="{:.1f}".format(t_1)+" days")
    ax[0].plot(D_2["x"],D_2["A"],c=clr2,marker=mrk2,markersize=6,linestyle='None',label="{:.1f}".format(t_2)+" d")
    ax[0].plot(D_3["x"],D_3["A"],c=clr3,marker=mrk3,markersize=9,linestyle='None',label="{:.1f}".format(t_3)+" d")
    # calculate and plot the analytical (theoretical) curves
    distx = [0]*21
    y1 = [1]*21
    y2 = [1]*21
    y3 = [1]*21
    N = 11
    for i in range(1,21):
        distx[i] = L*i/20
        y1[i] = c_c0_S(t_1,distx[i],0,N)
        y2[i] = c_c0_S(t_2,distx[i],0,N)
        y3[i] = c_c0_S(t_3,distx[i],0,N)
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
    d=rsl[(rsl["x"]==x076) & (rsl["z"]==0)].iloc[::incl_t,:]
    #ax[1].axis([0, 0.25, 0, 1.])
    ax[1].set_ylim([-0.1, 1.1])
    ax[1].plot(d['time'],d['A'],c=clr2,marker=mrk2,linestyle='None',markersize=4,label="{:.2f}".format(x076)+" m")
    ax[1].set_ylabel(yLabel)
    ax[1].set_xlabel(xLabel)
    ax[1].set_title("$x=$"+"{:.2f}".format(x076)+" m", x=0.7, y=0.7)
    # plot the analytical (theoretical) curve
    timex = [0]*101
    y1 = [0]*101
    X = 0.76 # distance along fracture (m)
    N = 11
    for i in range(1,101):
        timex[i] = 4*i/100
        y1[i] = c_c0_S(timex[i],X,0,N)
    ax[1].plot(timex,y1,c=clr2,marker="")

    txt = annotTxt()
    ax[1].annotate(txt, xy=(1.07, 1.), xycoords='axes fraction',
            fontsize=8, verticalalignment='top')

    plotName = "Tang_Fract"+".png"
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
        dMn=res[(res["time"]==d) & (res["x"]==xmn)]
        dMd=res[(res["time"]==d) & (res["x"]==xmd)]
        dMx=res[(res["time"]==d) & (res["x"]==xmx)]
        # ax[0].axis([0, 3, 0, 1.])
        ax[i].yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
        ax[i].plot(dMn["z"],dMn["A"],c=clr1,marker=mrk1,markersize=9,linestyle='None',label="x= {:.2f}".format(Delta_x/2.)+" m")
        ax[i].plot(dMd["z"],dMd["A"],c=clr2,marker=mrk2,markersize=6,linestyle='None',label="x= {:.2f}".format(xmd)+" m")
        ax[i].plot(dMx["z"],dMx["A"],c=clr3,marker=mrk3,markersize=9,linestyle='None',label="x= {:.2f}".format(xmx)+" m")
        # calculate and plot the analytical (theoretical) curves
        distz = [0]*21
        yXmn = [1]*21;  yXmd = [1]*21;  yXmx = [1]*21
        N = 11
        for j in range(0,21):
            distz[j]= z_tot_max*j/20
            yXmn[j] = c_c0_S(d,xmn,distz[j],N)
            yXmd[j] = c_c0_S(d,xmd,distz[j],N)
            yXmx[j] = c_c0_S(d,xmx,distz[j],N)
        ax[i].plot(distz,yXmn,c=clr1,marker="None")
        ax[i].plot(distz,yXmd,c=clr2,marker="None")
        ax[i].plot(distz,yXmx,c=clr3,marker="None")
        if i == 1: ax[i].set_xlabel(xLabel)
        ax[i].set_ylabel(yLabel)
        ax[i].set_title("$t=$"+"{:.2f}".format(d)+" days", x=0.7, y=0.7)

    ax[1].legend(fontsize='small', bbox_to_anchor=(1.07,0.), loc='lower left', borderaxespad=0)

    txt = annotTxt()
    ax[0].annotate(txt, xy=(1.07, 1.), xycoords='axes fraction',
            fontsize=8, verticalalignment='top')

    plotName = "Tang_Matrix"+".png"
    plt.savefig('plots/'+plotName)
    plt.show()
    
    plt.close()

# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------

res=pd.read_csv('Tang_1.chem.xyz.tsv',sep='\t',encoding='ANSI',decimal='.')
# remove blank space from column texts
res.rename(columns=lambda x: x.strip(), inplace=True)
# remove last column with 'NaN'
res.drop(res.columns[len(res.columns)-1],axis="columns", inplace=True)

Delta_x = res.loc[(res['z']==0) & (res['time']==60),'x'][1]  # 60 secs is the first output

# x-coordinate: get the minimum, maximum and mid distance
print("--- x-coordinate: finding min, mid and max")
xmn = 1e30
xmx = -1e20
distances = res.loc[(res['x']>=0) & (res['time']==60) & (res['z']<=0),'x']
for x in distances:
    if(x < xmn):
        xmn = x
    if(x > xmx):
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

print("--- x-coordinate: finding 0.76m")
# x-coordinate: get exact distance closest to 0.76 m
target = 0.76 # m to compare with the results in Fig.4b of Weatherill et al (2008)
dist = 1e20
x076 = -1
for x in distances:
    if(abs(x-target) < dist):
        x076 = x
        dist = abs(x-target)
print('    X-distance:  0.76 ='+str(x076)+' m')
if (x076 <0):
    print('Error, did not find 0.76m...')
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

# convert seconds to days
res['time']=res['time'].apply(lambda x: x/86400)

L = xmx # m, the length of the rock model
# the interstitial groundwater velocity is calculated in PHAST
# eqn.(D.3) in the user's guide
v = (K_x/phi_m) * (Delta_head/L)
print('    velocity = '+str(v)+' m/s;  '+str(v*86400)+' m/day')
# v = 0.75 /(24*60*60) # 0.75 m/day = 8.68056e-6 m/s (groundwater velocity)

# For the purposes of plotting, that is, to select the results at a given simulation time...
time_1 = 0.5; time_2 = 1.5 # days
print("--- finding results close to "+str(time_1)+", "+str(time_2)+" and max. days")
t_1 = -1; t_2 = -1; t_3 = -1
dist_1 = 1e20; dist_2 = 1e20; dist_3 = 1e20
days = res.loc[(res['z']==0) & (res['x']==0),'time']
for t in days:
    if(abs(t-time_1) < dist_1):
        t_1 = t
        dist_1 = abs(t-time_1)
    if(abs(t-time_2) < dist_2):
        t_2 = t
        dist_2 = abs(t-time_2)
    if(t > t_3):
        t_3 = t
print('    '+str(time_1)+' days = '+str(t_1)+'\n    '+str(time_2)+' days = '+str(t_2)+'\n    max days = '+str(t_3))
if (t_1 <0) | (t_2 <0) | (t_3 < 0):
    print('Error, one of these three values is negative...')
    raise SystemExit

n_times = days.shape[0]-1
print('    n_times='+str(n_times))

res0 = res[(res["time"]==t_1) & (res["x"]>=0) & (res["z"]==0)]
n_cells = res0.shape[0]-1
print('    n_cells='+str(n_cells))

Delta_x = L/n_cells
# estimated numerical spatial and temporal dispersivity
#             eqn.(D.7)        eqn.(D.10)
alpha_num = (Delta_x/2.) + (v*Delta_t/(2.)) 
print('    given dispersivity: '+str(alpha_L)+' m,  numerical dispersivity: '+"{:.4f}".format(alpha_num)+' m')

Pe = Delta_x / alpha_L # Peclet number, eqn.(D.5)
print('    Pe = '+str(Pe)+' (=Delta_x/alpha_L) Peclet nbr must be < 2')

Cr = v * Delta_t / Delta_x # Courant number, eqn.(D.8)
print('    Cr = '+str(Cr)+' (=v*Delta_t/Delta_x) Courant nbr should be <= Pe')

# 'include' is the frequency of markers on the plots.
# there should be around n_cells steps in the output file.
# but it is enough to plot around 15 markers
incl_x = max(1,round(n_cells/15))
incl_t = max(1,round(n_times/20))
plot(res,incl_x,incl_t)
plotMtrx(res)
