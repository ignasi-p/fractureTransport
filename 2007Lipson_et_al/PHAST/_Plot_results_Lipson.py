# -*- coding: utf-8 -*-
"""
Script to plot the results from a PHAST output file
Plots are saved in folder "plots"
@author: Ignasi
"""
import sys
import matplotlib.pyplot as plt
import pandas as pd
sys.path.insert(0, 'D:/OneDrive/Documents/Geokemi/Modelling/PhreeqC/1D_transport/Analytical_eqn/CraFlush')
from craflush import CraFlush

# the results from PHAST
tsv = 'Lipson_1.chem.xyz.tsv'

# Data used in the PHAST calculation
alpha_L = 0.001 # dispersivity
alpha_A = 0.001 # dispersivity used in the analytical model
D_0 = 1e-9 # m2/s
tau = 0.44  # this gives D_im = 4.4e-10 # m2/s, diffusion coefficient in the rock matrix ("im" for immobile)
D_im = D_0*tau # m2/s, diffusion coefficient in the rock matrix ("im" for immobile)
b = 0.0004  # m (fracture aperture)
spacing = 0.5 # m the distance between two fractures from Fig.2 in Weatherill et al (2008) and Watanabe and Kolditz (2015)
phi_m = 1 # porosity in the fracture
phi_im = 0.2 # (rock matrix porosity; "im" for immobile)
K_x = 0.00531915 # m/s
Delta_head = 0.01 # m
Delta_t = 86400 # s = 1 day

L = 3 # m, the length of the rock model
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
    alpha = alpha_A # alpha_num + alpha_L  # Fracture dispersivity (m)
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
    modelTxt = '5 stagnant layers of\n    increasing size$_{}$'
    annotateTxt = (
        modelTxt+"\n"+
        "$n\mathsf{_{cells}}$ = "+str(n_cells)+"\n"+
        "$\\alpha_L$ = "+"{:.3f}".format(alpha_L)+" m\n"+
        "$D_{im}$ = "+"{:.2e}".format(D_im)+" m$^2$/s\n")  # used in PHREEQC
    annotateTxt = annotateTxt +("Markers: PHAST$_{ }$\n"+
        "Lines: analytical model\n")
    annotateTxt = annotateTxt + "  of Sudicky et al.$^{ }$\n"
    annotateTxt = annotateTxt + "using: $\\alpha_{tot}$ = "+"{:.3f}".format(alpha_A)+" m\n"
    annotateTxt = annotateTxt + "$D'$ = "+"{:.2e}".format(D_im)+" m$^2$/s"
    return annotateTxt

def plot(rsl,incl_x):
    # Plot results into a single figure containing two plots.
    # The figure is saved in folder 'plots'
    # Not all points are plotted, the variable incl_x is used to select how many points are plotted.
    if incl_x < 1: incl_x = 1
    print("--- making plots along fracture and matrix.")

    # Create a figure containing two subfigures in two columns.
    # - in the left two plots ('axes') showing concentrations versus distance
    #   and into the matrix
    # - in the right the two legends (one for each plot)
    #   and a text box with the parameters for the simulation
    fig = plt.figure(layout='constrained',figsize=(5.12,3.84),dpi=100)
    subfigs = fig.subfigures(1,2, wspace = 0.04, width_ratios=[2, 1])

    # Left subfigure:
    # Two plots (axes) in a single column
    axsLeft = subfigs[0].subplots(2, 1)
    yLabel = 'C/C$_0$'

    # First plot with concentrations versus distance along the fracture for
    # three different simulation times (t_1, t_2 and t_3).
    D_1 = rsl[(rsl["time"]==t_1) & (rsl["x"]>=0) & (rsl["z"]==0)].iloc[::incl_x,:]
    D_2 = rsl[(rsl["time"]==t_2) & (rsl["x"]>=0) & (rsl["z"]==0)].iloc[::incl_x,:]
    D_3 = rsl[(rsl["time"]==t_3) & (rsl["x"]>=0) & (rsl["z"]==0)].iloc[::incl_x,:]

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
    axsLeft[0].set_ylim([-0.05, 1.05])
    axsLeft[0].plot(D_1["x"],D_1["A"],c=clr1,marker=mrk1,markersize=9,linestyle='None')
    axsLeft[0].plot(D_2["x"],D_2["A"],c=clr2,marker=mrk2,markersize=6,linestyle='None')
    axsLeft[0].plot(D_3["x"],D_3["A"],c=clr3,marker=mrk3,markersize=9,linestyle='None')
    axsLeft[0].set_ylabel(yLabel)
    axsLeft[0].set_xlabel(xLabel)
    # calculate and plot the analytical (theoretical) curves
    distx = [0]*21
    y1 = [1]*21
    y2 = [1]*21
    y3 = [1]*21
    for i in range(1,21):
        distx[i] = L*i/20
        N = 11
        y1[i] = c_c0_S(t_1,distx[i],0,N)
        y2[i] = c_c0_S(t_2,distx[i],0,N)
        y3[i] = c_c0_S(t_3,distx[i],0,N)
    axsLeft[0].plot(distx,y1,c=clr1,marker="None")
    axsLeft[0].plot(distx,y2,c=clr2,marker="None")
    axsLeft[0].plot(distx,y3,c=clr3,marker="None")
    
    # Second plot with concentrations into the fracture at three
    # distances (along the fracture) for a single simulation time (t_3).

    dMn=res[(res["time"]==t_3) & (res["x"]==xmn)]
    dMd=res[(res["time"]==t_3) & (res["x"]==xmd)]
    dMx=res[(res["time"]==t_3) & (res["x"]==xmx)]
    clr11 = "#AA4499" # mauve
    clr12 = "#CC6677"  # orange
    clr13 = "#FFB000"  # yellow
    xLabel = 'Distance z in matrix (m)'
    axsLeft[1].set_ylim([-0.05, 1.05])
    axsLeft[1].plot(dMn['z'],dMn['A'],c=clr11,marker=mrk1,linestyle='None',markersize=9)
    axsLeft[1].plot(dMd['z'],dMd['A'],c=clr12,marker=mrk2,linestyle='None',markersize=6)
    axsLeft[1].plot(dMx['z'],dMx['A'],c=clr13,marker=mrk3,linestyle='None',markersize=9)
    axsLeft[1].set_ylabel(yLabel)
    axsLeft[1].set_xlabel(xLabel)
    axsLeft[1].set_title("{:.1f}".format(t_3)+" days", x=0.7, y=0.7)
    # calculate and plot the analytical (theoretical) curves
    distz = [0]*21
    y1 = [1]*21
    y2 = [1]*21
    y3 = [1]*21
    for i in range(0,21):
        distz[i]=(spacing/2.)*i/20
        y1[i] = c_c0_S(t_3,xmn,distz[i],N)
        y2[i] = c_c0_S(t_3,xmd,distz[i],N)
        y3[i] = c_c0_S(t_3,xmx,distz[i],N)
    axsLeft[1].plot(distz,y1,c=clr11,marker="None")
    axsLeft[1].plot(distz,y2,c=clr12,marker="None")
    axsLeft[1].plot(distz,y3,c=clr13,marker="None")

    # Subfigure to the right:
    # Three subplots (axes) in a single column.
    axsRight = subfigs[1].subplots(3, 1)

    # do not show the axes
    axsRight[0].axis('off')
    axsRight[1].axis('off')
    axsRight[2].axis('off')

    # First 'plot': an empty plot defining the legend texts
    axsRight[0].plot([],[],c=clr1,marker=mrk1,linestyle='None',markersize=9,label="{:.1f}".format(t_1)+" days")
    axsRight[0].plot([],[],c=clr2,marker=mrk2,linestyle='None',markersize=6,label="{:.1f}".format(t_2)+" d")
    axsRight[0].plot([],[],c=clr3,marker=mrk3,linestyle='None',markersize=9,label="{:.1f}".format(t_3)+" d")
    axsRight[0].legend(fontsize='small', loc=(0,1), borderaxespad=0)

    # Second 'plot': a text box listing parameters used.
    txt=annotTxt()
    axsRight[1].annotate(txt, fontsize=8,
             xy=(0.06,0.75), xycoords='subfigure fraction',  verticalalignment='top',
             bbox=dict(boxstyle="square,pad=0.5",linewidth=1,facecolor='white'))

    # Third 'plot': an empty plot defining the legend texts
    axsRight[2].plot([],[],c=clr11,marker=mrk1,linestyle='None',markersize=9,label="{:.2f}".format(Delta_x/2.)+" m")
    axsRight[2].plot([],[],c=clr12,marker=mrk2,linestyle='None',markersize=6,label="{:.2f}".format(xmd)+" m")
    axsRight[2].plot([],[],c=clr13,marker=mrk3,linestyle='None',markersize=9,label="{:.2f}".format(xmx)+" m")
    axsRight[2].legend(fontsize='small', loc=(0,1), borderaxespad=0) #, bbox_to_anchor=(0.,2.)

    # save the figure
    plotName = "Lipson.png"
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

# x-coordinate: get the minimum, maximum and mid distance
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

# time: get the closest time to 10, 30 and 100 days (for plotting)
time_1 = 10; time_2 = 30; time_3 = 100
print("--- finding results close to "+str(time_1)+", "+str(time_2)+" and "+str(time_3)+" days")
t_1 = -1; t_2 = -1; t_3 = -1
dist_1 = 1e20; dist_2 = 1e20; dist_3 = 1e20
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
print('    '+str(time_1)+' years = '+str(t_1)+'\n    '+str(time_2)+' yr = '+str(t_2)+'\n    '+str(time_3)+' yr = '+str(t_3))
if (t_1 <0) | (t_2 <0) | (t_3 < 0):
    print('Error, one of these values is negative...')
    raise SystemExit

n_times = times.shape[0]-1
print('    n_times='+str(n_times))

res0 = res[(res["time"]==t_1) & (res["x"]>=0) & (res["z"]==0)]
n_cells = res0.shape[0]-1
print('    n_cells='+str(n_cells))
Delta_x = L/n_cells
print('    n_cells='+str(n_cells)+'   Delta_x={:.3f}'.format(Delta_x))

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
plot(res,incl_x)
