# -*- coding: utf-8 -*-
"""
Created on Thu May 11 20:39:17 2023

@author: Marc Mignard
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import femm #pip install pyfemm

femm.openfemm()
draw=False

##########################################################################
###    Utility functions to analyze manetorque rods with FEMM
###     
##########################################################################
mu0 = 1.257e-6 #permeability of vacuum, H/m or N/A^2

def addnode(x,y,draw=False):
    femm.mi_addnode(x,y) 
    if draw:
        plt.plot(x,y,'bo')
        
def addsegment(x1,y1,x2,y2,draw=False):
    femm.mi_addsegment(x1,y1,x2,y2)
    if draw:
        plt.plot([x1,x2],[y1,y2],'b')
        
def addblocklabel(x,y,matl,a1,a2,a3,a4,a5,a6):      
    femm.mi_addblocklabel(x,y)
    femm.mi_selectlabel(x,y)
    femm.mi_setblockprop(matl,a1,a2,a3,a4,a5,a6)
    femm.mi_clearselected()
    if draw:
        plt.plot(x,y,'go')
        
def makeSolenoid(cPts,wPts,draw=False):
    ''' 
    Create a solenoid using numpy arrays to define the points. Using cylindrically
    symmetric design, so only need to define half the solenoid (in the radial direction)
    With an antiperiodic boundary condition at the bottom in the axial direction the 
    the simulations would be twice as fast, but I never got around to doing that.
    
    cPts is the core boundary. For a cylindrical core, only need two points with 
        [[r,-h/2], [r,h/2]]
    because 0,-h/2 and 0,h/2 are assumed to be part of the core.
    
    wPts defines the winding boundary. Again only two points are needed with
        [[r+t,-h/2+bel], [r+t,h/2-bel]] 
    where 't' is the radial thickness of the winding, and 
    where 'bel' is the length of core at the end not covered by winding
    The top and bottom of the winding are assumed to be flat, and the winding
    extends to the left all the way to the core.
    '''
    
    #calculate total height (length) of the core
    h = np.max(cPts[:,1]) - np.min(cPts[:,1])
    
    #draw core
    addnode(0,-h/2,draw)        
    addnode(0,h/2,draw)
    for z in np.arange(0,cPts.shape[0]): #right edge of core
        addnode(cPts[z,0],cPts[z,1],draw)  
        if z>0:
            addsegment(cPts[z-1,0],cPts[z-1,1], \
                       cPts[z,0],cPts[z,1],draw)
    addsegment(0,-h/2,0,h/2,draw)  #center
    addsegment(0,-h/2,cPts[0,0],cPts[0,1],draw)  #bottom
    addsegment(0,h/2,cPts[-1,0],cPts[-1,1],draw) #top

    #draw winding
    for z in np.arange(0,wPts.shape[0]): #right edge of core
        addnode(wPts[z,0],wPts[z,1],draw)  
        if z>0:
            addsegment(wPts[z-1,0],wPts[z-1,1], \
                       wPts[z,0],wPts[z,1],draw)
    #bottom
    interpR = cPts[0,0]-(cPts[1,0]-cPts[0,0])/(cPts[1,1]-cPts[0,1])*(cPts[0,1]-wPts[0,1])
    if (np.abs(interpR-cPts[0,0])<1e-3) & (np.abs(wPts[0,1]-cPts[0,1])<1e-3):
        print('ignore interpolation')
        addsegment(cPts[0,0],cPts[0,1],wPts[0,0],wPts[0,1],draw)  
    else:
        addnode(interpR,wPts[0,1],draw) 
        addsegment(interpR,wPts[0,1],wPts[0,0],wPts[0,1],draw)
    #top
    n=cPts.shape[0]-1
    interpR = cPts[n,0]-(cPts[n,0]-cPts[n-1,0])/(cPts[n,1]-cPts[n-1,1])*(cPts[n,1]-wPts[-1,1])
    if (np.abs(interpR-cPts[n,0])<1e-3) & (np.abs(wPts[-1,1]-cPts[n,1])<1e-3):
        print('ignore interpolation')
        addsegment(cPts[n,0],cPts[n,1],wPts[-1,0],wPts[-1,1],draw)  #bottom
    else:
        addnode(interpR,wPts[-1,1],draw) 
        addsegment(interpR,wPts[-1,1],wPts[-1,0],wPts[-1,1],draw)  #bottom

def Bradial(R,L,M,mu0 = 1.257e-6):
    return mu0*M/(4*np.pi)/np.power(R**2+L**2/4,1.5)

def Baxial(R,L,M,mu0 = 1.257e-6):
    return mu0*M/(4*np.pi)*((R/L-1/2)/np.power(R**2-R*L+L**2/4,1.5) - \
                            (R/L-1/2)/np.power(R**2+R*L+L**2/4,1.5))
        
ErrFuncAxial=lambda M,p,fit: Baxial(p,h/100,M)-fit    
ErrFuncRadial=lambda M,p,fit: Bradial(p,h/100,M)-fit    
        
def analyzeSolenoid(cPts,wPts,coreMat,airMesh,wCurr,pos):
    '''
    Create the physical design of a solenoid, then calculate the magnetic flux
    in the axial and radial positions defined by numpy array 'pos'
    
    for definitions of cPts & wPts see function makeSolenoid() above
    coreMat is a magnetic material from the FEMM library
    airMesh controls the coarseness of the mesh, see documentation on femm.mi_setblockprop()
    wCurr is the current in the winding
    pos is a numpy array that defines where to measure the magnetic flux. This is confusing.
        I should have defined pos as a 2-d array, or as posR and posA for radial and axial directions.
        Instead, for the axial flux calculation, pos is assumed to be on a vertical line at radius=0.
        For the radial flux calculation, pos is assumed to be on a horizontal line at z=0.
        
    Usually I call analyzeSolenoid() many times to collect data with different currents.
    Could have made this more efficient by calling makeSolenoid() one time, and changing the
    winding current on that one design.
    '''
    femm.newdocument(0) #0 is magnetics problem
    femm.mi_probdef(0,'centimeters','axi',1e-8,0.5,30,0)
    #femm.mi_seteditmode('blocks') #don't seem to need this
    femm.mi_setgrid(1,'cart')
    femm.mi_addmaterial('Air',1,1,0,0,0,0,0,1,0,0,0)
    femm.mi_getmaterial(coreMat);
    femm.mi_getmaterial('26 AWG'); 
    femm.mi_addcircprop('winding',wCurr,1) #winding is 1000A-turns, and it is in series
    makeSolenoid(cPts,wPts,False)  
    ctr = int(cPts.shape[0]/2)
    addblocklabel(cPts[ctr,0]/2,0,coreMat,1,0,0,0,0,0)    
    addblocklabel((cPts[ctr,0]+wPts[ctr,0])/2,-2,'26 AWG',1,0,'winding',1,0,0)  
    if airMesh==0:
        addblocklabel(cPts[ctr,0]+2,-4,'Air',1,0,0,0,0,0)
    else:
        addblocklabel(cPts[ctr,0]+2,-4,'Air',0,airMesh,0,0,0,0) #airMesh=1 is fairly good. 0.5 is better but slower
        
    femm.mi_makeABC(7,120,0,0,0)
    
    femm.mi_saveas("temp.fem")
    femm.mi_createmesh()
    femm.mi_analyze()
    femm.mi_loadsolution()
    femm.mo_showdensityplot(1,0,1,0,'bmag')
    
    bz=np.zeros(pos.size)
    for i in np.arange(0,pos.size):
        bz[i] = femm.mo_getb(0,pos[i])[1]
    
    br=np.zeros(pos.size)
    for i in np.arange(0,pos.size):
        br[i]  = femm.mo_getb(pos[i],0)[1]
        
    femm.mi_close()
    return -br,bz
    
##########################################################################
###    Create and analyze a 300mm long magnetorquer rod with mu metal core
###     
##########################################################################

r=0.7
h=30
cPts = np.array([
    [r,-h/2],
    [r,h/2]])
wPts = np.array([
    [r+1,-h/2+1],
    [r+1,h/2-1]])

#core saturation
pos = np.linspace(0,h/2,50)
curr = np.array([3000,2000,1000,750,500])
bz = np.zeros([pos.size,curr.size])
for c in np.arange(0,curr.size):
    br,bz[:,c] = analyzeSolenoid(cPts,wPts,'Mu Metal',0,curr[c],pos)

plt.figure(figsize=(5,3.5),dpi=300)
plt.title('B axial in center of mu metal core')
for c in np.arange(0,curr.size):
    plt.plot(pos,bz[:,c],label='i={:.1f}kA-turn'.format(curr[c]/1000))
plt.annotate('Mu metal\nsaturates\nat 0.7T', xy=(6,0.69), xytext=(7,0.2), arrowprops=dict(facecolor='black', width=1, headwidth=5, headlength=8))

plt.xlabel('position (cm)')
plt.ylabel('B (T)')
plt.grid(True)
plt.xlim([pos[0],pos[-1]])
#plt.ylim([0,2])
plt.legend()
plt.savefig('BinCoreMuMetalSaturation.svg', bbox_inches='tight')
plt.show()

##########################################################################
###    Analyze a 300mm magnetorquer rod with Hiperco-50 core
###     
##########################################################################

r=0.7
h=30
cPts = np.array([
    [r,-h/2],
    [r,h/2]])
wPts = np.array([
    [r+1,-h/2+1],
    [r+1,h/2-1]])

#core saturation
pos = np.linspace(0,h/2,50)
curr = np.array([8000,5000,3000,2000,1000])
#curr = np.array([8000,5000,4000,3000,2000,1000,500,200])
bz = np.zeros([pos.size,curr.size])
for c in np.arange(0,curr.size):
    br,bz[:,c] = analyzeSolenoid(cPts,wPts,'Hiperco-50',0,curr[c],pos)
plt.figure(figsize=(5,3.5),dpi=300)
plt.title('B axial in center of Hiperco-50 core')
for c in np.arange(0,curr.size):
    plt.plot(pos,bz[:,c],label='i={:.1f}kA-turn'.format(curr[c]/1000))
plt.annotate('Hiperco-50\nsaturates\nat 2.3T', xy=(4,2.3), xytext=(2,1.1), 
             arrowprops=dict(facecolor='black', width=1, headwidth=5, headlength=8))
plt.xlabel('position (cm)')
plt.ylabel('B (T)')
plt.grid(True)
plt.xlim([pos[0],pos[-1]])
#plt.ylim([0,2])
plt.legend()
plt.savefig('BinCoreHP50Saturation.svg', bbox_inches='tight')
plt.show()

##########################################################################
###    Flux comparison to calculated values
###     
##########################################################################

#estimate magnetic moment
pos=np.linspace(50,100,101)
br,bz = analyzeSolenoid(cPts,wPts,'Mu Metal',0,1000,pos)

#In the axial direction, the simulation doesn't give accurate values until a little further
#from the end of the solenoid core, so begin the fitting later. To get the simulation to
#give more accurate values in the near field requires a finer mesh which takes longer to simulate.
#Estimate the magnet moment that best fits the simulation data
beg = 10; MeAxial,success = sp.optimize.leastsq(ErrFuncAxial,10,args=(pos[beg:-1]/100,bz[beg:-1]))
beg = 0; MeRadial,success = sp.optimize.leastsq(ErrFuncRadial,10,args=(pos[beg:-1]/100,br[beg:-1]))

plt.figure(figsize=(5,3.5),dpi=300)
plt.title('flux in axial and radial directions, Mu Metal core')
plt.plot(pos,bz*1e6,label='sim z')
plt.plot(pos,Baxial(pos/100,h/100,MeAxial)*1e6,':',label='CASJ eq z, M={:.1f}'.format(MeAxial[0]))
plt.plot(pos,br*1e6,label='sim r')
plt.plot(pos,Bradial(pos/100,h/100,MeRadial)*1e6,':',label='CASJ eq r, M={:.1f}'.format(np.abs(MeRadial[0])))
plt.xlabel('position (cm)')
plt.ylabel('B (uT)')
plt.grid(True)
plt.xlim([pos[0],pos[-1]])
#plt.ylim([0,2])
plt.legend()
plt.savefig('MuMetalFlux.svg', bbox_inches='tight')
plt.show()

##########################################################################
###   Comparison of magnetic moments with different materials
###     
##########################################################################

#magnetic moment as a function of current
pos=np.linspace(50,100,101)
curr = np.array([8000,5000,4000,3000,2000,1000,500,200])
bz = np.zeros([pos.size,curr.size])
br = np.zeros([pos.size,curr.size])
MeAxialMu = np.zeros(curr.size)
MeRadialMu = np.zeros(curr.size)
MeAxialVP = np.zeros(curr.size)
MeRadialVP = np.zeros(curr.size)
for c in np.arange(0,curr.size):
    br[:,c],bz[:,c] = analyzeSolenoid(cPts,wPts,'Mu Metal',0,curr[c],pos)
    beg = 10; MeAxialMu[c],success = sp.optimize.leastsq(ErrFuncAxial,10,args=(pos[beg:-1]/100,bz[beg:-1,c]))
    beg = 0; MeRadialMu[c],success = sp.optimize.leastsq(ErrFuncRadial,10,args=(pos[beg:-1]/100,br[beg:-1,c]))

    br[:,c],bz[:,c] = analyzeSolenoid(cPts,wPts,'Hiperco-50',0,curr[c],pos)
    beg = 10; MeAxialVP[c],success = sp.optimize.leastsq(ErrFuncAxial,10,args=(pos[beg:-1]/100,bz[beg:-1,c]))
    beg = 0; MeRadialVP[c],success = sp.optimize.leastsq(ErrFuncRadial,10,args=(pos[beg:-1]/100,br[beg:-1,c]))

plt.figure(figsize=(5,3.5),dpi=300)
plt.title('Comparison of magnetic moments\n(Mu metal & Hiperco-50')
plt.plot(curr/1000,MeRadialMu,'ro-',markersize=3,label='Mu radial')
plt.plot(curr/1000,MeAxialMu,'ro--',markersize=3,label='Mu axial')
plt.plot(curr/1000,MeRadialVP,'bo-',markersize=3,label='H50 radial')
plt.plot(curr/1000,MeAxialVP,'bo--',markersize=3,label='H50 axial')
plt.xlabel('current (kAˑturns)')
plt.ylabel('magnetic moment (Am^2)')
plt.grid(True)
plt.xlim([0,curr[0]/1000])
plt.ylim([0,100])
plt.legend()
plt.savefig('MomentComparison.svg', bbox_inches='tight')
plt.show()


