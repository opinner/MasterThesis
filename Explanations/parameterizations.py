import matplotlib.pyplot as plt
import numpy as np
#matplotlib preferences:
MINI_SIZE = 9
SMALL_SIZE = 10.95
MEDIUM_SIZE = 12
BIGGER_SIZE = 12
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE, titleweight = "bold")     # fontsize of the axes title
plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MINI_SIZE)    # legend fontsize
plt.rc('figure', titlesize=MEDIUM_SIZE, titleweight = "bold")  # fontsize of the figure title

#basic version
def basic_BB(Reb):
    Pr = 24.14

    if Reb < 0.18:
        return 0
        
    elif ((Reb >= 0.18) and (Reb < 96.56)):
        return 0.1 * Pr**(-0.5) * Reb**(0.5)

    elif ((Reb >= 96.56) and (Reb < 100)):
        return 0.2
    else:
        return 2*Reb**(-0.5)

#vectorized version (acts on every element of an array)
def BB(Reb):
    """
    vectorized version of the Osborn(1980) turbulence parametrizaion
    """
    vBB = np.vectorize(basic_BB,otypes=[float])
    return vBB(Reb)



def basic_Skif(Reb):
    if Reb < 7:
        return 0
        
    elif ((Reb >= 7) and (Reb<=100)):
        return 0.2
        
    else:
        return 2*Reb**(-0.5)

#vectorized version (acts on every element of an array)
def Skif(Reb):
    """
    vectorized version of the Shih(2005) turbulence parametrizaion
    """

    vSkif = np.vectorize(basic_Skif, otypes=[float]) 
    return vSkif(Reb)
     
        
def basic_Osborn(Reb):
    if Reb < 7:
        return 0
    else:
        return 0.2
        
#vectorized version (acts on every element of an array)
def Osborn(Reb):
    """
    vectorized version of the Osborn(1980) turbulence parametrizaion
    """
    vOsborn = np.vectorize(basic_Osborn, otypes=[float])    
    return vOsborn(Reb)


width = 6.2012
height = width / 1.618

#beamer figure sizes
#width = 1.7*4.252 #6.2012
#height = 1.5*3.7341 #* 4/3 #1.618

f1,axis1 = plt.subplots()

 
Reb = np.linspace(0,650,1000)

#bigger linewidths
linewidth = 4

axis1.plot(Reb,Osborn(Reb), linewidth=linewidth, label = "Osborn, 1980")
#axis1.plot(Reb,BB(Reb), linewidth=linewidth, label = "Bouffard & Boegman")
axis1.plot(Reb,Skif(Reb), "--", linewidth=linewidth ,label = "Shih et al., 2005")


axis1.set_ylim(-0.01,0.21) 
#axis1.set_xlim(-27,620)

axis1.set_xlabel(r"buoyancy Reynolds number $Re_b$")
axis1.set_ylabel(r"flux coefficient $\Gamma$")

axis1.legend(loc = "lower right")
axis1.set_title("Parametrizations of the flux coefficient ")

f1.set_size_inches(width,height)
f1.subplots_adjust(top=0.923,bottom=0.136,left=0.142,hspace=0.058,wspace=0.185)
f1.savefig("./parametrization",dpi = 600)
plt.show()      
