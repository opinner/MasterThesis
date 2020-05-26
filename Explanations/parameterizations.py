import matplotlib.pyplot as plt
import numpy as np

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
    return 0.2
        
#vectorized version (acts on every element of an array)
def Osborn(Reb):
    """
    vectorized version of the Osborn(1980) turbulence parametrizaion
    """
    vOsborn = np.vectorize(basic_Osborn, otypes=[float])    
    return vOsborn(Reb)


f1,axis1 = plt.subplots()

 
Reb = np.linspace(0,1000,1000)

#bigger linewidths
linewidth = 4

axis1.plot(Reb,Osborn(Reb), linewidth=linewidth, label = "Osborn")
#axis1.plot(Reb,BB(Reb), linewidth=linewidth, label = "Bouffard & Boegman")
axis1.plot(Reb,Skif(Reb), "--", linewidth=linewidth ,label = "Shih et al")


axis1.set_ylim(-0.01,0.21) 

axis1.set_xlabel(r"$Re_b$")
axis1.set_ylabel(r"$\Gamma$")

axis1.legend(loc = "lower right")
axis1.set_title("Turbulence parametrizations")

plt.tight_layout
plt.savefig("./parametrization",dpi = 300)
plt.show()      
