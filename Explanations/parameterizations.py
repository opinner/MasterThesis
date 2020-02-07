import matplotlib.pyplot as plt
import numpy as np

#basic version
def basic_BB(Reb):
    Pr = 700

    if Reb < 0.18:
        return 0
        
    elif ((Reb >= 0.18) and (Reb <= 96.56)):
        return 0.1 * Pr**(-0.5) * Reb**(0.5)

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


 
Reb = np.linspace(0,1000,1000)

#bigger linewidths
linewidth = 4

plt.plot(Reb,Osborn(Reb), linewidth=linewidth, label = "Osborn")
plt.plot(Reb,BB(Reb), linewidth=linewidth, label = "BB")
plt.plot(Reb,Skif(Reb), "--", linewidth=linewidth ,label = "Skif")


plt.legend()
plt.show()      
