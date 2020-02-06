import matplotlib.pyplot as plt
import numpy as np


def BB(Reb):
    Pr = 700

    if Reb < 0.18:
        return 0
        
    elif ((Reb >= 0.18) and (Reb <= 96.56)):
        return 0.1 * Pr**(-0.5) * Reb**(0.5)

    else:
        return 2*Reb**(-0.5)

def Skif(Reb):
    if Reb < 7:
        return 0
        
    elif ((Reb >= 7) and (Reb<=100)):
        return 0.2
        
    else:
        return 2*Reb**(-0.5)
        
def Osborn(Reb):
    return 0.2
        
vBB = np.vectorize(BB,otypes=[float])     
vSkif = np.vectorize(Skif, otypes=[float]) 
vOsborn = np.vectorize(Osborn, otypes=[float]) 

 
Reb = np.linspace(0,1000,1000)

linewidth = 4

plt.plot(Reb,vOsborn(Reb), linewidth=linewidth, label = "Osborn")
plt.plot(Reb,vBB(Reb), linewidth=linewidth, label = "BB")
plt.plot(Reb,vSkif(Reb), "--", linewidth=linewidth ,label = "Skif")


plt.legend()
plt.show()     
