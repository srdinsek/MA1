import math
import numpy as np
import matplotlib.pyplot as plt
from sympy import simplify
from random import random
from scipy.optimize import minimize
from scipy.optimize import leastsq
from scipy import linalg
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from numpy import linalg as LA
from scipy.linalg import solve
from timeit import default_timer as timer
from random import random
from secrets import randbelow
from random import SystemRandom
from scipy.stats import chisquare, norm, kstest, chi2, moment, ks_2samp
from scipy.special import erf
from mpl_toolkits.mplot3d import axes3d
from scipy.signal import gaussian








import seaborn as sns
sns.set_palette(sns.color_palette("autumn", 5))

#%% 
###########################################################################################################
#
#       1. NALOGA - KONSTANTNA GOSTOTA - (volumen, vztrajnostni)
#
###########################################################################################################

# to je pravilno!

def Konstantna(N):
    vsota=0
    vztrajnostni=0
    vztrajnostni2=0
    for i in range(N):
        x=1-2*random()
        y=1-2*random()
        z=1-2*random()
        if x*x+y*y <= 1 and y*y+z*z <= 1 and x*x+z*z <= 1:
            vsota+=1
            vztrajnostni+= x*x+y*y+z*z
            vztrajnostni2+= (x*x+y*y+z*z)*(x*x+y*y+z*z)
    return 'Masa znaša:',8*vsota/N,'+-',(8*vsota/N-8*(2-np.sqrt(2)))/(8*(2-np.sqrt(2))), 'in ocena napake:',8*np.sqrt(vsota/N-vsota*vsota/(N*N))/np.sqrt(N),'\t Vztrajnostni moment pa:',8*vztrajnostni/N,'in ocena napake:',8*np.sqrt(vztrajnostni2/N-vztrajnostni*vztrajnostni/(N*N))/np.sqrt(N)


def Natančnost(N):
    ana=[5**i for i in range(N)]
    vsota=[]
    iner=[]
    vsota2=[]
    for i in ana:
        vsota.append(Konstantna(i)[3])
        vsota2.append(Konstantna(i)[5])
        iner.append(Konstantna(i)[-1])
    plt.title('Natančnost mase in vztrajsnostnega momenta')
    plt.plot(ana,vsota,'r.', label='točna napaka za maso')
    plt.plot(ana,vsota2,'b.', label='ocena napake za maso')
    plt.plot(ana,iner,'g.', label='ocena napake za vztrajnostni moment')
    plt.xlabel('Število točk')
    plt.ylabel('Velikost napake')
    plt.yscale('log')
    plt.xscale('log')
    plt.legend()

#%%
###########################################################################################################
#
#       1. NALOGA - GOSTATA FUNKCIJA r - (volumen, vztrajnostni)
#
###########################################################################################################

# to je zaenkrat narobe!!

# ta spremeni porazdelitev in to sferično

def Potenca(N,p):
    masa=0
    volumen=0
    volumen2=0
    vztrajnostni=0
    vztrajnostni2=0
    koren=np.sqrt(3/2)
    for i in range(N):
        fi=2*np.pi*random()
        th=np.arccos(2*random()-1)
        r=koren*(random())**(1/(p+3))

        
        x=r*np.cos(fi)*np.sin(th)
        y=r*np.sin(fi)*np.sin(th)
        z=r*np.cos(th)
        l=x*x+y*y
        
        if x*x+y*y <= 1 and y*y+z*z <= 1 and x*x+z*z <= 1:
            volumen+=(koren/r)**p
            volumen2+=(koren/r)**(2*p)
            masa+=1
            vztrajnostni+=l
            vztrajnostni2+= l*l
    return 'Volumen znaša:',4*np.pi*koren*koren*koren/(p+3)*volumen/N, 'in ocena napake:',4*np.pi*koren*koren*koren/(p+3)*np.sqrt(volumen2/N-volumen*volumen/(N*N))/np.sqrt(N),'Masa znaša:',4*np.pi*koren*koren*koren/(p+3)*masa/N, 'in ocena napake:',4*np.pi*koren*koren*koren/(p+3)*np.sqrt(masa/N-masa*masa/(N*N))/np.sqrt(N),'\t Vztrajnostni moment pa:',4*np.pi*koren*koren*koren/(p+3)*vztrajnostni/N,'in ocena napake:',4*np.pi*koren*koren*koren/(p+3)*np.sqrt(vztrajnostni2/N-vztrajnostni*vztrajnostni/(N*N))/np.sqrt(N)


# ta le spremeni funkcijo porazdelitev ostane ista

def PotencaK(N,p):
    masa=0
    masa2=0
    volumen=0
    volumen2=0
    vztrajnostni=0
    vztrajnostni2=0
    koren=np.sqrt(3/2)
    for i in range(N):
        x=1-2*random()
        y=1-2*random()
        z=1-2*random()
        
        l=x*x+y*y
        
        if x*x+y*y <= 1 and y*y+z*z <= 1 and x*x+z*z <= 1:
            r=(x*x+y*y+z*z)
            volumen+=1
            volumen2+=1
            masa+=(np.sqrt(r)/koren)**p
            masa2+=(np.sqrt(r)/koren)**(2*p)
            vztrajnostni+=(l)*(np.sqrt(r)/koren)**p
            vztrajnostni2+= (l)*(l)*(np.sqrt(r)/koren)**(2*p)
    return 'Volumen znaša:',8*volumen/N, 'in ocena napake:',8*np.sqrt(volumen2/N-volumen*volumen/(N*N))/np.sqrt(N),'Masa znaša:',8*masa/N, 'in ocena napake:',8*np.sqrt(masa2/N-masa*masa/(N*N))/np.sqrt(N),'\t Vztrajnostni moment pa:',8*vztrajnostni/N,'in ocena napake:',8*np.sqrt(vztrajnostni2/N-vztrajnostni*vztrajnostni/(N*N))/np.sqrt(N)




# natančnost če spremenimo porazdelitev


def NatančnostP(N,p):
    ana=[5**i for i in range(N)]
    for j in range(p):
        volumen=[]
        masa=[]
        iner=[]
        for i in ana:
#            volumen.append(Potenca(i,j)[3])
            masa.append(Potenca(i,j)[7])
#            iner.append(Potenca(i,j)[-1])
        plt.title('Napaka pri različnih potencah p')
#        plt.plot(ana,volumen,'.', label='ocena napake za volumen p={}'.format(j))
        plt.plot(ana,masa,'o',label='ocena napake za maso p={}'.format(j))
#        plt.plot(ana,iner,'x', label='ocena napake za vztrajnostni moment p={}'.format(j))
    plt.xlabel('Število točk')
    plt.ylabel('Velikost napake')
    plt.yscale('log')
    plt.xscale('log')
    plt.legend()


# rišemo vrednosti, če spremenimo porazdelitev

def risanjeP(N,p):
    ana=[i for i in range(p)]
    volumen=[]
    masa=[]
    iner=[]
    for j in range(p):
        volumen.append(Potenca(N,j)[1])
        masa.append(Potenca(N,j)[5])
        iner.append(Potenca(N,j)[9])
    plt.title('Vrednosti pri različnih potencah p in N={}'.format(N))
    plt.plot(ana,volumen,'g:',label='volumen')
    plt.plot(ana,masa,'r:', label='masa')
    plt.plot(ana,iner,'b:', label='vztrajnostni moment')
    plt.ylabel('Vrednots volumna, mase oz. vztrajnostnega momenta')
    plt.xlabel('log(Vrednost p)')
    plt.yscale('log')
#    plt.xscale('log')
    plt.legend()


# rišemo porazdelitev če porazdelitev ostane ista

def risanjeK(N,p):
    ana=[i-p for i in range(2*p)]
    volumen=[]
    masa=[]
    iner=[]
    for j in range(2*p):
        volumen.append(PotencaK(N,j-p)[1])
        masa.append(PotencaK(N,j-p)[5])
        iner.append(PotencaK(N,j-p)[-3])
    plt.title('Vrednosti pri različnih potencah p in N={}'.format(N))
    plt.plot(ana,volumen,'g:',label='volumen')
    plt.plot(ana,masa,'r:', label='masa')
    plt.plot(ana,iner,'b:', label='vztrajnostni moment')
    plt.xlabel('Vrednots volumna, mase oz. vztrajnostnega momenta')
    plt.ylabel('log(Vrednost p)')
    plt.yscale('log')
#    plt.xscale('log')
    plt.legend()


#%%
###########################################################################################################
#
#       2. NALOGA - KROGLA KOT PORODNIŠNICA
#
###########################################################################################################


def normalno(N,mu):
    prehod=0
    for i in range(N):
        r=(random())**(1/3)
        th=np.arccos(2*random()-1)

        cos=np.cos(th)
        
        d=(-r*cos+np.sqrt(r*r*cos*cos-(r*r-1)))
        dtilde=-mu*np.log(1-random())
    
        if dtilde > d:
            prehod+=1
    return 'Ušlo je',prehod/N,'gama žarkov'

def nariši(N,mu):
    seznam=[normalno(N,i/20)[1] for i in range(mu)]
    t=[i/20 for i in range(mu)]
    plt.title('Verjetnost pobega v odvisnosti od razmerja $\mu$ N={}'.format(N))
    plt.plot(t,seznam,'r:')
    plt.xlabel('$\mu$')
    plt.ylabel('n/N')
    
#%%
#---------------------------------------------------------------------------------------------------------
#       2.NALOGA - odvisnost od radija
#---------------------------------------------------------------------------------------------------------    

def nono(N,mu,h):
    prehod=0
    koti=[]
    for i in range(N):
        r=(random())**(1/3)
        th=np.arccos(2*random()-1)

        cos=np.cos(th)
        
        d=(-r*cos+np.sqrt(r*r*cos*cos-(r*r-1)))
        dtilde=-mu*np.log(1-random())
        
        if r> h-0.05 and r < h+0.05:
            if dtilde > d:
                prehod+=1
                koti.append(th)
    plt.hist(koti,bins='fd')
    return 'Ušlo je',prehod/N,'gama žarkov'







#%%
###########################################################################################################
#
#       3. NALOGA - NEVTRONSKI REFLEKTOR
#
###########################################################################################################

#---------------------------------------------------------------------------------------------------------
#       tu ni napak!
#---------------------------------------------------------------------------------------------------------

def nevtron(N):
    kn=[]
    for i in range(N):
        n=0
        s=0
        while s >= 0 and s <= 1:
            krneki=random()
            s+=(1-2*krneki)/krneki*0.5*np.log(1-random())
            n+=1

        kn.append(n)
    
    plt.title('sipanje na plošči v ngativni smeri N={}'.format(N))
    plt.hist(kn,bins='fd')
    plt.xlabel('število korakov')
    plt.ylabel('število dogodkov pri teh korakih')
#    plt.yscale('log')

def nevtronL(N,mu):
    nn=0
    kp=0
    for i in range(N):
        s=0
        while s >= 0 and s <= 1:
            krneki=random()
            s+=(1-2*krneki)/krneki*mu*np.log(1-random())

        if s > 0:
            kp+=1
        if s < 0:
            nn+=1
            
    return kp/nn

def narisi(N,mu):
    totole=[nevtronL(N,i/10+0.1) for i in range(mu)]
    t=[i/10+0.1 for i in range(mu)]
    plt.title('Razmeje med sipanjem naprej in nazaj')
    plt.plot(t,totole,'r:')
    plt.xlabel('vrednost $\mu$')
    plt.ylabel('razmerje naprej/nazaj')


#%%
###########################################################################################################
#
#       3. NALOGA - NEVTRONSKI REFLEKTOR - IZOTROPNO
#
###########################################################################################################

def izotropno(N):
    kn=[]
    for i in range(N):
        n=0
        s=0
        while s >= 0 and s <= 1:
            krneki=random()
            s+=(1-2*krneki)/krneki*0.5*np.log(1-random())*(2*random()-1)
            n+=1

        kn.append(n)
    
    plt.title('sipanje na plošči v ngativni smeri N={}'.format(N))
    plt.hist(kn,bins='fd')
    plt.xlabel('število korakov')
    plt.ylabel('število dogodkov pri teh korakih')
#    plt.yscale('log')

def izotropnoL(N,mu):
    nn=0
    kp=0
    for i in range(N):
        s=0
        while s >= 0 and s <= 1:
            krneki=random()
            s+=(1-2*krneki)/krneki*mu*np.log(1-random())*(2*random()-1)

        if s > 0:
            kp+=1
        if s < 0:
            nn+=1
            
    return kp/nn

def izonarisi(N,mu):
    totole=[izotropnoL(N,i/10+0.1) for i in range(mu)]
    t=[i/10+0.1 for i in range(mu)]
    plt.title('Razmeje med sipanjem naprej in nazaj, če upoštevamo izotropnost')
    plt.plot(t,totole,'r:')
    plt.xlabel('vrednost $\mu$')
    plt.ylabel('razmerje naprej/nazaj')




