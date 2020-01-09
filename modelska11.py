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
from winsound import Beep
from scipy.stats import poisson





import seaborn as sns
sns.set_palette(sns.color_palette("autumn", 8))


#%%
###########################################################################################################
#
#       1. NALOGA - ENOSTAVNO EKSPONENTNO
#
###########################################################################################################

def umirajoča(n,beta,N):
    rešitev=[]
    for i in range(n):
        if N>0:
            rešitev.append(N)
            N=N-np.random.poisson(N*beta)
        else:
            rešitev.append(0)
    return [i for i in range(n)],rešitev

def rojevajoča(n,beta,N):
    rešitev=[]
    for i in range(n):
        if N>0:
            rešitev.append(N)
            N=N+np.random.poisson(N*beta)
        else:
            rešitev.append(0)
    return [i for i in range(n)],rešitev

def rojstvosmrt(n,betaR,betaS,t,N):
    rešitev=[]
    for i in range(n):
        if N>0:
            rešitev.append(N)
            N=N+np.random.poisson(N*betaR*t)-np.random.poisson(N*t*betaS)
        else:
            rešitev.append(0)
    return [i for i in range(n)],rešitev



#%%
#RIŠEM VEČ PRIMEROV UMIRAJOČIH
for i in range(6):
    reši=umirajoča(100,0.05*i+0.01,250)
    plt.plot(reši[0],reši[1],label='razmerje beta*dt t={0:.2f}'.format(0.05*i+0.01))
plt.title('Začetna populacija = 25. Model umirajočih')
plt.xlabel('Število časovnih korakov')
plt.ylabel('Število osebkov v populaciji N')
plt.legend(loc=0)

#%%
#RIŠEM VEČ ISTIH UMIRAJOČIH
for i in range(10):
    reši=umirajoča(600,0.01,250)
    plt.plot(reši[0],reši[1])
plt.title('Začetna populacija = 250. Model umirajočih beta=0.01')
plt.xlabel('Število časovnih korakov')
plt.ylabel('Število osebkov v populaciji N')
plt.legend(loc=0)

#%%
#RIŠEM VEČ PRIMEROV ROJEVAJOČIH
for i in range(3):
    reši=rojevajoča(100,0.05*i+0.01,25)
    plt.plot(reši[0],reši[1],label='razmerje beta*dt t={0:.2f}'.format(0.05*i+0.01))
plt.title('Začetna populacija = 25. Model rojevajočih')
plt.xlabel('Število časovnih korakov')
plt.ylabel('Število osebkov v populaciji N')
plt.legend(loc=0)

#%%
#RIŠEM VEČ PRIMEROV ROJEVAJOČIH IN UMIRAJOČIH
for i in range(6):
    reši=rojstvosmrt(600,4,5,0.05*i+0.01,250)
    plt.plot(reši[0],reši[1],label='razmerje beta*dt t={0:.2f}'.format(0.05*i+0.01))
plt.title('Začetna populacija = 250. Model rojevajočih in umirajočih')
plt.xlabel('Število časovnih korakov')
plt.ylabel('Število osebkov v populaciji N')
plt.legend(loc=0)

#%%
#RIŠEM VEČ ISTIH ROJEVAJOČIH IN UMIRAJOČIH
for i in range(10):
    reši=rojstvosmrt(600,4,5,0.01,250)
    plt.plot(reši[0],reši[1])
plt.title('Začetna populacija = 250. Model rojevajočih in umirajočih beta=0.01')
plt.xlabel('Število časovnih korakov')
plt.ylabel('Število osebkov v populaciji N')
plt.legend(loc=0)

#%%
###########################################################################################################
#
#       1. NALOGA - STATISTIKA ČASOV IZUMRTJA - ENOJNA
#
###########################################################################################################

def stat_umirajoča(n,beta,t,N):
    rešitev=[]
    for j in range(n):
        i=0
        število=N
        while število > 0:
            število=število-np.random.poisson(število*beta*t)
            i+=1
        rešitev.append(i*t)
    return [j for j in range(n)],rešitev

#%%
#RIŠEM HISTOGRAM ZA ČASE IZUMRTJA - UMIRAJOČI MODEL

for i in range(3):
    rešitev=stat_umirajoča(5000,(0.1*i+0.1),0.01,25)
    plt.hist(rešitev[1],bins=100,alpha=0.75,label='beta={0:.2f}'.format((0.1*i+0.1)))
plt.title('Statistika časov izumrtja N = 25. Model umirajočih, dt=0.01')
plt.xlabel('Čas izumrtja')
plt.ylabel('Število dogodkov')
plt.legend(loc=0)

#%%
#POVPREČNI ČAS IZUMRTJA V ODVISNOSTI OD BETA - UMIRAJOČI MODEL
povprečje=[]
N=20
for k in range(4):
    povprečje=[]
    for i in range(N):
        rešitev=stat_umirajoča(5000,0.05*i+0.1,(0.1*k+0.1)*3,250)
        povprečje.append(np.mean(rešitev[1]))
    plt.title('Povprečni čas izumrtja N = 250. Model umirajočih')
    plt.plot([i for i in range(N)],povprečje,'.',label='dt={0:.2f}'.format((0.1*k+0.1)*3))
    plt.plot([i for i in range(N)],povprečje,':')
plt.xlabel('vrednost beta')
plt.ylabel('Povprečni čas izumrtja')
plt.legend(loc=0)

#%%
#POVPREČNI ČAS IZUMRTJA V ODVISNOSTI OD dt - UMIRAJOČI MODEL
povprečje=[]
N=200
for i in range(4):
    povprečje=[]
    for k in range(N):
        rešitev=stat_umirajoča(5000,0.2*i+0.4,0.01*k+0.01,250)
        povprečje.append(np.mean(rešitev[1]))
    plt.title('Povprečni čas izumrtja N = 250. Model umirajočih')
    plt.plot([0.01*i+0.01 for i in range(N)],povprečje,'.',label='beta={0:.2f}'.format(0.2*i+0.4))
    plt.plot([0.01*i+0.01 for i in range(N)],povprečje,':')
plt.xlabel('vrednost dt')
plt.ylabel('Povprečni čas izumrtja')
plt.legend(loc=0)

#%%
###########################################################################################################
#
#       1. NALOGA - STATISTIKA ČASOV IZUMRTJA - DVOJNA
#
###########################################################################################################

def stat_umirajoča(n,beta,t,N):
    rešitev=[]
    for j in range(n):
        i=0
        število=N
        while število > 0:
            število=število-np.random.poisson(število*beta*t/2)-np.random.poisson(število*beta*t/2)
            i+=1
        rešitev.append(i*t)
    return [j for j in range(n)],rešitev

#%%
#RIŠEM HISTOGRAM ZA ČASE IZUMRTJA - UMIRAJOČI MODEL - DVOJNA

for i in range(3):
    rešitev=stat_umirajoča(5000,(0.1*i+0.1)*3,0.01,250)
    plt.hist(rešitev[1],bins=100,alpha=0.75,label='beta={0:.2f}'.format((0.1*i+0.1)*3))
plt.title('Statistika časov izumrtja N = 250. Model umirajočih, dt=0.01')
plt.xlabel('Čas izumrtja')
plt.ylabel('Število dogodkov')
plt.legend(loc=0)

#%%
###########################################################################################################
#
#       1. NALOGA - STATISTIKA ČASOV IZUMRTJA - ROJSTVO SMRT
#
###########################################################################################################

def stat_rojstvosmrt(n,betaR,betaS,t,N):
    rešitev=[]
    for j in range(n):
        i=0
        število=N
        while število > 0:
            število=število-np.random.poisson(število*betaS*t)+np.random.poisson(število*betaR*t)
            i+=1
        rešitev.append(i*t)
    return [j for j in range(n)],rešitev

#%%
#RIŠEM HISTOGRAM ZA ČASE IZUMRTJA - UMIRAJOČI MODEL

for i in range(3):
    rešitev=stat_rojstvosmrt(5000,4,5,0.01*i+0.01,250)
    plt.hist(rešitev[1],bins=100,alpha=0.5,label='beta={0:.2f}'.format(0.01*i+0.01))
plt.title('Statistika časov izumrtja N = 250. Model rojstvo-smrt, dt=0.01')
plt.xlabel('Čas izumrtja')
plt.ylabel('Število dogodkov')
plt.legend(loc=0)


#%%
#POVPREČNI ČAS IZUMRTJA V ODVISNOSTI OD dt - UMIRAJOČI MODEL
povprečje=[]
N=200
for i in range(N):
    rešitev=stat_rojstvosmrt(5000,4,5,0.01*i+0.01,25)
    povprečje.append(np.mean(rešitev[1]))
plt.title('Povprečni čas izumrtja N = 25. Model umirajočih')
plt.plot([0.01*i+0.01 for i in range(N)],povprečje,'.')
plt.plot([0.01*i+0.01 for i in range(N)],povprečje,':')
plt.xlabel('vrednost dt')
plt.ylabel('Povprečni čas izumrtja')
#plt.legend(loc=0)


#%%
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
#
#       2. NALOGA - MATRIKA PREHODOV ZA ZGORNJI MODEL
#
###########################################################################################################
###########################################################################################################
#      konstantni imirljivosti in smrtnosti
###########################################################################################################

def Matrika_prehodov(N,lamb,mu,dt):
    if (lamb+mu)*dt<1:
        return [[1-(lamb+mu)*dt if i==j else lamb*dt if i-1==j else mu*dt if i+1==j else 0 for j in range(N)] for i in range(N)]
    else:
        return False
def vektor(N):
    return [0 if i<N-1 else 1 for i in range(N)]

def vektor_n(N,lamb,mu,dt,n):
    b=np.array(vektor(N))
    A=np.array(Matrika_prehodov(N,lamb,mu,dt))
    for i in range(n):
        b=A.dot(b)
        vsota=sum(b)
        b=[b[i]/vsota for i in range(N)]
    return b
        
 
#zgled
#1-(lamb+mu)*dt         mu*dt              0            0
#   lamb*dt       1-(lamb+mu)*dt         mu*dt          0
#      0               lamb*dt      1-(lamb+mu)*dt    mu*dt
    
# mu so smrti
# lamb so rojstva
 
#%%
N=100
stanje=[i for i in range(N)]
for i in range(20):
        plt.plot(stanje,vektor_n(N,4,5,0.1,100*i+10),'.')
plt.title('Verjetnost za  pri dt=0.1, koraki 100*n+10')
plt.xlabel('populacijo s toliko osebki')
plt.ylabel('verjetnost za dano populacijo')

#zgled kako izgleda naša matrika za N=3, 4,5, dt=0.01
#[0.91, 0.05,  0  ]
#[0.04, 0.91, 0.05]
#[  0 , 0.04, 0.91]

#%%
###########################################################################################################
#      linearno odvisni imirljivosti in smrtnosti
###########################################################################################################

def Matrika_prehodovL(N,lamb,mu,dt):
    if (lamb+mu)*N*dt<1:
        return [[1-(lamb+mu)*(j)*dt if i==j else lamb*dt*(j) if i-1==j else mu*dt*(j) if i+1==j else 0 for j in range(N)] for i in range(N)]
    else:
        return False
    
def vektorL(N,š):
    return [0 if i!=š-1 else 1 for i in range(N)]

def vektor_nL(N,lamb,mu,dt,n,š):
    b=np.array(vektorL(N,š))
    A=np.array(Matrika_prehodovL(N,lamb,mu,dt))
    for i in range(n):
        b=A.dot(b)
        vsota=sum(b)
        b=[b[i]/vsota for i in range(N)]
    return b
        
 
#zgled
#1-(lamb+mu)*n*dt        mu*dt*n             0              0
#   lamb*dt*n       1-(lamb+mu)*dt*n      mu*dt*n           0
#      0               lamb*dt*n      1-(lamb+mu)*dt*n    mu*dt*n
    
# mu so smrti
# lamb so rojstva
 
#%%
    

import seaborn as sns
sns.set_palette(sns.color_palette("autumn", 10))


N=150
stanje=[i for i in range(N)]
for i in range(10):
#        plt.plot(stanje,vektor_nL(N,5,4,0.0001,500*i**2+100,100),',')
        plt.plot(stanje,vektor_nL(N,0,5,0.0001,150*i**2+100,100),'')
plt.plot([100,100],[0,0.15],'b:')
plt.title('Verjetnost pri samo umirajoči pop. dt=0.0001, koraki 50*n^2+100')
plt.xlabel('populacijo s toliko osebki')
plt.ylabel('verjetnost za dano populacijo')

#%%
###########################################################################################################
#      povprecna vrednost raznih porazdelitev
###########################################################################################################

import seaborn as sns
sns.set_palette(sns.color_palette("autumn", 10))

def vektor_nR(N,lamb,mu,dt,n,š):
    b=np.array(vektorL(N,š))
    A=np.array(Matrika_prehodovL(N,lamb,mu,dt))
    seznam=[]
    števec=[]
    for i in range(n):
        b=A.dot(b)
        vsota=sum(b)
        b=[b[i]/vsota for i in range(N)]
        if i%100:
            povpr=sum(b[i]*i for i in range(N))
            seznam.append(povpr)
            števec.append(i)
    return števec,seznam

N=150
n=10000

povprečje=[]



stanje=[i*0.0001 for i in range(n)]
resi=vektor_nR(N,4,10,0.0001,n,100)
plt.plot(resi[0],resi[1],'r-',label='umirajoča pop.= 10/4')
resi=vektor_nR(300,5,5,0.0001,n,100)
plt.plot(resi[0],resi[1],'g-',label='nevtralna pop.= 5/5')
resi=vektor_nR(500,10,4,0.0001,n,100)
plt.plot(resi[0],resi[1],'b-',label='rojevajoča pop.= 4/10')

plt.title('Povprečne vrednosti za razne primere dt=0.0001')
plt.xlabel('čas')
plt.ylabel('Povprečno število osebkov')
plt.legend()

#%%
###########################################################################################################
#      varianca raznih porazdelitev
###########################################################################################################

import seaborn as sns
sns.set_palette(sns.color_palette("autumn", 10))

def vektor_nV(N,lamb,mu,dt,n,š):
    b=np.array(vektorL(N,š))
    A=np.array(Matrika_prehodovL(N,lamb,mu,dt))
    seznam=[]
    števec=[]
    for i in range(n):
        b=A.dot(b)
        vsota=sum(b)
        b=[b[i]/vsota for i in range(N)]
        if i%100:
            povpr=sum(b[i]*i for i in range(N))
            var=sum((i-povpr)*(i-povpr)*b[i] for i in range(N))
            seznam.append(var)
            števec.append(i)
    return števec,seznam

N=150
n=10000

povprečje=[]



stanje=[i*0.0001 for i in range(n)]
resi=vektor_nV(N,4,10,0.0001,n,100)
plt.plot(resi[0],resi[1],'r-',label='umirajoča pop.= 10/4')
resi=vektor_nV(N,5,5,0.0001,n,100)
plt.plot(resi[0],resi[1],'g-',label='nevtralna pop.= 5/5')
resi=vektor_nV(N,10,4,0.0001,n,100)
plt.plot(resi[0],resi[1],'b-',label='rojevajoča pop.= 4/10')

plt.title('Varianca za razne primere dt=0.0001')
plt.xlabel('čas')
plt.ylabel('Varianca')
plt.legend()


#%%
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
#
#       2. NALOGA - MATRIKA PREHODOV ZA ZGORNJI MODEL
#
###########################################################################################################
###########################################################################################################
#      začetni model 5/4
###########################################################################################################




def zajci_lisica(n,a,b,Zo,Lo,dt):
    zajci=[]
    lisice=[]
    Z=Zo
    L=Lo
    for i in range(n):
        if Z>0 and L>0:
            zajci.append(Z)
            lisice.append(L)
            Z=Z+np.random.poisson(5*a*Z*dt)-np.random.poisson(4*a*Z*dt)-np.random.poisson(a/Lo*Z*L*dt)
            L=L+np.random.poisson(4*b*L*dt)-np.random.poisson(5*b*L*dt)+np.random.poisson(b/Zo*Z*L*dt)
        else:
            zajci.append(Z)
            lisice.append(L)
    return [i*dt for i in range(n)],zajci,lisice


#%%
#RIŠEM FAZNI DIAGRAM
    
import seaborn as sns
sns.set_palette(sns.color_palette("autumn", 2))


for i in range(1):
    reši=zajci_lisica(50000,1,2,200,50,0.01)
    
    plt.figure(0)
    plt.plot(reši[1],reši[2],'r')#label='razmerje beta*dt t={0:.2f}'.format(0.05*i+0.01))
    plt.title('Fazni diagram zajci-lisice')
    plt.xlabel('Število zajcev')
    plt.ylabel('Število lisic')
    
    plt.figure(1)
    plt.plot(reši[0],reši[1],'g:',label='Zajci')
    plt.plot(reši[0],reši[2],'r:',label='Lisice')
    plt.title('Odvisnost populacij od časa')
    plt.xlabel('Čas')
    plt.ylabel('Število osebkov')
    plt.legend(loc=0)

#%%
#ALGORITEM ZA RISANJE LEPIH

def zajci_lisicaHUD(n,a,b,Zo,Lo,dt):
    zajci=[]
    lisice=[]
    Z=Zo
    L=Lo
    for i in range(n):
        if Z>0 and L>0:
            zajci.append(Z)
            lisice.append(L)
            Zp=Z
            Z=Z+np.random.poisson(5*a*Z*dt)-np.random.poisson(4*a*Z*dt)-np.random.poisson(a/Lo*Z*L*dt)
            L=L+np.random.poisson(4*b*L*dt)-np.random.poisson(5*b*L*dt)+np.random.poisson(b/Zo*Zp*L*dt)
            j=i
        else:
            zajci.append(Z)
            lisice.append(L)
    return [i*dt for i in range(n)],zajci,lisice,j*dt
    
    
import seaborn as sns
sns.set_palette(sns.color_palette("autumn", 2))

reši=zajci_lisicaHUD(10000,1,2,200,50,0.01)
while reši[3]<60:
    reši=zajci_lisicaHUD(10000,1,2,200,50,0.01)

#%%
plt.figure(0)
plt.plot(reši[1],reši[2],'r')#label='razmerje beta*dt t={0:.2f}'.format(0.05*i+0.01))
plt.title('Fazni diagram zajci-lisice')
plt.xlabel('Število zajcev')
plt.ylabel('Število lisic')
    
plt.figure(1)
plt.plot(reši[0],reši[1],'g:',label='Zajci')
plt.plot(reši[0],reši[2],'r:',label='Lisice')
plt.title('Odvisnost populacij od časa')
plt.xlabel('Čas')
plt.ylabel('Število osebkov')
plt.legend(loc=0)


zajc=np.abs(np.fft.fft([reši[1][i] for i in range(5000)]))
lisc=np.abs(np.fft.fft([reši[2][i] for i in range(5000)]))
plt.figure(2)
plt.title('Fourijejeva transformacija za preimer zajcev')
plt.plot([i for i in range(len(zajc))],zajc)
plt.xlabel('Valovno število')
plt.ylabel('Amplituda')
plt.figure(3)
plt.title('Fourijejeva transformacija za preimer lisic')
plt.plot([i for i in range(len(lisc))],lisc)
plt.xlabel('Valovno število')
plt.ylabel('Amplituda')

#%%

vrhovi=[2.1,5.4,9.3,14.29,18,20.9,24.64,29.3,32.9,37,39.5,43.6,48.7,52.6,57.2]

razmaki=[(vrhovi[i+1]-vrhovi[i]) for i in range(len(vrhovi)-1)]

povprecen_razmak=sum(razmaki)/len(razmaki)
print(povprecen_razmak)

#%%
###########################################################################################################
#
#      3. NALOGA povprečna življenska doba običajne populacije 
#
###########################################################################################################

def zajci_lisicaHUD2(n,a,b,Zo,Lo,dt):
    Z=Zo
    L=Lo
    i=0
    while i<n and Z>0 and L>0:
        Zp=Z
        Z=Z+np.random.poisson(5*a*Z*dt)-np.random.poisson(4*a*Z*dt)-np.random.poisson(a/Lo*Z*L*dt)
        L=L+np.random.poisson(4*b*L*dt)-np.random.poisson(5*b*L*dt)+np.random.poisson(b/Zo*Zp*L*dt)
        j=i
        i=i+1
    return j*dt

#%%
N=10**4
sistema=[]
for i in range(N):
    d=zajci_lisicaHUD2(10000,1,2,200,50,0.01)
    sistema.append(d)
#%%

plt.title('Statistika življenskih dob')
plt.hist(sistema,150,alpha=0.75)
plt.xlabel('Življenska doba')
plt.ylabel('Število dogodkov')


print(sum(sistema)/len(sistema))



#%%

plt.figure(0)
plt.plot(b,c,'r')#label='razmerje beta*dt t={0:.2f}'.format(0.05*i+0.01))
plt.title('Fazni diagram zajci-lisice')
plt.xlabel('Število zajcev')
plt.ylabel('Število lisic')
    
plt.figure(1)
plt.plot(a,b,'g:',label='Zajci')
plt.plot(a,c,'r:',label='Lisice')
plt.title('Odvisnost populacij od časa')
plt.xlabel('Čas')
plt.ylabel('Število osebkov')
plt.legend(loc=0)


povb=sum(b[i] for i in range(len(b)))/len(b)
povc=sum(c[i] for i in range(len(c)))/len(c)
zajc=np.abs(np.fft.fft([b[i]-povb for i in range(5000)]))
lisc=np.abs(np.fft.fft([c[i]-povc for i in range(5000)]))
plt.figure(2)
plt.title('Fourijejeva transformacija za preimer zajcev')
plt.plot([i for i in range(len(zajc))],zajc)
plt.xlabel('Valovno število')
plt.ylabel('Amplituda')
plt.figure(3)
plt.title('Fourijejeva transformacija za preimer lisic')
plt.plot([i for i in range(len(lisc))],lisc)
plt.xlabel('Valovno število')
plt.ylabel('Amplituda')






