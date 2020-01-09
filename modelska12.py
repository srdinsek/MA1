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
import matplotlib.image as mpimg
from PIL import Image, ImageFilter




import seaborn as sns
sns.set_palette(sns.color_palette("autumn", 4))

#%%
###########################################################################################################
#
#       1. NALOGA - ZAPIŠEMO FOURIJEJEVO TRANSFORMACIJO
#
###########################################################################################################




def fourie(signal):
    d=len(signal)
    if d%2 == 0:
        return 'ni sodo število točk'
    else:
        return [sum(signal[i]*np.exp(2j*np.pi*k*i/d) for i in range(d)) for k in range(d)]

start=timer()    
print(np.abs(fourie([1,2,3,4,5,2,3,4,2,2,6,7,2])))
end=timer()
print("\n",end-start)

start=timer() 
print(np.abs(np.fft.fft([1,2,3,4,5,2,3,4,2,2,6,7,2])))
end=timer()
print("\n",end-start)

#VIDIMO, DA JE MOJA FUNKCIJA DOSTI POČASNEJŠA

#%%
#----------------------------------------------------------------------------------------------------------
#       1. NALOGA - PRENESEMO PODATKE
#----------------------------------------------------------------------------------------------------------
with open('val2.dat') as f:
    podaci=[l.strip().split("\t") for l in f]
signal2=[float(podaci[i][0]) for i in range(len(podaci))]

with open('val3.dat') as f:
    podaci=[l.strip().split("\t") for l in f]
signal3=[float(podaci[i][0]) for i in range(len(podaci))]


t=[i for i in range(len(signal2)) ]

plt.figure(0)
plt.title('signal val2')
plt.plot(t,signal2,'r-')
plt.xlabel('Čas')
plt.ylabel('Signal')

plt.figure(1)
plt.title('signal val3')
plt.plot(t,signal3,'r-')
plt.xlabel('Čas')
plt.ylabel('Signal')

#%%
#----------------------------------------------------------------------------------------------------------
#       1. NALOGA - FFT in razna okna
#----------------------------------------------------------------------------------------------------------

plt.figure(0)
plt.title('spekter val2')
plt.plot(t,abs(np.fft.fft(signal2)),'g-',label='čist')

d=len(signal2)

Bartlett2=[signal2[i]*(1-abs((i-d/2)/(d/2))) for i in range(d)]
plt.plot(t,abs(np.fft.fft(Bartlett2)),'b-',label='Bartlett')

Hann2=[signal2[i]*(1-np.cos(2*np.pi*i/d))/2 for i in range(d)]
plt.plot(t,abs(np.fft.fft(Hann2)),'r-',label='Hann')

Welch2=[signal2[i]*(1-((i-d/2)/(d/2))**2) for i in range(d)]
plt.plot(t,abs(np.fft.fft(Welch2)),'y-',label='Welch')

plt.xlabel('Frekvenca')
plt.ylabel('Amplituda')
plt.legend()

plt.figure(2)
plt.title('signal(spekter(signal)) za val2')
plt.plot(t,np.fft.fft(abs(np.fft.fft(signal2))),'g-',label='čist')#bolje je uporabiti np.fft.ifft()
plt.xlabel('Čas')
plt.ylabel('Signal')
plt.legend()
#------------------------------------------------------

plt.figure(1)
plt.title('spekter val3')
plt.plot(t,abs(np.fft.fft(signal3)),'g-',label='čist')

Bartlett3=[signal3[i]*(1-abs((i-d/2)/(d/2))) for i in range(d)]
plt.plot(t,abs(np.fft.fft(Bartlett3)),'b-',label='Bartlett')

Hann3=[signal3[i]*(1-np.cos(2*np.pi*i/d))/2 for i in range(d)]
plt.plot(t,abs(np.fft.fft(Hann3)),'r-',label='Hann')

Welch3=[signal3[i]*(1-((i-d/2)/(d/2))**2) for i in range(d)]
plt.plot(t,abs(np.fft.fft(Welch3)),'y-',label='Welch')

plt.xlabel('Frekvenca')
plt.ylabel('Amplituda')
plt.legend()

plt.figure(3)
plt.title('signal(spekter(signal)) za val3')
plt.plot(t,np.fft.fft(abs(np.fft.fft(signal3))),'g-',label='čist')
plt.xlabel('Čas')
plt.ylabel('Signal')
plt.legend()

#%%
#----------------------------------------------------------------------------------------------------------
#       1. NALOGA - RAZREDČEN SIGNAL - FFT in razna okna 
#----------------------------------------------------------------------------------------------------------

t=[i+1 for i in range(len([signal2[i] for i in range(128)])) ]

plt.figure(0)

plt.title('spekter val2')
plt.plot(t,abs(np.fft.fft([signal2[i] for i in range(128)])),'g-',label='čist')

d=len([signal2[i] for i in range(128)])

Bartlett2=[signal2[i]*(1-abs((i-d/2)/(d/2))) for i in range(d)]
plt.plot(t,abs(np.fft.fft(Bartlett2)),'b-',label='Bartlett')

Hann2=[signal2[i]*(1-np.cos(2*np.pi*i/d))/2 for i in range(d)]
plt.plot(t,abs(np.fft.fft(Hann2)),'r-',label='Hann')

Welch2=[signal2[i]*(1-((i-d/2)/(d/2))**2) for i in range(d)]
plt.plot(t,abs(np.fft.fft(Welch2)),'y-',label='Welch')

plt.xlabel('Frekvenca')
plt.ylabel('Amplituda')
plt.legend()

#------------------------------------------------------


plt.figure(1)
plt.title('spekter val3')
plt.plot(t,abs(np.fft.fft([signal3[i] for i in range(128)])),'g-',label='čist')

Bartlett3=[signal3[i]*(1-abs((i-d/2)/(d/2))) for i in range(d)]
plt.plot(t,abs(np.fft.fft(Bartlett3)),'b-',label='Bartlett')

Hann3=[signal3[i]*(1-np.cos(2*np.pi*i/d))/2 for i in range(d)]
plt.plot(t,abs(np.fft.fft(Hann3)),'r-',label='Hann')

Welch3=[signal3[i]*(1-((i-d/2)/(d/2))**2) for i in range(d)]
plt.plot(t,abs(np.fft.fft(Welch3)),'y-',label='Welch')

plt.xlabel('Frekvenca')
plt.ylabel('Amplituda')
plt.legend()

#%%
###########################################################################################################
#
#       definicija funkcije za spektralno gostoto moči
#
###########################################################################################################

def abi(seznam, d):
    return [(abs(seznam[i])**2+abs(seznam[-i])**2)/2 for i in range(int(d/2))]

#%%
#----------------------------------------------------------------------------------------------------------
#       1. NALOGA - FFT in razna okna P
#----------------------------------------------------------------------------------------------------------

d=len(signal2)
t=[i+1 for i in range(int(d/2)) ]

plt.figure(0)
plt.title('spekter val2')
plt.plot(t,abi(np.fft.fft(signal2),d),'g-',label='čist')



Bartlett2=[signal2[i]*(1-abs((i-d/2)/(d/2))) for i in range(d)]
plt.plot(t,abi(np.fft.fft(Bartlett2),d),'b-',label='Bartlett')

Hann2=[signal2[i]*(1-np.cos(2*np.pi*i/d))/2 for i in range(d)]
plt.plot(t,abi(np.fft.fft(Hann2),d),'r-',label='Hann')

Welch2=[signal2[i]*(1-((i-d/2)/(d/2))**2) for i in range(d)]
plt.plot(t,abi(np.fft.fft(Welch2),d),'y-',label='Welch')

plt.xlabel('Frekvenca')
plt.ylabel('Amplituda')
plt.legend()

plt.figure(2)
plt.title('signal(spekter(signal)) za val2')
plt.plot(t,np.fft.fft(abi(np.fft.fft(signal2),d)),'g-',label='čist')#bolje je uporabiti np.fft.ifft()
plt.xlabel('Čas')
plt.ylabel('Signal')
plt.legend()
#------------------------------------------------------

plt.figure(1)
plt.title('spekter val3')
plt.plot(t,abi(np.fft.fft(signal3),d),'g-',label='čist')

Bartlett3=[signal3[i]*(1-abs((i-d/2)/(d/2))) for i in range(d)]
plt.plot(t,abi(np.fft.fft(Bartlett3),d),'b-',label='Bartlett')

Hann3=[signal3[i]*(1-np.cos(2*np.pi*i/d))/2 for i in range(d)]
plt.plot(t,abi(np.fft.fft(Hann3),d),'r-',label='Hann')

Welch3=[signal3[i]*(1-((i-d/2)/(d/2))**2) for i in range(d)]
plt.plot(t,abi(np.fft.fft(Welch3),d),'y-',label='Welch')

plt.xlabel('Frekvenca')
plt.ylabel('Amplituda')
plt.legend()

plt.figure(3)
plt.title('signal(spekter(signal)) za val3')
plt.plot(t,np.fft.fft(abi(np.fft.fft(signal3),d)),'g-',label='čist')
plt.xlabel('Čas')
plt.ylabel('Signal')
plt.legend()

#%%
#----------------------------------------------------------------------------------------------------------
#       1. NALOGA - RAZREDČEN SIGNAL - FFT in razna okna - P - to je še potrebno napisat. uporabi abi, ampak verjetno je trebe še kaj.. :/
#----------------------------------------------------------------------------------------------------------

t=[i+1 for i in range(len([signal2[i] for i in range(32)])) ]

plt.figure(0)

plt.title('spekter val2 N=64')
plt.plot(t,abi(np.fft.fft([signal2[i] for i in range(64)]),64),'g-',label='čist')

d=len([signal2[i] for i in range(64)])

Bartlett2=[signal2[i]*(1-abs((i-d/2)/(d/2))) for i in range(d)]
plt.plot(t,abi(np.fft.fft(Bartlett2),d),'b-',label='Bartlett')

Hann2=[signal2[i]*(1-np.cos(2*np.pi*i/d))/2 for i in range(d)]
plt.plot(t,abi(np.fft.fft(Hann2),d),'r-',label='Hann')

Welch2=[signal2[i]*(1-((i-d/2)/(d/2))**2) for i in range(d)]
plt.plot(t,abi(np.fft.fft(Welch2),d),'y-',label='Welch')

plt.xlabel('Frekvenca')
plt.ylabel('Amplituda')
plt.legend()

#------------------------------------------------------


plt.figure(1)
plt.title('spekter val3 N=64')
plt.plot(t,abi(np.fft.fft([signal3[i] for i in range(64)]),d),'g-',label='čist')

Bartlett3=[signal3[i]*(1-abs((i-d/2)/(d/2))) for i in range(d)]
plt.plot(t,abi(np.fft.fft(Bartlett3),d),'b-',label='Bartlett')

Hann3=[signal3[i]*(1-np.cos(2*np.pi*i/d))/2 for i in range(d)]
plt.plot(t,abi(np.fft.fft(Hann3),d),'r-',label='Hann')

Welch3=[signal3[i]*(1-((i-d/2)/(d/2))**2) for i in range(d)]
plt.plot(t,abi(np.fft.fft(Welch3),d),'y-',label='Welch')

plt.xlabel('Frekvenca')
plt.ylabel('Amplituda')
plt.legend()
#%%
#%%
#%%
#%%
###########################################################################################################
#
#       2. NALOGA - PRENESEMO PODATKE
#
###########################################################################################################

import seaborn as sns
sns.set_palette(sns.color_palette("Set1", 4))

sigi=[]
for k in range(4):
    with open('signal{}.dat'.format(k)) as f:
        podaci=[l.strip().split("\t") for l in f]
    sigi.append([float(podaci[i][0]) for i in range(len(podaci))])

for i in range(4):
    plt.plot([k+1 for k in range(len(sigi[i]))],sigi[i],label='signal{}'.format(i))
plt.legend()

#%%
#----------------------------------------------------------------------------------------------------------
#       2. NALOGA - fft signala 0 z vsemi okni
#----------------------------------------------------------------------------------------------------------

d=len(sigi[0])

t=[i for i in range(d) ]
povprečje=sum(sigi[0][i] for i in range(d))/d


plt.figure(0)
plt.title('spekter signal0')
plt.plot(t,abs(np.fft.fft([sigi[0][i]-povprečje for i in range(d)])),label='čist')

Bartlett2=[sigi[0][i]*(1-abs((i-d/2)/(d/2))) for i in range(d)]
plt.plot(t,abs(np.fft.fft(Bartlett2)),label='Bartlett')

Hann2=[sigi[0][i]*(1-np.cos(2*np.pi*i/d))/2 for i in range(d)]
plt.plot(t,abs(np.fft.fft(Hann2)),label='Hann')

Welch2=[sigi[0][i]*(1-((i-d/2)/(d/2))**2) for i in range(d)]
plt.plot(t,abs(np.fft.fft(Welch2)),label='Welch')

plt.xlabel('Frekvenca')
plt.ylabel('Amplituda')
plt.legend()

#%%
#----------------------------------------------------------------------------------------------------------
#       2. NALOGA - luščenje signala  0
#----------------------------------------------------------------------------------------------------------
d=len(sigi[0])

def r(t,N):
    if t<=N/2:
        return 1/32*np.exp(-abs(t)/16)
    if t>N/2:
        return 1/32*np.exp(-abs(t-N)/16)

filt=[r(i,d) for i in range(d)]
R=np.fft.fft(filt)
spekter0=np.fft.fft(sigi[0])
deljenec=[spekter0[i]/R[i] for i in range(d)]
vhodni=np.fft.ifft(deljenec)

plt.title('Luščenje vhodnega signala iz izhodnega brez šuma (signal0)')
plt.plot([i for i in range(d)],sigi[0],label='Izhodni signal')
plt.plot([i for i in range(d)],vhodni,label='Izluščen vhodni signal')
plt.xlabel('Čas')
plt.ylabel('Signal')
plt.legend()


#%%
#----------------------------------------------------------------------------------------------------------
#       2. NALOGA - fft vseh signalov
#----------------------------------------------------------------------------------------------------------


d=len(sigi[0])

t=[i for i in range(d) ]
povprečje=sum(sigi[0][i] for i in range(d))/d


plt.figure(0)
for g in range(4):
    plt.title('spekter za vse signale brez okna'.format(g))
    plt.plot(t,abs(np.fft.fft([sigi[g][i]-povprečje for i in range(d)])),label='signal{}'.format(g))

plt.xlabel('Frekvenca')
plt.ylabel('Amplituda')
plt.legend()


#%%
###########################################################################################################
#
#       2. NALOGA - luščenje signala - sumljivi model
#
###########################################################################################################

d=len(sigi[2])
saga=np.fft.fft(sigi[2])
N=10

utilde=[saga[i]*((abs(saga[i])**2-N*N)/(abs(saga[i])**2))/R[i] for i in range(d)]

#%%


plt.title('Ocena za vhodni spekter (signal2)')
plt.plot([i for i in range(d)],sigi[2],label='Izhodni signal')
plt.plot([i for i in range(d)],np.fft.ifft(utilde),label='Izluščen vhodni signal')
plt.xlabel('Čas')
plt.ylabel('Signal')
plt.legend()


#%%
#----------------------------------------------------------------------------------------------------------
#       2. NALOGA - luščenje signala - težji model
#----------------------------------------------------------------------------------------------------------

d=len(sigi[3])
saga=np.fft.fft(sigi[1])
filt=[r(i,d) for i in range(d)]
R=np.fft.fft(filt)
saga2=np.fft.fft(sigi[2])

tisto=np.array([np.log10(abs(saga2[i+1])) for i in range(30)])#logaritem
t=np.array([i for i in range(30)])
k=np.polyfit(t,tisto,1)

S=[np.exp(k[0]*i+k[1]) if i <= d/2 else np.exp(k[0]*(d-i)+k[1]) for i in range(d)]
N=sum(abs(saga[i]) for i in [j+30 for j in range(170)])/170#logaritem


utilde=[saga[i]*(abs(S[i])**2/(abs(S[i])**2+abs(N)**2))/R[i] for i in range(d)]


#%%

plt.title('Ocena za vhodni spekter (signal3)')
plt.plot([i for i in range(d)],sigi[0],label='Izhodni signal')
plt.plot([i for i in range(d)],np.fft.ifft(utilde),label='Izluščen vhodni signal')
plt.xlabel('Čas')
plt.ylabel('Signal')
plt.legend()

#%%
#%%
#%%
#%%
###########################################################################################################
#
#       3. NALOGA - VRSTIČNO ODBIRANJE LINCOLNA
#
###########################################################################################################

#----------------------------------------------------------------------------------------------------------
#       3. NALOGA - začetniški fail
#----------------------------------------------------------------------------------------------------------

figi=[]
vmes=[]
for k in range(5):
    with open('lincoln_L30_N{}0.pgm'.format(k)) as f:
        podaci=[l.strip().split(" ") for l in f]
    podaci=[[int(podaci[i][j]) for j in range(len(podaci[i])) if i > 0] for i in range(len(podaci))]
    podaci=[podaci[i][j]  for i in range(len(podaci)) for j in range(len(podaci[i]))]
    vmes.append(podaci)
    slika=np.array([[podaci[i+j]/256 for j in range(256)] for i in range(313)])
    slika=[[slika[i][j] if slika[i][j]<=1 else 0 for j in range(len(slika[3]))] for i in range(len(slika))]
    #figi.append([[int(podaci[i][0])] for i in range(len(podaci))])
    figi.append(slika)


imgplot = plt.imshow(figi[0])


#img=mpimg.imread('lincoln_L30_N00.pgm')
plt.figure(10)
plt.plot([i for i in range(len(vmes[0]))],vmes[0])


#%%
#----------------------------------------------------------------------------------------------------------
#       3. NALOGA - legit importane slike v .png formatu
#----------------------------------------------------------------------------------------------------------

figi=[]
slikce=[]
for i in range(5):
    img=mpimg.imread("lincoln_L30_N{}0.png".format(i))
    x=Image.open( "lincoln_L30_N{}0.png".format(i) )
    figi.append(img)
    slikce.append(x)
plt.figure(10)
imgplot = plt.imshow(figi[3],cmap='copper')
#%%
# definiramo novo funkcijo

def r2(t):
    return 1/30*np.exp(-t/30) # ta bi mogla bit ta prava...

#%%
# mali tester
    
plt.plot([i for i in range(len(R))],filt)

#%%
#----------------------------------------------------------------------------------------------------------
#       3. NALOGA - STOLPCI ta se je izkazala za uspešno
#----------------------------------------------------------------------------------------------------------


filt=[r2(i) for i in range(256)]
R=np.fft.fft(filt)

a=[]
for i in range(313):
    mati=figi[0]
    trololo=[mati[j][i] for j in range(256)]
    bebe=np.fft.fft(trololo)
    a.append(abs(np.fft.ifft([bebe[k]/R[k] for k in range(len(R))])))
    
m=np.max(a)
a=np.array([[a[j][i]/m if j>0 and i>0 else 0 for j in range(len(a))] for i in range(len(a[2]))])
#a=np.array(a)

plt.figure(3)
imgplot = plt.imshow(a,cmap='Greys_r')


#plt.figure(2)
#plt.plot([i for i in range(256)],np.real(np.fft.ifft([bebe[k]/R[k] for k in range(len(R))])))
#%%
#----------------------------------------------------------------------------------------------------------
#       3. NALOGA - OBOJE ta ni zares smiselna
#----------------------------------------------------------------------------------------------------------

filt=[r2(i) for i in range(256)]
R=np.fft.fft(filt)

filt2=[r2(i) for i in range(313)]
R2=np.fft.fft(filt2)

a=[]
for i in range(313):
    mati=figi[0]
    trololo=[mati[j][i] for j in range(256)]
    bebe=np.fft.fft(trololo)
    a.append(abs(np.fft.ifft([bebe[k]/R[k] for k in range(len(R))])))

a=np.array([[a[j][i] if j>0 and i>0 else 0 for j in range(len(a))] for i in range(len(a[2]))])
A=[]
for i in range(256):
    trololo=[a[i][j] for j in range(313)]
    bebe=np.fft.fft(trololo)
    A.append(abs(np.fft.ifft([bebe[k]/R2[k] for k in range(len(R2))])))
    
#m=np.max(a)
A=np.array([[A[j][i] if j>0 and i>0 else 0 for j in range(len(A))] for i in range(len(A[2]))])
#a=np.array(a)


imgplot = plt.imshow(A,cmap='Greys_r')
#%%
#----------------------------------------------------------------------------------------------------------
#       3. NALOGA - STOLPCI za zahtevnejše primere tudi N(f)
#----------------------------------------------------------------------------------------------------------


filt=[r2(i) for i in range(256)]
R=np.fft.fft(filt)

a=[]
for i in range(313):
    
    saga=np.array(np.fft.fft([figi[4][j][i] for j in range(256)]))
    
    
    titi=np.array([np.log10(abs(np.fft.fft([figi[4][i][150] for i in range(256)])))[i+2] for i in range(20)])# logaritem
    t=np.array([i+2 for i in range(20)])
    k=np.polyfit(t,titi,1)
    
    titi2=np.array([np.log10(abs(np.fft.fft([figi[4][i][150]-figi[0][i][150] for i in range(256)])))[i+50] for i in range(70)])# logaritem
    t2=np.array([i+50 for i in range(70)])
    k2=np.polyfit(t2,titi2,1)
    
    
    S=[np.exp(k[0]*i+k[1])  for i in range(256)]
    N=[np.exp(k2[0]*i+k2[1]) if i <= 256/2 else np.exp(k2[0]*(256-i)+k2[1]) for i in range(256)]
#    N=sum(abs(saga[i]) for i in [j+50 for j in range(50)])/50
    
    
    
    utilde=[saga[i]*(abs(S[i])**2/(abs(S[i])**2+abs(N[i])**2))/R[i] for i in range(256)]
    a.append(abs(np.fft.ifft(utilde)))
    
plt.figure(8)
titi=[abs(np.fft.fft([figi[2][i][10] for i in range(256)]))[i+2] for i in range(254)]
plt.plot([i for i in range(254)],titi)
    

a=np.array([[a[j][i] if a[j][i]<1 else 1 for j in range(len(a))] for i in range(len(a[2]))])
#a=np.array(a)

plt.figure(5)
imgplot = plt.imshow(a,cmap='copper')


#%%
plt.figure(100)
plt.plot([i for i in range(256)],[a[i][200] for i in range(256)])
#if i <= 256/2 else np.exp(k[0]*(256-i)+k[1])


#%%

plt.figure(8)
plt.title('Spektri za stolpec 10')
for j in [2]:
    titi=[abs(np.fft.fft([figi[j][i][10] for i in range(256)]))[i+2] for i in range(254)]
    plt.plot([i for i in range(254)],titi,label='slika{}'.format(j))
plt.xlabel('Frekvenca')
plt.ylabel('Spekter')
plt.legend()

#%%
#----------------------------------------------------------------------------------------------------------
#       3. NALOGA - STOLPCI vsak posebej tudi N(f)
#----------------------------------------------------------------------------------------------------------


filt=[r2(i) for i in range(256)]
R=np.fft.fft(filt)

a=[]
for i in range(313):
    
    saga=np.array(np.fft.fft([figi[0][j][i] for j in range(256)]))
    
    
    titi=np.array([np.log10(abs(np.fft.fft([figi[0][j][i] for j in range(256)])))[j+2] for j in range(50)])# logaritem
    t=np.array([i+2 for i in range(50)])
    k=np.polyfit(t,titi,1)
    
    titi2=np.array([np.log10(abs(np.fft.fft([figi[0][j][i] for j in range(256)])))[j+70] for j in range(50)])# logaritem
    t2=np.array([i+70 for i in range(50)])
    k2=np.polyfit(t2,titi2,1)
    
    
    S=[np.exp(k[0]*i+k[1]) if i <= 256/2 else np.exp(k[0]*(256-i)+k[1]) for i in range(256)]
    N=[np.exp(k2[0]*i+k2[1]) if i <= 256/2 else np.exp(k2[0]*(256-i)+k2[1]) for i in range(256)]
#    N=sum(abs(saga[i]) for i in [j+50 for j in range(50)])/50
    
    
    
    utilde=[saga[i]*(abs(S[i])**2/(abs(S[i])**2+abs(N[i])**2))/R[i] for i in range(256)]
    a.append(abs(np.fft.ifft(utilde)))
    

titi=[abs(np.fft.fft([figi[0][i][10] for i in range(256)]))[i+2] for i in range(254)]
plt.plot([i for i in range(254)],titi)
    

a=np.array([[a[j][i] if i%253>3 else 1 for j in range(len(a))] for i in range(len(a[2]))])
#a=np.array(a)

plt.figure(4)
imgplot = plt.imshow(a,cmap='Greys_r')


#%%



















