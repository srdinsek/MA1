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
from scipy.linalg import solve_toeplitz
from scipy.optimize import root


#%%
#%%
#%%
###########################################################################################################
#
#       1. NALOGA - ZAPIŠEMO MATRIKO AVTOKORELACIJ IN REŠIMO ZA A in Z
#
###########################################################################################################


def R(i,seznam,N):
    return sum(seznam[n]*seznam[n+i]/(i-N) for n in range(N-2-i))
    
def M(seznam,p,N):
    return np.array([[R(abs(i-j),seznam,N) for i in range(p)] for j in range(p)])
def b(seznam,p,N):
    return np.array([-R(i+1,seznam,N) for i in range(p)])

def a(seznam,p):
    N=len(seznam)
    return np.linalg.solve(M(seznam,p,N), b(seznam,p,N))

def a_h(seznam,p):
    N=len(seznam)
    c=[R(i,seznam,N) for i in range(p)]
    b=[-R(i+1,seznam,N) for i in range(p)]
    return linalg.solve_toeplitz(c, b)

def z(ak):
    p=len(ak)
    return (np.roots([ak[p-1-i] if i<p else 1 for i in range(p+1)]))**-1

def zabs(zk,ak):
    p=len(ak)
    z=[(zk[i]/abs(zk[i]))**-1 for i in range(p)]
    rešitev=np.poly1d(z, True)
    return [rešitev[i+1]/rešitev[0] for i in range(p)]

def P(ni,ak):
    N=256
    p=len(ak)
    return 1/(abs(1+sum(ak[i]*np.exp(-1j*ni*2*np.pi*(i+1)/N) for i in range(p)))**2)

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

with open('co22.txt') as f:
    podaci=[l.strip().split(" ") for l in f]
co=[[float(podaci[i][0]),float(podaci[i][1])] for i in range(len(podaci)) if float(podaci[i][1])!=-99.99]
co2=[float(podaci[i][1]) for i in range(len(podaci)) if float(podaci[i][1])!=-99.99] 

with open('borza.dat') as f:
    podaci=[l.strip().split("\t") for l in f]
borza=[float(podaci[i][0]) for i in range(len(podaci))]

with open('luna1.txt') as f:
    podaci=[l.strip().split(" ") for l in f]
luna=[[float(podaci[i][1]),float(podaci[i][2])] for i in range(len(podaci))]
luna=[luna[i][1] for i in range(len(luna))]

with open('Wolf_number.dat') as f:
    podaci=[l.strip().split(" ") for l in f]
wolf=[[float(podaci[j][i]) for i in range(len(podaci[j])) if podaci[j][i]!=''] for j in range(len(podaci))]
wolf=[wolf[i][2] for i in range(len(wolf))]

#%%
#----------------------------------------------------------------------------------------------------------
#       1. NALOGA - REŠIMO MATRIČNO ENAČBO - primerjamo hitrosti metod
#----------------------------------------------------------------------------------------------------------

start=timer()
print(a(signal2,6))
end=timer()
print('čas',end-start)# ta traja 5.14 sekund (250)

start=timer()
print('\n\n',a_h(signal2,8))
end=timer()
print('čas',end-start)# ta traja 0.040 sekund (250)


#%%
#----------------------------------------------------------------------------------------------------------
#       1. NALOGA - REŠIMO SISTEM
#----------------------------------------------------------------------------------------------------------


ali2=a_h(signal2,5)
zk2=z(ali2)
alin2=zabs(zk2,ali2)

ali3=a_h(signal3,5)
zk3=z(ali3)
alin3=zabs(zk3,ali3)

aliCO2=a_h(co2,3)
zkCO2=z(aliCO2)
alinCO2=zabs(zkCO2,aliCO2)

#%%
###########################################################################################################
#       1. NALOGA - RIŠEMO KROG Z REŠITVAMI
###########################################################################################################

circle1 = plt.Circle((0, 0), 1, color='black',fill=False)
fig,ax = plt.subplots()

ax.set_aspect("equal")
for i in [5,6,8,20]:
    ali2=a_h(signal2,i)
    zk2=z(ali2)
    ax.scatter(zk2.real,zk2.imag,label='p={}'.format(i))
plt.plot()
ax.add_artist(circle1)

plt.title('Poli za val2')
plt.xlabel('Realni del')
plt.ylabel('Imaginarni del')
plt.legend()



#%%
circle1 = plt.Circle((0, 0), 1, color='black',fill=False)
fig,ax = plt.subplots()

ax.set_aspect("equal")
for i in [5,6,8,20]:
    ali3=a_h(signal3,i)
    zk3=z(ali3)
    ax.scatter(zk3.real,zk3.imag,label='p={}'.format(i))
plt.plot()
ax.add_artist(circle1)

plt.title('Poli za val3')
plt.xlabel('Realni del')
plt.ylabel('Imaginarni del')
plt.legend()

#%%

circle1 = plt.Circle((0, 0), 1, color='black',fill=False)
fig,ax = plt.subplots()

ax.set_aspect("equal")
ax.scatter(zkCO2.real,zkCO2.imag)
plt.plot()
ax.add_artist(circle1)

plt.title('Poli za p=250 CO2')
plt.xlabel('Realni del')
plt.ylabel('Imaginarni del')

#%%
###########################################################################################################
#       1. NALOGA - RIŠEMO SPEKTER
###########################################################################################################

# val2

plt.title('Spekter P($\omega$) val2')
for i in range(1):
    ali2=a_h(signal2,(i+1)*4)
    zk2=z(ali2)
    alin2=zabs(zk2,ali2)
    plt.plot([i for i in range(int(512/2))],[P(i,ali2) for i in range(int(512/2))],label='P($\omega$) za p={}'.format((i+1)*4))
    plt.plot([i for i in range(int(512/2))],[P(i,alin2) for i in range(int(512/2))],'g',label='P($\omega$) za p={} normirano'.format((i+1)*4))

speki=abs(np.fft.fft(signal2))
plt.plot([i for i in range(int(512/2))],[(speki[i]**2+speki[-i]**2)/2 for i in range(int(512/2))],'b',alpha=0.4,label='P(val2)')
plt.xlabel('$\omega$')
plt.ylabel('log[P($\omega$)]')
plt.legend()
plt.yscale('symlog')


#%%

# val3

plt.title('Spekter P($\omega$) val3')
for i in [13]:
    ali3=a_h(signal3,i)
    zk3=z(ali3)
    alin3=zabs(zk3,ali3)
    plt.plot([i for i in range(int(512/2))],[100*P(i,ali3) for i in range(int(512/2))],label='100*P($\omega$) za p={}'.format(i))
    plt.plot([i for i in range(int(512/2))],[100*P(i,alin3) for i in range(int(512/2))],'g',label='100*P($\omega$) za p={} normirano'.format(i))


#plt.plot([i for i in range(int(512/2))],[10000*P(i,ali3) for i in range(int(512/2))],label='P($\omega$)')
speki=abs(np.fft.fft(signal3))
plt.plot([i for i in range(int(512/2))],[(speki[i]**2+speki[-i]**2)/2 for i in range(int(512/2))],'b',label='P(val3)')
plt.xlabel('$\omega$')
plt.ylabel('P($\omega$)')
plt.legend()
plt.yscale('symlog')

#%%
###########################################################################################################
#       1. NALOGA - SUPER TEST
###########################################################################################################

sigolino=[np.sin(5*2*np.pi*k/256)+np.sin(10*2*np.pi*k/256)+np.sin(15*2*np.pi*k/256) for k in range(256)]

speki=abs(np.fft.fft(sigolino))
Piki=[(speki[i]**2+speki[-i]**2)/2 for i in range(int(len(speki)/2))]
plt.plot([i for i in range(len(Piki))],Piki)

aligo=a_h(sigolino,100)
zkgo=z(aligo)
alingo=zabs(zkgo,aligo)

plt.plot([i for i in range(int(len(sigolino)/2))],[P(i,aligo) for i in range(int(len(sigolino)/2))],label='100*P($\omega$) za p={}'.format(5))
plt.yscale('symlog')


#%%
###########################################################################################################
#       1. NALOGA - CO2
###########################################################################################################
# co2

d=len(co2)
tisto=np.array(co2)
t=np.array([i for i in range(d)])
k=np.polyfit(t,tisto,1)
Sol=[co2[i]-(k[0]*i+k[1]) for i in range(d)]

plt.title('Spekter P($\omega$) val3')
for i in [30,90]:
    aliCO2=a_h(Sol,i)
    zkCO2=z(aliCO2)
    alinCO2=zabs(zkCO2,aliCO2)
#    plt.plot([i for i in range(int(512/2))],[100*P(i,aliCO2) for i in range(int(512/2))],label='100*P($\omega$) za p={}'.format((i+1)*9))
    plt.plot([i for i in range(int(512/2))],[10*P(i,aliCO2) for i in range(int(512/2))],label='10*P($\omega$) za p={}'.format(i))


d=len(co2)
speki=abs(np.fft.fft([co2[i]*(1-np.cos(2*np.pi*i/d))/2 for i in range(d)]))
plt.plot([i for i in range(int(512/2))],[(speki[i]**2+speki[-i]**2)/2 for i in range(int(512/2))],'b',alpha=0.4,label='P(CO2)')
plt.title('P($\omega$) za CO2 Hann')
plt.xlabel('$\omega$')
plt.ylabel('P($\omega$)')
plt.legend()
plt.yscale('symlog')

#%%
#----------------------------------------------------------------------------------------------------------
#       1. NALOGA - OKNA IN FFT
#----------------------------------------------------------------------------------------------------------


#Bartlett2=[co2[i]*(1-abs((i-d/2)/(d/2))) for i in range(d)]
#Hann2=[signal2[i]*(1-np.cos(2*np.pi*i/d))/2 for i in range(d)]
#Welch2=[signal2[i]*(1-((i-d/2)/(d/2))**2) for i in range(d)]

speki=abs(np.fft.fft([Sol[i]*(1-abs((i-d/2)/(d/2))) for i in range(d)]))
plt.plot([i for i in range(int(512/2))],[(speki[i]**2+speki[-i]**2)/2 for i in range(int(512/2))],'b',label='Bartlett')
speki=abs(np.fft.fft([Sol[i]*(1-np.cos(2*np.pi*i/d))/2 for i in range(d)]))
plt.plot([i for i in range(int(512/2))],[(speki[i]**2+speki[-i]**2)/2 for i in range(int(512/2))],'r',label='Hann')
speki=abs(np.fft.fft([Sol[i]*(1-((i-d/2)/(d/2))**2) for i in range(d)]))
plt.plot([i for i in range(int(512/2))],[(speki[i]**2+speki[-i]**2)/2 for i in range(int(512/2))],'g',label='Welch')
plt.title('P($\omega$) za CO2')
plt.xlabel('$\omega$')
plt.ylabel('P($\omega$)')
plt.legend()
plt.yscale('symlog')

#%%
#----------------------------------------------------------------------------------------------------------
#       1. NALOGA - OD CO2 PODATKOV ODŠTEJEMO TREND
#----------------------------------------------------------------------------------------------------------

d=len(co2)
tisto=np.array(co2)
t=np.array([i for i in range(d)])
k=np.polyfit(t,tisto,1)
Sol=[co2[i]-(k[0]*i+k[1]) for i in range(d)]

plt.title('CO2')
plt.plot(t,[co2[i]-300 for i in range(d)],'r',label='CO2 - 300')
plt.plot(t,Sol,'b',label='CO2 - linearni fit')
plt.xlabel('Korak')
plt.ylabel('Koncentracija')
plt.legend()

#%%
#%%
#%%
#%%
###########################################################################################################
#
#       2. NALOGA - LINEARNA NAPOVED
#
###########################################################################################################
d=len(co2)
tisto=np.array(co2)
t=np.array([i for i in range(d)])
k=np.polyfit(t,tisto,1)
Sol=[co2[i]-(k[0]*i+k[1]) for i in range(d)]

#%%
#----------------------------------------------------------------------------------------------------------
#       2. NALOGA - DEFINIRAM FUNKCIJO
#----------------------------------------------------------------------------------------------------------

def S(ak, signal):
    Stil=[signal[i] for i in range(int(len(signal)/2))]
    N=len(Stil)
    p=len(ak)
    for n in range(int(len(signal)/2)-p):
        Stil.append(sum(-ak[i]*Stil[n+N-i-1] for i in range(p)))
    return Stil


#%%
#----------------------------------------------------------------------------------------------------------
#       2. NALOGA - RIŠEM NAPOVED ZA CO2
#----------------------------------------------------------------------------------------------------------    

d=len(Sol)
    
sns.set_palette(sns.color_palette("autumn",2))


plt.title('Linearna napoved brez normiranja')

for g in [60]:
    aliCO2=a_h([Sol[i] for i in range(int(d/2))],g)
    zkCO2=z(aliCO2)
    alinCO2=zabs(zkCO2,aliCO2)

    sigi=S(aliCO2,Sol)
    N=len(sigi)
    plt.plot([i for i in range(N)],sigi,'g',label='p={}'.format(g))

plt.plot([i for i in range(N)],[Sol[i] for i in range(N)],'b',[302,302],[-6,6],'b',alpha=0.4)
#plt.yscale('symlog')
plt.xlabel('Korak')
plt.ylabel('Koncentracija')
plt.legend()

#%%
#----------------------------------------------------------------------------------------------------------
#       2. NALOGA - RIŠEM PODATKE
#----------------------------------------------------------------------------------------------------------
plt.figure(1)
plt.title('Borza')
plt.plot([i for i in range(len(borza))],borza,color='green')
plt.xlabel('Korak')
plt.ylabel('Vrednost')
plt.figure(2)
plt.title('Luna')
plt.plot([i for i in range(len(luna))],luna,color='blue')
plt.xlabel('Korak')
plt.ylabel('Jakost')
plt.figure(3)
plt.title('Sončeve pege')
plt.plot([i for i in range(len(wolf))],wolf,color='red')
plt.xlabel('Korak')
plt.ylabel('Jakost')

#%%
#----------------------------------------------------------------------------------------------------------
#       2. NALOGA - RIŠEM NAPOVED ZA BORZO
#----------------------------------------------------------------------------------------------------------    

d=len(borza)
    
sns.set_palette(sns.color_palette("autumn",3))


plt.title('Linearna napoved z normiranjem')
for g in [10,30,60]:
    aliCO2=a_h([borza[i] for i in range(int(d/2))],g)
    zkCO2=z(aliCO2)
    alinCO2=zabs(zkCO2,aliCO2)

    sigi=S(aliCO2,borza)
    N=len(sigi)
    plt.plot([i for i in range(N)],sigi,label='p={}'.format(g),alpha=0.7)

plt.plot([i for i in range(N)],[borza[i] for i in range(N)],'b',[d/2,d/2],[0,16000],'b',alpha=0.4)
#plt.yscale('symlog')
plt.xlabel('Korak')
plt.ylabel('Koncentracija')
plt.legend()

#%%
#----------------------------------------------------------------------------------------------------------
#       2. NALOGA - RIŠEM NAPOVED ZA LUNA
#----------------------------------------------------------------------------------------------------------    

d=len(luna)
    



plt.title('Linearna napoved z normiranjem')
for g in [5,20,40,60]:
    aliCO2=a_h([luna[i] for i in range(int(d/2))],g)
    zkCO2=z(aliCO2)
    alinCO2=zabs(zkCO2,aliCO2)

    sigi=S(alinCO2,luna)
    N=len(sigi)
    t=1
    if g==60:
        t=0.3
    plt.plot([i for i in range(N)],sigi,label='p={}'.format(g),alpha=t)

plt.plot([i for i in range(N)],[luna[i] for i in range(N)],'y',[d/2,d/2],[-25,25],'b',alpha=0.4)
#plt.yscale('symlog')
plt.xlabel('Korak')
plt.ylabel('Koncentracija')
plt.legend()

#%%
#----------------------------------------------------------------------------------------------------------
#       2. NALOGA - RIŠEM NAPOVED ZA SONCE
#----------------------------------------------------------------------------------------------------------    

d=len(wolf)

plt.title('Linearna napoved brez normiranja')
for g in [600,710]:
    aliCO2=a_h([wolf[i] for i in range(int(d/2))],g)
    zkCO2=z(aliCO2)
    alinCO2=zabs(zkCO2,aliCO2)

    sigi=S(aliCO2,wolf)
    N=len(sigi)
    t=1
    plt.plot([i for i in range(N)],sigi,label='p={}'.format(g),alpha=t)

plt.plot([i for i in range(N)],[wolf[i] for i in range(N)],'y',[d/2,d/2],[-25,250],'b',alpha=0.4)
#plt.yscale('symlog')
plt.xlabel('Korak')
plt.ylabel('Koncentracija')
plt.legend()

#%%

















