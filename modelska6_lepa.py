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
from numpy.polynomial.chebyshev import chebvander2d
from numpy.polynomial.chebyshev import chebfit
from numpy.polynomial.legendre import legvander2d
from numpy.polynomial.hermite import hermvander2d
from scipy.stats import norm
from scipy.special import mathieu_cem, mathieu_sem



import seaborn as sns
sns.set_palette(sns.color_palette("cool", 10))

#%%# tukaj vnesem podatke in jih preuredim v matriko A in vektor b, ker rašim matrično enačbo. x predstavlja rešitev prve naloge
###########################################################################################################
#
#       1. NALOGA - vnesem podatke (x je rešitev)
#
###########################################################################################################

with open('farmakoloski1.txt') as f:
    podaci=[l.strip().split("\t") for l in f]   
seznam=[[float(podaci[i][j]) for j in range(len(podaci[i]))] for i in range(len(podaci)) if i >> 0]
dolzina=len(seznam)
fi=[[1/seznam[i][0] if j==1 else 1 for j in [0,1]] for i in range(dolzina)]
A=[[sum((seznam[i][1]**2)*fi[i][j]*fi[i][k]/3 for i in range(dolzina)) for j in range(len(fi[1]))] for k in range(len(fi[1]))]
b=[sum((1/seznam[i][1])*fi[i][k]*seznam[i][1]**2/3 for i in range(dolzina)) for k in [0,1]]
x=solve(A,b)



a=sum( (seznam[i][1]**2/3)**2 for i in range(dolzina))
Ax=sum( fi[i][0]*(seznam[i][1]**2/3)**2 for i in range(dolzina))
Ay=sum( fi[i][1]*(seznam[i][1]**2/3)**2 for i in range(dolzina))

a_1= (A[0][0]*Ay-Ax*A[0][1])/(A[0][0]*a-Ax**2)
a_2= (a*A[0][1]-Ax*Ay)/(A[0][0]*a-Ax**2)



print('\t \t')
print('Imamo matriko:\t', A)
print('In vektor b:\t\t', b)
print('Ko tako enačbo rešimo, dobimo rešitev:\t', x, 'pri čemer, so konstante \n\n a =', x[1]/x[0], 'in y0 =', 1/x[0])
print('\n\n Če pa jih rešimo analitično dobimo za parametre vrednsoti a_1;', a_1 ,'\t in za a_2:\t', a_2)

#%% tukaj rišem graf of prve naloge
#----------------------------------------------------------------------------------------------
#       1. NALOGA - rišem normalno rešitev z napakami
#----------------------------------------------------------------------------------------------


#časovna os
t2 = np.arange(0.1, 1000, 0.05)

#ostale osi in arrayi potrebni za grafe
iks=[seznam[i][0] for i in range(dolzina)]
ipsilon=[seznam[i][1] for i in range(dolzina)]
e=[3 for i in range(dolzina)]

# \chi^2
ksi=sum((seznam[i][1]-(seznam[i][0]/x[0])/(seznam[i][0]+x[1]/x[0]))**2/9 for i in range(dolzina))/5


#izrišemo
plt.plot(iks,ipsilon,'ro',label="$\chi^2/(m-n)$: {0:.3f}".format(ksi))#ta je samo zato, da napiše v legendo \chi^2
plt.plot(t2,(t2/x[0])/(t2+x[1]/x[0]),'b--',label='prilagoditvena funkcija')#zvezna prilagoditvena funkcija
plt.errorbar(iks,ipsilon,yerr=e,color='gray',fmt='o',label='meritve')#točke meritve
plt.title('Prilagoditvena funkcija za dane meritve')
plt.ylabel('y')
plt.xlabel('x')
plt.legend()

#%%
#----------------------------------------------------------------------------------------------
#       1. NALOGA - rišem linearizirano rešitev z lineariziranimi napakami
#----------------------------------------------------------------------------------------------


#linearizirani arrayi
noviiks=[1/t2[i] for i in range(len(t2))]
noviipsilon=[(1/((t2[i]/x[0])/(t2[i]+x[1]/x[0]))) for i in range(len(t2))]
e=[(3/seznam[i][1]**2) for i in range(dolzina)]
iks=[1/iks[i] for i in range(len(iks))]
ipsilon=[(1/ipsilon[i]) for i in range(len(ipsilon))]


#izrišemo
plt.plot(noviiks,noviipsilon,'b--',label='prilagoditvena funkcija')
plt.errorbar(iks,ipsilon,yerr=e,color='gray',fmt='o',label='meritve')
plt.title('Linearizirana prilagoditvena funkcija za dane meritve 2')
plt.ylabel('log(1/y)')
plt.xlabel('1/x')
plt.yscale('log')
plt.legend()

#neuspeli poskus
#e=[math.log10(3/seznam[i][1]**2) for i in range(dolzina)]
#plt.plot([1/t2[i] for i in range(len(t2))],[math.log10(1/((t2[i]/x[0])/(t2[i]+x[1]/x[0]))) for i in range(len(t2))],'b--',label='prilagoditvena funkcija')
#plt.errorbar([1/iks[i] for i in range(len(iks))],[math.log10(1/ipsilon[i]) for i in range(len(ipsilon))],yerr=e,color='gray',fmt='o',label='meritve')
#plt.legend()
#plt.ylabel('log(1/y)')
#plt.xlabel('1/x')
#plt.title('Linearizirana prilagoditvena funkcija za dane meritve 3')
#plt.show()



#%%
#%%
#%%
###########################################################################################################
#
#       2. NALOGA - vnesem podatke
#
###########################################################################################################

with open('thtg-xfp-thfp.dat') as f:
    vmesni=[l.strip().split(" ") for l in f]
    podaci=[[float(vmesni[j][i]) for i in range(len(vmesni[j])) if vmesni[j][i] != ""] for j in range(len(vmesni))]
dolzina=len(podaci)

#%%
#----------------------------------------------------------------------------------------------
#       2. NALOGA - definiram model(metoda njamanjših kvadratov) - LINEARNO (to sem opustil - brez normiranja)
#----------------------------------------------------------------------------------------------

th=np.array([podaci[i][1] for i in range(dolzina)])
x=np.array([podaci[i][2] for i in range(dolzina)])
y=np.array([podaci[i][0] for i in range(dolzina)])
    
def model(n):    
    def residuals(p, y, x, th, n):
        err = (y-pval(th,x,n,p))/0.0573
        return(err)
        
    def pval(th, x, n, p):
        return sum((x**j)*sum(p[i+j]*th**i for i in range(n)) for j in range(n))
    
    p0 = np.array([1.1 for i in range(2*n)])
    
    plsq = leastsq(residuals, p0, args=(y, th, x, n), maxfev=2000)
    ksi=sum((residuals(plsq[0],y[i],x[i],th[i],n)/0.0573)**2 for i in range(dolzina))/(dolzina-len(p0))
    
    return ksi
    
#%%

plt.plot([i+1 for i in range(20)],[model(i+1) for i in range(20)])
plt.yscale('log')
#%%
#----------------------------------------------------------------------------------------------
#       2. NALOGA - definiram model(metoda SVD) - LINEARNO (brez normiranja)
#----------------------------------------------------------------------------------------------
th=np.array([podaci[i][1] for i in range(dolzina)])
x=np.array([podaci[i][2] for i in range(dolzina)])
y=np.array([podaci[i][0] for i in range(dolzina)])

def pow(baza,potenca):
    rezultat=1
    for i in range(potenca):
        rezultat=rezultat*baza
    return rezultat

def A(n):
    return [[(pow(x[k],i))*(pow(th[k],j))/(0.0573) for i in range(n) for j in range(n)] for k in range(dolzina)]

def hitra(n):
    start=timer()
    celksi=[]
    for j in range(n):
        start1=timer()
        AS=A(j+1)
        sol,_,_,_= np.linalg.lstsq(AS,y)
        vrednost=np.dot(AS,sol)
        ksi=sum(((y[i]-vrednost[i])/(0.0573))**2 for i in range(dolzina))/(dolzina-len(sol))
        celksi.append(ksi)
        end1=timer()
        print('čas', j, '-tega koraka:', end1-start1,'\n \chi tega koraka pa je:',ksi)
    iksos=[i for i in range(n)]
    mini=np.argmin(celksi)
    plt.plot(iksos,celksi,label='minimum: ({0:.0f},  {1:.2f})'.format(mini,celksi[mini]))
    plt.title('$\chi^2$ pri linearni enečbi (SVD) - ni normirano')
    plt.xlabel('red polinoma n/2')
    plt.ylabel('$\chi^2$')
    plt.legend()
    end=timer()
    print('čas, ki ga je porabil algoritem za izračun te funkije je:', end-start)   

#%%
#----------------------------------------------------------------------------------------------
#       2. NALOGA - histogram residualov - LINEARNO (brez normiranja)
#----------------------------------------------------------------------------------------------


AS=A(8)
sol,_,_,_= np.linalg.lstsq(AS,y)
vrednost=np.dot(AS,sol)

kiklop=[]
for j in range(dolzina):
    kiklop.append((y[j]-vrednost[j])/(0.0573))
mu, std = norm.fit(kiklop)

plt.title('Histogram $ | \, \delta y \, | $ za linearno - brez normiranja (n=8)')
plt.hist(kiklop,500,normed=1, facecolor='green', alpha=0.75)
n, bins, patches = plt.hist(kiklop,500,normed=1, facecolor='green', alpha=0.75)
fitnormalne = mlab.normpdf( bins, mu, std)#prilagajamo normalno funkcijo
l = plt.plot(bins, fitnormalne, 'r--', linewidth=1,label='N( {0:.3f}, $\sigma$ = {1:.2f})'.format(mu,std))#prilagajamo normalno funkcijo
plt.ylabel('Število residualov $ | \, \delta y \, | $ z vrednostjo na danem intervalu')
plt.xlabel('Velikost residualov $ | \, \delta y \, | $')
plt.legend()
plt.grid(True)



#%%
#----------------------------------------------------------------------------------------------
#       2. NALOGA - definiram model(metoda SVD) - LINEARNO (normirano)
#----------------------------------------------------------------------------------------------

#spremenimo interval kjer s epojavljajo x in th zato, da jih lahko vstavimo v čebiševe polinome
th=np.array([podaci[i][1]/164.102 for i in range(dolzina)])
x=np.array([podaci[i][2]/164.102 for i in range(dolzina)])
y=np.array([podaci[i][0]/164.102 for i in range(dolzina)])

def pow(baza,potenca):
    rezultat=1
    for i in range(potenca):
        rezultat=rezultat*baza
    return rezultat

def A(n):
    return [[(pow(x[k],i))*(pow(th[k],j))/(0.0573/164.102) for i in range(n) for j in range(n)] for k in range(dolzina)]

def hitra(n):
    start=timer()
    celksi=[]
    for j in range(n):
        start1=timer()
        AS=A(j+1)
        sol,_,_,_= np.linalg.lstsq(AS,y)
        vrednost=np.dot(AS,sol)
        ksi=sum(((y[i]-vrednost[i])/(0.0573/164.102))**2 for i in range(dolzina))/(dolzina-len(sol))
        celksi.append(ksi)
        end1=timer()
        print('čas', j, '-tega koraka:', end1-start1,'\n \chi tega koraka pa je:',ksi)
    iksos=[i for i in range(n)]
    mini=np.argmin(celksi)
    plt.plot(iksos,celksi,label='minimum: ({0:.0f},  {1:.2f})'.format(mini,celksi[mini]))
    plt.title('$\chi^2$ pri linearni enečbi (SVD)')
    plt.xlabel('red polinoma n/2')
    plt.ylabel('$\chi^2$')
    plt.legend()
    end=timer()
    print('čas, ki ga je porabil algoritem za izračun te funkije je:', end-start)   

#%%
#----------------------------------------------------------------------------------------------
#       2. NALOGA - histogram residualov - LINEARNO (normirano)
#----------------------------------------------------------------------------------------------


AS=A(15)
sol,_,_,_= np.linalg.lstsq(AS,y)
vrednost=np.dot(AS,sol)

kiklop=[]
for j in range(dolzina):
    kiklop.append((y[j]-vrednost[j])/(0.0573/164.102))
mu, std = norm.fit(kiklop)

plt.title('Histogram $ | \, \delta y \, | $ za linearno (n=24)')
plt.hist(kiklop,500,normed=1, facecolor='green', alpha=0.75)
n, bins, patches = plt.hist(kiklop,500,normed=1, facecolor='green', alpha=0.75)
fitnormalne = mlab.normpdf( bins, mu, std)#prilagajamo normalno funkcijo
l = plt.plot(bins, fitnormalne, 'r--', linewidth=1,label='N( {0:.3f}, $\sigma$ = {1:.2f})'.format(mu,std))#prilagajamo normalno funkcijo
plt.ylabel('Število residualov $ | \, \delta y \, | $ z vrednostjo na danem intervalu')
plt.xlabel('Velikost residualov $ | \, \delta y \, | $')
plt.legend()
plt.grid(True)

#%%
#----------------------------------------------------------------------------------------------
#       2. NALOGA - graf meritev in fit - POLINOMI ČEBIŠEVA
#----------------------------------------------------------------------------------------------

iksos=np.arange(-1, 1, 0.05)
thetos=np.arange(-1, 1, 0.05)
def splosno(n):
    return [[(pow(iksos[k],i))*(pow(thetos[k],j))/(0.0573/164.102) for i in range(n) for j in range(n)] for k in range(len(iksos))]
fitos=np.dot(splosno(15),sol)



plt.title('prilagoditvena funkcija in podatki za linearno (n=15)')

t=[i for i in range(dolzina)]
l = plt.plot(x, y, 'r,', linewidth=1)
plt.plot(thetos,fitos,'b:',linewidth=0.5)
plt.ylabel('meritve na tarči (normirano)')
plt.xlabel('kot na detektorju (normirano)')
plt.legend()
plt.grid(True)

#%%
#----------------------------------------------------------------------------------------------
#       2. NALOGA - definiram model(metoda SVD) - POLINOMI ČEBIŠEVA
#----------------------------------------------------------------------------------------------

#spremenimo interval kjer s epojavljajo x in th zato, da jih lahko vstavimo v čebiševe polinome
th=np.array([podaci[i][1]/164.102 for i in range(dolzina)])
x=np.array([podaci[i][2]/164.102 for i in range(dolzina)])
y=np.array([podaci[i][0]/164.102 for i in range(dolzina)])



def A(n):
    return chebvander2d(x,th,[n,n])/(0.0573/164.102)



def chebi(n):
    start=timer()
    celksi=[]
    for j in range(n):
        start1=timer()
        AS=A(j)
        sol,c,g,l= np.linalg.lstsq(AS,y)
        vrednost=np.dot(AS,sol)
        ksi=sum(((y[i]-vrednost[i])/(0.0573/164.102))**2 for i in range(dolzina))/(dolzina-len(sol))
        celksi.append(ksi)
        end1=timer()
        print('čas', j, '-tega koraka:', end1-start1,'\n \chi tega koraka pa je:',ksi)
    iksos=[i for i in range(n)]
    mini=np.argmin(celksi)
    plt.plot(iksos,celksi,label='minimum: ({0:.0f},  {1:.2f})'.format(mini,celksi[mini]))
    plt.title('$\chi^2$ pri vrsti Čebiševih polinomov (SVD)')
    plt.xlabel('red polinoma n/2')
    plt.ylabel('$\chi^2$')
    plt.legend()
    end=timer()
    print('čas, ki ga je porabil algoritem za izračun te funkije je:', end-start)
    
#%%
#----------------------------------------------------------------------------------------------
#       2. NALOGA - histogram residualov - POLINOMI ČEBIŠEVA
#----------------------------------------------------------------------------------------------


AS=A(20)
sol,c,g,l= np.linalg.lstsq(AS,y)
vrednost=np.dot(AS,sol)

kiklop=[]
for j in range(dolzina):
    kiklop.append((y[j]-vrednost[j])/(0.0573/164.102))
mu, std = norm.fit(kiklop)

plt.title('Histogram $ | \, \delta y \, | $ za Čebiševe polinome (n=29)')
n, bins, patches = plt.hist(kiklop,500,normed=1, facecolor='green', alpha=0.75)
fitnormalne = mlab.normpdf( bins, mu, std)#prilagajamo normalno funkcijo
l = plt.plot(bins, fitnormalne, 'r--', linewidth=1,label='N( {0:.3f}, $\sigma$ = {1:.2f})'.format(mu,std))#prilagajamo normalno funkcijo
plt.ylabel('Število residualov $ | \, \delta y \, | $ z vrednostjo na danem intervalu')
plt.xlabel('Velikost residualov $ | \, \delta y \, | $')
plt.legend()
plt.grid(True)

#%%
#----------------------------------------------------------------------------------------------
#       2. NALOGA - graf residualov - POLINOMI ČEBIŠEVA
#----------------------------------------------------------------------------------------------

plt.title('Graf residualov  $ | \, \delta y \, | $ za Čebiševe polinome (n=29)')

t=[i for i in range(dolzina)]
l = plt.plot(t, kiklop, 'r:', linewidth=1)
plt.ylabel('Velikost residualov $ | \, \delta y \, | $')
plt.xlabel('število meritve')
plt.legend()
plt.grid(True)

#%%
#----------------------------------------------------------------------------------------------
#       2. NALOGA - graf meritev in fit - POLINOMI ČEBIŠEVA
#----------------------------------------------------------------------------------------------
iksos=np.arange(-1, 1, 0.05)
thetos=np.arange(-1, 1, 0.05)
def splosno(n):
    return chebvander2d(iksos,thetos,[n,n])/(0.0573/164.102)
fitos=np.dot(splosno(20),sol)



plt.title('prilagoditvena funkcija in podatki za Čebiševe polinome (n=29)')

t=[i for i in range(dolzina)]
l = plt.plot(x, y, 'r,', linewidth=1)
plt.plot(thetos,fitos,'b:',linewidth=0.5)
plt.ylabel('meritve na tarči (normirano)')
plt.xlabel('kot na detektorju (normirano)')
plt.legend()
plt.grid(True)

#%%
#----------------------------------------------------------------------------------------------
#       2. NALOGA - definiram model(metoda SVD) - POLINOMI ČEBIŠEVA na roke!
#----------------------------------------------------------------------------------------------


#def čebišev(n):
#    start=timer()
#    celksi=[]
#    for j in range(n):
#        start1=timer()
#        AS=A(j)
        U, s, V = np.linalg.svd(AS, full_matrices=True)
        sol=[np.multiply(np.divide(np.dot(np.transpose(U)[i],y),s[i]),V[i]) for i in range(len(V))] # edina sprememba je tu, kjer se barva ne spremeni
        sol=np.dot(sol,[1 for i in range(len(sol))])
        vrednost=np.dot(AS,sol)
        ksi=sum(((y[i]-vrednost[i])/0.0573)**2 for i in range(dolzina))/(dolzina-len(sol))
#        celksi.append(ksi)
#        end1=timer()
#        print('čas', len(sol), '-tega koraka:', end1-start1,'\n \chi tega koraka pa je:',ksi)
#    iksos=[i for i in range(n)]
#    plt.plot(iksos,celksi)
#    plt.title('$\chi^2$ pri vrsti Čebiševih polinomov (SVD)')
#    plt.xlabel('red polinoma n/2')
#    plt.ylabel('$\chi^2/(m-n)$')
#    end=timer()
#    print('čas, ki ga je porabil algoritem za izračun te funkije je:', end-start)




#%%
#----------------------------------------------------------------------------------------------
#       2. NALOGA - definiram model(metoda SVD) - LEGENDROVI POLINOMI
#----------------------------------------------------------------------------------------------

#spremenimo interval kjer s epojavljajo x in th zato, da jih lahko vstavimo v čebiševe polinome
th=np.array([podaci[i][1]/164.102 for i in range(dolzina)])
x=np.array([podaci[i][2]/164.102 for i in range(dolzina)])
y=np.array([podaci[i][0]/164.102 for i in range(dolzina)])



def A(n):
    return legvander2d(x,th,[n,n])/(0.0573/164.102)



def legendi(n):
    start=timer()
    celksi=[]
    for j in range(n):
        start1=timer()
        AS=A(j)
        sol,c,g,l= np.linalg.lstsq(AS,y)
        vrednost=np.dot(AS,sol)
        ksi=sum(((y[i]-vrednost[i])/(0.0573/164.102))**2 for i in range(dolzina))/(dolzina-len(sol))
        celksi.append(ksi)
        end1=timer()
        print('čas', j, '-tega koraka:', end1-start1,'\n \chi tega koraka pa je:',ksi)
    iksos=[i for i in range(n)]
    mini=np.argmin(celksi)
    plt.plot(iksos,celksi,label='minimum: ({0:.0f},  {1:.2f})'.format(mini,celksi[mini]))
    plt.title('$\chi^2$ pri vrsti legendrovih polinomov (SVD)')
    plt.xlabel('red polinoma n/2')
    plt.ylabel('$\chi^2$')
    plt.legend()
    end=timer()
    print('čas, ki ga je porabil algoritem za izračun te funkije je:', end-start)

#%%
#----------------------------------------------------------------------------------------------
#       2. NALOGA - histogram residualov - LEGENDROVI POLINOMI
#----------------------------------------------------------------------------------------------


AS=A(29)
sol,_,_,_= np.linalg.lstsq(AS,y)
vrednost=np.dot(AS,sol)

kiklop=[]
for j in range(dolzina):
    kiklop.append((y[j]-vrednost[j])/(0.0573/164.102))
mu, std = norm.fit(kiklop)

plt.title('Histogram $ | \, \delta y \, | $ za Legendrove polinome (n=29)')
n, bins, patches = plt.hist(kiklop,500,normed=1, facecolor='green', alpha=0.75)
fitnormalne = mlab.normpdf( bins, mu, std)#prilagajamo normalno funkcijo
l = plt.plot(bins, fitnormalne, 'r--', linewidth=1,label='N( {0:.3f}, $\sigma$ = {1:.2f})'.format(mu,std))#prilagajamo normalno funkcijo
plt.ylabel('Število residualov $ | \, \delta y \, | $ z vrednostjo na danem intervalu')
plt.xlabel('Velikost residualov $ | \, \delta y \, | $')
plt.legend()
plt.grid(True)

#%%
#----------------------------------------------------------------------------------------------
#       2. NALOGA - definiram model(metoda SVD) - HERMITOVI POLINOMI
#----------------------------------------------------------------------------------------------

#spremenimo interval kjer s epojavljajo x in th zato, da jih lahko vstavimo v čebiševe polinome
th=np.array([podaci[i][1]/164.102 for i in range(dolzina)])
x=np.array([podaci[i][2]/164.102 for i in range(dolzina)])
y=np.array([podaci[i][0]/164.102 for i in range(dolzina)])



def A(n):
    return hermvander2d(x,th,[n,n])/(0.0573/164.102)



def hermit(n):
    start=timer()
    celksi=[]
    for j in range(n):
        start1=timer()
        AS=A(j)
        sol,c,g,l= np.linalg.lstsq(AS,y)
        vrednost=np.dot(AS,sol)
        ksi=sum(((y[i]-vrednost[i])/(0.0573/164.102))**2 for i in range(dolzina))/(dolzina-len(sol))
        celksi.append(ksi)
        end1=timer()
        print('čas', j, '-tega koraka:', end1-start1,'\n \chi tega koraka pa je:',ksi)
    iksos=[i for i in range(n)]
    mini=np.argmin(celksi)
    plt.plot(iksos,celksi,label='minimum: ({0:.0f},  {1:.2f})'.format(mini,celksi[mini]))
    plt.title('$\chi^2$ pri vrsti Hermitovih polinomov (SVD)')
    plt.xlabel('red polinoma n/2')
    plt.ylabel('$\chi^2$')
    plt.legend()
    end=timer()
    print('čas, ki ga je porabil algoritem za izračun te funkije je:', end-start)

#%%
#----------------------------------------------------------------------------------------------
#       2. NALOGA - histogram residualov - HERMITOVI POLINOMI
#----------------------------------------------------------------------------------------------


AS=A(11)
sol,_,_,_= np.linalg.lstsq(AS,y)
vrednost=np.dot(AS,sol)

kiklop=[]
for j in range(dolzina):
    kiklop.append((y[j]-vrednost[j])/(0.0573/164.102))
mu, std = norm.fit(kiklop)

plt.title('Histogram $ | \, \delta y \, | $ za Hermitove polinome (n=11)')
n, bins, patches = plt.hist(kiklop,500,normed=1, facecolor='green', alpha=0.75)
fitnormalne = mlab.normpdf( bins, mu, std)#prilagajamo normalno funkcijo
l = plt.plot(bins, fitnormalne, 'r--', linewidth=1,label='N( {0:.3f}, $\sigma$ = {1:.2f})'.format(mu,std))#prilagajamo normalno funkcijo
plt.ylabel('Število residualov $ | \, \delta y \, | $ z vrednostjo na danem intervalu')
plt.xlabel('Velikost residualov $ | \, \delta y \, | $')
plt.legend()
plt.grid(True)


#%%
#%%
#%%
###########################################################################################################
#
#       3. NALOGA - vnesem podatke - RENDGENSKI ABSORBCIJSKI ROBOVI
#
###########################################################################################################
# This file contains normalized mu(E) from:
# listi_cell_wall_epi.xmu (column 2)
# listi_cell_wall_mezo.xmu (column 3)
# CdSO4_GSH.xmu (column 4)
# CdSO4_pectin_tr.xmu (column 5)
with open('zadnjanaloga.norm') as f:
    vmesni=[l.strip().split(" ") for l in f]
    podaci=[[float(vmesni[j][i]) for i in range(len(vmesni[j])) if vmesni[j][i] != ""] for j in range(len(vmesni)) if j > 7]
dolzina=len(podaci)

E= [podaci[i][0] for i in range(dolzina)]
A = [[podaci[i][3], podaci[i][4]] for i in range(dolzina)]
yepi = [podaci[i][1] for i in range(dolzina)]
ymezo = [podaci[i][2] for i in range(dolzina)]

solepi,_,_,_= np.linalg.lstsq(A,yepi)
vrednostepi=np.dot(A,solepi)
ksiepi=sum((yepi[i]-vrednostepi[i])**2 for i in range(dolzina))/(dolzina-len(solepi))

solmezo,_,_,_= np.linalg.lstsq(A,ymezo)
vrednostmezo=np.dot(A,solmezo)
ksimezo=sum((ymezo[i]-vrednostmezo[i])**2 for i in range(dolzina))/(dolzina-len(solmezo))

print(solepi,'in napaka:',ksiepi,'\n\n',solmezo,'in napaka:',ksimezo)

#%%
#----------------------------------------------------------------------------------------------
#       3. NALOGA - rešitev za prvo rastlino - RENDGENSKI ABSORBCIJSKI ROBOVI
#----------------------------------------------------------------------------------------------


fitepi=np.dot(A,solepi)



plt.title('prilagoditvena funkcija in podatki za prvo rastlino')

plt.plot(E,fitepi,'y.',linewidth=1,label='delež GSH:{0:.2f} in pektina:{1:.2f}'.format(solepi[0],solepi[1]))
plt.plot(E,fitepi,'g.',linewidth=1,label='reducirani $\chi^2$:{0:.4f}'.format(ksiepi))
plt.plot(E,[podaci[i][3]*solepi[0] for i in range(dolzina)],'k.',linewidth=1,label='delež GSH')
plt.plot(E,[podaci[i][4]*solepi[1] for i in range(dolzina)],'c.',linewidth=1,label='delež pektina')
plt.plot(E, yepi, 'b.', linewidth=1,label='meritve')
plt.plot(E,fitepi,'r.',linewidth=1,label='naša rešitev')
plt.ylabel('absorbcija')
plt.xlabel('energija')
plt.legend(loc=2)
plt.grid(True)

#%%
#----------------------------------------------------------------------------------------------
#       3. NALOGA - rešitev za drugo rastlino - RENDGENSKI ABSORBCIJSKI ROBOVI
#----------------------------------------------------------------------------------------------


fitmezo=np.dot(A,solmezo)



plt.title('prilagoditvena funkcija in podatki za drugo rastlino')

plt.plot(E,fitmezo,'y.',linewidth=1,label='delež GSH:{0:.2f} in pektina:{1:.2f}'.format(solmezo[0],solmezo[1]))
plt.plot(E,fitmezo,'g.',linewidth=1,label='reducirani $\chi^2$:{0:.4f}'.format(ksimezo))
plt.plot(E,[podaci[i][3]*solmezo[0] for i in range(dolzina)],'k.',linewidth=1,label='delež GSH')
plt.plot(E,[podaci[i][4]*solmezo[1] for i in range(dolzina)],'c.',linewidth=1,label='delež pektina')
plt.plot(E, ymezo, 'b.', linewidth=1,label='meritve')
plt.plot(E,fitmezo,'r.',linewidth=1,label='naša rešitev')
plt.ylabel('absorbcija')
plt.xlabel('energija')
plt.legend(loc=2)
plt.grid(True)


#%%
#----------------------------------------------------------------------------------------------
#       3. NALOGA - vnesem podatke - RENDGENSKI ABSORBCIJSKI ROBOVI (nova napaka)
#----------------------------------------------------------------------------------------------
# This file contains normalized mu(E) from:
# listi_cell_wall_epi.xmu (column 2)
# listi_cell_wall_mezo.xmu (column 3)
# CdSO4_GSH.xmu (column 4)
# CdSO4_pectin_tr.xmu (column 5)
with open('zadnjanaloga.norm') as f:
    vmesni=[l.strip().split(" ") for l in f]
    podaci=[[float(vmesni[j][i]) for i in range(len(vmesni[j])) if vmesni[j][i] != ""] for j in range(len(vmesni)) if j > 7]
dolzina=len(podaci)



E= [podaci[i][0] for i in range(dolzina)]
Aepi = [[podaci[i][3]/math.sqrt(0.000177294784184), podaci[i][4]/math.sqrt(0.000177294784184)] for i in range(dolzina)]
Amezo = [[podaci[i][3]/math.sqrt(0.000355617427079), podaci[i][4]/math.sqrt(0.000355617427079)] for i in range(dolzina)]
yepi = [podaci[i][1]/math.sqrt(0.000177294784184) for i in range(dolzina)]
ymezo = [podaci[i][2]/math.sqrt(0.000355617427079) for i in range(dolzina)]

solepi,_,_,_= np.linalg.lstsq(Aepi,yepi)
vrednostepi=np.dot(Aepi,solepi)
ksiepi=sum((yepi[i]-vrednostepi[i])**2 for i in range(dolzina))/(dolzina-len(solepi))

solmezo,_,_,_= np.linalg.lstsq(Amezo,ymezo)
vrednostmezo=np.dot(Amezo,solmezo)
ksimezo=sum((ymezo[i]-vrednostmezo[i])**2 for i in range(dolzina))/(dolzina-len(solmezo))

print(solepi,'in napaka:',ksiepi,'\n\n',solmezo,'in napaka:',ksimezo)

#%%
#----------------------------------------------------------------------------------------------
#       3. NALOGA - rešitev za prvo rastlino - RENDGENSKI ABSORBCIJSKI ROBOVI (nova napaka)
#----------------------------------------------------------------------------------------------


fitepi=np.dot(Aepi,solepi)



plt.title('prilagoditvena funkcija in podatki za prvo rastlino')

plt.plot(E,fitepi,'y.',linewidth=1,label='delež GSH:{0:.2f} in pektina:{1:.2f}'.format(solepi[0],solepi[1]))
plt.plot(E,fitepi,'g.',linewidth=1,label='reducirani $\chi^2$:{0:.4f}'.format(ksiepi))
plt.plot(E,[podaci[i][3]*solepi[0]/math.sqrt(0.000177294784184) for i in range(dolzina)],'k.',linewidth=1,label='delež GSH')
plt.plot(E,[podaci[i][4]*solepi[1]/math.sqrt(0.000177294784184) for i in range(dolzina)],'c.',linewidth=1,label='delež pektina')
plt.plot(E, yepi, 'b.', linewidth=1,label='meritve')
plt.plot(E,fitepi,'r.',linewidth=1,label='naša rešitev')
plt.ylabel('absorbcija')
plt.xlabel('energija')
plt.legend(loc=2)
plt.grid(True)

#%%
#----------------------------------------------------------------------------------------------
#       3. NALOGA - rešitev za drugo rastlino - RENDGENSKI ABSORBCIJSKI ROBOVI (nova napaka)
#----------------------------------------------------------------------------------------------


fitmezo=np.dot(Amezo,solmezo)



plt.title('prilagoditvena funkcija in podatki za drugo rastlino')

plt.plot(E,fitmezo,'y.',linewidth=1,label='delež GSH:{0:.2f} in pektina:{1:.2f}'.format(solmezo[0],solmezo[1]))
plt.plot(E,fitmezo,'g.',linewidth=1,label='reducirani $\chi^2$:{0:.4f}'.format(ksimezo))
plt.plot(E,[podaci[i][3]*solmezo[0]/math.sqrt(0.000355617427079) for i in range(dolzina)],'k.',linewidth=1,label='delež GSH')
plt.plot(E,[podaci[i][4]*solmezo[1]/math.sqrt(0.000355617427079) for i in range(dolzina)],'c.',linewidth=1,label='delež pektina')
plt.plot(E, ymezo, 'b.', linewidth=1,label='meritve')
plt.plot(E,fitmezo,'r.',linewidth=1,label='naša rešitev')
plt.ylabel('absorbcija')
plt.xlabel('energija')
plt.legend(loc=2)
plt.grid(True)

#%%
















