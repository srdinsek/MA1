import math
import numpy as np
import matplotlib.pyplot as plt
from sympy import simplify
from random import random
from scipy.optimize import minimize
from scipy.optimize import leastsq
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from numpy import linalg as LA
from scipy.linalg import solve




import seaborn as sns
sns.set_palette(sns.color_palette("cool", 10))


#PRVA NALOGA
if 0==1:
    with open('farmakoloski1.txt') as f:
        podaci=[l.strip().split("\t") for l in f]
    seznam=[[float(podaci[i][j]) for j in range(len(podaci[i]))] for i in range(len(podaci)) if i >> 0]
    dolzina=len(seznam)
    fi=[[1/seznam[i][0] if j==1 else 1 for j in [0,1]] for i in range(dolzina)]
    x=np.array([seznam[i][0] for i in range(dolzina)])
    y=np.array([seznam[i][1] for i in range(dolzina)])
    
    
    def residuals(p, y, x):
        err = (y-pval(x,p))/3
        return(err)
    
    def pval(x, p):
        return(p[0]*x/(x+p[1]))
    
    
    p0 = np.array([33 , 6.5, 1.3])
    e=[3 for i in range(dolzina)]
    
    
    
##    plt.errorbar(x,y,yerr=e,color='gray',fmt='o',label='meritve')
    plsq = leastsq(residuals, p0, args=(y, x), maxfev=2000)
    zanimivo2=[np.abs(residuals(plsq[0],y[i],x[i])) for i in range(dolzina)]
    zanimivod2=[i for i in range(len(zanimivo))]
###    slika=np.abs(np.fft.rfft([residuals(plsq[0],y[i],x[i]) for i in range(dolzina)]))
    print(plsq)
    plt.plot(zanimivod,zanimivo,color='red',marker='.', label='p različen od 1')
    plt.fill_between(zanimivod,zanimivo,interpolate=True,color='green',alpha=0.1)
    plt.plot(zanimivod2,zanimivo2,color='blue',marker='.', label='p enak 1')
    plt.fill_between(zanimivod2,zanimivo2,interpolate=True,color='brown',alpha=0.1)
    #ksi=sum((seznam[i][1]-(plsq[0][0]*seznam[i][0]**plsq[0][2]/(seznam[i][0]**plsq[0][2]+plsq[0][1]**plsq[0][2])))**2/9 for i in range(dolzina))
    ksi=sum((residuals(plsq[0],y[i],x[i]))**2 for i in range(dolzina))
##    plt.plot(x,y,'ro',label="$\chi^2$: {0:.2f}".format(ksi))
    t=np.linspace(0,1200,1000)
###    plt.bar([0,1,2,3,4],slika,color='red', alpha=0.5)
##    plt.plot(t,(t*109.299529422)/(t+26.604825283),'r--',label='linearna prilagoditvena funkcija')
##    plt.plot(t,pval(t,plsq[0]),'b--',label='nelinearna prilagoditvena funkcija')
    plt.title('Absolutne vrednosti odstopanja')
    plt.xlabel('korak')
    plt.ylabel('odstopanje $|\delta\,\, y|$')
    plt.legend(loc=1)
    

   

    
    #linearizirano ampak so napake napačne, niso lineariziriane!
    #e=[3 for i in range(dolzina)]
    #plt.errorbar([1/x[i] for i in range(dolzina)],[1/y[i] for i in range(dolzina)],yerr=[1/e[i] for i in range(dolzina)],color='gray',fmt='o',label='meritve')
    #ksi=sum((residuals(plsq[0],y[i],x[i]))**2/9 for i in range(dolzina))
    #plt.plot([1/x[i] for i in range(dolzina)],[1/y[i] for i in range(dolzina)],'ro',label="$\chi^2$: {0:.2f}".format(ksi))
    #t=np.linspace(-10,1000,200000)
    #resitu=pval(t,plsq[0])
    #plt.plot([1/t[i] for i in range(len(t))],[1/resitu[i] for i in range(len(t))],'b--',label='prilagoditvena funkcija')
    #plt.title('Nelinearna prilagoditvena funkcija za dane meritve')
    #plt.xlabel('1/x')
    #plt.ylabel('1/y')
    #plt.yscale('log')
    #plt.legend(loc=4)
    
    
    ###########################################################################################################################################
    
    A=[[sum((seznam[i][1]**2)*fi[i][j]*fi[i][k]/3 for i in range(dolzina)) for j in range(len(fi[1]))] for k in range(len(fi[1]))]
    b=[sum((1/seznam[i][1])*fi[i][k]*seznam[i][1]**2/3 for i in range(dolzina)) for k in [0,1]]
    x=solve(A,b)
    
    print('\t \t')
    print('Imamo matriko:\t', seznam)
    print('In rešitev:\t\t', b)
    print('Ko tako enačbo rešimo, dobimo rešitev:\t', x, 'pri čemer, so konstante a =', x[1]/x[0], 'in y0 =', 1/x[0])

#DRUGA NALOGA
if 0==0:

    with open('ledvice.dat') as f:
        podaci=[l.strip().split("\t") for l in f]
    seznam=[[float(podaci[i][j]) for j in range(len(podaci[i]))] for i in range(len(podaci)) if i >> 0] #prvi stolpec je čas, drugi stolpec je število sunkov na detektorju
    dolzina=len(seznam)
    print(dolzina)
    x=np.array([seznam[i][0] for i in range(dolzina)])
    y=np.array([seznam[i][1] for i in range(dolzina)])
    print(seznam,'\n\n',x,'\n\n',y)
    
    def residuals2(c, y, x):
        err = (y-pval2(x,c))
        return(err)
    
    def pval2(x, c):
        return c[0]*math.exp(1)**(-x*c[1])+c[2]*math.exp(1)**(-x*c[3])+c[4]*math.exp(1)**(-x*c[5])+c[6]*math.exp(1)**(-x*c[7])
    
    p0 = np.array([1,0.1,1,0.2,1,0.5,1,0.2])
    
    plsq2 = leastsq(residuals2, p0, args=(y, x), maxfev=2000)
    
#    print(plsq2[0])
    ksi2=sum((residuals2(plsq2[0],y[i],x[i]))**2 for i in range(dolzina))/(len(x)-len(p0))
    print(ksi2)
    print(2)
    
##    slika=np.abs(np.fft.rfft([np.abs(residuals2(plsq2[0],y[i],x[i])) for i in range(dolzina)]))
##    plt.bar([i for i in range(len(slika))],slika,color='red', alpha=0.5)

    if 4==0:    
        zanimivo5=[np.abs(residuals2(plsq2[0],y[i],x[i])) for i in range(dolzina)]
        zanimivod5=[i for i in range(len(zanimivo))]
        
        plt.plot(zanimivod,zanimivo,color='red',marker='.',label='navaden')
        plt.fill_between(zanimivod,zanimivo,interpolate=True,color='green',alpha=0.1)
        
        plt.plot(zanimivod2,zanimivo2,color='blue',marker='.', label='navaden + konstanta')
        plt.fill_between(zanimivod2,zanimivo2,interpolate=True,color='brown',alpha=0.1)
        
        plt.plot(zanimivod3,zanimivo3,color='yellow',marker='.', label='keronsko')
        plt.fill_between(zanimivod3,zanimivo3,interpolate=True,color='red',alpha=0.1)
        
        plt.plot(zanimivod4,zanimivo4,color='green',marker='.', label='korensko + konstanta')
        plt.fill_between(zanimivod4,zanimivo4,interpolate=True,color='gold',alpha=0.1)
        
    #    plt.plot(zanimivod5,zanimivo5,color='black',marker='.', label='dva kompartmenta')
    #    plt.fill_between(zanimivod5,zanimivo5,interpolate=True,color='pink',alpha=0.1)
    

    
    
    
    
    e=[0 for i in range(dolzina)]
    plt.errorbar(x,y,yerr=e,color='gray',fmt='o',label='meritve')
    plt.plot(x,y,'ro',label="$\chi^2/(m-n)$: {0:.2f}".format(ksi2))
#    plt.plot([i for i in range(14)],resitev)
    t=np.linspace(0,2200,1000)
    tau=np.linspace(0,500,1000)
    plt.plot(t,pval2(t,plsq2[0]),'b--',label='štirje eksponentni členi')
#    plt.plot(t,plsq2[0][0]*math.exp(1)**(-t*plsq2[0][1]),'r--',label='$Ae^{-t*\lambda_{A}}$')
#    plt.plot(tau,plsq2[0][2]*math.exp(1)**(-tau*plsq2[0][3]),'g--',label='$Be^{-t*\lambda_{B}}$')
    #plt.plot(t,resitev_koren[0][0]*math.exp(1)**(-t**(1/2)*resitev_koren[0][1]),'r--',label='$c[0]e^{-\sqrt{t}*\lambda}$')
    #plt.plot(t,resitev_konstanta[0][0]*math.exp(1)**(-t*resitev_konstanta[0][1])+resitev_konstanta[0][2],'g--',label='$c[0]e^{-t*\lambda}+konstanta$')
    #plt.plot(t,resitev_brez[0][0]*math.exp(1)**(-t*resitev_brez[0][1]),'y--',label='$c[0]e^{-t*\lambda}$')
    
    plt.title('prilagoditvena funkcija za štiri čene')
    plt.xlabel('napetost U')
    plt.ylabel('logaritem tok I')
    plt.legend(loc=1)
    plt.yscale('symlog')




if 0==1:
    
    with open('korozija.txt') as f:
        vmesni=[l.strip().split(" ") for l in f]
        seznam=[[float(vmesni[j][i]) for i in range(len(vmesni[j]))] for j in range(len(vmesni))] #prvi stolpec je čas, drugi stolpec je število sunkov na detektorju
    print(podaci)
    dolzina=len(seznam)
    print(dolzina)
    y=np.array([seznam[i][0] for i in range(dolzina)])
    x=np.array([seznam[i][1] for i in range(dolzina)])
    print(x,'\n\n\n',y)
    
    def residuals3(I, y, x):
        err = (y-pval3(x,I))
        return(err)
        
    def pval3(x, I):
        return I[0]*x+I[1]*x**2+I[2]*x**3#+I[3]*x**4+I[4]*x**5+I[5]*x**6
    
    #I[0]*(math.exp(1)**((x-I[3])/I[1])-math.exp(1)**(-(x-I[3])/I[2]))
        
    p0 = np.array([100 , 10,5])
    plsq3 = leastsq(residuals3, p0, args=(y, x), maxfev=2000)
    print('\n\nkonstante so: \t\t\t\t',plsq3[0])
    ksi3=sum((residuals3(plsq3[0],y[i],x[i]))**2 for i in range(dolzina))/(dolzina-len(p0))
    print('$\chi^2$ je:\t\t\t',ksi3)
        
w###    slika=np.abs(np.fft.rfft([residuals3(plsq3[0],y[i],x[i]) for i in range(dolzina)]))
###    plt.bar([i for i in range(len(slika))],slika,color='red', alpha=0.5)
    
    
    
    zanimivo2=[np.abs(residuals3(plsq3[0],y[i],x[i])) for i in range(dolzina)]
    zanimivod2=[i for i in range(len(zanimivo2))]
    
    
####    plt.plot(zanimivod,zanimivo,color='red',marker='.', label='brez zamika')
####    plt.fill_between(zanimivod,zanimivo,interpolate=True,color='green',alpha=0.1)
####    plt.plot(zanimivod2,zanimivo2,color='blue',marker='.', label='z zamikom')
####    plt.fill_between(zanimivod2,zanimivo2,interpolate=True,color='brown',alpha=0.1)
    
    
    e=[0 for i in range(dolzina)]
    plt.errorbar(x,y,yerr=e,color='gray',fmt='o',label='meritve')
    plt.plot(x,y,'bo',label="$\chi^2/(m-n)$: {0:.2f}".format(ksi3))
    t=np.linspace(-0.01,0.006,1000)
    plt.plot(t,pval3(t,plsq3[0]),'b--',label='$A U+B U^2 + C U^3$')
 #   plt.plot(t,plsq3[0][0]*(math.exp(1)**((t-plsq3[0][3])/plsq3[0][1])),'r--',label='$I_{0}*exp((U-U_{0})/U_{a})$')
 #   plt.plot(t,-plsq3[0][0]*math.exp(1)**(-(t-plsq3[0][3])/plsq3[0][2]),'g--',label='$-I_{0}*exp(-(U-U_{0})/U_{c}) $ ')
    plt.title('Prilagoditvena funkcija')
    plt.xlabel('napetost U')
    plt.ylabel('tok I')
    plt.legend(loc=2)
    plt.show()

