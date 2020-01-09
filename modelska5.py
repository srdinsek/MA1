import math
import numpy as np
import matplotlib.pyplot as plt
from sympy import simplify
from random import random
from scipy.optimize import minimize
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm


import seaborn as sns
sns.set_palette(sns.color_palette("cool", 10))

#def binarna(y,t,p,q,r):
#    a,ao,c,b=y
#    dydt=[-p*a**2+q*a*ao,p*a**2-q*a*ao-r*ao,r*ao,r*ao]
#    return dydt


#def binarnaapprox(y,t,p,q,r):
#    a,c,b=y
#    dydt=[-p*a**2+q*a*p*a**2/(q*a+r),r*p*a**2/(q*a+r),r*p*a**2/(q*a+r)]
#    return dydt


def drugareakcija(y,t,p,q,r,s,T):
    u,v,z,x,y=y
    dydt=[s*x*y-u*z*r,q*z**2-v*p-T*v*y,s*x*y-r*z*u+T*y*v-q*z**2+v*p,r*z*u-s*x*y+T*y*v,r*z*u-s*x*y-T*y*v]
    return dydt


#def ure(y,t,p1,p2):
#    x,y,w,a=y
#    dydt=[-p1*x*y,2*p2*a*w-2*p1*x*y,p1*x*y-p2*a*w,-2*p2*a*w]
#    return dydt




if 0==0:
    
    for i in np.linspace(0,10,1):
        y0=[10,15,22,0.01,20]
        t=np.linspace(0,0.05,20000)
                
        sol = odeint(drugareakcija, y0, t, args=(2,2,3,4,5))
        q=i
        p=3000
        r=10
     #   ao=[p*a**2/(q*a+r) for a in sol[:,0]]        
        neki=[(((sol[i+1,3]-sol[i,3])/(0.05/20000))/sol[i,0])**2/sol[i,1] for i in range(len(sol[:,3])-1)]
        neki.append(0)
#        plt.plot(t,sol[:,0]+sol[:,3]+sol[:,1]+sol[:,2], color='black',label='vsota')
#        plt.plot(t,sol[:,0],color='r',label='u')
     #   plt.plot(sol[:,1],sol[:,3],label='u(0)={}'.format(y0[0]))
#        plt.plot(t,sol[:,2],color='gold',label='z')
#        plt.plot(t,[sol[i,3]**2/sol[i,0] for i in range(len(sol[:,3]))],color='blue',label='x')
        plt.plot(t,neki,color='blue',label='x')
#        plt.plot(t,sol[:,3],color='blue',label='x')
#        plt.plot(t,sol[:,4],color='violet',label='y')
#        plt.plot(t,sol[:,1],color='green',label='v')
        plt.ylabel('koncentracije')
        plt.xlabel('čas')
        plt.legend()
        plt.title('konstanta k/m  v=15 in x=0,01')
        plt.grid()
        plt.show()
    
else:    
    import heapq
    maks=[]
    gji=[]
    for g in np.linspace(50,71,3): 
        
        zacetniI=[]
        i=0.1
        j=0
        cas=[]
                        
        while j<220:
            y0=[10,g,0,j]
            t = np.linspace(0,i,100)
                        
            sol = odeint(ure, y0, t, args=(2,200))#(ime funkcije, začetni pogoji, čas, argumenti)
            vsota=sol[:,0]+sol[:,3]+sol[:,1]+sol[:,2] #(funkcija sol vrne [[],[],[],[]] in zato seštejem vse, da dobim lepo krivuljo)
                       
            if (vsota[-1]-vsota[-3]) > -0.001 and (vsota[-1]-vsota[-3]) <0.001: #ta pogoj mi pove, kdaj naj neham povečevati i (dolžino časovnega intervala)
                zacetniI.append(j)
                cas.append(i)
                i=0.1 #ponastavim čas na začetno vrednost
                j=j+1 #tu potem povečam začetno koncentracijo tiosulfata in začnem z novim ciklom
            else:
                i=i+0.2   #to je zanka v zanki z njo povečujem čas
                
#        tra=heapq.nlargest(2, range(len(cas)), key=cas.__getitem__)
#        maks.append(np.argmax(cas))
#        gji.append(g)
                
            plt.plot(zacetniI,cas)
            plt.ylabel('najdaljše trajanje reakcije')
            plt.xlabel('začetna vrednost persulfata')
            plt.grid()
            plt.title('razmerje je 100')
            plt.legend()
            plt.show()