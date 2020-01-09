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
sns.set_palette(sns.color_palette("cool", 6))

#%%
#----------------------------------------------------------------------------------------------
#       PREIZKUŠANJE RAZNIH GENERATORJEV
#----------------------------------------------------------------------------------------------


start=timer()
t=random() # od njega je čas ko enakomerno zafila histogram z natančnostjo 500 okoli 10 000 000 korakov
end=timer()
print('Število:',t,'čas:',end-start,'\n\n\n')


start=timer()
t=randbelow(1000000)/1000000 # od njega je čas ko enakomerno zafila histogram z natančnostjo 500 okoli 10 000 000 korakov, ampak rata po 1 000 000 koraku zelo počasen
end=timer()
print('Število:',t,'čas:',end-start,'\n\n\n')


start=timer()
t=SystemRandom().random() # od njega je čas ko enakomerno zafila histogram z natančnostjo 500 okoli 10 000 000 korakov, ampak rata po 1 000 000 koraku zelo počasen
end=timer()
print('Število:',t,'čas:',end-start,'\n\n\n')


def casrand(n):
    zahtevnostrand=[]
    zahtevnostbelow=[]
    zahtevnostsist=[]
    korak=[i*1000 for i in range(n)]
    for i in korak:
        start=timer()
        razporeditev=[randbelow(1000000)/1000000 for k in range(i)]
        end=timer()
        zahtevnostbelow.append(end-start)
    for i in korak:
        start=timer()
        razporeditev=[SystemRandom().random() for k in range(i)]
        end=timer()
        zahtevnostsist.append(end-start)
    for i in korak:
        start=timer()
        razporeditev=[random() for k in range(i)]
        end=timer()
        zahtevnostrand.append(end-start) 
    
    plt.title('Časovni zahtevnosti klicev naključnih števil')
    plt.plot(korak,zahtevnostrand,'r:',label='ranom()')
    plt.plot(korak,zahtevnostbelow,'b:',label='randbelow()')
    plt.plot(korak,zahtevnostsist,'g:',label='SystemRandom().random()')
    plt.xlabel('število števil')
    plt.ylabel('čas generiranja števil')
    plt.legend()

#%% 
###########################################################################################################
#
#       1. NALOGA - časovna zahtevnost - GAUSSOVA PORAZDELITEV
#
###########################################################################################################


    
def cas(n):
    zahtevnostkon=[]
    zahtevnostmull=[]
    zahtevnostpol=[]
    korak=[i*100 for i in range(n)]
    for i in korak:
        start=timer()
        razporeditev=[np.sqrt(-2*np.log(random()))*np.cos(2*math.pi*random()) for k in range(i)]
        end=timer()
        zahtevnostmull.append(end-start)
    for i in korak:
        start=timer()
        razporeditev=[(random()+random()+random()+random()+random()+random()-(random()+random()+random()+random()+random()+random())) for k in range(i)]
        end=timer()
        zahtevnostkon.append(end-start)
    for i in korak:
        start=timer()
        u=random()
        v=random()
        s=u*u+v*v
        razporeditev=[u*np.sqrt(-2*np.log(s)/s) for k in range(i)]
        end=timer()
        zahtevnostpol.append(end-start) 
    
    plt.title('Časovni zahtevnosti generatorjev')
    plt.plot(korak,zahtevnostkon,'r:',label='Konvolucijski generator')
    plt.plot(korak,zahtevnostmull,'b:',label='Box-Mullerjev generator')
    plt.plot(korak,zahtevnostpol,'g:',label='Marsaglijev polarni generator')
    plt.xlabel('število števil')
    plt.ylabel('čas generiranja števil')
    plt.legend()



#%%
#%%# 
###########################################################################################################
#
#       1. NALOGA - uvedem funkcije ki modelirajo in testorajo porazdelitev - GAUSSOVA PORAZDELITEV
#
###########################################################################################################

koren=np.sqrt(2)
#tu smo le preverjali periodo
def perioda(n):
    randic=[SystemRandom().random() for i in range(n)]
    plt.hist(randic,500,normed=1, facecolor='green', alpha=0.75)

#sedaj pa uporabimo Muller-box ovo funkcijo
def mullerbox(n):
    razporeditev=[np.sqrt(-2*np.log(random()))*np.cos(2*math.pi*random()) for i in range(n)]
    g = np.histogram(razporeditev,15)
    fitnormalne=np.array([1/2*erf(g[1][i+1]/koren)-1/2*erf(g[1][i]/koren) for i in range(len(g[1])-1)])
    fitnormalne=np.multiply(fitnormalne,n)
    chisq,pstat=chisquare(g[0],fitnormalne)
    g = plt.hist(razporeditev,15,normed = True, facecolor='green', alpha=0.75)
    Dstat,pkol=kstest(razporeditev,'norm')
    t=np.linspace(-3,3,200)
    plt.plot(t,1/(np.sqrt(2*math.pi))*np.exp(-t*t/2),'r--',label='$D$ = {0:.3f} in $p$ = {1:.3f}'.format(Dstat,pkol))
    plt.plot(t,1/(np.sqrt(2*math.pi))*np.exp(-t*t/2),'r--',label='$\chi^2$ = {0:.3f} in $p$ = {1:.3f}'.format(chisq,pstat))
    plt.plot(t,1/(np.sqrt(2*math.pi))*np.exp(-t*t/2),'r--',label='N( $\mu$={0:.2f} , $\sigma$={1:.2f})'.format(0,1))
    plt.title('Box-Mullerjev generator za n={}'.format(n))
    plt.xlabel('Intervali')
    plt.ylabel('Število naključnih števil na intervalu (normirano)')
    plt.show()
    plt.legend()

# 9*random()-9*random() 
def konvolucija(n):
    razporeditev=[(random()+random()+random()+random()+random()+random()-(random()+random()+random()+random()+random()+random())) for i in range(n)]
    g = np.histogram(razporeditev,15)
    fitnormalne=np.array([1/2*erf(g[1][i+1]/koren)-1/2*erf(g[1][i]/koren) for i in range(len(g[1])-1)])
    fitnormalne=np.multiply(fitnormalne,n)
    chisq,pstat=chisquare(g[0],fitnormalne)
    g=plt.hist(razporeditev,15,normed=True, facecolor='green', alpha=0.75)
    Dstat,pkol=kstest(razporeditev,'norm')
    t=np.linspace(-3,3,200)
    plt.plot(t,1/(np.sqrt(2*math.pi))*np.exp(-t*t/2),'r--',label='$D$ = {0:.3f} in $p$ = {1:.3f}'.format(Dstat,pkol))
    plt.plot(t,1/(np.sqrt(2*math.pi))*np.exp(-t*t/2),'r--',label='$\chi^2$ = {0:.3f} in $p$ = {1:.3f}'.format(chisq,pstat))
    plt.plot(t,1/(np.sqrt(2*math.pi))*np.exp(-t*t/2),'r--',label='N( $\mu$={0:.2f} , $\sigma$={1:.2f})'.format(0,1))
    plt.title('Konvolucijski generator za n={}'.format(n))
    plt.xlabel('Intervali')
    plt.ylabel('Število naključnih števil na intervalu (normirano)')
    plt.legend()
    print(pkol)
    


#Z_0 = R \cos(\Theta) =\sqrt{-2 \ln U_1} \cos(2 \pi U_2)\,


#%%
#----------------------------------------------------------------------------------------------
#       1. NALOGA - IZRIŠEMO SEDAJ GRAFE IN JIH SHRANIMO V MAPO
#----------------------------------------------------------------------------------------------



for i in [100,1000,10000,100000,1000000]:
    f=plt.figure(i)
    mullerbox(i)
    f.savefig("mullerbox{}nova.pdf".format(i), bbox_inches='tight')

for i in [100,1000,10000,100000,1000000]:
    f=plt.figure(i+1)
    konvolucija(i)
    f.savefig("konvolucija{}nova.pdf".format(i), bbox_inches='tight')
plt.show()


#%%# 
###########################################################################################################
#
#       1. NALOGA - Rišem porazdelitve statistik - GAUSSOVA PORAZDELITEV
#
###########################################################################################################

def statistikeBMchi(for_zanka,n):
    chi=[]
    for j in range(for_zanka):
        razporeditev=[np.sqrt(-2*np.log(random()))*np.cos(2*math.pi*random()) for i in range(n)]
        g=np.histogram(razporeditev,15)
        fitnormalne=[1/2*erf(g[1][i+1]/koren)-1/2*erf(g[1][i]/koren) for i in range(len(g[1])-1)]
        fitnormalne=np.multiply(fitnormalne,n)
        chisq,pstat=chisquare(g[0],fitnormalne)
        chi.append(chisq)
    plt.hist(chi, bins='fd', normed=True,facecolor='orange',alpha=0.75)
    t=np.arange(0,200,0.1)
    plt.plot(t,chi2.pdf(t,15-1),'black', alpha=0.6)
    plt.title('Porazdelitev $\chi^2$ za n={0} in št. korakov ={1} (Box-Muller)'.format(n,for_zanka))
    plt.xlabel('Intervali - velikost $\chi^2$')
    plt.ylabel('Število $\chi^2$ na nekem intervalu')
    plt.legend()

def statistikeKONchi(for_zanka,n):
    chi=[]
    for j in range(for_zanka):
        razporeditev=[(random()+random()+random()+random()+random()+random()-(random()+random()+random()+random()+random()+random())) for i in range(n)]
        g=np.histogram(razporeditev,15)
        fitnormalne=[1/2*erf(g[1][i+1]/koren)-1/2*erf(g[1][i]/koren) for i in range(len(g[1])-1)]
        fitnormalne=np.multiply(fitnormalne,n)
        chisq,pstat=chisquare(g[0],fitnormalne)
        chi.append(chisq)
    plt.hist(chi, bins='fd',normed=True,facecolor='orange',alpha=0.75)
    t=np.arange(0,200,0.1)
    plt.plot(t,chi2.pdf(t,15-1),'black', alpha=0.6)
    plt.title('Porazdelitev $\chi^2$ za n={0} in št. korakov ={1} (Konvolucija)'.format(n,for_zanka))
    plt.xlabel('Intervali - velikost $\chi^2$')
    plt.ylabel('Število $\chi^2$ na nekem intervalu')
    plt.legend()


def statistikeBMkol(for_zanka,n):
    chi=[]
    for j in range(for_zanka):
        razporeditev=[np.sqrt(-2*np.log(random()))*np.cos(2*math.pi*random()) for i in range(n)]
        Dstat,pkol=Dstat,pkol=kstest(razporeditev,'norm')
        chi.append(Dstat)
    plt.hist(chi, bins='fd', normed=1,facecolor='orange',alpha=0.75)
    plt.title('Porazdelitev Kolmogorova za n={0} in št. korakov ={1} (Box-Muller)'.format(n,for_zanka))
    plt.xlabel('Intervali - velikost $\chi^2$')
    plt.ylabel('Število $\chi^2$ na nekem intervalu')


def statistikeKONkol(for_zanka,n):
    chi=[]
    for j in range(for_zanka):
        razporeditev=[(random()+random()+random()+random()+random()+random()-(random()+random()+random()+random()+random()+random())) for i in range(n)]
        Dstat,pkol=Dstat,pkol=kstest(razporeditev,'norm')
        chi.append(Dstat)
    plt.hist(chi, bins='fd', normed=1,facecolor='orange',alpha=0.75)
    plt.title('Porazdelitev Kolmogorova za n={0} in št. korakov ={1} (Konvolucija)'.format(n,for_zanka))
    plt.xlabel('Intervali - velikost $\chi^2$')
    plt.ylabel('Število $\chi^2$ na nekem intervalu')

#%%
#----------------------------------------------------------------------------------
#            1. NALOGA - rišemo zdAJ grafe in jih shranimo za porazdelitve statistik
#----------------------------------------------------------------------------------

f=plt.figure(0)
statistikeBMchi(10000,1000)
f.savefig("statistikeBMchi.pdf".format(i), bbox_inches='tight')
f=plt.figure(1)
statistikeBMkol(10000,1000)
f.savefig("statistikeBMkol.pdf".format(i), bbox_inches='tight')
f=plt.figure(2)
statistikeKONchi(10000,1000)
f.savefig("statistikeKONchi.pdf".format(i), bbox_inches='tight')
f=plt.figure(3)
statistikeKONkol(10000,1000)
f.savefig("statistikeKONkol.pdf".format(i), bbox_inches='tight')
plt.show()


#%%
#%% 
###########################################################################################################
#
#       2. NALOGA - PORAZDELITEV FOTONOV PO PROSTORSKEM KOTU
#
###########################################################################################################

#Tu pač kličem funkcije
sfera(1000)
sfera(100)



#%%
#------------------------------------------------------------------------------
#           2.NALOGA - funkcije za risanje krogle
#------------------------------------------------------------------------------

def set_axes_equal(ax):

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])
    
    
    
def sfera(stevilo_kroglic):
    N=stevilo_kroglic
    
    # PLOT
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    
    
    # sfera
    u_s = np.linspace(0, 2*np.pi, 100)
    v_s = np.linspace(0, np.pi, 100)
    x_s = np.outer(np.cos(u_s), np.sin(v_s))
    y_s = np.outer(np.sin(u_s), np.sin(v_s))
    z_s = np.outer(np.ones(np.size(u_s)), np.cos(v_s))
    
    ax.plot_surface(x_s, y_s, z_s, color='grey',lw=0., alpha=0.10)
    
    
    # (podani fi, theta izmenično)        
    t1 = [np.arccos(random()*2-1) for k in range(N)]
    f1 = [random()*2*np.pi for k in range(N)]
    x1 = np.cos(f1)*np.sin(t1)
    y1 = np.sin(f1)*np.sin(t1)
    z1 = np.cos(t1)
    ax.scatter(x1, y1, z1,color="r",s=3)
    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    x = np.cos(u)*np.sin(v)
    y = np.sin(u)*np.sin(v)
    z = np.cos(v)
    ax.plot_wireframe(x, y, z, color="black",linewidth=0.7,alpha=0.4)

    
    
    
    ax.set_aspect("equal")
    set_axes_equal(ax)
    fig.tight_layout()
    plt.title('Naključno izžrebane smeri za n={}'.format(N))
    plt.axis('off')

def sferasin(stevilo_kroglic):
    N=stevilo_kroglic
    
    # PLOT
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    
    
    # sfera
    u_s = np.linspace(0, 2*np.pi, 100)
    v_s = np.linspace(0, np.pi, 100)
    x_s = np.outer(np.cos(u_s), np.sin(v_s))
    y_s = np.outer(np.sin(u_s), np.sin(v_s))
    z_s = np.outer(np.ones(np.size(u_s)), np.cos(v_s))
    
    ax.plot_surface(x_s, y_s, z_s, color='grey',lw=0., alpha=0.10)
    
    
    # (podani fi, theta izmenično)        
    t1 = [math.pi*np.sin(random()*math.pi)**2 for k in range(N)]
    f1 = [random()*2*np.pi for k in range(N)]
    x1 = np.cos(f1)*np.sin(t1)
    y1 = np.sin(f1)*np.sin(t1)
    z1 = np.cos(t1)
    ax.scatter(x1, y1, z1,color="r",s=1)
    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    x = np.cos(u)*np.sin(v)
    y = np.sin(u)*np.sin(v)
    z = np.cos(v)
    ax.plot_wireframe(x, y, z, color="black",linewidth=0.7,alpha=0.4)

    
    
    
    ax.set_aspect("equal")
    set_axes_equal(ax)
    fig.tight_layout()
    plt.title('Naključno izžrebane smeri za n={} porazdelitev $ \sin()^2$'.format(N))
    plt.axis('off')
#%%
#------------------------------------------------------------------------------
#           2.NALOGA - funkcije za risanje dipola na krogli
#------------------------------------------------------------------------------
    
    
def F(x):
    return (3/4)*(np.cos(x)*np.cos(x)*np.cos(x)/3-np.cos(x)+2/3)


def sferadip(stevilo_kroglic):
    N=stevilo_kroglic
    
    # PLOT
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    
    
    # sfera
    u_s = np.linspace(0, 2*np.pi, 100)
    v_s = np.linspace(0, np.pi, 100)
    x_s = np.outer(np.cos(u_s), np.sin(v_s))
    y_s = np.outer(np.sin(u_s), np.sin(v_s))
    z_s = np.outer(np.ones(np.size(u_s)), np.cos(v_s))
    
    ax.plot_surface(x_s, y_s, z_s, color='grey',lw=0., alpha=0.10)
    
    
    # (podani fi, theta izmenično)        
    t1=[]
    
    for i in range(N):
        u=10**20
        x=1
        while u > (3/4)*np.sin(x)**3/(2):
            u=random()*math.pi
            x=random()*math.pi
        t1.append(x)
            
    f1 = [random()*2*np.pi for k in range(N)]
    x1 = np.cos(f1)*np.sin(t1)
    y1 = np.sin(f1)*np.sin(t1)
    z1 = np.cos(t1)
    ax.scatter(x1, y1, z1,color="r",s=1)
    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    x = np.cos(u)*np.sin(v)
    y = np.sin(u)*np.sin(v)
    z = np.cos(v)
    ax.plot_wireframe(x, y, z, color="black",linewidth=0.7,alpha=0.4)

    
    
    
    ax.set_aspect("equal")
    set_axes_equal(ax)
    fig.tight_layout()
    plt.title('Naključno izžrebane smeri za n={} porazdelitev $ \sin()^2$'.format(N))
    plt.axis('off')


#%%
#------------------------------------------------------------------------------
#           2.NALOGA - funkcije za verjetnostne porazdelitve
#------------------------------------------------------------------------------

    
def razporeditev(N):
    t1 = [np.arccos(random()*2-1) for k in range(N)]
    plt.title('verjetnostna porazdelitev po kotih enakomerno po smereh')
    plt.hist(t1,bins='fd',normed=1,label='število števil n={}'.format(N))
    plt.xlabel('kot v radianih')
    plt.ylabel('število zadetkov na intervalu (normirano)')
    plt.legend()  
    
    
    


def razporeditevsin2(N):
    t1=[]
    for i in range(N):
        u=10**20
        x=1
        while u > (3/4)*np.sin(x)**3/(2):
            u=random()*math.pi
            x=random()*math.pi
        t1.append(x)
    plt.title('verjetnostna porazdelitev po kotih za dipolno sevanje $\sin()^2$')
    plt.hist(t1,bins='fd',normed=1,label='število števil n={}'.format(N))
    t=np.arange(0,math.pi,0.1)
    plt.plot(t,3/4*np.sin(t)**3,'g--',label='$3/4 sin()^3$')
    plt.plot(t,F(t),'r:',label='F()')
    plt.xlabel('kot v radianih')
    plt.ylabel('število zadetkov na intervalu (normirano)')
    plt.legend()

#%%
#------------------------------------------------------------------------------
#           2.NALOGA - funkcije za risanje momentov
#------------------------------------------------------------------------------




def dipolna(N):
    momi=[]
    for j in range(N):
        t1=[]
        for i in range(j):
            u=10**20
            x=1
            while u > (3/4)*np.sin(x)**3/(2):
                u=random()*math.pi
                x=random()*math.pi
            t1.append(x)
        momi.append(moment(np.cos(t1),moment=[1,2,3,4,5,6]))
    t=[i for i in range(N)]
    plt.title('razni momenti pri dipolni porazdelitvi')
    plt.plot(t,[momi[i][0] for i in range(N)],label='<cos()>')
    plt.plot(t,[momi[i][1] for i in range(N)],label='<cos()^2>')
    plt.plot(t,[momi[i][2] for i in range(N)],label='<cos()^3>')
    plt.plot(t,[momi[i][3] for i in range(N)],label='<cos()^4>')
    plt.plot(t,[momi[i][4] for i in range(N)],label='<cos()^5>')
    plt.plot(t,[momi[i][5] for i in range(N)],label='<cos()^6>')
    plt.plot(t,[0 for i in range(N)],'r:',label='teoretično vsi lihi momenti=0')
    plt.plot(t,[1/5 for i in range(N)],'r:',label='teoretično <cos()^2>=1/5')
    plt.plot(t,[3/35 for i in range(N)],'y:',label='teoretično <cos()^4>=3/35')
    plt.plot(t,[3/63 for i in range(N)],'g:',label='teoretično <cos()^6>=3/63')
#    plt.plot(t,povprecje2,'b:',label='<cos()^2>')
    plt.xlabel('dolžina korakov')
    plt.ylabel('vrednost momentov')
    plt.legend(loc=1)
    
    
def navadna(N):
    momi=[]
    for j in range(N):
        t1 = [np.arccos(random()*2-1) for k in range(N)]
        momi.append(moment(np.cos(t1),moment=[1,2,3,4,5,6]))
    t=[i for i in range(N)]
    plt.title('razni momenti pri enakomerni porazdelitvi')
    plt.plot(t,[momi[i][0] for i in range(N)],label='<cos()>')
    plt.plot(t,[momi[i][1] for i in range(N)],label='<cos()^2>')
    plt.plot(t,[momi[i][2] for i in range(N)],label='<cos()^3>')
    plt.plot(t,[momi[i][3] for i in range(N)],label='<cos()^4>')
    plt.plot(t,[momi[i][4] for i in range(N)],label='<cos()^5>')
    plt.plot(t,[momi[i][5] for i in range(N)],label='<cos()^6>')
    plt.plot(t,[0 for i in range(N)],'r:',label='teoretično vsi lihi momenti=0')
    plt.plot(t,[1/3 for i in range(N)],'r:',label='teoretično <cos()^2>=1/3')
    plt.plot(t,[1/5 for i in range(N)],'y:',label='teoretično <cos()^4>=1/5')
    plt.plot(t,[1/7 for i in range(N)],'g:',label='teoretično <cos()^6>=1/7')
    plt.xlabel('dolžina korakov')
    plt.ylabel('vrednost momentov')
    plt.legend()



#%%
#%% 
###########################################################################################################
#
#       3. NALOGA - ODDAJANJE NALOG
#
###########################################################################################################

letnik14=[]
for j in range(13):
    j=j+101
    string='mod_tm14_{}.dat'.format(j)
    with open(string) as f:
        podaci=[l.strip().split(":") for l in f]
        letnik14.append(podaci)
letnik14=[[24*60*int(letnik14[i][j][0])+60*int(letnik14[i][j][1])+int(letnik14[i][j][2]) for j in range(len(letnik14[i]))] for i in range(len(letnik14))]


letnik13=[]
for j in range(13):
    j=j+101
    string='mod_tm13_{}.dat'.format(j)
    with open(string) as f:
        podaci=[l.strip().split(":") for l in f]
        letnik13.append(podaci)
letnik13=np.asarray([[24*60*float(letnik13[i][j][0])+60*float(letnik13[i][j][1])+float(letnik13[i][j][2]) for j in range(len(letnik13[i]))] for i in range(len(letnik13))])



letnik11=[]
for j in range(12):
    j=j+101
    string='mod_tm11_{}.dat'.format(j)
    with open(string) as f:
        podaci=[l.strip().split(":") for l in f]
        letnik11.append(podaci)
letnik11=[[24*60*float(letnik11[i][j][0])+60*float(letnik11[i][j][1])+float(letnik11[i][j][2]) for j in range(len(letnik11[i]))] for i in range(len(letnik11))]


letnik10=[]
for j in range(13):
    j=j+101
    string='mod_tm10_{}.dat'.format(j)
    with open(string) as f:
        podaci=[l.strip().split(":") for l in f]
        letnik10.append(podaci)
letnik10=[[24*60*float(letnik10[i][j][0])+60*float(letnik10[i][j][1])+float(letnik10[i][j][2]) for j in range(len(letnik10[i]))] for i in range(len(letnik10))]

#iks.append(24*60*float(i[0])+60*float(i[1])+float(i[2]))


neki=[]
for i in range(len(letnik14)):
    neki = neki+letnik14[i]
for i in range(len(letnik13)):
    neki = neki+letnik13[i]
for i in range(len(letnik11)):
    neki = neki+letnik11[i]
for i in range(len(letnik10)):
    neki = neki+letnik10[i]
nekiure=[neki[i]/(60*24) for i in range(len(neki))]

#%%
#------------------------------------------------------------------------------
#           3.NALOGA - splošno
#------------------------------------------------------------------------------

plt.title('Oddajanje nalog v odvisnosti od časa')
plt.hist(nekiure, bins='fd',facecolor='green',alpha=0.55)
plt.plot([0 ,0],[0,110],'r:')
plt.xlabel('čas $[dnevi]$')
plt.ylabel('število oddanih nalog')
plt.show()

#%%
#------------------------------------------------------------------------------------------------------------------------------------------------------------
#           3.NALOGA - oddajanje nalog v odvisnoti od časa na intervalu [-dan,dan]
#------------------------------------------------------------------------------------------------------------------------------------------------------------



seznam=[letnik14[i]+letnik13[i]+letnik11[i]+letnik10[i] for i in range(12)]
tisto=[np.histogram(seznam[i],bins=(-24*60,24*60))[0][0] for i in range(12)]
Amezo=[[1,i+1] for i in range(12)]
solmezo,_,_,_= np.linalg.lstsq(Amezo,tisto)
resitev=np.dot(Amezo,solmezo)
ksimezo=sum((tisto[i]-resitev[i])**2 for i in range(12))/(12-2)
t=[i+1 for i in range(12)]

plt.title('Število oddanih nalog v intervau $[-dan,dan]$ okoli roka za oddajo')
plt.plot(t,tisto,'r.',label='podatki')
plt.plot(t,resitev,'b:',label='reducirani $\chi^2$={0:.2f}'.format(ksimezo))
plt.plot(t,resitev,'b:',label='prilagoditvena linearna funkcija\n k={0:.2f} in c={1:.2f}'.format(solmezo[1],solmezo[0]))
plt.xlabel('Številka naloge')
plt.ylabel('Števil oddanih nalog na intervalu $[-dan,dan]$')
plt.legend()
#%%
#------------------------------------------------------------------------------------------------------------------------------------------------------------
#           3.NALOGA - oddajanje nalog v odvisnoti od časa na intervalu [-dan,dan] za vsak letnik posebej
#------------------------------------------------------------------------------------------------------------------------------------------------------------


seznam14=[letnik14[i] for i in range(12)]
seznam13=[letnik13[i] for i in range(12)]
seznam11=[letnik11[i] for i in range(12)]
seznam10=[letnik10[i] for i in range(12)]

tisto14=[np.histogram(seznam14[i],bins=(-24*60,24*60))[0][0] for i in range(12)]
tisto13=[np.histogram(seznam13[i],bins=(-24*60,24*60))[0][0] for i in range(12)]
tisto11=[np.histogram(seznam11[i],bins=(-24*60,24*60))[0][0] for i in range(12)]
tisto10=[np.histogram(seznam10[i],bins=(-24*60,24*60))[0][0] for i in range(12)]

Amezo=[[1,i+1] for i in range(12)]

solmezo14,_,_,_= np.linalg.lstsq(Amezo,tisto14)
solmezo13,_,_,_= np.linalg.lstsq(Amezo,tisto13)
solmezo11,_,_,_= np.linalg.lstsq(Amezo,tisto11)
solmezo10,_,_,_= np.linalg.lstsq(Amezo,tisto10)

resitev14=np.dot(Amezo,solmezo14)
resitev13=np.dot(Amezo,solmezo13)
resitev11=np.dot(Amezo,solmezo11)
resitev10=np.dot(Amezo,solmezo10)

ksimezo14=sum((tisto14[i]-resitev14[i])**2 for i in range(12))/(12-2)
ksimezo13=sum((tisto13[i]-resitev13[i])**2 for i in range(12))/(12-2)
ksimezo11=sum((tisto11[i]-resitev11[i])**2 for i in range(12))/(12-2)
ksimezo10=sum((tisto10[i]-resitev10[i])**2 for i in range(12))/(12-2)

t=[i+1 for i in range(12)]

plt.title('Število oddanih nalog v intervau $[-dan,dan]$ okoli roka za oddajo')

plt.plot(t,tisto14,'.',color='red',label='podatki 2014')
plt.plot(t,resitev14,color='red',linestyle=':',label='reducirani $\chi^2$={0:.2f}'.format(ksimezo14))
plt.plot(t,resitev14,color='red',linestyle=':',label='prilagoditvena linearna funkcija\n k={0:.2f} in c={1:.2f}'.format(solmezo14[1],solmezo14[0]))

plt.plot(t,tisto13,'.',color='gold',label='podatki 2013')
plt.plot(t,resitev13,color='gold',linestyle=':',label='reducirani $\chi^2$={0:.2f}'.format(ksimezo13))
plt.plot(t,resitev13,color='gold',linestyle=':',label='prilagoditvena linearna funkcija\n k={0:.2f} in c={1:.2f}'.format(solmezo13[1],solmezo13[0]))

plt.plot(t,tisto11,'.',color='blue',label='podatki 2011')
plt.plot(t,resitev11,color='blue',linestyle=':',label='reducirani $\chi^2$={0:.2f}'.format(ksimezo11))
plt.plot(t,resitev11,color='blue',linestyle=':',label='prilagoditvena linearna funkcija\n k={0:.2f} in c={1:.2f}'.format(solmezo11[1],solmezo11[0]))

plt.plot(t,tisto10,'.',color='green',label='podatki 2010')
plt.plot(t,resitev10,color='green',linestyle=':',label='reducirani $\chi^2$={0:.2f}'.format(ksimezo10))
plt.plot(t,resitev10,color='green',linestyle=':',label='prilagoditvena linearna funkcija\n k={0:.2f} in c={1:.2f}'.format(solmezo10[1],solmezo10[0]))

plt.xlabel('Številka naloge')
plt.ylabel('Števil oddanih nalog na intervalu $[-dan,dan]$')
plt.legend()

#%%
###########################################################################################################
#
#       3. NALOGA - KOlMOGOROV
#
###########################################################################################################

#-----------------------------------------------------------------------------------------------------------------------
# verjetnost p
#-----------------------------------------------------------------------------------------------------------------------

m=[[ks_2samp(letnik10[i],letnik11[i])[1],ks_2samp(letnik10[i],letnik13[i])[1],ks_2samp(letnik10[i],letnik14[i])[1],ks_2samp(letnik11[i],letnik13[i])[1],ks_2samp(letnik11[i],letnik14[i])[1],ks_2samp(letnik13[i],letnik14[i])[1]] for i in range(12)]
m=[[m[i][j] for i in range(len(m))]for j in range(len(m[0]))]

matrix = np.matrix(m)
matrix
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.set_aspect('equal')
plt.title('Podobnost med letniki - verjetnost p')
plt.imshow(matrix, interpolation='nearest', cmap=plt.cm.YlGn,
    extent=(0.5,np.shape(matrix)[1]+0.5,0.5,np.shape(matrix)[0]+0.5))
plt.colorbar()
plt.xlabel('Številka naloge')
plt.ylabel('letnik z letnikom')
plt.show()

#%%

#-----------------------------------------------------------------------------------------------------------------------
# statistika D
#-----------------------------------------------------------------------------------------------------------------------

m=[[ks_2samp(letnik10[i],letnik11[i])[0],ks_2samp(letnik10[i],letnik13[i])[0],ks_2samp(letnik10[i],letnik14[i])[0],ks_2samp(letnik11[i],letnik13[i])[0],ks_2samp(letnik11[i],letnik14[i])[0],ks_2samp(letnik13[i],letnik14[i])[0]] for i in range(12)]
m=[[m[i][j] for i in range(len(m))]for j in range(len(m[0]))]

matrix = np.matrix(m)
matrix
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.set_aspect('equal')
plt.title('Podobnost med letniki - vrednost D')
plt.imshow(matrix, interpolation='nearest', cmap=plt.cm.YlGn,
    extent=(0.5,np.shape(matrix)[1]+0.5,0.5,np.shape(matrix)[0]+0.5))
plt.colorbar()
plt.xlabel('Številka naloge')
plt.ylabel('letnik z letnikom')
plt.show()


#%%

#-----------------------------------------------------------------------------------------------------------------------
# naloge med sabo - statistika D
#-----------------------------------------------------------------------------------------------------------------------

seznam=[letnik11[i]+letnik13[i]+letnik14[i] for i in range(12)]
m=[[ks_2samp(seznam[i],seznam[j])[0] if j>i else 0 for j in range(11)] for i in range(11)]
m=[[m[10-i][j] for j in range(11)] for i in range(11)]


matrix = np.matrix(m)
matrix
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.set_aspect('equal')
plt.title('Podobnost med nalogami - vrednost D')
plt.imshow(matrix, interpolation='nearest', cmap=plt.cm.YlGn,
    extent=(0.5,np.shape(matrix)[1]+0.5,0.5,np.shape(matrix)[0]+0.5))
plt.colorbar()
plt.xlabel('naloga s katero primerjamo')
plt.ylabel('naloga ki jo primerjamo z nalogami')
plt.show()


#%%

#-----------------------------------------------------------------------------------------------------------------------
# naloge med sabo - verjetnost p
#-----------------------------------------------------------------------------------------------------------------------

seznam=[letnik11[i]+letnik13[i]+letnik14[i] for i in range(12)]
m=[[ks_2samp(seznam[i],seznam[j])[1] if j>i else 0 for j in range(11)] for i in range(11)]
m=[[m[10-i][j] for j in range(11)] for i in range(11)]


matrix = np.matrix(m)
matrix
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.set_aspect('equal')
plt.title('Podobnost med nalogami - vrednost p')
plt.imshow(matrix, interpolation='nearest', cmap=plt.cm.YlGn,
    extent=(0.5,np.shape(matrix)[1]+0.5,0.5,np.shape(matrix)[0]+0.5))
plt.colorbar()
plt.xlabel('naloga s katero primerjamo')
plt.ylabel('naloga ki jo primerjamo z nalogami')
plt.show()










