# Code-for-High_Tc-dataimport numpy as np
from matplotlib import pyplot as plt
import matplotlib
from scipy.odr import *
from scipy.stats import chisquare
import scipy as sp
from scipy.constants import *
%matplotlib nbagg

def Read(Dateiname, nValues=1,X_then_dX=True, Skiprows=1):   #X, dX, Y, dY = StandardRead(test.txt, 2)
    Values=[]
    if X_then_dX==True:
        nValues=2*nValues
    for i in range(nValues):
        Values.append(np.array(np.loadtxt(Dateiname, skiprows=1, usecols=(i,))))
    return Values
    
    Energy, de2_250fs = Read("de2_250fs_digitized.txt", 2, False)

Energy_reduced=[]
for i in Energy:
    if i<2.6:
        Energy_reduced.append(i)
Energy_reduced =np.array(Energy_reduced)
# e2_0_reduced=np.array([e2_0[i] for i in range(len(Energy_reduced))])
de2_250fs_reduced = np.array([de2_250fs[i] for i in range(len(Energy_reduced))])
# e2_10ps_reduced=np.array([e2_10ps[i] for i in range(len(Energy_reduced))])
# de2_350fs_reduced=np.array([e2_0[i] for i in range(len(Energy_reduced))])
# de2_10ps_reduced=np.array([e2_0[i] for i in range(len(Energy_reduced))])


dEnergy_reduced = np.array([0.01 for i in range(len(Energy_reduced))])

dde2_250fs_reduced  = dEnergy_reduced
def Fit(FitModel, x, y, dx, dy, beta, Name="testFittiings", Title="Forgot the Title", 
                Label="Data", X_Label= "X", Y_Label = "Y", Alpha=0.3, LS="none", color="black", 
                ): 
    model = Model(FitModel)
    
    fig, ax = plt.subplots()
    ax.errorbar(x, y, xerr = dx, yerr = dy, marker = ".", ms = 3, label=Label, alpha = Alpha, ls=LS)
    fitBeta=[]
    fitSDBeta=[]
    
    
        

        
    
    data = RealData(x, y, sx = dx, sy = dy)
    #print(data)
    odr = ODR(data, model, beta0 = beta) 
    fit = odr.run() 
    fit.pprint()
    G = fit.beta
    if (0 not in dx and 0 not in dy):
        print('Chi-square =',fit.res_var,'\n')    #chi^2 macht keinen Sinn, wenn die Fehler Null sind!
    fitBeta.append(fit.beta)
    fitSDBeta.append(fit.sd_beta)
        
    x_plot = np.linspace(x[0],x[len(x)-1],1000)
    y_plot = FitModel(G, x_plot)
        
    ax.plot(x_plot, y_plot, label = "Fit für " + Label, color = color)
      

    
    
    ax.legend(loc='best', fontsize=9)
    ax.set_xlabel(X_Label, fontsize=12)
    ax.set_ylabel(Y_Label, fontsize=12)
    ax.set_title(Title, fontsize=13)
    ax.grid()
    ax.axis('tight')
    plt.savefig(Name+".pdf")
    fig.tight_layout();
    return fitBeta, fitSDBeta
    
    fig, ax = plt.subplots()
# ax.plot(np.linspace(0.1,0.62,1000), Drude([0.64054891, 0.54802584], np.linspace(0.1,0.62,1000)))
ax.plot(Energy_reduced, de2_250fs_reduced, color = 'g')


print(np.linspace(0.00,0.6,100), Drude([20, 0.5], np.linspace(0.00,0.6,100))) #this will print the data
#This will save simulated data as txt file
asdf = np.linspace(0.00,0.6,100)
sdfg = Drude([20, 0.5], np.linspace(0.00,0.6,100))
with open("test1.txt", "w") as f:
    f.write("Energy\tEpsilon\n")
    for i in range (100):
        f.write(str(asdf[i])+"\t"+str(sdfg[i])+"\n")
        
try:
    sdfgdf=1/1
except Exception as ex:
    print(ex)
    if ex == "division by zero": #not fully functional(does not compare correctly), errors
        print("0")
    else:
        print("fuck")
finally:
    print("always")
#------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------
#scipy curve fit method
from scipy.optimize import curve_fit
from scipy.misc import imsave

def Drude1Lorentz2_real(x, wpD, tauD, wpl1, wl1, taul1, wpl2, wl2, taul2, wpl3, wl3, taul3):
    return 5 - wpD**2/(x**2+tauD**(-2))+ (wpl1**2*(wl1**2-x**2))/((wl1**2-x**2)**2+x**2/taul1**2)+(wpl2**2*(wl2**2-x**2))/((wl2**2-x**2)**2+x**2/taul2**2)+ (wpl3**2*(wl3**2-x**2))/((wl3**2-x**2)**2+x**2/taul3**2)

def Drude1Lorentz2_imag(x, wpD, tauD, wpl1, wl1, taul1, wpl2, wl2, taul2, wpl3, wl3, taul3):
    return (1/(x*tauD))*wpD**2/(x**2+tauD**(-2))+ (wpl1**2*x/taul1)/((wl1**2-x**2)**2+x**2/taul1**2)+ (wpl2**2*x/taul2)/((wl2**2-x**2)**2+x**2/taul2**2) + (wpl3**2*x/taul3)/((wl3**2-x**2)**2+x**2/taul3**2)

def Drude1Lorentz2_real_fix(x, wpD, tauD, taul1, taul2):
    wl1=0.6
    wl2=2.3
    wpl1=0.2
    wpl2=1.8
    return 5 - wpD**2/(x**2+tauD**(-2))+ (wpl1**2*(wl1**2-x**2))/((wl1**2-x**2)**2+x**2/taul1**2)+(wpl2**2*(wl2**2-x**2))/((wl2**2-x**2)**2+x**2/taul2**2)
            

def Drude1Lorentz2_imag_fix(x, wpD, tauD, taul1, taul2):
    wl1=0.6
    wl2=2.3
    wpl1=0.2
    wpl2=1.8
    return (1/(x*tauD))*wpD**2/(x**2+tauD**(-2))+ (wpl1**2*x/taul1)/((wl1**2-x**2)**2+x**2/taul1**2)+ (wpl2**2*x/taul2)/((wl2**2-x**2)**2+x**2/taul2**2)
    
   x_data = Energy
x_fit = np.linspace(0.2,10,1000)
y_data_real = e1_100mw
y_data_imag = e2_100mw
yerr_d = de1_100mw_reduced

print(x_data.shape,y_data.shape, yerr_d.shape)
fit_func = Drude1Lorentz2_real_fix
test_func = Drude1Lorentz2_imag_fix
parlables = 'wpD, tauD, taul1, taul2'.split(', ')

wpD_est = 0.5
wpD_tol = .1
tauD_est = 3.16527054
tauD_tol = .1
taul1_est = 4.66922027
taul1_tol = .1
taul2_est = 0.3993293
taul2_tol = .1

# guess = [1.38259923, 0.16527054, 0.23112961, 0.65622163, 4.66922027, 1.90463712, 2.4314406,  0.93993293]
# bounds = ((-np.inf,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf),
#            (np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf))
guess = [wpD_est, tauD_est, taul1_est, taul2_est]
bounds = ((-np.inf,-np.inf,-np.inf,-np.inf),
           (np.inf, np.inf, np.inf, np.inf))
bounds = ((wpD_est*(1-wpD_tol),tauD_est*(1-tauD_tol),taul1_est*(1-taul1_tol),taul2_est*(1-taul2_tol)),
           (wpD_est*(1+wpD_tol), tauD_est*(1+tauD_tol), taul1_est*(1+taul1_tol), taul2_est*(1+taul2_tol)))

#----------------------------------------------------------------------

popt, pcov = curve_fit(fit_func,x_data,y_data_real,p0=guess,bounds=bounds,sigma=yerr_d)
perr = np.sqrt(np.diag(pcov))
f, ax = plt.subplots(1,2,figsize=(9,4))
ax[0].errorbar(x_data,y_data_real,yerr=yerr_d)
ax[0].plot(x_fit,(fit_func(x_fit,*popt)),'--')
ax[0].set_title('Real part')
ax[0].set_xlabel('Energy [eV]')
ax[1].errorbar(x_data,y_data_imag,yerr=yerr_d)
ax[1].plot(x_fit,(test_func(x_fit,*popt)),'--')
ax[1].set_title('imaginary part')
ax[1].set_xlabel('Energy [eV]')

for axis in ax:
    axis.axvline(2.2, linestyle='--',alpha=0.5)
#     axis.axvline(0.6, linestyle='-',alpha=0.5)
    axis.axvline(popt[2], linestyle='--',alpha=0.5)
    axis.axvline(popt[3], linestyle='--',alpha=0.5)
    axis.set_xlim(0,5)

print('par| popt | std')
for i in range(len(popt)):
    print('p{}: {:.3f} | {:.3f}'.format(parlables[i],popt[i],perr[i]))
#-------------------------------------------------------------------------------------------------
def Drude1Lorentz2_real(x, wpD, tauD, wpl1, wl1, taul1, wpl2, wl2, taul2, wpl3, wl3, taul3):
    return 5 - wpD**2/(x**2+tauD**(-2))+ (wpl1**2*(wl1**2-x**2))/((wl1**2-x**2)**2+x**2/taul1**2)+(wpl2**2*(wl2**2-x**2))/((wl2**2-x**2)**2+x**2/taul2**2)+ (wpl3**2*(wl3**2-x**2))/((wl3**2-x**2)**2+x**2/taul3**2)

def Drude1Lorentz2_imag(x, wpD, tauD, wpl1, wl1, taul1, wpl2, wl2, taul2, wpl3, wl3, taul3):
    return (1/(x*tauD))*wpD**2/(x**2+tauD**(-2))+ (wpl1**2*x/taul1)/((wl1**2-x**2)**2+x**2/taul1**2)+ (wpl2**2*x/taul2)/((wl2**2-x**2)**2+x**2/taul2**2) + (wpl3**2*x/taul3)/((wl3**2-x**2)**2+x**2/taul3**2)
    
    popt , pcov = curve_fit(fit_func,x_data,y_data_real,p0=guess,bounds=bounds,sigma=yerr_d)
perr = np.sqrt(np.diag(pcov))
f, ax = plt.subplots(1,2,figsize=(14,7))
ax[0].plot(x_data,y_data_real, 'o')
ax[0].plot(x_fit,(fit_func(x_fit,*popt)),'--')
ax[0].set_title('Real part', 'r')
ax[0].set_xlabel('Energy [eV]', 'r')
ax[1].plot(x_data,y_data_imag, 'o')
ax[1].plot(x_fit,(test_func(x_fit,*popt)),'--')
ax[1].set_title('imaginary part', 'r')
ax[1].set_xlabel('Energy [eV]', 'r')

for axis in ax:
    axis.axvline(2.2, linestyle='--',alpha=0.5)
    axis.axvline(1.3, linestyle='--',alpha=0.5)
    axis.axvline(0.6, linestyle='--',alpha=0.5)
#     axis.axvline(popt[1], linestyle='--',alpha=0.5)
#     axis.axvline(popt[2], linestyle='--',alpha=0.5)
#     axis.axvline(popt[3], linestyle='--',alpha=0.5)
    axis.set_xlim(0,5)

print('par| popt | std')
for i in range(len(popt)):
    print('p{}: {:.3f} | {:.3f}'.format(parlables[i],popt[i],perr[i]))
f.savefig('e1_e2_best fit.png')
# imsave('rgb_gradient.png', )
    
    x_data = Energy
x_fit = np.linspace(0.2,10,1000)
y_data_real = e1_60mw
y_data_imag = e2_60mw
# yerr_d = de1_100mw_reduced

print(x_data.shape,y_data.shape, yerr_d.shape)
fit_func = Drude1Lorentz2_real
test_func = Drude1Lorentz2_imag
parlables = 'wpD, tauD, wpl1, wl1, taul1, wpl2, wl2, taul2, wpl3, wl3, taul3'.split(', ')

wpD_est = 0.55
wpD_tol = 0.1
tauD_est = 1.76527054
tauD_tol = .1
wpl1_est = 0.33
wpl1_tol =0.1
wl1_est = 0.59
wl1_tol = 0.1
taul1_est = 2.96922027
taul1_tol = .1
wpl2_est = 0.45
wpl2_tol = 0.1
wl2_est = 1.4
wl2_tol= 0.1
taul2_est = 1.793293
taul2_tol = .1
wpl3_est = 1.3
wpl3_tol = 0.1
wl3_est = 2.2
wl3_tol = 0.1
taul3_est = 1.8
taul3_tol = 0.1


# guess = [1.38259923, 0.16527054, 0.23112961, 0.65622163, 4.66922027, 1.90463712, 2.4314406,  0.93993293]
# bounds = ((-np.inf,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf),
#            (np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf))
guess = [wpD_est, tauD_est, wpl1_est, wl1_est, taul1_est, wpl2_est, wl2_est, taul2_est, wpl3_est, wl3_est, taul3_est]
bounds = ((-np.inf,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf),
           (np.inf, np.inf, np.inf, np.inf, np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf))
bounds = ((wpD_est*(1-wpD_tol),tauD_est*(1-tauD_tol),wpl1_est*(1-wpl1_tol), wl1_est*(1-wl1_tol),taul1_est*(1-taul1_tol), wpl2_est*(1-wpl2_tol), wl2_est*(1-wl2_tol),taul2_est*(1-taul2_tol), wpl3_est*(1-wpl3_tol), wl3_est*(1-wl3_tol), taul3_est*(1-taul3_tol)),
           (wpD_est*(1+wpD_tol),tauD_est*(1+tauD_tol),wpl1_est*(1+wpl1_tol),wl1_est*(1+wl1_tol), taul1_est*(1+taul1_tol), wpl2_est*(1+wpl2_tol), wl2_est*(1+wl2_tol),taul2_est*(1+taul2_tol), wpl3_est*(1+wpl3_tol), wl3_est*(1+wl3_tol), taul3_est*(1+taul3_tol)))


#---------------------------------------------------------------------------------------------------------------------------------------
#Systematic fitting

import numpy as np
from matplotlib.pyplot import cm
from matplotlib import pyplot as plt
import matplotlib
from scipy.odr import *
from scipy.stats import chisquare
import scipy as sp
from scipy.optimize import curve_fit, leastsq, least_squares
from scipy.misc import imsave
from scipy.constants import *
%matplotlib nbagg


def Read(Dateiname, nValues=1, X_then_dX=True, Skiprows=1):   #X, dX, Y, dY = StandardRead(test.txt, 2)
    Values=[]
    if X_then_dX==True:
        nValues=2*nValues
    for i in range(nValues):
        Values.append(np.array(np.loadtxt(Dateiname, skiprows=1, usecols=(i,))))
    return Values
def to_complex(real,imag):
    y_comp = []
    for i in range(len(real)):
        y_comp.append(np.complex(real[i],imag[i])) 
    return np.array(y_comp)
    
 energy_axis, *data_e1 = Read("e1_values_2ev.txt", 10, False)
_, *data_e2 = Read("e2_values_2ev.txt", 10, False)
data_e2 = np.array(data_e2)
data_e1 = np.array(data_e1)
data_labels = [0,1,5,10,20,40,60,80,100]
print(data_e1.shape, data_e2.shape)

f, ax = plt.subplots(1,2,figsize=(10,5))
for i in range(data_e1.shape[0]):
    ax[0].plot(energy_axis,data_e1[i,:], '-x', label=data_labels[i])
    ax[1].plot(energy_axis,data_e2[i,:], '-x', label=data_labels[i])
    ax[0].legend()
    ax[1].legend()
    
    
    def Drude1Lorentz2_real(x, wpD, wpl1, taul1,wpl2, taul2):
#     wpD=0.31
    tauD=8
    wl1=0.64
    wl2=2.1
    return 5 - wpD**2/(x**2+tauD**(-2))+ (wpl1**2*(wl1**2-x**2))/((wl1**2-x**2)**2+x**2/taul1**2)+(wpl2**2*(wl2**2-x**2))/((wl2**2-x**2)**2+x**2/taul2**2)


def Drude1Lorentz2_imag(x,wpD, wpl1, taul1,wpl2, taul2):
#     wpD=0.31
    tauD=8
    wl1=0.64
    wl2=2.1
    return (1/(x*tauD))*wpD**2/(x**2+tauD**(-2))+ (wpl1**2*x/taul1)/((wl1**2-x**2)**2+x**2/taul1**2)+ (wpl2**2*x/taul2)/((wl2**2-x**2)**2+x**2/taul2**2)

def residuals(p,y,x):
    real = Drude1Lorentz2_real(x,*p)
    imag = Drude1Lorentz2_imag(x,*p)
    c = to_complex(real,imag)
    a = y - c
    return a.real ** 2 + a.imag ** 2
    
    
   parlables = ['wpD', 'wpl1', 'taul1','wpl2', 'taul2']#, 'wpl3', 'taul3']

guess = [0.0, .26, 3.17, 1.179, 2.47]
multiguess = [[0, 0.26, 3.17, 1.179, 2.47],
             [0.34, 0.26, 3.17, 1.179, 2.47],
             [0.34, 0.26, 3.17, 1.179, 2.47],
             [0.34, 0.26, 3.17, 1.179, 2.47],
             [0.34, 0.26, 3.17, 1.179, 2.47],
             [0.41, 0.26, 3.17, 1.179, 2.47],
             [0.41, 0.26, 3.17, 1.179, 2.47],
             [0.41, 0.26, 3.17, 1.179, 2.47],
             [0.41, 0.26, 3.17, 1.179, 2.47]]

bounds = [((-0.01,-np.inf,-np.inf,-np.inf,-np.inf),(0.3, np.inf, np.inf, np.inf, np.inf)),
          ((-np.inf,-np.inf,-np.inf,-np.inf,-np.inf),(np.inf, np.inf, np.inf, np.inf, np.inf)),
          ((-np.inf,-np.inf,-np.inf,-np.inf,-np.inf),(np.inf, np.inf, np.inf, np.inf, np.inf)),
          ((-np.inf,-np.inf,-np.inf,-np.inf,-np.inf),(np.inf, np.inf, np.inf, np.inf, np.inf)),
          ((-np.inf,-np.inf,-np.inf,-np.inf,-np.inf),(np.inf, np.inf, np.inf, np.inf, np.inf)),
          ((.40,-np.inf,-np.inf,-np.inf,-np.inf),(np.inf, np.inf, np.inf, np.inf, np.inf)),
          ((.40,-np.inf,-np.inf,-np.inf,-np.inf),(np.inf, np.inf, np.inf, np.inf, np.inf)),
          ((.40,-np.inf,-np.inf,-np.inf,-np.inf),(np.inf, np.inf, np.inf, np.inf, np.inf)),
          ((.40,-np.inf,-np.inf,-np.inf,-np.inf),(np.inf, np.inf, np.inf, np.inf, np.inf))]
# bounds = ((wpD_est*(1-wpD_tol),wpl1_est*(1-wpl1_tol),taul1_est*(1-taul1_tol), wpl2_est*(1-wpl2_tol),taul2_est*(1-taul2_tol)),
#            (wpD_est*(1+wpD_tol),wpl1_est*(1+wpl1_tol), taul1_est*(1+taul1_tol), wpl2_est*(1+wpl2_tol),taul2_est*(1+taul2_tol)))


n = 8
f, ax = plt.subplots(1,2,figsize=(10,5))
x_data = energy_axis
x_fit = np.linspace(.3,3,100)

fit_data_e1 = np.ndarray((100,n+1))
fit_data_e2 = np.ndarray((100,n+1))
fit_data_e1[:,0] = x_fit
fit_data_e2[:,0] = x_fit

fit_pars = {}
for par in parlables:
    fit_pars[par] = []

popt = guess
#variable n should be number of curves to plot (I skipped this earlier thinking that it is obvious when looking at picture - sorry my bad mistake xD): n=len(array_of_curves_to_plot)
#version 1:

color=cm.rainbow(np.linspace(0,1,n))
for i,c in zip(range(n),color):

    y_real = data_e1[i,:]
    y_imag = data_e2[i,:]
    z = to_complex(y_real,y_imag)
    ax[0].plot(x_data,y_real, 'o',c=c)
    ax[1].plot(x_data,y_imag, 'o',c=c)
    ax[0].set_ylabel('Epsilon_1')
    ax[0].set_xlabel('Energy [eV]')
    ax[1].set_ylabel('Epsilon_2')
    ax[1].set_xlabel('Energy [eV]')
    
    p0 = popt#multiguess[i]
    b = bounds[i]
    
#     popt, cov_x, infodict, mesg, ier = leastsq(residuals, p0, args=(z,x_data), full_output=1)
    res = least_squares(residuals, p0, args=(z,x_data), bounds=b)
    popt = res['x']
    for j, par in enumerate(popt):
        fit_pars[parlables[j]].append(par)
    
    fit_curve_e1 = (Drude1Lorentz2_real(x_fit,*popt))
    fit_curve_e2 = (Drude1Lorentz2_imag(x_fit,*popt))
    fit_data_e1[:,i+1] = fit_curve_e1
    fit_data_e2[:,i+1] = fit_curve_e2
    
    ax[0].plot(x_fit,fit_curve_e1,'-', c=c, label='{} mW'.format(data_labels[i]))
    ax[1].plot(x_fit,fit_curve_e2,'-', c=c, label='{} mW'.format(data_labels[i]))
    ax[0].legend()
    ax[1].legend()
#     res_percent_real = 100*(y_real-(Drude1Lorentz2_real(x_data,*popt)))/y_real
#     res_percent_imag = 100*(y_imag-(Drude1Lorentz2_imag(x_data,*popt)))/y_imag
#     ax[1,0].plot(x_data,res_percent_real,'-', c=c)
#     ax[1,1].plot(x_data,res_percent_imag,'-', c=c)    
#     ax[1,0].set_ylabel('%')
#     ax[1,1].set_ylabel('%')
# guess = []
# res = leastsq(residuals, p0, args=(z,Energy), full_output=1)
# popt, cov_x, infodict, mesg, ier  = res


f.savefig('e1_e2_values_best fit_till2_2ev.png')
np.savetxt('asd_e1.txt',fit_data_e1)
np.savetxt('asd_e2.txt',fit_data_e2)

f, ax = plt.subplots(2,3,figsize=(10,6))
i=0
for key, val in fit_pars.items():
    a,b = i//3,i%3
    ax[a,b].plot(data_labels[:len(val)],val, 'yo', c=c)
    ax[a,b].set_ylabel(key)
    ax[a,b].set_xlabel('Power [mW]')
    i+=1
f.savefig('e1_e2_fit_parameters_till2_2ev.png')

#---------------------------------------------------------------------------------------------------------------------------------------
def lorentz_real(w,s,w0,t):
    num = s**2 * (w0**2 - w**2)
    den = (w0**2 - w**2)**2 + w**2/t**-2
    return num/den
def lorentz_imag(w,s,w0,t):
    num = s**2 * w / t
    den = (w0**2 - w**2)**2 + w**2/t**-2
    return num/den
def de_lorentz_w0_real(w,s,w0,t):
    num = s**2*t**2*(-2*w0**5*t**2+4*w0**3*w**2*t**2+2*w0*w**2-2*w0*w**4*t**2)
    den = t**2*((w0**2-w**2)**2 + w**2)
    return num/den
def de_lorentz_w0_imag(w,s,w0,t):
    num = -(4*s**2*w*t**3*w0*(w0**2-w**2))
    den = (t**2*(w0**2 - w**2)**2 +w**2)**2
    return num/den
def de_lorentz_t_real(w,s,w0,t):
    num = 2*s**2*w**2*t(w0**2-w**2)
    den = (t**2*(w0**2-w**2)**2 + w**2)**2
    return num/den
def de_lorentz_t_imag(w,s,w0,t):
    num = s**2*w*(-t**2*w0**4+2*w**2*t**2*w0**2-w**4*y**2+w**2)
    den = (t**2*(w0**2-w**2)**2 + w**2)**2
    return num/den
    
    
 def delta_e1(w,C,Sd,Smir,Wct,Tct):
    Td = 8
    Wmir = 0.64
    Tmir = 'fix'
    Sct = fix
    Wct0 = 'fix from fit'
    Tct0 = 'fix from fit'
    eq = C
    eq += lorentz_real(w,Sd,0,Td)
    eq += lorentz_real(w,Smir,Wmir,Tmir)
    eq += de_lorentz_real(w,Sct,Wct,Tct)*(Wct-Wct0)
    eq += de_lorentz_real(w,Sct,Wct,Tct)*(Tct-Tct0)
    return eq
def delta_e2(w,C,Sd,Smir,Wct,Tct):
    Td = 'from literature'
    Wmir = fix
    Tmir = fix
    Sct = fix
    Wct0 = 'fix from fit'
    Tct0 = 'fix from fit'
    eq = 0
    eq += lorentz_imag(w,Sd,0,Td)
    eq += lorentz_imag(w,Smir,Wmir,Tmir)
    eq += de_lorentz_imag(w,Sct,Wct,Tct)*(Wct-Wct0)
    eq += de_lorentz_imag(w,Sct,Wct,Tct)*(Tct-Tct0)
    return eq

def residuals_delta(p,y,x):
    real = delta_e1(x,*p)
    imag = delta_e2(x,*p)
    c = to_complex(real,imag)
    a = y - c
    return a.real ** 2 + a.imag ** 2
#-------------------------------------------------------------------------------------------------
def lorentz_real(x,s,w0,t):
    num = s**2 * (w0**2 - x**2)
    den = (w0**2 - x**2)**2 + x**2/t**-2
    return num/den
def lorentz_imag(x,s,w0,t):
    num = s**2 * x / t
    den = (w0**2 - x**2)**2 + x**2/t**-2
    return num/den
def de_lorentz_w0_real(x,s,w0,t):
    num = s**2*t**2*(-2*w0**5*t**2+4*w0**3*x**2*t**2+2*w0*x**2-2*w0*x**4*t**2)
    den = t**2*((w0**2-x**2)**2 + x**2)
    return num/den
def de_lorentz_w0_imag(x,s,w0,t):
    num = -(4*s**2*x*t**3*w0*(w0**2-x**2))
    den = (t**2*(w0**2 - x**2)**2 +x**2)**2
    return num/den
def de_lorentz_t_real(x,s,w0,t):
    num = 2*s**2*x**2*t*(w0**2-x**2)
    den = (t**2*(w0**2-x**2)**2 + x**2)**2
    return num/den
def de_lorentz_t_imag(x,s,w0,t):
    num = s**2*x*(-t**2*w0**4+2*x**2*t**2*w0**2-x**4*t**2+x**2)
    den = (t**2*(w0**2-x**2)**2 + x**2)**2
    return num/den
def drude_real(x, s, t):
    num= s**2
    den = x**2 + 1/t**2
    return num/den
def drude_imag(x, s, t):
    num = s**2/(x*t)
    den = x**2 + 1/t**2
    return num/den


def delta_e1(x,Sd,Smir,Wct,Tct):
    Td = 8
    Wmir = 0.64
    Tmir = 3.1
    Sct = 0.53
    Wct0 = 2.3
    Tct0 = 2
    eq = 8
    eq = drude_real(x,Sd,Td)
    eq += lorentz_real(x,Smir,Wmir,Tmir)
    eq += de_lorentz_w0_real(x,Sct,Wct,Tct)*(Wct-Wct)
    eq += de_lorentz_t_real(x,Sct,Wct,Tct)*(Tct-Tct0)
    return eq
def delta_e2(x,Sd,Smir,Wct,Tct):
    Td = 8
    Wmir = 0.64
    Tmir = 3.1
    Sct = 0.53
    Wct0 = 2.3
    Tct0 = 2
    eq = 0
    eq += drude_imag(x,Sd,Td)
    eq += lorentz_imag(x,Smir,Wmir,Tmir)
    eq += de_lorentz_w0_imag(x,Sct,Wct,Tct)*(Wct-Wct0)
    eq += de_lorentz_t_imag(x,Sct,Wct,Tct)*(Tct-Tct0)
    return eq
    
    
    def residuals_delta(p,y,x):
    real = delta_e1(x,*p)
    imag = delta_e2(x,*p)
    c = to_complex(real,imag)
    a = y - c
    return a.real ** 2 + a.imag ** 2
    
    
    parlables = ['Sd', 'Smir', 'Wct','Tct']#, 'wpl3', 'taul3']

guess = [ 0.15757571,  1.81161059,  2.24411459,  0.97826626]
multiguess = [[ 0.15757571,  1.81161059,  2.24411459,  0.97826626],
             [ 0.15757571,  1.81161059,  2.24411459,  0.97826626],
             [ 0.15757571,  1.81161059,  2.24411459,  0.97826626],
             [ 0.15757571,  1.81161059,  2.24411459,  0.97826626],
             [ 0.15757571,  1.81161059,  2.24411459,  0.97826626],
             [ 0.15757571,  1.81161059,  2.24411459,  0.97826626],
             [ 0.15757571,  1.81161059,  2.24411459,  0.97826626]]
             

bounds = [((-np.inf,-np.inf,-np.inf,-np.inf),(np.inf, np.inf, np.inf, np.inf, )),
          ((-np.inf,-np.inf,-np.inf,-np.inf),(np.inf, np.inf, np.inf, np.inf)),
          ((-np.inf,-np.inf,-np.inf,-np.inf),(np.inf, np.inf, np.inf, np.inf)),
          ((-np.inf,-np.inf,-np.inf,-np.inf),(np.inf, np.inf, np.inf, np.inf)),
          ((-np.inf,-np.inf,-np.inf,-np.inf),(np.inf, np.inf, np.inf, np.inf)),
          ((-np.inf,-np.inf,-np.inf,-np.inf),(np.inf, np.inf, np.inf, np.inf)),
          ((-np.inf,-np.inf,-np.inf,-np.inf),(np.inf, np.inf, np.inf, np.inf))]
# bounds = ((wpD_est*(1-wpD_tol),wpl1_est*(1-wpl1_tol),taul1_est*(1-taul1_tol), wpl2_est*(1-wpl2_tol),taul2_est*(1-taul2_tol)),
#            (wpD_est*(1+wpD_tol),wpl1_est*(1+wpl1_tol), taul1_est*(1+taul1_tol), wpl2_est*(1+wpl2_tol),taul2_est*(1+taul2_tol)))


n = 4
f, ax = plt.subplots(1,2,figsize=(10,5))
x_data = energy_axis
x_fit = np.linspace(.3,3,100)

fit_data_de1 = np.ndarray((100,n+1))
fit_data_de2 = np.ndarray((100,n+1))
fit_data_de1[:,0] = x_fit
fit_data_de2[:,0] = x_fit

fit_pars = {}
for par in parlables:
    fit_pars[par] = []

popt = guess

#variable n should be number of curves to plot (I skipped this earlier thinking that it is obvious when looking at picture - sorry my bad mistake xD): n=len(array_of_curves_to_plot)
#version 1:

color=cm.rainbow(np.linspace(0,1,n))
for i,c in zip(range(n),color):

    y_real = data_de1[i,:]
    y_imag = data_de2[i,:]
    z = to_complex(y_real,y_imag)
    ax[0].plot(x_data,y_real, 'o',c=c)
    ax[1].plot(x_data,y_imag, 'o',c=c)
    ax[0].set_ylabel('Epsilon_1')
    ax[0].set_xlabel('Energy [eV]')
    ax[1].set_ylabel('Epsilon_2')
    ax[1].set_xlabel('Energy [eV]')
    
    p0 = multiguess[i]
    b = bounds[i]
    
#     popt, cov_x, infodict, mesg, ier = leastsq(residuals, p0, args=(z,x_data), full_output=1)
    res = least_squares(residuals_delta, p0, args=(z,x_data), bounds=b)
    popt = res['x']
    for j, par in enumerate(popt):
        fit_pars[parlables[j]].append(par)
    
    fit_curve_de1 = (delta_e1(x_fit,*popt))
    fit_curve_de2 = (delta_e2(x_fit,*popt))
    fit_data_de1[:,i+1] = fit_curve_de1
    fit_data_de2[:,i+1] = fit_curve_de2
    
    ax[0].plot(x_fit,fit_curve_de1,'-', c=c, label='{} mW'.format(data_labels[i]))
    ax[1].plot(x_fit,fit_curve_de2,'-', c=c, label='{} mW'.format(data_labels[i]))
    ax[0].legend()
    ax[1].legend()
#     res_percent_real = 100*(y_real-(Drude1Lorentz2_real(x_data,*popt)))/y_real
#     res_percent_imag = 100*(y_imag-(Drude1Lorentz2_imag(x_data,*popt)))/y_imag
#     ax[1,0].plot(x_data,res_percent_real,'-', c=c)
#     ax[1,1].plot(x_data,res_percent_imag,'-', c=c)    
#     ax[1,0].set_ylabel('%')
#     ax[1,1].set_ylabel('%')
# guess = []
# res = leastsq(residuals, p0, args=(z,Energy), full_output=1)
# popt, cov_x, infodict, mesg, ier  = res
print(popt)


f, ax = plt.subplots(2,3,figsize=(10,6))
i=0
for key, val in fit_pars.items():
    a,b = i//3,i%3
    ax[a,b].plot(data_labels[:len(val)],val, 'yo', c=c)
    ax[a,b].set_ylabel(key)
    ax[a,b].set_xlabel('Power [mW]')
    i+=1
#f.savefig('e1_e2_fit_parameters_till2_2ev.png')
