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




