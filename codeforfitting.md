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
