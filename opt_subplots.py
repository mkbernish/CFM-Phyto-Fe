'''
From a800_05_12_07_SameQc.py
'''

import sys
sys.path.append("C:\\Users\\mathp\\OneDrive\\Desktop\\FE_opt_Spyd\\Optimization")
sys.path.append("C:\\Users\\mathp\\OneDrive\\Desktop\\FE_opt_Spyd")
from pylab import *
from sf import *
from sf_opt import *
import numpy
from af001_energy_calculation import *
from Solver_2D import *
from Solver_3D import *
from Solver_2D_O import *
from Qc_essential_computation import *
from FigSetting2 import *
import random
import time
import pandas

#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
# Parameter sets (preparation)
#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

m=3.79146798299876E-19         #(mol C s-1 cell-1) maintenance carbonhydrate consumption (idea from 172-7)
Pmax=0.00320513285659728

OT=0.00863364097132997
# Cnbiosynth=4.34728279914354E-10        #(molN cell-1 s) Constant for varible part of biosynthesis protein nitrogen (193-37)
Nstore_max=2.91679384515998E-15         #(molN cell-1) Constant protein pool in nitrogen (193-25)
# Cnrna_variable=6212.59249917364        #(s) Constant for Variable part of RNA (193-26)
Ypthylakoid_chl=0.0281633095303638        #((molP cell-1)/(molC chl cell-1)) the shoichiometric ratio for cell phosphorus in thylakoid membrane to chlorophyll (193-26)
Pconst_other=5.44534485638617E-17               #(molP cell-1) Constant part of phosphorus (193-26) * This includes ATP ADP, Phospholipid and DNA RNA
Qp_max=25.26/(3.097e16)                                                                              #(molP cell-1) total phosphorus content in the cell (193-26)
Cessential=1.51786753491048E-15          #(molC cell-1) essential carbon (lipid membrane, etc.) *8.33e-14/10 is 10%
Nconst_protein = 4.45336898828389E-15   #(molN cell-1) Constant protein pool in nitrogen (193-25)



E3=evalue()
E=E3.E
Qc=1.00*10**(-12)/12      #(molC/cell) biomass C per cell (196-18)(average of N and P limited cases from Healey 1985)
YchlN_C=4/55

#Conversion parameters================
CNprotein=4.49   #(molC molN) the ratio of C to N in protein (derived from Brown 1991) calculation in "13 Amino acid composition of different phytoplankton.xlsx"
YcyanoC_N=2                             #(molC molN) C/N molar ratio of cyanophycin
YpgC_P=40                           #(molC molP) C/P molar ratio of PG: Phosphatidyl glycerol (assuming C16 fatty acids (since actually mostly C16 (Huflejt et al., 1990)
Nunit=1/Qc#*14*10**6/(12*10**3)          #((ug N / mgC)/(molN cell-1) unit conversion term (164-20)
Punit=1/Qc#*30.97*10**6/(12*10**3)       #((ug P / mgC)/(molP cell-1) unit conversion term (164-20)

#=========================================
#DNA and RNA-> Kei 193-28
#=========================================              

Molar_mass_DNA_AT_average=307.47        #(g mol-1) Molar mass AT average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")
Molar_mass_DNA_CG_average=307.97        #(g mol-1) Molar mass CG average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")
Molar_mass_RNA_AT_average=316.47        #(g mol-1) Molar mass AT average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")
Molar_mass_RNA_CG_average=323.97        #(g mol-1) Molar mass CG average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")

#------------------------------------
#E coli
#------------------------------------
CG_Ecoli=0.506          #(dimensionless) from [http://www.ncbi.nlm.nih.gov/genome/167 (accessed 06/18/2016)]
AT_Ecoli=1-CG_Ecoli     #(dimensionless) 

Molar_mass_DNA_Ecoli=Molar_mass_DNA_AT_average*CG_Ecoli+Molar_mass_DNA_CG_average*AT_Ecoli     #(g mol-1) Molar mass of DNA unit
Molar_mass_RNA_Ecoli=Molar_mass_RNA_AT_average*CG_Ecoli+Molar_mass_RNA_CG_average*AT_Ecoli     #(g mol-1) Molar mass of RNA unit

RNA_DNA_mass_ratio=17.844/6.5239  #(ug/ug) from values ad D=0 "07 Bremer and Dennis 1996 data plot.xlsx"
RNA_DNA_molar_ratio=RNA_DNA_mass_ratio/Molar_mass_RNA_Ecoli*Molar_mass_DNA_Ecoli    #(mol mol-1)

#---------------------------------------------
#Stoichiometric parameters for DNA and RNA
#---------------------------------------------

CG=0.563                   #GC% not CG but I started with CG so I stick with it; it does not matter as "AT GC".   [http://www.ncbi.nlm.nih.gov/genome/13522 (accessed 06/18/2016)]
YnucacidP_N=1/(3.5*(1-CG)+4*CG)               #(molP molN-1) P/N molar ratio of RNA (193-26) values (193-28) excel file "08 N to P ratio in DNA and RNA.xlsx"

YdnaC_N=3.5*(1-CG)+2.5*CG       #(molC molN-1) C/N molar ratio of dna (based on "10 Amino acids and nucleic acids stoichiometry.xlsx)
YrnaC_N=3.25*(1-CG)+2.5*CG      #(molC molN-1) C/N molar ratio of rna (based on "10 Amino acids and nucleic acids stoichiometry.xlsx)

DNAmb=2.1269                   #(Mb) Megabase pair of synechococcus DNA in mega (million) base pairs [http://www.ncbi.nlm.nih.gov/genome/13522 (accessed 06/18/2016)]
Avogadro=6.022*10**23           #(molecules mol-1) Avogadro constant
Pdna_const=DNAmb*2*10**6/Avogadro                #(molP cell-1) Constant part of DNA in phosphorus 
Prna_const=Pdna_const*RNA_DNA_molar_ratio       #(molP cell-1) Constant part of RNA in phosphorus
#* Make sure to multiply by 2 as they are base PAIRs"
Ndna_const=Pdna_const/YnucacidP_N      #(molN cell-1) Constant part of DNA in nitrogen
Nrna_const=Ndna_const*RNA_DNA_molar_ratio   #(molN cell-1) Constatn part of RNA in phosphorus
Ndna=Ndna_const    #(molN cell-1) DNA in nitrogen (here assuming constant)
#Ynphoto_chl=3.56099164557551 
YphotoFe_N=0.001636364  #(molFe molN-1) Fe/N ratio in photosystem iron (0.001636364 in Inomura et al., 2022)

def function(aFe,Ynphoto_chl,Fe,I_s,author_name):
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #Main calculation 199-21~
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    Pchl=Pmax*(1-exp(-OT*I_s)) #(C mol s-1 Chl mol-1) Carbohydrate fixation rate per chlorophyll (167-1)(193-25)
    if author_name == "Jabre (2020) 1" or author_name == "Jabre (2020) 3" or author_name == "Jabre (2020) 6":
        ld_cycle = 1
    else:
        ld_cycle =     0.71428571428
        
    Pchl = Pchl*ld_cycle
    Vfe=aFe*Fe    #(mol P cell-1 s-1) iron uptake per cell
    # A_chl = (1+E)/Pchl
    # B_chl = (m/Pchl)
    # Apho_fe = 5.83e-3
    A=((1+E)*Qc*Ynphoto_chl)/Pchl+Cnbiosynth
    B=Nconst_protein+(m*Ynphoto_chl)/Pchl
    L = ((1 + E)*Qc*Ypthylakoid_chl)/Pchl
    M = (m*Ypthylakoid_chl)/Pchl
    #========================
    # Fe limitation related
    #========================
    
    R=((1+E)*Qc*Ynphoto_chl*YphotoFe_N)/Pchl
    S=(m*Qc*Ynphoto_chl*YphotoFe_N)/Pchl
    
    aFe = R
    bFe = S
    cFe = -Vfe
    
    DFe=DSolver(aFe,bFe,cFe)
    D=DFe.rQ
    # DFe=solver_2D(aFe,bFe,cFe)

    Chl_const = m/Pchl                              # (molC chl cell-1) cN[i]hlrophyll concentration (193-25)
    Chl_D = (1 + E)*Qc/Pchl                           # (molC chl cell-1) cN[i]hlrophyll concentration (193-25)
    
    #Nitrogen related
    Nchl_D = Chl_D*YchlN_C                          # (molN chl cell-1) Chlorophyll N concentration
    Nphoto_D = Chl_D*Ynphoto_chl                    # (molN cell-1) Photosynthesis related protein nitrogen (193-25)
    Nchl_const = Chl_const*YchlN_C                  # (molN chl cell-1) Chlorophyll N concentration
    Nphoto_const = Chl_const*Ynphoto_chl            # (molN cell-1) Photosynthesis related protein nitrogen (193-25)
    Nbiosynth_D = Cnbiosynth                        # (molN cell-1) various part of biosynthesis related protein in N (193-37)
    Nprotein_D = Nphoto_D + Nbiosynth_D             # (molN cell-1) All the proteins in N (193-26)
    Nprotein_const = Nphoto_const + Nconst_protein  # (molN cell-1) All the proteins in N (193-26)
    Nrna_D = Nprotein_const*Cnrna_variable          # (molN cell-1) variable part of nitrogen in RNA (193-26)(193-37)
    Nrna_D2 = Nprotein_D*Cnrna_variable             # (molN cell-1) variable part of nitrogen in RNA (193-26)(193-37)
    
    #Constant carbon parameters------------------
    Cconst_protein = Nconst_protein*CNprotein  #(molC cell-1) carbon in other protein assumed constant (195-16)
    Cdna_const = Ndna_const*YdnaC_N      #(molC cell-1) carbon in constant part of DNA (195-16)   
    
    #Calculating factors and Qc_essential-------------
    
    Qc_D2 = A*Cnrna_variable*YrnaC_N
    Qc_D = (1+E)*Qc/Pchl + A*CNprotein + L*YpgC_P + B*Cnrna_variable*YrnaC_N
    Qc_const = m/Pchl + B*CNprotein + M*YpgC_P + Nrna_const*YrnaC_N + Ndna_const*YdnaC_N + Cessential
    
    Qc_essential = Qc_essential_computation(D,Chl_const, Chl_D, Nchl_const, Nchl_D, Nphoto_const, Nphoto_D, Nbiosynth_D,\
                                                  Nprotein_const, Nprotein_D, Nrna_const, Nrna_D, Nrna_D2, Ypthylakoid_chl, YrnaC_N,\
                                                  CNprotein, YpgC_P, Cdna_const, Cconst_protein, Cessential)  
    
    i0 = Qc_essential >= Qc   #conditions where Mu max applies
    i1 = isnan(Qc_essential)  #Including this avoide when Qc_essential is nan with high Vn
    i = i0 + i1
    D[i] = 2*(Qc - Qc_const)/(Qc_D + sqrt(Qc_D*Qc_D + 4*Qc_D2*(Qc - Qc_const)))
    #D[i] = 2*(Qc - Qc_const)/(Qc_D + sqrt(Qc_D*Qc_D + 4*Qc_D2*(Qc - Qc_const)))
    #DSolver(Qc_D2,Qc_D,Qc_const)(0)
    return D
# function ending here

#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
# AW method: Largely copied from J001 00 00
#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO


data = pandas.read_csv("C:\\Users\\mathp\\OneDrive\\Desktop\\FE_opt_Spyd\\Optimization\\fe_all_new.csv",delimiter=',', index_col ="author")
# "C:\Users\mathp\OneDrive\Desktop\FE_opt_Spyd\Optimization\fe_all_new.csv"
# Iron data is converted to micromolar [dissolved Fe]
df = pandas.DataFrame(data) 
df["Fe (uM)"] = df["Fe (uM)"].astype(float)

namesforloop = ["Sunda (1995) hux","Sunda (1995) cal","Sunda (1997) mic", "Sunda (1997) min","Sunda (1995) oce","Sunda (1997) pse","Sunda (1997) wei",\
                "Jabre (2020) 6","Jabre (2020) 3","Jabre (2020) 1"]
    
# namesforloop = ["Sunda (1995) hux","Sunda (1995) cal","Sunda (1997) mic", "Sunda (1997) min","Sunda (1995) oce","Sunda (1997) pse","Sunda (1997) wei",\
                    # "Jabre (2020) 6"]

    
fig, axs = plt.subplots(nrows=2, ncols=5, figsize=(20, 8))
# fig, axs = plt.subplots(nrows=2, ncols=4, figsize=(16, 8))
# plt.subplots_adjust(hspace=0.11,wspace=0.05)
plt.subplots_adjust(left=0.14, bottom=0.12, right=0.15, top=None,wspace=None, hspace=None)

axs = axs.ravel()

afe_best_array = np.zeros([np.size(namesforloop),1])
yn_best_array = np.zeros([np.size(namesforloop),1])
station_fe = np.zeros([np.size(namesforloop),1])
for j in range(10):
    author_name = namesforloop[j]
    Ea = 70000 #activation energy
    R = 8.3
    A =Ea/R
    #Tt - ambient temperature that you are reading in
    if author_name == "Jabre (2020) 1":
        Tt = 1 + 273.15
    elif author_name == "Jabre (2020) 3":
        Tt = 3 + 273.15
    elif author_name == "Jabre (2020) 6":
        Tt = 6 + 273.15
    else:
        Tt = 20 + 273.15

    Tref = 293
    Arr = exp(-A*((1/Tt)-(1/Tref)))
    
    Cnbiosynth=4.34728279914354E-10/Arr        #(molN cell-1 s) Constant for varible part of biosynthesis protein nitrogen (193-37)
    Cnrna_variable=6212.59249917364/Arr       #(s) Constant for Variable part of RNA (193-26)

    if author_name == "Sunda (1995) hux" or author_name == "Sunda (1995) oce" or author_name == "Sunda (1995) cal" or author_name == "Sunda (1997) pse" or \
        author_name == "Sunda (1997) wei":
        I = 500
    elif author_name == "Sunda (1997) mic" or author_name == "Sunda (1997) min" or author_name == "Jabre (2020) 6" or author_name == "Jabre (2020) 3" or \
        author_name == "Jabre (2020) 1" or author_name == "Hudson (1990) br":
        I = 50
    elif author_name == "Hudson (1990)":
        I = 100
    elif author_name == "Sunda (1995) hux 175" or author_name == "Sunda (1995) oce 175":
        I = 175
        

    rows = df.loc[author_name] 
    Fe =rows[["Fe (uM)"]].values
    Mu_d = rows[["mu"]].values
    
    #Removing negative values
    Fe[Mu_d<0] = nan
    Mu_d[Mu_d<0] = nan
    
    #Standard deviation of samples
    sigma = nanstd(Mu_d)*sqrt(size(Mu_d)-sum(isnan(Mu_d)))/sqrt(size(Mu_d)-sum(isnan(Mu_d))-1)
    
    
    #print(sigma)
    
    # Initial values
    if author_name == "Sunda (1995) hux 175":
        aFe = 1*Qc/86400
    else:
        aFe = 0.1*Qc/86400  #aFe #(mol N cell-1 s-1 / umolN L-1) affinity of Fe (initial value from d007_06_00)
    Yn = 2
    repeat_times = arange(100000+1)
    n = zeros(size(repeat_times))*nan
    aFe_array = copy(n)
    Yn_array = copy(n)
    Nc_array = copy(n)
    
    #change ratio per iteration
    a = 0.5
    a_aFe = a*aFe
    a_Yn = a*Yn
    
    #For step 0
    Mu0 = function(aFe,Yn,Fe,I,author_name)*86400
    
    X2 = nansum((Mu_d - Mu0)**2/(2*sigma**2))
    
    P = exp(-X2)
    Pbest = P
    aFebest = aFe
    Ynbest = Yn
    
    #aPbest= 6.247661717797246e-19
    t0 = time.time()
    for i in repeat_times:
    
        #checking if the parameter value come within the realistic range
        while True:
            aFe1 = aFe + a_aFe*random.uniform(-1,1)
            if aFe1>0 and aFe1<1e-16:
                break
        
        while True:
            Yn1 = Yn + a_Yn*random.uniform(-1,1)
            if Yn1>0 and Yn1<70:
                break
        
        #Candidate value calculation
        Mu1 = function(aFe1,Yn1,Fe,I,author_name)*86400
        X2 = nansum((Mu_d - Mu1)**2/(2*sigma**2))
        P1 = exp(-X2)
        Pratio = P1/P
        r01 = random.uniform(0,1)
    
        #testing the quality of Pratio
        if r01 < Pratio:
            P = P1
            aFe = aFe1
            Yn = Yn1
            
            #Pbest update
            if P>Pbest:
                Pbest = P
                aFebest = aFe
                Ynbest = Yn
        
        #Recording values
        aFe_array[i] = aFe
        Yn_array[i] = Yn
        
        #Time printing
        if mod(i,10000) == 0:
            print(i,round(time.time()-t0,2),'(s)')
            t0 = time.time()
    print('Pbest=', Pbest,', aFebest=', aFebest, ', Ynbest=', Ynbest)
    fe_max = numpy.nanmax(Fe) + 1e-4
    afe_best_array[j] = aFebest
    yn_best_array[j] = Ynbest
    Fe_high_res = arange(1e-7,fe_max,1e-7) 
    Numbertoarray=ones(size(Fe_high_res))     
    Mu = function(aFebest,Ynbest,Fe_high_res,I,author_name)
    # Mu = function(1e-18*.7,1,Fe_high_res)
    D = Mu

    axs[j].plot(Fe_high_res*1000,D*86400,'k')
    axs[j].plot(Fe*1000,Mu_d,'ro')
    axs[j].set_xlim(left=0,right=fe_max*1000)


    def max_number_index(input_list):
        idx = 0
        max_number = input_list[0]
        for i in range(1, len(input_list)):
            if input_list[i] > max_number:
                max_number = input_list[i]
                idx = i
        return idx
    idx_max = max_number_index(D)
    station_fe[j]  = Fe_high_res[idx_max]

df1=pandas.DataFrame(afe_best_array)
df2=pandas.DataFrame(yn_best_array)
# df1.to_csv("C:\\Users\\mathp\\OneDrive\\Desktop\\FE_opt_Spyd\\Optimization\\afebest_light_dark.csv",index=False,header=False)
# df2.to_csv("C:\\Users\\mathp\\OneDrive\\Desktop\\FE_opt_Spyd\\Optimization\\ynbest_light_dark.csv",index=False,header=False)

fss = 24
fig.text(0.5, 0.005, 'Fe (nmol L$^{-1}$)', ha='center')
fig.text(0.001, 0.5, '$\mathit{\mu}$ (d$^{-1}$)', va='center', rotation='vertical')
# for 10 subplots
fig.text(0.16,0.62,'a',fontsize=fss)
fig.text(0.36,0.62,'b',fontsize=fss)
fig.text(0.562,0.62,'c',fontsize=fss)
fig.text(0.754,0.62,'d',fontsize=fss)
fig.text(0.956,0.62,'e',fontsize=fss)
fig.text(0.16,0.131,'f',fontsize=fss)
fig.text(0.36,0.131,'g',fontsize=fss)
fig.text(0.562,0.131,'h',fontsize=fss)
fig.text(0.754,0.131,'i',fontsize=fss)
fig.text(0.956,0.131,'j',fontsize=fss)
# fig.text(0.709,0.094,'K.',fontsize=fss)
# fig.text(0.956,0.094,'L.',fontsize=fss)

# for 8 subplots
# fig.text(0.2,0.61,'a',fontsize=fss)
# fig.text(0.452,0.61,'b',fontsize=fss)
# fig.text(0.695,0.61,'c',fontsize=fss)
# fig.text(0.94,0.61,'d',fontsize=fss)
# fig.text(0.2,0.131,'e',fontsize=fss)
# fig.text(0.452,0.131,'f',fontsize=fss)
# fig.text(0.695,0.131,'g',fontsize=fss)
# fig.text(0.94,0.131,'h',fontsize=fss)
# fig.text(0.754,0.131,'I.',fontsize=fss)
# fig.text(0.956,0.131,'J.',fontsize=fss)
# fig.text(0.709,0.094,'K.',fontsize=fss)
# fig.text(0.956,0.094,'L.',fontsize=fss)
sf_opt("subplots_opt_6_18")


