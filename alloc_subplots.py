# -*- coding: utf-8 -*-
"""
Created on Wed Dec 20 12:31:25 2023

@author: mathp
"""

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
Cnbiosynth=4.34728279914354E-10        #(molN cell-1 s) Constant for varible part of biosynthesis protein nitrogen (193-37)
Nstore_max=2.91679384515998E-15         #(molN cell-1) Constant protein pool in nitrogen (193-25)
Cnrna_variable=6212.59249917364        #(s) Constant for Variable part of RNA (193-26)
Ypthylakoid_chl=0.0281633095303638        #((molP cell-1)/(molC chl cell-1)) the shoichiometric ratio for cell phosphorus in thylakoid membrane to chlorophyll (193-26)
Pconst_other=5.44534485638617E-17               #(molP cell-1) Constant part of phosphorus (193-26) * This includes ATP ADP, Phospholipid and DNA RNA
Qp_max=25.26/(3.097e16)                                                                              #(molP cell-1) total phosphorus content in the cell (193-26)
Cessential=1.51786753491048E-15          #(molC cell-1) essential carbon (lipid membrane, etc.) *8.33e-14/10 is 10%
Nconst_protein = 4.45336898828389E-15   #(molN cell-1) Constant protein pool in nitrogen (193-25)

#==============================

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
Feunit = 1/Qc
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


Color_DNA_const='orange'
Color_RNA_const='yellow'
Color_protein_const='blue'
Color_photo='yellow'
Color_RNA_variable='red'
Color_DNA_variable='blue'
Color_protein_biosynthesis='pink'
Color_chl='#66FF66'
Color_other='#006400'
Color_other='#008000'
Color_P_const='#CCFFFF'
Color_DNA='black'
Color_RNA='red'
Color_Nstore='purple'
Color_Pstore='#D0CECE'
Color_Cessential='brown'
Color_thylakoid='#FFD966'

data = pandas.read_csv("C:\\Users\\mathp\\OneDrive\\Desktop\\FE_opt_Spyd\\Optimization\\fe_all_new.csv",delimiter=',', index_col ="author")
df = pandas.DataFrame(data)
df["Fe (uM)"] = df["Fe (uM)"].astype(float)

# namesforloop = ["Sunda (1995) hux", "Sunda (1995) oce","Sunda (1995) cal","Sunda (1997) pse","Sunda (1997) wei",\
#                     "Sunda (1997) mic","Sunda (1997) min","Jabre (2020) 6","Jabre (2020) 3","Jabre (2020) 1"]

namesforloop = ["Sunda (1995) hux","Sunda (1995) cal","Sunda (1997) mic", "Sunda (1997) min","Sunda (1995) oce","Sunda (1997) pse","Sunda (1997) wei",\
                "Jabre (2020) 6","Jabre (2020) 3","Jabre (2020) 1"]
    
# namesforloop = ["Sunda (1995) hux","Sunda (1995) cal","Sunda (1997) mic", "Sunda (1997) min","Sunda (1995) oce","Sunda (1997) pse","Sunda (1997) wei",\
#                     "Jabre (2020) 6"]


    
fig, axs = plt.subplots(nrows=2, ncols=5, figsize=(21, 8))
# fig, axs = plt.subplots(nrows=2, ncols=4, figsize=(16, 8))
plt.subplots_adjust(left=0.12, bottom=0.12, right=None, top=None,wspace=None, hspace=None)
afebest = pandas.read_csv("C:\\Users\\mathp\\OneDrive\\Desktop\\FE_opt_Spyd\\Optimization\\afebest_light_dark.csv",header=None)
ynbest = pandas.read_csv("C:\\Users\\mathp\\OneDrive\\Desktop\\FE_opt_Spyd\\Optimization\\ynbest_light_dark.csv",header=None)
    
axs = axs.ravel()
for j in range(np.size(namesforloop)):
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
    rows = df.loc[author_name] 
    Pchl=Pmax*(1-exp(-OT*I))
    Fe =rows[["Fe (uM)"]].values
    if author_name == "Jabre (2020) 1" or author_name == "Jabre (2020) 3" or author_name == "Jabre (2020) 6":
        ld_cycle = 1
    else:
        ld_cycle =     0.71428571428
        
    Pchl = Pchl*ld_cycle

    aFe = afebest.loc[j,0]
    Ynphoto_chl = ynbest.loc[j,0]
    fe_max = numpy.nanmax(Fe) + 1e-4
    Fe_high_res = arange(1e-7,fe_max,1e-6) 
    Numbertoarray=ones(size(Fe_high_res))     
    Vfe=aFe*Fe_high_res 
    Pchl=Pmax*(1-exp(-OT*I))
    A=((1+E)*Qc*Ynphoto_chl)/Pchl+Cnbiosynth
    B=Nconst_protein+(m*Ynphoto_chl)/Pchl
    L = ((1 + E)*Qc*Ypthylakoid_chl)/Pchl
    M = (m*Ypthylakoid_chl)/Pchl
    R=((1+E)*Qc*Ynphoto_chl*YphotoFe_N)/Pchl
    S=(m*Ynphoto_chl*YphotoFe_N)/Pchl
    
    aFe = R
    bFe = S
    cFe = -Vfe
    
    DFe=DSolver(aFe,bFe,cFe)
    D23=DFe.rQ
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
    
    Qc_essential = Qc_essential_computation(D23,Chl_const, Chl_D, Nchl_const, Nchl_D, Nphoto_const, Nphoto_D, Nbiosynth_D,\
                                                  Nprotein_const, Nprotein_D, Nrna_const, Nrna_D, Nrna_D2, Ypthylakoid_chl, YrnaC_N,\
                                                  CNprotein, YpgC_P, Cdna_const, Cconst_protein, Cessential)  
    
    i0 = Qc_essential >= Qc   #conditions where Mu max applies
    i1 = isnan(Qc_essential)  #Including this avoide when Qc_essential is nan with high Vn
    i = i0 + i1
    D23[i] = 2*(Qc - Qc_const)/(Qc_D + sqrt(Qc_D*Qc_D + 4*Qc_D2*(Qc - Qc_const)))
    inp=D23
    Numtoarray=np.ones(np.shape(inp)) 
    Chl       = Chl_const      + Chl_D*inp# masked aray, ~4.4 e-11
    Nchl      = Nchl_const     + Nchl_D*inp #masked array, ~3.2 e-12             #(molN chl cell-1) Chlorophyll N concentration
    Nphoto    = Nphoto_const   + Nphoto_D*inp #masked array, ~5.6 e-10    # Nphoto_const should be like e-16  #(molN cell-1) Photosynthesis related protein nitrogen (193-25)6
    Nbiosynth =                  Nbiosynth_D*inp# masked array, ~1.7 e-10
    Nprotein  = Nprotein_const + Nprotein_D*inp  #masked array, ~7.4 e-10     #(mol N cell-1) all the protein
    Nrna      = Nrna_const     + Nrna_D*inp + Nrna_D2*inp*inp#masked array, 1.8 ~e-6
    Nessential= Nchl + Nphoto + Nbiosynth + Nconst_protein + Nrna + Ndna #masked array, ~1.8 e -6
     # Qn_max    = Nessential + Nstore_max #masked array, ~1.8 e -6
    
      #Pdna=Ndna*YnucacidP_N #7.06e-18  #(mol P cell-1) DNA in phosphorus
    Pthylakoid=Chl*Ypthylakoid_chl #masked array, 1.2 e-12         #(molP cell-1) Phosphorus in thylakoid membranes: phospholipid, etc. (193-26)
      #Prna=Nrna*YnucacidP_N #masked array, 4.9 e-7    #(molP cell-1) Phosphorus in RNA 
      #Pessential=Pthylakoid+Prna+Pdna+Pconst_other #masked array, 4.9 e-7   #(molP cell-1)
     
    Cchl = Chl   #masked array, 4.4 e-11                 #(molC cell-1) carbon in chlorophyll (195-16)
    Cphoto = Nphoto*CNprotein #masked array, 2.5 e-9 #Nphoto should be like e-15     #(molC cell-1) carbon in photosystem protein (195-16)
    Cbiosynth = Nbiosynth*CNprotein #masked array, 7.99 e-10  #(molC cell-1) carbon in biosynthesis protein (195-16)
    Crna = Nrna*YrnaC_N #masked array, 5.3 e-6     #(molC cell-1) carbon RNA (195-16)
    CthylakoidPG = Pthylakoid*YpgC_P  #masked array, 5.0 e-11         #(molC cell-1) carbon in PG (phosphatidyl glycerol) in thylakoid membranes
    Cnstore = np.zeros(np.shape(inp))#Nstore*YcyanoC_N
      
    Cother = Qc - Cphoto - Cbiosynth - Cconst_protein - Cchl\
              - Crna - Cdna_const\
              - Cessential - CthylakoidPG - Cnstore #masked array, -5.2 e-6
    percentorratio=100       #100: percent, 1:ratio
    Cphoto_plot=Cphoto/Qc*percentorratio   #masked array, 3 e+6  #Cphoto should be like e -14/e -15       
    Cbiosynth_plot=Cbiosynth/Qc*percentorratio #masked array, 9e+5
    Cconst_protein_plot=Cconst_protein/Qc*percentorratio*Numtoarray #not masked array, 23 (all same value)
    Cchl_plot=Cchl/Qc*percentorratio #masked array, 5 e+4
    Crna_plot=Crna/Qc*percentorratio*Numtoarray #masked array, 6.4 e+9
    Cdna_const_plot=Cdna_const/Qc*percentorratio*Numtoarray #not masked array, 0.09 (all same)
    Cother_plot=Cother/Qc*percentorratio #masked array, -6.4 e+9
    Cessential_plot=Cessential/Qc*percentorratio*Numtoarray  #not masked array, 1.82 (all same)
    CthylakoidPG_plot=CthylakoidPG/Qc*percentorratio #masked array, 6 e+4
    Cnstore_plot=Cnstore/Qc*percentorratio #not masked array, all zeros
    C_Photo = Cphoto_plot  + CthylakoidPG_plot  + Cchl_plot 
    C_Bio = Crna_plot + Cbiosynth_plot
    C_Other = Cessential_plot + Cconst_protein_plot + Cdna_const_plot
    
    #======================
    Color_Photo='#ffdd74'   # bright orange 
    Color_Bio='#b3ffb3'     #bright green
    Color_RNA='#99FFCC'
    Color_Other='#85ccff'  # dark blue #change
    Color_Nstore='#FFCCCC'  # pink (doesn't show)
    Color_other='#FFFF99' #light yellow
    y1 = Cessential_plot+Cconst_protein_plot+Cdna_const_plot
    y2 = Cphoto_plot+CthylakoidPG_plot+Cchl_plot
    y3 = Crna_plot+Cbiosynth_plot
    y4 = Cother_plot
    StackPlotColorsC2=(Color_Other,Color_Photo,Color_Bio,Color_other)
    # axs[j].stackplot(Fe_high_res*1000,Cessential_plot+Cconst_protein_plot+Cdna_const_plot,\
    #                   Cphoto_plot+CthylakoidPG_plot+Cchl_plot,Crna_plot+Cbiosynth_plot,Cnstore_plot,Cother_plot,colors=StackPlotColorsC2)
    axs[j].stackplot(Fe_high_res*1000,y1,y2,y3,y4,colors=StackPlotColorsC2)
    
    A_FePho = 5.83*10e-3
    Q_Fe = (A_FePho*Chl_D*inp)+(A_FePho*Chl_const)

    # axs[j].stackplot(Fe_high_res*1000,Q_Fe*Feunit,colors=StackPlotColorsC2)
    
    # axs[j].plot(Fe_high_res*1000,y1,'#85ccff')
    # axs[j].plot(Fe_high_res*1000,y2,'#ffdd74' )
    # axs[j].plot(Fe_high_res*1000,y3,'#b3ffb3')
    # axs[j].plot(Fe_high_res*1000,y4,'#FFCCCC')

    axs[j].set_xlim(left=0,right=fe_max*1000)
    
    axs[j].set_ylim(bottom=0,top=100)

    
fss = 24
fig.text(0.5, 0.005, 'Fe (nmol L$^{-1}$)', ha='center')
fig.text(0.005, 0.5, 'C allocation ($\%$)', va='center', rotation='vertical')
# fig.text(0.001, 0.5, 'Fe:C (mol mol$^{-1}$)', va='center', rotation='vertical')
# fig.text("P:C (mol mol$^{-1}$)")

# for 10 subplots
fig.text(0.162,0.62,'a',fontsize=fss)
fig.text(0.36,0.62,'b',fontsize=fss)
fig.text(0.562,0.62,'c',fontsize=fss)
fig.text(0.76,0.62,'d',fontsize=fss)
fig.text(0.956,0.62,'e',fontsize=fss)
fig.text(0.162,0.131,'f',fontsize=fss)
fig.text(0.36,0.131,'g',fontsize=fss)
fig.text(0.562,0.131,'h',fontsize=fss)
fig.text(0.76,0.131,'i',fontsize=fss)
fig.text(0.956,0.131,'j',fontsize=fss)


# for 8 subplots
# fig.text(0.2,0.61,'A.',fontsize=fss)
# fig.text(0.452,0.61,'B.',fontsize=fss)
# fig.text(0.695,0.61,'C.',fontsize=fss)
# fig.text(0.94,0.61,'D.',fontsize=fss)
# fig.text(0.2,0.131,'E.',fontsize=fss)
# fig.text(0.452,0.131,'F.',fontsize=fss)
# fig.text(0.695,0.131,'G.',fontsize=fss)
# fig.text(0.94,0.131,'H.',fontsize=fss)

# fig.legend(['Essential', 'Photosynthesis', 'Biosynthesis','Storage'],loc="lower left",ncol=4,fontsize="x-small")
# fig.legend(loc="upper left")
# sf_opt("subplots_allocation_legend")
# sf_opt("all_figs_alloc_light_new_xlims")
sf_opt("fe_allocation_plot_620")

