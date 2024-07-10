# -*- coding: utf-8 -*-
"""
Created on Wed Dec 27 13:03:31 2023

@author: mathp
"""
import sys
sys.path.append("C:\\Users\\mathp\\OneDrive\\Desktop\\FE_opt_Spyd\\Optimization")
sys.path.append("C:\\Users\\mathp\\OneDrive\\Desktop\\FE_opt_Spyd")
from pylab import *
from sf import *
from sf_opt import *
import numpy as np
from af001_energy_calculation import *
from Solver_2D import *
from Solver_3D import *
from Solver_2D_O import *
from Qc_essential_computation import *
from FigSetting2 import *
import random
import time
import pandas
import matplotlib.pyplot as plt
import netCDF4 as nc
import cartopy
import cartopy.crs as ccrs
import cartopy.mpl.ticker as cticker
import matplotlib.cm as cm
from light_with_chla import *

m=3.79146798299876E-19         #(mol C s-1 cell-1) maintenance carbonhydrate consumption (idea from 172-7)
Pmax=0.00320513285659728

OT=0.00863364097132997
Nstore_max=2.91679384515998E-15         #(molN cell-1) Constant protein pool in nitrogen (193-25)
Ypthylakoid_chl=0.0281633095303638        #((molP cell-1)/(molC chl cell-1)) the shoichiometric ratio for cell phosphorus in thylakoid membrane to chlorophyll (193-26)
Pconst_other=5.44534485638617E-17               #(molP cell-1) Constant part of phosphorus (193-26) * This includes ATP ADP, Phospholipid and DNA RNA
Qp_max=25.26/(3.097e16)                                                                              #(molP cell-1) total phosphorus content in the cell (193-26)
Cessential=1.51786753491048E-15          #(molC cell-1) essential carbon (lipid membrane, etc.) *8.33e-14/10 is 10%
Nconst_protein = 4.45336898828389E-15   #(molN cell-1) Constant protein pool in nitrogen (193-25)
E3=evalue()
E=E3.E
Qc=1.00*10**(-12)/12      #(molC/cell) biomass C per cell (196-18)(average of N and P limited cases from Healey 1985)
YchlN_C=4/55
depth_to_z = 5
#Conversion parameters================
CNprotein=4.49   #(molC molN) the ratio of C to N in protein (derived from Brown 1991) calculation in "13 Amino acid composition of different phytoplankton.xlsx"
YcyanoC_N=2                             #(molC molN) C/N molar ratio of cyanophycin
YpgC_P=40                           #(molC molP) C/P molar ratio of PG: Phosphatidyl glycerol (assuming C16 fatty acids (since actually mostly C16 (Huflejt et al., 1990)
Nunit=1/Qc#*14*10**6/(12*10**3)          #((ug N / mgC)/(molN cell-1) unit conversion term (164-20)
Punit=1/Qc#*30.97*10**6/(12*10**3)       #((ug P / mgC)/(molP cell-1) unit conversion term (164-20)
Molar_mass_DNA_AT_average=307.47        #(g mol-1) Molar mass AT average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")
Molar_mass_DNA_CG_average=307.97        #(g mol-1) Molar mass CG average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")
Molar_mass_RNA_AT_average=316.47        #(g mol-1) Molar mass AT average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")
Molar_mass_RNA_CG_average=323.97        #(g mol-1) Molar mass CG average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")
CG_Ecoli=0.506          #(dimensionless) from [http://www.ncbi.nlm.nih.gov/genome/167 (accessed 06/18/2016)]
AT_Ecoli=1-CG_Ecoli     #(dimensionless) 

Molar_mass_DNA_Ecoli=Molar_mass_DNA_AT_average*CG_Ecoli+Molar_mass_DNA_CG_average*AT_Ecoli     #(g mol-1) Molar mass of DNA unit
Molar_mass_RNA_Ecoli=Molar_mass_RNA_AT_average*CG_Ecoli+Molar_mass_RNA_CG_average*AT_Ecoli     #(g mol-1) Molar mass of RNA unit

RNA_DNA_mass_ratio=17.844/6.5239  #(ug/ug) from values ad D=0 "07 Bremer and Dennis 1996 data plot.xlsx"
RNA_DNA_molar_ratio=RNA_DNA_mass_ratio/Molar_mass_RNA_Ecoli*Molar_mass_DNA_Ecoli    #(mol mol-1)
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
YphotoFe_N=0.001636364  #(molFe molN-1) Fe/N ratio in photosystem iron (0.001636364 in Inomura et al., 2022)

###### GRIDDED TEMPERATURE DATA ##############
Cq_temp = scipy.io.loadmat('C:\\Users\\mathp\\OneDrive\\Desktop\\Short_Course\\temp_gridded.mat')
# "C:\Users\mathp\OneDrive\Desktop\Short_Course\temp_gridded.mat"
temp = Cq_temp["ptemp2"]
t1 = temp[0:45,:,:]
t2 = temp[46:91,:,:]
woa_temp = np.concatenate((t1,t2),0)
Ea = 70000 #activation energy
R = 8.3
A =Ea/R
#Tt - ambient temperature that you are reading in
ab = np.zeros(shape=(24,90,180))
Tref = 20 + 273.15
Tt = woa_temp + 273.15 
Arr = exp(-A*((1/Tt)-(1/Tref)))

Cnbiosynth_temp=4.34728279914354E-10/Arr        #(molN cell-1 s) Constant for varible part of biosynthesis protein nitrogen (193-37)
Cnrna_variable_temp=6212.59249917364/Arr       #(s) Constant for Variable part of RNA (193-26)

Cnbiosynth = np.transpose(np.nanmean(Cnbiosynth_temp[:,:,0:depth_to_z],2))
Cnrna_variable = np.transpose(np.nanmean(Cnrna_variable_temp[:,:,0:depth_to_z],2))


afebest = pandas.read_csv("C:\\Users\\mathp\\OneDrive\\Desktop\\FE_opt_Spyd\\Optimization\\afebest.csv",header=None)
ynbest = pandas.read_csv("C:\\Users\\mathp\\OneDrive\\Desktop\\FE_opt_Spyd\\Optimization\\ynbest.csv",header=None)

# importing iron data
ds1 = nc.Dataset('C:\\Users\\mathp\\OneDrive\\Desktop\\FE_opt_Spyd\\Optimization\\DFe_Pasquier_Holzer_2018.nc')
lat = ma.getdata(ds1.variables["lat"][:])
lon = ma.getdata(ds1.variables["lon"][:])
depth = ds1.variables["depth"][:]
iron = ma.getdata(ds1.variables["DFe_TYP"][:]/1000)
iron_surf = np.nanmean(iron[0:depth_to_z,:,:],0)
Fe = iron_surf
light = ma.getdata(Ical())
surf_light = np.nanmean(light[0:depth_to_z,:,:],0)
# surf_light[surf_light<0]=nan

# plotting stuff
cmap1 = cm.get_cmap("RdYlBu_r",lut=20)
fig, axs = plt.subplots(nrows=2, ncols=5, figsize=(30, 6),subplot_kw={'projection':ccrs.PlateCarree()},layout='constrained')
# fig, axs = plt.subplots(nrows=2, ncols=4, figsize=(34, 8),subplot_kw={'projection':ccrs.PlateCarree()},layout='constrained')

axs = axs.ravel()

for j in range(10):
   aFe = afebest.loc[j,0]
   Ynphoto_chl = ynbest.loc[j,0]
   Pchl=Pmax*(1-exp(-OT*surf_light)) #(C mol s-1 Chl mol-1) Carbohydrate fixation rate per chlorophyll (167-1)(193-25)
   Vfe=aFe*Fe    #(mol P cell-1 s-1) iron uptake per cell

   A=((1+E)*Qc*Ynphoto_chl)/Pchl+Cnbiosynth
   B=Nconst_protein+(m*Ynphoto_chl)/Pchl
   L = ((1 + E)*Qc*Ypthylakoid_chl)/Pchl
   M = (m*Ypthylakoid_chl)/Pchl
   #========================
   # Fe limitation related
   #========================
   R=((1+E)*Qc*Ynphoto_chl*YphotoFe_N)/Pchl
   S=(m*Ynphoto_chl*YphotoFe_N)/Pchl
   
   aFe = R
   bFe = S
   cFe = -Vfe
   
   # DFe=DSolver(aFe,bFe,cFe)
   # D23=DFe.rQ
   DFe=solver_2D(aFe,bFe,cFe)
   D23=DFe.X
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
   
   msk = (Qc_essential >= Qc) | (np.isnan(Qc_essential))
   uhh = np.nonzero(msk) 
   xidx = uhh[0]
   yidx = uhh[1]
   for k in range(len(uhh[0])):
       x1 = xidx[k]
       y1 = yidx[k]
       D23[x1,y1] = 2*(Qc - Qc_const[x1,y1])/(Qc_D[x1,y1] + sqrt(Qc_D[x1,y1]*Qc_D[x1,y1] + 4*Qc_D2[x1,y1]*(Qc - Qc_const[x1,y1])))


   im=axs[j].pcolormesh(lon,lat, np.transpose(D23*86400),cmap=cmap1,vmin=0,vmax=3)
   lon_formatter = cticker.LongitudeFormatter()
   lat_formatter = cticker.LatitudeFormatter()
   if  j==1 or j==2 or j == 3 or j ==4:
          axs[j].set_yticks([])
          axs[j].set_xticks([])
   elif  j==6 or j==7 or j == 8 or j ==9:
          lon_formatter = cticker.LongitudeFormatter()
          axs[j].xaxis.set_major_formatter(lon_formatter)
          axs[j].set_xticks(np.arange(-120,121,60), crs=ccrs.PlateCarree())
   elif j == 5:
          lat_formatter = cticker.LatitudeFormatter()
          axs[j].yaxis.set_major_formatter(lat_formatter)
          axs[j].set_yticks(np.arange(-90,91,30), crs=ccrs.PlateCarree())
       
          lon_formatter = cticker.LongitudeFormatter()
          axs[j].xaxis.set_major_formatter(lon_formatter)
          axs[j].set_xticks(np.arange(-120,121,60), crs=ccrs.PlateCarree())
   elif j==0:
          lat_formatter = cticker.LatitudeFormatter()
          axs[j].yaxis.set_major_formatter(lat_formatter)
          axs[j].set_yticks(np.arange(-90,91,30), crs=ccrs.PlateCarree())
         
   # if  j==1 or j==2 or j == 3:
   #             axs[j].set_yticks([])
   #             axs[j].set_xticks([])
   # elif  j==5 or j==6 or j == 7:
   #             lon_formatter = cticker.LongitudeFormatter()
   #             axs[j].xaxis.set_major_formatter(lon_formatter)
   #             axs[j].set_xticks(np.arange(-120,121,60), crs=ccrs.PlateCarree())
   # elif j == 4:
   #             lat_formatter = cticker.LatitudeFormatter()
   #             axs[j].yaxis.set_major_formatter(lat_formatter)
   #             axs[j].set_yticks(np.arange(-90,91,30), crs=ccrs.PlateCarree())
             
   #             lon_formatter = cticker.LongitudeFormatter()
   #             axs[j].xaxis.set_major_formatter(lon_formatter)
   #             axs[j].set_xticks(np.arange(-120,121,60), crs=ccrs.PlateCarree())
   # elif j==0:
   #             lat_formatter = cticker.LatitudeFormatter()
   #             axs[j].yaxis.set_major_formatter(lat_formatter)
   #             axs[j].set_yticks(np.arange(-90,91,30), crs=ccrs.PlateCarree())

   # axs[j].yaxis.set_major_formatter(lat_formatter)
   axs[j].add_feature(cartopy.feature.COASTLINE, edgecolor='black')
   axs[j].set_extent([-180, 180, -90, 85],ccrs.PlateCarree())


# fig.colorbar(im,ax=axs.ravel().tolist(),location='bottom',fraction=0.046, pad=0)
# plt.colorbar(im, ax=axs.ravel().tolist(),fraction=0.046, pad=0.04)
# cb_ax = fig.add_axes([0.15, 0.94, 0.8, 0.03])
# cbar = fig.colorbar(im, ax=axs.ravel().tolist(),cmap=cmap1,orientation='horizontal')
# cbar.set_label('$\mathit{\mu}$ (d$^{-1}$)')
# cbaxes = fig.add_axes([0.1, 0.1, 0.02, 0.8])  
 
# position for the colorbar
cb = plt.colorbar(im, ax=axs.ravel().tolist(),cmap=cmap1,orientation='vertical',location="right",shrink=0.9,pad=0.017)
cb.set_label('$\mathit{\mu}$ (d$^{-1}$)')
fss = 32
fig.text(0.176,0.565,'a',fontsize=fss)
fig.text(0.355,0.565,'b',fontsize=fss)
fig.text(0.54,0.565,'c',fontsize=fss)
fig.text(0.73,0.565,'d',fontsize=fss)
fig.text(0.913,0.565,'e',fontsize=fss)
fig.text(0.176,0.1,'f',fontsize=fss)
fig.text(0.355,0.1,'g',fontsize=fss)
fig.text(0.54,0.1,'h',fontsize=fss)
fig.text(0.73,0.1,'i',fontsize=fss)
fig.text(0.913,0.112,'j',fontsize=fss)


# fig.text(0.21,0.545,'A.',fontsize=fss)
# fig.text(0.448,0.545,'B.',fontsize=fss)
# fig.text(0.684,0.545,'C.',fontsize=fss)
# fig.text(0.915,0.545,'D.',fontsize=fss)
# fig.text(0.21,0.082,'E.',fontsize=fss)
# fig.text(0.448,0.082,'F.',fontsize=fss)
# fig.text(0.684,0.082,'G.',fontsize=fss)
# fig.text(0.915,0.082,'H.',fontsize=fss)
# fig.text(0.20,0.0001,'J.',fontsize=fss,color="white")
sf_opt("subplots_mu_627")
# plt.savefig("C:\\Users\\mathp\\OneDrive\\Desktop\\FE_opt_Spyd\\Optimization\\Figures\\Opt\\new_species\\mu_save2")
# sf_opt("all_maps_mu_withk_newhuxetc_NEWIRON_upper200")




# fig, ax = plt.subplots(figsize = [12,8],subplot_kw={'projection':ccrs.PlateCarree()})
# x,y = np.meshgrid(lon, lat)
# im = ax.pcolormesh(lon, lat, np.transpose(Pmax-Pchl),cmap=cmap1,vmin=0)
# cbar = plt.colorbar(im, ax = ax, orientation = 'vertical', fraction = 0.02, pad = 0.02)
# cbar.set_label('Pchl')
# ax.set_xticks(np.arange(-180,181,60), crs=ccrs.PlateCarree())
# lon_formatter = cticker.LongitudeFormatter()
# ax.xaxis.set_major_formatter(lon_formatter)
# ax.set_yticks(np.arange(-90,91,30), crs=ccrs.PlateCarree())
# lat_formatter = cticker.LatitudeFormatter()
# ax.yaxis.set_major_formatter(lat_formatter)
# ax.coastlines()