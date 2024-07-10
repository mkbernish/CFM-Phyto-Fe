
'''
Created on May 18, 2014
This one reads dpi
@author: Keisuke
'''
from pylab import * 

def sf_opt(figName):
    First_part="C:\\Users\\mathp\\OneDrive\\Desktop\\FE_opt_Spyd\\Optimization\\Figures\\Opt\\new_species\\"
    Figure_name=str(figName)
    Last_part=".png"
    savefig(First_part+Figure_name+Last_part,dpi=450)
