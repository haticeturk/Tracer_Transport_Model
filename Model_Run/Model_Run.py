# %%
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as dates
import datetime as dt
import numpy.testing as npt
import numpy as np
import matplotlib.pyplot as plt
import logging
import importlib
import scipy as sp
import warnings
import seaborn as sns   
#import psutil
#import statsmodels 
from pandas import DataFrame
from numpy import vstack, add, float64, multiply, exp, zeros, nan_to_num, vstack, where, power
from numpy import add, float64, multiply, exp, zeros, nan_to_num, vstack, where, power
from scipy.stats import rv_continuous
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import scipy.stats
from scipy.stats import gamma, beta
from numpy import add, float64, multiply, exp, zeros, nan_to_num, vstack, where, power
from numpy import sin, arange, isclose
from scipy.optimize import differential_evolution, Bounds, shgo, dual_annealing, minimize
from scipy import signal
from importlib import reload
warnings.filterwarnings("ignore")
#plot the results
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import os
from mpl_toolkits import mplot3d
from matplotlib import cm, colors
#from statsmodels.graphics.tsaplots import plot_acf

#Set your  directory
os.chdir(". /Tracer_Transport_Model/Model_Run")
saveDir = './racer_Transport_Model/Output/Plots'
import tracer_models as tracer_models
from tracer_models import ObjectiveFunction
from tracer_models import Tracer_Mod_Wetness

# %%
## Preperation of the data 
#Read Hyrological Data 
meas_df =  pd.read_csv("Data/Data_test.csv", parse_dates=['Date'],index_col='Date',).squeeze("columns")
#Eliminate Duplicate Date values
meas_df = meas_df[~meas_df.index.duplicated(keep='first')]
#interpolote NAN values meas_df PET
meas_df.PET = meas_df.PET.interpolate(method='linear', limit_direction='both')
# Meas_df  select data from date 
# Example: meas_df = meas_df.loc['2009-10-01':'2013-09-30'] 
meas_df = meas_df.loc['2009-10-01':'2013-09-30']

# Define the warm up period and calibration period
wp =365    ## warm up period corresoond to 1 year
nt = meas_df.Qmmd.size ## calibration period end index as length of the data 


## lets prepare the data to numpy array input format to run the model
evap =meas_df.PET[0:nt].to_numpy()
#apend zero to evap < 0
evap = np.where(evap < 0, 0, evap)
prec = meas_df.Pmmd[0:nt].to_numpy()
temp = meas_df.T_mean[0:nt].to_numpy()
flow = meas_df.Qmmd[0:nt].to_numpy()

## Measured tracer data "Deuterium" data P is precipitation and Q is flow
ConP_H= meas_df.P_H[0:nt ].to_numpy()
ConQ_H= meas_df.Q_H[0:nt ].to_numpy()

#Lets  calculate 3 month average to get Runoff Coefficient

df = pd.DataFrame({"flow":meas_df.Qmmd[wp:nt] , "prec": meas_df.Pmmd[wp:nt]})
df_n  = df.groupby(np.arange(len(df))//90).mean()
df_n['QperP'] = df_n.flow/df_n.prec


##  Prepare the inital storage Tracer composition for  unsaturated zone Su and Passive Groundwater Ssp
valid_indices_H = ~np.isnan(ConQ_H)
filtered_Qo_H = flow[valid_indices_H]
filtered_ConcQ_H = ConQ_H[valid_indices_H]
percentile_5_H = np.percentile(filtered_Qo_H, 5)
filtered_ConcQ_H = filtered_ConcQ_H[filtered_Qo_H <= percentile_5_H]
ConcSsp_avg_H = np.nanmean(filtered_ConcQ_H)
ConcSu_avg_H = np.nansum(ConP_H* prec) / np.nansum(prec)

# %%
############################## Model calibration ############################## 
#Set the case to "Cal" to inform the model that is calibration case
########If the Case is "Cal" the model will retuns ########
######### 1) Total Flux matrix  =C_Q*Q or C_rm*RM................ such that
# Flux_df = pd.DataFrame({'Flux_RM': Flux_RM, 'Flux_Pe': Flux_Pe, 'Flux_Ei': Flux_Ei, 'Flux_rf': Flux_rf, 'Flux_ea': Flux_ea, 'Flux_rs': Flux_rs, 'Flux_rfs': Flux_rfs, 'Flux_rfn': Flux_rfn, 'Flux_Qo': Flux_Qo, 'Flux_Qf': Flux_Qf, 
# Flux_Qof': Flux_Qof, 'Flux_Qstot': Flux_Qstot, 'Flux_Ql': Flux_Ql})
########## 2) Age matix of hydrological fluxes  = Q_i*TTD or RM_i*TTD.......... such that each matrix is a 2D matrix with rows as age and columns as time step ########
# Age_MR = [RM_Age_MR, Pe_Age_MR, rf_Age_MR, rfs_Age_MR, rfn_Age_MR, Qf_Age_MR, Qof_Age_MR, Qstot_Age_MR, Ql_Age_MR, Qo_Age_MR, rs_Age_MR, Ei_Age_MR, ea_Age_MR]
#  
########### 3) df_WB ## will return the hydrological model water balance components and time variable alphat ########
# df_WB = ["Qtot", "Qs", "Ql", "Qf", "Qof", "Qo", "Ss", "Si", "Sr", "Ssa", "Ssp", "Sf", "Rs","Rsg", "Rmelt", "Rff", "Rfn", "Rfs", "Ei", "Ea", "Etot", "Pe",
#Ps", "M", "Cr","Cn", "alphat"]
############################## ----Define the Model Parameter boudries----#################################
parameters = pd.DataFrame(columns=["initial", "pmin", "pmax"])
parameters.loc["srmax"] = (300.0, 100, 500)          #0  Root-zone storage capacity
parameters.loc["lp"] =    (0.26, 1e-5, 1)            #1  Transpiration water stress factor
parameters.loc["Rs\max"] =  (1.05, 0.0, 1.2)         #2  Maximum percolation rate
parameters.loc["beta"] = (0.17, 1e-5, 5.0)           #3  Parameter determining the nonlinearity of outflow / recharge. Shape factor
parameters.loc["Kp"] =  (1e-8, 1e-10, 1e-5)          #4  Loss coefficient
parameters.loc["Imax"] = (4.97, 1.2, 5.0)            #5  Interception capacity
parameters.loc["Ka"] = (0.18, 0.01, 0.2)             #6  Storage coefficient of the slow-responding reservoir
parameters.loc["tt"] = (-2.0, -4.0, 5.0)             #7  Threshold temperature for snow melt
parameters.loc["fmelt"] =  (2.81, 1.0, 5.0)          #8  Melt factor
parameters.loc["sfmax"] = (1.013, 1.0, 20.0)         #9  Fast response storage capacity
parameters.loc["kf"] = (0.9, 0.01, 2.0)              #10 fast responce coefficient
parameters.loc["cp"] = (0.5, 0.1, 1.0)               #11 Division parameter for fast groundwater recharge (1-cp) =Rfs
parameters.loc["cn"] = (0.6, 0.1, 1.0)               #12 Division parameter for fraction of overland flow
parameters.loc["Bf"] = (1e-8, 1e-10,1e-5)            #13 Saturation-excess overland flow coefficient
parameters.loc["SSp"] = (3000.0, 500.0, 10000.0)     #14 Passive storage capacity
parameters.loc["Ptres"] = (6.0, 2.0, 20.0)           #15 Threshold precipitation for overland flow
parameters.loc["SU_Alpha"] = (0.1, 0.01, 1.0)        #16 SAS alpha shape parameter for root zone 
parameters.loc["SG_Alpha"] = (0.99, 0.98, 1.0)       #17 SAS alpha shape parameter for groundwater

##############################----Get inital parameter to initialize the model----#################################
par =parameters.initial.values

#######---- Define fittenss function to minimize----#######
def fitness(par):
    par =par
    Flux_TB, Agn_MR , df_WB = Tracer_Mod_Wetness.Tracer_Mod().get_tracer_balance(prec, evap, temp, ConP_H, ConcSsp_avg_H, ConcSu_avg_H,  p = par,  Case ='Cal')
    Flux_Qtot  = Flux_TB.Flux_Qstot+ Flux_TB.Flux_Qo + Flux_TB.Flux_Qof + Flux_TB.Flux_Qf
    # get Qtot _Age  = Qf_Age_MR, Qof_Age_MR, Qstot_Age_MR,
    Qtot_Age  = Agn_MR[5] + Agn_MR[6] + Agn_MR[7]+ Agn_MR[9]
    #Sum Qtot_Age columns
    Qtot = Qtot_Age.sum(axis=0)
    Flux_Qtot = Flux_Qtot[:-1] # remove the last value as it is extra
    Conc_sim  = Flux_Qtot/Qtot[:-1]
    C_obs = ConQ_H[wp:nt]
    C_sim  = Conc_sim[wp:nt].values
    sim_Q = df_WB.Qtot[wp:nt].values
    obs_Q = flow[wp:nt]
    QP= df_n.QperP
    Q_Avarage_sim  = pd.DataFrame({"flow":df_WB.Qtot[wp:nt] , "prec":prec[wp:nt]})
    df_OP = Q_Avarage_sim.groupby(np.arange(len(Q_Avarage_sim))//90).mean()
    OP_sim = df_OP.flow/df_OP.prec
    Ei =  ObjectiveFunction.nashsutcliffe(obs_Q, sim_Q)
    TB_nse = ObjectiveFunction.nashsutcliffe(C_obs,C_sim)
    Eii = ObjectiveFunction.lognashsutcliffe(obs_Q, sim_Q)
    Eiii =  ObjectiveFunction.nse_FDC(obs_Q, sim_Q)
    Ebase = ObjectiveFunction.nashsutcliffe(QP,OP_sim)
    #Eve = ObjectiveFunction.volume_error(QP, OP_sim)
    
    DE = np.sqrt(0.5*(((0.25*(1-Ei)**2 + 0.25*(1-Ebase)**2 + 0.25*(1-Eii)**2 +0.25*(1-Eiii)**2)) + (1-TB_nse)**2))
    return DE

#######---- Define the parameter bounds----#######
bounds = Bounds(lb= parameters.pmin.values, ub=parameters.pmax.values, keep_feasible=False)

#######----Initialize a list to hold solutions that meet the criterion ----#######
satisfactory_solutions = []

# Maximum number of solutions to keep
max_solutions = 1000  # Adjust as needed
def callback_efficient(xk, convergence):
    global satisfactory_solutions
    # Check if the solution meets the criterion
    current_fitness = fitness(xk)
    if current_fitness < 1:  # Criterion: fitness value should be less than 1
        satisfactory_solutions.append((xk, current_fitness))
        # Sort the list based on fitness (for minimizing)
        satisfactory_solutions.sort(key=lambda sol: sol[1])
        # Ensure the list does not exceed the maximum length
        satisfactory_solutions = satisfactory_solutions[:max_solutions]
# %%
# Run the differential evolution optimizer set the workers as nedeed
# Change the updating parameter to 'differ" to change the updating strategy
#fitness is your fittnes function valu of "DE" 
result = differential_evolution(fitness, bounds,  callback=callback_efficient, updating='immediate', workers=8, maxiter=500)

print("optimization is done!!!!")

# Save Result .x as data frame 
df_result = pd.DataFrame({"par":result.x}) 
df_result.to_csv(".Set_Directory_to save/Cal_Results/Cal_Parameters.csv")
## Save  satisfactory_solutions
satSol = pd.DataFrame(satisfactory_solutions, columns=["par", "fitness"])
import pickle
satSol.to_pickle("..Set_Directory_to save/Cal_Results/Satisfactory_solutions.pkl")

# %%
############Tracer Model with calibrated parameters################

Parameter = pd.read_csv("Cal_Results/Cal_Parameters.csv")
parOpt = Parameter.par.values
################Run the model with calibrated parameters ################
# Set the case as "notCal" to retun the Storage age matrix and Flux matrix as well as Tracer balance components
# Storage_Age_MR =  [SS_Age_MR, Si_Age_MR, Sr_Age_MR, Sf_Age_MR, SSa_Age_MR]
# Storage_Flux_MR = [SS_Flux_MR, Si_Flux_MR, Sr_Flux_MR, Sf_Flux_MR, SSa_Flux_MR]

Agn_MR, Flux_dfn, Storage_Age_MR, Storage_Flux_MR, TrB, WB_test =Tracer_Mod_Wetness.Tracer_Mod().get_tracer_balance(prec, evap, temp, ConP_H, ConcSsp_avg_H, ConcSu_avg_H, p = parOpt,  Case ='notCal')
# analyze the results 
Sim_df =  WB_test # hydrological fluexes
Sim_df.index  = meas_df.index[0:nt]
Sim_df['date'] = Sim_df.index.map(lambda x: x.strftime('%Y-%m-%d'))
#Check goodnes of fit 
Evaluation =pd.DataFrame(ObjectiveFunction.calculate_all_functions(flow[wp:nt], Sim_df.Qtot[wp:nt]))
# Get Qtot Age matrix as sum of  Qf_Age_MR, Qof_Age_MR, Qstot_Age_MR,
Qtot_Age  = Agn_MR[5] + Agn_MR[6]  + Agn_MR[7] + Agn_MR[9]
#Get Etot Age matrix as sum of Ei_Age_MR, ea_Age_MR
Etot_Age = Agn_MR[12]  + Agn_MR[11]
#Sum Qtot_Age by column to get total Q for each time step
Qtot = Qtot_Age.sum(axis=0)
#Sum Etot_Age by column to get total E for each time step
Etot = Etot_Age.sum(axis=0)
WB_Qtot = WB_test.Qtot[wp:nt].values
#Get total fluxes for streamflow and evaporation
Flux_Qtot  = Flux_dfn.Flux_Qstot+ Flux_dfn.Flux_Qo + Flux_dfn.Flux_Qof + Flux_dfn.Flux_Qf
Flux_Etot = Flux_dfn.Flux_Ei + Flux_dfn.Flux_ea
Conc_sim  = Flux_Qtot/Qtot
Conc_Esim = Flux_Etot/Etot
C_obs =ConQ_H[wp:nt]
C_sim  = Conc_sim[wp:nt].values
#Get the evaluation metrics for Tracer simulation
TB_nse = ObjectiveFunction.calculate_all_functions(C_obs,C_sim)

# %% 
#################Lets get backward TTD results  ############################
## Example of computing backward Transit Time Distribution (TTD) for streamflow and Evaporation ##
#  This can be extended to other fluxes such as recharge, runoff, etc

# Define traking length 
#### The TTD  can be computed by be Age_matrix diviwed by the sum of the age matrix for each column. Keep in mind that the row are 
# Age ranked from 0 to n. For example the last row is the oldest age and the first row is the youngest age. 
#If you want to sum age fraction T<90 than Tr_lwngth is = 90 .
Tr_length = nt # Traking lenght 
TTD_Qtot=np.zeros((Tr_length+1,Tr_length+1),  dtype=float64) 
for i in range(0,Tr_length):
    TTD_Qtot[:,i] = Qtot_Age[:,i]/np.sum(Qtot_Age[:,i])
TTD_Qtot[:,Tr_length] = Qtot_Age[:,Tr_length-1]/np.sum(Qtot_Age[:,Tr_length-1])
TTD_Evap=np.zeros((Tr_length+1,Tr_length+1),  dtype=float64)
for i in range(0,Tr_length):
    TTD_Evap[:,i] = Etot_Age[:,i]/np.sum(Etot_Age[:,i])
TTD_Evap[:,Tr_length] = Etot_Age[:,Tr_length-1]/np.sum(Etot_Age[:,Tr_length-1])
#Remove TTD_Qtot last column as its defined extra one column more in the model
TTD_Qtot = TTD_Qtot[:, :-1]
Qtot_Age = Qtot_Age[:, :-1]
TTD_Evap = TTD_Evap[:, :-1]
Etot_Age = Etot_Age[:, :-1]

# %%
##Plot culmulative distribution function bacward TTD
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(25, 10))
for i in range(1000, Tr_length):
     ax1.plot(np.cumsum(TTD_Qtot[0:1000,i]), linewidth=0.5, color="#7C9D96")
ax1.set_xlabel('Age[d]')
for i in range(1000, Tr_length):
     ax2.plot(np.cumsum(TTD_Evap[0:1000,i]), linewidth=0.5, color="#7C9D96")
ax1.set_xlabel('Age[d]')
#Set x axes from zero 
ax1.set_xlim(1, 1000)
ax2.set_xlim(1, 1000)
#Set X label Age
ax2.set_xlabel('Age', fontsize=20)
#Set y lim ax2
ax2.set_ylim(0, 1.5)
#scale xaxes log
ax1.set_ylabel('$eCDF\ [-]$')
ax2.set_ylabel('$eCDF\ [-]$')
#Set x axes 0 to 100
ax2.set_xscale('log')
ax1.set_xscale('log')
plot_title = 'Cumulative Distribution Function of backward TTD'
plt.show()
