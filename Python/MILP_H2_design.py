
"""

Citation notice:

If you use the MILP optimization framework, please cite:
F. Superchi and A. Bianchini, Renewable Energy, vol. 256, p. 124470, 2026.
DOI: 10.1016/j.renene.2025.124470
    
If you use the market analysis results implemented in this code, please cite:
F. Superchi et al., Renewable Energy, vol. 245, p. 122813, 2025.
DOI: 10.1016/j.renene.2025.122813

"""

import numpy as np 
import pandas as pd
import gurobipy as gp
from gurobipy import GRB

MILP_design = True

kWh_factor = 1    # hourly simulation


#%%
'Remuneration scheme'
P_wind_band         = 137.5   / 1000  #€/kWh
P_pv_band           = 112.574 / 1000  #€/kWh
P_wind_over         = 110     / 1000  #€/kWh
P_pv_over           = 60.148  / 1000  #€/kWh
P_bess              = 165     / 1000  #€/kWh

'MILP penalties'
#penalties to overshooting and undershooting energy
C_fix = 0.165/kWh_factor  #€/kW   0.165 €/kWh and a time span of 1 minute: time_factor = 60
C_pen_p = C_fix * 1.5 # * 10000
C_pen_n = C_fix * 1.5 # * 10000

#penalties applied to components utilization 
C_batt = 0.1     #cost of utilizing the BESS 
C_el   = 0.5     #cost of utilizing the electrolyzer
C_fc   = 1     #cost of utilizing the fuel cell
C_excess = 5

#%%
'grid power'
# C_grid = 165 #€/MWh puchase price from the grid
C_grid = 336 / 1000 #https://www.europarl.europa.eu/doceo/document/E-8-2017-006974_EN.html?
# C_grid = 600 /1000 #10.1016/j.egypro.2018.12.065

P_grid_max = 400 #maximum power absorption from the grid

#%%

"""
USER INPUT REQUIRED: hourly time series (1 year)

This framework expects three *hourly* time series over the same period:
- load: electrical demand (e.g., kW or MW)
- PV: PV electrical power (uncurtailed or delivered, same unit as load)
- WT: wind turbine electrical power (uncurtailed or delivered, same unit as load)

How to provide the data (choose ONE option):
1) Provide three CSV files (one per series), each containing:
   - a timestamp column named "time"
   - a value column (see VALUE_COL_* below)

2) Provide three pandas DataFrames already loaded in memory, with the same structure.

Important assumptions:
- Time resolution must be exactly 1 hour (hourly).
- The three series must share the same timestamps.
- Timestamps should be timezone-consistent (preferably timezone-naive or all in the same timezone).
- Expected length is usually 8760 hours (non-leap year) or 8784 (leap year).
"""

'RES and load dataset'

date = pd.date_range('xxx-xx-xx 00:00:00', 'xxxx-xx-xx 00:00:00', freq='h')

df_hh = pd.DataFrame()
df_hh['time'] = date
df_hh['load'] = pd.read_csv('Load.csv', sep = ';')
df_hh['PV'] = pd.read_csv('PV.csv', sep = ';')
df_hh['WT'] = pd.read_csv('WT.csv', sep = ';')



#%%
date = pd.date_range('2022-11-01 00:00:00', '2023-11-01 00:00:00', freq='h')

df=pd.DataFrame() #creating the DataFrame
df['time_stamp']=date
df['month']=df.time_stamp.dt.month
df['day']=df.time_stamp.dt.day
df['hour']=df.time_stamp.dt.hour
df['minute']=df.time_stamp.dt.minute
df['load'] = df_hh['load']
df['PV_power'] = df_hh['PV']
df['wind_power'] = df_hh['WT']
df = df.fillna(value = 0)



#%%

df['SOC']=list([50.000]*len(df.time_stamp))
df['BESS_power']=list([0]*len(df.time_stamp))
df['power_out_var']=list([0]*len(df.time_stamp))
df['power_out_tot']=list([0]*len(df.time_stamp))
df['earnings']=list([0]*len(df.time_stamp))
df['curt']=list([0]*len(df.time_stamp))
df['gap']=list([0.00]*len(df.time_stamp))

#%%
'cost of components'
#market analysis form 10.1016/j.renene.2025.122813

prices = pd.read_excel('prices_excel.xlsx', sheet_name = ['Li-BESS', 'ALK EL', 'PEM FC', 'H2 Tank', 'PV', 'Onshore WT'],
                       usecols = 'V:Y', skiprows= [0,1], nrows = 3 )

year = 2020   #time horizon for prices

if year == 2020:
    k = 0
elif year == 2030:
    k = 1
elif year == 2050:
    k = 2

EL_cost      = prices['ALK EL']['avg'][k]           # €/kW
FC_cost      = prices['PEM FC']['avg'][k]           # €/kW     
HP_tank_cost = prices['H2 Tank']['avg'][k]           # €/kg_h2
LP_tank_cost = prices['H2 Tank']['avg'][k]           # €/kg_h2
bess_cost    = prices['Li-BESS']['avg'][k]           # €/kWh
WT_cost      = prices['Onshore WT']['avg'][k]          # €/kW
PV_cost      = prices['PV']['avg'][k]           # €/kWp

comp_cost = 60000    #€/unit

components = {}

components['EL'] = {'total installation costs': EL_cost,        # €/
                    'OeM': 0.0275*EL_cost,                      # €/kW/y
                    'lifetime': 10, 'relpacement': 0.4}

components['FC'] = {'total installation costs': FC_cost,           # €/kW - ref. file 'Costi.xls' 460 €/kg
                    'OeM': 0.0275*FC_cost,                          # €/kW/h
                    'lifetime': 10, 'relpacement': 0.4}

components['BESS'] = {'total installation costs': bess_cost,       # €/MWh
                      'OeM': 0.025*bess_cost,                      # €/MWh/y
                      'lifetime': 10, 'relpacement': 0.8}

components['HP_tank'] = {'total installation costs': HP_tank_cost,      # €/kg
                        'OeM': 0.01*HP_tank_cost,                       # €/kg/y
                        'lifetime': 25, 'relpacement': 0}

components['LP_tank'] = {'total installation costs': LP_tank_cost,      # €/kg
                        'OeM': 0.01*LP_tank_cost,                       # €/kg/y
                        'lifetime': 25, 'relpacement': 0}

components['WT'] = {'total installation costs': WT_cost,           # €/kW - ref. file 'Costi.xls' 460 €/kg
                     'OeM': 0.025*WT_cost,                                    # €/kW/h
                     'lifetime': 25, 'relpacement': 0}

components['PV'] = {'total installation costs': PV_cost,           # €/kW - ref. file 'Costi.xls' 460 €/kg
                    'OeM': 0.025*PV_cost,                          # €/kW/h
                    'lifetime': 25, 'relpacement': 0}

components['compressor'] = {'total installation costs': comp_cost,           # €/kW - ref. file 'Costi.xls' 460 €/kg
                            'OeM': 0.025*comp_cost,                                    # €/kW/h
                            'lifetime': 25, 'relpacement': 0, 'size' : 1}

#%%
'''Components parameters'''
'BESS parameters'
eta_c = 0.95       
eta_d= 0.95

C_rate_charge    = 1   # limit to C-rate charge of 1C
C_rate_discharge = 3   # limit to C-rate discharge of 3C

'EL parameters'
phi_el = 18  /1000/kWh_factor # kg/MWh to kg/kWh for this timestep - kg of H2 produced per MWh electrical energy input

'FC parameters'
phi_fc = 59  /1000/kWh_factor #kg/MWh to kg/kWh for this timestep - kg of H2 consumed per MWh electrical energy output

'compressor power'
E_comp = 1.18  #kWh/kg   energy requets to compress 1 kg of H2
P_comp = 70  #kW   power absorption to 

eta_comp = 0.965  # compression efficiency

'MILP parameters'
Time_Limit = 5*24*3600
MIPgap = 0.01

'''MILP optimization'''
end_of_sim= 24 * kWh_factor #length of time horizon: 1 day (24 h * 60)
M = np.arange(1,13)  # months
Ms = np.arange(4,10).tolist()    # during summer months overshootings are allowed and remunerated

'fixed design'
if MILP_design == False:
    BESS_capacity = 7 * 1000  # MWh to kWh
    P_el_nom = 250                # kW nominal power of electrolyzer
    P_fc_nom =  890               #kW nominal power of fuel cell
    tank_capacity = 6000 # kg capacity of the H2 tank
    PV_scale = 5 # x160kWp   - 1810/160 for 100% SSd

#%%

#####
'whole year sim'
i = 0
hours = len(df) 
#####
'(for cycle)'
# i = 464             #winter demo
# i = 6576            #summer demo 1
# i = 6976            #summer demo 2
# i = 6076            #summer demo 3

# hours = 24*30

# i = 5000
# hours = 24*60

#####

#%%

WT_pot = np.array(df.loc[i:i+(hours*kWh_factor - 1),'wind_power'])
PV_pot0 = np.array(df.loc[i:i+(hours*kWh_factor - 1),'PV_power']) #* PV_scale
load = np.array(df.loc[i:i+(hours*kWh_factor - 1),'load'])

n = len(WT_pot)                               #len of time horizon
length = np.arange(n)

'Gurobi Model'

m = gp.Model("mip1")                  #model selection
m.setParam('TimeLimit', Time_Limit)   #Limits the total time expended (in seconds).
m.setParam('MIPGap', MIPgap)          #The MIP solver will terminate when the gap between the lower and upper objective bound is less than MIPGap times the absolute value of the incumbent objective value.

m.setParam('MIPFocus', 1) 
m.setParam('Presolve', 2) # 0 to turn it off
m.setParam('Heuristics', 0.3) # 0 to turn it off

# m.setParam("Threads", 18)  # limits threads to preserve cpu

'control variables'
#variables to be optimized are the quantities trend during the following day: n-sized vectors
P_fix      = m.addVars(length,    vtype=GRB.INTEGER,    name= "P_fix" )   #integer because it must lay on one of the declared power levels
y          = m.addVars(length, 3, vtype=GRB.BINARY,     name='power')     #3d auxiliary binary variable to constraint power levels

P_var_p    = m.addVars(length,    vtype=GRB.CONTINUOUS, name="P_var_p")   #variable power can be continuous 
P_var_n    = m.addVars(length,    vtype=GRB.CONTINUOUS, name="P_var_n")

WT_band    = m.addVars(length,    vtype=GRB.CONTINUOUS, name="P_wind_band")
PV_band    = m.addVars(length,    vtype=GRB.CONTINUOUS, name="P_pv_band")
WT_over    = m.addVars(length,    vtype=GRB.CONTINUOUS, name="P_wind_over")
PV_over    = m.addVars(length,    vtype=GRB.CONTINUOUS, name="P_pv_over")
WT_excess  = m.addVars(length,    vtype=GRB.CONTINUOUS, name="P_wind_excess")
PV_excess  = m.addVars(length,    vtype=GRB.CONTINUOUS, name="P_pv_excess")
WT_load    = m.addVars(length,    vtype=GRB.CONTINUOUS, name="P_wind_load")   # wind energy that covers the load
PV_load    = m.addVars(length,    vtype=GRB.CONTINUOUS, name="P_pv_load")     # solar energy that covers the load

WT_bess    = m.addVars(length,    vtype=GRB.CONTINUOUS, name="P_wind_load")   # wind energy that charges the bess
PV_bess    = m.addVars(length,    vtype=GRB.CONTINUOUS, name="P_pv_load")     # solar energy that charges the bess
WT_el    = m.addVars(length,    vtype=GRB.CONTINUOUS, name="P_wind_load")   # wind energy that feeds the electrolyzer
PV_el    = m.addVars(length,    vtype=GRB.CONTINUOUS, name="P_pv_load")     # solar energy that feeds the electrolyzer

P_BESS_tot = m.addVars(length,    vtype=GRB.CONTINUOUS, name="P_bess_tot")
P_BESS_p   = m.addVars(length,    vtype=GRB.CONTINUOUS, name="P_bess_p")    #power coming from the BESS (discharging)

P_BESS_p_load   = m.addVars(length,    vtype=GRB.CONTINUOUS, name="P_bess_load")    #power coming from the BESS (discharging)
P_BESS_p_grid   = m.addVars(length,    vtype=GRB.CONTINUOUS, name="P_bess_grid")    #power coming from the BESS (discharging)

P_BESS_n   = m.addVars(length,    vtype=GRB.CONTINUOUS, name="P_bess_n")    #power going to the BESS (charging)
h          = m.addVars(length, 2, vtype=GRB.BINARY,     name='battery')   #auxiliary variable to prevent simulatenous charge-discharge
P_BESS_abs = m.addVars(length,    vtype=GRB.CONTINUOUS, name="P_bess_abs")

H2_tank  = m.addVars(length,    vtype=GRB.CONTINUOUS, name="H2 soc",  lb = 0 ) #limitations to the SOC of H2 tank
En_BESS  = m.addVars(length,    vtype=GRB.CONTINUOUS, name="BESS soc")            #limitations to the SOC of BESS

P_el      = m.addVars(length,    vtype=GRB.CONTINUOUS, name="P_el", lb = 0)
P_fc      = m.addVars(length,    vtype=GRB.CONTINUOUS, name="P_fc", lb = 0)
k         = m.addVars(length, 2,   vtype=GRB.BINARY,     name='hydrogen')   #auxiliary variable to prevent simulatenous produiction and consumption of H2

P_excess        = m.addVars(length,    vtype=GRB.CONTINUOUS, name="P_excess", lb = 0)
P_exp        = m.addVars(length,    vtype=GRB.CONTINUOUS, name="P_exp", lb = 0)

P_grid          = m.addVars(length,    vtype=GRB.CONTINUOUS, name="P_grid", lb = 0)        # power taken from the grid to cover the load
P_load_covered  = m.addVars(length,    vtype=GRB.CONTINUOUS, name="P_load_covered", lb = 0)        # power taken from the grid to cover the load

PV_pot  = m.addVars(length,    vtype=GRB.CONTINUOUS, name="P_pv")        # power taken from the grid to cover the load
# PV_pot = PV_pot0


'design variables'
if MILP_design == True:
    BESS_capacity = m.addVar(vtype=GRB.INTEGER, name="BESS_capacity",  lb = 0 , ub = 20000)  # kWh     #2880
    P_el_nom      = m.addVar(vtype=GRB.INTEGER, name="P_el" , lb = 0 , ub = 1000)             # kW
    P_fc_nom      = m.addVar(vtype=GRB.INTEGER, name="P_fc", lb = 0 , ub = 1000)              # kW
    tank_capacity = m.addVar(vtype=GRB.INTEGER, name="Tank_capacity", lb = 0 , ub = 10000)    # kg
    PV_scale      = m.addVar(vtype=GRB.INTEGER, name="PV_scale" , lb = 1 , ub = 1)            # x160 kW
    # PV_scale = 1

'cost variables'
years    = 21   #lifetime of the project
r        = 0.05
d_vector = [1/(1+r)**x for x in np.arange(0,years)]     #depreciation vector

R_t    = m.addVars(years,    vtype=GRB.CONTINUOUS, name="R_t", lb=-GRB.INFINITY)    # revenues
OEM_t  = m.addVars(years,    vtype=GRB.CONTINUOUS, name="OEM_t")  # operation and maintenance
C_t    = m.addVars(years,    vtype=GRB.CONTINUOUS, name="C_t")    # substitution costs

I0     = m.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name="I0")    #investment cost
NPV    = m.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name="NPV")  #net present value
C_0    = m.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name="C0")      #subtstitution cost
OEM_0  = m.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name="O&M0")  #yearly O&M

R_0 = m.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name="R_0")    

'''control constraints'''
m.addConstrs( PV_pot[z] == PV_pot0[z] * PV_scale for z in range(n))

'Power in allowed bands'
m.addConstrs( (P_fix[z] == 0 * y[z,0] + 200 * y[z,1] + 400 * y[z,2] for z in range(n)) )   #fixed power must be equal to one of the defined power levels in each timestep
m.addConstrs((y[z,0] + y[z,1] + y[z,2] == 1 for z in range(n)))     #fixed power must lay only in one of the fixed power levels simultaneously

# m.addConstrs( P_fix[z] == P_fix[z+i] for i in range(kWh_factor) for z in range(0,n-kWh_factor,kWh_factor))  #mantain the same fixed power output for one-hour

'H2 constraints'
m.addConstrs(P_el[z] <= P_el_nom     * k[z,0] for z in range(n))       #max el power limitation
m.addConstrs(P_el[z] >= P_el_nom*0.2 * k[z,0] for z in range(n))       #min el power limitation
m.addConstrs(P_fc[z] <= P_fc_nom     * k[z,1] for z in range(n))       #max fc power limitation
m.addConstrs((k[z,0] + k[z,1] == 1 for z in range(n))) #avoid simultaneous production and consumption of H2 

m.addConstr( H2_tank[0]   == tank_capacity * 0.4)     #takes the soc of the previous simulation
m.addConstr( H2_tank[n-1] == tank_capacity * 0.4)     #takes the soc of the previous simulation

m.addConstrs(H2_tank[z+1] == H2_tank[z] + P_el[z]*phi_el*eta_comp - P_fc[z]*phi_fc for z in range(n-1)) #soc update considering chage and discharge efficiency
m.addConstrs(H2_tank[z]   <= tank_capacity for z in range(n))

'BESS constraints'
m.addConstr( En_BESS[0]== BESS_capacity*0.5)     #takes the soc of the previous simulation
# m.addConstr( En_BESS[n-1] >= BESS_capacity*0.3)      #tries to reach SOC = 50% at the end of the second day

m.addConstrs( (En_BESS[z+1] == En_BESS[z] + (P_BESS_p[z]/eta_d - P_BESS_n[z]*eta_c) / kWh_factor for z in range(n-1))) #soc update considering chage and discharge efficiency
m.addConstrs( En_BESS[z] <= BESS_capacity*0.95 for z in range(n)) #soc update considering chage and discharge efficiency
m.addConstrs( En_BESS[z] >= BESS_capacity*0.15 for z in range(n)) #soc update considering chage and discharge efficiency

m.addConstrs(P_BESS_p[z]   <= BESS_capacity * C_rate_discharge     * h[z,0] for z in range(n))   #discharge power limitation
m.addConstrs(P_BESS_n[z]  <= BESS_capacity * C_rate_charge  * h[z,1] for z in range(n))      #charge power limitation
m.addConstrs((h[z,0] + h[z,1] == 1 for z in range(n)))   #cannot charge and discharge simultaneously 

m.addConstrs ((P_BESS_tot[z] == P_BESS_p[z] + P_BESS_n[z]  for z in range(n)))

# for z in range(n):
#     m.addGenConstrAbs(P_BESS_abs[z], P_BESS_tot[z], "absconstr")

'global energy balance'
# m.addConstrs( (WT_pot[z] + PV_pot[z] + P_fc[z] + P_BESS_p[z] + P_grid[z] == load[z] + P_el[z] + P_BESS_n[z] + P_excess [z] for z in range(n)), 'energy balance')

m.addConstrs( (WT_pot[z] + PV_pot[z] + P_fc[z] + P_BESS_p[z] + P_grid[z] ==     #+ P_var_n[z]
               load[z] + P_el[z] + P_BESS_n[z] + P_exp[z] + P_excess [z] 
               for z in range(n)), 'energy balance')

'load request'
m.addConstrs( P_grid[z]  <= P_grid_max for z in range(n))  #limit power absorption from the grid

m.addConstrs( (WT_load[z] + PV_load[z] + P_fc[z] + P_BESS_p_load[z] + P_grid[z] == load[z] for z in range(n)), 'load energy balance')

m.addConstrs( (P_load_covered[z] == load[z] - P_grid[z] for z in range(n)), 'load energy balance')

'export - energy balances'
m.addConstrs( (WT_band[z] + WT_over[z] + WT_load[z] + WT_excess [z] + WT_bess[z] + WT_el[z] == WT_pot[z] for z in range(n)), 'wind energy balance')  #non è uguale, può andare anche dentro alla bess o fc
m.addConstrs( (PV_band[z] + PV_over[z] + PV_load[z] + PV_excess[z] + PV_bess[z] + PV_el[z] == PV_pot[z] for z in range(n)), 'pv energy balance')

m.addConstrs( (WT_bess[z] + PV_bess[z] == P_BESS_n[z] for z in range(n)), 'BESS input balance')
m.addConstrs( (WT_el[z] + PV_el[z] == P_el[z] for z in range(n)), 'EL input balance')

m.addConstrs( (WT_band[z] + WT_over[z] + PV_band[z] + PV_over[z] + P_BESS_p_grid[z] == P_exp[z] for z in range(n)), 'export energy balance')
m.addConstrs( (WT_excess[z] + PV_excess[z] == P_excess[z] for z in range(n)), 'excess energy balance')

m.addConstrs( P_exp[z]  <= P_grid_max for z in range(n))  #limit power export to the grid


m.addConstrs( (P_BESS_p_load[z] + P_BESS_p_grid[z] == P_BESS_p[z]  for z in range(n)), 'bess discharge balance')

'fix band - grid request'
m.addConstrs( (P_fix[z] == WT_band[z] + PV_band[z] + P_BESS_p_grid[z]   for z in range(n)), 'grid energy balance')   #+ P_var_n[z]

m.addConstr( R_0 ==   gp.quicksum(WT_band[z]       * P_wind_band for z in range(n)) 
                    + gp.quicksum(PV_band[z]       * P_pv_band   for z in range(n))   
                    + gp.quicksum(WT_over[z]       * P_wind_over for z in range(n)) 
                    + gp.quicksum(PV_over[z]       * P_pv_over   for z in range(n))   
                    + gp.quicksum(P_BESS_p_grid[z] * P_bess      for z in range(n))      
                    
                    # minimize undershootings
                    # - gp.quicksum(P_var_n[z] * C_pen_n for z in range(n)) 
                    
                    # minimize components usage
                    # - gp.quicksum(P_el[z] * C_el for z in range(n)) 
                    # - gp.quicksum(P_fc[z] * C_fc for z in range(n)) 
                    # - gp.quicksum(P_BESS_abs[z] * C_batt for z in range(n)) 
                    # - gp.quicksum(P_excess[z] * C_excess for z in range(n)) 
                    
                    # savings from covering the load using RES instead of grid power
                    + gp.quicksum(P_load_covered[z] * C_grid for z in range(n))
                    )

m.addConstrs( R_t[0] == 0 for j in range(1,years))
m.addConstrs( R_t[j] == R_0 for j in range(1,years))

m.addConstrs( OEM_t[0] == 0 for j in range(1,years))
m.addConstrs( OEM_t[j] == OEM_0 for j in range(1,years))

m.addConstrs( C_t[10] == C_0 for j in range(1,years))
m.addConstrs( C_t[j] == 0 for j in range(0,9))
m.addConstrs( C_t[j] == 0 for j in range(11,years))

m.addConstr( I0 ==  BESS_capacity * bess_cost
                    + P_el_nom * EL_cost 
                    + P_fc_nom * FC_cost 
                    + tank_capacity * HP_tank_cost 
                    # + PV_scale * 160 * PV_cost 
                    # +  800 * WT_cost
            )

m.addConstr( OEM_0 == I0 * 0.08)

m.addConstr( NPV ==   gp.quicksum(R_t[j]   *  d_vector[j] for j in range(years)) 
                    - gp.quicksum(OEM_t[j] *  d_vector[j] for j in range(years))
                    # - gp.quicksum(C_t[j]   *  d_vector[j] for j in range(years))
                    - I0
                    )

'''NPV Objective function'''  
m.setObjective(NPV, GRB.MAXIMIZE)
# m.setObjective(R_0, GRB.MAXIMIZE)

'optimization'
m.optimize()   

#%%
'hourly flows outpout'
df_h = pd.DataFrame()

df_h['WT_prod [kW]']  = WT_pot
df_h['PV_prod [kW]']  = [PV_pot[z].x for z in range(hours)]
df_h['RES_prod [kW]']  = df_h['WT_prod [kW]'] + df_h['PV_prod [kW]']
df_h['load [kW]'] = load

df_h['p_var_p']   = [P_var_p[z].x for z in range(hours)]
df_h['p_var_n']    = [P_var_n[z].x for z in range(hours)]
df_h['p_fix']      = [P_fix[z].x for z in range(hours)]
df_h['h_tot']      = [h[z,0].x for z in range(hours)]
   
df_h['e_wind_band'] = [WT_band[z].x for z in range(hours)]
df_h['e_pv_band']   = [PV_band[z].x for z in range(hours)]
df_h['e_wind_over'] = [WT_over[z].x for z in range(hours)]
df_h['e_pv_over']   = [PV_over[z].x for z in range(hours)]

df_h['soc_BESS']  = [En_BESS[z].x for z in range(hours)]
df_h['soc_tank']  = [H2_tank[z].x for z in range(hours)]
df_h['P_el']      = [P_el[z].x for z in range(hours)]
df_h['P_fc']      = [P_fc[z].x for z in range(hours)]

df_h['P_excess']  = [P_excess[z].x for z in range(hours)]
df_h['P_grid']  = [P_grid[z].x for z in range(hours)]

df_h['P_BESS_tot']  = [P_BESS_tot[z].x for z in range(hours)]
df_h['P_BESS_p']   = [P_BESS_p[z].x for z in range(hours)]
df_h['P_BESS_n']    = [P_BESS_n[z].x for z in range(hours)]
df_h['P_BESS_abs']  = [P_BESS_abs[z].x for z in range(hours)]

df_h['P_BESS_grid']   = [P_BESS_p_grid[z].x for z in range(hours)]

df_h['P_load_covered']   = [P_load_covered[z].x for z in range(hours)]
df_h['WT_load']   = [WT_load[z].x for z in range(hours)]
df_h['PV_load']   = [PV_load[z].x for z in range(hours)]
df_h['BESS_load']   = [P_BESS_p_load[z].x for z in range(hours)]

df_h['WT_bess']   = [WT_bess[z].x for z in range(hours)]
df_h['PV_bess']   = [PV_bess[z].x for z in range(hours)]
df_h['WT_el']   = [WT_el[z].x for z in range(hours)]
df_h['PV_el']   = [PV_el[z].x for z in range(hours)]

#%%
#output sizing

if MILP_design == True:
    df_sizes = pd.DataFrame()
    df_sizes['BESS_capacity'] = [BESS_capacity.x]
    df_sizes['P_el_nom'] = [P_el_nom.x]
    df_sizes['P_fc_nom'] = [P_fc_nom.x]
    df_sizes['tank_capacity'] = [tank_capacity.x]
    df_sizes['PV_scale'] = [PV_scale.x]
    
    df_sizes['NPV'] = [NPV.x]
    df_sizes['I0'] = [I0.x]
    
    df_NPV = pd.DataFrame()
    df_NPV['R_t'] = [R_t[i].x for i in range(len(R_t))]
    df_NPV['OEM_t'] = [OEM_t[i].x for i in range(len(OEM_t))]
    df_NPV['C_t'] = [C_t[i].x for i in range(len(OEM_t))]
    
    df_NPV['R_t_d'] = df_NPV['R_t'] * d_vector
    df_NPV['OEM_t_d'] = df_NPV['OEM_t'] * d_vector
    df_NPV['C_t_d'] = df_NPV['C_t'] * d_vector
    
    I0_test = (df_sizes['BESS_capacity'] * bess_cost 
                + df_sizes['P_el_nom'] * EL_cost 
                + df_sizes['P_fc_nom'] * FC_cost 
                + df_sizes['tank_capacity'] * HP_tank_cost
                + df_sizes['PV_scale'] * 160 * PV_cost 
                + 800 * WT_cost
                + comp_cost
                )
    
    NPV_test = sum(df_NPV['R_t_d']) - sum(df_NPV['OEM_t_d']) - sum(df_NPV['C_t_d']) - I0_test
    
    NPV_list = []
    
    NPV_y = -I0_test[0]
    
    for y in range(years): 
        NPV_y = NPV_y + df_NPV['R_t_d'][y] + df_NPV['OEM_t_d'][y] + df_NPV['C_t_d'][y]
        NPV_list.append(NPV_y)
        
    df_NPV['NPV'] = NPV_list

#%%
df_results = pd.DataFrame()
df_results['Earings [€/y]'] = [R_0.x]

df_results['e_wind_band'] = [sum(df_h['e_wind_band'])]
df_results['e_pv_band']   = [sum(df_h['e_pv_band'])]
df_results['e_wind_over'] = [sum(df_h['e_wind_over'])]
df_results['e_pv_over']   = [sum(df_h['e_pv_over'])]

df_results['soc_BESS']  = [sum(df_h['soc_BESS'])]
df_results['soc_tank']  = [sum(df_h['soc_tank'])]
df_results['e_el']      = [sum(df_h['P_el'])]
df_results['e_fc']      = [sum(df_h['P_fc'])]

df_results['e_excess'] = [sum(df_h['P_excess'])] 
df_results['e_grid']   = [sum(df_h['P_grid'])]

df_results['e_BESS_tot']  = [sum(df_h['P_BESS_tot'])]
df_results['e_BESS_p']    = [sum(df_h['P_BESS_p'])]
df_results['e_BESS_n']    = [sum(df_h['P_BESS_n'])]
df_results['e_BESS_abs']  = [sum(df_h['P_BESS_abs'])]

df_results['e_BESS_grid']  = [sum(df_h['P_BESS_grid'])]
df_results['e_load_covered'] = [sum(df_h['P_load_covered'])]

df_results['R_wind_band'] = df_results['e_wind_band'] * P_wind_band
df_results['R_pv_band']   = df_results['e_pv_band'] * P_pv_band
df_results['R_wind_over'] = df_results['e_wind_over'] * P_wind_over
df_results['R_pv_over']   = df_results['e_pv_over'] * P_pv_over
df_results['R_BESS_grid'] = df_results['e_BESS_grid'] * P_bess 
df_results['R_load_covered'] = df_results['e_load_covered'] * C_grid

#%%
df_NPV.to_csv('NPV.csv', sep = ';')
df_results.to_csv('results.csv', sep = ';')
df_h.to_csv('hourly.csv', sep = ';')
df_sizes.to_csv('sizes.csv', sep = ';')





