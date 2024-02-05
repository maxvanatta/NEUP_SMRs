## V24 Including steam release (inequaltiy in the capacity constraint)
# Also preventing market runaway by limitng the the peak demand

## V23 Streamlining for Greatlakes

## V22 Improving the day-by-day DAM suboptimizations into a whole year rather than full year.
# Adds in initial states into the formulation (lb, ub # mods, reactor gen, elec gen, one/off state)
# Takes a first pass with reduced hours to find the number of modules needed 
# Then fixes the number of modules for each daily optimization (need to add a lower bound of number lb = ub)
# 

## V21 Adding in subsets of the year start date and number of hours

## V20 Adding in the hourly profiles. will do this using the naics code

## V19 Making the script more able to handle the per facility format. Demands in particular are now not a csv which is brough in, but a single hourly value. 06/13/2022

## V18 Dropping the Avoided cost of ff in the objective function. Add in more post process for comparing and analizing revenue breakdown 06/09/2022

## V17 More cleanup

## V16 Cleaning up some constraintes based upon the writing. Making the total charge and discharge for the TES for TES scale not for the modules.

## V15 Making TES paramters optimization -- Lied, did not do that. Just updated some things

## V14 Removing Fast solve as it was far slower. And inputing TES without optimization

## V13 This is a short-cut for resolving results more quickly, by brute force

## V12 Adding the fast solve functionality with solely a float value number of modules option. 
# Adding this possibility on top of current layout to provide layered effect.
# Can run fast, choose best 3, then run full analysis.

## V11 Adding in thermal penalties for coolant to working fluid transfer. 01/21/22
# This is primarily done in the Generator spreadsheet 'sites.csv'
# The Temperature coefficient is 0.75 for all systems currently and are adjusted to occur prior to importing.
# Energy loss to the heat transfer are adjusted to the maximum thermal capacity BUT the electric capacity already accounts for this in the efficiceny therefore it is only applied to the pTCap value.

## V10 Fixing the MSL for thermal baseline issue.
# This means that the vTGen is now the total thermal generation from the 'boiler' aspect of the cycle.
# vGen remains unchanged
# The thermal demand for the facility is satisfied by (vTGen-vGen/pEff)
# This allows for the system to have a more specific link between the thermal requirements of the nuclear reactor itself which has a distinct ramping and range
# and for the electricicty to behave a certain way without tryign to directly couple them in vTGenabovemin which would provide complications if the 
# vTGen were simply the heat to industry.
# new value is put into the results vHeatGen which is the MWht for the supplied industrial heating.
# this opens the question is the pTVOC is solely the thermal cycle and if the pVOC is just to operate the turbine and generate electricty outside of the thermal behavior.

import pandas as pd
import numpy as np
import cplex
import os

import matplotlib.pyplot as plt
import seaborn as sns
import math

import plotly.graph_objects as go
import plotly.io as pio
import plotly.express as px

import pyomo.environ as pe
import pyomo.opt
from pyomo.opt import SolverStatus

from datetime import datetime
s = datetime.now()
print(s,'start')

import argparse


parser = argparse.ArgumentParser(description='Import the input variables for the model run')
parser.add_argument('--Start', type=int, help='Starting facility index',required=True)
parser.add_argument('--End', type=int, help='Ending facility index',required=True)
parser.add_argument('--Region', type=str, help='SPP or ERCOT',required=True)
parser.add_argument('--Year', type=int, help='Year',required=True)
parser.add_argument('--Zero', type=int, help='no electricity market',required=False)
parser.add_argument('--NoHeat', type=bool, help='no Heat demand',required=False)
parser.add_argument('--NoDAM_EACH', type=bool, help='no DAM, purely looking at the wholesale possibilities in each market',required=False)
parser.add_argument('--NoDAM', type=bool, help='no DAM, purely looking at the wholesale possibilities in each market',required=False)
parser.add_argument('--Sens', type=str, help='Which sensitivity to run',required=False)
parser.add_argument('--LimitedMod', type=bool, help='Limits Module maximum to peak demand +10%, Default is to limit',required=False)
# sensitivities can inlcude: Market Year, TES, Hourly Variation, Market Foresight


args = parser.parse_args()



import Facility_Processing_2015 as FP

class EconomicAssessment_OneGen:
    def __init__(self):
        
        self.Generator_IDs = None             # Names of the egenrators being evaluated
        
        #### sets ####
        self.G = None                         # Generator indexes [0,1,2....]
        self.T = None                         # Timestep indexes [0,1,2....8759]
        self.M = None #list(range(6))              # Maximum number of modules for any generator
        #self.YR = None                        # Number of years
        
        #### paramteres ####
        
        # Environmental paramters/System
        self.pLMP = None                      # LMP electric $MWhe [t] or [t,yr]
        #self.pTLMP = None                     # LMP heat $/MWht [g,t] or [t,yr]
        self.pTempdemanded = None             # Minimum temperature which the reactor must supply for the given TLMP single value
        self.pTEnergyDemand = None            # Thermal energy required per hour in MWt [t]
        
        # Engineering parameters for generators
        self.pOnoroffinitial = None           # initial state, default, 1 = On [g]
        self.pStartupfixedcost = None         # $/MW [g]
        self.pVOC = None                      # Variable O&M ($/MWhe) [g]
        self.pTVOC = None                     # Variable O&M ($/MWht) [g]
        self.pCap = None                      # Maximum Electric capacity MWe == pTCap*ThermalEff [g]
        self.pTCap = None                     # Maximum Thermal capacity MWt [g]
        self.pThermaleff = None               # Thermal efficiency [g]
        self.pThermaTransfereff = None        # V11 Thermal transsfer efficicnecy between coolant and workign fluid.
        self.pMSL = None                      # Minimum stable thermal load (MWt) [g]
        self.pMSLTurb = None                  # Minimum stable trubine load minumum core + turbine power [g]
        self.pMDT = None                      # Minimum down time (hours) [g]
        self.pGenabovemininitial = None       # Initial generation above the min (MWe) [g]
        self.pRGenInitial = None              # Initial generation of reactor (MWth) [g]
        self.pRamprate = None                 # Maximum ramp (up and down) (MW/hr) [g]
        
        self.pCAPEX = None                    # Capital costs of the generator $/MWe [g]
        self.pFOPEX = None                    # Fixed O&M costs of the generator $/MWe [g]
        self.pModulesize = None               # Module size in MWe [g]
        self.pOutlettemp = None               # SMR outlet temperature [g]
        self.pWorkingTemp = None              # SMR working fluid temp 
        self.pMaxMod = None                   # maximum number of modules [g]
        self.pMinMod = None                   # minimum number of modules [g] V222 71522    
        
        self.pIR = 0.07                       # Discount/ Interest Rate 7% (OMB Circular)
        self.pLifetime = 30                   # Lifetime, years
        
        # V14 Additions MV 3/1/2022
        ## TES Parameters ####
        self.pTESCapacity = None              # Storage thermal capacity (MWh-thermal)
        self.pTESChargeEff = None             # Charging Efficiency (%) 
        self.pTESDischargeEff = None          # Discharging Efficicency (%) 
        self.pTESLossRate = None              # Rate at which the stored thermal energy is lost (%/hr)
        self.pTESPowerCap = None              # Power in or out maximum (MWt)
        
        ## TES Variables #### 
        self.vTESSOC = None                   # TES State of Charge (MWht)
        self.vTESSOC_Charge = None            # TES Charge amount (MWht)
        self.vTESSOC_Discharge = None         # TES Discharge amount (MWht)
        # V15
        self.vTESCapacity = None              # DV of the total storage (MWht)
        self.vTESPowerCap = None              # DV of CHarge/Discharge Power (MWt)
        # V22 
        self.pTESInitSOC = None               # Start SOC of TES (MWht)
        self.pTES_Ch0 = None
        self.pTES_Dis0 = None
        #### variables ####
        
        # Float variables
        self.vGen = None                      # Electric generation (MWhe) [g,t] --> [t]
        self.vTGen = None                     # Heat generation (MWht) [g,t] --> [t]
        self.vGenabovemin = None              # Electric generation above the MSLTurb (MWhe) [m,t]
        self.vTGenabovemin = None             # Thermal Generation above MSL (MWht) [m,t]
        
        # Binary variables
        self.vTurnon = None                   # Does the generator turn on (1) or not (0) [m,t] --> [t]
        self.vTurnoff = None                  # Gnerator turns off(1) or not (0) at time t [m,t] --> [t]
        self.vOnoroff = None                  # Generator state at time t [m,t] --> [t]
        self.vEOnoroff = None                 # Electricity generation on (1) or off (0) [m,t] --> [t]
        self.vTOnoroff = None                 # Heat generation on (1) or off (0) [m,t] --> [t]
        
        # results/objective variables  --These are not explicit model components and ar calculated from other variable values
        self.rEProfit = None                  # Profits from vGen*LMP [m,t]
        self.rTProfit = None                  # Profits from vTgen*TLMP [m,t]
        self.rProfit = None                   # Hourly profits ($) [m,t]
        
        self.vTotalProfits = None             # Total Profits ($) [-] single value
        self.rFeasible = None                 # is this configurations feasible     
        
        self.vGenerator = None                # Index of the generator chosen for the
        
        #V18 
        self.rACC = None                      # Annualized capital cost
        self.rFOMC = None                     # Fixed cost
        self.rVC = None                       # Variable cost
        self.rFC = None                       # Fuel cost
        self.rEP = None                       # Electricity profit
        self.rAVC = None                      # Avoided cost of natural gas
        self.rSUC = None                      # Startup costs
        self.rTES = None                      # TES CAPEX
        self.RampCost = 1.0
        
        
    def BuildModel(self, DAM = False):
        pTS = []
        for i in range(len(self.G)):
            if self.pWorkingTemp[i] >= self.pTempdemanded:
                pTS.append(1)
            else:
                pTS.append(0)
            
        self.pThermalSat = pTS
        
        if not DAM:
            print('Building Model')
        model = pe.ConcreteModel()
        
        #### SETS ####     
        
        model.G = pe.Set(initialize = self.G, doc = "Generator IDs")
        model.T = pe.Set(initialize = self.T, doc = "Timestep")
        model.M = pe.Set(initialize = self.M, doc = "Modules")
        #model.YR = pe.Set(initialize = self.YR, doc = 'Years of model runs-- CURRENTLY NOT USED')
        
        
        #### PARAMETERS ####
        
        model.pLMP = pe.Param(model.T,initialize = {self.T[i]: self.pLMP[i] for i in range(len(self.T))} , doc = "LMP Electric $/MWhe")
        #model.pTLMP = pe.Param(model.G, model.T,initialize = self.pTLMP, doc = "LMP Thermal $/MWht")
        model.pTEnergyDemand = pe.Param(model.T,initialize = {self.T[i]: self.pTEnergyDemand[i] for i in range(len(self.T))} , doc = "Energy demand MWht")
        
        model.pOnoroffinitial = pe.Param(model.M,initialize = {self.M[i]: self.pOnoroffinitial[i] for i in range(len(self.M))}, doc = "Inital state, 1 on, 0 off")
        model.pStartupfixedcost = pe.Param(model.G,initialize = {self.G[i]: self.pStartupfixedcost[i] for i in range(len(self.G))}, doc = "Startup costs, $500 default")
        #model.pVOC = pe.Param(model.G,initialize = {self.G[i]: self.pVOC[i] for i in range(len(self.G))}, doc = "Variable O&M Electric $/MWhe")
        model.pTVOC = pe.Param(model.G,initialize = {self.G[i]: self.pTVOC[i] for i in range(len(self.G))}, doc = "Variable O&M Thermal $/MWhe")
        model.pCap = pe.Param(model.G,initialize = {self.G[i]: self.pCap[i] for i in range(len(self.G))} , doc = "Maximum Electric Capacity MWe")
        model.pTCap = pe.Param(model.G,initialize = {self.G[i]:self.pTCap[i] for i in range(len(self.G))}, doc = "Maximum Thermal Capacity MWt")
        model.pThermaleff = pe.Param(model.G,initialize = {self.G[i]: self.pThermaleff[i] for i in range(len(self.G))}, doc = "Thermal efficiency (MWe/MWt)")
        model.pThermalTransfereff = pe.Param(model.G,initialize = {self.G[i]: self.pThermaTransfereff[i] for i in range(len(self.G))}, doc = "Thermal transfer efficiency (MW -coolant/MWt- working fluid)")
        model.pMSL = pe.Param(model.G,initialize = {self.G[i]: self.pMSL[i] for i in range(len(self.G))}, doc = "Minimum Stable Load Thermal MWt")
        model.pMSLTurb = pe.Param(model.G,initialize = {self.G[i]: self.pMSLTurb[i] for i in range(len(self.G))}, doc = "Minimum Stable Load of Turbine MWe")
        model.pMDT = pe.Param(model.G,initialize = {self.G[i]: self.pMDT[i] for i in range(len(self.G))} , doc = "Mimimum Down time (hours)")
        #model.pGenabovemininitial = pe.Param(model.M,initialize = {self.M[i]: self.pGenabovemininitial[i] for i in range(len(self.M))}, doc = "Initial Generation above the Minimun MWe")
        model.pRGenInitial = pe.Param(model.M,initialize = {self.M[i]: self.pRGenInitial[i] for i in range(len(self.M))}, doc = "Initial Generation MWt")
        model.pRamprate = pe.Param(model.G,initialize = {self.G[i]: self.pRamprate[i] for i in range(len(self.G))}, doc = "Ramp rate up and down MW/hr")        
        model.pFOPEX = pe.Param(model.G,initialize = {self.G[i]: self.pFOPEX[i] for i in range(len(self.G))}, doc = "FOPEX pparamter $/MWe")
        model.pCAPEX = pe.Param(model.G,initialize = {self.G[i]: self.pCAPEX[i] for i in range(len(self.G))}, doc = "CAPEX values$/MWe")
        model.pWorkingTemp = pe.Param(model.G,initialize = {self.G[i]: self.pWorkingTemp[i] for i in range(len(self.G))}, doc = "working temperature which the reactor can supply in Degrees C")
        model.pThermalSat  = pe.Param(model.G,initialize = {self.G[i]: self.pThermalSat[i] for i in range(len(self.G))}, doc = "SMR can supply the temperature")
        model.pMaxMod = pe.Param(model.G, initialize = {self.G[i]: self.pMaxMod[i] for i in range(len(self.G))}, doc = "maximum number of modules")
        #V22 MV 7/15/2022
        model.pMinMod = pe.Param(model.G, initialize = {self.G[i]: self.pMinMod[i] for i in range(len(self.G))}, doc = "maximum number of modules")
        model.pTESInitSOC = pe.Param(model.G, initialize = {self.G[i]: self.pTESInitSOC[i] for i in range(len(self.G))}, doc = "initial SOC")
        if DAM:
            model.pTES_Ch0 = pe.Param(model.G, initialize = {self.G[i]: self.pTES_Ch0[i] for i in range(len(self.G))}, doc = "initial charge DAM")
            model.pTES_Dis0 = pe.Param(model.G, initialize = {self.G[i]: self.pTES_Dis0[i] for i in range(len(self.G))}, doc = "initial discharge DAM")

        
        # V14 Additions MV 3/1/2022
        model.pTESCapacity = pe.Param(model.G,initialize = {self.G[i]: self.pTESCapacity[i] for i in range(len(self.G))}, doc = "TES max storage")
        model.pTESChargeEff = pe.Param(model.G,initialize = {self.G[i]: self.pTESChargeEff[i] for i in range(len(self.G))}, doc = "TES Charge Efficiency")
        model.pTESDischargeEff = pe.Param(model.G,initialize = {self.G[i]: self.pTESDischargeEff[i] for i in range(len(self.G))}, doc = "TES Discharge Efficiency")
        model.pTESLossRate = pe.Param(model.G,initialize = {self.G[i]: self.pTESLossRate[i] for i in range(len(self.G))}, doc = "TES Loss Rate if just sitting")
        model.pTESPowerCap = pe.Param(model.G,initialize = {self.G[i]: self.pTESPowerCap[i] for i in range(len(self.G))}, doc = "Maximum power in or out of stroage")
        
        #### VARIABLES ####
        model.vGenerator = pe.Var(model.G,within=pe.Binary,doc = "Choice of generator")
        model.vModuleS = pe.Var(model.G,model.M,within = pe.Binary, doc = "choice of modules for full optimization")
        
        model.vGen = pe.Var(model.G,model.M,model.T,within=pe.NonNegativeReals, doc = "Electric generation (MWhe)")
        model.vTGen = pe.Var(model.G,model.M,model.T,within=pe.NonNegativeReals, doc = " Heat generation (MWht) ")
        model.vRGen = pe.Var(model.G,model.M,model.T,within=pe.NonNegativeReals, doc = " Total Reactor generation (MWht) ")
        
        model.vGenabovemin = pe.Var(model.G,model.M,model.T,within=pe.NonNegativeReals, doc = "Electric generation above the MSLTurb (MWhe)")
        model.vRGenabovemin = pe.Var(model.G,model.M,model.T,within=pe.NonNegativeReals, doc = "Thermal Generation above MSL (MWht)")
        
        model.vTurnon = pe.Var(model.G,model.M,model.T,within=pe.Binary, doc = "Does the generator turn on (1) or not (0)")
        model.vTurnoff = pe.Var(model.G, model.M,model.T,within=pe.Binary, doc = "Generator turns off(1) or not (0) at time t")
        model.vOnoroff = pe.Var(model.G,model.M,model.T,within=pe.Binary, doc = "Generator state at time t")
        model.vEOnoroff = pe.Var(model.G,model.M,model.T,within=pe.Binary, doc = "Electricity generation on (1) or off (0)")
        model.vTOnoroff = pe.Var(model.G,model.M,model.T,within=pe.Binary, doc = "Heat generation on (1) or off (0)")
        
        # V14 Additions MV 3/1/2022
        model.vTESSOC = pe.Var(model.G,model.T,within=pe.NonNegativeReals, doc = "TES State of Charge")
        model.vTESSOC_Charge = pe.Var(model.G,model.M,model.T,within=pe.NonNegativeReals, doc = "TES SOC charge")
        model.vTESSOC_Discharge = pe.Var(model.G,model.T,within=pe.NonNegativeReals, doc = "TES SOC charge")
        
        # V15 
        #model.vTESCapacity = pe.Var(model.G,within=pe.NonNegativeReals, doc = "DV of TES Capacity")
        #model.vTESPowerCap = pe.Var(model.G,within=pe.NonNegativeReals, doc = "DV of TES Power Charge/Discharge")
        
        if not DAM:
            print('Params & Vars established')
        
        
        #### OBJECTIVES ####
        
        def ERevenues(model):
            return sum(sum(sum(model.vGen[g,m,t]*model.pLMP[t]for g in model.G) for m in model.M) for t in model.T)
        
        #def TRevenues(model):
        #    return sum(sum((sum(model.vTGen[g,m,t] for m in model.M)+ model.vTESSOC_Discharge[g,t]*model.pTESDischargeEff[g])*model.pTLMP[g,t] for g in model.G) for t in model.T)
        
        def RampCost(model):
            return sum(sum(sum((((model.vRGenabovemin[g,m,model.T.prevw(t)])-(model.vRGenabovemin[g,m,t]))**2)*self.RampCost for t in model.T) for m in model.M) for g in model.G)
        
        def Costs(model):
            return (sum(model.pStartupfixedcost[g]*sum(sum(model.vTurnon[g,m,t] for t in model.T) for m in model.M) for g in model.G) 
                    + sum(sum(sum((model.vRGen[g,m,t]*model.pTVOC[g]) for g in model.G) for m in model.M) for t in model.T) 
                    + sum(sum(model.pCap[g]*(model.pCAPEX[g]*model.vModuleS[g,m]*(self.pIR/(1-((1+self.pIR)**(-1*self.pLifetime))))+model.pFOPEX[g]*model.vModuleS[g,m]) for m in model.M) for g in model.G)*(len(self.T)/8760)
                    + sum(model.vGenerator[g]*model.pTESCapacity[g]*self.pTESCapex*(self.pIR/(1-((1+self.pIR)**(-1*self.pLifetime))))*(len(self.T)/8760) for g in model.G))
        
        
        def Obj_Profit(model):
            return ERevenues(model)-Costs(model)-RampCost(model)  # +TRevenues(model)
        model.Obj_Profit = pe.Objective(rule = Obj_Profit, sense=pe.maximize, doc = "Maximize the profits by balancing thermal and electric generation")
                

        #### CONSTRAINTS ####
        
        model.GeneratorChooser = pe.Constraint(expr = sum(model.vGenerator[g] for g in model.G)<=1)
        model.ModChooserS = pe.Constraint(expr = sum(sum(model.vModuleS[g,m] for g in model.G)for m in model.M)>=1)

        def ThermalEquality(model,t):
            return sum(sum(model.vTGen[g,m,t] for m in model.M) + model.vTESSOC_Discharge[g,t]*model.pTESDischargeEff[g] for g in model.G) == model.pTEnergyDemand[t]
        model.ThermalEquality = pe.Constraint(model.T, rule = ThermalEquality, doc = "Limits the thermal output (and therefore the 'revenue from heat')")
        
        def ModuleLimit(model,g):
            return sum(model.vModuleS[g,m] for m in model.M) <= model.pMaxMod[g]
        model.ModuleLimit = pe.Constraint(model.G,rule = ModuleLimit, doc = "limits the number of available modules")
        
        def ModuleLimitLower(model,g):
            return sum(model.vModuleS[g,m] for m in model.M) >= model.pMinMod[g]
        model.ModuleLimitLower = pe.Constraint(model.G,rule = ModuleLimitLower, doc = "limits the lower number of available modules")
        
        def ModGenS(model,g,m):
            return model.vModuleS[g,m] <= model.vGenerator[g]
        model.ModGenS = pe.Constraint(model.G,model.M,rule = ModGenS, doc = "cannot build modules for the non-optimal gens Slow")
        
        def ModStatus(model,g,m,t):
            return model.vOnoroff[g,m,t] <= model.vModuleS[g,m]
        model.ModStatus = pe.Constraint(model.G,model.M, model.T, rule = ModStatus, doc = "The status of the plant must be the previous state plus the turnon/off")
        
        def ECapacityLimit(model,g,m,t):
            return model.vGen[g,m,t] <= model.pCap[g]
        model.ECapacityLimit = pe.Constraint(model.G,model.M,model.T,rule = ECapacityLimit, doc = "The sum of the electric generation must not surpass maximum generation")
        
        def RGenEquality(model,g,m,t):
            return (model.vTGen[g,m,t] + model.vGen[g,m,t]/model.pThermaleff[g] + 
                    model.vTESSOC_Charge[g,m,t]/model.pTESChargeEff[g]) <= model.vRGen[g,m,t]*model.pThermalTransfereff[g] # Updated V14 MV 03012022 # Updated to inequaltiy in V24 10/26/2022
        model.RGenEquality = pe.Constraint(model.G,model.M,model.T,rule = RGenEquality, doc = "electricity and thermal must equal ttoal")
        
        def EStatus(model,g,m,t): 
            return model.vEOnoroff[g,m,t] <= model.vOnoroff[g,m,t]
        model.Estatus = pe.Constraint(model.G,model.M,model.T, rule = EStatus, doc = "The electric generation cannot be on when the overall status is off")
        
        def Estatus2(model,g,m,t):
            return model.vGen[g,m,t] <= model.pCap[g]*model.vEOnoroff[g,m,t]
        model.Estatus2 = pe.Constraint(model.G,model.M,model.T,rule = Estatus2, doc = "The sum of the electric generation must not surpass maximum generation")
        
        def RStatus(model,g,m,t):
            return model.vRGen[g,m,t]*model.pThermaleff[g]<= model.pCap[g]*model.vOnoroff[g,m,t]
        model.RStatus = pe.Constraint(model.G,model.M,model.T,rule = RStatus, doc = "The sum of the electric and thermal generation must not surpass maximum generation")
        
        def TAllow(model,g,m,t):
            return model.vOnoroff[g,m,t] <= model.pThermalSat[g]
        model.TAllow = pe.Constraint(model.G,model.M,model.T, rule = TAllow, doc = "Must satisfty heat constraint")
        
        def RGenAbove(model,g,m,t):
            return model.vRGenabovemin[g,m,t]  == ((model.vRGen[g,m,t])-(((model.pMSL[g])/model.pThermaleff[g])*model.vOnoroff[g,m,t]))
        model.RGenAbove = pe.Constraint(model.G,model.M,model.T, rule = RGenAbove, doc = "This defines the generation above the minimum stable load")
        
        def EGenAbove(model,g,m,t):
            return model.vGenabovemin[g,m,t]  == (model.vGen[g,m,t]-(model.pMSLTurb[g]*model.vEOnoroff[g,m,t]))
        model.EGenAbove = pe.Constraint(model.G,model.M,model.T, rule = EGenAbove, doc = "This defines the generation above the minimum stable load")
        
        def RRamprateUp(model,g,m,t):
            if t == model.T[1]: 
                return (model.pRGenInitial[m])-(model.vRGenabovemin[g,m,t]) >= -1*model.pRamprate[g]/(model.pThermaleff[g])
                #return pe.Constraint.Feasible # V11C THis means that we can start iur reactor at any output.
            return (model.vRGenabovemin[g,m,model.T.prev(t)])-(model.vRGenabovemin[g,m,t]) >= -1*model.pRamprate[g]/(model.pThermaleff[g])
        model.RRamprateUp = pe.Constraint(model.G,model.M,model.T, rule = RRamprateUp, doc = "The hourly change must not exceed the ramprate between time steps")
        
        def RRamprateDn(model,g,m,t):
            if t == model.T[1]:
                return ((model.pRGenInitial[m])-(model.vRGenabovemin[g,m,t]))*(model.pThermaleff[g]) <= model.pRamprate[g]
                #return pe.Constraint.Feasible # V11C THis means that we can start iur reactor at any output.
            return (model.vRGenabovemin[g,m,model.T.prev(t)]-model.vRGenabovemin[g,m,t])*(model.pThermaleff[g]) <= model.pRamprate[g]
        model.RRamprateDn = pe.Constraint(model.G,model.M,model.T, rule = RRamprateDn, doc = "The hourly change must not exceed the ramprate between time steps")
        ''''''
        def Status(model,g,m,t):
            if t == model.T[1]:
                #return pe.Constraint.Feasible
                return model.vOnoroff[g,m,t] == (model.pOnoroffinitial[m]+ model.vTurnon[g,m,t] - model.vTurnoff[g,m,t])
            return model.vOnoroff[g,m,t] == (model.vOnoroff[g,m,model.T.prev(t)] + model.vTurnon[g,m,t] - model.vTurnoff[g,m,t])
        model.Status = pe.Constraint(model.G,model.M, model.T, rule = Status, doc = "The status of the plant must be the previous state plus the turnon/off")
        
        def SingleGenStatus(model,g,m,t):
            return model.vOnoroff[g,m,t] <= model.vGenerator[g]
        model.SingleGenStatus = pe.Constraint(model.G,model.M, model.T, rule = SingleGenStatus, doc = "The status of the plant must be the previous state plus the turnon/off")
        '''
        def TurnOn(model,g,m,t):
            if t == model.T[1]:
                return model.vTurnon[g,m,t] == (model.vOnoroff[g,m,t]-model.pOnoroffinitial[g] + model.vTurnoff[g,m,t])
            return model.vTurnon[g,m,t] == (model.vOnoroff[g,m,t]-model.vOnoroff[g,m,model.T.prev(t)] + model.vTurnoff[g,m,t])
        model.TurnOn = pe.Constraint(model.G,model.M,model.T, rule = TurnOn, doc = "Ensure turn on variable goes from off to on states")
        
        def TurnOff(model,g,m,t):
            if t == model.T[1]:
                return model.vTurnoff[g,m,t] == (model.pOnoroffinitial[g] - model.vOnoroff[g,m,t] + model.vTurnon[g,m,t])
            return model.vTurnoff[g,m,t] == model.vOnoroff[g,m,model.T.prev(t)] - model.vOnoroff[g,m,t] + model.vTurnon[g,m,t]
        model.TurnOff = pe.Constraint(model.G,model.M,model.T, rule = TurnOff, doc = "Ensure turn off variable goes from on to off states")
        '''
        def MinDownTime(model,g,m,t):
            if t == model.T[1]:
                return pe.Constraint.Feasible
            if t <= model.pMDT[g]:
                return sum((1-model.vOnoroff[g,m,model.T.prev(t,td+1)]) for td in range(t)) >= t*model.vTurnon[g,m,t]
            return sum((1-model.vOnoroff[g,m,model.T.prev(t,td+1)]) for td in range(model.pMDT[g])) >= model.pMDT[g]*model.vTurnon[g,m,t]
        model.MinDownTime = pe.Constraint(model.G,model.M,model.T, rule = MinDownTime, doc = "We must enforce minimum down time")
        
        # V14 Updates
        
        def TES_SOC(model,g,t):
            if t == model.T[1]:
                if DAM == True:
                    return model.vTESSOC[g,t] == (model.pTESInitSOC[g]*(1-model.pTESLossRate[g]) 
                                      - model.pTES_Dis0[g] + model.pTES_Ch0[g])
                else:
                    return model.vTESSOC[g,t] == model.pTESInitSOC[g]
            return model.vTESSOC[g,t] == (model.vTESSOC[g,model.T.prev(t)]*(1-model.pTESLossRate[g]) 
                                          - model.vTESSOC_Discharge[g,model.T.prev(t)]
                                          + sum(model.vTESSOC_Charge[g,m,model.T.prev(t)] for m in model.M))
        model.TES_SOC = pe.Constraint(model.G,model.T, rule = TES_SOC, doc = "SOC is equal to the in/out plus the previous step with losses.")
        
        def TES_SOC_Max(model,g,t):
            return model.vTESSOC[g,t] <= model.pTESCapacity[g]
        model.TES_SOC_Max = pe.Constraint(model.G,model.T, rule = TES_SOC_Max, doc = "SOC has to be below the capacity limit")
        
        def TES_SOC_Power(model,g,t):
            return model.vTESSOC_Discharge[g,t] <= model.pTESPowerCap[g]
        model.TES_SOC_Power = pe.Constraint(model.G, model.T, rule = TES_SOC_Power, doc = "TES cannot surpass the power out max")
        
        def TES_SOC_Power2(model,g,t):
            return sum(model.vTESSOC_Charge[g,m,t]for m in model.M) <= model.pTESPowerCap[g]
        model.TES_SOC_Power2 = pe.Constraint(model.G, model.T, rule = TES_SOC_Power2, doc = "TES cannot surpass the power out max")
        
        self.model = model
        if not DAM:
            print('Model Built')
        
        
    def SolveModel(self, solver='cplex', DAM = False, init_run = False):
        self.BuildModel(DAM = DAM)
        if not DAM:
            print('Solving...')
        opt = pyomo.opt.SolverFactory(solver,tee = True) #iterlim = 1,
        opt.options.mipgap = 0.1
        

        results = opt.solve(self.model, tee = False, logfile='CPLEX_test_log_S'+str(args.Start)+'E'+str(args.End)+str(args.Region)+str(args.Year)+'3.log')
        # print('>>Solver status is {} and solver termination condition is {}'.format(results.solver.status,results.solver.termination_condition))
        
        if results.solver.termination_condition =='optimal':
            #print(pe.value(self.model.Obj_Profit))
            self.vTotalProfits = pe.value(self.model.Obj_Profit)
            self.rFeasible = True
            self.ResultData(DAM = DAM)
            self.ResultOutput(DAM = DAM)
            
        else:
            self.rFeasible = False
            self.vTotalProfits= None
            self.rACC = 0
            self.rFOMC = 0
            self.rVC = 0
            self.rFC = 0
            self.rEP = 0
            self.rAVC = 0
            self.rSUC = 0
            self.rTES = 0
            self.ModCount = 0
            
    def ResultData(self, DAM = False):
        self.vGen = np.array([[[pe.value(self.model.vGen[g,m,t]) for g in self.G] for m in self.M] for t in self.T])
        self.vTGen = np.array([[[pe.value(self.model.vTGen[g,m,t]) for g in self.G] for m in self.M] for t in self.T])
        self.vRGen = np.array([[[pe.value(self.model.vRGen[g,m,t]) for g in self.G] for m in self.M] for t in self.T])
        
        self.vTESSOC_Discharge = np.array([[pe.value(self.model.vTESSOC_Discharge[g,t]) for g in self.G] for t in self.T])
        self.vTESSOC_Charge = np.array([[[pe.value(self.model.vTESSOC_Charge[g,m,t]) for g in self.G] for m in self.M] for t in self.T])
        self.vTESSOC = np.array([[pe.value(self.model.vTESSOC[g,t]) for g in self.G] for t in self.T])
        
        # Binary variables
        self.vTurnon = np.array([[[pe.value(self.model.vTurnon[g,m,t]) for g in self.G] for m in self.M] for t in self.T])
        self.vTurnoff = np.array([[[pe.value(self.model.vTurnoff[g,m,t]) for g in self.G] for m in self.M] for t in self.T])
        self.vOnoroff = np.array([[[pe.value(self.model.vOnoroff[g,m,t]) for g in self.G] for m in self.M] for t in self.T])
        self.vEOnoroff = np.array([[[pe.value(self.model.vEOnoroff[g,m,t]) for g in self.G] for m in self.M] for t in self.T])
        #self.vTOnoroff = np.array([[[pe.value(self.model.vTOnoroff[g,m,t]) for g in self.G] for m in self.M] for t in self.T])
        self.vGenabovemin = np.array([[[pe.value(self.model.vGenabovemin[g,m,t]) for g in self.G] for m in self.M] for t in self.T])
        self.vRGenabovemin = np.array([[[pe.value(self.model.vRGenabovemin[g,m,t]) for g in self.G] for m in self.M] for t in self.T])
        
        self.vTotalProfits = pe.value(self.model.Obj_Profit)
        self.vGenerator = [pe.value(self.model.vGenerator[g]) for g in self.G]
        self.vModuleS = np.array([[pe.value(self.model.vModuleS[g,m]) for g in self.G] for m in self.M])
        self.vModDAM = [i[0] for i in self.vModuleS]
        self.ModCount = int(self.vModuleS.sum())
        
    def ResultOutput(self, DAM = False):
        self.vRGen_PD = pd.DataFrame(self.vRGen[:,:,self.vGenerator.index(1)])
        self.vGen_PD = pd.DataFrame(self.vGen[:,:,self.vGenerator.index(1)])
        self.vTGen_PD = pd.DataFrame(self.vTGen[:,:,self.vGenerator.index(1)])
        self.vTurnon_PD = pd.DataFrame(self.vTurnon[:,:,self.vGenerator.index(1)])
        self.vOnoroff_PD = pd.DataFrame(self.vOnoroff[:,:,self.vGenerator.index(1)])
        
        self.vTESSOC_Discharge = pd.DataFrame(self.vTESSOC_Discharge[:,self.vGenerator.index(1)])
        self.vTESSOC_Charge = pd.DataFrame(self.vTESSOC_Charge[:,:,self.vGenerator.index(1)])
        
        self.rHeatDelivered = self.vTGen_PD  - (self.vGen_PD/self.pThermaleff[self.vGenerator.index(1)])
        
        self.rERevenue = self.vGen_PD.multiply(self.pLMP, axis = 0)
        if not DAM:
            print(self.Generator_IDs[self.vGenerator.index(1)],"was chosen for operation")
            print(sum(sum(self.vModuleS)),"modules used")
        
        self.Output = pd.DataFrame()
        m = [i for i, x in enumerate(self.vModuleS[:,0].tolist()) if x == 1]
        for x in m:
            self.Output['Reactor_Gen [MWht]_'+str(x)] = self.vRGen_PD[x]
            self.Output['Thermal_Gen [MWht]_'+str(x)] = self.vTGen_PD[x]
            self.Output['Elec_Gen [MWhe]_'+str(x)] = self.vGen_PD[x]
            self.Output['Elec_Gen [MWht]_'+str(x)] = self.Output['Elec_Gen [MWhe]_'+str(x)]/float(self.pThermaleff[0])
            self.Output['TES_Charge_'+str(x)] = self.vTESSOC_Charge[x]
        self.Output['TES_Disharge'] = self.vTESSOC_Discharge[0]
        self.Output['SOC'] = self.vTESSOC[:,self.vGenerator.index(1)]
        self.Output
        
        self.rACC = float(sum(self.pCap)*sum(self.pCAPEX)*sum(self.vModuleS)*(self.pIR/(1-((1+self.pIR)**(-1*self.pLifetime))))*(len(self.T)/8760))
        self.rFOMC = float(sum(self.pCap)*sum(self.pFOPEX)*sum(self.vModuleS)*(len(self.T)/8760))
        self.rVC = float(self.vRGen_PD.values.sum()*self.pVOC.values.sum()*self.pThermaleff[self.vGenerator.index(1)])
        self.rFC = float(self.vRGen_PD.values.sum()*self.FuelCost.values.sum()*sum(self.pThermaleff))
        self.rEP = float(self.rERevenue.values.sum())
        self.rAVC = float((self.vTGen_PD*self.NGprice).values.sum())
        self.rSUC = float(sum(self.pStartupfixedcost)* sum(sum(sum(self.vTurnon))))
        self.rTES = float(sum(self.pTESCapacity)*self.pTESCapex*(self.pIR/(1-((1+self.pIR)**(-1*self.pLifetime))))*(len(self.T)/8760))
        
        
def ParamsVarsPD(base,sites,lmps, dem, TES, naics, EP, heatLMP = 4, minTemp = 500, MaxMod = 2, hrProfiles = None,startHr = 0, hrcount=8760, DAM = False, DAM_data = None, init_run = False):
    

    # Environmental paramters/System
    base.G = list(range(len(sites['Sites'])))
    base.T = list(range(hrcount))
    if init_run:
        base.T = list(range(len(base.T[::6])))
        
    base.M = list(range(MaxMod)) #will need to be fixed eventually!!!!!!!!!!!!!!!!!!!!!!
    base.Generator_IDs = sites['Sites'].tolist()

    base.NGprice = heatLMP
    
    lmpList = lmps['Settlement Point Price'].tolist()
    if init_run:
        base.pLMP = lmpList[::6]
    else:
        base.pLMP = lmpList[startHr:startHr+hrcount]
    if args.Zero == 0:
        base.pLMP = len(base.pLMP)*[0.000000001]

    base.pTempdemanded = minTemp

    if args.NoHeat:
        base.pTempdemanded = 0
        base.pTEnergyDemand = [0]*hrcount

    
    else:
        if naics== None:
            base.pTEnergyDemand = [dem]*len(base.pLMP)

    
        if naics != None:
            ###########################################################3333
            TotalHeat = dem*8760
            TotalHours = sum(hrProfiles[:8760])
            CorrectiveFactor = TotalHeat/TotalHours
            NormalizedHour = [hrProfiles[i] * CorrectiveFactor for i in range(len(hrProfiles))]
            base.pTEnergyDemand = NormalizedHour[startHr:startHr+hrcount]
            #print(CorrectiveFactor)
        if init_run:
            base.pTEnergyDemand = base.pTEnergyDemand[::6]
        #print(max(base.pTEnergyDemand))

    
    
    # Engineering parameters for generators
    base.pCAPEX = (sites['CAPEX $/kWe']*1000).tolist()   
    base.pFOPEX = (sites['FOPEX $/kWe']*1000).tolist()   
    base.pOutlettemp = sites['Outlet Temp (C)'].tolist()
    base.pWorkingTemp = sites['Working Fluid Temp (C)'].tolist()
        
    base.pStartupfixedcost = sites['Startupfixedcost in $'].tolist()
    
    
    
    base.pCap = sites['Power in MWe'].tolist() 
    base.pTCap = (sites['Power in MWt']).tolist() 
    base.pThermaleff = (sites['Thermal Efficiency']).tolist() 
    base.pThermaTransfereff = sites['Thermal Transfer Efficiency'].tolist() #V11
    base.pMSL = sites['MSL in MWe'].tolist()
    if args.Sens == 'CF95':
        base.pMSL = (sites['Power in MWe']*0.95).tolist() #V25
    base.pMSLTurb = sites['MSL_turb in MWe'].tolist()
    
        
    
    #base.pVOC = sites['Effective VOC in $/MWh-e'].tolist() # V18  
    #base.pTVOC = sites['Effective HVOC in $/MWh-t'].tolist() # V18
    
    base.pVOC = sites['VOM in $/MWh-e']
    base.FuelCost = sites['FC in $/MWh-e']
    
    base.pTVOC = ((base.pVOC+base.FuelCost)*base.pThermaleff).tolist()
    
    base.pMDT = sites['MDT in hours'].tolist()
    diff = [base.pCap[i] - base.pMSL[i] for i in range(len(base.pMSL))]
    base.pGenabovemininitial = [diff[i]/1.2 for i in range(len(diff))] # This value is currently not being used as it was too often limiting in the conditions. 
    base.pRamprate = sites['Ramp Rate (MW/hr)'].tolist()
    
    CapValue = TES
    base.pTESCapacity = [CapValue]*len(base.G)
    base.pTESChargeEff = [0.99]*len(base.G) # MIT report - >98% RTE
    base.pTESDischargeEff = [0.99]*len(base.G) # MIT report - >98% RTE
    base.pTESLossRate = [0.01]*len(base.G)
    base.pTESPowerCap = [CapValue/EP]*len(base.G)
    base.pTESCapex = 20000 # $/MWh-th generally from table 3 of https://aip.scitation.org/doi/pdf/10.1063/1.4984433 + MIT report For Nitrate Salt Storage

    base.pTESInitSOC = [0]*len(base.G)
    
    base.YR = None
    
    # V22
    if DAM:
        base.pMaxMod = [DAM_data[0]]*len(base.G)
        base.pMinMod = [DAM_data[0]]*len(base.G) 
        base.pTESInitSOC = [DAM_data[3]]*len(base.G)
        base.pRGenInitial = []
        for n in range(DAM_data[0]):
            d = DAM_data[1][n]-base.pMSL[0]
            if d > 0:
                base.pRGenInitial.append(d)
            else:
                base.pRGenInitial.append(0)
        base.pOnoroffinitial = DAM_data[2]
        
        base.pTES_Ch0 = [DAM_data[4]]*len(base.G)
        base.pTES_Dis0 = [DAM_data[5]]*len(base.G)
        
    else:
        base.pMaxMod = [MaxMod]*len(base.G) # THis can and likely will need to be changed, but for now it will work nicely.
        base.pMinMod = [0]*len(base.G) 
        base.pRGenInitial = [0]*MaxMod
        base.pOnoroffinitial = [0]*len(base.M) #sites['Onoroffinitial'].tolist()
        
def ImportData(sites , lmps , hrcount = 8760, startHr = 0):
    Sites_pd = pd.read_csv(sites)
    LMP_pd = pd.read_csv(lmps)
    
    LMP_pd = LMP_pd.iloc[startHr:startHr+hrcount]    
    return Sites_pd, LMP_pd

def NGTempCostCurve(Temp,NG_Cost = 3.0, AHF_Coeffs = [0,-0.00038,0.90556]):
    HHV = 54 # MJ/kg
    Density = 0.68 # kg/m3
    cfTom3 = 35.31 # Unit conversion
    AHF = AHF_Coeffs[0]*(Temp^2) + AHF_Coeffs[1]*(Temp) + AHF_Coeffs[2] # avaialble Heat fraction - Deep Patel Equation
    
    HHV = HHV*Density*(1/cfTom3)*(1/1000000)*(1000) # returns  TJ/thousand cf 
    Cost = NG_Cost*(1/HHV)*(1/AHF)*(1/277.778) # returns the Cost in $/MWh
    return Cost

def QuickGraphOutputs(Output,LMP,name,hr = 72,demand = False, TES = True):
    
    labels = Output.columns.tolist()
    labels_2 = []
    
    if TES:
        for i in labels:
            if 'Reactor' in i:
                pass
            elif 'SOC' in i:
                pass
            elif 'MWhe' in i:
                pass
            else:
                labels_2.append(i)
        OG = Output[labels_2]

        t = 0
        td = hr
        dashes = []
        colors = ['xkcd:orangered','xkcd:cyan','pink']
        pal = []
        l = len(OG.columns.tolist())
        d1 = 3
        d2 = 0
        while len(dashes)<l-1:
            dashes+=[(d1, d2),(d1, d2),(d1, d2)]
            pal+=colors
            d2+=2
        dashes+=[(1,2)]
        pal+=['xkcd:black']
        sns.set_style("whitegrid")
        ticks_24 = list(range(0,td+1,24))

        fig = plt.figure()
        plt.figure(figsize = (8,2))
        g = sns.lineplot(data=OG.loc[t:td,:], palette=pal, dashes=dashes)
        g.set_ylabel('Hourly Power [MWt]', fontsize=16)
        #g.lines[0].set_linestyle("--")
        plt.xticks(ticks_24)
        plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left', borderaxespad=0)
        g.figure.savefig(str(name)+'_Reactors.png', dpi = 1000, bbox_inches = "tight")

        fig2 = plt.figure()
        plt.figure(figsize = (8,2))
        g2 = sns.lineplot(data = LMP.loc[t:td+1,'Settlement Point Price'])
        g2.set_ylabel('LMP [$/MWh]', fontsize=16)
        g2.set_xlabel('Hour', fontsize=16)
        plt.xticks(ticks_24)
        g2.figure.savefig(str(name)+'_LMPS.png', dpi = 1000, bbox_inches = "tight")
    else:
        for i in labels:
            if 'Reactor' in i:
                pass
            elif 'SOC' in i:
                pass
            elif 'MWhe' in i:
                pass
            elif 'TES' in i:
                pass
            else:
                labels_2.append(i)
        OG = Output[labels_2]

        t = 0
        td = hr
        dashes = []
        colors = ['xkcd:orangered','xkcd:cyan']
        pal = []
        l = len(OG.columns.tolist())
        d1 = 3
        d2 = 0
        while len(dashes)<l-1:
            dashes+=[(d1, d2),(d1, d2)]
            pal+=colors
            d2+=2
        sns.set_style("whitegrid")
        ticks_24 = list(range(0,td+1,24))


        fig = plt.figure()
        plt.figure(figsize = (16,4))
        g = sns.lineplot(data=OG.loc[t:td,:], palette=pal, dashes=dashes)
        g.set_ylabel('Hourly Power [MWt]', fontsize=16)
        #g.lines[0].set_linestyle("--")
        plt.xticks(ticks_24)
        plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left', borderaxespad=0)
        g.figure.savefig(str(name)+'_Reactors.png')

        fig2 = plt.figure()
        plt.figure(figsize = (16,4))
        g2 = sns.lineplot(data=LMP.loc[t:td+1,'Settlement Point Price'])
        g2.set_ylabel('LMP [$/MWh]', fontsize=16)
        g2.set_xlabel('Hour', fontsize=16)
        plt.xticks(ticks_24)
        g2.figure.savefig(str(name)+'_LMPS.png')
    
def WaterfallPlot(Data,name, gtitle,hours = 8760):
    C = 8760/hours
    if isinstance(Data,list):
        y = [-1*Data[0]*C,-1*Data[1]*C,-1*Data[2]*C,-1*Data[3]*C,-1*Data[4]*C,Data[5]*C, ((Data[5])-((Data[0]+Data[1])+Data[2]+Data[3]+Data[4]))*C]
    else:
        y = [-1*A.rACC*C,-1*A.rFOMC*C,-1*A.rVC*C,-1*A.rFC*C,-1*A.rSUC*C,A.rEP*C, ((A.rEP)-((A.rACC+A.rFOMC)+A.rVC+A.rFC+A.rSUC))*C] # ,A.rAVC*C
    
    fig = go.Figure(go.Waterfall(
        orientation = "v",
        measure = ["relative", "relative", "relative", "relative", "relative", "relative", "total"],
        x = ["Annualized Capital Cost", "Fixed Cost", "Variable Cost", "Fuel Cost", "Start-Up Cost", "Electric Wholesale Revenues", "Net Revenue"],
        textposition = "outside",
        y = y,
        connector = {"line":{"color":"rgb(63, 63, 63)"}}
    )) #"TES Annualized Capital",  , "relative",    , "Natural Gas Cost"

    fig.update_layout(
            title = gtitle,
            showlegend = True,
            width=400, height=400,
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)'
    )

    fig.write_html(name+".html")
    fig.write_image(name+".png")


def main(Temp = 300, NGCost = 4, hrcount = 8760, dem = 100, lmps = 'lmps.csv',sites = 'sites_3.csv', genInd = None, TES = 0, MaxMod = 2, naics = None, year = 2019, EP = 10, startHr = 0, DAM = False, DAM_data = None, init_run = False):
    if isinstance(lmps,(str)):
        S,L= ImportData(lmps = lmps, sites = sites)
    else:
        S,L = sites,lmps
    Assessment = EconomicAssessment_OneGen()
    heatPrice = NGTempCostCurve(Temp,NG_Cost = NGCost)
    if genInd == None:
        pass
    else:
        S = S.loc[S.index == genInd]
    
    if naics == None:
        hrProfiles = [1]*8760
    if args.NoHeat:
        hrProfiles = [0]*8760
    else:
        hrProfiles = FP.QuickProfile(year, naics)
    s = datetime.now()
    # print(s,'Profiles complete')
    ParamsVarsPD(Assessment,S,L,dem, TES, naics, EP, heatLMP = heatPrice, minTemp = Temp, MaxMod =MaxMod, hrProfiles = hrProfiles, startHr = startHr, hrcount=hrcount, DAM = DAM, DAM_data = DAM_data, init_run = init_run)
    Assessment.SolveModel(DAM = DAM, init_run = init_run)
    return Assessment

def DAMmain(Temp = 300, NGCost = 4, dem = 100, lmps = 'SPP_SOUTH_2021.csv',sites = 'Base4_MSR.csv', genInd = None, TES = 0, InitMaxMod = 2, naics = 311221, year = 2021, EP = 10, startHr = 0, days = 1):
    # do the rough version
    S,L= ImportData(lmps = lmps, sites = sites)
    L2 = L.iloc[::12,:].reset_index()
    print('Running initial pass')
    hr = len(L2.index)
    A = main(Temp = Temp, NGCost = NGCost, hrcount = hr, dem = dem, lmps = L2,sites = S, genInd = genInd, TES = TES, MaxMod = InitMaxMod, naics = naics, year = year, EP = EP, startHr = startHr, DAM = False, DAM_data = None, init_run = True)
    
    MaxMod = A.ModCount
    if MaxMod ==0:
        MaxMod = 1
    if args.LimitedMod == True:    # V24 10/26/2022
        hrProfiles = FP.QuickProfile(year, naics)
        TotalHeat = dem*8760
        TotalHours = sum(hrProfiles[:8760])
        CorrectiveFactor = TotalHeat/TotalHours
        NormalizedHour = [hrProfiles[i] * CorrectiveFactor for i in range(len(hrProfiles))]
        peak = max(NormalizedHour)

        S = S.loc[S.index == genInd]
        cap = float(S['Power in MWt'])
        print(S)
        if int(math.ceil((peak/cap)*1.1)) >MaxMod:
            pass
        else:
            MaxMod = int(math.ceil((peak/cap)*1.1))
        print(MaxMod)
        
    InitGen = [0]*MaxMod
    InitStat = [0]*MaxMod
    InitSOC = 0
    Ch0 = 0
    Dis0 = 0
    Output = pd.DataFrame()
    hP = []
    hACC = []
    hFOMC = []
    hVC = []
    hFC = []
    hSUC = []
    hrTES = []
    hEP = []
    F = []
    if A.ModCount > 0:
        try: 
            for x in range(days):
                if x == days-1:
                    hrc = 24
                else:
                    hrc = 48
                print("____________________________{}___________________".format(x))
                A = main(Temp = Temp, NGCost = NGCost, hrcount = hrc, dem = dem, lmps = lmps,sites = sites, genInd = genInd, TES = TES, MaxMod = MaxMod, naics = naics, year = year, EP = EP, startHr = 24*x, DAM = True, DAM_data = [MaxMod,InitGen,InitStat, InitSOC, Ch0, Dis0], init_run = False)
                print('0')
                InitGen = A.vRGen_PD.iloc[23,:].tolist()
                InitStat = A.vOnoroff_PD.iloc[-1,:].tolist()
                Output = pd.concat((Output,A.Output.iloc[:24,:]))
                InitSOC = float(A.vTESSOC.sum(axis=1)[23])
                Ch0 = float(A.vTESSOC_Charge.sum(axis = 1)[23])
                Dis0 = float(A.vTESSOC_Discharge.sum(axis = 1)[23])
                print('a')


                vRGen = A.vRGen.sum(axis = 1)
                vGen = A.vGen.sum(axis = 1)
                TO = A.vTurnon.sum(axis = 1)[:24]
                EP_List = [A.pLMP[k]*A.vGen[k] for k in range(24)]
                print('b')
                rVC = float(sum(vRGen[:24])*A.pVOC*A.pThermaleff)
                rFC = float(sum(vRGen[:24])*A.FuelCost*A.pThermaleff)
                rSUC = float(sum(TO)*A.pStartupfixedcost)
                rEP = float(sum(sum(EP_List)))        
                print('c')
                hACC.append(A.rACC*(24/hrc))
                hFOMC.append(A.rFOMC*(24/hrc))
                hVC.append(rVC)
                hFC.append(rFC)
                hSUC.append(rSUC)
                hrTES.append(A.rTES*(24/hrc))
                hEP.append(rEP)
                hP.append(rEP-(A.rACC*(24/hrc)+A.rFOMC*(24/hrc)+A.rTES*(24/hrc)+rVC+rFC+rSUC))
                F.append(0)
                print('d')
        except:
            print('Failure - DAM Market')
            print(dem)
            F.append(2)
            hACC.append(0)
            hFOMC.append(0)
            hVC.append(0)
            hFC.append(0)
            hSUC.append(0)
            hrTES.append(0)
            hEP.append(0)
            hP.append(0)
    else:
        F.append(1)
        hACC.append(0)
        hFOMC.append(0)
        hVC.append(0)
        hFC.append(0)
        hSUC.append(0)
        hrTES.append(0)
        hEP.append(0)
        hP.append(0)
                
    df = {'Profit':hP,"Annualized Capital Cost":hACC, "Fixed Cost":hFOMC, "Variable Cost":hVC, "Fuel Cost":hFC, 
          "Start-Up Cost":hSUC, "TES Annualized Capital":hrTES, "Electric Wholesale Revenues":hEP, "Feasibility":F}


    Finances = pd.DataFrame(df)
    
    Output = Output.reset_index(drop= True)
    
    Fl = Finances.sum().tolist()
    WF = [Fl[1],Fl[2],Fl[3],Fl[4],Fl[5],Fl[6],Fl[7],Fl[8]] # This is missing the TES portion right now with 
    
    return Output, L, Finances, WF


def PerFacilityRuns(facilities,region,year,outDir, TESval = 0, MaxMod = 2, startHr = 0,days = 365):

    try: 
        os.mkdir(outDir) 
    except OSError as error: 
        print('Directory Already Exists, Overwriting')

    if region =='ERCOT':
        try: 
            LMPs = ['ERCOT_'+str(year)+'_HubAVG.csv']
        except:
            print("Do not have that vintage of LMP, can we check the back from something else?")
        if args.Year == 0:
            print('Will fail until new data is entered')
            LMPs = ['ERCOT_2018_HubAVG.csv','ERCOT_2019_HubAVG.csv','ERCOT_2020_HubAVG.csv','ERCOT_2021_HubAVG.csv']

    if region =='SPP':
        try: 
            LMPs = ['SPP_'+str(year)+'_HubAVG.csv']
        except:
            print("Do not have that vintage of LMP, can we check the back from something else?")
        if args.Year == 0:
            print('Will fail until new data is entered')
            LMPs = ['SPP_2018_HubAVG.csv','SPP_2019_HubAVG.csv','SPP_2020_HubAVG.csv','SPP_2021_HubAVG.csv']

    if args.Sens == 'MarketYear':
        LMPs = ['SPP_2020_HubAVG.csv','SPP_2021_HubAVG.csv','ERCOT_2018_HubAVG.csv','ERCOT_2019_HubAVG.csv','ERCOT_2020_HubAVG.csv','ERCOT_2021_HubAVG.csv']

    TES_caps = [0]
    if args.Sens == 'TES':
        TES_caps = [0,0.5,1,2]

    pdf2 = pd.DataFrame()
    pdfBestPer = pd.DataFrame()
    

    for TES_c in TES_caps:

        for LMP_Val in LMPs:
            T = facilities['Temp_degC'].tolist()
        
            facID = facilities['FACILITY_ID'].tolist()

            Ds = facilities['Thermal MWh/hr'].tolist()
        
            NAICS = facilities['FINAL_NAICS_CODE'].tolist()

            S,L = ImportData(lmps = LMP_Val, sites = 'Base5.csv')
            j = args.Start
            if args.NoHeat:
                Ds = [0]
                Ts = [0]
        
            while (j < len(T))& (j <= args.End):
                print(facID[j])
                print("####### ",str(j),"of ",str(len(T)),' ##########')
                G = []
                P = []
                LMP_List = []
                TES_list = []
                D_list = []
                T_list = []
                F_list = []
                Feas_list = []
                ACC = []
                FOMC = []
                VC = []
                FC = []
                SUC = []
                rTES = []
                EP = []
                AVC = []
                t = T[j]
                D = Ds[j]
                NAICS_j = NAICS[j]
                if args.Sens == 'Flat':
                    NAICS_j = 0
                N = []
                DEMTYPE = []
                GI = 0
                while GI < len(S.index):
                    print('Checkpoint 2')
                
                    D_list.append(D)
                    T_list.append(t)
                    F_list.append(facID[j])
                    A,L,F,WF = DAMmain(days = days, sites = 'Base5.csv', lmps = LMP_Val,dem = D, genInd = GI, Temp = t, TES = TES_c*D, InitMaxMod = MaxMod, startHr = 0, year = year, naics = NAICS_j)
                    LMP_List.append(LMP_Val)
                    TES_list.append(TES_c)
                    G.append(S.loc[S.index==GI]['Sites'].tolist()[0])
                    #print(WF)
                    TotalP = WF[6]-(WF[0]+WF[1]+WF[2]+WF[3]+WF[4]+WF[5])
                    P.append(TotalP)
                    if TotalP != None:
                        A.to_csv(outDir+'/'+str(facID[j])+'_'+str(t)+'_Op_'+str(S.loc[S.index==GI]['Sites'].tolist()[0])+'.csv') 
                        F.to_csv(outDir+'/'+str(facID[j])+'_'+str(t)+'_Fin_'+str(S.loc[S.index==GI]['Sites'].tolist()[0])+'.csv') 
                        N.append((len(A.columns.tolist())-2)/5)
                    else:
                        N.append(0)
                    Feas_list.append(WF[7])
                    ACC.append(-1*WF[0])
                    FOMC.append(-1*WF[1])
                    VC.append(-1*WF[2])
                    FC.append(-1*WF[3])
                    SUC.append(-1*WF[4])
                    rTES.append(-1*WF[5])
                    EP.append(WF[6])

                    DEMTYPE.append('hourly')
                
                
                    GI+=1 
                
                df = {'LMPType':LMP_List,'TES VAlue': TES_list,'Facility ID':F_list,'Temperature Req':T_list,'Generator':G,'Modules':N,
                      'Facility Demand':D_list,'Profit':P,"Annualized Capital Cost":ACC, 
                      "Fixed Cost":FOMC, "Variable Cost":VC, "Fuel Cost":FC, "Start-Up Cost":SUC,"TES Annualized Capital":rTES, 
                      "Electric Wholesale Revenues":EP,'Demand type':DEMTYPE, 'Feasibility':Feas_list} #  "Natural Gas Cost":AVC} 

                
                j+=1
                pdf = pd.DataFrame(df)
            
                prof_best = pdf['Profit'].max()
            
                pdfBestPer = pd.concat((pdfBestPer,pdf.loc[pdf['Profit']==prof_best]))
            
                pdf2 = pd.concat((pdf2,pdf))
                pdfBestPer.to_csv(outDir+'/BestPer.csv')   
                pdf2.to_csv(outDir+'/PerFacilityTest.csv')
                e = datetime.now()
                print(e,'end')
                d = e-s
                print(d)
        pdfBestPer.to_csv(outDir+'/BestPer.csv')   
        pdf2.to_csv(outDir+'/PerFacilityTest.csv')
    return pdfBestPer


def PerFacilityRunsNoDAM(facilities,region,year,outDir, TESval = 0, MaxMod = 2, startHr = 0,days = 365):

    try: 
        os.mkdir(outDir) 
    except OSError as error: 
        print('Directory Already Exists, Overwriting')

    if region =='ERCOT':
        try: 
            LMPs = ['ERCOT_'+str(year)+'_HubAVG.csv']
        except:
            print("Do not have that vintage of LMP, can we check the back from something else?")
        if args.Year == 0:
            print('Will fail until new data is entered')
            LMPs = ['ERCOT_2018_HubAVG.csv','ERCOT_2019_HubAVG.csv','ERCOT_2020_HubAVG.csv','ERCOT_2021_HubAVG.csv']

    if region =='SPP':
        try: 
            LMPs = ['SPP_'+str(year)+'_HubAVG.csv']
        except:
            print("Do not have that vintage of LMP, can we check the back from something else?")
        if args.Year == 0:
            print('Will fail until new data is entered')
            LMPs = ['SPP_2018_HubAVG.csv','SPP_2019_HubAVG.csv','SPP_2020_HubAVG.csv','SPP_2021_HubAVG.csv']

    pdf2 = pd.DataFrame()
    pdfBestPer = pd.DataFrame()
    for LMP_Val in LMPs:
        T = facilities['Temp_degC'].tolist()
        
        facID = facilities['FACILITY_ID'].tolist()

        Ds = facilities['Thermal MWh/hr'].tolist()
        
        NAICS = facilities['FINAL_NAICS_CODE'].tolist()

        S,L = ImportData(lmps = LMP_Val, sites = 'Base5.csv')
        j = args.Start
        if args.NoHeat:
            Ds = [0]
            Ts = [0]
        
        while (j < len(T))& (j <= args.End):
            print("####### ",str(j),"of ",str(len(T)),' ##########')
            G = []
            P = []
            LMP_List = []
            D_list = []
            T_list = []
            F_list = []
            Feas_list = []
            ACC = []
            FOMC = []
            VC = []
            FC = []
            SUC = []
            rTES = []
            EP = []
            AVC = []
            t = T[j]
            D = Ds[j]
            NAICS_j = NAICS[j]
            N = []
            DEMTYPE = []
            GI = 0
            h = 8
            while GI < len(S.index):
                print('Checkpoint 2')
                
                D_list.append(D)
                T_list.append(t)
                F_list.append(facID[j])
                A,L,F,WF = DAMmain(days = days, sites = 'Base5.csv', lmps = LMP_Val,dem = D, genInd = GI, Temp = t, TES = TESval, InitMaxMod = MaxMod, startHr = 0, year = year, naics = NAICS[j])

                print('Checkpoint 2')
                h = 8760
                A = main(Temp = t, hrcount = h, dem = D, lmps = LMP_Val,sites = 'Base5.csv', genInd = GI, TES = 0, MaxMod = 12, year = args.Year, EP = 10, startHr = 0, DAM = False, DAM_data = None, init_run = False)
                G.append(S.loc[S.index==GI]['Sites'].tolist()[0])

                vRGen = A.vRGen.sum(axis = 1)
                vGen = A.vGen.sum(axis = 1)
                TO = A.vTurnon.sum(axis = 1)
                EP_List = [A.pLMP[k]*A.vGen[k] for k in range(h)]

                rVC = float(sum(vRGen)*A.pVOC*A.pThermaleff)
                rFC = float(sum(vRGen)*A.FuelCost*A.pThermaleff)
                rSUC = float(sum(TO)*A.pStartupfixedcost)
                rEP = float(sum(sum(EP_List)))        
                LMP_List.append(LMP_Val)
    
                WF = [A.rACC,A.rFOMC,rVC,rFC,rSUC,rEP,0] # This is missing the TES portion right now with Fl[6]



                print(WF)
                TotalP = WF[5]-(WF[0]+WF[1]+WF[2]+WF[3]+WF[4])
                P.append(TotalP)
                if TotalP != None:
                    A.Output.to_csv(outDir+'/'+'_Op_'+str(S.loc[S.index==GI]['Sites'].tolist()[0])+'.csv') 
                    N.append((len(A.Output.columns.tolist())-2)/5)
                else:
                    N.append(0)
                Feas_list.append(WF[6])
                ACC.append(-1*WF[0])
                FOMC.append(-1*WF[1])
                VC.append(-1*WF[2])
                FC.append(-1*WF[3])
                SUC.append(-1*WF[4])
                #rTES.append(-1*WF[5])
                EP.append(WF[5])
                
            df = {'LMPType':LMP_List,'Facility ID':F_list,'Temperature Req':T_list,'Generator':G,'Modules':N,
                  'Facility Demand':D_list,'Profit':P,"Annualized Capital Cost":ACC, 
                  "Fixed Cost":FOMC, "Variable Cost":VC, "Fuel Cost":FC, "Start-Up Cost":SUC, 
                  "Electric Wholesale Revenues":EP,'Demand type':DEMTYPE, 'Feasibility':Feas_list} #  "Natural Gas Cost":AVC}"TES Annualized Capital":rTES, 

            print(df)
            j+=1
            pdf = pd.DataFrame(df)
            
            prof_best = pdf['Profit'].max()
            
            pdfBestPer = pd.concat((pdfBestPer,pdf.loc[pdf['Profit']==prof_best]))
            
            pdf2 = pd.concat((pdf2,pdf))
            pdfBestPer.to_csv(outDir+'/BestPer.csv')   
            pdf2.to_csv(outDir+'/PerFacilityTest.csv')
            e = datetime.now()
            print(e,'end')
            d = e-s
            print(d)
    pdfBestPer.to_csv(outDir+'/BestPer.csv')   
    pdf2.to_csv(outDir+'/PerFacilityTest.csv')
    return pdfBestPer


def PerMarketRunsEACH(region,year,outDir, TESval = 0, MaxMod = 12, startHr = 0,days = 365):

    try: 
        os.mkdir(outDir) 
    except OSError as error: 
        print('Directory Already Exists, Overwriting')

    if region =='ERCOT':
        try: 
            LMPs = ['ERCOT_'+str(year)+'_HubAVG.csv']
        except:
            print("Do not have that vintage of LMP, can we check the back from something else?")
        if args.Year == 0:
            print('Will fail until new data is entered')
            LMPs = ['ERCOT_2018_HubAVG.csv','ERCOT_2019_HubAVG.csv','ERCOT_2020_HubAVG.csv','ERCOT_2021_HubAVG.csv']

    if region =='SPP':
        try: 
            LMPs = ['SPP_'+str(year)+'_HubAVG.csv']
        except:
            print("Do not have that vintage of LMP, can we check the back from something else?")
        if args.Year == 0:
            print('Will fail until new data is entered')
            LMPs = ['SPP_2018_HubAVG.csv','SPP_2019_HubAVG.csv','SPP_2020_HubAVG.csv','SPP_2021_HubAVG.csv']

    pdf2 = pd.DataFrame()
    pdfBestPer = pd.DataFrame()
    for LMP_Val in LMPs:

        S,L = ImportData(lmps = LMP_Val, sites = 'Base5.csv')

        print("####### ",str(LMP_Val),' ##########')
        G = []
        P = []
        LMP_List = []
        D_list = []
        T_list = []
        F_list = []
        Feas_list = []
        ACC = []
        FOMC = []
        VC = []
        FC = []
        SUC = []
        rTES = []
        EP = []
        AVC = []
        N = []
        DEMTYPE = []
        GI = 0
        while GI < len(S.index):
            print('Checkpoint 2')
            h = 8760
            A = main(Temp = 0, hrcount = h, dem = 0, lmps = LMP_Val,sites = 'Base5.csv', genInd = GI, TES = 0, MaxMod = 12, year = 2020, EP = 10, startHr = 0, DAM = False, DAM_data = None, init_run = False)
            G.append(S.loc[S.index==GI]['Sites'].tolist()[0])

            vRGen = A.vRGen.sum(axis = 1)
            vGen = A.vGen.sum(axis = 1)
            TO = A.vTurnon.sum(axis = 1)
            EP_List = [A.pLMP[k]*A.vGen[k] for k in range(h)]

            rVC = float(sum(vRGen)*A.pVOC*A.pThermaleff)
            rFC = float(sum(vRGen)*A.FuelCost*A.pThermaleff)
            rSUC = float(sum(TO)*A.pStartupfixedcost)
            rEP = float(sum(sum(EP_List)))        
            LMP_List.append(LMP_Val)
    
            WF = [A.rACC,A.rFOMC,rVC,rFC,rSUC,rEP,0] # This is missing the TES portion right now with Fl[6]



            print(WF)
            TotalP = WF[5]-(WF[0]+WF[1]+WF[2]+WF[3]+WF[4])
            P.append(TotalP)
            if TotalP != None:
                A.Output.to_csv(outDir+'/'+'_Op_'+str(S.loc[S.index==GI]['Sites'].tolist()[0])+'.csv') 
                N.append((len(A.Output.columns.tolist())-2)/5)
            else:
                N.append(0)
            Feas_list.append(WF[6])
            ACC.append(-1*WF[0])
            FOMC.append(-1*WF[1])
            VC.append(-1*WF[2])
            FC.append(-1*WF[3])
            SUC.append(-1*WF[4])
            #rTES.append(-1*WF[5])
            EP.append(WF[5])

            DEMTYPE.append('hourly')
                
                
            GI+=1 
                
        df = {'LMPType':LMP_List,'Generator':G,'Profit':P,"Annualized Capital Cost":ACC, 
                "Fixed Cost":FOMC, "Variable Cost":VC, "Fuel Cost":FC, "Start-Up Cost":SUC, 
                "Electric Wholesale Revenues":EP,'Feasibility':Feas_list} 

        print(df)
        pdf = pd.DataFrame(df)
            
        pdf2 = pd.concat((pdf2,pdf))
        pdf2.to_csv(outDir+'/PerFacilityTest.csv')
        e = datetime.now()
        print(e,'end')
        d = e-s
        print(d)
    return pdfBestPer

if args.Region == 'ERCOT':
    Facilities = pd.read_csv('ERCOT_NG_Corr_01182024.csv')
else:
    Facilities = pd.read_csv('SPP_NG_Corr_01182024.csv')

if isinstance(args.Sens,str):
    Facilities = pd.read_csv('50_sensitivity.csv')

if args.NoDAM_EACH:
    PerMarketRunsEACH(args.Region,args.Year, args.Region+str(args.Year)+'_NOHeat', TESval = 0, MaxMod = 12, startHr = 0,days = 365)

elif args.NoDAM:
    PerFacilityRunsNoDAM(Facilities,args.Region,args.Year, args.Region+str(args.Year)+'_NoDAM', TESval = 0, MaxMod = 12, startHr = 0,days = 365)
    
elif isinstance(args.Sens,str):
    if args.Sens == 'NoDAM': # V25
        PerFacilityRunsNoDAM(Facilities,args.Region,args.Year, args.Region+str(args.Year)+'_NoDAM', TESval = 0, MaxMod = 12, startHr = 0,days = 365)
    else:
        PerFacilityRuns(Facilities,args.Region,args.Year,args.Region+str(args.Year)+'_NG_'+str(args.Sens), TESval = 0, MaxMod = 12, days = 365)
else:
    Results = PerFacilityRuns(Facilities,args.Region,args.Year,args.Region+str(args.Year)+'_NG_S'+str(args.Start)+'E'+str(args.End), TESval = 0, MaxMod = 12, days = 365)
