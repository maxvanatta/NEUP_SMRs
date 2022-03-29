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

import matplotlib.pyplot as plt
import seaborn as sns
import math

import pyomo.environ as pe
import pyomo.opt
from pyomo.opt import SolverStatus

from datetime import datetime
s = datetime.now()
print(s,'start')

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
        self.pTLMP = None                     # LMP heat $/MWht [g,t] or [t,yr]
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
        self.pRamprate = None                 # Maximum ramp (up and down) (MW/hr) [g]
        
        self.pCAPEX = None                    # Capital costs of the generator $/MWe [g]
        self.pFOPEX = None                    # Fixed O&M costs of the generator $/MWe [g]
        self.pModulesize = None               # Module size in MWe [g]
        self.pOutlettemp = None               # SMR outlet temperature [g]
        self.pWorkingTemp = None              # SMR working fluid temp 
        self.pMaxMod = None                   # maximum number of modules [g]
        self.FastMod = None                   # V12 Fast mdouel helper       
        
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
        
        
        
        
    def BuildModel(self):
        pTS = []
        for i in range(len(self.G)):
            if self.pOutlettemp[i] >= self.pTempdemanded:
                pTS.append(1)
            else:
                pTS.append(0)
            
        self.pThermalSat = pTS
        
        print('Building Model')
        model = pe.ConcreteModel()
        
        #### SETS ####     
        
        model.G = pe.Set(initialize = self.G, doc = "Generator IDs")
        model.T = pe.Set(initialize = self.T, doc = "Timestep")
        model.M = pe.Set(initialize = self.M, doc = "Modules")
        #model.YR = pe.Set(initialize = self.YR, doc = 'Years of model runs-- CURRENTLY NOT USED')
        
        
        #### PARAMETERS ####
        
        model.pLMP = pe.Param(model.T,initialize = {self.T[i]: self.pLMP[i] for i in range(len(self.T))} , doc = "LMP Electric $/MWhe")
        model.pTLMP = pe.Param(model.G, model.T,initialize = self.pTLMP, doc = "LMP Thermal $/MWht")
        model.pTEnergyDemand = pe.Param(model.T,initialize = {self.T[i]: self.pTEnergyDemand[i] for i in range(len(self.T))} , doc = "Energy demand MWht")
        
        model.pOnoroffinitial = pe.Param(model.G,initialize = {self.G[i]: self.pOnoroffinitial[i] for i in range(len(self.G))}, doc = "Inital state, 1 on, 0 off")
        model.pStartupfixedcost = pe.Param(model.G,initialize = {self.G[i]: self.pStartupfixedcost[i] for i in range(len(self.G))}, doc = "Startup costs, $500 default")
        model.pVOC = pe.Param(model.G,initialize = {self.G[i]: self.pVOC[i] for i in range(len(self.G))}, doc = "Variable O&M Electric $/MWhe")
        model.pTVOC = pe.Param(model.G,initialize = {self.G[i]: self.pTVOC[i] for i in range(len(self.G))}, doc = "Variable O&M Thermal $/MWhe")
        model.pCap = pe.Param(model.G,initialize = {self.G[i]: self.pCap[i] for i in range(len(self.G))} , doc = "Maximum Electric Capacity MWe")
        model.pTCap = pe.Param(model.G,initialize = {self.G[i]:self.pTCap[i] for i in range(len(self.G))}, doc = "Maximum Thermal Capacity MWt")
        model.pThermaleff = pe.Param(model.G,initialize = {self.G[i]: self.pThermaleff[i] for i in range(len(self.G))}, doc = "Thermal efficiency (MWe/MWt)")
        model.pThermalTransfereff = pe.Param(model.G,initialize = {self.G[i]: self.pThermaTransfereff[i] for i in range(len(self.G))}, doc = "Thermal transfer efficiency (MW -coolant/MWt- working fluid)")
        model.pMSL = pe.Param(model.G,initialize = {self.G[i]: self.pMSL[i] for i in range(len(self.G))}, doc = "Minimum Stable Load Thermal MWt")
        model.pMSLTurb = pe.Param(model.G,initialize = {self.G[i]: self.pMSLTurb[i] for i in range(len(self.G))}, doc = "Minimum Stable Load of Turbine MWe")
        model.pMDT = pe.Param(model.G,initialize = {self.G[i]: self.pMDT[i] for i in range(len(self.G))} , doc = "Mimimum Down time (hours)")
        model.pGenabovemininitial = pe.Param(model.G,initialize = {self.G[i]: self.pGenabovemininitial[i] for i in range(len(self.G))}, doc = "Initial Generation above the Minimun MWe")
        model.pRamprate = pe.Param(model.G,initialize = {self.G[i]: self.pRamprate[i] for i in range(len(self.G))}, doc = "Ramp rate up and down MW/hr")        
        model.pFOPEX = pe.Param(model.G,initialize = {self.G[i]: self.pFOPEX[i] for i in range(len(self.G))}, doc = "FOPEX pparamter $/MWe")
        model.pCAPEX = pe.Param(model.G,initialize = {self.G[i]: self.pCAPEX[i] for i in range(len(self.G))}, doc = "CAPEX values$/MWe")
        model.pWorkingTemp = pe.Param(model.G,initialize = {self.G[i]: self.pWorkingTemp[i] for i in range(len(self.G))}, doc = "working temperature which the reactor can supply in Degrees C")
        model.pThermalSat  = pe.Param(model.G,initialize = {self.G[i]: self.pThermalSat[i] for i in range(len(self.G))}, doc = "SMR can supply the temperature")
        model.pMaxMod = pe.Param(model.G, initialize = {self.G[i]: self.pMaxMod[i] for i in range(len(self.G))}, doc = "maximum number of modules")
        
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
        

        print('Params & Vars established')
        
        #### OBJECTIVES ####
        
        def ERevenues(model):
            return sum(sum(sum(model.vGen[g,m,t]*model.pLMP[t]for g in model.G) for m in model.M) for t in model.T)
        
        def TRevenues(model):
            return sum(sum((sum(model.vTGen[g,m,t] for m in model.M)+ model.vTESSOC_Discharge[g,t]*model.pTESDischargeEff[g])*model.pTLMP[g,t] for g in model.G) for t in model.T)
        
        def Costs(model):
            return (sum(model.pStartupfixedcost[g]*sum(sum(model.vTurnon[g,m,t] for t in model.T) for m in model.M) for g in model.G) 
                    + sum(sum(sum((model.vRGen[g,m,t]*model.pTVOC[g]) for g in model.G) for m in model.M) for t in model.T) 
                    + sum(sum(model.pCap[g]*(model.pCAPEX[g]*model.vModuleS[g,m]*(self.pIR/(1-((1+self.pIR)**(-1*self.pLifetime))))+model.pFOPEX[g]*model.vModuleS[g,m]) for m in model.M) for g in model.G)*(len(self.T)/8760)) # currently the switch to telectricty is cost-free
        #
        
        # V15 
        # Will need to input the cost portion of this function for the TES system. 
        
        
        def Obj_Profit(model):
            return ERevenues(model)+TRevenues(model)-Costs(model)
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
                    model.vTESSOC_Charge[g,m,t]/model.pTESChargeEff[g]) == model.vRGen[g,m,t]*model.pThermalTransfereff[g] # Updated V14 MV 03012022
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
                return pe.Constraint.Feasible # V11C THis means that we can start iur reactor at any output.
            return (model.vRGenabovemin[g,m,model.T.prev(t)])-(model.vRGenabovemin[g,m,t]) >= -1*model.pRamprate[g]/(model.pThermaleff[g])
        model.RRamprateUp = pe.Constraint(model.G,model.M,model.T, rule = RRamprateUp, doc = "The hourly change must not exceed the ramprate between time steps")
        
        def RRamprateDn(model,g,m,t):
            if t == model.T[1]:
                return pe.Constraint.Feasible # V11C THis means that we can start iur reactor at any output.
            return (model.vRGenabovemin[g,m,model.T.prev(t)]-model.vRGenabovemin[g,m,t])*(model.pThermaleff[g]) <= model.pRamprate[g]
        model.RRamprateDn = pe.Constraint(model.G,model.M,model.T, rule = RRamprateDn, doc = "The hourly change must not exceed the ramprate between time steps")
                
        def Status(model,g,m,t):
            if t == model.T[1]:
                return model.vOnoroff[g,m,t] == (model.pOnoroffinitial[g]+ model.vTurnon[g,m,t] - model.vTurnoff[g,m,t])
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
                return model.vTESSOC[g,t] == model.pTESCapacity[g]/2
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
        #model.pprint()
        print('Model Built')
        
        
    def SolveModel(self, solver='cplex'):
        self.BuildModel()
        print('Solving...')
        
        opt = pyomo.opt.SolverFactory(solver,tee = True)
        opt.options.mipgap = 0.01
        
        results = opt.solve(self.model, tee = True, logfile="CPLEX_test_log.log")
        
        print('>>Solver status is {} and solver termination condition is {}'.format(results.solver.status,results.solver.termination_condition))
        if results.solver.termination_condition =='optimal':
            print(pe.value(self.model.Obj_Profit))
            self.vTotalProfits = pe.value(self.model.Obj_Profit)
            self.rFeasible = True
            self.ResultData()
            #print(self.vTotalProfits)
            self.ResultOutput()
            
        else:
            self.rFeasible = False
            self.vTotalProfits= None        
            
    def ResultData(self):
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
        '''
        print(self.vGenerator)
        print(self.vModule)
        print(self.vTGen[:,:,self.vGenerator.index(1)])
        '''
        
    def ResultOutput(self):
        self.vRGen_PD = pd.DataFrame(self.vRGen[:,:,self.vGenerator.index(1)])
        self.vGen_PD = pd.DataFrame(self.vGen[:,:,self.vGenerator.index(1)])
        #self.vGen_PD.columns = self.Generator_IDs
        self.vTGen_PD = pd.DataFrame(self.vTGen[:,:,self.vGenerator.index(1)])
        #self.vTGen_PD.columns = self.Generator_IDs
        self.vTurnon_PD = pd.DataFrame(self.vTurnon[:,:,self.vGenerator.index(1)])
        #self.vTurnon_PD.columns = self.Generator_IDs
        
        self.vTESSOC_Discharge = pd.DataFrame(self.vTESSOC_Discharge[:,self.vGenerator.index(1)])
        self.vTESSOC_Charge = pd.DataFrame(self.vTESSOC_Charge[:,:,self.vGenerator.index(1)])
        
        self.rHeatDelivered = self.vTGen_PD  - (self.vGen_PD/self.pThermaleff[self.vGenerator.index(1)])
        
        self.rERevenue = self.vGen_PD.multiply(self.pLMP, axis = 0)
        self.rTRevenue = self.vTGen_PD*self.pTLMP_PD  
        
        #self.rStartupcosts = pd.DataFrame(self.vTurnon[:,:,self.vGenerator.index(1)]).multiply(self.pStartupfixedcost, axis = 1)
        #self.rStartupcosts.columns = self.Generator_IDs

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
        self.Output['SOC'] = self.vTESSOC
        self.Output

def ParamsVarsPD(base,sites,lmps, dems, TES, heatLMP = None, minTemp = 500, facility = 'Demanded Energy MWht', MaxMod = 2):
    # Environmental paramters/System
    base.G = list(range(len(sites['Sites'])))
    base.T = list(range(len(lmps['Settlement Point Price'])))
    base.M = list(range(MaxMod)) #will need to be fixed eventually!!!!!!!!!!!!!!!!!!!!!!
    base.Generator_IDs = sites['Sites'].tolist()

    base.FastMod = {(base.G[i],base.M[j]): np.ones([len(base.M),len(base.G)],dtype = int) for i in range(len(base.G)) for j in range(len(base.M))}
    
    if heatLMP == None:
        #lmps['Settlement Point Price'][:10]= -10
        
        HLMP = pd.DataFrame()
        for x in base.G:
            HLMP[x] = [sites['HLMP in $/MWh-t'][x]]*len(base.T)
        base.pTLMP_PD = HLMP
        HLMP = np.array(HLMP)
        HLMP_Dict ={(base.G[i],base.T[j]):HLMP[j,i] for i in range(len(base.G)) for j in range(len(base.T))}
        base.pTLMP = HLMP_Dict
    elif isinstance(heatLMP,(int,float)):
        HLMP = pd.DataFrame()
        for x in base.G:
            HLMP[x] = [heatLMP]*len(base.T)
        base.pTLMP_PD = HLMP
        HLMP = np.array(HLMP)
        HLMP_Dict ={(base.G[i],base.T[j]):HLMP[j,i] for i in range(len(base.G)) for j in range(len(base.T))}
        base.pTLMP = HLMP_Dict
    
    elif isinstance(heatLMP,(np.ndarray, np.generic)):
        if heatLMP.shape[0]==len(base.T):
            if heatLMP.shape[1]==len(base.G):
                HLMP = pd.DataFrame(heatLMP)
            else:
                print('FAILED DUE TO heatLMP SHAPE, axis 1')
        else:
            print('FAILED DUE TO heatLMP SHAPE, axis 0')
        base.pTLMP_PD = HLMP
        HLMP = np.array(HLMP)
        HLMP_Dict ={(base.G[i],base.T[j]):HLMP[j,i] for i in range(len(base.G)) for j in range(len(base.T))}
        base.pTLMP = HLMP_Dict
        
    elif isinstance(heatLMP,(list)):
        HLMP = pd.DataFrame(heatLMP)
        if len(heatLMP) == len(base.T):
            for x in base.G:
                HLMP[x] = heatLMP
        else:
            print('FAILED DUE TO heatLMP LIST LENGTH MISMATCH')
        base.pTLMP_PD = HLMP
        HLMP = np.array(HLMP)
        HLMP_Dict ={(base.G[i],base.T[j]):HLMP[j,i] for i in range(len(base.G)) for j in range(len(base.T))}
        base.pTLMP = HLMP_Dict
    else:
        print('heatLMP data types are float, int, np.darray, or list.  What did you put in?')
         
    base.pLMP = lmps['Settlement Point Price'].tolist()
    base.pTempdemanded = minTemp
    base.pTEnergyDemand = dems[facility].tolist()
    
    # Engineering parameters for generators
    base.pCAPEX = (sites['CAPEX $/kWe']*1000).tolist()   
    base.pFOPEX = (sites['FOPEX $/kWe']*1000).tolist()   
    base.pOutlettemp = sites['Outlet Temp (C)'].tolist()
    base.pWorkingTemp = sites['Working Fluid Temp (C)'].tolist()
    base.pOnoroffinitial = [0]*len(base.G) #sites['Onoroffinitial'].tolist()    
    base.pStartupfixedcost = sites['Startupfixedcost in $'].tolist()
    base.pVOC = sites['Effective VOC in $/MWh-e'].tolist()
    base.pTVOC = sites['Effective HVOC in $/MWh-t'].tolist()
    base.pCap = sites['Power in MWe'].tolist() 
    base.pTCap = (sites['Power in MWt']).tolist() 
    base.pThermaleff = (sites['Thermal Efficiency']).tolist() 
    base.pThermaTransfereff = sites['Thermal Transfer Efficiency'].tolist() #V11
    base.pMSL = sites['MSL in MWe'].tolist()
    base.pMSLTurb = sites['MSL_turb in MWe'].tolist()
    base.pMaxMod = [MaxMod]*len(base.G) # THis can and likely will need to be changed, but for now it will work nicely.
    
    base.pMDT = sites['MDT in hours'].tolist()
    diff = [base.pCap[i] - base.pMSL[i] for i in range(len(base.pMSL))]
    base.pGenabovemininitial = [diff[i]/1.2 for i in range(len(diff))] # This value is currently not being used as it was too often limiting in the conditions. 
    base.pRamprate = sites['Ramp Rate (MW/hr)'].tolist()
    
    CapValue = TES
    base.pTESCapacity = [CapValue]*len(base.G)
    base.pTESChargeEff = [0.85]*len(base.G)
    base.pTESDischargeEff = [0.85]*len(base.G)
    base.pTESLossRate = [0.01]*len(base.G)
    base.pTESPowerCap = [CapValue/10]*len(base.G)
    
    base.YR = None
        
def ImportData(sites = 'sites.csv', lmps = 'lmps.csv', dems = 'demands.csv', hrcount = 8760):
    Sites_pd = pd.read_csv(sites)
    LMP_pd = pd.read_csv(lmps)
    dem_pd = pd.read_csv(dems)
    
    LMP_pd = LMP_pd.iloc[:hrcount]    
    dem_pd = dem_pd.iloc[:hrcount]
    return Sites_pd, LMP_pd, dem_pd

def NGTempCostCurve(Temp,NG_Cost = 3.0, AHF_Coeffs = [0,-0.00038,0.90556]):
    HHV = 54 # MJ/kg
    Density = 0.68 # kg/m3
    cfTom3 = 35.31 # Unit conversion
    AHF = AHF_Coeffs[0]*(Temp^2) + AHF_Coeffs[1]*(Temp) + AHF_Coeffs[2] # avaialble Heat fraction - Deep Patel Equation
    
    HHV = HHV*Density*(1/cfTom3)*(1/1000000)*(1000) # returns  TJ/thousand cf 
    Cost = NG_Cost*(1/HHV)*(1/AHF)*(1/277.778) # returns the Cost in $/MWh
    return Cost

def QuickGraphOutputs(base,name,hr = 72):
    labels = A.Output.columns.tolist()
    labels_2 = []
    for i in labels:
        if 'Reactor' in i:
            pass
        elif 'SOC' in i:
            pass
        elif 'MWhe' in i:
            pass
        else:
            labels_2.append(i)
    OG = A.Output[labels_2]

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
    plt.figure(figsize = (16,4))
    g = sns.lineplot(data=OG.loc[t:td,:], palette=pal, dashes=dashes)
    g.set_ylabel('Hourly Power [MWt]', fontsize=16)
    #g.lines[0].set_linestyle("--")
    plt.xticks(ticks_24)
    plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left', borderaxespad=0)
    g.figure.savefig('NuScale_Powers_200'+str(name)+'.png')

    fig2 = plt.figure()
    plt.figure(figsize = (16,4))
    g2 = sns.lineplot(data=A.pLMP[t:td+1])
    g2.set_ylabel('LMP [$/MWh]', fontsize=16)
    g2.set_xlabel('Hour', fontsize=16)
    plt.xticks(ticks_24)
    g2.figure.savefig('NuScale_LMPs_200'+str(name)+'.png')


def main(Temp = 300, NGCost = 4, hrcount = 8760, dems = 'demands.csv', lmps = 'lmps.csv',sites = 'sites_3.csv', genInd = None, facility = 'Demanded Energy MWht', TES = 0, MaxMod = 2):
    if isinstance(dems,(str)):
        S,L,D = ImportData(dems = dems, lmps = lmps, sites = sites, hrcount= hrcount)
    else:
        S,L,D = sites,lmps, dems
    Assessment = EconomicAssessment_OneGen()
    heatPrice = NGTempCostCurve(Temp,NG_Cost = NGCost)

    if genInd == None:
        pass
    else:
        S = S.loc[S.index == genInd]
    ParamsVarsPD(Assessment,S,L,D, TES, heatLMP = heatPrice, minTemp = Temp, facility = facility, MaxMod =MaxMod)
    Assessment.SolveModel( solver='cplex')
    return Assessment

A = main(hrcount = int(8760),sites = 'sites_nuscale.csv',lmps = 'lmps_TX_camb.csv',NGCost = 4, MaxMod=2)
e = datetime.now()
print(e,'end')
d = e-s
print(d)
