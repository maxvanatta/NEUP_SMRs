# NEUP_SMRs

This repo is the code for work on "Deployment Pathways of Small Modular Reactors in Electric Power and Industrial Markets to Achieve Cost Reductions and Widespread Use" as part of a NEUP grant at the University of Michigan, SEAS, ASSET Lab under Dr. Michael Craig. 

The purpose of this research is to evaluate the economic feasibility of SMRs at providing process heat to industrial facilities while also gaining profits from participation in the electricity wholesale markets. More information will be available with the upcoming publication of this research. 

# Code Directions

There are two options for running this code: 1) using the Jupyter notebook included here or 2) the python script directly. 
For either, the version will be deonted with a V## filename. Early versions may still include bugs and limited functionality, therefore the most recent will be the most appropriate for use. 

**To run this code, first prepare your environment with:**

pip install pyomo
pip install cplex 
pip install docplex
pip install pandas
pip install seaborn
pip install plotly
pip install altair

- The reactors are included in the $base5.csv$ file. Documentation on sourcing these parameters will be included in the publication. 
- Energy markets included here are from ERCOT 2020 and SPP 2021. All years from 2018-2021 will be included soon.
- Facility demands are derived from NREL code here: https://data.nrel.gov/submissions/118
- This yearly demand is then dived into hourly heat loads using either an even distirbution of load or through distibution into healt load curves from the same link when available. 

**There are multiple options for running the code:**

**A)** run a single facility/SMR combination for the optimal yearly behavior with perfect foresight: $main() $

**B)** run a single facility/SMR combination for optimal behavior in a Day Ahead Market condition (daily 48 hour horizons): $DAMmain()$

**C)** run a series of all facilities within a market and all reactors using a DAM condition: $PerFacilityRuns()$

  PerFacilityRuns(facilities,region,year,outDir, TESval = 0, MaxMod = 2, startHr = 0,days = 365)
  - _facilites_ is the string of the csv file for facilities
  - _region_ is string for market (SPP or ERCOT)
  - _year_ is the int for the year tested
  - _outDir_ is string folder name for output files to be placed
  - _TESval_ is TES capacity in MWhth
  - _MasMod_ is initial maximum number of modules which can be built onsite
  - _startHr_ is int hour within the year to start optimization (keep at 0 is the days are 365)
  - _days_ is int of days to run model across

Ultimately this is still a work in progress. If you are using it, and if you cannot get it functioning or want more email: mvanatta@umich.edu
