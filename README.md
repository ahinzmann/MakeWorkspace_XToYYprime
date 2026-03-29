This project is used to make the XToYYprime analysis datacard and workspace, and it also has function to get the limit and fit plots. You should first use the pepper\_XToYYprime framework(https://github.com/TaozheYu/pepper\_XToYYprime) to produce the 3D histogram and then run this project.

You can follow these steps to procude the datacard and workspace:

# Setup the combine tool

We use the Combine v10 in CMSSW 14\_1\_X runs on el9

```bash
cmsrel CMSSW_14_1_0_pre4
cd CMSSW_14_1_0_pre4/src
cmsenv
git -c advice.detachedHead=false clone --depth 1 --branch v10.5.1 https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
scramv1 b clean; scramv1 b -j$(nproc --ignore=2) # always make a clean build, with n - 2 cores on the system
```
# Creat the local working area

Copy the project in the location

```bash
cd $CMSSW_BASE/src
git clone git@github.com:TaozheYu/MakeWorkspace_XToYYprime.git 
cd MakeWorkspace_XToYYprime
```
# Save the worksapces and datacards
```bash
python3 SaveXToYYprime_workspaces.py -y <year> -t <topology> -i <the input file path>
```

`-y` means the year of Run2, the choice is `2016`, `2016APV`, `2017` or `2018`

`-t` means the the topology of Yprime, the choice is `resolved` or `boosted`

`-i` means the 3D histogram of input path 

# Use the combine tool to get limit or do fit

```bash 
python3 run_combine.py -y <year> -t <topology> -rt <the run model>
```
`-rt` has two choices, `limit` means using asymptotic method to get limit; `fit` means diagnostics fit  

## Get the limit plots
```bash
python3 plotter_combineLimit.py -y <year> -t <topology> 
```

## Get the pre-fit and post-fit plots
```bash
python3 plot_prefit_and_postfit.py -y <year> -t <topology> -ft <fit type>
```
`-ft` has two choices, `prefit` means producing prefit plots, `postfit` means producing postfit plots




