#!/usr/bin/env python

from analysis import analysis
from stats import  stats
from helpers import par_options
from helpers import scan_options

def toinit():
    indir = "/data/user/zzhang1/fit_checkmanuel/input/"
    #indir = "/data/user/zzhang1/fit_reprohans_try/input/"
    # create analysis.
    wrapper = analysis()
    # need to specifically create the analysis.
    # this function contains all the analysis specific details and can be edited in ./src/python/analysis.py
    wrapper.create(indir)
    
    # get occasional prints of minimization progress from inside likelihood function
    
    # use stats class to analyze NuFit::analysis
    mini = stats(wrapper)
    
    astro_norm = par_options("astro_norm", 1.8, 0.01, 0., 10.)
    astro_index = par_options("astro_index", 2.3, 0.01, 0, 10.0)
    muon_norm = par_options("muon_norm", 1.5, 0.01, 0.0, 10.0)
    muon_norm1 = par_options("muon_norm_mlb", 1.0, 0.01, 0.0, 10.0)
    
    conv_norm = par_options("conv_norm", 0.9, 0.01, 0.0, 10.0)
    prompt_norm = par_options("prompt_norm", 0.2, 0.1, 0.0, 100.0)
    
    options = {
            astro_norm.name : astro_norm,
            astro_index.name : astro_index,
            conv_norm.name : conv_norm,
            muon_norm.name : muon_norm,
            muon_norm1.name : muon_norm1,
            prompt_norm.name: prompt_norm,
            }
    
    mini.set_options(options)
    return mini
#ana = mini.analysis["cascade_all"]
#mini.analysis[0].model.update_hists([1.47,1.05,0,1,1.59,2.51])
#print(mini.analysis[0].mcsum.value)
#print(mini.analysis[0].data.hist.value)
#mini.set_tolerance(1)
#mini.fit(False) # if true -> get profile LLH errors after minimization
