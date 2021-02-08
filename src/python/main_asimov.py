#!/usr/bin/env python

from analysis import analysis
from stats import  stats
from helpers import par_options
from helpers import scan_options
# snowstorm indir
indir = "/data/user/zzhang1/GlobalFit-cascade/fit_input_snowstorm/"
#indir = "/data/user/zzhang1/fit_checkmanuel/input_manuelsim/"
#indir = "/data/user/zzhang1/fit_checkmanuel/input/"
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
prompt_norm = par_options("prompt_norm", 0.2, 0.1, 0.0, 0.0)
#prompt_norm = par_options("prompt_norm", 0, 0.1, 0.0, 0.0)
scattering = par_options("scattering", 1, 0.01, 0.9, 1.1)
absorption = par_options("absorption", 1, 0.01, 0.9, 1.1)
anisotropy_scale = par_options("anisotropy_scale", 1, 0.1, 0, 2)
dom_efficiency = par_options("dom_efficiency", 1, 0.01, 0.9, 1.1)
holeiceforward_unified1 = par_options("holeiceforward_unified1", 0, 0.1, -1, 1)
holeiceforward_unified2 = par_options("holeiceforward_unified2", 0, 0.01, -0.2, 0.2)
'''
# to fix sys parameters
scattering = par_options("scattering", 1, 0.01, 1, 1)
absorption = par_options("absorption", 1, 0.01, 1, 1)
anisotropy_scale = par_options("anisotropy_scale", 1, 0.1, 1, 1)
dom_efficiency = par_options("dom_efficiency", 1, 0.01, 1, 1)
holeiceforward_unified1 = par_options("holeiceforward_unified1", 0, 0.1, 0, 0)
holeiceforward_unified2 = par_options("holeiceforward_unified2", 0, 0.01, 0, 0)
'''

options = {
        astro_norm.name : astro_norm,
        astro_index.name : astro_index,
        conv_norm.name : conv_norm,
        muon_norm.name : muon_norm,
        muon_norm1.name : muon_norm1,
        prompt_norm.name: prompt_norm,
        scattering.name : scattering,
        absorption.name : absorption,
        anisotropy_scale.name : anisotropy_scale,
        dom_efficiency.name : dom_efficiency,
        holeiceforward_unified1.name : holeiceforward_unified1,
        holeiceforward_unified2.name : holeiceforward_unified2
        }

mini.set_options(options)
#mini.set_tolerance(1)
mini.set_tolerance(0.0000001)

inject_pars = {
        astro_norm.name: 1.66,
        astro_index.name: 2.53,
        conv_norm.name: 1.07,
        muon_norm.name: 1.45,
        muon_norm1.name: 1.0,
        prompt_norm.name: 0.0,
        scattering.name: 1.02,
        absorption.name: 1.03,
        anisotropy_scale.name: 1.2,
        dom_efficiency.name: 0.94,
        holeiceforward_unified1.name: 0.5,
        holeiceforward_unified2.name: 0.1
        }

mini.asimov_fit("injected_hist233.h5", inject_pars) # if true -> get profile LLH errors after minimization
#mini.inject_data("injected_hist233.h5", inject_pars) # if true -> get profile LLH errors after minimization
