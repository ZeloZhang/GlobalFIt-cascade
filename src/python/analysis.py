#!/usr/bin/env python

import numpy as np

from hists import hists
from models.astro_model_single_plaw import astro_model_single_plaw
from model_base_snowstorm import model_base_snowstorm

class analysis:
    def __init__(self):
        self.verbose_llh = False
        self.n_llh_evals = 0
    
    def create(self,defdir):
        print("... start wrapping analyses")

        # here goes the analysis code (parsing of data and mc)

        # define cascade signal sample
        binsx = np.linspace(2.6, 7, 23) # logE
        binsy = [-2, 0.2, 0.6, 2] # cosz (3 bins: Northern Sky, middle, Southern Sky)
        binsz = [-10, 10] # ra (1 bin)

        dir_myc = defdir+"myc/"
        dir_mlb = defdir+"mlb/"
        dir_as = defdir+"as/"

        dir_myc_HKKMS06 = dir_myc + "HI_30cm/"
        dir_mlb_HKKMS06 = dir_mlb + "HI_30cm/"

        analysis1 = hists("cascade_all", binsx, binsy, binsz)
        name_nue = "neutrino/nue_cscd_all_neutrino.txt"
        name_numu = "neutrino/numu_cscd_all_neutrino.txt"
        name_nutau = "neutrino/nutau_cscd_all_neutrino.txt"
        name_mu = "HI_50cm/neutrino/mgun_cscd_all_neutrino.txt"
        name_data = "data/txt/data_cascade_dustcorr.txt"
        analysis1.read(dir_myc_HKKMS06+name_nue, dir_myc_HKKMS06+name_numu, dir_myc_HKKMS06+name_nutau, dir_myc+name_mu, dir_myc+name_data)

        # define muon background sample
        binsx2 = [2.6, 4.778] # logE
        binsy2 = [-2, 2] # cosz (1 bin: all-sky)
        binsz2 = [-10, 10] # ra (1 bin)

        analysis2 = hists("muon", binsx2, binsy2, binsz2)
        name_nue2 = "muon/muon_nue.txt"
        name_numu2 = "muon/muon_numu.txt"
        name_nutau2 = "muon/muon_nutau.txt"
        name_mu2 = "HI_50cm/muon/muon_mgun.txt"
        name_data2 = "data/txt/data_muon.txt"
        analysis2.read(dir_myc_HKKMS06+name_nue2, dir_myc_HKKMS06+name_numu2, dir_myc_HKKMS06+name_nutau2, dir_myc+name_mu2, dir_myc+name_data2)

        # define numu control sample
        binsx3 = np.linspace(2.6, 4.8, 12)
        analysis3 = hists("hybrid", binsx3, binsy2, binsz)
        name_nue3 = "neutrino/nue_hybrid_neutrino.txt"
        name_numu3 = "neutrino/numu_hybrid_neutrino.txt"
        name_nutau3 = "neutrino/nutau_hybrid_neutrino.txt"
        name_mu3 = "HI_50cm/neutrino/mgun_hybrid_neutrino.txt"
        name_data3 = "data/txt/data_hybrid.txt"
        analysis3.read(dir_myc_HKKMS06+name_nue3, dir_myc_HKKMS06+name_numu3, dir_myc_HKKMS06+name_nutau3, dir_myc+name_mu3, dir_myc+name_data3)

        # define mlb sample
        binsx4 = np.linspace(4, 7, 16)
        analysis4 = hists("cascade_mlb", binsx4, binsy, binsz)
        name_nue4 = "txt/nue_mlb_neutrino.txt"
        name_numu4 = "txt/numu_mlb_neutrino.txt"
        name_nutau4 = "txt/nutau_mlb_neutrino.txt"
        name_mu4 = "HI_50cm/txt/cors_mlb.txt"
        name_data4 = "data/txt/data_mlb_dustcorr.txt"
        #analysis4.read(dir_mlb_HKKMS06+name_nue4, dir_mlb_HKKMS06+name_numu4, dir_mlb_HKKMS06+name_nutau4, dir_mlb+name_mu4, dir_mlb+name_data4);

        # define as sample
        binsx_as = [4.53 + i*0.352857 for i in range(8)]
        binsy_as = [-2,2]
        analysis5 = hists("cascade_as", binsx_as, binsy_as, binsz)
        name_nue5 = "as_nue.txt"
        name_numu5 = "as_numu.txt"
        name_nutau5 = "as_nutau.txt"
        name_mu5 = "as_cors.txt"
        name_data5 = "as_data_SPE.txt"
        #analysis5.read(dir_as+name_nue5, dir_as+name_numu5, dir_as+name_nutau5, dir_as+name_mu5, dir_as+name_data5)

        #analyses = [analysis1]
        analyses = [analysis1,analysis2, analysis3]
        #analyses = [analysis1, analysis2, analysis3, analysis4]
        #analyses = [analysis1, analysis2, analysis3, analysis4, analysis5]

        # create astro model. then create base model
        astro = astro_model_single_plaw()
        mymodel = model_base_snowstorm(analyses, astro)
        self.model = mymodel

        print("... wrapping done")

    def change_astro_model(self, astro):
        self.model.change_astro_model(astro)

    def get_likelihood(self, pars):
        neglogl = self.model.likelihood(pars)
        self.n_llh_evals += 1
        # factor of 2 for Wilk's theorem
        return 2 * neglogl

    def get_likelihood_gof(self, pars):
        # factor of 2 for Wilk's theorem
        return 2 * self.model.likelihood_gof(pars)

    def get_likelihood_abs(self, pars):
        # factor of 2 for Wilk's theorem
        return 2 * self.model.likelihood_abs(pars)

    def get_n_llh_evals(self):
        return self.n_llh_evals

    def get_histograms(self, outfile, pars):
        self.model.get_histograms(outfile, pars)

    def get_hist_mcsum(self, pars):
        self.model.get_hist_mcsum(pars)

    def get_llh_evals(self):
        return self.n_llh_evals

    def get_npars(self):
        return self.model.get_npars()

    def get_par_names(self):
        return self.model.get_par_names()
       
    def get_histograms(self, outfile, **pars):
        self.model.get_histograms(outfile, **pars)

    def get_hist_mcsum(self, **pars):
        return self.model.get_hist_mcsum(**pars)

    def reset_n_llh_evals(self):
        self.n_llh_evals = 0

    def get_analysis_names(self):
        return self.model.get_analysis_names()

    def set_hist(self, analysis, hist):
        self.model.set_hist(analysis, hist)
        
    def cache_data_hists(self):
        self.model.cache_data_hists()

    def restore_data_hists(self):
        self.model.restore_data_hists()

    def update_auxillary_data(self,**pars):
        self.model.update_auxillary_data(**pars)

    def reset_auxillary_data(self):
        self.model.reset_auxillary_data()
