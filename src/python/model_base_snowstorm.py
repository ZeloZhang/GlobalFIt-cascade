import numpy as np
import pandas as pd
import copy
import sys
import scipy
from scipy.stats import norm
from scipy.optimize import minimize
from scipy import integrate

class model_base_snowstorm:
    def __init__(self, fitdata, astro):
        self.input_hists = fitdata
        self.astro_model = astro
        par0 = "muon_norm"
        par1 = "conv_norm"
        par2 = "prompt_norm"
        par3 = "muon_norm_mlb"
        par4 = "scattering"
        par5 = "absorption"
        par6 = "anisotropy_scale"
        par7 = "dom_efficiency"
        par8 = "holeiceforward_unified1"
        par9 = "holeiceforward_unified2"
        self.par_names = [par0, par1, par2, par3, par4, par5, par6, par7, par8, par9]
        self.store_parameters()

        # default set for snowstorm reweighting
        self.sigmas = {
            "scattering": 0.04, 
            "absorption": 0.04,
            "anisotropy_scale": 0.4,
            "dom_efficiency": 0.04,
            "holeiceforward_unified1": 0.4,
            "holeiceforward_unified2": 0.08,
            }
        self.snowstorm_pars_range = {
            "scattering": (0.9,1.1), 
            "absorption": (0.9,1.1),
            "anisotropy_scale": (0,2),
            "dom_efficiency": (0.9,1.1),
            "holeiceforward_unified1": (-1,1),
            "holeiceforward_unified2": (-0.2,0.2),
            }
        print("... jointly analyze {} event selections.".format(len(fitdata)))

    def __del__(self):
        del self.astro_model
        print("... end of analysis. all models cleaned from memory.")

    def store_parameters(self):
        # store the ordering in a dictionary
        self.parameters = {}
        for i in range(len(self.par_names)):
            self.parameters[self.par_names[i]]=i

        self.npars_base = len(self.parameters)
        self.npars_astro = self.astro_model.get_npars()
        self.npars = self.npars_base + self.npars_astro

        # add names for astro parameters
        self.parameters.update(self.astro_model.get_par_names())
        self.par_names += list(self.astro_model.get_par_names().keys())

        #print(self.parameters)
        #print(self.par_names)
        # have to increase parameter indices for astro parameters by npars_base
        # to ensure that base parameters are always in the front of the containers
        for i in range(self.npars_base, self.npars):
            previous_indice = self.parameters[self.par_names[i]]
            self.parameters[self.par_names[i]]=previous_indice+self.npars_base

        # proceed with book keeping
        self.ndatasets = len(self.input_hists) 

        self.input_indices = {}
        for i in range(self.ndatasets):
            dataset = self.input_hists[i]
            self.input_indices[dataset.name] = i
        
        # and force caching of log factorials for subsequent computations
        self.cache_logfactorials()

        #print(self.parameters)

    def cache_logfactorials(self):
        current_val = 0
        self.log_factorials = {0:0}
        self.log_factorials.update({1:0})
        for i in range(2,101):
            current_val += np.log(i)
            self.log_factorials.update({i: current_val})

        self.stirling_constant = 0.5*np.log(2*np.pi)

    def log_factorial(self, value):
        if value <= 100:
            return self.log_factorials[value]
        else:
            return (value + 0.5) * np.log(value) - value + self.stirling_constant

    def likelihood(self, pars):
        # calculates a joint poisson-likelihood over all bins
        # ignores terms that are independent of the model parameters

        neglogl = 0

        # set all histogram to current parameters
        self.update_hists(pars)

        # need to loop over all eventselections and all bins
        for i in range(0, self.ndatasets):
            dataset = self.input_hists[i]

            observed = 0
            expected = 0
            expected = dataset.mcsum.value
            expected = [[list(map(lambda x:10**(-20) if x < 10**(-20) else x, d1)) for d1 in d2] for d2 in expected]
            #df_expected = pd.DataFrame(expected)
            #df_expected = df_expected.applymap(lambda x:10**(-20) if x < 10**(-20) else x)
            expected = np.array(expected)
            observed = np.array(dataset.data.hist.value)
            if ~expected.all():
                print("0 in expected")
                print(expected)
            llh_in_bin = (expected-observed*np.log(expected))
            neglogl += np.sum(llh_in_bin)
        if np.isnan(neglogl):
            print(" !!!!! FATAL ERROR !!!! - Likelihood evaluated to NaN. ")
            print(pars)
        print(pars)
        print(neglogl)
        return neglogl

    def likelihood_gof(self, pars):
        # calculates full likelihood (including data-only dependent terms)
        neglogl = self.likelihood(pars)
        neglogl = self.likelihood_gof_parameter_independent_part(neglogl) # adds the parameter independent terms
        return neglogl

    def likelihood_gof_parameter_independent_part(self, neglogl):
        # need to loop over all eventselections and all bins
        # adds terms to likelihood that are parameter independent
        # idential to the saturated gof (Baker, Cousins 1984)
        # note that delta-llh values cancel these terms again

        for i in range(0, self.ndatasets):
            dataset = self.input_hists[i]
            observed = 0
            observed = dataset.data.hist.value
            def get_likelihood_in_bin(observed):
                if observed:
                    llh_in_bin = (observed * np.log(observed) - observed)
                else:
                    llh_in_bin = 0
                return llh_in_bin
            llh_in_bin = [[get_likelihood_in_bin(i) for i in j] for j in observed] 
            neglogl += np.sum(llh_in_bin)

        return neglogl

    def likelihood_abs(self, pars):
        # calculates full likelihood (including data-only dependent terms)
        neglogl = self.likelihood(pars)
        neglogl = self.likelihood_abs_parameter_independent_part(neglogl) # adds the parameter independent terms
        return neglogl

    def likelihood_abs_parameter_independent_part(self, neglogl):
        # need to loop over all eventselections and all bins
        # adds terms to likelihood that are parameter independent

        for i in range(0, self.ndatasets):
            dataset = self.input_hists[i]
            observed = 0
            observed = dataset.data.hist.value
            llh_in_bin = [[self.log_factorial(int(i)) for i in j] for j in observed]
            neglogl += np.sum(llh_in_bin)

        return neglogl

    def get_par_names(self):
        return self.parameters

    def update_hists(self, pars):
        astro_pars = pars[self.npars_base:]
        self.update_astro(astro_pars,pars)
        self.update_atmospherics(pars)
        self.update_sum()

    def truncated_normal_distribution(self, xs, mu, sigma, variable_range):
        # input list of x avoid calculat integral multiple times
        def my_norm_distribution(x, loc, scale):
            return np.exp(-((x-loc)/scale)**2/2)/np.sqrt(2*np.pi)/scale

        f = scipy.stats.norm(loc=mu, scale=sigma)
        normalization = integrate.quad(f.pdf, variable_range[0], variable_range[1])
        return my_norm_distribution(xs,loc=mu,scale=sigma)/normalization[0] # assuming the original distribution of systematic parameters are uniform

    def get_snowstorm_reweight(self, pars, component):
        # to get snowstorm reweight factor
        # hard code order of parameters here, need to be able to fit with models automatically

        reweight_factor = 1
        # reweight_factor = truncated_normal_distribution()*(range_max-range_min)
        reweight_factor = reweight_factor * self.truncated_normal_distribution(xs=component.scattering, mu=pars[4], sigma=self.sigmas["scattering"], variable_range=self.snowstorm_pars_range["scattering"]) * (self.snowstorm_pars_range["scattering"][1]-self.snowstorm_pars_range["scattering"][0])
        reweight_factor = reweight_factor * self.truncated_normal_distribution(xs=component.absorption, mu=pars[5], sigma=self.sigmas["absorption"], variable_range=self.snowstorm_pars_range["absorption"]) * (self.snowstorm_pars_range["absorption"][1]-self.snowstorm_pars_range["absorption"][0])
        reweight_factor = reweight_factor * self.truncated_normal_distribution(xs=component.anisotropy_scale, mu=pars[6], sigma=self.sigmas["anisotropy_scale"], variable_range=self.snowstorm_pars_range["anisotropy_scale"]) * (self.snowstorm_pars_range["anisotropy_scale"][1]-self.snowstorm_pars_range["anisotropy_scale"][0])
        reweight_factor = reweight_factor * self.truncated_normal_distribution(xs=component.dom_efficiency, mu=pars[7], sigma=self.sigmas["dom_efficiency"], variable_range=self.snowstorm_pars_range["dom_efficiency"]) * (self.snowstorm_pars_range["dom_efficiency"][1]-self.snowstorm_pars_range["dom_efficiency"][0])
        reweight_factor = reweight_factor * self.truncated_normal_distribution(xs=component.holeiceforward_unified1, mu=pars[8], sigma=self.sigmas["holeiceforward_unified1"], variable_range=self.snowstorm_pars_range["holeiceforward_unified1"]) * (self.snowstorm_pars_range["holeiceforward_unified1"][1]-self.snowstorm_pars_range["holeiceforward_unified1"][0])
        reweight_factor = reweight_factor * self.truncated_normal_distribution(xs=component.holeiceforward_unified2, mu=pars[9], sigma=self.sigmas["holeiceforward_unified2"], variable_range=self.snowstorm_pars_range["holeiceforward_unified2"]) * (self.snowstorm_pars_range["holeiceforward_unified2"][1]-self.snowstorm_pars_range["holeiceforward_unified2"][0])
        return reweight_factor

    def update_astro(self, astro_pars, all_pars):
        # adjust histograms according to parameter values
        # this needs to be performed on each dataset that enters the likelihood

        for i in range(self.ndatasets):
            dataset = self.input_hists[i]
            # first reset all histograms
            dataset.nue.astro.Reset()
            dataset.numu.astro.Reset()
            dataset.nutau.astro.Reset()

            # adjust astro histograms
            flux = self.astro_model.get_flux(astro_pars, dataset.nue.energy_prim, dataset.nue.coszenith_prim, dataset.nue.ra_prim, dataset.nue.ptype)
            snowstorm_reweight = self.get_snowstorm_reweight(all_pars, dataset.nue)
            dataset.nue.astro.Fill(dataset.nue.logenergy_rec, dataset.nue.coszenith_rec, dataset.nue.ra_rec, [rw * aw * f for rw, aw,f in zip(snowstorm_reweight, dataset.nue.astro_weight, flux)])

            flux = self.astro_model.get_flux(astro_pars, dataset.numu.energy_prim, dataset.numu.coszenith_prim, dataset.numu.ra_prim, dataset.numu.ptype)
            snowstorm_reweight = self.get_snowstorm_reweight(all_pars, dataset.numu)
            dataset.numu.astro.Fill(dataset.numu.logenergy_rec, dataset.numu.coszenith_rec, dataset.numu.ra_rec, [rw * aw * f for rw, aw,f in zip (snowstorm_reweight, dataset.numu.astro_weight, flux)])

            flux = self.astro_model.get_flux(astro_pars, dataset.nutau.energy_prim, dataset.nutau.coszenith_prim, dataset.nutau.ra_prim, dataset.nutau.ptype)
            snowstorm_reweight = self.get_snowstorm_reweight(all_pars, dataset.nutau)
            dataset.nutau.astro.Fill(dataset.nutau.logenergy_rec, dataset.nutau.coszenith_rec, dataset.nutau.ra_rec, [rw * aw * f for rw,aw,f in zip(snowstorm_reweight, dataset.nutau.astro_weight, flux)])

    def update_sum(self):
        for i in range(self.ndatasets):
            dataset = self.input_hists[i]
            dataset.astro.Reset()
            dataset.astro.Add(dataset.nue.astro)
            dataset.astro.Add(dataset.numu.astro)
            dataset.astro.Add(dataset.nutau.astro)

            dataset.mcsum.Reset()
            dataset.mcsum.Add(dataset.astro)
            dataset.mcsum.Add(dataset.atm_conv)
            dataset.mcsum.Add(dataset.atm_prompt)
            dataset.mcsum.Add(dataset.muon.hist)

    def update_atmospherics(self,pars):
        for i in range(self.ndatasets):
            dataset = self.input_hists[i]
            if dataset.name=="cascade_mlb":
                dataset.muon.hist.Reset()
                dataset.muon.hist = copy.deepcopy(dataset.muon.hist_orig)
                dataset.muon.hist.Scale(pars[3]) # muon norm (no shape change)
            elif dataset.name=="cascade_as":
                dataset.muon.hist.Reset()
                dataset.muon.hist = copy.deepcopy(dataset.muon.hist_orig)
                dataset.muon.hist.Scale(pars[4]) # muon norm (no shape change)
            else:
                dataset.muon.hist.Reset()
                dataset.muon.hist = copy.deepcopy(dataset.muon.hist_orig)
                dataset.muon.hist.Scale(pars[0])

            dataset.atm_conv.Reset()
            dataset.atm_prompt.Reset()
            
            dataset.nue.conv.Reset()
            dataset.nue.prompt.Reset()
            dataset.numu.conv.Reset()
            dataset.numu.prompt.Reset()
            dataset.nutau.prompt.Reset()

            snowstorm_reweight = self.get_snowstorm_reweight(pars, dataset.nue)
            dataset.nue.conv.Fill(dataset.nue.logenergy_rec, dataset.nue.coszenith_rec, dataset.nue.ra_rec, [rw*cw for rw, cw in zip(snowstorm_reweight, dataset.nue.conv_weight)])
            dataset.nue.prompt.Fill(dataset.nue.logenergy_rec, dataset.nue.coszenith_rec, dataset.nue.ra_rec, [rw*pw for rw, pw in zip(snowstorm_reweight, dataset.nue.prompt_weight)])

            snowstorm_reweight = self.get_snowstorm_reweight(pars, dataset.numu)
            dataset.numu.conv.Fill(dataset.numu.logenergy_rec, dataset.numu.coszenith_rec, dataset.numu.ra_rec, [rw*cw for rw, cw in zip(snowstorm_reweight, dataset.numu.conv_weight)])
            dataset.numu.prompt.Fill(dataset.numu.logenergy_rec, dataset.numu.coszenith_rec, dataset.numu.ra_rec, [rw*pw for rw, pw in zip(snowstorm_reweight, dataset.numu.prompt_weight)])
            
            snowstorm_reweight = self.get_snowstorm_reweight(pars, dataset.nutau)
            dataset.nutau.prompt.Fill(dataset.nutau.logenergy_rec, dataset.nutau.coszenith_rec, dataset.nutau.ra_rec, [rw*pw for rw, pw in zip(snowstorm_reweight, dataset.nutau.prompt_weight)])

            dataset.atm_conv.Add(dataset.nue.conv)
            dataset.atm_conv.Add(dataset.numu.conv)
            dataset.atm_conv.Scale(pars[1]) # conv norm (no shape change)

            dataset.atm_prompt.Add(dataset.nue.prompt)
            dataset.atm_prompt.Add(dataset.numu.prompt)
            dataset.atm_prompt.Add(dataset.nutau.prompt)
            dataset.atm_prompt.Scale(pars[2]) # prompt norm (no shape change)

    def get_npars_base(self):
        return self.npars_base

    def get_npars(self):
        return self.npars

    def get_histograms(self, outfile, pars_map):
        # order parameters as required by model. check for parameter names
        pars = []
        pars = self.fill_parameters(pars_map)
        
        self.update_hists(pars)

        print("... writing analysis histograms")

        for single_hists in self.input_hists:
            single_hists.atm_conv.write_hist(outfile)
            single_hists.atm_prompt.write_hist(outfile)
            single_hists.astro.write_hist(outfile)
            single_hists.mcsum.write_hist(outfile)

            single_hists.nue.conv.write_hist(outfile)
            single_hists.nue.prompt.write_hist(outfile)
            single_hists.nue.astro.write_hist(outfile)
            single_hists.nutau.conv.write_hist(outfile)
            single_hists.nutau.prompt.write_hist(outfile)
            single_hists.nutau.astro.write_hist(outfile)
            single_hists.numu.conv.write_hist(outfile)
            single_hists.numu.prompt.write_hist(outfile)
            single_hists.numu.astro.write_hist(outfile)
            
            single_hists.muon.hist.write_hist(outfile)
            single_hists.data.hist.write_hist(outfile)

        print("... done")



    def get_hist_mcsum(self, pars_map):
        pars = []
        pars = self.fill_parameters(pars_map)
        self.update_hists(pars)
        hists_mcsum = []
        for i in range(self.ndatasets):
            dataset = self.input_hists[i]
            hists_mcsum.append(copy.deepcopy(dataset.mcsum))
        return hists_mcsum

    def fill_parameters(self, pars_map):
        # check if we have received the expected number of parameter values
        pars = np.zeros(self.npars)
        if not len(pars_map)==self.npars:
            print("!!! FATAL: received unexpected number of parameters: {} instead of: {}".format(len(pars_map), self.npars))
            print("... exiting")
            sys.exit(1)

        for it in self.parameters:
            # does parameter exist?
            if it not in pars_map:
                print("!!! FATAL: parameter {} not found. Did you specify it correctly?".format(it))
                print("!!! only found the following parameters:")
                print(pars_map)
                print("... exiting")
                sys.exit(1)
            # parameter found
            pars[self.parameters[it]]=pars_map[it]
        return pars

    def get_analysis_names(self):
        names = []
        for i in range(self.ndatasets):
            print(self.input_hists[i].name)
            names+=[self.input_hists[i].name]
        return names

    def cache_data_hists(self):
        self.data_hist_cache = [copy.deepcopy(item.data.hist) for item in self.input_hists]

    def restore_data_hists(self):
        for i in range(self.ndatasets):
            self.input_hists[i].data.hist = copy.deepcopy(self.data_hist_cache[i])
        self.data_hist_cache = []

    def set_hist(self, analysis, hist):
        index = self.input_indices[analysis]
        self.input_hists[index].data.hist = copy.deepcopy(hist)

    def gaussian_prior_penalty(self, x, mean, sigma):
        return 1/(2*sigma*sigma)*(x-mean)*(x-mean)

    def bivariate_prior_penalty(abso, scat, mean_abs, mean_scat):
        # need inverse covariance matrix. inverse is symmetric.
        # corresponding to sigma = 0.07, rho = -0.1

        cov_inf11 = 206.143
        cov_inf12 = 20.614

        delta_1 = abso - mean_abs
        delta_2 = scat - mean_scat

        return 0.5*(delta_1 * cov_inf11 *delta_1 + delta_2 * conv_inf11 * delta_2 + 2 * (delta_1 * cov_inf12 * delta_2))

    def update_auxillary_date(point):
        # use only during toy fits
        # make sure means are consistent with what parameter values have been simulated as toy
        return 













