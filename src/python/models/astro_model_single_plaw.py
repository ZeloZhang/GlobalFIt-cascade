import numpy as np

from astro_model_base import astro_model_base

class astro_model_single_plaw(astro_model_base):
    def __init__(self):
        # model name
        self.model_name = "single_powerlaw"

        # specify names and ordering (indices) of model parameters
        par0 = "astro_norm"
        par1 = "astro_index"
        self.par_names = [par0, par1]
        
        # store the ordering in a map
        self.parameters = {}
        for i in range(len(self.par_names)):
            self.parameters[self.par_names[i]]=i

        self.npars = len(self.parameters)

    def get_flux(self, pars, energy, coszen, ra, ptype):
        # assumes pars is ordered according to names provided above
        return [pars[0]*10**(-18)*np.power(en/(10**5), -1.0*pars[1]) for en in energy]
