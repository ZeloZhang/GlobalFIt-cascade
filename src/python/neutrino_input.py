import pandas as pd
import sys

from hist import hist

class neutrino_input:
    def __init__(self, name, bins_x, bins_y, bins_z):
        self.conv = hist(name+"_conv", bins_x, bins_y, bins_z) # histogram of observables. changes during fitting.
        self.prompt = hist(name+"_prompt", bins_x, bins_y, bins_z) # histogram of observables. changes during fitting.
        self.astro = hist(name+"_astro", bins_x, bins_y, bins_z) # histogram of observables. changes during fitting.

        self.binsx = bins_x
        self.nbinsx = len(bins_x)-1
        self.binsy = bins_y
        self.nbinsy = len(bins_y)-1
        self.binsz = bins_z
        self.nbinsz = len(bins_z)-1

    def read(self, infile):
        # file infile needs to be ASCII containing the following 11 columns
        # run_id, event_id, Enu, theta_nu, azimuth_nu, Erec, theta_rec, azimuth_rec, conv_weight, prompt_weight, astro_weight

        ncols = 18
        try:
            self.dataset = pd.read_csv(infile, sep="\s+")
            self.logenergy_rec = self.dataset.logEnergyRec
            self.coszenith_rec = self.dataset.cosZenithRec
            self.ra_rec = self.dataset.AziRec
            self.conv_weight = self.dataset.conv_weight
            self.prompt_weight = self.dataset.prompt_weight
            self.astro_weight = self.dataset.astro_weight
            self.energy_prim = self.dataset.EnergyPrim
            self.coszenith_prim = self.dataset.cosZenithPrim
            self.ra_prim = self.dataset.AziPrim
            self.ptype = self.dataset.TypePrim
            self.scattering = self.dataset.snowstorm_Parameters_Scattering
            self.absorption = self.dataset.snowstorm_Parameters_Absorption
            self.anisotropy_scale = self.dataset.snowstorm_Parameters_AnisotropyScale
            self.dom_efficiency = self.dataset.snowstorm_Parameters_DOMEfficiency
            self.holeiceforward_unified1 = self.dataset.snowstorm_Parameters_HoleIceForward_Unified1
            self.holeiceforward_unified2 = self.dataset.snowstorm_Parameters_HoleIceForward_Unified2
        except IOError:
            print("FATAL! unable to open file: ")
            print("... exiting")
            sys.exit(1)
                       
        if not len(self.dataset.keys())==ncols:
            print("FATAL! unexpected number of columns!")
            print("... exiting")
            sys.exit(1)
        self.size = len(self.dataset)

    def get_size(self):
        return self.size

    def get_nbinsx(self):
        return self.nbinsx

    def get_nbinsy(self):
        return self.nbinsy

    def get_nbinsz(self):
        return self.nbinsz

    def clear(self):
        self.dataset = 0

