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

        ncols = 12
        try:
            self.dataset = pd.read_csv(infile, sep="\s+")
            self.logenergy_rec = self.dataset.logEnergyRec.tolist()
            self.coszenith_rec = self.dataset.cosZenithRec.tolist()
            self.ra_rec = self.dataset.AziRec.tolist()
            self.conv_weight = self.dataset.conv_weight.tolist()
            self.prompt_weight = self.dataset.prompt_weight.tolist()
            self.astro_weight = self.dataset.astro_weight.tolist()
            self.energy_prim = self.dataset.EnergyPrim.tolist()
            self.coszenith_prim = self.dataset.cosZenithPrim.tolist()
            self.ra_prim = self.dataset.AziPrim.tolist()
            self.ptype = self.dataset.TypePrim.tolist()
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

