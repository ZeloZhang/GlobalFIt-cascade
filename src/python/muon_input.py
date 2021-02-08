import pandas as pd
import sys

from hist import hist

class muon_input:
    def __init__(self, name, bins_x, bins_y, bins_z):
        self.hist = hist(name,bins_x,bins_y,bins_z)

    def read(self,infile):
        # file infile needs to be ASCII containing the following 9 columns
        # run_id, event_id, Enu, theta_nu, azimuth_nu, Erec, theta_rec, azimuth_rec, muon_weight

        ncols = 9
        try:
            self.dataset = pd.read_csv(infile, sep="\s+")
            self.logenergy_rec = self.dataset.logEnergyRec
            self.coszenith_rec = self.dataset.cosZenithRec
            self.ra_rec = self.dataset.AziRec
            self.muon_weight = self.dataset.weight
            '''
            self.scattering = self.dataset.snowstorm_Parameters_Scattering
            self.absorption = self.dataset.snowstorm_Parameters_Absorption
            self.anisotropy_scale = self.dataset.snowstorm_Parameters_AnisotropyScale
            self.dom_efficiency = self.dataset.snowstorm_Parameters_DOMEfficiency
            self.holeiceforward_unified1 = self.dataset.snowstorm_Parameters_HoleIceForward_Unified1
            self.holeiceforward_unified2 = self.dataset.snowstorm_Parameters_HoleIceForward_Unified2
            '''
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


