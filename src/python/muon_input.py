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
            self.logenergy_rec = self.dataset.logEnergyRec.tolist()
            self.coszenith_rec = self.dataset.cosZenithRec.tolist()
            self.ra_rec = self.dataset.AziRec.tolist()
            self.muon_weight = self.dataset.weight.tolist()
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


