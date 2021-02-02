import pandas as pd

from hist import hist

class data_input:
    def __init__(self, name, bins_x, bins_y, bins_z):
        self.hist = hist(name, bins_x, bins_y, bins_z)

    def read(self,infile):
        # file infile needs to be ASCII containing the following 9 columns
        # run_id, event_id, Erec, theta_rec, azimuth_rec

        ncols = 5
        self.mydf = pd.read_csv(infile, sep = "\s+")
        self.run = self.mydf.Run.tolist()
        self.event = self.mydf.Event.tolist()
        self.logenergy_rec = self.mydf.logEnergyRec.tolist()
        self.coszenith_rec = self.mydf.cosZenithRec.tolist()
        self.ra_rec = self.mydf.AziRec.tolist()
        self.size = len(self.mydf)

    def get_size(self):
        return self.size
