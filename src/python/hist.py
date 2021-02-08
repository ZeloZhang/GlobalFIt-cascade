import numpy as np
import pandas as pd

class hist:
    def __init__(self, name, userbins_x, userbins_y, userbins_z):
        self.name = name
        self.userbins_x = userbins_x
        self.userbins_y = userbins_y
        self.userbins_z = userbins_z
        nbx = len(userbins_x)-1
        nby = len(userbins_y)-1
        nbz = len(userbins_z)-1
        self.value = np.zeros([nbx,nby,nbz])

    def Add(self, hist_component):
        self.value += hist_component.value

    def Reset(self):
        nbx = len(self.userbins_x)-1
        nby = len(self.userbins_y)-1
        nbz = len(self.userbins_z)-1
        self.value = np.zeros([nbx,nby,nbz])

    def Fill(self, logenergy_rec, coszenith_rec, ra_rec, weight=1):
        if weight is 1:
            weight = np.ones(len(logenergy_rec))
        self.value = np.histogramdd([logenergy_rec,coszenith_rec,ra_rec],bins=[self.userbins_x,self.userbins_y,self.userbins_z],weights=weight)[0]

    def Scale(self, norm):
        self.value *= norm

    def write_hist(self, outfile):
        index = pd.MultiIndex.from_product([range(s) for s in self.value.shape], names=["logenergy_rec","coszenith_rec","ra_rec"])
        df = pd.DataFrame({self.name:self.value.flatten()}, index=index)
        df.to_hdf(outfile,self.name)

    def set_hist(self, infile, key):
        nbx = len(self.userbins_x)-1
        nby = len(self.userbins_y)-1
        nbz = len(self.userbins_z)-1
        self.value = pd.read_hdf(infile,key).values.reshape(nbx,nby,nbz)

