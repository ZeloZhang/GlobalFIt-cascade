import copy

from hist import hist
from neutrino_input import neutrino_input
from muon_input import muon_input
from data_input import data_input

class hists:
    def __init__(self, name, userbins_x, userbins_y, userbins_z):
        self.atm_conv = hist(name+"_all_conv", userbins_x, userbins_y, userbins_z)
        self.atm_prompt = hist(name+"_all_prompt", userbins_x, userbins_y, userbins_z)
        self.astro = hist(name+"_all_astro", userbins_x, userbins_y, userbins_z)
        self.mcsum = hist(name+"_mcsum", userbins_x, userbins_y, userbins_z)

        self.nue = neutrino_input(name+"_nue", userbins_x, userbins_y, userbins_z)
        self.numu = neutrino_input(name+"_numu", userbins_x, userbins_y, userbins_z)
        self.nutau = neutrino_input(name+"_nutau", userbins_x, userbins_y, userbins_z)
        self.muon = muon_input(name+"_muon", userbins_x, userbins_y, userbins_z)
        self.data = data_input(name+"_data", userbins_x, userbins_y, userbins_z)

        # keep a copy of input bins
        self.binsx = userbins_x
        self.nbinsx = len(userbins_x)-1
        self.binsy = userbins_y
        self.nbinsy = len(userbins_y)-1
        self.binsz = userbins_z
        self.nbinsz = len(userbins_z)-1

        self.name = name

    def read(self, f_nue, f_numu, f_nutau, f_muon, f_data):
        # read input from all neutrinos
        # and initialize data histogram

        print("... preparing event selection: ", self.name)
        print("... reading nue")
        self.nue.read(f_nue)
        print("... done")
        print("... reading numu")
        self.numu.read(f_numu)
        print("... done")
        print("... reading nutau")
        self.nutau.read(f_nutau)
        print("... done")
        print("... reading muon")
        self.muon.read(f_muon)
        print("... done")
        print("... reading data")
        self.data.read(f_data)
        print("... done")

        # initialize data histograms
        print("... filling histograms")

        # nue
        self.nue.conv.Fill(self.nue.logenergy_rec, self.nue.coszenith_rec, self.nue.ra_rec, self.nue.conv_weight)
        self.nue.prompt.Fill(self.nue.logenergy_rec, self.nue.coszenith_rec, self.nue.ra_rec, self.nue.prompt_weight)
        # numu
        self.numu.conv.Fill(self.numu.logenergy_rec, self.numu.coszenith_rec, self.numu.ra_rec, self.numu.conv_weight)
        self.numu.prompt.Fill(self.numu.logenergy_rec, self.numu.coszenith_rec, self.numu.ra_rec, self.numu.prompt_weight)
        # nutau
        self.nutau.prompt.Fill(self.nutau.logenergy_rec, self.nutau.coszenith_rec, self.nutau.ra_rec, self.nutau.prompt_weight)

        # add nue, numu, nutau histograms
        self.atm_conv.Add(self.nue.conv)
        self.atm_conv.Add(self.numu.conv)
        # copy the histogram
        self.atm_conv_orig = copy.deepcopy(self.atm_conv)

        self.atm_prompt.Add(self.nue.prompt)
        self.atm_prompt.Add(self.numu.prompt)
        self.atm_prompt.Add(self.nutau.prompt)
        # copy the histogram
        self.atm_prompt_orig = copy.deepcopy(self.atm_prompt)

        # muon
        self.muon.hist.Fill(self.muon.logenergy_rec, self.muon.coszenith_rec, self.muon.ra_rec, self.muon.muon_weight)
        # copy the histogram
        self.muon.hist_orig = copy.deepcopy(self.muon.hist)

        # data
        self.data.hist.Fill(self.data.logenergy_rec, self.data.coszenith_rec, self.data.ra_rec)

        print("... done")


