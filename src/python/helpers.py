class par_options:
    def __init__(self, name_="name", seed_=0, stepsize_=0.01, limit_low_=-10**5, limit_high_=10**5):
        self.name = name_
        self.seed = seed_
        self.stepsize = stepsize_
        self.limit_low = limit_low_
        self.limit_high = limit_high_

    def set_options(self, seed_, stepsize_, limit_low_, limit_high_):
        self.name = name_
        self.seed = seed_
        self.stepsize = stepsize_
        self.limit_low = limit_low_
        self.limit_high = limit_high_

class scan_options:
    def __init__(self, name_="name", nsteps_=10, range_low_=-10**5, range_high_=10**5):
        self.name = name_
        self.nsteps = nsteps_
        self.range_low = range_low_
        self.range_high = range_high_

    def set_options(self, name_, nsteps_, range_low_, range_high_):
        self.name = name_
        self.nsteps = nsteps_
        self.range_low = range_low_
        self.range_high = range_high_

