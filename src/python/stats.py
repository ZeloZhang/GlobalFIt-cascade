import copy
import sys
import time
from scipy.stats import norm
from scipy.optimize import minimize
from scipy import integrate

class stats:
    def __init__(self, analysis_):
        self.analysis = copy.deepcopy(analysis_)
        self.target_func = self.analysis.get_likelihood
        self.tolerance = 0.0000001 # default, can be changed using set_tolerance()
        self.maxiter = 100000 # default
        self.flush_rate = 20 # default
        self.seeds = []

    def set_options(self,opts):
        self.seeds = []
        self.stepsizes = []
        self.limits_low = []
        self.limits_high = []

        npars = self.analysis.get_npars()

        # check whether we have received options for all parameter
        if not len(opts)==npars:
            print("!!! FATAL: received unexpected number of parameters: {} instead of: {}".format(len(opts), npars))
            print("... exiting")
            sys.exit(1)

        # get model parameters and set user options
        names = list(self.analysis.get_par_names().keys())

        for i in range(npars):
            # check wether parameter name exists in user options
            if names[i] not in opts:
                print("!!! FATAL: parameter {} not found. Did you specify it correctly?".format(names[i]))
                print("!!! only found the following parameters:")
                print(opts)
                print("... exitting")
                sys.exit(1)

            # name exists
            self.seeds += [opts[names[i]].seed]
            self.stepsizes += [opts[names[i]].stepsize]
            self.limits_low += [opts[names[i]].limit_low]
            self.limits_high += [opts[names[i]].limit_high]

    def set_tolerance(self, tol):
        self.tolerance = tol

    def fit(self, werrors):
        npars = self.analysis.get_npars()
        self.bestpars = []
        par_names = list(self.analysis.get_par_names().keys())

        print("... commencing numerical minimization ({} parameters)".format(npars))
    
        # time for full likelihood evaluation
        start = time.time()
    
        bounds = [(i,j) for i,j in zip(self.limits_low, self.limits_high)]
        print(bounds)
        local_target_func = self.target_func
        res = minimize(local_target_func, self.seeds, bounds=bounds, tol=self.tolerance)
        print("... done")
        print("RESULTS (fitstatus: {})".format(res.status))
        #print(par_names)
        #print(res.x)
        [print(pn,result) for pn,result in zip(par_names,res.x)]

        # set seed to best fit (for subsequent minimizations)
    
        self.seeds = []
        self.bestpars = copy.deepcopy(res.x)
        self.seeds = copy.deepcopy(res.x)
    
        duration = time.time() - start
        print("likelihood value (gof): {}".format(self.analysis.get_likelihood_gof(res.x))) 
        print("likelihood value (abs): {}".format(self.analysis.get_likelihood_abs(res.x))) 
        print("time to find minimium and evaluate exact likelihood: {}s".format(duration))
        print("number of llh evaluations: {}".format(self.analysis.get_n_llh_evals()))
    
        if werrors:
            print("need to finish stats.py werror part")

    def get_bestpars(self,**pars):
        # check if number of parameters requested matches model parameters
        npars = self.analysis.get_npars()
        if not len(pars) == npars:
            print("!!! received npars={} but model has: {} parameters".format(len(par), npars))
            print("exiting ...")
            sys.exit(1)

        parameters = self.analysis.get_par_names()

        for it in parameters.items():
            # check if parameter exists
            if it[0] not in pars:
                print("!!! can not find model parameter: {}".format(it[0]))
                print("!!! the following parameters are found: ")
                print(pars)
                print("... exiting")
                sys.exit(1)
            pars[it[0]] = copy.deepcopy(self.bestpars[it[1]])

    def change_astro_model(self,astro):
        self.analysis.change_astro_model(astro)
        # need to update functor to update total number of parameter dimensions
        self.target_func = self.analysis.get_likelihood

    def get_profile_llh_value(self, point, verbose=True, scan=False):
        if verbose:
            npars = self.analysis.get_naprs()
            # this is a single profile_llh evaluation. not a scan!
            # be more verbose and do more error checking on user input
            # also since this is not a scan, we need to setup the minimizer here

            print("... requested profile llh")

            if npars<=len(point):
                # user input can not be profiled
                print("!!! FATAL: can not calculate profile likelihood for dim(point) >= number of model parameters. For dim(point) == number of model parameters use stats::get_llh_value(std::map<std::string, double> &point) instead!")
                print("... exiting")
                sys.exit(1)

        par_names = self.analysis.get_par_names()
        parameters = self.analysis.get_par_names()

        # figure out what parameters need to be fixed to a point
        indices = []
        for it in point.items():
            if it[0] not in parameters:
                # point specified by user concerns a parameter that does not exist :-(
                print("!!! FATAL: parameter {} not found. Are you sure this parameter is in the model?".format(it[0]))
                print("!!! the model has the following parameters:")
                print(parameters)
                print("... exiting")
                sys.exit(1)

            # user provided a valid point
            indices.append(parameters.get[it[0]])

        for i in range(npars):
            # fix requested variables
            if i in indices:
                print("... treating variable: {} with index {} as constant.".format(par_names[i], i))
                self.limits_low[i] = point.get(par_names[i])
                self.limits_high[i] = point.get(par_names[i])
               
        print("... commencing numerical minimization ( {} parameters)".format(npars-len(point)))

        # this be executed for both cases:
        # single point evaluation by user as well as single point evaluation as part of a llh scan
        
        bounds = [(i,j) for i,j in zip(self.limits_low, self.limits_high)]
        res = minimize(self.target_func, self.seeds, bounds=bounds, tol=self.tolerance)
        xs = res.x
        profile_llh = self.analysis.get_likelihood_gof(xs)

        npars = self.analysis.get_npars()
        par_names = self.analysis.get_par_names()

        if verbose:
            # this is not a scan
            # print some information for user
            print("... done")
            print("RESULTS (fitstatus: {})".format(res.status))
            print(res.x)
        else:
            # this is a scan. or toyfit.  need to update containers:
            for i in range(npars):
                point[par_names[i]] = xs[i]
            point["fit_status"] = res.status
            point["llh"] = profile_llh

            # set next seed point
            if scan:
                # the next scan point will be in neighborhood of this solution.
                # set seed to current best fit for next iteration, also update point with current solution
                if res.status==0:
                    self.seeds = copy.deepcopy(xs)
                else:
                    print(" ... fit failed. seed next scan point from global best fit.")
                    self.seeds = self.bestpars

    def get_llh_value(self, **point):
        npars = self.analysis.get_npars()
        if not npars==len(point):
            print("!!! FATAL: can not calculate profile likelihood for dim(point) != number pf model parameters")
            print("... exiting")
            sys.exit(1)

        # need to order parameters to pass to likelihood
        parameters = analysis.get_par_names()
        pars = []
