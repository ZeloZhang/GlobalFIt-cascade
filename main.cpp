#include <cstdlib>
#include <string>
#include <map>
#include <TMath.h>
#include "include/stats.h"
#include "include/analysis.h"
#include "include/helpers.h"
#include "include/models/astro_model_single_plaw_wcutoff.h"
#include "include/bootstrap/toymc.h"



int main(int argc, char **argv)
{

	//std::string outdir=std::string("/data/user/zzhang1/fit_checkmanuel/output/");

	if (argc!=5) {
		std::cout << "expect four integer arguments: jobid and njobs indir outdir quitting!" << std::endl;
		return 1;
	}

	// get arguments relevant for profile scan
	int jobid = std::stoi(argv[1]);
        int random_seed = jobid;
	int njobs = std::stoi(argv[2]);
        std::string indir = argv[3];
        std::string outdir = argv[4];
//        int njobs = 1;

	/** create analysis. */
	NuFit::analysis wrapper;

	/** need to specifically create the analysis. */
	/** this function contains all the analysis specific details and can be edited in ./src/analysis.cpp */

	wrapper.create(indir); 

	/** get occasional prints of minimization progress from inside likelihood function */
	//wrapper.set_verbosity(true);

	/** use stats class to analyze NuFit::analysis */
	NuFit::stats min(wrapper);
	/** option to change model (NuFit::analysis uses single powerlaw as default) */
	//NuFit::astro_model_single_plaw_wcutoff *new_astro_model = new NuFit::astro_model_single_plaw_wcutoff();
	//min.change_astro_model(new_astro_model);

	
	/** specify seed, stepsize, limits for each parameter */
	/** parameter names need to match specification in model class, but parameter ordering does not matter! */
	/** e.g. ./src/astro_model_single_plaw.cpp + ./src/model_base.cpp */
	/** order of arguments: name, seed, stepsize, limit_low, limit_high */
	NuFit::helpers::par_options astro_norm("astro_norm", 1.8, 0.01, 0., 10.);
	NuFit::helpers::par_options astro_index("astro_index", 2.3, 0.01, 0, 10.0);
	//NuFit::helpers::par_options energy_cut("energy_cut", 6, 0.01, 3.0, 7.0);	
        NuFit::helpers::par_options muon_norm("muon_norm", 1.5, 0.01, 0.0, 10.0);
	NuFit::helpers::par_options muon_norm1("muon_norm_mlb", 1.0, 0.01, 0.0, 10.0);
	//NuFit::helpers::par_options muon_norm1("muon_norm_mlb", 1.0, 0.0, 0.0, 10.0);

	NuFit::helpers::par_options conv_norm("conv_norm", 0.9, 0.01, 0.0, 10.0);
	NuFit::helpers::par_options prompt_norm("prompt_norm", 0.2, 0.1, 0.0, 100.0);
        
        /*
        //original setting
        NuFit::helpers::par_options dom_eff("dom_efficiency", 1.0, 0.01, 0.0, 10.0);
	
        NuFit::helpers::par_options scattering("scattering", 1.0, 0.01, 0.0, 10.0);
	NuFit::helpers::par_options absorption("absorption", 1.0, 0.01, 0.0, 10.0);
	NuFit::helpers::par_options cr_index("delta_cr", 0.0, 0.01, -0.5, 0.5);
	NuFit::helpers::par_options holeice("zholeice_scattering", 1.3, 0.01, 0.0, 10.0);
        */

        /*
        // 4 analysis setting
        NuFit::helpers::par_options dom_eff("dom_efficiency", 1.025, 0.00, 1.025, 1.025);
        NuFit::helpers::par_options scattering("scattering", 1.02, 0.00, 1.02, 1.02);
	NuFit::helpers::par_options absorption("absorption", 1.033, 0.00, 1.033, 1.033);
	NuFit::helpers::par_options cr_index("delta_cr", 0.017, 0.00, 0.017, 0.017);
	NuFit::helpers::par_options holeice("zholeice_scattering", 1.72, 0.00, 1.72, 1.72);
        */
        
        // cscd only analysis setting
        /*
        NuFit::helpers::par_options dom_eff("dom_efficiency", 1.02, 0.0, 1.02, 1.02);
        NuFit::helpers::par_options scattering("scattering", 1.01, 0.0, 1.01, 1.01);
	NuFit::helpers::par_options absorption("absorption", 1.04, 0.00, 1.04, 1.04);
	NuFit::helpers::par_options cr_index("delta_cr", 0.33, 0.00, 0.33, 0.33);
	NuFit::helpers::par_options holeice("zholeice_scattering", 1.85, 0.0, 1.85, 1.85);
        */

	/** .. package everything */
	std::map<std::string, NuFit::helpers::par_options> options;
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_norm.name, astro_norm));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(astro_index.name, astro_index));
	//options.insert(std::pair<std::string, NuFit::helpers::par_options>(energy_cut.name, energy_cut));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(conv_norm.name, conv_norm));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(muon_norm.name, muon_norm));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(muon_norm1.name, muon_norm1));
	options.insert(std::pair<std::string, NuFit::helpers::par_options>(prompt_norm.name, prompt_norm));	
	
        //options.insert(std::pair<std::string, NuFit::helpers::par_options>(dom_eff.name, dom_eff));
	
        //options.insert(std::pair<std::string, NuFit::helpers::par_options>(scattering.name, scattering));
	//options.insert(std::pair<std::string, NuFit::helpers::par_options>(absorption.name, absorption));
	//options.insert(std::pair<std::string, NuFit::helpers::par_options>(cr_index.name, cr_index));
	//options.insert(std::pair<std::string, NuFit::helpers::par_options>(holeice.name, holeice));
        
	/** .. and send to minimizer */
	min.set_options(options);

	/** write seed histograms */
	std::map<std::string, double> pars;
	pars.insert(std::pair<std::string, double>(muon_norm.name, muon_norm.seed));
	pars.insert(std::pair<std::string, double>(muon_norm1.name, muon_norm1.seed));
	pars.insert(std::pair<std::string, double>(conv_norm.name, conv_norm.seed));
	pars.insert(std::pair<std::string, double>(prompt_norm.name, prompt_norm.seed));
	pars.insert(std::pair<std::string, double>(astro_norm.name, astro_norm.seed));
	pars.insert(std::pair<std::string, double>(astro_index.name, astro_index.seed));
	//pars.insert(std::pair<std::string, double>(energy_cut.name, energy_cut.seed));
	
        //pars.insert(std::pair<std::string, double>(dom_eff.name, dom_eff.seed));
	
        //pars.insert(std::pair<std::string, double>(scattering.name, scattering.seed));
	//pars.insert(std::pair<std::string, double>(absorption.name, absorption.seed));
	//pars.insert(std::pair<std::string, double>(cr_index.name, cr_index.seed));
	//pars.insert(std::pair<std::string, double>(holeice.name, holeice.seed));
        
        min.set_tolerance(1);

        min.fit(false); // ** if true -> get profile LLH errors after minimization from ROOT Minuit2
        //min.fit(true);
	
	//std::string outfile("./hists_seed.root");
	//wrapper.get_histograms(outfile, pars);

        //min.get_bestpars(pars);

        /** write best fit histograms */
	//std::string outfile = outdir+std::string("hists_fit.root");
        //wrapper.get_histograms(outfile, pars);
        /** example of how to do a profile llh scan */
        /*	
        int nsteps_index = 30;
        double index_min = 2.2;
        double index_max = 3;
        double ds = (index_max-index_min) / (nsteps_index-1);

        if (nsteps_index % njobs != 0) {
                std::cout << "please choose njobs to be an integer fraction of total number of steps" << std::endl;
                return 1;
        }
        int nsteps_job = nsteps_index / njobs;
        double index_min_tj = index_min + jobid * nsteps_job * ds;
        double index_max_tj = index_min_tj + (nsteps_job-1) * ds;
        NuFit::helpers::scan_options s_astro_norm("astro_norm", 30, 0.5, 2.5);
        NuFit::helpers::scan_options s_astro_index("astro_index", nsteps_job, index_min_tj, index_max_tj);
        

        //NuFit::helpers::scan_options s_astro_norm("astro_norm", 10, 1.5, 1.7);
        //NuFit::helpers::scan_options s_astro_index("astro_index", 10, 2.3, 2.6);
        std::map<std::string, NuFit::helpers::scan_options> scan;
        scan.insert(std::pair<std::string, NuFit::helpers::scan_options>(s_astro_norm.name, s_astro_norm));
        scan.insert(std::pair<std::string, NuFit::helpers::scan_options>(s_astro_index.name, s_astro_index));
        min.set_flush_rate(1); // set frequency of file dumps of scan results during scan
        min.scan_profile_llh(outdir+std::string("outfile_profile_llh_2d_part")+std::to_string(jobid)+std::string(".txt"), scan);
        */
        return 0;
}
