#include "../include/analysis.h"

NuFit::analysis::analysis() {

	verbose_llh = false;	
	n_llh_evals = 0;

}

NuFit::analysis::~analysis() {
	delete model;
}

void NuFit::analysis::create(std::string defdir)
{

        std::cout << "... start wrapping analyses" << std::endl;
        std::cout << std::endl;

        // here goes the analysis code (parsing of data and mc)

        // define cascade signal sample
        std::vector<double> binsx; // logE
	for (unsigned int i=0; i<23; ++i)	
		binsx.push_back(2.6 + 0.2 * i);

        std::vector<double> binsy; // cosz (2 bins: Northern Sky, Southern Sky)
        binsy.push_back(-2.0);
        binsy.push_back(0.2);
	binsy.push_back(0.6);
        binsy.push_back(2.0);

        std::vector<double> binsz; // ra (1 bin)
        binsz.push_back(-10.);
        binsz.push_back(10.);

	//std::string defdir("/data/user/zzhang1/fit_checkmanuel/input_hanssim/");
	std::string dir_myc=defdir+std::string("myc/");
	std::string dir_mlb=defdir+std::string("mlb/");
	std::string dir_as=defdir+std::string("as/");

	// flasher holeice p2=0.0
	std::string dir_myc_HKKMS06=dir_myc+std::string("HI_30cm/");
	std::string dir_mlb_HKKMS06=dir_mlb+std::string("HI_30cm/");
	
	//SIBYLL2.3c 
        //std::string dir30=defdir+std::string("mceq/SIBYLL2.3c/holeice30/");	
	
	//DPMJETIII
	//std::string dir30=defdir+std::string("mceq/DPMJETIII/holeice30/");
	
	//EPOSLHC
	//std::string dir30=defdir+std::string("mceq/EPOSLHC/holeice30/");
	
	//QGSJETII04
	//std::string dir30=defdir+std::string("mceq/QGSJETII04/holeice30/");


	
        NuFit::hists analysis1("cascade_all", binsx, binsy, binsz); 
        std::string name_nue("neutrino/nue_cscd_all_neutrino.txt");
        std::string name_numu("neutrino/numu_cscd_all_neutrino.txt");
        std::string name_nutau("neutrino/nutau_cscd_all_neutrino.txt");
	std::string name_mu("HI_50cm/neutrino/mgun_cscd_all_neutrino.txt"); 
	std::string name_data("data/txt/data_cascade_dustcorr.txt");

        analysis1.read(dir_myc_HKKMS06+name_nue, dir_myc_HKKMS06+name_numu, dir_myc_HKKMS06+name_nutau, dir_myc+name_mu, dir_myc+name_data);
        
        // define muon background sample
        std::vector<double> binsx2; // logE
	binsx2.push_back(2.6);
	binsx2.push_back(4.778);

        std::vector<double> binsy2; // cosz (1 bin: all-sky)
        binsy2.push_back(-2.0);
        binsy2.push_back(2.0);

        std::vector<double> binsz2; // ra (1 bin)
        binsz2.push_back(-10.);
        binsz2.push_back(10.);

        NuFit::hists analysis2("muon", binsx2, binsy2, binsz2); 
        std::string name_nue2("muon/muon_nue.txt");
        std::string name_numu2("muon/muon_numu.txt");
        std::string name_nutau2("muon/muon_nutau.txt");
        std::string name_mu2("HI_50cm/muon/muon_mgun.txt");
        std::string name_data2("data/txt/data_muon.txt");
        analysis2.read(dir_myc_HKKMS06+name_nue2, dir_myc_HKKMS06+name_numu2, dir_myc_HKKMS06+name_nutau2, dir_myc+name_mu2, dir_myc+name_data2);
        
	// define numu control sample
	std::vector<double> binsx3;
        for (unsigned int i=0; i<12; ++i)
             binsx3.push_back(2.6 + 0.2 * i);

        /*
        std::vector<double> binsy2; // cosz (1 bin: all-sky)
        binsy2.push_back(-2.0);
        binsy2.push_back(2.0);
        */

        NuFit::hists analysis3("hybrid", binsx3, binsy2, binsz);
        std::string name_nue3("neutrino/nue_hybrid_neutrino.txt");
        std::string name_numu3("neutrino/numu_hybrid_neutrino.txt");
        std::string name_nutau3("neutrino/nutau_hybrid_neutrino.txt");
        std::string name_mu3("HI_50cm/neutrino/mgun_hybrid_neutrino.txt");
        std::string name_data3("data/txt/data_hybrid.txt");
        analysis3.read(dir_myc_HKKMS06+name_nue3, dir_myc_HKKMS06+name_numu3, dir_myc_HKKMS06+name_nutau3, dir_myc+name_mu3, dir_myc+name_data3);


	// define mlb sample
        std::vector<double> binsx4; // logE
        for (unsigned int i=0; i<16; ++i)
                binsx4.push_back(4.0 + 0.2 * i);

	NuFit::hists analysis4("cascade_mlb", binsx4, binsy, binsz);
        std::string name_nue4("txt/nue_mlb_neutrino.txt");
        std::string name_numu4("txt/numu_mlb_neutrino.txt");
        std::string name_nutau4("txt/nutau_mlb_neutrino.txt");
        std::string name_mu4("HI_50cm/txt/cors_mlb.txt");
        std::string name_data4("data/txt/data_mlb_dustcorr.txt");

        //analysis4.read(dir_mlb_HKKMS06+name_nue4, dir_mlb_HKKMS06+name_numu4, dir_mlb_HKKMS06+name_nutau4, dir_mlb+name_mu4, dir_mlb+name_data4);
	/*
	// define as sample
	
        std::vector<double> binsx_as;
        for (int i=0; i<8; i++)
                binsx_as.push_back(4.53 + 0.352857 * i);

        std::vector<double> binsy_as;
        binsy_as.push_back(-2.0);
        binsy_as.push_back(2.0);

        NuFit::hists analysis5("cascade_as", binsx_as, binsy_as, binsz);
        std::string name_nue5("as_nue.txt");
        std::string name_numu5("as_numu.txt");
        std::string name_nutau5("as_nutau.txt");
        std::string name_mu5("as_cors.txt"); 
        std::string name_data5("as_data_SPE.txt");

        //analysis5.read(dir_as+name_nue5, dir_as+name_numu5, dir_as+name_nutau5, dir_as+name_mu5, dir_as+name_data5);
	*/
        
        
        std::vector<NuFit::hists> analyses;
        analyses.push_back(analysis1);
        analyses.push_back(analysis2);
        analyses.push_back(analysis3);
	//analyses.push_back(analysis4);
	//analyses.push_back(analysis5);

        // create astro model. then create base model
        NuFit::astro_model_single_plaw *astro = new astro_model_single_plaw();
        //NuFit::astro_model_plaw_singlep *astro = new astro_model_plaw_singlep();
        //base_model needs to know the input data as well as the astro model
        NuFit::model_base *mymodel = new model_base(analyses, astro);
	model = mymodel;

	
	// comment out code below if you don't need systematics		
	/*
         	 
	std::vector<std::string> analysis_names;
	std::map<std::string, NuFit::hists*> map_analyses;
	for (unsigned int i=0; i<analyses.size(); ++i) 
	{
		analysis_names.push_back(analyses[i].name);
		map_analyses.insert(std::pair<std::string, NuFit::hists*>(analyses[i].name, &(analyses[i])));
	}

	std::map<std::string, NuFit::interpolated_par *> systematics;

	
	// start with dom efficiency
	std::string sysname_eff("dom_efficiency");
	NuFit::interpolated_par_domeff *domeff = new interpolated_par_domeff(sysname_eff, analysis_names, map_analyses);

	
	// scattering
	std::string sysname_scat("scattering");
	NuFit::interpolated_par_scattering *scattering = new interpolated_par_scattering(sysname_scat, analysis_names, map_analyses);

	// absorption
	std::string sysname_abs("absorption");
	NuFit::interpolated_par_absorption *absorption = new interpolated_par_absorption(sysname_abs, analysis_names, map_analyses);
	
	std::string sysname_hi("zholeice_scattering");
	NuFit::interpolated_par_holeice *holeice = new interpolated_par_holeice(sysname_hi, analysis_names, map_analyses);
        
	// all objects in this map will be destroyed by model_base_sys later on
	systematics.insert(std::pair<std::string, NuFit::interpolated_par *> (sysname_eff, domeff));
	systematics.insert(std::pair<std::string, NuFit::interpolated_par *> (sysname_scat, scattering));
	systematics.insert(std::pair<std::string, NuFit::interpolated_par *> (sysname_abs, absorption));
	systematics.insert(std::pair<std::string, NuFit::interpolated_par *> (sysname_hi, holeice));
	
	
	NuFit::model_base_sys *mymodel = new model_base_sys(analyses, astro, systematics, 10);
      	model = mymodel; 
	*/
	
        // end of analysis code

        std::cout << std::endl;
        std::cout << "... wrapping done" << std::endl; 
}

// the functions below interface to the model class

void NuFit::analysis::change_astro_model(NuFit::astro_model_base *astro)
{
	model->change_astro_model(astro);
	return;
}

double NuFit::analysis::get_likelihood(const double *pars)
{
	// this is being minimized when stats class is used (interface to ROOT Minuit2)
	double neglogl = model->likelihood(pars);
	++n_llh_evals;

	// add verbosity if requested
	if(verbose_llh) {
		if(n_llh_evals%50==0 || n_llh_evals==1) {
			std::cout << "... evaluated llh " << n_llh_evals << " times already." << std::endl;
			std::cout << "f( ";
			for(unsigned int i=0; i<get_npars(); ++i)
				std::cout << pars[i] << " ";
			std::cout << ") = " << neglogl << std::endl;
			std::cout << std::endl;
		}
		
	}

	// factor of 2 for Wilk's theorem
	return 2.0 * neglogl;
}

double NuFit::analysis::get_likelihood_gof(const double *pars)
{
	// factor of 2 for Wilk's theorem
	return 2.0 * model->likelihood_gof(pars);
}

double NuFit::analysis::get_likelihood_abs(const double *pars)
{
	// factor of 2 for Wilk's theorem
	return 2.0 * model->likelihood_abs(pars);
}

double NuFit::analysis::get_lnprob(boost::python::numeric::array pars_)
{
	// assume that pars is of length npars
	unsigned int npars = get_npars();	
	double pars[npars];
	for (unsigned int i=0; i<npars; ++i)
		pars[i] = boost::python::extract<double>(pars_[i]);

	// make log-likelihood positive
	return (-1.0) * model->likelihood_abs(pars); 
}

unsigned int NuFit::analysis::get_n_llh_evals() 
{
	return n_llh_evals;
}

unsigned int NuFit::analysis::get_npars()
{
	return model->get_npars();
}

void NuFit::analysis::get_par_names(std::vector<std::string> &names)
{
	model->get_par_names(names);
	return;
}

void NuFit::analysis::get_par_names(std::map<std::string, unsigned int> &names)
{
        model->get_par_names(names);
        return;
}

void NuFit::analysis::get_histograms(std::string outfile, std::map<std::string, double> &pars)
{
	model->get_histograms(outfile, pars);
	return;
}

std::vector<TH3D*> NuFit::analysis::get_hist_mcsum(std::map<std::string, double> &pars)
{ 
        return model->get_hist_mcsum(pars);
}

void NuFit::analysis::reset_n_llh_evals() 
{
	n_llh_evals=0;
	return;
}

void NuFit::analysis::set_verbosity(bool flag=false) {
	verbose_llh=flag;
	return;
}

std::vector<std::string> NuFit::analysis::get_analysis_names()
{
	return model->get_analysis_names();
}

void NuFit::analysis::set_hist(std::string analysis, TH3D* hist)
{
	model -> set_hist(analysis, hist);
	return;
}

void NuFit::analysis::cache_data_hists()
{
	model -> cache_data_hists();
	return;
}

void NuFit::analysis::restore_data_hists()
{
	model -> restore_data_hists();
	return;
}

void NuFit::analysis::update_auxillary_data(std::map<std::string, double> &pars)
{
        model -> update_auxillary_data(pars);
        return;
}

void NuFit::analysis::reset_auxillary_data()
{
        model -> reset_auxillary_data();
        return;
}


#include <boost/python.hpp>
using namespace boost::python;
BOOST_PYTHON_MODULE(analysis)
{
        // expose functions create
	boost::python::numeric::array::set_module_and_type("numpy", "ndarray");
        class_<NuFit::analysis>("analysis", init<>())	
                .def("create", &NuFit::analysis::create)
		.def("get_lnprob", &NuFit::analysis::get_lnprob)
        ;
}
