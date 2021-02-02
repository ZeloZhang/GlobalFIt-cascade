class astro_model_base:
    def get_par_names(self):
        return self.parameters

    def get_npars(self):
        return self.npars

    def get_model_name(self):
        return self.model_name
