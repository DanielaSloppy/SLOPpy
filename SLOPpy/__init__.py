#!/usr/bin/env python
# -*- coding: utf-8 -*-
from SLOPpy.sloppy_run import sloppy_run
from SLOPpy.subroutines.io_subroutines import yaml_parser, pars_input
#from SLOPpy.subroutines.interpol import interpolate1d_grid_nocheck
from SLOPpy.prepare_datasets import *
from SLOPpy.sky_correction import *
from SLOPpy.differential_refraction_preparation import *
from SLOPpy.differential_refraction import *
from SLOPpy.check_differential_refraction import *
from SLOPpy.telluric_template import *
from SLOPpy.telluric_molecfit_v1_preparation import *
from SLOPpy.telluric_molecfit_v1_coadd import *
from SLOPpy.telluric_molecfit_preparation import *
from SLOPpy.telluric_molecfit_coadd import *
from SLOPpy.telluric_template_alternative import *
from SLOPpy.telluric_airmass_stellarRF import *
from SLOPpy.telluric_airmass_observerRF import *
from SLOPpy.telluric_airmass_observerRF_chunks import *
from SLOPpy.telluric_observerRF_skycalc import *
from SLOPpy.interstellar_lines import *
from SLOPpy.master_out import *
from SLOPpy.transmission_spectrum_preparation import *
from SLOPpy.transmission_spectrum import *
from SLOPpy.transmission_spectrum_average import *
from SLOPpy.transmission_spectrum_shortcuts import *
from SLOPpy.second_telluric_correction_on_transmission import *
#from SLOPpy.clv_rm_modelling import *
from SLOPpy.compare_clv_rm_effects import *
from SLOPpy.spectra_lightcurve import *
from SLOPpy.spectra_lightcurve_average import *
from SLOPpy.transmission_lightcurve import *
from SLOPpy.transmission_lightcurve_average import *
from SLOPpy.write_output_spectra import *
from SLOPpy.write_output_transmission import *

from SLOPpy.quick_transmission import *

# NEW
from SLOPpy.clv_rm_models import *
from SLOPpy.clv_rm_models_doubleprecision import *
from SLOPpy.clv_rm_models_lines import *
from SLOPpy.transmission_mcmc import *
from SLOPpy.transmission_binned_mcmc import *

from SLOPpy.pca_preparation import *
from SLOPpy.sysrem_correction import *

__version__ = "1.5.1"
