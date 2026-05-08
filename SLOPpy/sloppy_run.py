import SLOPpy
import argparse
import os
import sys
import collections

def sloppy_run(file_conf=None ):

    print()
    print('SLOPpy v{0}'.format(SLOPpy.__version__))
    print()
    print('Python version in use:')
    print(sys.version)
    #if sys.version_info[0] == 3 and sys.version_info[1] > 7:
    #    print('WARNING MESSAGES SUPPRESSED!')
    #print()

    if file_conf is None:
        parser = argparse.ArgumentParser(prog='SLOPpy_Run', description='SLOPpy runner')
        parser.add_argument('config_file', type=str, nargs=1, help='config file')
    
        args = parser.parse_args()
        file_conf = args.config_file[0]

    config_in = SLOPpy.yaml_parser(file_conf)

    SLOPpy.pars_input(config_in)

    #return()

    print()
    """ creation of the pickle files """
    SLOPpy.prepare_datasets(config_in)

    #""" Retrieving the dictionary with the pipeline recipes """
    #pipeline = config_in['pipeline']

    """ Recipes must be performed in a given order... that's why we must use and ordered dictionary"""
    """ Each of the following recipes has to be performed on the whole spectrum """
    pipeline_common_routines = collections.OrderedDict()

    pipeline_common_routines['sky_correction'] = SLOPpy.compute_sky_correction

    #pipeline_common_routines['PCA_test01'] = SLOPpy.PCA_test01

    pipeline_common_routines['pca_preparation'] = SLOPpy.compute_pca_preparation
    pipeline_common_routines['sysrem_correction'] = SLOPpy.compute_sysrem_correction


    # molecfit version 1.5
    pipeline_common_routines['telluric_molecfit_v1_preparation'] = SLOPpy.compute_telluric_molecfit_v1_preparation
    pipeline_common_routines['telluric_molecfit_v1'] = SLOPpy.compute_telluric_molecfit_v1
    pipeline_common_routines['telluric_molecfit_v1_coadd'] = SLOPpy.compute_telluric_molecfit_v1_coadd

    # molecfit new version
    pipeline_common_routines['telluric_molecfit_preparation'] = SLOPpy.compute_telluric_molecfit_preparation
    pipeline_common_routines['telluric_molecfit'] = SLOPpy.compute_telluric_molecfit
    pipeline_common_routines['telluric_molecfit_coadd'] = SLOPpy.compute_telluric_molecfit_coadd

    pipeline_common_routines['telluric_template'] = SLOPpy.compute_telluric_template
    pipeline_common_routines['telluric_template_reference'] = SLOPpy.compute_telluric_template_reference
    pipeline_common_routines['telluric_template_alternative'] = SLOPpy.compute_telluric_template_alternative

    pipeline_common_routines['telluric_airmass_stellarRF'] = SLOPpy.compute_telluric_airmass_stellarRF
    pipeline_common_routines['telluric_airmass_reference_stellarRF'] = SLOPpy.compute_telluric_airmass_reference_stellarRF

    pipeline_common_routines['telluric_airmass_observerRF'] = SLOPpy.compute_telluric_airmass_observerRF
    pipeline_common_routines['telluric_airmass_berv_observerRF'] = SLOPpy.compute_telluric_airmass_berv_observerRF
    pipeline_common_routines['telluric_airmass_reference_observerRF'] = SLOPpy.compute_telluric_airmass_reference_observerRF
    pipeline_common_routines['telluric_airmass_berv_reference_observerRF'] = SLOPpy.compute_telluric_airmass_berv_reference_observerRF


    #pipeline_routines['telluric_obsolete_wyttenbach'] = SLOPpy.compute_telluric_obsolete_wyttenbach
    #pipeline_routines['telluric_airmass_observerRF_chunks'] = SLOPpy.compute_telluric_airmass_observerRF_chunks
    pipeline_common_routines['telluric_observerRF_skycalc'] = SLOPpy.compute_telluric_observerRF_skycalc

    pipeline_common_routines['differential_refraction'] = SLOPpy.compute_differential_refraction
    pipeline_common_routines['differential_refraction_update'] = SLOPpy.compute_differential_refraction_update

    pipeline_common_routines['check_differential_refraction'] = SLOPpy.check_differential_refraction
    pipeline_common_routines['write_differential_refraction'] = SLOPpy.write_differential_refraction

    pipeline_common_routines['interstellar_lines'] = SLOPpy.compute_interstellar_lines

    pipeline_common_routines['master_out'] = SLOPpy.compute_master_out

    pipeline_common_routines['clv_rm_models'] = SLOPpy.compute_clv_rm_models
    pipeline_common_routines['clv_rm_models_doubleprecision'] = SLOPpy.compute_clv_rm_models_doubleprecision
    pipeline_common_routines['transmission_spectrum_preparation'] = SLOPpy.compute_transmission_spectrum_preparation
    pipeline_common_routines['wiggle_correction'] = SLOPpy.compute_wiggle_correction

    pipeline_common_routines['write_output_transmission'] = SLOPpy.write_output_transmission
    pipeline_common_routines['write_output_transmission_stellarRF'] = SLOPpy.write_output_transmission_stellarRF
    pipeline_common_routines['write_output_transmission_planetRF'] = SLOPpy.write_output_transmission_planetRF
    pipeline_common_routines['write_output_transmission_observerRF'] = SLOPpy.write_output_transmission_observerRF


    """ Legacy routines for testing purposes """
    pipeline_routines = collections.OrderedDict()


    #pipeline_routines['transmission_spectrum_planetRF'] = SLOPpy.compute_transmission_spectrum_planetRF
    #pipeline_routines['transmission_spectrum_observerRF'] = SLOPpy.compute_transmission_spectrum_observerRF
    #pipeline_routines['transmission_spectrum_stellarRF'] = SLOPpy.compute_transmission_spectrum_stellarRF
    #pipeline_routines['transmission_spectrum'] = SLOPpy.compute_transmission_spectrum

    #pipeline_routines['second_telluric_correction_on_transmission'] = SLOPpy.compute_second_telluric_correction_on_transmission

    #pipeline_routines['transmission_map'] = SLOPpy.compute_transmission_map
    #pipeline_routines['transmission_clv_rm_map'] = SLOPpy.compute_transmission_clv_rm_map


    """ Each of the following recipes has to be performed independently on each
    set of spectral lines """

    pipeline_lines_routines = collections.OrderedDict()

    """
    pipeline_lines_routines['transmission_spectrum_planetRF'] = SLOPpy.compute_transmission_spectrum_planetRF
    pipeline_lines_routines['transmission_spectrum_observerRF'] = SLOPpy.compute_transmission_spectrum_observerRF
    pipeline_lines_routines['transmission_spectrum_stellarRF'] = SLOPpy.compute_transmission_spectrum_stellarRF
    pipeline_lines_routines['transmission_spectrum'] = SLOPpy.compute_transmission_spectrum

    pipeline_lines_routines['second_telluric_correction_on_transmission'] = SLOPpy.compute_second_telluric_correction_on_transmission


    pipeline_lines_routines['transmission_map'] = SLOPpy.compute_transmission_map
    pipeline_lines_routines['transmission_clv_rm_map'] = SLOPpy.compute_transmission_clv_rm_map

    #pipeline_lines_routines['clv_rm_modelling'] = SLOPpy.compute_clv_rm_modelling
    """

    # ! NEW
    pipeline_lines_routines['quick_transmission'] = SLOPpy.compute_quick_transmission
    pipeline_lines_routines['clv_rm_models_lines'] = SLOPpy.compute_clv_rm_models_lines
    pipeline_lines_routines['transmission_mcmc'] = SLOPpy.compute_transmission_mcmc
    pipeline_lines_routines['transmission_mcmc_iterative'] = SLOPpy.compute_transmission_mcmc_iterative
    pipeline_lines_routines['transmission_binned_mcmc'] = SLOPpy.compute_transmission_binned_mcmc
    pipeline_lines_routines['transmission_binned_mcmc_iterative'] = SLOPpy.compute_transmission_binned_mcmc_iterative

    pipeline_lines_routines['transmission_spectrum_planetRF'] = SLOPpy.compute_transmission_spectrum_planetRF
    pipeline_lines_routines['transmission_spectrum_observerRF'] = SLOPpy.compute_transmission_spectrum_observerRF
    pipeline_lines_routines['transmission_spectrum_stellarRF'] = SLOPpy.compute_transmission_spectrum_stellarRF
    pipeline_lines_routines['transmission_spectrum'] = SLOPpy.compute_transmission_spectrum

    pipeline_lines_routines['transmission_spectrum_planetRF_iterative'] = SLOPpy.compute_transmission_spectrum_planetRF_iterative
    pipeline_lines_routines['transmission_spectrum_observerRF_iterative'] = SLOPpy.compute_transmission_spectrum_observerRF_iterative
    pipeline_lines_routines['transmission_spectrum_stellarRF_iterative'] = SLOPpy.compute_transmission_spectrum_stellarRF_iterative
    pipeline_lines_routines['transmission_spectrum_iterative'] = SLOPpy.compute_transmission_spectrum_iterative


    pipeline_lines_routines['transmission_spectrum_average_planetRF'] = SLOPpy.compute_transmission_spectrum_average_planetRF
    pipeline_lines_routines['transmission_spectrum_average_observerRF'] = SLOPpy.compute_transmission_spectrum_average_observerRF
    pipeline_lines_routines['transmission_spectrum_average_stellarRF'] = SLOPpy.compute_transmission_spectrum_average_stellarRF
    pipeline_lines_routines['transmission_spectrum_average'] = SLOPpy.compute_transmission_spectrum_average

    pipeline_lines_routines['transmission_spectrum_average_planetRF_iterative'] = SLOPpy.compute_transmission_spectrum_average_planetRF
    pipeline_lines_routines['transmission_spectrum_average_observerRF_iterative'] = SLOPpy.compute_transmission_spectrum_average_observerRF
    pipeline_lines_routines['transmission_spectrum_average_stellarRF_iterative'] = SLOPpy.compute_transmission_spectrum_average_stellarRF
    pipeline_lines_routines['transmission_spectrum_average_iterative'] = SLOPpy.compute_transmission_spectrum_average


    pipeline_lines_routines['transmission_lightcurve'] = SLOPpy.compute_transmission_lightcurve
    pipeline_lines_routines['transmission_lightcurve_average'] = SLOPpy.compute_transmission_lightcurve_average

    pipeline_lines_routines['spectra_lightcurve'] = SLOPpy.compute_spectra_lightcurve
    pipeline_lines_routines['spectra_lightcurve_average'] = SLOPpy.compute_spectra_lightcurve_average


    # TODO: to be updated to support single line(s) set
    """
    pipeline_lines_routines['transmission_clv_rm_correction_planetRF'] = SLOPpy.compute_transmission_clv_rm_correction_planetRF
    pipeline_lines_routines['transmission_clv_rm_correction_observerRF'] = SLOPpy.compute_transmission_clv_rm_correction_observerRF
    pipeline_lines_routines['transmission_clv_rm_correction_stellarRF'] = SLOPpy.compute_transmission_clv_rm_correction_stellarRF
    pipeline_lines_routines['transmission_clv_rm_correction'] = SLOPpy.compute_transmission_clv_rm_correction

    pipeline_lines_routines['transmission_clv_rm_average_planetRF'] = SLOPpy.compute_transmission_clv_rm_average_planetRF
    pipeline_lines_routines['transmission_clv_rm_average_observerRF'] = SLOPpy.compute_transmission_clv_rm_average_observerRF
    pipeline_lines_routines['transmission_clv_rm_average_stellarRF'] = SLOPpy.compute_transmission_clv_rm_average_stellarRF
    pipeline_lines_routines['transmission_clv_rm_average'] = SLOPpy.compute_transmission_clv_rm_average



    pipeline_lines_routines['spectra_lightcurve'] = SLOPpy.compute_spectra_lightcurve
    pipeline_lines_routines['spectra_lightcurve_average'] = SLOPpy.compute_spectra_lightcurve_average

    pipeline_lines_routines['excess_lightcurve'] = SLOPpy.compute_spectra_lightcurve
    pipeline_lines_routines['excess_lightcurve_average'] = SLOPpy.compute_spectra_lightcurve_average

    pipeline_lines_routines['spectra_lightcurve_clv_rm_correction'] = SLOPpy.compute_spectra_lightcurve_clv_rm_correction
    pipeline_lines_routines['spectra_lightcurve_average_clv_rm_correction'] = SLOPpy.compute_spectra_lightcurve_average_clv_rm_correction

    pipeline_lines_routines['excess_lightcurve_clv_rm_correction'] = SLOPpy.compute_spectra_lightcurve_clv_rm_correction
    pipeline_lines_routines['excess_lightcurve_average_clv_rm_correction'] = SLOPpy.compute_spectra_lightcurve_average_clv_rm_correction



    pipeline_lines_routines['transmission_lightcurve_planetRF'] = SLOPpy.compute_transmission_lightcurve_planetRF
    pipeline_lines_routines['transmission_lightcurve_observerRF'] = SLOPpy.compute_transmission_lightcurve_observerRF
    pipeline_lines_routines['transmission_lightcurve_stellarRF'] = SLOPpy.compute_transmission_lightcurve_stellarRF
    pipeline_lines_routines['transmission_lightcurve'] = SLOPpy.compute_transmission_lightcurve

    pipeline_lines_routines['write_output_spectra'] = SLOPpy.write_output_spectra

    pipeline_lines_routines['transmission_lightcurve_average_planetRF'] = SLOPpy.compute_transmission_lightcurve_average_planetRF
    pipeline_lines_routines['transmission_lightcurve_average_observerRF'] = SLOPpy.compute_transmission_lightcurve_average_observerRF
    pipeline_lines_routines['transmission_lightcurve_average_stellarRF'] = SLOPpy.compute_transmission_lightcurve_average_stellarRF
    pipeline_lines_routines['transmission_lightcurve_average'] = SLOPpy.compute_transmission_lightcurve_average
    """


    plot_preparation_routines = collections.OrderedDict()
    #plot_preparation_routines['clv_rm_modelling'] = SLOPpy.plot_clv_rm_modelling

    # ! New
    plot_preparation_routines['clv_rm_models'] = SLOPpy.plot_clv_rm_models
    plot_preparation_routines['clv_rm_models_doubleprecision'] = SLOPpy.plot_clv_rm_models


    plot_routines = collections.OrderedDict()

    plot_routines['plot_dataset'] = SLOPpy.plot_dataset
    plot_routines['prepare_dataset'] = SLOPpy.plot_dataset
    plot_routines['dataset'] = SLOPpy.plot_dataset
    plot_routines['sky_correction'] = SLOPpy.plot_sky_correction

    plot_routines['differential_refraction'] = SLOPpy.plot_differential_refraction
    plot_routines['differential_refraction_update'] = SLOPpy.plot_differential_refraction_update

    plot_routines['check_differential_refraction'] = SLOPpy.plot_check_differential_refraction
    #plot_routines['write_differential_refraction'] = SLOPpy.write_differential_refraction

    #plot_routines['PCA_test01'] = SLOPpy.PCA_test01

    # molecfit version 1.5
    plot_routines['telluric_molecfit_v1'] = SLOPpy.plot_telluric_molecfit_v1
    plot_routines['telluric_molecfit_v1_coadd'] = SLOPpy.plot_telluric_molecfit_v1_coadd

    # molecfit new version
    plot_routines['telluric_molecfit'] = SLOPpy.plot_telluric_molecfit
    plot_routines['telluric_molecfit_coadd'] = SLOPpy.plot_telluric_molecfit_coadd

    plot_routines['telluric_template'] = SLOPpy.plot_telluric_template
    plot_routines['telluric_template_reference'] = SLOPpy.plot_telluric_template_reference
    plot_routines['telluric_template_alternative'] = SLOPpy.plot_telluric_template_alternative

    plot_routines['telluric_airmass_stellarRF'] = SLOPpy.plot_telluric_airmass_stellarRF
    plot_routines['telluric_airmass_reference_stellarRF'] = SLOPpy.plot_telluric_airmass_reference_stellarRF

    plot_routines['telluric_airmass_observerRF'] = SLOPpy.plot_telluric_airmass_observerRF
    plot_routines['telluric_airmass_berv_observerRF'] = SLOPpy.plot_telluric_airmass_berv_observerRF
    plot_routines['telluric_airmass_reference_observerRF'] = SLOPpy.plot_telluric_airmass_reference_observerRF
    plot_routines['telluric_airmass_berv_reference_observerRF'] = SLOPpy.plot_telluric_airmass_berv_reference_observerRF

    #plot_routines['telluric_obsolete_wyttenbach'] = SLOPpy.plot_telluric_obsolete_wyttenbach
    #plot_routines['telluric_airmass_observerRF_chunks'] = SLOPpy.plot_telluric_airmass_observerRF_chunks
    plot_routines['telluric_observerRF_skycalc'] = SLOPpy.plot_telluric_observerRF_skycalc

    plot_routines['interstellar_lines'] = SLOPpy.plot_interstellar_lines

    plot_routines['master_out'] = SLOPpy.plot_master_out

    plot_routines['telluric_molecfit_preparation'] = SLOPpy.plot_telluric_molecfit_preparation

    # ! NEW
    plot_routines['transmission_spectrum_preparation'] = SLOPpy.plot_transmission_spectrum_preparation
    plot_routines['wiggle_correction'] = SLOPpy.plot_wiggle_correction

    """
    plot_routines['transmission_spectrum_planetRF'] = SLOPpy.plot_transmission_spectrum_planetRF
    plot_routines['transmission_spectrum_observerRF'] = SLOPpy.plot_transmission_spectrum_observerRF
    plot_routines['transmission_spectrum_stellarRF'] = SLOPpy.plot_transmission_spectrum_stellarRF
    plot_routines['transmission_spectrum'] = SLOPpy.plot_transmission_spectrum

    plot_routines['second_telluric_correction_on_transmission'] = SLOPpy.plot_second_telluric_correction_on_transmission


    plot_routines['transmission_clv_rm_correction_planetRF'] = SLOPpy.plot_transmission_clv_rm_correction_planetRF
    plot_routines['transmission_clv_rm_correction_observerRF'] = SLOPpy.plot_transmission_clv_rm_correction_observerRF
    plot_routines['transmission_clv_rm_correction_stellarRF'] = SLOPpy.plot_transmission_clv_rm_correction_stellarRF
    plot_routines['transmission_clv_rm_correction'] = SLOPpy.plot_transmission_clv_rm_correction


    plot_routines['spectra_lightcurve'] = SLOPpy.plot_spectra_lightcurve
    plot_routines['excess_lightcurve'] = SLOPpy.plot_spectra_lightcurve

    plot_routines['spectra_lightcurve_clv_rm_correction'] = SLOPpy.plot_spectra_lightcurve_clv_rm_correction
    plot_routines['excess_lightcurve_clv_rm_correction'] = SLOPpy.plot_spectra_lightcurve_clv_rm_correction


    plot_routines['transmission_lightcurve_planetRF'] = SLOPpy.plot_transmission_lightcurve_planetRF
    plot_routines['transmission_lightcurve_observerRF'] = SLOPpy.plot_transmission_lightcurve_observerRF
    plot_routines['transmission_lightcurve_stellarRF'] = SLOPpy.plot_transmission_lightcurve_stellarRF
    plot_routines['transmission_lightcurve'] = SLOPpy.plot_transmission_lightcurve

    plot_routines['transmission_map'] = SLOPpy.plot_transmission_map
    plot_routines['transmission_clv_rm_map'] = SLOPpy.plot_transmission_clv_rm_map
    """

    plot_lines_routines = collections.OrderedDict()
    plot_lines_routines['clv_rm_models_lines'] = SLOPpy.plot_clv_rm_models_lines


    plot_lines_routines['transmission_binned_mcmc'] = SLOPpy.plot_transmission_binned_mcmc

    plot_lines_routines['transmission_spectrum_planetRF'] = SLOPpy.plot_transmission_spectrum_planetRF
    plot_lines_routines['transmission_spectrum_observerRF'] = SLOPpy.plot_transmission_spectrum_observerRF
    plot_lines_routines['transmission_spectrum_stellarRF'] = SLOPpy.plot_transmission_spectrum_stellarRF
    plot_lines_routines['transmission_spectrum'] = SLOPpy.plot_transmission_spectrum

    plot_lines_routines['transmission_spectrum_planetRF_iterative'] = SLOPpy.plot_transmission_spectrum_planetRF_iterative
    plot_lines_routines['transmission_spectrum_observerRF_iterative'] = SLOPpy.plot_transmission_spectrum_observerRF_iterative
    plot_lines_routines['transmission_spectrum_stellarRF_iterative'] = SLOPpy.plot_transmission_spectrum_stellarRF_iterative
    plot_lines_routines['transmission_spectrum_iterative'] = SLOPpy.plot_transmission_spectrum_iterative





    plot_average_routines = collections.OrderedDict()
    plot_average_routines['compare_master_out'] = SLOPpy.plot_compare_master_out


    plot_lines_average_routines = collections.OrderedDict()

    # ! These should be removed and performed line by line !

    plot_lines_average_routines['transmission_binned_mcmc'] = SLOPpy.plot_transmission_binned_mcmc
    plot_lines_average_routines['transmission_spectrum_average_planetRF'] = SLOPpy.plot_transmission_spectrum_average_planetRF
    plot_lines_average_routines['transmission_spectrum_average_observerRF'] = SLOPpy.plot_transmission_spectrum_average_observerRF
    plot_lines_average_routines['transmission_spectrum_average_stellarRF'] = SLOPpy.plot_transmission_spectrum_average_stellarRF
    plot_lines_average_routines['transmission_spectrum_average'] = SLOPpy.plot_transmission_spectrum_average


    """
    plot_average_routines['excess_lightcurve_average'] = SLOPpy.plot_spectra_lightcurve_average
    plot_average_routines['spectra_lightcurve_average'] = SLOPpy.plot_spectra_lightcurve_average

    plot_average_routines['spectra_lightcurve_average_clv_rm_correction'] = \
        SLOPpy.plot_spectra_lightcurve_average_clv_rm_correction
    plot_average_routines['excess_lightcurve_average_clv_rm_correction'] = \
        SLOPpy.plot_spectra_lightcurve_average_clv_rm_correction

    plot_average_routines['transmission_average_planetRF'] = SLOPpy.plot_transmission_average_planetRF
    plot_average_routines['transmission_average_observerRF'] = SLOPpy.plot_transmission_average_observerRF
    plot_average_routines['transmission_average_stellarRF'] = SLOPpy.plot_transmission_average_stellarRF
    plot_average_routines['transmission_average'] = SLOPpy.plot_transmission_average

    plot_average_routines['transmission_lightcurve_average_planetRF'] = SLOPpy.plot_transmission_lightcurve_average_planetRF
    plot_average_routines['transmission_lightcurve_average_observerRF'] = SLOPpy.plot_transmission_lightcurve_average_observerRF
    plot_average_routines['transmission_lightcurve_average_stellarRF'] = SLOPpy.plot_transmission_lightcurve_average_stellarRF
    plot_average_routines['transmission_lightcurve_average'] = SLOPpy.plot_transmission_lightcurve_average

    plot_average_routines['transmission_clv_rm_average_planetRF'] = SLOPpy.plot_transmission_clv_rm_average_planetRF
    plot_average_routines['transmission_clv_rm_average_observerRF'] = SLOPpy.plot_transmission_clv_rm_average_observerRF
    plot_average_routines['transmission_clv_rm_average_stellarRF'] = SLOPpy.plot_transmission_clv_rm_average_stellarRF
    plot_average_routines['transmission_clv_rm_average'] = SLOPpy.plot_transmission_clv_rm_average

    plot_average_routines['compare_clv_rm_effects_planetRF'] = SLOPpy.plot_compare_clv_rm_effects_planetRF
    plot_average_routines['compare_clv_rm_effects_observerRF'] = SLOPpy.plot_compare_clv_rm_effects_observerRF
    plot_average_routines['compare_clv_rm_effects_stellarRF'] = SLOPpy.plot_compare_clv_rm_effects_stellarRF
    plot_average_routines['compare_clv_rm_effects'] = SLOPpy.plot_compare_clv_rm_effects

    plot_average_routines['transmission_map_average'] = SLOPpy.plot_transmission_map_average
    plot_average_routines['transmission_clv_rm_map_average'] = SLOPpy.plot_transmission_clv_rm_map_average

    """

    """
    Execution of subroutines
    """

    # ! NEW !
    print()
    print("*** Data preparation analysis ***")

    try:
        pipeline = config_in['pipeline']
        has_plots = len(pipeline)
    except (KeyError, TypeError):
        pipeline = {}


    for key in pipeline:
        if key in pipeline_common_routines:
            print()
            pipeline_common_routines[key](config_in)


    # ! Kept here for legacy purposes !
    for key in pipeline:
        if key in pipeline_routines:
            print()
            pipeline_routines[key](config_in)

    # ! NEW !
    print()
    print("*** Spectral lines analysis  ***")


    try:
        spectral_lines = config_in['spectral_lines']
        has_plots = len(spectral_lines)
    except (KeyError, TypeError):
        pipeline = {}

    for lines_label in spectral_lines:
        for key in pipeline:
            if key in pipeline_lines_routines:
                print()
                pipeline_lines_routines[key](config_in, lines_label)

    #for key, func in pipeline_routines.items():
    #    if key in pipeline: func(config_in)

    #TODO: must be updated to be performed on a single set of spectral lines

    try:
        plots = config_in['plots']
        has_plots = len(plots)
    except (KeyError, TypeError):
        return

    print()
    print("*** Plot Subroutines ***")
    print()

    plots = config_in['plots']
    nights = config_in['nights']

    for key in plots:
        if key in plot_preparation_routines:
            plot_preparation_routines[key](config_in)
            print()

    for key in plots:
        if key in plot_routines:
            plot_routines[key](config_in)
            print()

    for lines_label in config_in['spectral_lines']:

        for key in plots:

            for night in nights:
                    if key in plot_lines_routines:
                        plot_lines_routines[key](config_in, lines_label, night)
                        print()

            if key in plot_lines_average_routines:
                plot_lines_average_routines[key](config_in, lines_label)
                print()

        #for key, func in plot_preparation_routines.items():
        #    if key in plots: func(config_in)
        #
        #for night in nights:
        #    for key, func in plot_routines.items():
        #        if key in plots: func(config_in, night)
        #
        #for key, func in plot_average_routines.items():
        #    if key in plots: func(config_in)
