from heeps.util.download_from_gdrive import extract_zip
import astropy.units as u
import os
import numpy as np
import proper
proper.print_it = False


def read_config(verbose=False, **update_conf):

    conf = dict(

    # =============================================================================
    #           Console and file management 
    # =============================================================================
    cpu_count = 10,                     # 1 = single core; None = max-1 cores
    send_to = None,                     # user's email, for notifications
    prefix = '',                        # for saved files: e.g. 'test_'
    headless = False,                   # true if running on a headless server
    gdrive_id = '1wj3onWQ9GVW-l8X58JMgAj-9TNqalKb-', # Google Drive ID
    # required directories for data (e.g. fits files)
    dir_current = '$HOME/heeps_metis',
    dir_input = 'input_files',
    dir_output = 'output_files',
    dir_temp = 'temp_files',

    # =============================================================================
    #           Parameters for telescope
    # =============================================================================
    focal = 658.6,                      # focal distance in m
    pupil_img_size = 39.9988,           # pupil image in m (for PROPER)
    diam_nominal = 38.542,              # nominal diameter (for LS oversizing)
    diam_ext = 36.905,                  # effective outer circular aperture in m
    diam_int = 11.213,                  # effective central obscuration in m
    file_pupil = 'pupil/ELT_allglass.fits',# entrance pupil file
    spi_width = 0.54,                   # spider width in m
    spi_angles = [0, 60, 120],          # spider angles in deg
    # if no valid pupil file, pupil will be created with the following params:
    seg_width = 1.45,                   # segment width in m
    seg_gap = 0.004,                    # gap between segments in m
    seg_rms = 0,                        # rms of the reflectivity of all segments
    select_petal = None,
    npetals = 6,                        # number of petals
    # number of hexagonal segments per column (from left to right)
    seg_ny = np.array([10, 13, 16, 19, 22, 23, 24, 25, 26, 27, 28, 29, \
                      30, 31, 30, 31, 30, 31, 30, 31, 30, 31, 30, 29, \
                      28, 27, 26, 25, 24, 23, 22, 19, 16, 13, 10]),
    # coordinates of missing segments
    seg_missing = [],
#    seg_missing = [(-2,5),(-2,6),(-3,4),(-3,5),(-3,6),(-4,5),(-4,6)],


    # =============================================================================
    #           Parameters for observing modes and bands
    # =============================================================================
    # Types of various modes include (e.g. for METIS):
    #    1. mode = 'RAVC' for Ring Apodized Vortex Coronagraph
    #    2. mode = 'CVC'  for Classical Vortex Coronagraph
    #    3. mode = 'APP'  for Apodizing Phase Plate
    #    4. mode = 'CLC'  for Classical Lyot Coronagraph
    #    5. mode = 'ELT'  for no coronagraph (only telescope)
    # Default mode: L-band Ring Apodized Vortex
    mode = 'RAVC',                      # HCI mode
    band = 'L',                         # spectral band
    lam = 3.81e-6,                      # wavelength in m
    pscale = 5.47,                      # pixel scale in mas/pix
    ngrid = 1024,                       # number of pixels of the wavefront array
    npupil = 285,                       # number of pixels of the pupil
    ndet = 365,                         # size of the detector plane array
    hfov = 1.1,                         # (optional) half FOV in arcsec (updates ndet)
    add_bckg = False,                   # true means background flux and photon noise are added 
    mag = 5,                            # star magnitude at selected band
    mag_ref = 0,                        # reference magnitude for star and background fluxes
    flux_star = 8.999e+10,              # [e-/s] HCI-L long, mag 0 (Jan 21, 2020)
    flux_bckg = 8.878e+04,              # [e-/s/pix]
    cube_duration = 3600,               # cube duration in seconds
    lat = -24.59,                       # telescope latitude in deg (Armazones=-24.59 ,Paranal -24.63)
    dec = -5,                           # star declination in deg (e.g. 51 Eri -2.47)
    file_lyot_stop = 'pupil/ls_ravc_allglass_285.fits', # lyot stop file
    ls_dRext = 0.0282,                  # LS Rext undersize (% diam ext)
    ls_dRint = 0.0282,                  # LS Rint oversize (% diam ext)
    ls_dRspi = 0.037,                   # LS spider oversize (% diam ext)
    ls_misalign = [0,0,0,0,0,0],        # Lyot stop misalignment
    vc_charge = 2,                      # vortex topological charge
    vc_zoffset = 0,                     # vortex defocus in m (z axis)
    vc_chrom_leak = 2e-3                # vortex chromatic leakage
    ravc_calc = True,                   # calculate RA params (Mawet2013)
    ravc_t = 0.7903,                    # (calc=False) mean-M1 RA trans
    ravc_r = 0.6033,                    # (calc=False) mean-M1 RA radius wrt allglass
    ravc_misalign = [0,0,0,0,0,0],      # RA misalignment
    clc_diam = 80,                      # CLC occulter diam in mas
    file_vc_trans = 'optics/agpm_trans.fits', # vortex transmittance
    file_app_trans = 'optics/metis_gvapp_tx.fits', # APP transmittance
    file_app_amp = 'optics/APP_stop_L_285_v2.fits', # APP amplitude
    file_app_phase = 'optics/vAPP_Dshape_Lband_asymmetric.fits', # APP phase
    app_strehl = 0.64,                   # APP Strehl ratio
    app_single_psf = 0.48,               # APP single PSF (4% leakage)
    student_distrib = True,              # use Student's distribution instead of Gaussian
    # Multiple spectral bands
    bands = ['L', 'M', 'N1', 'N2'],
    band_specs = {  
        'L': {'lam': 3.82e-6,
            'pscale': 5.47,
            'flux_star': 8.999e+10,                 # HCI-L long
            'flux_bckg': 8.878e+04,
            'ls_dRspi': 0.037},
        'M': {'lam': 4.8e-6,
            'pscale': 5.47,
            'flux_star': 2.452e+10,                 # CO ref
            'flux_bckg': 6.714e+05,
            'ls_dRspi': 0.037},
        'N1': {'lam': 8.65e-6,
            'pscale': 6.79,
            'flux_star': 3.684e+10,                 # GeoSnap N1
            'flux_bckg': 4.725e+07,
            'ls_dRspi': 0.037},
        'N2': {'lam': 11.25e-6, 
            'pscale': 6.79,
            'flux_star': 3.695e+10,                 # GeoSnap N2
            'flux_bckg': 1.122e+08,
            'ls_dRspi': 0.037},
        'N1a': {'lam': 8.65e-6, 
            'pscale': 10.78,
            'flux_star': 2.979e+10,                 # Aquarius N1
            'flux_bckg': 9.630e+07,
            'ls_dRspi': 0.037},
        'N2a': {'lam': 11.25e-6, 
            'pscale': 10.78,
            'flux_star': 2.823e+10,                 # Aquarius N2
            'flux_bckg': 2.142e+08,
            'ls_dRspi': 0.037}
        },
    # Multiple HCI modes
    modes = ['RAVC', 'CVC', 'CLC', 'APP', 'ELT'],
    mode_specs = {
        'RAVC': {'ls_dRint': 0.0282},
        'CVC': {'ls_dRint': 0.05},
        'CLC': {'ls_dRint': 0.10}
        },

    # =============================================================================
    #           Parameters for wavefront
    # ============================================================================
    nframes = 10,                       # number of frames to crop the input data
    nstep = 1,                          # take 1 frame every nstep (cubesize = nframes/nstep)

    add_phase = True,                   # phase screens (SCAO residuals, NCPA, petal piston)
    file_phase = 'wavefront/COMPASS_201810_RandomWind_100screens_meters.fits',
    add_amp = False,                    # amplitude screens (Talbot effect)
    file_amp = 'wavefront/Talbot_LM_20201120_IMGP_meridian_allglass.fits',

    rms_phase_sta = 35.9,               # static (nm)
    rms_phase_qlsf = 20,                # quasistatic low spatial freq (nm)
    rms_phase_qhsf = 20,                # quasistatic high spatial freq (nm)
    rms_phase_dyn = 40,                 # dynamic (nm)

    add_point_err = False,              # pointing errors
    file_point_err = 'wavefront/point_all_3600s_300ms.fits',
    rms_point_qsta = 0.4,               # quasistatic (mas)
    rms_point_dyn = 2,                  # dynamic (mas)

    add_apo_drift = False,              # apodizer drift
    ptv_drift = 0.02,                   # (%)
    
    add_vort_chrom_leak = False,        # add chromatic leakage in the vortex plane
    add_det_chrom_leak = False          # add chromatic leakage in the detector plane

    )                                   # end of default conf dict
 
    # =============================================================================
    #           Perform initialization actions
    # =============================================================================
    
    # update conf dictionary
    conf.update(**update_conf)    
    if verbose is True:
        print('Default config: band=%s, mode=%s'%(conf['band'], conf['mode']))
        print('\u203e'*15)
        print('   npupil=%s, pscale=%s mas, lam=%3.4E m'\
            %(conf['npupil'], conf['pscale'], conf['lam']))
        hfov = conf['ndet']/2*conf['pscale']/1e3
        hfov_lamD = hfov*u.arcsec.to('rad')/(conf['lam']/conf['diam_ext'])
        print('   ndet=%s (-> hfov=%s arcsec, %s lam/D)\n'%(conf['ndet'], \
            round(hfov,2), round(hfov_lamD,2)))
    
    # create directories
    conf['dir_current'] = os.path.normpath(os.path.expandvars(conf['dir_current']))
    conf['dir_input'] = os.path.join(conf['dir_current'], conf['dir_input'], '')
    conf['dir_output'] = os.path.join(conf['dir_current'], conf['dir_output'], '')
    conf['dir_temp'] = os.path.join(conf['dir_current'], conf['dir_temp'], '')
    os.makedirs(conf['dir_input'], exist_ok=True)
    os.makedirs(conf['dir_output'], exist_ok=True)
    os.makedirs(conf['dir_temp'], exist_ok=True)
    
    # create paths to fits files
    for filename in ['file_pupil', 'file_phase', 'file_amp', 'file_point_err', \
            'file_lyot_stop', 'file_vc_trans', 'file_app_trans', \
            'file_app_amp', 'file_app_phase']:
        conf[filename] = os.path.join(conf['dir_input'], conf[filename])
    
    # downloading input files from Google Drive
    verif_file = 'wavefront/COMPASS_201810_RandomWind_100screens_meters.fits'
    if not os.path.isfile(os.path.join(conf['dir_input'], verif_file)):
        print("Downloading input files from Google Drive to \n'%s'\n"%conf['dir_input'])
        extract_zip(conf['gdrive_id'], conf['dir_input'])
    
    # disable matplotlib display to run on a headless server
    if conf['headless'] is True:
        import matplotlib; matplotlib.use('agg')

    # sort alphabetically
    conf = {k: v for k, v in sorted(conf.items())}

    return conf
