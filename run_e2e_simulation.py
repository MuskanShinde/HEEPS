
# ===================== #
# End-to-end simulation #
# ===================== #
#scp mshinde@fenrir.astro.ulg.ac.be:~/HEEPS/clc_optim/circular_lyot_stop/10_min_SCAO/cc_adi_bckg0_nominal_pup_Circular_Lyot_stop_best_L_CLC.fits E:\IISER\5th_year\Fenrir\lyot_stop_designs_for_clc_diam_85\circular_lyot_stop\contrast_curves\10_min_SCAO
import heeps
import sys

def run_e2e_simulation(savename='ADI_contrast_curve_CLC_ELT_fullM1.png'):

    # 1. Create a config dictionary with your simulation parameters in read_config
    conf = dict(
         band = 'L',
         mode = 'CVC',
         file_pupil = 'pupil/ELT_fullM1.fits',
         #file_phase = 'WFerrors/cube_Cbasic_20201130_3600s_300ms_0piston_meters_scao_only_285.fits',
         file_lyot_stop = '',
         add_phase = True,
         add_amp = True,
         add_bckg=True,
         dir_output = 'post_processing/adi',
         cpu_count = None,
         nframes = 20,
         nstep = 1,
    )
    conf = heeps.config.read_config(verbose=False, **conf)

    # 2. Update config parameters. The following parameters 
    conf = heeps.config.update_config(saveconf=True, verbose=True, **conf)
    conf['ls_dRext'] = 0.0291
    conf['ls_dRint'] = 0.08
    conf['ls_dRspi'] = 0.0317

    # 3. Load entrance pupil, and create 'wavefront' object
    wf = heeps.pupil.pupil(savefits=True, verbose=True, **conf)

    # 4. Load wavefront errors
    phase_screens, amp_screens, tiptilts, misaligns = heeps.wavefront.load_errors(verbose=True, **conf)

    # 5. Propagate one frame of offaxis psf (i.e. planet)
    psf = heeps.wavefront.propagate_one(wf, onaxis=False, savefits=True, verbose=True, **conf)

    # 6. Propagate cube of onaxis psfs (i.e. star)
    psfs = heeps.wavefront.propagate_cube(wf, phase_screens=phase_screens, \
        amp_screens=amp_screens, tiptilts=tiptilts, misaligns=misaligns, onaxis=True, \
        savefits=True, verbose=True, **conf)

    # 7. Produce 5-sigma sensitivity (contrast) curves for each set of PSFs (modes, bands)
    sep, sen = heeps.contrast.adi_one(tag='amp1', savepsf=True, savefits=True, verbose=True, **conf)

    # 8. Create a figure 
    if savename != '':
        import matplotlib.pyplot as plt
        from matplotlib.pylab import ScalarFormatter
        plt.figure(figsize=(12,4))
        plt.grid(True), plt.grid(which='minor', linestyle=':')
        plt.loglog(), plt.gca().xaxis.set_major_formatter(ScalarFormatter())
        plt.xlabel('Angular separation $[arcsec]$')
        plt.ylabel('5-$\sigma$ sensitivity (contrast)')
        label = '%s-band %s'%(conf['band'], conf['mode'])
        plt.plot(sep, sen, label=label, marker='d', markevery=0.12, markersize=4)
        plt.legend()
        plt.xticks([0.02, 0.1, 0.5, 1])
        plt.xlim(0.02,1)
        plt.ylim(1e-6,1e-2)
        plt.savefig(savename, dpi=300, transparent=True)

if __name__ == "__main__":
    '''
    Terminal command line example
    > python run_e2e_simulation.py test.png
    '''
    if len(sys.argv) > 1:
        run_e2e_simulation(sys.argv[1])
    else:
        run_e2e_simulation()