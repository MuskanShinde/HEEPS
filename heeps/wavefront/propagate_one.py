from heeps.optics import apodizer, fp_mask, lyot_stop, detector
from heeps.util.img_processing import resize_img, pad_img
from copy import deepcopy
import numpy as np
import proper

def propagate_one(wf, conf, phase_screen=None, misalign=None, zernike=None, verbose=False):
            
    """ Create a function to propagate one single wavefront """

    # keep a copy of the input wavefront
    wf1 = deepcopy(wf)

    # apply phase screen (scao residuals, ncpa, petal piston)
    if phase_screen is not None:
        phase_screen = np.nan_to_num(phase_screen)# in microns
        phase_screen = resize_img(phase_screen, conf['npupil'])
        phase_screen = pad_img(phase_screen, conf['ngrid'])        
        proper.prop_add_phase(wf1, phase_screen)

    # apply tip-tilt (Zernike 2,3)
    if zernike is not None:
        # convert the tip/tilt (mas) to RMS phase errors (m)
        # theta ~ tan(theta) = (2*RMS) / pupil radius => RMS = theta*D/4
        zernike = np.array(zernike, ndmin=1)*u.mas.to('rad')*conf['diam_ext']/4
        zern_inds = np.array(zern_inds, ndmin=1)
        assert zern_inds.size == zernike.size, \
            "Zernike values and indices must be arrays/lists of same length."
        proper.prop_zernikes(wf1, zern_inds, zernike)

    # RA misalignment
    conf['RAVC_misalign'] = misalign

    # pupil-plane apodization
    wf1, apo_amp, apo_phase = apodizer(wf1, get_amp=True, verbose=verbose, **conf)
    # focal-plane mask
    wf1 = fp_mask(wf1, conf, verbose=verbose)
    # Lyot-stop or APP
    wf1, ls_amp, ls_phase = lyot_stop(wf1, get_amp=True, verbose=verbose, **conf)
    # detector
    psf = detector(wf1, conf, verbose=verbose)

    return psf
