import heeps
from heeps.optics import lyot_stop, lens
from heeps.pupil import pupil
from heeps.wavefront.add_errors import add_errors
import proper
from astropy.io import fits
import os
import numpy as np

def detector(wf, phase_screen=None, amp_screen=None, tiptilt=None, misalign=[0,0,0,0,0,0], 
    ngrid=1024, ndet=365, onaxis=False, dir_output='output_files', savefits=False, 
    verbose=False, **conf):
    
    assert(ngrid >= ndet), 'Error: final image is bigger than initial grid size'
    
    # propagate to detector
    lens(wf, **conf)
    
    # add chromatic leakage
    if conf['add_det_chrom_leak'] is True and onaxis == True:
        if verbose == True:
            print('   Add chromatic leakage in detector plane')
        wf_cl = pupil(savefits=False, verbose=False, **conf)
        wf_cl = add_errors(wf_cl, phase_screen=phase_screen, amp_screen=amp_screen, \
                tiptilt=tiptilt, misalign=misalign, **conf)
        wf_cl = lyot_stop(wf_cl, verbose=False, **conf)
        lens(wf_cl, **conf)
        proper.prop_multiply(wf_cl, np.sqrt(conf['vc_chrom_leak']))
        wf._wfarr += np.transpose(wf_cl._wfarr)
        
    # get intensity (A^2)
    (psf, _) = proper.prop_end(wf, NOABS = False)
    # crop to detector size
    start = int(ngrid/2 - ndet/2) + 1
    end = int(ngrid/2 + ndet/2) + 1
    psf = psf[start:end, start:end]

    if verbose is True:
        print("   extract PSF on the detector: ndet=%s"%ndet)

    # save psf as fits file
    if savefits == True:
        os.makedirs(dir_output, exist_ok=True)
        filename = os.path.join(dir_output, 'PSF_IMG_%s.fits'%conf['band'])
        fits.writeto(filename, np.float32(psf), overwrite=True)

    return psf