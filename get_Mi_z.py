import numpy as np
from SED_Model import lrt_model

def get_Mi_z(mags, mag_errs, zs, Mi_catalog=None, z_catalog=None):

    obj = lrt_model()

    obj.jy = np.zeros(mags.shape)
    obj.ejy = np.zeros(mags.shape)
    obj.jyuse = np.ones(mags.shape, dtype=np.int32)
    for j, mag in enumerate(mags):
        if np.isnan(mag_errs[j]):
            obj.jyuse[j] = 0
        else:
            obj.jy[j] = 3631.*np.exp(-0.4*mags[j])
            obj.ejy[j] = 0.4*np.log(10.)*obj.jy[j]*mag_errs[j]

    norm = 0
    if Mi_catalog is not None:   
        obj.zspec = z_catalog
        obj.kc_fit()
        norm = obj.abs_mag[0] - Mi_catalog

    Mi_z = np.zeros(zs.shape)
    for k,z in enumerate(zs):
        obj.zspec = z
        obj.kc_fit()
        Mi_z[k] = obj.abs_mag[0] - norm

    return Mi_z