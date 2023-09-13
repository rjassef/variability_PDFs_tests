import numpy as np
from astropy.table import Table
from SED_Model import lrt_model

tab = Table.read("20random_sources.dat", format='ascii.csv')

bands = ["u","g","r","i","z"]
for i in range(len(tab)):
    obj = lrt_model()
    obj.jy = np.zeros(len(bands))
    obj.ejy = np.zeros(len(bands))
    obj.jyuse = np.ones(len(bands), dtype=np.int32)
    for j, band in enumerate(bands):
        if hasattr(tab['psMagErr_{}'.format(band)],'mask') and tab['psMagErr_{}'.format(band)].mask[i]:
            obj.jyuse[j] = 0
        else:
            obj.jy[j] = 3631.*np.exp(-0.4*tab['psMag_{}'.format(band)][i])
            obj.ejy[j] = 0.4*np.log(10.)*obj.jy[j]*tab['psMagErr_{}'.format(band)][i]

    obj.ejy = np.where(obj.ejy<0.1*obj.jy, 0.1*obj.jy, obj.ejy)
    obj.chi2zop = 1
    obj.pz_fit()
