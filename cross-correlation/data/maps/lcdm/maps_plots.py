import numpy as np
import matplotlib.pyplot as plt
import healpy as hp

band1_nomask=hp.read_map("contrast_2mpz_nside32_5deg_n0p05_nomask_band1.fits")
band2_nomask=hp.read_map("contrast_2mpz_nside32_5deg_n0p05_nomask_band2.fits")
band3_nomask=hp.read_map("contrast_2mpz_nside32_5deg_n0p05_nomask_band3.fits")
band4_nomask=hp.read_map("contrast_2mpz_nside32_5deg_n0p05_nomask_band4.fits")

no_mask=[band1_nomask, band2_nomask, band3_nomask, band4_nomask]

mask=hp.read_map("mask_combined_xsc_kq85_ns32.fits")

for i in range(4):
    plt.figure()
    hp.mollview(no_mask[i], title="Band {} with no mask".format(i+1), unit="mK", min=-1, max=1, cmap='RdYlBu_r')
    hp.graticule()
    plt.savefig("band{}_nomask.png".format(i+1))
    plt.close()



