import numpy as np
import matplotlib.pyplot as plt
import healpy as hp

#mpz_nomask=[]
#mpz_wmask=[]
xsc_nomask=[]
xsc_wmask=[]

#mpz_pattern="contrast_2mpz_nside32_5deg_n0p05_"
xsc_pattern="../masking/contrast_xsc_nside32_5deg_n0p05_"

for i in range(4):
    #mpz_nmaski=hp.read_map("{0}nomask_band{1}.fits".format(mpz_pattern, i+1))    
    #mpz_wmaski=hp.read_map("{0}wmask2_band{1}.fits".format(mpz_pattern, i+1))
    xsc_nmaski=hp.read_map("{0}nomask_band{1}.fits".format(xsc_pattern, i+1))
    xsc_wmaski=hp.read_map("{0}wmask_band{1}.fits".format(xsc_pattern, i+1))

    #mpz_nomask.append(mpz_nmaski)
    #mpz_wmask.append(mpz_wmaski)
    xsc_nomask.append(xsc_nmaski)
    xsc_wmask.append(xsc_wmaski)

#wmap_ilc=hp.read_map("wmap_ilc_9yr_v5_beam4p9degs_ns32_2muK.fits")
wmapQ=hp.read_map("../masking/wmap_band_forered_imap_r9_9yr_Q_v5_beam5deg_ns32_2muK.fits")
wmapV=hp.read_map("../masking/wmap_band_forered_imap_r9_9yr_V_v5_beam5deg_ns32_2muK.fits")
wmapW=hp.read_map("../masking/wmap_band_forered_imap_r9_9yr_W_v5_beam5deg_ns32_2muK.fits")

wmapQ_wmask=hp.read_map("../masking/wmap_band_forered_imap_r9_9yr_Q_v5_beam5deg_ns32_2muK_wmask.fits")
wmapV_wmask=hp.read_map("../masking/wmap_band_forered_imap_r9_9yr_V_v5_beam5deg_ns32_2muK_wmask.fits")
wmapW_wmask=hp.read_map("../masking/wmap_band_forered_imap_r9_9yr_W_v5_beam5deg_ns32_2muK_wmask.fits")

wmap_nomask=[wmapQ, wmapV, wmapW]
wmap_wmask=[wmapQ_wmask, wmapV_wmask, wmapW_wmask]
titles_wmap=["WMAP (Q channel)", "WMAP (V channel)", "WMAP (W channel)"]
fnames_wmap=["wmap_Q", "wmap_V", "wmap_W"]

mask=hp.read_map("mask_combined_xsc_kq85_ns32.fits")
plt.figure()
hp.mollview(mask, title="Mask", cmap='RdYlBu_r')
plt.savefig("mask.png")
plt.close()

for i in range(4):
    ''' template for 2mpz if needed
    plt.figure()
    hp.mollview(mpz_nomask[i], title="2MASS contrast - band {} with no mask (2mpz)".format(i+1), min=-1, max=1, cmap='jet')
    hp.graticule()
    plt.savefig("band{}_nomask_2mpz.png".format(i+1))
    plt.close()
    
    plt.figure()
    hp.mollview(mpz_wmask[i], title="2MASS contrast - band {} with mask (2mpz)".format(i+1), min=-1, max=1, cmap='jet')
    hp.graticule()
    plt.savefig("band{}_wmask_2mpz.png".format(i+1))
    plt.close()
    '''

    plt.figure()
    hp.mollview(xsc_nomask[i], title="2MASS contrast - band {} with no mask (xsc)".format(i+1), min=-1, max=1, cmap='RdYlBu_r')
    hp.graticule()
    plt.savefig("band{}_nomask_xsc.png".format(i+1))
    plt.close()

    plt.figure()
    hp.mollview(xsc_wmask[i], title="2MASS contrast - band {} with mask (xsc)".format(i+1), min=-1, max=1, cmap='RdYlBu_r')
    hp.graticule()
    plt.savefig("band{}_wmask_xsc.png".format(i+1))
    plt.close()

    if i<3:
        plt.figure()
        hp.mollview(wmap_nomask[i], title=titles_wmap[i], unit="mK", cmap='RdYlBu_r')#min=-0.15, max=0.15
        hp.graticule()
        plt.savefig("{}_nomask.png".format(fnames_wmap[i]))
        plt.close()

        plt.figure()
        hp.mollview(wmap_wmask[i], title=titles_wmap[i], unit="mK", cmap='RdYlBu_r')
        hp.graticule()
        plt.savefig("{}_wmask.png".format(fnames_wmap[i]))
        plt.close()
