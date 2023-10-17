import numpy as np
import matplotlib.pyplot as plt
import healpy as hp

mpz_nomask=[]
mpz_wmask=[]
xsc_nomask=[]
xsc_wmask=[]

mpz_pattern="contrast_2mpz_nside32_5deg_n0p05_"
xsc_pattern="contrast_xsc_nside32_5deg_n0p05_"

for i in range(4):
    mpz_nmaski=hp.read_map("{0}nomask_band{1}.fits".format(mpz_pattern, i+1))    
    mpz_wmaski=hp.read_map("{0}wmask_band{1}.fits".format(mpz_pattern, i+1))
    xsc_nmaski=hp.read_map("{0}nomask_band{1}.fits".format(xsc_pattern, i+1))
    xsc_wmaski=hp.read_map("{0}wmask_band{1}.fits".format(xsc_pattern, i+1))

    mpz_nomask.append(mpz_nmaski)
    mpz_wmask.append(mpz_wmaski)
    xsc_nomask.append(xsc_nmaski)
    xsc_wmask.append(xsc_wmaski)

wmap_ilc=hp.read_map("wmap_ilc_9yr_v5_beam4p9degs_ns32_2muK.fits")
wmapQ=hp.read_map("wmap_band_forered_imap_r9_9yr_Q_v5_beam5deg_ns32_2muK.fits")
wmapV=hp.read_map("wmap_band_forered_imap_r9_9yr_V_v5_beam5deg_ns32_2muK.fits")
wmapW=hp.read_map("wmap_band_forered_imap_r9_9yr_W_v5_beam5deg_ns32_2muK.fits")

wmap_nomask=[wmap_ilc, wmapQ, wmapV, wmapW]
titles_wmap=["Ctt WMAP (ilc)", "Ctt WMAP (Q channel)", "Ctt WMAP (V channel", "Ctt WMAP (W channel)"]
fnames_wmap=["wmap_ilc.png", "wmap_Q.png", "wmap_V.png", "wmap_W.png"]

mask=hp.read_map("mask_combined_xsc_kq85_ns32.fits")
plt.figure()
hp.mollview(mask, title="Mask", cmap="jet")
plt.savefig("mask.png")
plt.close()

for i in range(4):
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

    plt.figure()
    hp.mollview(xsc_nomask[i], title="2MASS contrast - band {} with no mask (xsc)".format(i+1), min=-1, max=1, cmap='jet')
    hp.graticule()
    plt.savefig("band{}_nomask_xsc.png".format(i+1))
    plt.close()

    plt.figure()
    hp.mollview(xsc_wmask[i], title="2MASS contrast - band {} with mask (xsc)".format(i+1), min=-1, max=1, cmap='jet')
    hp.graticule()
    plt.savefig("band{}_wmask_xsc.png".format(i+1))
    plt.close()

    plt.figure()
    hp.mollview(wmap_nomask[i], title=titles_wmap[i], unit="mK", min=-0.2, max=0.2, cmap="jet")
    hp.graticule()
    plt.savefig(fnames_wmap[i])
    plt.close()
