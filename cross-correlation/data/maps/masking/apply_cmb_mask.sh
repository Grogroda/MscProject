#!/usr/bin/bash

echo "Usage: ./apply_cmb_mask.sh [band (Q,V or W)]"

prefix=wmap_band_forered_imap_r9_9yr_
echo "Band: $1"
suffix=_v5_beam5deg_ns32_2muK

/home/arthurdmeirelles/powertools/v1.6/bin/apply_mask ${prefix}${1}${suffix}.fits mask_combined_xsc_kq85_ns32.fits ${prefix}${1}${suffix}_wmask.fits
