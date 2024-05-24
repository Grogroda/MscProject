#!/usr/bin/bash

echo "Usage: ./apply_contrast_mask.sh [band (1,2,3 or 4)]"

prefix=contrast_xsc_nside32_5deg_n0p05_
echo "Band: $1"

/home/arthurdmeirelles/powertools/v1.6/bin/apply_mask ${prefix}nomask_band${1}.fits mask_combined_xsc_kq85_ns32.fits ${prefix}wmask_band${1}.fits
