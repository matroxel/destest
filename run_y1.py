import numpy as np
import os
import fitsio as fio
import healpy as hp
import src.catalog as catalog
import src.config as config
import src.fig as fig
import src.lin as lin
import src.sys_split as sys_split
import src.corr as corr
import src.field as field
import src.pz as pz
import src.cosmo as cosmo
import src.y1shearcat as y1
import src.txt as txt

#400 seconds to read mcal
#250 seconds to read i3

NSIDE_DEFAULT = 4096


goldfile  = '/global/cscratch1/sd/troxel/finaly1cats/y1a1-gold-mof-badregion.fits'
i3file    = '/global/cscratch1/sd/troxel/finaly1cats/y1a1-im3shape_v5-unblind_v2_matched_v3.fits'
mcalfile  = '/global/cscratch1/sd/troxel/finaly1cats/mcal-y1a1-combined-riz-unblind-v4-matched.fits'
i3pickle  = '/global/cscratch1/sd/troxel/finaly1cats/i3.cpickle'
mcalpickle= '/global/cscratch1/sd/troxel/finaly1cats/mcal.cpickle'
bpzfile   = '/global/cscratch1/sd/troxel/finaly1cats/y1a1-gold-mof-badregion_BPZ.fits'
bpz0file  = '/global/cscratch1/sd/troxel/finaly1cats/mcal-y1a1-combined-griz-blind-v3-matched_BPZ.fits'
i3epochdir= '/project/projectdirs/des/wl/desdata/wlpipe/im3shape_y1a1_v5/bord/epoch/'
mcalepochdir = '/global/cscratch1/sd/troxel/finaly1cats/'
i3epochpickle  = '/global/cscratch1/sd/troxel/finaly1cats/i3epoch.cpickle'
mcalepochpickle= '/global/cscratch1/sd/troxel/finaly1cats/mcalepoch.cpickle'
rmfile    = '/global/cscratch1/sd/troxel/redmagicv6.4.11/5bins_hidens_hilum_higherlum_jointmask_0.15-0.9_magauto_mof_combo_removedupes_spt_fwhmi_exptimei_cut_badpix_sample_weighted2sig_sorted.fits.gz'
psfdir    = '/global/cscratch1/sd/troxel/finaly1cats/psf_cats/'
psfexpinf = '/global/cscratch1/sd/troxel/finaly1cats/exposure_info_y1a1-v02.fits'
imagefile = '/global/cscratch1/sd/troxel/finaly1cats/y1a1_image_id.fits'
special_points_file = '/global/cscratch1/sd/zuntz/y1a1_special_field_points.fits'

psfpickle = '/global/cscratch1/sd/troxel/finaly1cats/psf.cpickle'
rm_maskfile = '/global/cscratch1/sd/troxel/finaly1cats/5bins_hidens_hilum_higherlum_jointmask_0.15-0.9_magauto_mof_combo_removedupes_spt_fwhmi_exptimei_cut_mask.fits.gz'
rm_randomsfile = '/global/cscratch1/sd/troxel/finaly1cats/5bins_hidens_hilum_higherlum_jointmask_0.15-0.9_magauto_mof_combo_removedupes_spt_fwhmi_exptimei_cut_randoms.fits.gz'

def build_special(special_points_file, rm_maskfile, rm_randomsfile):
    special=catalog.CatalogStore("special", catfile=special_points_file, cattype='gal', cols=-1, 
        cutfunc=catalog.CatalogMethods.final_null_cuts_ra(), ranfile=rm_randomsfile)
    special.add_pixel(4096, nest=False)

    rm_mask_hpix = fio.FITS(rm_maskfile)[1]['HPIX'][:]
    inmask=np.in1d(special.pix, rm_mask_hpix)
    catalog.CatalogMethods.match_cat(special,inmask)
    return special

special = build_special(special_points_file, rm_maskfile, rm_randomsfile)


# Load various data - Probably don't do this all in same job, lots of memory if you're not loading the pickles, and even then...
i3,mcal  = y1.y1.load_data(i3file,mcalfile,goldfile,bpzfile,bpz0file,i3pickle=i3pickle,mcalpickle=mcalpickle)

print "(Not loading PSF)"
#psf = y1.y1.load_psf_data(psfdir,psfpickle)
mcalepoch,i3epoch = y1.y1.load_epoch_data(i3epochdir,mcalepochdir,imagefile,i3pickle,mcalpickle,i3epochpickle,mcalepochpickle)

# rm = catalog.CatalogStore('rm',cutfunc=catalog.CatalogMethods.final_null_cuts_ra(),cattype='gal',catfile=rmfile,cols=['coadd','ra','dec'])
# special=catalog.CatalogStore("special", catfile=special_points_file, cattype='gal', cols=-1, cutfunc=catalog.CatalogMethods.final_null_cuts_ra())


# def add_i3_cal_to_epoch(epoch, i3):
#     x=np.in1d(epoch.coadd,i3.coadd,assume_unique=False)
#     catalog.CatalogMethods.match_cat(epoch,x)
#     mask=np.argsort(epoch.coadd)
#     catalog.CatalogMethods.match_cat(epoch,mask)

#     diff=np.diff(epoch.coadd)
#     diff=np.where(diff!=0)[0]+1
#     # diff = np.concatenate([[0], diff, [None]])
#     diff=np.append([0],diff)
#     diff=np.append(diff,[None])

#     epoch.m1=np.zeros(len(epoch.coadd))
#     epoch.m2=np.zeros(len(epoch.coadd))
#     epoch.c1=np.zeros(len(epoch.coadd))
#     epoch.c2=np.zeros(len(epoch.coadd))
#     epoch.w=np.zeros(len(epoch.coadd))

#     for i in xrange(len(diff)-1):
#         epoch.m1[diff[i]:diff[i+1]]=i3.m1[i]
#         epoch.m2[diff[i]:diff[i+1]]=i3.m2[i]
#         epoch.c1[diff[i]:diff[i+1]]=i3.c1[i]
#         epoch.c2[diff[i]:diff[i+1]]=i3.c2[i]
#         epoch.w[diff[i]:diff[i+1]]=i3.w[i]
    

# def make_randoms(name, cat, rans_per_object):
#     if os.path.exists(name+"_random.fits.gz"):
#         print "Using existing random", name+"_random.fits.gz"
#     else:
#         print "Making new random (one off)"
#         catalog.CatalogMethods.create_random_cat_from_cat(cat, nran=len(cat.e1)*rans_per_object, label=name+"_")
#     f = fitsio.FITS(name+"_random.fits.gz")
#     data = f[1].read_columns(['ra', 'dec'])
#     cat.ran_ra = data['ra']
#     cat.ran_dec = data['dec']
#     f.close()



# add_randoms(i3, 1)
# add_randoms(metacal, 1)


#txt.write_methods.heading('---------------',mcal,label='y1_paper',create=True)



# Making all the plots for the paper (should probably not try to run all these at once in same job)

y1.y1_plots.tangential_shear_plot(i3,mcal,special)


# Normal catalog plots
y1.y1_plots.mean_e(i3,mcal)
y1.y1_plots.footprint_plot(mcal)


# single-epoch plots
y1.y1_plots.mean_e_row_plot(mcalepoch)
y1.y1_plots.mean_e_fov(mcalepoch,i3epoch)

print "PSF plots disabled for now"
import sys
sys.exit(1)
# psf catalog plots - from Mike's psfex files
y1.y1_plots.psf_star_dist(psf)
y1.y1_plots.psf_mag_res_plot(psf)
y1.y1_plots.psf_star_fwhm_dist(psf,psfexpinf)
y1.y1_plots.psf_e_fov(psf)

