import numpy as np
import src.catalog as catalog
import src.config as config
import src.fig as fig
import src.lin as lin
import src.sys_split as sys_split
import src.corr as corr
import src.field as field
import src.pz as pz
import src.cosmosis as cosmosis
import sys

# Shear tests

def summary_tests(cat,mask):
  lin.summary_stats.i3_flags_dist(cat,mask=mask)
  lin.summary_stats.e1_psf_stats(cat,mask=mask)
  return

def linear_tests(cat,cols,flagcols,mask):
  lin.hist.hist_tests(cols,cat,mask=mask)
  lin.hist.hist_2D_tests(cols,cols,cat,mask=mask)
  lin.footprint.hexbin_tests(cols,cat,mask=mask)
  lin.footprint.footprint_tests(cat,[],mask=mask,bins=100,label='All')
  lin.footprint.footprint_tests(cat,flagcols,mask=mask,bins=100,label='')
  return

def tile_tests(cat,vals,mask):
  lin.summary_stats.tile_stats(cat,vals,mask=mask)
  lin.hist.tile_tests(vals,cat)
  lin.hist.tile_tests_2D(vals,vals,cat)
  lin.footprint.tile_tests(vals,cat,mask=mask)
  return

def field_tests(epochcat,mask):
  field.field.footprint(epochcat,mask=mask)
  field.field.whisker(epochcat,mask=mask)
  field.field.whisker_chip(i3epoch,mask=epochmask)
  return

def split_tests(cat,lenscat,cols,mask):
  sys_split.split.cat_splits_lin(cols,cat,mask=mask)
  sys_split.split.cat_splits_2pt(cols,cat,lenscat,mask=mask)
  return

def corr_tests(cat,mask):
  corr.xi_2pt.xi_2pt(cat,corr='GG',maska=mask,plot=True)
  corr.xi_2pt.xi_2pt(cat,corr='GG',ga='psf',maska=mask,plot=True)
  corr.xi_2pt.xi_2pt(cat,cat,corr='GG',gb='psf',maska=mask,plot=True)
  return

def single_tests(epochcat,epochmask,cat,mask,lenscat):
  vals=['e1','e2','psf1','psf2','psffwhm']
  summary_tests(cat,mask)
  tile_tests(cat,vals,mask)
  cols=['psffwhm','rgp','rad','psf1','psf2','snr','e1','e2','evals','iter','ra_off','dec_off','flux','chi2pix','invfluxfrac','resmin','resmax','pos','psfpos','airmass','fwhm','maglimit','skybrite','skysigma','dpsf']
  flagcols=['info','error']
  linear_tests(cat,cols,flagcols,mask)
  cols=['psffwhm','rgp','rad','psf1','psf2','snr','e1','e2','evals','iter','ra_off','dec_off','flux','chi2pix','invfluxfrac','resmin','resmax','psfpos','airmass','fwhm','maglimit','skybrite','skysigma','dpsf']
  split_tests(cat,lenscat,cols,mask)
  field_tests(epochcat,epochmask)
  corr_tests(cat,mask)
  return

# def pair_tests(cat,cat2,mask,mask2,match_col):
#   vals=['e1','e2','psf1','psf2','psffwhm']
#   lin.hist.tile_tests_2D(vals,vals,cat,cat2=cat2)
#   cols=['psffwhm','rgp','rad','psf1','psf2','snr','e1','e2','evals','iter','ra_off','dec_off','flux','chi2pix','invfluxfrac','resmin','resmax','psfpos','airmass','fwhm','maglimit','skybrite','skysigma','dpsf']
#   cols=['psffwhm','rgp','skysigma','dpsf']
#   lin.hist.hist_2D_tests(cols,cols,cat,cat2=cat2,mask=mask,mask2=mask2,match_col=match_col)
#   return


if __name__ == '__main__':


  ###   Shear testing Examples    ###


  #Load shear catalog

  #Select a subset of tiles (useful for testing). Set tiles to None for all tiles.
  tiles=['DES0428-5205', 'DES0428-5622', 'DES0428-5914', 'DES0429-4414','DES0429-4457', 'DES0429-5248', 'DES0429-5705', 'DES0430-5331','DES0430-5414', 'DES0430-5957', 'DES0431-5457', 'DES0431-5748','DES0431-6039', 'DES0432-4540', 'DES0432-4623', 'DES0432-4706','DES0432-4748', 'DES0432-4831', 'DES0432-4914', 'DES0432-4957']

  #Select a set of columns to read in from catalog files. A dictionary in config.py translates these shortnames to the actual column names in the files.
  cols=['modelmu','bamp','cov12','cov11','modelsig','stamp','damp','nexp','chi2pix','psffwhm','coadd','info','like','cov22','bflux','mu','flux','rsnr','rgp','dec','evals','rad','dec_off','ra_off','dflux','fluxfrac','psf1','psf2','hsmpsf1','hsmpsf2','modmax','modmin','sig','ra','resmax','tile','maskfrac','error','snr','resmin','e1','e2','flagi','flagr','iter','pz','r','g','i','z']

  #Read in files and build catalog class. See catalog.py for further details.
  i3=catalog.CatalogStore('y1_i3_sv_v1',cutfunc=None,cattype='i3',cols=cols,catdir='/share/des/disc2/y1/im3shape/single_band/r/y1v1/tmp/',release='y1',tiles=tiles)

  #Load in systematics maps and match to galaxy positions.
  sys_split.split_methods.load_maps(i3)


  # Select columns to read in from epoch catalog files. A dictionary in config.py translates these shortnames to the actual column names in the files.
  cols=['coadd','expnum','xoff','yoff','psf1_exp','psf2_exp','ccd','row','col','e1','e2']

  #Read in files and build catalog class. See catalog.py for further details.
  i3epoch=catalog.CatalogStore('y1_i3_sv_epoch_v1',cutfunc=None,cattype='i3',cols=cols,catdir='/share/des/disc2/y1/im3shape/single_band/r/y1v1/spte_sv/epoch_cats/',release='y1',tiles=tiles)

  #Match epoch entries to main catalog entries.
  epochmask=np.in1d(i3epoch.coadd,i3.coadd[i3.info==0],assume_unique=False)

  #Read in lens catalog for tangential shear.
  rm10=catalog.CatalogStore('y1_rm10',cutfunc=catalog.CatalogMethods.default_rm_cuts(),cattype='gal',cols=['coadd','ra','dec','zp'],catfile=config.redmagicdir+'sva1_gold_1.0.2_run_redmapper_v6.3.3_redmagic_0.5-10.fit',ranfile=config.redmagicdir+'random.npy',release='y1')

  #Run full set of shear tests on shear catalog. Separate tests are available for comparing two catalogs, but untested in this version - commented out above.
  single_tests(i3epoch,epochmask,i3,i3.info==0,rm10)


  ###   Photo-z testing Examples    ###

  # Load photo-z catalog of pdfs - h5 version

  sn=catalog.PZStore('skynet',setup=True,pztype='SKYNET',filetype='h5',file='sv_skynet_final.h5')

  # Build photo-z bins from PZStore object

  pz.pz_methods.build_nofz_bins(sn,0.3,1.3,cat=None,bins=3,split='mean',pzmask=None,catmask=None)

  # Load photo-z n(z)'s and bootstrap samples from standard dict file used for spec validation.
  
  sn=catalog.PZStore('skynet',setup=True,pztype='SKYNET',filetype='dict',file='WL_test_3_bins.pickle')

  # Plot comparison of photo-z vs spec.

  fig.plot_methods.plot_nofz(sn,'pz_test')

  # Write cosmosis n(z) files from PZStore object.

  cosmosis.make.nofz(sn,'pz_test')

  # Submit cosmosis runs for spec validation

  cosmosis.run.submit_pz_spec_test(sn,'pz_test',boot=False,cosmo=False)
  cosmosis.run.submit_pz_spec_test(sn,'pz_test',boot=True,cosmo=False)
  # First two jobs must finish writing simulated data before running second two commands
  cosmosis.run.submit_pz_spec_test(sn,'pz_test',boot=False,cosmo=True)
  cosmosis.run.submit_pz_spec_test(sn,'pz_test',boot=True,cosmo=True)

  # When cosmosis jobs finished, calculate output figures and stats

  fig.plot_methods.plot_pz_sig8('pz_test',boot=True,tomobins=3)

  