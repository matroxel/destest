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

def summary_tests(cat,mask):
  lin.summary_stats.i3_flags_dist(cat,mask=mask)
  lin.summary_stats.e1_psf_stats(cat,mask=mask)
  return

def linear_tests(cat,cols,flagcols,mask):
  lin.hist.hist_tests(cols,cat,mask=mask)
  lin.hist.hist_2D_tests(cols,cols,cat,mask=mask)
  lin.footprint.hexbin_tests(cols,cat,mask=mask)
  lin.footprint.footprint_tests(cat,[],mask=mask,bins=100,label='All')
  #lin.footprint.footprint_tests(cat,flagcols,mask=mask,bins=100,label='')
  return

def tile_tests(cat,vals,mask):
  lin.summary_stats.tile_stats(cat,vals,mask=mask)
  lin.hist.tile_tests(vals,cat)
  lin.hist.tile_tests_2D(vals,vals,cat)
  lin.footprint.tile_tests(vals,cat,mask=mask)
  return

def field_tests(cat,mask):
  lin.field.footprint(cat,mask=mask)
  lin.field.whisker(cat,mask=mask)
  return

def split_tests(cat,cols,mask):
  sys_split.split.cat_splits_lin(cols,cat,mask=mask)
  sys_split.split.cat_splits_2pt(cols,cat,mask=mask)
  return

def corr_tests(cat,mask):
  corr.xi_2pt.xi_2pt(cat,corr='GG',maska=mask,plot=True)
  corr.xi_2pt.xi_2pt(cat,corr='GG',ga='psf',maska=mask,plot=True)
  corr.xi_2pt.xi_2pt(cat,cat,corr='GG',gb='psf',maska=mask,plot=True)
  return

def single_tests(cat,mask):
  vals=['e1','e2','psf1','psf2','psffwhm']
  summary_tests(cat,mask)
  tile_tests(cat,mask)
  cols=['psffwhm','rgp','rad','psf1','psf2','snr','e1','e2','evals','iter','ra_off','dec_off','flux','chi2pix','invfluxfrac','resmin','resmax','pos','psfpos','airmass','exptime','fwhm','maglimit','skybrite','skysigma','dpsf']
  flagcols=['info','error']
  linear_tests(cat,cols,flagcols,mask)
  cols=['psffwhm','rgp','rad','psf1','psf2','snr','e1','e2','evals','iter','ra_off','dec_off','flux','chi2pix','invfluxfrac','resmin','resmax','psfpos','airmass','exptime','fwhm','maglimit','skybrite','skysigma','dpsf']
  split_tests(cat,cols,mask)
  corr_tests(cat,mask)
  return

def pair_tests(cat,cat2,mask,mask2,match_col):

  vals=['e1','e2','psf1','psf2','psffwhm']
  lin.hist.tile_tests_2D(vals,vals,cat,cat2=cat2)

  cols=['psffwhm','rgp','rad','psf1','psf2','snr','e1','e2','evals','iter','ra_off','dec_off','flux','chi2pix','invfluxfrac','resmin','resmax','psfpos','airmass','exptime','fwhm','maglimit','skybrite','skysigma','dpsf']
  cols=['psffwhm','rgp','skysigma','dpsf']
  lin.hist.hist_2D_tests(cols,cols,cat,cat2=cat2,mask=mask,mask2=mask2,match_col=match_col)

  return


if __name__ == '__main__':

  cols=['modelmu','bamp','cov12','cov11','modelsig','stamp','damp','nexp','chi2pix','psffwhm','coadd','info','like','cov22','bflux','mu','flux','rsnr','rgp','dec','evals','rad','dec_off','ra_off','dflux','fluxfrac','psf1','psf2','hsmpsf1','hsmpsf2','modmax','modmin','sig','ra','resmax','tile','maskfrac','error','snr','resmin','e1','e2','flagi','flagr','iter','pz','r','g','i','z']
  i3=catalog.CatalogStore('y1_i3_sv_v1',cutfunc=None,cattype='i3',cols=cols,catdir='/share/des/disc2/y1/im3shape/single_band/r/y1v1/tmp/',release='y1')

  