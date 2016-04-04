# Select columns to read from catalog (no cols option indicates to read all columns from appropriate dict in config.py)
cols=['stamp','nexp','chi2pix','psffwhm','coadd','info','error','like','flux','rgp','dec','evals','rad','dec_off','ra_off','dflux','fluxfrac','psf1','psf2','hsmpsf1','hsmpsf2','modmax','modmin','ra','resmax','tile','maskfrac','snr','resmin','e1','e2','iter','bflux','dflux','flux','cov11','cov22','zp','g','r','i','z','bfrac']
# Build CatalogStore object for shape catalog - could include cutfunc option that points to a function to, for example, not read objects with error flags set. We want to look at the error distribution, so we don't do this here.
i3=catalog.CatalogStore('y1_i3_sv_v1',cattype='i3',cols=cols,catdir='/share/des/disc2/y1/im3shape/single_band/r/y1v1/complete/main/',release='y1',maxiter=500)
# Remove duplicate unique ids (if necessary)
catalog.CatalogMethods.remove_duplicates(i3)

# Read in galaxy catalog
rm=catalog.CatalogStore('y1_rm_highdens',cattype='gal',cols=['coadd','ra','dec','zp'],catfile=config.redmagicdir+'y1a1_gold_1.0.2b-full_redmapper_v6.4.11_redmagic_highdens_0.5-10.fit',release='y1',ranfile=config.redmagicdir+'y1a1_gold_1.0.2b-full_redmapper_v6.4.11_redmagic_highdens_0.5-10_randoms.fit')

nbc=np.load('/home/troxel/destest/i3nbcv1.npy')
a=np.argsort(nbc[:,0])
mask=np.diff(nbc[a,0])
mask=mask==0
mask=~mask
mask=a[mask]
nbc=nbc[mask]

x,y=catalog.CatalogMethods.sort2(nbc[:,0],i3.coadd)
catalog.CatalogMethods.match_cat(i3,y)
i3.m=nbc[x,1]
i3.c1=nbc[x,2]
i3.c2=nbc[x,3]
i3.w=nbc[x,4]
i3.w/=np.sum(i3.w)

# Pre-cut summary of error properties
lin.summary_stats.i3_flags_vals_check(i3)
# Plots location of flagged failures
lin.footprint.footprint_tests(i3,['info','error'])

# Cut bad objects from catalog
mask=(i3.info==0)&(i3.rgp>1.13)&(i3.snr>12)&(i3.snr<200)&(i3.rgp<3)&(~(np.isnan(i3.psf1)|np.isnan(i3.psf2)|np.isnan(i3.snr)|np.isnan(i3.psffwhm)))&(i3.g<99)&(i3.r<99)&(i3.i<99)&(i3.z<99)
catalog.CatalogMethods.match_cat(i3,mask)

# Post-cut summary of error properties
lin.summary_stats.i3_flags_vals_check(i3)

# Optionally specify sub-set of columns for the following functions
#cols=['chi2pix','psffwhm','nlike','flux','rgp','dec','evals','rad','dec_off','ra_off','dflux','invfluxfrac','psf1','psf2','hsmpsf1','hsmpsf2','modmax','modmin','ra','resmax','maskfrac','snr','resmin','e1','e2','iter','bflux','dflux','flux','zp','g','r','i','z','pos','e','psfe','dpsf','psfpos','hsmpsfe','hsmdpsf','hsmpsfpos']

# Check summary statistics of catalog values (default all)
i3.wt=False
i3.bs=False
lin.summary_stats.val_stats(i3)

# Check summary statistics of catalog values (default all)
i3.wt=True
i3.bs=True
lin.summary_stats.val_stats(i3)

# Compile summary statistics per tile of catalog values (default all)
lin.summary_stats.tile_stats(i3)

# Produce histograms of catalog columns (default all)
lin.hist.hist_tests(i3)
lin.hist.hist_2D_tests(i3)

# Produce histograms of catalog columns per tile (default all)
lin.hist.tile_tests(i3)
lin.hist.tile_tests_2D(i3)

# Produce plots of mean values of columns across survey footprint and by tile (default all)
lin.footprint.hexbin_tests(i3)
lin.footprint.tile_tests(i3)

# Produce object density plot across survey footprint (default all)
lin.footprint.footprint_tests(cat,[],label='All')

i3.wt=False
i3.bs=False
# Produce plots and statistics on linear correlation of ellipticity with catalog columns
sys_split.split.cat_splits_lin_e(i3,cols=cols)
# Produce plots and statistics on linear correlation of non-ellipticity catalog columns
sys_split.split.cat_splits_lin_full(i3,cols=cols)

# Show the relative disagreement between 2pt statistics when using only galaxies split into each half of a given quantity
sys_split.split.cat_splits_2pt(i3,rm,cols=cols)

i3.wt=True
i3.bs=True
# Produce plots and statistics on linear correlation of ellipticity with catalog columns
sys_split.split.cat_splits_lin_e(i3,cols=cols)
# Produce plots and statistics on linear correlation of non-ellipticity catalog columns
sys_split.split.cat_splits_lin_full(i3,cols=cols)

# Show the relative disagreement between 2pt statistics when using only galaxies split into each half of a given quantity
sys_split.split.cat_splits_2pt(i3,rm,cols=cols)


cols=['coadd','expnum','xoff','yoff','psf1_exp','psf2_exp','ccd','row','col','e1','e2']
i3epoch=catalog.CatalogStore('y1_i3_sv_epoch_v1',cutfunc=None,cattype='i3',cols=cols,catdir='/share/des/disc2/y1/im3shape/single_band/r/y1v1/complete/epoch/',release='y1',maxiter=500)


