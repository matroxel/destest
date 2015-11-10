import numpy as np
import fitsio as fio

import catalog
import config
import fig
import txt
import lin
import corr

class field(object):

  @staticmethod
  def whisker(cat,mask=None):
    """
    Calculate whisker plot for e and psf e over field of view.
    """

    mask=catalog.CatalogMethods.check_mask(cat.coadd,mask)

    #x,y=field_methods.get_field_pos(cat)

    cx=field_methods.ccd_centres()[:,1]
    cy=field_methods.ccd_centres()[:,0]

    dc=2048./4.

    x0=[]
    y0=[]
    pos0=[]
    psfpos0=[]
    e0=[]
    psf0=[]
    for i in range(len(cx)):
      if (i==1)|(i==30)|(i==60):
        continue
      print 'chip',i
      #pos1=2.*(cat.pos[mask&(cat.ccd==i)]-np.pi/2.)
      e1=cat.e1[mask&(cat.ccd==i)]
      e2=cat.e2[mask&(cat.ccd==i)]
      psf1=cat.psf1_exp[mask&(cat.ccd==i)]
      psf2=cat.psf2_exp[mask&(cat.ccd==i)]
      x1=cat.row[mask&(cat.ccd==i)]
      y1=cat.col[mask&(cat.ccd==i)]
      for j in xrange(4):
        for k in xrange(8):
          x0=np.append(x0,cx[i]-field_methods.ccdx/2.+(j+.5)*field_methods.ccdx/8.)
          y0=np.append(y0,cy[i]-field_methods.ccdy/2.+(k+.5)*field_methods.ccdy/4.)
          mask1=(x1>k*dc)&(x1<=(k+1.)*dc)&(y1>j*dc)&(y1<=(j+1.)*dc)
          pos0=np.append(pos0,0.5*np.arctan2(np.mean(e2[mask1]),np.mean(e1[mask1])))
          psfpos0=np.append(psfpos0,0.5*np.arctan2(np.mean(psf2[mask1]),np.mean(psf1[mask1])))
          e0=np.append(e0,np.sqrt(np.mean(e1[mask1])**2.+np.mean(e2[mask1])**2.))
          psf0=np.append(psf0,np.sqrt(np.mean(psf1[mask1])**2.+np.mean(psf2[mask1])**2.))

    fig.plot_methods.plot_whisker(y0,x0,np.sin(pos0)*e0,np.cos(pos0)*e0,name=cat.name,label='shear',scale=.01,key=r'$\langle e\rangle$')
    fig.plot_methods.plot_whisker(y0,x0,np.sin(psfpos0)*psf0,np.cos(psfpos0)*psf0,name=cat.name,label='psf',scale=.01,key=r'$\langle$ PSF $e\rangle$')

    return

  @staticmethod
  def whisker_chip(cat,mask=None):
    """
    Calculate whisker plot for e and psf e over each chip.
    """

    mask=catalog.CatalogMethods.check_mask(cat.coadd,mask)

    #x,y=field_methods.get_field_pos(cat)

    cx=field_methods.ccd_centres()[:,1]
    cy=field_methods.ccd_centres()[:,0]

    dc=2048./20.

    e1=cat.e1[mask]
    e2=cat.e2[mask]
    psf1=cat.psf1_exp[mask]
    psf2=cat.psf2_exp[mask]
    x1=cat.row[mask]
    y1=cat.col[mask]
    for i in xrange(len(cx)):
      x0=[]
      y0=[]
      pos0=[]
      psfpos0=[]
      e0=[]
      psf0=[]
      for j in xrange(20):
        for k in xrange(40):
          x0=np.append(x0,j*dc)
          y0=np.append(y0,k*dc)
          mask1=(x1>k*dc)&(x1<=(k+1.)*dc)&(y1>j*dc)&(y1<=(j+1.)*dc)
          pos0=np.append(pos0,0.5*np.arctan2(np.mean(e2[mask1]),np.mean(e1[mask1])))
          psfpos0=np.append(psfpos0,0.5*np.arctan2(np.mean(psf2[mask1]),np.mean(psf1[mask1])))
          e0=np.append(e0,np.sqrt(np.mean(e1[mask1])**2.+np.mean(e2[mask1])**2.))
          psf0=np.append(psf0,np.sqrt(np.mean(psf1[mask1])**2.+np.mean(psf2[mask1])**2.))

      fig.plot_methods.plot_whisker(y0,x0,np.sin(pos0)*e0,np.cos(pos0)*e0,name=cat.name,label='chip_'+str(i)+'_shear',scale=.01,key=r'$\langle e\rangle$')
      fig.plot_methods.plot_whisker(y0,x0,np.sin(psfpos0)*psf0,np.cos(psfpos0)*psf0,name=cat.name,label='chip_'+str(i)+'_psf',scale=.01,key=r'$\langle$ PSF $e\rangle$')

    return

  @staticmethod
  def footprint(cat,mask=None):
    """
    Calculate and plot object number density over field of view.
    """

    mask=catalog.CatalogMethods.check_mask(cat.coadd,mask)

    cat.ra,cat.dec=field_methods.get_field_pos(cat)
    fig.plot_methods.plot_footprint(cat,mask=mask,label='field',bins=2)

    return

  @staticmethod
  def build_special_points(chunk):
    """
    Used to build parts of catalog of special points.
    """

    dchunk=int(fio.FITS(config.wcsfile)[-1].get_nrows())/40
    ia=dchunk*chunk
    print ia
    ib=dchunk*(chunk+1)
    if chunk==39:
      ib=int(fio.FITS(config.wcsfile)[-1].get_nrows()) 
    print ib

    tmp=fio.FITS(config.wcsfile)[-1][ia:ib]
    tiles=fio.FITS(config.coaddtiles)[-1].read()

    image=np.empty(tmp.shape, dtype=tmp.dtype.descr + [('naxis1',int)]+[('naxis2',int)])
    for name in tmp.dtype.names:
      image[name]=tmp[name]    

    image['naxis1']=np.ones(len(image))*2048
    image['naxis2']=np.ones(len(image))*4096

    for i in range(ib-ia):
      print i
      print str(image['expnum'][i])+' '+str(image['ccdnum'][i])
      line=str(i)+' '+str(image['expnum'][i])+' '+str(image['ccdnum'][i])+' '
      radec=field_methods.translate_to_wcs([1024,2048],image[i])
      if field_methods.get_coadd_tile(radec[0],radec[1],tiles=tiles) in image['tilename'][i]:
        line+=str(radec[0])+' '+str(radec[1])+' '
      else:
        line+=str(999)+' '+str(999)+' '        
      radec=field_methods.translate_to_wcs([0,0],image[i])
      if field_methods.get_coadd_tile(radec[0],radec[1],tiles=tiles) in image['tilename'][i]:
        line+=str(radec[0])+' '+str(radec[1])+' '
      else:
        line+=str(999)+' '+str(999)+' '        
      radec=field_methods.translate_to_wcs([2048,0],image[i])
      if field_methods.get_coadd_tile(radec[0],radec[1],tiles=tiles) in image['tilename'][i]:
        line+=str(radec[0])+' '+str(radec[1])+' '
      else:
        line+=str(999)+' '+str(999)+' '        
      radec=field_methods.translate_to_wcs([0,4096],image[i])
      if field_methods.get_coadd_tile(radec[0],radec[1],tiles=tiles) in image['tilename'][i]:
        line+=str(radec[0])+' '+str(radec[1])+' '
      else:
        line+=str(999)+' '+str(999)+' '        
      radec=field_methods.translate_to_wcs([2048,4096],image[i])
      if field_methods.get_coadd_tile(radec[0],radec[1],tiles=tiles) in image['tilename'][i]:
        line+=str(radec[0])+' '+str(radec[1])+' '
      else:
        line+=str(999)+' '+str(999)+' '        
      with open('y1a1_special_points_'+str(chunk)+'.txt','a') as f:
        f.write(line+'\n')

    return

  @staticmethod
  def build_special_points_fits(): 
    """
    Combines parts of special points catalog into single fits catalog.
    """

    import fitsio as fio

    tmp=fio.FITS(config.wcsfile)[-1].read()
    a=np.sort(np.unique(tmp['expnum']))
    b=np.sort(np.unique(tmp['ccdnum']))-1
    store=np.empty((len(a)*(len(b)+1)),dtype=[('exposure',int)]+[('ccd',int)]+[('type',int)]+[('ra','f8')]+[('dec','f8')])
    for i in range(len(a)):
      store['exposure'][i*len(b):(i+1)*len(b)]=a[i]
      for j in range(len(b)):
        store['ccd'][i*len(b)+j]=b[j]

    for i in range(40):
      print i
      tmp=np.genfromtxt('y1a1_special_points_'+str(i)+'.txt',names=['index','exposure','ccd','racenter','deccenter','rall','decll','raul','decul','ralr','declr','raur','decur'])
      for j in range(len(tmp)):
        if j%1000==0:
          print j
        mask=(store['exposure']==tmp['exposure'][j])&(store['ccd']==tmp['ccd'][j])
        if tmp['racenter'][j]!=999:
          store['type'][mask]=0
          store['ra'][mask]=tmp['racenter'][j]
          store['dec'][mask]=tmp['deccenter'][j]
        if tmp['rall'][j]!=999:
          store['type'][mask]=1
          store['ra'][mask]=tmp['rall'][j]
          store['dec'][mask]=tmp['decll'][j]
        if tmp['raul'][j]!=999:
          store['type'][mask]=2
          store['ra'][mask]=tmp['raul'][j]
          store['dec'][mask]=tmp['decul'][j]
        if tmp['ralr'][j]!=999:
          store['type'][mask]=3
          store['ra'][mask]=tmp['ralr'][j]
          store['dec'][mask]=tmp['declr'][j]
        if tmp['raur'][j]!=999:
          store['type'][mask]=4
          store['ra'][mask]=tmp['raur'][j]
          store['dec'][mask]=tmp['decur'][j]

    for i in range(len(a)):
      if i%1000==0:
        print i
      store['exposure'][len(a)*len(b)+i]=a[i]
      store['ccd'][len(a)*len(b)+i]=-1
      store['type'][len(a)*len(b)+i]=-1
      mask=(store['exposure']==store['exposure'][len(a)*len(b)+i])&(store['type']==0)&((store['ccd']==27)|(store['ccd']==34))
      store['ra'][len(a)*len(b)+i]=np.mean(store['ra'][mask])
      store['dec'][len(a)*len(b)+i]=np.mean(store['dec'][mask])

    fio.write(config.spointsfile,array,clobber=True)

    return

  @staticmethod
  def corr_points(cat,mask):
    """
    Calculate and plot tangential shear and mean shear around special points in catalog.
    """

    import fitsio as fio

    mask=catalog.CatalogMethods.check_mask(cat.coadd,mask)

    fits=fio.FITS(config.spointsfile)
    tmp=fits[-1].read()
    pointings=catalog.CatalogStore('field_points',setup=False)
    pointings.coadd=np.arange(len(tmp))
    pointings.ra=tmp['ra']
    pointings.dec=tmp['dec']
    pointings.tbins=50
    pointings.sep=np.array([.1,500.])

    mask1=tmp['ccd']==-1
    theta,out,err,chi2=corr.xi_2pt.xi_2pt(pointings,catb=cat,k=None,ga=None,gb=None,corr='NG',maska=mask1,maskb=mask,wa=None,wb=None,ran=False,mock=False,erron=True,jkmask=None,label0='',plot=False)
    fig.plot_methods.plot_field_corr(cat,theta,out,err,'centre')
    theta,out,err,chi2=corr.xi_2pt.xi_2pt(pointings,catb=cat,k='e1',ga=None,gb=None,corr='NK',maska=mask1,maskb=mask,wa=None,wb=None,ran=False,mock=False,erron=True,jkmask=None,label0='',plot=False)
    theta,out2,err2,chi2=corr.xi_2pt.xi_2pt(pointings,catb=cat,k='e2',ga=None,gb=None,corr='NK',maska=mask1,maskb=mask,wa=None,wb=None,ran=False,mock=False,erron=True,jkmask=None,label0='',plot=False)
    fig.plot_methods.plot_field_corr2(cat,theta,out,err,out2,err2,'centre')

    mask1=tmp['type']==0
    theta,out,err,chi2=corr.xi_2pt.xi_2pt(pointings,catb=cat,k=None,ga=None,gb=None,corr='NG',maska=mask1,maskb=mask,wa=None,wb=None,ran=False,mock=False,erron=True,jkmask=None,label0='',plot=False)
    fig.plot_methods.plot_field_corr(cat,theta,out,err,'chip-centre')
    theta,out,err,chi2=corr.xi_2pt.xi_2pt(pointings,catb=cat,k='e1',ga=None,gb=None,corr='NK',maska=mask1,maskb=mask,wa=None,wb=None,ran=False,mock=False,erron=True,jkmask=None,label0='',plot=False)
    theta,out2,err2,chi2=corr.xi_2pt.xi_2pt(pointings,catb=cat,k='e2',ga=None,gb=None,corr='NK',maska=mask1,maskb=mask,wa=None,wb=None,ran=False,mock=False,erron=True,jkmask=None,label0='',plot=False)
    fig.plot_methods.plot_field_corr2(cat,theta,out,err,out2,err2,'chip-centre')

    mask1=tmp['type']==1
    theta,out,err,chi2=corr.xi_2pt.xi_2pt(pointings,catb=cat,k=None,ga=None,gb=None,corr='NG',maska=mask1,maskb=mask,wa=None,wb=None,ran=False,mock=False,erron=True,jkmask=None,label0='',plot=False)
    fig.plot_methods.plot_field_corr(cat,theta,out,err,'chip-cornera')
    theta,out,err,chi2=corr.xi_2pt.xi_2pt(pointings,catb=cat,k='e1',ga=None,gb=None,corr='NK',maska=mask1,maskb=mask,wa=None,wb=None,ran=False,mock=False,erron=True,jkmask=None,label0='',plot=False)
    theta,out2,err2,chi2=corr.xi_2pt.xi_2pt(pointings,catb=cat,k='e2',ga=None,gb=None,corr='NK',maska=mask1,maskb=mask,wa=None,wb=None,ran=False,mock=False,erron=True,jkmask=None,label0='',plot=False)
    fig.plot_methods.plot_field_corr2(cat,theta,out,err,out2,err2,'chip-cornera')

    mask1=tmp['type']==2
    theta,out,err,chi2=corr.xi_2pt.xi_2pt(pointings,catb=cat,k=None,ga=None,gb=None,corr='NG',maska=mask1,maskb=mask,wa=None,wb=None,ran=False,mock=False,erron=True,jkmask=None,label0='',plot=False)
    fig.plot_methods.plot_field_corr(cat,theta,out,err,'chip-cornerb')
    theta,out,err,chi2=corr.xi_2pt.xi_2pt(pointings,catb=cat,k='e1',ga=None,gb=None,corr='NK',maska=mask1,maskb=mask,wa=None,wb=None,ran=False,mock=False,erron=True,jkmask=None,label0='',plot=False)
    theta,out2,err2,chi2=corr.xi_2pt.xi_2pt(pointings,catb=cat,k='e2',ga=None,gb=None,corr='NK',maska=mask1,maskb=mask,wa=None,wb=None,ran=False,mock=False,erron=True,jkmask=None,label0='',plot=False)
    fig.plot_methods.plot_field_corr2(cat,theta,out,err,out2,err2,'chip-cornerb')

    mask1=tmp['type']==3
    theta,out,err,chi2=corr.xi_2pt.xi_2pt(pointings,catb=cat,k=None,ga=None,gb=None,corr='NG',maska=mask1,maskb=mask,wa=None,wb=None,ran=False,mock=False,erron=True,jkmask=None,label0='',plot=False)
    fig.plot_methods.plot_field_corr(cat,theta,out,err,'chip-cornerc')
    theta,out,err,chi2=corr.xi_2pt.xi_2pt(pointings,catb=cat,k='e1',ga=None,gb=None,corr='NK',maska=mask1,maskb=mask,wa=None,wb=None,ran=False,mock=False,erron=True,jkmask=None,label0='',plot=False)
    theta,out2,err2,chi2=corr.xi_2pt.xi_2pt(pointings,catb=cat,k='e2',ga=None,gb=None,corr='NK',maska=mask1,maskb=mask,wa=None,wb=None,ran=False,mock=False,erron=True,jkmask=None,label0='',plot=False)
    fig.plot_methods.plot_field_corr2(cat,theta,out,err,out2,err2,'chip-cornerc')

    mask1=tmp['type']==4
    theta,out,err,chi2=corr.xi_2pt.xi_2pt(pointings,catb=cat,k=None,ga=None,gb=None,corr='NG',maska=mask1,maskb=mask,wa=None,wb=None,ran=False,mock=False,erron=True,jkmask=None,label0='',plot=False)
    fig.plot_methods.plot_field_corr(cat,theta,out,err,'chip-cornerd')
    theta,out,err,chi2=corr.xi_2pt.xi_2pt(pointings,catb=cat,k='e1',ga=None,gb=None,corr='NK',maska=mask1,maskb=mask,wa=None,wb=None,ran=False,mock=False,erron=True,jkmask=None,label0='',plot=False)
    theta,out2,err2,chi2=corr.xi_2pt.xi_2pt(pointings,catb=cat,k='e2',ga=None,gb=None,corr='NK',maska=mask1,maskb=mask,wa=None,wb=None,ran=False,mock=False,erron=True,jkmask=None,label0='',plot=False)
    fig.plot_methods.plot_field_corr2(cat,theta,out,err,out2,err2,'chip-cornerd')

    mask1=tmp['type']>0
    theta,out,err,chi2=corr.xi_2pt.xi_2pt(pointings,catb=cat,k=None,ga=None,gb=None,corr='NG',maska=mask1,maskb=mask,wa=None,wb=None,ran=False,mock=False,erron=True,jkmask=None,label0='',plot=False)
    fig.plot_methods.plot_field_corr(cat,theta,out,err,'chip-corner')
    theta,out,err,chi2=corr.xi_2pt.xi_2pt(pointings,catb=cat,k='e1',ga=None,gb=None,corr='NK',maska=mask1,maskb=mask,wa=None,wb=None,ran=False,mock=False,erron=True,jkmask=None,label0='',plot=False)
    theta,out2,err2,chi2=corr.xi_2pt.xi_2pt(pointings,catb=cat,k='e2',ga=None,gb=None,corr='NK',maska=mask1,maskb=mask,wa=None,wb=None,ran=False,mock=False,erron=True,jkmask=None,label0='',plot=False)
    fig.plot_methods.plot_field_corr2(cat,theta,out,err,out2,err2,'chip-corner')

    return theta,out,err,chi2    


class field_methods(object):
  """
  Utilities for doing pixel and chip calculations.
  """

  chip_centres = {

  'N7':[16.908,191.670],
  'N6':[16.908,127.780],
  'N5':[16.908,63.890],
  'N4':[16.908,0.],
  'N3':[16.908,-63.890],
  'N2':[16.908,-127.780],
  'N1':[16.908,-191.670],
  'N13':[50.724,159.725],
  'N12':[50.724,95.835],
  'N11':[50.724,31.945],
  'N10':[50.724,-31.945],
  'N9':[50.724,-95.835],
  'N8':[50.724,-159.725],
  'N19':[84.540,159.725],
  'N18':[84.540,95.835],
  'N17':[84.540,31.945],
  'N16':[84.540,-31.945],
  'N15':[84.540,-95.835],
  'N14':[84.540,-159.725],
  'N24':[118.356,127.780],
  'N23':[118.356,63.890],
  'N22':[118.356,0.],
  'N21':[118.356,-63.890],
  'N20':[118.356,-127.780],
  'N28':[152.172,95.835],
  'N27':[152.172,31.945],
  'N26':[152.172,-31.945],
  'N25':[152.172,-95.835],
  'N31':[185.988,63.890],
  'N30':[185.988,0.],
  'N29':[185.988,-63.890],
  'S7':[-16.908,191.670],
  'S6':[-16.908,127.780],
  'S5':[-16.908,63.890],
  'S4':[-16.908,0.],
  'S3':[-16.908,-63.890],
  'S2':[-16.908,-127.780],
  'S1':[-16.908,-191.670],
  'S13':[-50.724,159.725],
  'S12':[-50.724,95.835],
  'S11':[-50.724,31.945],
  'S10':[-50.724,-31.945],
  'S9':[-50.724,-95.835],
  'S8':[-50.724,-159.725],
  'S19':[-84.540,159.725],
  'S18':[-84.540,95.835],
  'S17':[-84.540,31.945],
  'S16':[-84.540,-31.945],
  'S15':[-84.540,-95.835],
  'S14':[-84.540,-159.725],
  'S24':[-118.356,127.780],
  'S23':[-118.356,63.890],
  'S22':[-118.356,0.],
  'S21':[-118.356,-63.890],
  'S20':[-118.356,-127.780],
  'S28':[-152.172,95.835],
  'S27':[-152.172,31.945],
  'S26':[-152.172,-31.945],
  'S25':[-152.172,-95.835],
  'S31':[-185.988,63.890],
  'S30':[-185.988,0.],
  'S29':[-185.988,-63.890]

  }

  ccdid=['S29','S30','S31','S25','S26','S27','S28','S20','S21','S22','S23','S24','S14','S15','S16','S17','S18','S19','S8','S9','S10','S11','S12','S13','S1','S2','S3','S4','S5','S6','S7','N1','N2','N3','N4','N5','N6','N7','N8','N9','N10','N11','N12','N13','N14','N15','N16','N17','N18','N19','N20','N21','N22','N23','N24','N25','N26','N27','N28','N29','N30','N31']

  ccdx=4096.*15.e-6*1000.
  ccdy=2048.*15.e-6*1000.

  @staticmethod
  def ccd_centres():

    centrex=[]
    centrey=[]
    for i,x in enumerate(field_methods.ccdid):
      centrex=np.append(centrex,field_methods.chip_centres.get(x,None)[1])
      centrey=np.append(centrey,field_methods.chip_centres.get(x,None)[0])

    return np.vstack((centrex,centrey)).T

  @staticmethod
  def ccd_corners():

    cornersx=[]
    cornersy=[]
    for i,x in enumerate(field_methods.ccdid):
      centrex=np.append(centrex,field_methods.chip_centres.get(x,None)[0]-field_methods.ccdx/2.) # lower left
      centrex=np.append(centrex,field_methods.chip_centres.get(x,None)[0]-field_methods.ccdx/2.) # lower right
      centrex=np.append(centrex,field_methods.chip_centres.get(x,None)[0]+field_methods.ccdx/2.) # upper left
      centrex=np.append(centrex,field_methods.chip_centres.get(x,None)[0]+field_methods.ccdx/2.) # upper right

      centrey=np.append(centrey,field_methods.chip_centres.get(x,None)[1]-field_methods.ccdy/2.)
      centrey=np.append(centrey,field_methods.chip_centres.get(x,None)[1]+field_methods.ccdy/2.)
      centrey=np.append(centrey,field_methods.chip_centres.get(x,None)[1]-field_methods.ccdy/2.)
      centrey=np.append(centrey,field_methods.chip_centres.get(x,None)[1]+field_methods.ccdy/2.)

    return np.vstack((centrex,centrey)).T

  @staticmethod
  def ccd_to_field(ccd,ccdx,ccdy):

    centre=field_methods.ccd_centres()

    centrex=(centre[:,0])[[ccd]]
    centrey=(centre[:,1])[[ccd]]

    return ccdx*15e-6*1000+centrex,ccdy*15e-6*1000+centrey 

  @staticmethod
  def get_field_pos(cat):

    x,y=field_methods.ccd_to_field(cat.ccd,cat.row-field_methods.ccdx/2.,cat.col-field_methods.ccdy/2.)

    return x,y 

  @staticmethod
  def translate_to_wcs(pos,image):

    from esutil import wcsutil
    
    wcs=wcsutil.WCS(image, longpole=180.0, latpole=90.0, theta0=90.0)
    ra,dec=wcs.image2sky(pos[0],pos[1])

    return ra,dec

  @staticmethod
  def get_coadd_tile(ra,dec,tiles=None):

    if tiles is None:
      tiles=fio.FITS(config.coaddtiles)[-1].read()

    tmp=tiles['TILENAME'][(ra<tiles['URAUR'])&(dec<tiles['UDECUR'])&(ra>tiles['URALL'])&(dec>tiles['UDECLL'])]
    if len(tmp)==0:
      tmp=tiles['TILENAME'][((ra+360)<tiles['URAUR'])&(dec<tiles['UDECUR'])&((ra+360)>tiles['URALL'])&(dec>tiles['UDECLL'])]

    return tmp[0].rstrip()
