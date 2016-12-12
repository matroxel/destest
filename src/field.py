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
  def loop_epoch(nbc,val='e',key='e',scale=0.01,catdir='/share/des/disc2/y1/im3shape/single_band/r/y1v1/complete/epoch/',catname='i3',cattype='i3epoch'):
    """
    """
    import glob

    # whisker store
    y=[]
    x=[]
    mw=[]
    e=[]
    e1=[]
    e2=[]

    for ii in range(len(glob.glob(catdir))):

      epoch=catalog.CatalogStore(catname,cutfunc=None,cattype=cattype,cols=['coadd','e1','e2','ccd','col','row','psf1','psf2'],catdir=catdir,release='y1',maxrows=10000000,maxiter=50,exiter=ii)
      epoch.wt=True
      epoch.bs=True
      x=np.in1d(epoch.coadd,nbc[:,0],assume_unique=False)
      catalog.CatalogMethods.match_cat(epoch,x)
      mask=np.argsort(epoch.coadd)
      catalog.CatalogMethods.match_cat(epoch,mask)
      x=np.in1d(nbc[:,0],np.unique(epoch.coadd),assume_unique=False)
      nbc=nbc[x,:]
      diff=np.diff(epoch.coadd)
      diff=np.where(diff!=0)[0]+1
      diff=np.append([0],diff)
      diff=np.append(diff,[None])

      epoch.m=np.zeros(len(epoch.coadd))
      epoch.c1=np.zeros(len(epoch.coadd))
      epoch.c2=np.zeros(len(epoch.coadd))
      epoch.w=np.zeros(len(epoch.coadd))
      for i in range(len(diff)-1):
        # if i%1000==0:
        #   print i
        epoch.m1[diff[i]:diff[i+1]]=nbc[i,1]
        epoch.m2[diff[i]:diff[i+1]]=nbc[i,2]
        epoch.c1[diff[i]:diff[i+1]]=nbc[i,3]
        epoch.c2[diff[i]:diff[i+1]]=nbc[i,4]
        epoch.w[diff[i]:diff[i+1]]=nbc[i,5]

      tmp=[y,x,mw,e1,e2,e]
      for i,x in enumerate(field.whisker_loop(epoch,tmp)):
        print 'nums',len(tmp[i]),x,tmp[i],x
        if ii==0:
          tmp[i]=x
        else:
          tmp[i]+=x

    field.loop_epoch_finalise(cat,val,key,scale,[x for x in tmp])

    return 

  @staticmethod
  def whisker_calc(cat,col='e'):

    st=[[],[],[],[],[],[]]
    for i,x in enumerate(field.whisker_loop(cat,col=col)):
      print 'nums',len(st[i]),x,st[i],x
      if ii==0:
        st[i]=x
      else:
        st[i]+=x

    return tmp

  @staticmethod
  def loop_epoch_finalise(cat,val,key,scale,y,x,mw,e1,e2,e0):
    """
    """

    pos0=0.5*np.arctan2(e2/mw,e1/mw)
    e0/=mw
    for i in range(len(x)):
      x[i,:,:]+=field_methods.ccd_centres()[i,1]-field_methods.ccdx/2.
      y[i,:,:]+=field_methods.ccd_centres()[i,0]-field_methods.ccdy/2.
    fig.plot_methods.plot_whisker(y,x,np.sin(pos0)*e0,np.cos(pos0)*e0,name=cat.name,label=val,scale=scale,key=r'$\langle '+key+'\rangle$')

    return

  @staticmethod
  def whisker_loop(cat,col='e',nx=8,ny=4,label='',plot=False):
    """
    Calculate whisker plot for e and psf e over field of view.
    """

    cx=field_methods.ccd_centres()[:,1]
    cy=field_methods.ccd_centres()[:,0]

    dc=2048./4.

    x=np.zeros((len(cx),nx,ny))
    y=np.zeros_like(x)
    e0=np.zeros_like(x)
    mw=np.zeros_like(x)
    e1=np.zeros_like(x)
    e2=np.zeros_like(x)
    for i in range(len(cx)):
      if (i==1)|(i==30)|(i==60):
        continue
      print 'chip',i
      #pos1=2.*(cat.pos[mask&(cat.ccd==i)]-np.pi/2.)
      mask=(cat.ccd==i)
      e1=getattr(cat,col+'1')[mask]
      e2=getattr(cat,col+'2')[mask]
      if cat.bs:
        e1-=cat.c1[mask]
        e2-=cat.c2[mask]
        m=cat.m[mask]
      else:
        m=np.ones(np.sum(mask))
      if cat.wt:
        w=cat.w[mask0]
      else:
        w=np.ones(np.sum(mask))
      e1,x0,y0=np.histogram2d(cat.row[mask],cat.col[mask],bins=[nx,ny],weights=e1*w)
      e2,x0,y0=np.histogram2d(cat.row[mask],cat.col[mask],bins=[nx,ny],weights=e2*w)
      e0,x0,y0=np.histogram2d(cat.row[mask],cat.col[mask],bins=[nx,ny],weights=np.sqrt(e1**2+e2**2)*w)
      mw,x0,y0=np.histogram2d(cat.row[mask],cat.col[mask],bins=[nx,ny],weights=m*w)

      x[i,:,:]=(x0[1:]+x0[:-1])/2
      y[i,:,:]=(y0[1:]+y0[:-1])/2

    return y,x,mw,e1,e2,e

  @staticmethod
  def whisker(cat,mask=None,label='',plot=False):
    """
    Calculate whisker plot for e and psf e over field of view.
    """

    mask=catalog.CatalogMethods.check_mask(cat.coadd,mask)

    if not hasattr(cat, 'ra'):
      cat.ra,cat.dec=field_methods.get_field_pos(cat)

    #x,y=field_methods.get_field_pos(cat)

    cx=field_methods.ccd_centres()[:,1]
    cy=field_methods.ccd_centres()[:,0]

    dc=2048./4.

    x0=[]
    y0=[]
    pos0=[]
    psfpos0=[]
    e0=[]
    m0=[]
    psf0=[]
    pos1=[]
    psfpos1=[]
    e1=[]
    psf1=[]
    for i in range(len(cx)):
      if (i==1)|(i==30)|(i==60):
        continue
      print 'chip',i
      #pos1=2.*(cat.pos[mask&(cat.ccd==i)]-np.pi/2.)
      mask0=mask&(cat.ccd==i)
      e1=cat.e1[mask0]
      e2=cat.e2[mask0]
      if cat.bs:
        e1-=cat.c1[mask0]
        e2-=cat.c2[mask0]
        m=cat.m[mask0]
      if cat.wt:
        w=cat.w[mask0]
      else:
        w=np.ones(np.sum(mask0))
      psf1=cat.psf1[mask0]
      psf2=cat.psf2[mask0]
      x1=cat.row[mask0]
      y1=cat.col[mask0]
      # pos1=np.append(pos1,np.mean(0.5*np.arctan2(e2,e1)))
      # psfpos1=np.append(psfpos1,np.mean(0.5*np.arctan2(psf2,psf1)))
      # e1=np.append(e1,np.mean(np.sqrt(e1**2.+e2**2.)))
      # psf1=np.append(psf1,np.mean(np.sqrt(psf1**2.+psf2**2.)))
      for j in xrange(4):
        for k in xrange(8):
          x0=np.append(x0,cx[i]-field_methods.ccdx/2.+(j+.5)*field_methods.ccdx/8.)
          y0=np.append(y0,cy[i]-field_methods.ccdy/2.+(k+.5)*field_methods.ccdy/4.)
          mask1=(x1>k*dc)&(x1<=(k+1.)*dc)&(y1>j*dc)&(y1<=(j+1.)*dc)
          pos0=np.append(pos0,0.5*np.arctan2(np.average(e2[mask1],weights=w[mask1]),np.average(e1[mask1],weights=w[mask1])))
          psfpos0=np.append(psfpos0,0.5*np.arctan2(np.average(psf2[mask1],weights=w[mask1]),np.average(psf1[mask1],weights=w[mask1])))
          e0=np.append(e0,np.average(np.sqrt(e1[mask1]**2.+e2[mask1]**2.),weights=w[mask1]))
          if cat.bs:
            m0=np.append(m0,np.average(1.+m[mask1]))
          else:
            m0=np.append(m0,1.)
          psf0=np.append(psf0,np.average(np.sqrt(psf1[mask1]**2.+psf2[mask1]**2.),weights=w[mask1]))

    # fig.plot_methods.plot_whisker(cy,cx,np.sin(pos1)*e1,np.cos(pos1)*e1,name=cat.name,label='shear2',scale=.01,key=r'$\langle e\rangle$')
    # fig.plot_methods.plot_whisker(cy,cx,np.sin(psfpos1)*psf1,np.cos(psfpos1)*psf1,name=cat.name,label='psf2',scale=.01,key=r'$\langle$ PSF $e\rangle$')
    if plot:
      fig.plot_methods.plot_whisker(y0,x0,np.sin(pos0)*e0/m0,np.cos(pos0)*e0/m0,name=cat.name,label='shear'+label,scale=.01,key=r'$\langle e\rangle$')
      fig.plot_methods.plot_whisker(y0,x0,np.sin(psfpos0)*psf0,np.cos(psfpos0)*psf0,name=cat.name,label='psf'+label,scale=.01,key=r'$\langle$ PSF $e\rangle$')

    return y0,x0,m0,pos0,e0,psfpos0,psf0

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
      print 'chip', i
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

      fig.plot_methods.plot_whisker(y0,x0,np.sin(pos0)*e0,np.cos(pos0)*e0,name=cat.name,label='chip_'+str(i)+'_shear',scale=.01,key=r'$\langle e\rangle$',chip=True)
      fig.plot_methods.plot_whisker(y0,x0,np.sin(psfpos0)*psf0,np.cos(psfpos0)*psf0,name=cat.name,label='chip_'+str(i)+'_psf',scale=.01,key=r'$\langle$ PSF $e\rangle$',chip=True)

    return

  @staticmethod
  def footprint(cat,mask=None):
    """
    Calculate and plot object number density over field of view.
    """

    mask=catalog.CatalogMethods.check_mask(cat.coadd,mask)

    if not hasattr(cat, 'fx'):
      cat.fx,cat.fy=field_methods.get_field_pos(cat)

    fig.plot_methods.plot_field_footprint(cat.fx,cat.fy,cat.name,label='field',bins=1)

    return

  @staticmethod
  def loop_submit_sp():

    import mpi4py.MPI
    from mpi_pool import MPIPool
    import sys
    comm = mpi4py.MPI.COMM_WORLD
    pool = MPIPool(comm=comm,debug=True)
    
    if pool.is_master():
        commands = range(config.nchunk)
        sys.stdout.flush()
        print
        print
    else:
        commands = None
    comm.Barrier()
    pool.map(field.field.build_special_points, commands)
    pool.close()
    comm.Barrier()

    return


  @staticmethod
  def build_special_points_fits(sp=None): 
    """
    Combines parts of special points catalog into single fits catalog.
    """

    import fitsio as fio

    name=['center','ll','ul','lr','ur','tb1','tb2','tb3','tb4','tb5','tb6']
    #tiles=fio.FITS(config.coaddtiles)[-1].read()
    wcs=fio.FITS(config.wcsfile)[-1].read()
    a=np.sort(np.unique(wcs['expnum']))
    b=np.sort(np.unique(wcs['ccdnum']))-1
    store=np.empty((len(a)*len(b)*12),dtype=[('exposure',int)]+[('ccd',int)]+[('band','S1')]+[('type',int)]+[('ra','f8')]+[('dec','f8')])
    print len(store)
    # for i in range(len(a)):
    #   store['exposure'][i*len(b):(i+1)*len(b)]=a[i]
    #   for j in range(len(b)):
    #     store['ccd'][i*len(b)+j]=b[j]

    if sp is None:
      for i in range(config.nchunk):
        print i
        tmp=np.genfromtxt('y1a1_special_points_'+str(i)+'.txt',names=['index','exposure','ccd','band','racenter','deccenter','rall','decll','raul','decul','ralr','declr','raur','decur','ratb1','dectb1','ratb2','dectb2','ratb3','dectb3','ratb4','dectb4','ratb5','dectb5','ratb6','dectb6'],dtype=None)
        if i==0:
          sp=tmp
        else:
          sp=np.append(sp,tmp)
        a=np.argsort(sp,order=('exposure','ccd'))
        sp=sp[a]

    ind0=-1
    for i in range(len(sp)):
      if i%10000==0:
        print i
      if i>0:
        if (sp['ccd'][i]==sp['ccd'][i-1])&(sp['exposure'][i]==sp['exposure'][i-1]):
          continue
      ind=[]
      for j in range(100):
        if i+j+1==len(sp):
          break
        if (sp['ccd'][i]==sp['ccd'][i+j+1])&(sp['exposure'][i]==sp['exposure'][i+j+1]):
          ind.append(j)
        else:
          break
      for k in range(11):
        ind0+=1
        if ind0==len(store):
          print i,len(sp)
          break
        store['exposure'][ind0]=sp['exposure'][i]
        store['ccd'][ind0]=sp['ccd'][i]
        store['band'][ind0]=sp['band'][i]
        store['type'][ind0]=k
        store['ra'][ind0]=np.mean(sp['ra'+name[k]][i:i+j+1])
        store['dec'][ind0]=np.mean(sp['dec'+name[k]][i:i+j+1])

    # FOV centers from two center chips
    mask = ((store['ccd']==28)|(store['ccd']==33))&(store['type']==0)
    ra=store['ra'][mask]
    dec=store['dec'][mask]
    store['exposure'][ind0+1:ind0+1+len(ra)/2]=(store['exposure'][mask])[::2]
    store['ccd'][ind0+1:ind0+1+len(ra)/2]=-1
    store['band'][ind0+1:ind0+1+len(ra)/2]=(store['band'][mask])[::2]
    store['type'][ind0+1:ind0+1+len(ra)/2]=-1
    store['ra'][ind0+1:ind0+1+len(ra)/2]=(ra[::2]+ra[1::2])/2
    store['dec'][ind0+1:ind0+1+len(ra)/2]=(dec[::2]+dec[1::2])/2

    store=store[:ind0+1+len(ra)/2]

    fio.write(config.spointsfile,store,clobber=True)

    return

  @staticmethod
  def mean_shear(cat,mask):
    """
    Calculate and plot mean shear as a function of pixel row/column.
    """

    if not hasattr(cat, 'ra'):
      cat.ra,cat.dec=field_methods.get_field_pos(cat)

    cat.lbins=50
    split.split_methods.split_gals_lin_along(cat,'ra',mask=mask,plot=True,label='field-row',fit=True)
    split.split_methods.split_gals_lin_along(cat,'ra',mask=mask,plot=True,label='field-col',fit=True)


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
    pointings.type=tmp['type']
    pointings.ccd=tmp['ccd']
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

  @staticmethod
  def get_radec_coadd_tiles(tiles=None,tiles0=None,file=config.coaddtiles):

    if tiles is None:
      tiles=fio.FITS(file)[-1].read()

    if tiles0 is None:
      mask=np.ones(len(tiles)).astype(bool)
    else:
      mask=np.in1d(np.core.defchararray.strip(tiles['TILENAME']),tiles0,assume_unique=False)

    return tiles,np.vstack(((tiles['URAUR'][mask]+tiles['URALL'][mask])/2.,(tiles['UDECUR'][mask]+tiles['UDECLL'][mask])/2.)).T


def build_special_points(chunk):
  """
  Used to build parts of catalog of special points.
  """

  import re

  dchunk=int(fio.FITS(config.wcsfile)[-1].get_nrows())/config.nchunk
  ia=dchunk*chunk
  print ia
  ib=dchunk*(chunk+1)
  if chunk==config.nchunk-1:
    ib=int(fio.FITS(config.wcsfile)[-1].get_nrows())
  print ib

  with open(config.y1blacklist) as f:
    lines = f.readlines()
  blexp=[]
  blccd=[]
  for line in lines:
    blexp=np.append(blexp,int(re.compile('\w+').findall(line)[1][6:]))
    blccd=np.append(blccd,int(re.compile('\w+').findall(line)[2]))

  tmp=fio.FITS(config.wcsfile)[-1][ia:ib]
  image=np.empty(tmp.shape, dtype=tmp.dtype.descr + [('naxis1',int)]+[('naxis2',int)])
  for name in tmp.dtype.names:
    image[name]=tmp[name]

  image['naxis1']=np.ones(len(image))*2048
  image['naxis2']=np.ones(len(image))*4096

  tb = np.genfromtxt('../tape_bumps.txt',names=['ccd','t','l','b','r'],delimiter=',')

  for i in range(ib-ia):
    if image['expnum'][i] in blexp:
      if image['ccdnum'][i] in blccd[blexp==image['expnum'][i]]:
        continue
    # print i,str(image['expnum'][i])+' '+str(image['ccdnum'][i])
    line=str(i)+' '+str(image['expnum'][i])+' '+str(image['ccdnum'][i])+' '+str(image['band'][i])+' '
    rapos=[1024,0,2048,0,2048]
    decpos=[2048,0,4096,0,4096]
    tbmask=tb['ccd']==image['ccdnum'][i]
    for j in range(6):
      decpos.append(int((tb[tbmask]['t'][j]+tb[tbmask]['b'][j])/2))
      rapos.append(int((tb[tbmask]['l'][j]+tb[tbmask]['r'][j])/2))
    radec=field_methods.translate_to_wcs([rapos,decpos],image[i])
    # if field_methods.get_coadd_tile(radec[0],radec[1],tiles=tiles) in image['tilename'][i]:
    for j in range(11):
      line+=str(radec[0][j])+' '+str(radec[1][j])+' '

    with open('y1a1_special_points_'+str(chunk)+'.txt','a') as f:
      f.write(line+'\n')

  return
