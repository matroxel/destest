import numpy as np
import os

import catalog
import fig
import txt

vary = {
  
  'omega_m' : False,
  'h0' : False,
  'omega_b' : False,
  'sigma8_input' : False,
  'tau' : False,
  'n_s' : False,
  'A_s' : False,
  'omega_k' : False,
  'w' : False,
  'wa' : False,
  'omnuh2' : False,
  'massless_nu' : False,
  'massive_nu' : False,
  'bias_1' : False,
  'bias_2' : False,
  'bias_3' : False,
  'bias_4' : False,
  'bias_5' : False,
  'bias_6' : False,
  'm1' : False,
  'm2' : False,
  'm3' : False,
  'm4' : False,
  'm5' : False,
  'm6' : False,
  'A' : False,
  'alpha' : False,
  'a_planck' : False

}

prior = {
  
  'omega_m' : False,
  'h0' : False,
  'omega_b' : False,
  'sigma8_input' : False,
  'tau' : False,
  'n_s' : False,
  'A_s' : False,
  'omega_k' : False,
  'w' : False,
  'wa' : False,
  'omnuh2' : False,
  'massless_nu' : False,
  'massive_nu' : False,
  'bias_1' : False,
  'bias_2' : False,
  'bias_3' : False,
  'bias_4' : False,
  'bias_5' : False,
  'bias_6' : False,
  'm1' : False,
  'm2' : False,
  'm3' : False,
  'm4' : False,
  'm5' : False,
  'm6' : False,
  'A' : False,
  'alpha' : False,
  'a_planck' : False

}


dc1_params = {
  
  'omega_m' : (0.05, 0.3156, 0.6),
  'h0' : (0.4, 0.6727, 0.9),
  'omega_b' : (0.02, 0.0491685, 0.09),
  'sigma8_input' : (0.5, 0.831, 1.2),
  'tau' : (0.01, 0.08, 0.8),
  'n_s' : (0.84, 0.9645, 1.06),
  'A_s' : (1e-9, 2.1e-9, 3e-9),
  'omega_k' : (-1.0, 0.0, 1.0),
  'w' : (-2.1, -1.0, 0.0),
  'wa' : (-1.0, 0.0, 1.0),
  'omnuh2' : (0.0, 0.00065, 0.001),
  'massless_nu' : (1.0, 2.046, 4.0),
  'massive_nu' : (0, 1, 2),
  'bias_1' : (-0.3, 0., 0.3),
  'bias_2' : (-0.3, 0., 0.3),
  'bias_3' : (-0.3, 0., 0.3),
  'bias_4' : (-0.3, 0., 0.3),
  'bias_5' : (-0.3, 0., 0.3),
  'bias_6' : (-0.3, 0., 0.3),
  'm1' : (-0.2, 0., 0.2),
  'm2' : (-0.2, 0., 0.2),
  'm3' : (-0.2, 0., 0.2),
  'm4' : (-0.2, 0., 0.2),
  'm5' : (-0.2, 0., 0.2),
  'm6' : (-0.2, 0., 0.2),
  'A' : (-5., 1., 5.0),
  'alpha' : (-5., 0., 5.0),
  'a_planck' : (.5,1.,1.5)

}


dc1_priors = {
  
  'omega_m' : (0.3156, 0.5),
  'h0' : (0.6726, 0.2),
  'omega_b' : (0.0491685, 0.05),
  'sigma8_input' : (0.831, 0.25),
  'tau' : (0.08, 0.5),
  'n_s' : (0.9645, 0.2),
  'A_s' : (2.215e-9, 1e-9),
  'omega_k' : (0.0, 0.1),
  'w' : (-1.0, 1.),
  'wa' : (0.0, 1.),
  'omnuh2' : (0.00065, .0001),
  'massless_nu' : (2.046, 1.),
  'massive_nu' : (1, 1),
  'bias_1' : (0., 0.02),
  'bias_2' : (0., 0.02),
  'bias_3' : (0., 0.02),
  'bias_4' : (0., 0.02),
  'bias_5' : (0., 0.02),
  'bias_6' : (0., 0.02),
  'm1' : (0., 0.02),
  'm2' : (0., 0.02),
  'm3' : (0., 0.02),
  'm4' : (0., 0.02),
  'm5' : (0., 0.02),
  'm6' : (0., 0.02),
  'A' : (1., 1.0),
  'alpha' : (0., 1.),
  'a_planck' : (1.,.1)

}

class run(object):

  @staticmethod
  def loop_submit():

    run.submit(label='dc1_sig8_Om_03',nodes=1,procs=32,hr=48,pts=1000,mneff=0.5,mntol=0.5,ia=False,pz=False,mbias=False,params=dc1_params,priors=dc1_priors,data='data/datavector_cosmosis.txt',cov='data/cov.npy',nofz='data/n_of_zs.txt')
    run.submit(label='dc1_sig8_Om_03',nodes=1,procs=32,hr=48,pts=1000,mneff=0.5,mntol=0.5,ia=False,pz=True,mbias=False,params=dc1_params,priors=dc1_priors,data='data/datavector_cosmosis.txt',cov='data/cov.npy',nofz='data/n_of_zs.txt')
    run.submit(label='dc1_sig8_Om_03',nodes=1,procs=32,hr=48,pts=1000,mneff=0.5,mntol=0.5,ia=False,pz=False,mbias=True,params=dc1_params,priors=dc1_priors,data='data/datavector_cosmosis.txt',cov='data/cov.npy',nofz='data/n_of_zs.txt')
    run.submit(label='dc1_sig8_Om_03',nodes=1,procs=32,hr=48,pts=1000,mneff=0.5,mntol=0.5,ia=False,pz=True,mbias=True,params=dc1_params,priors=dc1_priors,data='data/datavector_cosmosis.txt',cov='data/cov.npy',nofz='data/n_of_zs.txt')
    run.submit(label='dc1_sig8_Om_03',nodes=1,procs=32,hr=48,pts=1000,mneff=0.5,mntol=0.5,ia=True,pz=True,mbias=True,params=dc1_params,priors=dc1_priors,data='data/datavector_cosmosis.txt',cov='data/cov.npy',nofz='data/n_of_zs.txt')

    return

  @staticmethod
  def submit(label='dc1',nodes=1,procs=32,hr=48,pts=1000,mneff=0.5,mntol=0.5,ia=False,pz=False,mbias=False,planck=False,tomobins=3,params=dc1_params,priors=dc1_priors,data='data/datavector_cosmosis.txt',cov='data/cov.npy',nofz='data/n_of_zs.txt',cldir=''):

    from popen2 import popen2
    import subprocess as sp
    import time

    vary['sigma8_input']=True
    vary['omega_m']=True
    vary['h0']=True
    vary['ns']=True
    vary['tau']=True
    vary['omega_b']=True
    vary['w']=True
    vary['wa']=False

    name=label
    sout=label+'_pz-'+str(pz)+'_m-'+str(mbias)+'_ia-'+str(ia)+'_planck-'+str(planck)
    outfile='out/'+sout+'.txt' #'multinest-'+str(pts)+'_'+str(mneff)+'_'+str(mntol)+
    spriors=config.cosmosisrootdir+sout+'_priors.ini'
    sparams=config.cosmosisrootdir+sout+'_values.ini'


    mods=r' '
    like=r'wl '
    if pz:
      mods+=r'photoz_bias '
      for i in xrange(tomobins):
        vary['bias_'+str(i+1)]=True
        prior['bias_'+str(i+1)]=True
    if ia:
      sia=r'T'
      mods+=r'linear_alignment shear_shear add_intrinsic '
      vary['A']=True
    else:
      sia=r'F'
      mods+=r'shear_shear '
    if mbias:
      mods+=r'shear_m_bias '
      for i in xrange(tomobins):
        vary['m'+str(i+1)]=True
        prior['m'+str(i+1)]=True
    if planck:
      # vary['a_planck']=True
      # prior['a_planck']=True
      mods+=r'planck'
      like+=r'planck2015 '

    print vary
    print prior

    make.values(params,vary,ia,pz,mbias,planck,tomobins,sparams)
    if make.priors(priors,prior,ia,pz,mbias,planck,tomobins,spriors)==0:
      spriors=''    

    p = sp.Popen('qsub', shell=True, bufsize=1, stdin=sp.PIPE, stdout=sp.PIPE, close_fds=True, cwd=config.cosmosisrootdir)
    output,input = p.stdout, p.stdin

    job_string = """#!/bin/bash
    #PBS -l nodes=%s:ppn=%s
    #PBS -l walltime=%s:00:00
    #PBS -N %s
    #PBS -o %s.log
    #PBS -j oe
    #PBS -m abe 
    #PBS -M michael.troxel@manchester.ac.uk
    module use /home/zuntz/modules/module-files
    module load python
    module use /etc/modulefiles/
    cd /home/troxel/cosmosis/
    source my-source
    cd %s
    cd $PBS_O_WORKDIR
    export POINTS="%s"
    export MNEFF="%s"
    export MNTOL="%s"
    export OUTFILE="%s"
    export DATA="%s"
    export COV="%s"
    export NOFZ="%s"
    export IA="%s"
    export MODS="%s"
    export LIKE="%s"
    export PRIORS="%s"
    export PARAMS="%s"
    export CLDIR="%s"
    mpirun -n %s cosmosis --mpi data_in.ini
    postprocess -o plots -p %s %s""" % (str(nodes),str(procs),str(hr),name,sout,config.cosmosisrootdir,str(pts),str(mneff),str(mntol),outfile,data,cov,nofz,sia,mods,like,spriors,sparams,cldir,str(procs),sout,outfile)    

    output,outputerr=p.communicate(input=job_string)

    print job_string
    print output

    time.sleep(0.1)

    return

  @staticmethod
  def submit_pz_spec_test(pz0,test,boot=False,cosmo=False,nodes=1,procs=32,hr=48,params=dc1_params,priors=dc1_priors):

    from popen2 import popen2
    import subprocess as sp
    import time

    if pz0.pztype+'.txt' not in os.listdir(config.pzrootdir+test+'/nofz'):
      print 'Missing '+pz0.pztype+'.txt'
    if 'notomo_'+pz0.pztype+'.txt' not in os.listdir(config.pzrootdir+test+'/nofz'):
      print 'Missing '+'notomo_'+pz0.pztype+'.txt'
    if 'spec_'+pz0.pztype+'.txt' not in os.listdir(config.pzrootdir+test+'/nofz'):
      print 'Missing '+'spec_'+pz0.pztype+'.txt'
    if 'notomo_spec_'+pz0.pztype+'.txt' not in os.listdir(config.pzrootdir+test+'/nofz'):
      print 'Missing '+'notomo_spec_'+pz0.pztype+'.txt'
    if boot&hasattr(pz0,'boot'):
      for i in xrange(pz0.boot):
        if pz0.pztype+'_'+str(i)+'.txt' not in os.listdir(config.pzrootdir+test+'/nofz'):
          print 'Missing '+pz0.pztype+'_'+str(i)+'.txt'
        if 'notomo_'+pz0.pztype+'_'+str(i)+'.txt' not in os.listdir(config.pzrootdir+test+'/nofz'):
          print 'Missing '+'notomo_'+pz0.pztype+'_'+str(i)+'.txt'  

    vary['sigma8_input']=True
    params['sigma8_input']=(0.85, 1., 1.15)

    spriors=config.pzrootdir+test+'/priors.ini'
    sparams=config.pzrootdir+test+'/values.ini'

    make.values(params,False,False,False,sparams)

    jobstring0="""#!/bin/bash
    #PBS -l nodes=%s:ppn=%s
    #PBS -l walltime=%s:00:00
    #PBS -N %s
    #PBS -o %s.log
    #PBS -j oe
    #PBS -m abe 
    #PBS -M michael.troxel@manchester.ac.uk
    module use /home/zuntz/modules/module-files
    module load python
    module use /etc/modulefiles/
    cd /home/troxel/cosmosis/
    source my-source
    cd %s
    cd $PBS_O_WORKDIR
    export DIR="%s"
    """ % (str(nodes),str(procs),str(hr),test,test,config.pzrootdir,config.pzrootdir+test)

    if boot:

      jobstring3="""cosmosis %sdata_in.ini
      """ % (config.pzrootdir)

      if cosmo:
        jobstring3="""mpirun -n 32 cosmosis --mpi %sdata_in.ini
        """ % (config.pzrootdir)

        jobstring=jobstring0
        p = sp.Popen('qsub', shell=True, bufsize=1, stdin=sp.PIPE, stdout=sp.PIPE, close_fds=True, cwd=config.pzrootdir+test)
        output,input = p.stdout, p.stdin

        jobstring3="""cosmosis %sdata_in.ini
        """ % (config.pzrootdir)

        jobstring1="""export SAMP="grid"
        export LIKE="xipm"
        """
        jobstring1="""export SAMP="grid"
        export LIKE="wl"
        """
        ii=0
        for nofz in os.listdir(config.pzrootdir+test+'/nofz'):
          if (pz0.pztype+'.txt' not in nofz):
            if 'notomo' in nofz:
              continue
            ii+=1
            if ii in [10,20,30,40]:
              print jobstring
              output,outputerr=p.communicate(input=jobstring)
              time.sleep(0.1)
              jobstring=jobstring0
              p = sp.Popen('qsub', shell=True, bufsize=1, stdin=sp.PIPE, stdout=sp.PIPE, close_fds=True, cwd=config.pzrootdir+test)
              output,input = p.stdout, p.stdin
            #jobstring2="""export SAVEXI="save_xi"
            jobstring2="""export SAVEXI="cl_like"
            export POFZ="%s"
            export POFZ2="%s"
            """ % ('spec_'+pz0.pztype,nofz[:-4])
            jobstring+=jobstring1+jobstring2+jobstring3

            jobstring4="""postprocess %s/out/spec_%s_%s.txt -o %s/out/sim_data_%s
            """ % (config.pzrootdir+test,pz0.pztype,nofz[:-4],dir+test,nofz[:-4])
            jobstring+=jobstring4


        print jobstring
        output,outputerr=p.communicate(input=jobstring)
        time.sleep(0.1)

      else:

        jobstring=jobstring0
        p = sp.Popen('qsub', shell=True, bufsize=1, stdin=sp.PIPE, stdout=sp.PIPE, close_fds=True, cwd=config.pzrootdir+test)
        output,input = p.stdout, p.stdin

        jobstring1="""export SAMP="test"
        export LIKE=""
        """

        for nofz in os.listdir(config.pzrootdir+test+'/nofz'):
          if pz0.pztype+'.txt' not in nofz:
            #jobstring2="""export SAVEXI="save_xi"
            jobstring2="""export SAVEXI="generate_dataset"
            export POFZ="%s"
            export POFZ2="%s"
            """ % (nofz[:-4],nofz[:-4])
            jobstring+=jobstring1+jobstring2+jobstring3

        output,outputerr=p.communicate(input=jobstring)
        time.sleep(0.1)

    else:

      jobstring3="""cosmosis %sdata_in.ini
      """ % (config.pzrootdir)
 
      if cosmo:

        jobstring3="""mpirun -n 32 cosmosis --mpi %sdata_in.ini
        """ % (config.pzrootdir)

        jobstring=jobstring0
        p = sp.Popen('qsub', shell=True, bufsize=1, stdin=sp.PIPE, stdout=sp.PIPE, close_fds=True, cwd=config.pzrootdir+test)
        output,input = p.stdout, p.stdin
        jobstring3="""cosmosis %sdata_in.ini
        """ % (config.pzrootdir)

        jobstring1="""export SAMP="grid"
        export LIKE="xipm"
        """
        jobstring1="""export SAMP="grid"
        export LIKE="wl"
        """

        for nofz in os.listdir(config.pzrootdir+test+'/nofz'):
          if ('spec' not in nofz)&(pz0.pztype+'.txt' in nofz):
            if 'notomo' in nofz:
              #jobstring2="""export SAVEXI="xipm_2d"
              jobstring2="""export SAVEXI="cl_like"
              export POFZ="%s"
              export POFZ2="%s"
              """ % ('notomo_spec_'+pz0.pztype,'notomo_'+pz0.pztype)
              jobstring4="""postprocess %s/out/%s_%s.txt -o %s/out/sim_data_%s
              """ % (config.pzrootdir+test,'notomo_spec_'+pz0.pztype,'notomo_'+pz0.pztype,config.pzrootdir+test,'notomo_'+pz0.pztype)
            else:
              #jobstring2="""export SAVEXI="save_xi"
              jobstring2="""export SAVEXI="cl_like"
              export POFZ="%s"
              export POFZ2="%s"
              """ % ('spec_'+pz0.pztype,pz0.pztype)
              jobstring4="""postprocess %s/out/spec_%s_%s.txt -o %s/out/sim_data_%s
              """ % (config.pzrootdir+test,pz0.pztype,nofz[:-4],config.pzrootdir+test,nofz[:-4])
            jobstring+=jobstring1+jobstring2+jobstring3

            jobstring+=jobstring4 

        #print jobstring           

        output,outputerr=p.communicate(input=jobstring)
        time.sleep(0.1)

      else:

        p = sp.Popen('qsub', shell=True, bufsize=1, stdin=sp.PIPE, stdout=sp.PIPE, close_fds=True, cwd=config.pzrootdir+test)
        output,input = p.stdout, p.stdin
        jobstring=jobstring0

        jobstring1="""export SAMP="test"
        export LIKE=""
        """

        for nofz in os.listdir(config.pzrootdir+test+'/nofz'):
          if pz0.pztype+'.txt' in nofz:
            if 'notomo' in nofz:
              #jobstring2="""export SAVEXI="save_xi_2d"
              jobstring2="""export SAVEXI="generate_dataset"
              export POFZ="%s"
              export POFZ2="%s"
              """ % (nofz[:-4],nofz[:-4])
            else:
              #jobstring2="""export SAVEXI="save_xi"
              jobstring2="""export SAVEXI="generate_dataset"
              export POFZ="%s"
              export POFZ2="%s"
              """ % (nofz[:-4],nofz[:-4])
            jobstring+=jobstring1+jobstring2+jobstring3

        output,outputerr=p.communicate(input=jobstring)
        time.sleep(0.1)

    time.sleep(0.1)

    return

class make(object):

  @staticmethod
  def values(params,vary,ia,pz,mbias,planck,tomobins,sfile):

    n='\n'

    with open(sfile,'w') as f:
      f.write('[cosmological_parameters]'+n)
      for x in ['omega_m','h0','omega_b','sigma8_input','tau','n_s','A_s','omega_k','w','wa']:#,'omnuh2','massless_nu','massive_nu'
        if (vary.get(x)):
          f.write(x+' = '+str(params.get(x)[0])+' '+str(params.get(x)[1])+' '+str(params.get(x)[2])+n)
        else:
          f.write(x+' = '+str(params.get(x)[1])+n)
      if pz:
        f.write('\n[wl_photoz_errors]'+n)
        for x in ['bias_1','bias_2','bias_3','bias_4','bias_5','bias_6'][:tomobins]:
          if (vary.get(x)):
            f.write(x+' = '+str(params.get(x)[0])+' '+str(params.get(x)[1])+' '+str(params.get(x)[2])+n)
          else:
            f.write(x+' = '+str(params.get(x)[1])+n)
      if mbias:
        f.write('\n[shear_calibration_parameters]'+n)
        for x in ['m1','m2','m3','m4','m5','m6'][:tomobins]:
          if (vary.get(x)):
            f.write(x+' = '+str(params.get(x)[0])+' '+str(params.get(x)[1])+' '+str(params.get(x)[2])+n)
          else:
            f.write(x+' = '+str(params.get(x)[1])+n)
      if ia:
        f.write('\n[intrinsic_alignment_parameters]'+n)
        for x in ['A']:
          if (vary.get(x)):
            f.write(x+' = '+str(params.get(x)[0])+' '+str(params.get(x)[1])+' '+str(params.get(x)[2])+n)
          else:
            f.write(x+' = '+str(params.get(x)[1])+n)
        # f.write('\n[ia_z_field]'+n)
        # for x in ['alpha']:
        #   if (vary.get(x)):
        #     f.write(x+' = '+str(params.get(x)[0])+' '+str(params.get(x)[1])+' '+str(params.get(x)[2])+n)
        #   else:
        #     f.write(x+' = '+str(params.get(x)[1])+n)
      if planck:
        f.write('\n[planck]'+n)
        for x in ['a_planck']:
          if (vary.get(x)):
            f.write(x+' = '+str(params.get(x)[0])+' '+str(params.get(x)[1])+' '+str(params.get(x)[2])+n)
          else:
            f.write(x+' = '+str(params.get(x)[1])+n)

    return

  @staticmethod
  def priors(params,prior,ia,pz,mbias,planck,tomobins,sfile):

    n='\n'
    cnt=0

    with open(sfile,'w') as f:
      f.write('[cosmological_parameters]'+n)
      for x in ['omega_m','h0','omega_b','sigma8_input','tau','n_s','A_s','omega_k','w','wa']:#,'omnuh2','massless_nu','massive_nu'
        if (prior.get(x)):
          cnt+=1
          f.write(x+' = gaussian '+str(params.get(x)[0])+' '+str(params.get(x)[1])+n)
      if pz:
        f.write('\n[wl_photoz_errors]'+n)
        for x in ['bias_1','bias_2','bias_3','bias_4','bias_5','bias_6'][:tomobins]:
          if (prior.get(x)):
            cnt+=1
            f.write(x+' = gaussian '+str(params.get(x)[0])+' '+str(params.get(x)[1])+n)
      if mbias:
        f.write('\n[shear_calibration_parameters]'+n)
        for x in ['m1','m2','m3','m4','m5','m6'][:tomobins]:
          if (prior.get(x)):
            cnt+=1
            f.write(x+' = gaussian '+str(params.get(x)[0])+' '+str(params.get(x)[1])+n)
      if ia:
        f.write('\n[intrinsic_alignment_parameters]'+n)
        for x in ['A']:
          if (prior.get(x)):
            cnt+=1
            f.write(x+' = gaussian '+str(params.get(x)[0])+' '+str(params.get(x)[1])+n)
        # f.write('\n[ia_z_field]'+n)
        # for x in ['alpha']:
        #   if (prior.get(x)):
        #     cnt+=1
        #     f.write(x+' = gaussian '+str(params.get(x)[0])+' '+str(params.get(x)[1])+n)
      if planck:
        f.write('\n[planck]'+n)
        for x in ['a_planck']:
          if (prior.get(x)):
            cnt+=1
            f.write(x+' = gaussian '+str(params.get(x)[0])+' '+str(params.get(x)[1])+n)

    return cnt


  @staticmethod
  def nofz(pz0,test):

    if not os.path.exists(config.pzrootdir+test):
      os.makedirs(config.pzrootdir+test)
    if not os.path.exists(config.pzrootdir+test+'/nofz'):
      os.makedirs(config.pzrootdir+test+'/nofz')
    if not os.path.exists(config.pzrootdir+test+'/out'):
      os.makedirs(config.pzrootdir+test+'/out')

    out=np.vstack((pz0.bin,pz0.pz[0,:]))
    out=np.row_stack(([0.,0.],out.T,[pz0.bin[-1]+(pz0.bin[1]-pz0.bin[0]),0.]))
    np.savetxt(config.pzrootdir+test+'/nofz/notomo_'+pz0.pztype+'.txt',out)

    out=np.vstack((pz0.bin,[pz0.pz[i+1,:] for i in range(pz0.tomo-1)]))
    out=np.row_stack(([0. for i in range(pz0.tomo)],out.T,np.append(pz0.bin[-1]+(pz0.bin[1]-pz0.bin[0]),[0. for i in range(pz0.tomo-1)])))
    np.savetxt(config.pzrootdir+test+'/nofz/'+pz0.pztype+'.txt',out)

    if hasattr(pz0, 'spec'):

      out=np.vstack((pz0.bin,pz0.spec[0,:]))
      out=np.row_stack(([0.,0.],out.T,[pz0.bin[-1]+(pz0.bin[1]-pz0.bin[0]),0.]))
      np.savetxt(config.pzrootdir+test+'/nofz/notomo_spec_'+pz0.pztype+'.txt',out)

      out=np.vstack((pz0.bin,[pz0.spec[i+1,:] for i in range(pz0.tomo-1)]))
      out=np.row_stack(([0. for i in range(pz0.tomo)],out.T,np.append(pz0.bin[-1]+(pz0.bin[1]-pz0.bin[0]),[0. for i in range(pz0.tomo-1)])))
      np.savetxt(config.pzrootdir+test+'/nofz/spec_'+pz0.pztype+'.txt',out)

    if hasattr(pz0,'boot'):

      for j in xrange(pz0.boot):

        out=np.vstack((pz0.bin,pz0.bootspec[0,j,:]))
        out=np.row_stack(([0.,0.],out.T,[pz0.bin[-1]+(pz0.bin[1]-pz0.bin[0]),0.]))
        np.savetxt(config.pzrootdir+test+'/nofz/notomo_'+pz0.pztype+'_'+str(j)+'.txt',out)

        out=np.vstack((pz0.bin,[pz0.bootspec[i+1,j,:] for i in range(pz0.tomo-1)]))
        out=np.row_stack(([0. for i in range(pz0.tomo)],out.T,np.append(pz0.bin[-1]+(pz0.bin[1]-pz0.bin[0]),[0. for i in range(pz0.tomo-1)])))
        np.savetxt(config.pzrootdir+test+'/nofz/'+pz0.pztype+'_'+str(j)+'.txt',out)

    return

# class stat_methods(object):

#   @staticmethod
#   def loop_output2():

#     import matplotlib
#     matplotlib.use ('agg')
#     import matplotlib.pyplot as plt
#     plt.style.use('/home/troxel/SVA1/SVA1StyleSheet.mplstyle')

#     x=[0.,0.01,0.03,0.05]
#     frac=np.zeros((4,5,4))
#     width2=np.zeros((4,5,4))
#     frac[:,:,0]=1.
#     width=[0.01144225,0.06967396,0.04952199,0.01559282]
#     for i in xrange(5):
#       width2[:,i,0]=width
#     print width
#     ind=-1
#     for a in ['True','False']:
#       for b in ['True','False']:
#         for c in ['True','False']:
#           if (a=='False')|(b=='False'):
#             if c=='True':
#               continue
#           y=[1.]
#           y2=[0.01144225]
#           ind+=1
#           for i,ii in enumerate(['_01','_03','_05']):
#             s='S50,f8,f8'
#             label='dc1_sig8'
#             tmp=np.genfromtxt('/home/troxel/cosmosis/cosmosis-des-library/wl/y1prep/plots/'+label+ii+'_pz-'+a+'_m-'+b+'_ia-'+c+'_means.txt',dtype=s)
#             for j in range(len(tmp)):
#               for k,param in enumerate(['sigma8']):
#                 if param in tmp[j][0]:
#                   print label,k,param,ii,a,b,c,np.around(tmp[j][1],3),np.around(tmp[j][2],3),np.around(tmp[j][2]/tmp[j][1],3),np.around(tmp[j][2]/width[0],3)
#                   y=np.append(y,tmp[j][2]/width[0])
#                   y2=np.append(y2,tmp[j][2])
#                   print ind
#                   frac[0,ind,i+1]=tmp[j][2]/width[0]
#                   width2[0,ind,i+1]=tmp[j][2]
#                   # print label,k,param,ii,a,b,c,tmp[j][2]

#     for k,param in enumerate(['sigma8','omega_m','s8']):
#       ind=-1
#       for a in ['True','False']:
#         for b in ['True','False']:
#           for c in ['True','False']:
#             if (a=='False')|(b=='False'):
#               if c=='True':
#                 continue
#             ind+=1
#             # if (a=='True')&(b=='True')&(c=='True'):
#             #   frac[k+1,ind,:]=frac[k+1,ind,0]
#             #   width2[k+1,ind,:]=width2[k+1,ind,0]
#             #   continue
#             for i,ii in enumerate(['_01','_03','_05']):
#               s='S50,f8,f8'
#               label='dc1_sig8_Om'
#               tmp=np.genfromtxt('/home/troxel/cosmosis/cosmosis-des-library/wl/y1prep/plots/'+label+ii+'_pz-'+a+'_m-'+b+'_ia-'+c+'_means.txt',dtype=s)
#               print '/home/troxel/cosmosis/cosmosis-des-library/wl/y1prep/plots/'+label+ii+'_pz-'+a+'_m-'+b+'_ia-'+c+'_means.txt'
#               for j in range(len(tmp)):
#                 print j,i,a,b,c,param,tmp[j][0]
#                 if param in tmp[j][0]:
#                   print label,k,param,ii,a,b,c,np.around(tmp[j][1],3),np.around(tmp[j][2],3),np.around(tmp[j][2]/tmp[j][1],3),np.around(tmp[j][2]/width[k+1],3)
#                   frac[k+1,ind,i+1]=tmp[j][2]/width[k+1]
#                   width2[k+1,ind,i+1]=tmp[j][2]

#                   # print label,k,param,ii,a,b,c,tmp[j][2]
#                   #print label,k,param,ii,a,b,c,tmp[j][2]
#           #         width[k,i+1]=tmp[j][2]
#           #         frac[k,i+1]=tmp[j][2]/width[k+1,0]
#           # print width,frac

#           # plt.plot(x,sig868)
#           # plt.savefig(label+i+'_pz-'+a+'_m-'+b+'_ia-'+c+'_sig8_68.png')
#           # plt.close()
#           # plt.plot(x,fsig868)
#           # plt.savefig(label+i+'_pz-'+a+'_m-'+b+'_ia-'+c+'_sig8_f68.png')
#           # plt.close()
#           # plt.plot(x,ratsig868)
#           # if (a=='True')&(b=='True')&(c=='True'):
#           #   plt.title('Photo-z + shear bias + IA')
#           # if (a=='True')&(b=='True')&(c=='False'):
#           #   plt.title('Photo-z + shear bias')
#           # if (a=='False')&(b=='True')&(c=='False'):
#           #   plt.title('shear bias')
#           # if (a=='True')&(b=='False')&(c=='False'):
#           #   plt.title('Photo-z bias')
#           # plt.xlabel('Percent prior in pz, m bias')
#           # plt.ylabel('Ratio to width with no systematics')
#           # plt.ylim([1.,3.])
#           # plt.xlim([0.,0.06])
#           # plt.savefig(label+i+'_pz-'+a+'_m-'+b+'_ia-'+c+'_sig8_rat68.png')
#           # plt.close()
#     print frac
#     print width2

#     for i,param in enumerate(['sigma8','sigma8b','omegam','s8']):
#       plt.plot(x,frac[i,0,:],label='Photo-z+shear+IA')
#       plt.plot(x,frac[i,1,:],label='Photo-z+shear')
#       plt.plot(x,frac[i,2,:],label='Photo-z')
#       plt.plot(x,frac[i,3,:],label='Shear')
#       plt.plot(x,frac[i,4,:],label='None')
#       plt.axhline(y=1.5,linestyle='--')
#       plt.title(param)
#       plt.xlabel('Percent prior')
#       plt.ylabel('Fractional increase in 1-sigma width')
#       plt.legend(loc='upper left',ncol=1, frameon=True,prop={'size':12})
#       plt.savefig('req_frac_'+param+'.png')
#       plt.close()
#     for i,param in enumerate(['sigma8','sigma8b','omegam','s8']):
#       plt.plot(x,width2[i,0,:],label='Photo-z+shear+IA')
#       plt.plot(x,width2[i,1,:],label='Photo-z+shear')
#       plt.plot(x,width2[i,2,:],label='Photo-z')
#       plt.plot(x,width2[i,3,:],label='Shear')
#       plt.plot(x,width2[i,4,:],label='None')
#       plt.axhline(y=1.5*width2[i,0,0],linestyle='--')
#       plt.title(param)
#       plt.xlabel('Percent prior')
#       plt.ylabel('Width of 1-sigma constraint')
#       plt.legend(loc='upper left',ncol=1, frameon=True,prop={'size':12})
#       plt.savefig('req_width_'+param+'.png')
#       plt.close()

#     return

#   @staticmethod
#   def loop_output(label):

#     import matplotlib
#     matplotlib.use ('agg')
#     import matplotlib.pyplot as plt
#     plt.style.use('/home/troxel/SVA1/SVA1StyleSheet.mplstyle')

#     for a in ['True','False']:
#       for b in ['True','False']:
#         for c in ['True','False']:
#           sig868=[]#[0.011]
#           fsig868=[]#[0.014]
#           ratsig868=[]#[1.]
#           # if (a=='True')&(b=='True')&(c=='True'):
#           #   continue
#           if (a=='False')|(b=='False'):
#             if c=='True':
#               continue
#           for i in ['_01','_03','_05']:#,'_02','_03','_05']:
#             s='S50,f8,f8'
#             tmp=np.genfromtxt('/home/troxel/cosmosis/cosmosis-des-library/wl/y1prep/plots/'+label+i+'_pz-'+a+'_m-'+b+'_ia-'+c+'_means.txt',dtype=s)
#             for j in range(len(tmp)):
#               if 'sigma8' in tmp[j][0]:
#                 print i,a,b,c,np.around(tmp[j][1],3),np.around(tmp[j][2],3),np.around(tmp[j][2]/tmp[j][1],3),np.around(tmp[j][2]/0.08074807,3)
#                 sig868=np.append(sig868,tmp[j][2])
#                 fsig868=np.append(fsig868,tmp[j][2]/tmp[j][1])
#                 ratsig868=np.append(ratsig868,tmp[j][2]/0.08074807)#01160847
#                 #print tmp[j][1],tmp[j][2]

#           x=[.01, .03,.05]#,.02,.03,.05]
#           #print label+i+'_pz-'+a+'_m-'+b+'_ia-'+c,x,sig868,fsig868,ratsig868

#           plt.plot(x,sig868)
#           plt.savefig(label+i+'_pz-'+a+'_m-'+b+'_ia-'+c+'_sig8_68.png')
#           plt.close()
#           plt.plot(x,fsig868)
#           plt.savefig(label+i+'_pz-'+a+'_m-'+b+'_ia-'+c+'_sig8_f68.png')
#           plt.close()
#           plt.plot(x,ratsig868)
#           if (a=='True')&(b=='True')&(c=='True'):
#             plt.title('Photo-z + shear bias + IA')
#           if (a=='True')&(b=='True')&(c=='False'):
#             plt.title('Photo-z + shear bias')
#           if (a=='False')&(b=='True')&(c=='False'):
#             plt.title('shear bias')
#           if (a=='True')&(b=='False')&(c=='False'):
#             plt.title('Photo-z bias')
#           plt.xlabel('Percent prior in pz, m bias')
#           plt.ylabel('Ratio to width with no systematics')
#           plt.ylim([1.,3.])
#           plt.xlim([0.,0.06])
#           plt.savefig(label+i+'_pz-'+a+'_m-'+b+'_ia-'+c+'_sig8_rat68.png')
#           plt.close()

#     return

#   @staticmethod
#   def loop_pp():

#     from popen2 import popen2
#     import subprocess as sp
#     import time

#     for a in ['True','False']:
#       for b in ['True','False']:
#         for c in ['True','False']:
#           if (a=='False')|(b=='False'):
#             if c=='True':
#               continue
#           s='postprocess -o plots --no-fill --derive derive_s8.py -p comp_s8'+'_pz-'+a+'_m-'+b+'_ia-'+c
#           print s
#           continue
#           for i in ['_05','_03','_01']:

#             s+='out/dc1_sig8'+i+'_pz-'+a+'_m-'+b+'_ia-'+c+'.txt '

#           print s

#     return


#   @staticmethod
#   def get_bounds(file,nparams):


#     def smooth_likelihood(x):

#         import kde

#         #Interpolate using KDE from cosmosis
#         kde = kde.KDE(x, factor=2.0)
#         x_axis, like = kde.grid_evaluate(100, (x.min(), x.max()) )
#         return x_axis, like      

#     tmp=np.loadtxt(file)

#     for i in range(nparams):

#       l0=tmp[:,i+1]
#       x0=tmp[:,0]

#       # x,l=smooth_likelihood(x0)

#       # level1, level2=stat_methods.find_grid_contours(l, 0.68, 0.95)
#       # print level1,level2

#       # above = np.where(l>level1)[0]
#       # left = above[0]
#       # right = above[-1]
#       # x1a = x[left]
#       # x2a = x[right]

#       # above = np.where(l>level2)[0]
#       # left = above[0]
#       # right = above[-1]
#       # x1b = x[left]
#       # x2b = x[right]

#       # print i,x1a,x2a,x1b,x2b

#       order=np.argsort(x0)
#       x0=x0[order]
#       l0=l0[order]

#       x = np.linspace(x0[0], x0[-1], 10*len(x0))
#       l = np.interp(x, x0, l0)
#       l -= l.max()

#       level1, level2=stat_methods.find_grid_contours(l, 0.68, 0.95)
#       #print level1,level2

#       above = np.where(l>level1)[0]
#       left = above[0]
#       right = above[-1]
#       x1a = x[left]
#       x2a = x[right]

#       above = np.where(l>level2)[0]
#       left = above[0]
#       right = above[-1]
#       x1b = x[left]
#       x2b = x[right]

#       #print i,x1a,x2a,x1b,x2b

#     return x1a,(x2a+x1a)/2.,x2a


#   @staticmethod
#   def find_grid_contours(log_like, contour1, contour2):
#     import scipy.optimize

#     #from cosmosis

#     like = np.exp(log_like)
#     like_total = like.sum()
#     def objective(limit, target):
#       w = np.where(like>limit)
#       return like[w].sum() - target
#     target1 = like_total*contour1
#     target2 = like_total*contour2
#     level1 = scipy.optimize.bisect(objective, like.min(), like.max(), args=(target1,))
#     level2 = scipy.optimize.bisect(objective, like.min(), like.max(), args=(target2,))
#     level1 = np.log(level1)
#     level2 = np.log(level2)

#     return level1, level2
