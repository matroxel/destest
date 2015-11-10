import numpy as np
import os

import catalog
import config
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
  'bias_1' : (-0.1, 0., 0.1),
  'bias_2' : (-0.1, 0., 0.1),
  'bias_3' : (-0.1, 0., 0.1),
  'bias_4' : (-0.1, 0., 0.1),
  'bias_5' : (-0.1, 0., 0.1),
  'bias_6' : (-0.1, 0., 0.1),
  'm1' : (-0.1, 0., 0.1),
  'm2' : (-0.1, 0., 0.1),
  'm3' : (-0.1, 0., 0.1),
  'm4' : (-0.1, 0., 0.1),
  'm5' : (-0.1, 0., 0.1),
  'm6' : (-0.1, 0., 0.1),
  'A' : (-3., 1., 5.0),
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
  def submit(label='dc1',nodes=1,procs=32,hr=48,pts=200,mneff=0.8,mntol=0.5,ia=False,pz=False,mbias=False,planck=False,tomobins=3,params=dc1_params,priors=dc1_priors,data='data/datavector_cosmosis.txt',cov='data/cov.npy',nofz='data/n_of_zs.txt',cldir='',resume=False,submit=True):
    """
    A wrapper to submit cosmosis runs. Currently works for my PBS environment. Will add alternate option to just print necessary commands to a file to run as desired.

    Use:

    ....


    """

    from popen2 import popen2
    import subprocess as sp
    import time

    vary['sigma8_input']=True
    vary['omega_m']=True
    vary['h0']=True
    vary['n_s']=True
    vary['tau']=True
    vary['omega_b']=True
    vary['w']=True
    vary['wa']=False

    name=label
    sout=label+'_pz-'+str(pz)+'_m-'+str(mbias)+'_ia-'+str(ia)+'_planck-'+str(planck)
    outfile='out/'+sout+'.txt' #'multinest-'+str(pts)+'_'+str(mneff)+'_'+str(mntol)+
    mnoutfile='out/'+sout+'.multinest'
    spriors=config.cosomsiscosmodir+sout+'_priors.ini'
    sparams=config.cosomsiscosmodir+sout+'_values.ini'


    mods=r' '
    like=r'xipm '
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

    p = sp.Popen('qsub', shell=True, bufsize=1, stdin=sp.PIPE, stdout=sp.PIPE, close_fds=True, cwd=config.cosomsiscosmodir)
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
    export MNOUT="%s"
    export MNRESUME="%s"
    mpirun -n %s cosmosis --mpi data_in.ini
    postprocess -o plots -p %s %s""" % (str(nodes),str(procs),str(hr),name,sout,config.cosomsiscosmodir,str(pts),str(mneff),str(mntol),outfile,data,cov,nofz,sia,mods,like,spriors,sparams,cldir,mnoutfile,str(resume),str(procs),sout,outfile)    

    output,outputerr=p.communicate(input=job_string)

    print job_string
    print output

    time.sleep(0.1)

    return

  @staticmethod
  def submit_pz_spec_test(pz0,test,boot=False,cosmo=False,nodes=1,procs=32,hr=48,params=dc1_params,priors=dc1_priors,submit=True):
    """
    A wrapper to submit cosmosis runs specifically for photo-z spec validation. Currently works for my PBS environment. Will add alternate option to just print necessary commands to a file to run as desired. Needs work. Currently works with Cls due to need of synthetic covariance, but will switch back to xi if used in WL analysis once covariances are available. Could also make it optional which (xi vs cl) to use.

    Use:

    ....


    """

    if submit:
      from popen2 import popen2
      import subprocess as sp
      import time

    if pz0.pztype+'.txt' not in os.listdir(config.pztestdir+test+'/nofz'):
      print 'Missing '+pz0.pztype+'.txt'
    if 'notomo_'+pz0.pztype+'.txt' not in os.listdir(config.pztestdir+test+'/nofz'):
      print 'Missing '+'notomo_'+pz0.pztype+'.txt'
    if 'spec_'+pz0.pztype+'.txt' not in os.listdir(config.pztestdir+test+'/nofz'):
      print 'Missing '+'spec_'+pz0.pztype+'.txt'
    if 'notomo_spec_'+pz0.pztype+'.txt' not in os.listdir(config.pztestdir+test+'/nofz'):
      print 'Missing '+'notomo_spec_'+pz0.pztype+'.txt'
    if boot&hasattr(pz0,'boot'):
      for i in xrange(pz0.boot):
        if pz0.pztype+'_'+str(i)+'.txt' not in os.listdir(config.pztestdir+test+'/nofz'):
          print 'Missing '+pz0.pztype+'_'+str(i)+'.txt'
        if 'notomo_'+pz0.pztype+'_'+str(i)+'.txt' not in os.listdir(config.pztestdir+test+'/nofz'):
          print 'Missing '+'notomo_'+pz0.pztype+'_'+str(i)+'.txt'  

    vary['sigma8_input']=True
    params['sigma8_input']=(0.85, 1., 1.15)

    spriors=config.pztestdir+test+'/priors.ini'
    sparams=config.pztestdir+test+'/values.ini'

    make.values(params,vary,False,False,False,False,0,sparams)

    if submit:
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
      """ % (str(nodes),str(procs),str(hr),test,test,config.pztestdir,config.pztestdir+test)
    else:
      jobstring0="""#!/bin/bash
      cd %s
      %s
      cd %s
      """ % (config.cosmosisrootdir,config.cosmosissource,config.pztestdir)


    if boot:

      jobstring3="""cosmosis %sdata_in.ini
      """ % (config.pztestdir)

      if cosmo:
        jobstring3="""mpirun -n 32 cosmosis --mpi %sdata_in.ini
        """ % (config.pztestdir)

        jobstring=jobstring0
        if submit:
          p = sp.Popen('qsub', shell=True, bufsize=1, stdin=sp.PIPE, stdout=sp.PIPE, close_fds=True, cwd=config.pztestdir+test)
          output,input = p.stdout, p.stdin

        jobstring3="""cosmosis %sdata_in.ini
        """ % (config.pztestdir)

        jobstring1="""export SAMP="grid"
        export LIKE="xipm"
        """
        jobstring1="""export SAMP="grid"
        export LIKE="wl"
        """
        ii=0
        for nofz in os.listdir(config.pztestdir+test+'/nofz'):
          if (pz0.pztype+'.txt' not in nofz):
            if 'notomo' in nofz:
              continue
            ii+=1
            if ii in [10,20,30,40]:
              print jobstring
              output,outputerr=p.communicate(input=jobstring)
              time.sleep(0.1)
              jobstring=jobstring0
              p = sp.Popen('qsub', shell=True, bufsize=1, stdin=sp.PIPE, stdout=sp.PIPE, close_fds=True, cwd=config.pztestdir+test)
              output,input = p.stdout, p.stdin
            #jobstring2="""export SAVEXI="save_xi"
            jobstring2="""export SAVEXI="cl_like"
            export POFZ="%s"
            export POFZ2="%s"
            """ % ('spec_'+pz0.pztype,nofz[:-4])
            jobstring+=jobstring1+jobstring2+jobstring3

            jobstring4="""postprocess %s/out/spec_%s_%s.txt -o %s/out/sim_data_%s
            """ % (config.pztestdir+test,pz0.pztype,nofz[:-4],dir+test,nofz[:-4])
            jobstring+=jobstring4

        if submit:
          output,outputerr=p.communicate(input=jobstring)
          time.sleep(0.1)
        else:
          with open('cosmosis_pz_boot-'+str(boot)+'_cosmo-'+str(cosmo)+'.submit','w') as f:
            f.write(jobstring)

      else:

        jobstring=jobstring0
        if submit:
          p = sp.Popen('qsub', shell=True, bufsize=1, stdin=sp.PIPE, stdout=sp.PIPE, close_fds=True, cwd=config.pztestdir+test)
          output,input = p.stdout, p.stdin

        jobstring1="""export SAMP="test"
        export LIKE=""
        """

        for nofz in os.listdir(config.pztestdir+test+'/nofz'):
          if pz0.pztype+'.txt' not in nofz:
            #jobstring2="""export SAVEXI="save_xi"
            jobstring2="""export SAVEXI="generate_dataset"
            export POFZ="%s"
            export POFZ2="%s"
            """ % (nofz[:-4],nofz[:-4])
            jobstring+=jobstring1+jobstring2+jobstring3

        if submit:
          output,outputerr=p.communicate(input=jobstring)
          time.sleep(0.1)
        else:
          with open('cosmosis_pz_boot-'+str(boot)+'_cosmo-'+str(cosmo)+'.submit','w') as f:
            f.write(jobstring)        

    else:

      jobstring3="""cosmosis %sdata_in.ini
      """ % (config.pztestdir)
 
      if cosmo:

        jobstring3="""mpirun -n 32 cosmosis --mpi %sdata_in.ini
        """ % (config.pztestdir)

        jobstring=jobstring0
        if submit:
          p = sp.Popen('qsub', shell=True, bufsize=1, stdin=sp.PIPE, stdout=sp.PIPE, close_fds=True, cwd=config.pztestdir+test)
          output,input = p.stdout, p.stdin

        jobstring3="""cosmosis %sdata_in.ini
        """ % (config.pztestdir)

        jobstring1="""export SAMP="grid"
        export LIKE="xipm"
        """
        jobstring1="""export SAMP="grid"
        export LIKE="wl"
        """

        for nofz in os.listdir(config.pztestdir+test+'/nofz'):
          if ('spec' not in nofz)&(pz0.pztype+'.txt' in nofz):
            if 'notomo' in nofz:
              #jobstring2="""export SAVEXI="xipm_2d"
              jobstring2="""export SAVEXI="cl_like"
              export POFZ="%s"
              export POFZ2="%s"
              """ % ('notomo_spec_'+pz0.pztype,'notomo_'+pz0.pztype)
              jobstring4="""postprocess %s/out/%s_%s.txt -o %s/out/sim_data_%s
              """ % (config.pztestdir+test,'notomo_spec_'+pz0.pztype,'notomo_'+pz0.pztype,config.pztestdir+test,'notomo_'+pz0.pztype)
            else:
              #jobstring2="""export SAVEXI="save_xi"
              jobstring2="""export SAVEXI="cl_like"
              export POFZ="%s"
              export POFZ2="%s"
              """ % ('spec_'+pz0.pztype,pz0.pztype)
              jobstring4="""postprocess %s/out/spec_%s_%s.txt -o %s/out/sim_data_%s
              """ % (config.pztestdir+test,pz0.pztype,nofz[:-4],config.pztestdir+test,nofz[:-4])
            jobstring+=jobstring1+jobstring2+jobstring3

            jobstring+=jobstring4 


        if submit:
          output,outputerr=p.communicate(input=jobstring)
          time.sleep(0.1)
        else:
          with open('cosmosis_pz_boot-'+str(boot)+'_cosmo-'+str(cosmo)+'.submit','w') as f:
            f.write(jobstring)          

      else:

        jobstring=jobstring0
        if submit:
          p = sp.Popen('qsub', shell=True, bufsize=1, stdin=sp.PIPE, stdout=sp.PIPE, close_fds=True, cwd=config.pztestdir+test)
          output,input = p.stdout, p.stdin

        jobstring1="""export SAMP="test"
        export LIKE=""
        """

        for nofz in os.listdir(config.pztestdir+test+'/nofz'):
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

        if submit:
          output,outputerr=p.communicate(input=jobstring)
          time.sleep(0.1)
        else:
          with open('cosmosis_pz_boot-'+str(boot)+'_cosmo-'+str(cosmo)+'.submit','w') as f:
            f.write(jobstring)          

    return

class make(object):

  @staticmethod
  def values(params,vary,ia,pz,mbias,planck,tomobins,sfile):
    """
    Writes values.ini files for cosmosis run submitted via run class.
    """

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
    """
    Writes priors.ini files for cosmosis run submitted via run class.
    """

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
    """
    Writes n(z) files in cosmosis format for photo-z spec validation testing.
    """

    if not os.path.exists(config.pztestdir+test):
      os.makedirs(config.pztestdir+test)
    if not os.path.exists(config.pztestdir+test+'/nofz'):
      os.makedirs(config.pztestdir+test+'/nofz')
    if not os.path.exists(config.pztestdir+test+'/out'):
      os.makedirs(config.pztestdir+test+'/out')

    out=np.vstack((pz0.bin,pz0.pz[0,:]))
    out=np.row_stack(([0.,0.],out.T,[pz0.bin[-1]+(pz0.bin[1]-pz0.bin[0]),0.]))
    np.savetxt(config.pztestdir+test+'/nofz/notomo_'+pz0.pztype+'.txt',out)

    out=np.vstack((pz0.bin,[pz0.pz[i+1,:] for i in range(pz0.tomo-1)]))
    out=np.row_stack(([0. for i in range(pz0.tomo)],out.T,np.append(pz0.bin[-1]+(pz0.bin[1]-pz0.bin[0]),[0. for i in range(pz0.tomo-1)])))
    np.savetxt(config.pztestdir+test+'/nofz/'+pz0.pztype+'.txt',out)

    if hasattr(pz0, 'spec'):

      out=np.vstack((pz0.bin,pz0.spec[0,:]))
      out=np.row_stack(([0.,0.],out.T,[pz0.bin[-1]+(pz0.bin[1]-pz0.bin[0]),0.]))
      np.savetxt(config.pztestdir+test+'/nofz/notomo_spec_'+pz0.pztype+'.txt',out)

      out=np.vstack((pz0.bin,[pz0.spec[i+1,:] for i in range(pz0.tomo-1)]))
      out=np.row_stack(([0. for i in range(pz0.tomo)],out.T,np.append(pz0.bin[-1]+(pz0.bin[1]-pz0.bin[0]),[0. for i in range(pz0.tomo-1)])))
      np.savetxt(config.pztestdir+test+'/nofz/spec_'+pz0.pztype+'.txt',out)

    if hasattr(pz0,'boot'):

      for j in xrange(pz0.boot):

        out=np.vstack((pz0.bin,pz0.bootspec[0,j,:]))
        out=np.row_stack(([0.,0.],out.T,[pz0.bin[-1]+(pz0.bin[1]-pz0.bin[0]),0.]))
        np.savetxt(config.pztestdir+test+'/nofz/notomo_'+pz0.pztype+'_'+str(j)+'.txt',out)

        out=np.vstack((pz0.bin,[pz0.bootspec[i+1,j,:] for i in range(pz0.tomo-1)]))
        out=np.row_stack(([0. for i in range(pz0.tomo)],out.T,np.append(pz0.bin[-1]+(pz0.bin[1]-pz0.bin[0]),[0. for i in range(pz0.tomo-1)])))
        np.savetxt(config.pztestdir+test+'/nofz/'+pz0.pztype+'_'+str(j)+'.txt',out)

    return
