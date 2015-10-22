import numpy as np

import catalog

class write_methods(object):

  @staticmethod
  def write_append(line,cat,label='',create=False,newline=True):

    if newline:
      n='\n'
      print line
    else:
      n=''
      print line,

    if create:
      with open('text/'+label+'_'+cat.name+'.txt','w') as f:
        f.write(line+n)
    else:
      with open('text/'+label+'_'+cat.name+'.txt','a') as f:
        f.write(line+n)

    return


  @staticmethod
  def heading(line,cat,label='',create=False):

    write_methods.write_append('',cat,label=label,create=create)
    write_methods.write_append('------------------',cat,label=label,create=create)
    write_methods.write_append('   '+line+'   ',cat,label=label,create=False)
    write_methods.write_append('------------------',cat,label=label,create=False)
    write_methods.write_append('',cat,label=label,create=create)

    return

  @staticmethod
  def get_file(cat,label=''):

    return 'text/'+label+'_'+cat.name+'.txt'


