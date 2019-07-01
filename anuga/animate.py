"""
A module to allow interactive plotting in a Jupyter notebook of quantities and mesh 
associated with an ANUGA domain and SWW file.
"""

import numpy as np
import matplotlib.pyplot as plt
import os

class Domain_plotter:
  """
  A class to wrap ANUGA domain centroid values for stage, height, elevation
  xmomentunm and ymomentum, and triangulation information.
  """
  
  def __init__(self, domain, plot_dir = '_plot'):
    
    self.plot_dir = plot_dir
    self.make_plot_dir()
    
    import matplotlib.tri as tri
    self.nodes = domain.nodes
    self.triangles = domain.triangles
    self.x = domain.nodes[:,0]
    self.y = domain.nodes[:,1]
    
    self.xc = domain.centroid_coordinates[:,0]
    self.yc = domain.centroid_coordinates[:,1]

    self.xllcorner = domain.geo_reference.xllcorner
    self.yllcorner = domain.geo_reference.yllcorner
    self.zone = domain.geo_reference.zone
    
    self.triang = tri.Triangulation(self.x, self.y, self.triangles)
    
    self.elev  = domain.quantities['elevation'].centroid_values
    self.depth = domain.quantities['height'].centroid_values
    self.stage = domain.quantities['stage'].centroid_values
    self.xmom  = domain.quantities['xmomentum'].centroid_values
    self.ymom  = domain.quantities['ymomentum'].centroid_values
    self.domain = domain
    
  def _depth_frame(self, figsize, dpi):
 
    name = self.domain.get_name()
    time = self.domain.get_time() 

    fig = plt.figure(figsize=figsize, dpi=dpi)

    plt.title('Time {0:0>4}'.format(time))
    
    self.triang.set_mask(self.depth>0.01)
    plt.tripcolor(self.triang, 
              facecolors = self.elev,
              cmap='Greys_r')
    
    self.triang.set_mask(self.depth<0.01)
    plt.tripcolor(self.triang, 
              facecolors = self.depth,
              cmap='viridis')

    plt.colorbar()
    
    return    
    
  def save_depth_frame(self):

    figsize=(10,6)
    dpi = 80
    plot_dir = self.plot_dir
    name = self.domain.get_name()
    time = self.domain.get_time()

    self._depth_frame(figsize,dpi);
    
    if plot_dir is None:
        plt.savefig(name+'_{0:0>10}.png'.format(int(time)))
    else:
        plt.savefig(os.path.join(plot_dir, name+'_{0:0>10}.png'.format(int(time))))
    plt.close()
    
    return    

  def plot_depth_frame(self):
  
    figsize=(5,3)
    dpi = 80
    
    self._depth_frame(figsize,dpi)
    
    plt.show()
    
    return

  def make_depth_animation(self):
    import numpy as np
    import glob
    from matplotlib import image, animation
    from matplotlib import pyplot as plt

    plot_dir = self.plot_dir
    name = self.domain.get_name()
    time = self.domain.get_time() 
    
    if plot_dir is None:
        expression = name+'_*.png'
    else:
        expression = os.path.join(plot_dir, name+'_*.png')
    img_files = sorted(glob.glob(expression))

    figsize=(10,6)

    fig = plt.figure(figsize=figsize, dpi=80)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis('off')  # so there's not a second set of axes
    im = plt.imshow(image.imread(img_files[0]))

    def init():
      im.set_data(image.imread(img_files[0]))
      return im,

    def animate(i):
      image_i=image.imread(img_files[i])
      im.set_data(image_i)
      return im,

    anim = animation.FuncAnimation(fig, animate, init_func=init,
                            frames=len(img_files), interval=200, blit=True)

    plt.close()
  
    return anim

  def make_plot_dir(self, clobber=True):
    """
    Utility function to create a directory for storing a sequence of plot
    files, or if the directory already exists, clear out any old plots.  
    If clobber==False then it will abort instead of deleting existing files.
    """

    plot_dir = self.plot_dir
    if plot_dir is None:
      return
    else:
      import os
      if os.path.isdir(plot_dir):
          if clobber:
              os.system("rm %s/*" % plot_dir)
          else:
              raise IOError('*** Cannot clobber existing directory %s' % plot_dir)
      else:
          os.system("mkdir %s" % plot_dir)
      print "Figure files for each frame will be stored in " + plot_dir


      
      
   
  
      
      
class SWW_plotter:
  """
  A class to wrap ANUGA swwfile centroid values for stage, height, elevation
  xmomentunm and ymomentum, and triangulation information.
  """
  
  def __init__(self, swwfile = 'domain.sww', plot_dir = '_plot', minimum_allowed_depth=1.0e-03):
    
    self.plot_dir = plot_dir
    self.make_plot_dir()
    
    import matplotlib.tri as tri
    import numpy as np
 
    import os
    self.name = os.path.splitext(swwfile)[0]

    from anuga.file.netcdf import NetCDFFile
    p = NetCDFFile(swwfile)
    

    self.x = np.array(p.variables['x'])
    self.y = np.array(p.variables['y'])
    self.triangles = np.array(p.variables['volumes'])

    vols0=self.triangles[:,0]
    vols1=self.triangles[:,1]
    vols2=self.triangles[:,2]
    
    self.triang = tri.Triangulation(self.x, self.y, self.triangles)
    
    
    self.xc=(self.x[vols0]+self.x[vols1]+self.x[vols2])/3.0
    self.yc=(self.y[vols0]+self.y[vols1]+self.y[vols2])/3.0
    
    
    self.xllcorner = p.xllcorner
    self.yllcorner = p.yllcorner
    self.zone = p.zone
    
    self.elev  = np.array(p.variables['elevation_c'])
    self.stage = np.array(p.variables['stage_c'])
    self.xmom  = np.array(p.variables['xmomentum_c'])
    self.ymom  = np.array(p.variables['ymomentum_c'])
    
    self.depth = np.zeros_like(self.stage)
    if(len(self.elev.shape)==2):
      self.depth = self.stage - self.elev
    else:
      for i in range(self.depth.shape[0]):
        self.depth[i,:] = self.stage[i,:]-self.elev
    
    
    self.xvel  = np.where(self.depth > minimum_allowed_depth, self.xmom / self.depth, 0.0)
    self.yvel  = np.where(self.depth > minimum_allowed_depth, self.ymom / self.depth, 0.0)
    
    self.speed = np.sqrt(self.xvel**2 + self.yvel**2)
    
    self.time  = np.array(p.variables['time'])
    
    
  def _depth_frame(self, figsize, dpi, frame):
 
    name = self.name
    time = self.time[frame] 
    depth = self.depth[frame,:]
    try:
      elev  = self.elev[frame,:]
    except:
      elev  = self.elev

    ims = []
    
    fig = plt.figure(figsize=figsize, dpi=dpi)

    plt.title('Depth: Time {0:0>4}'.format(time))
    
    self.triang.set_mask(depth>0.01)
    plt.tripcolor(self.triang, 
              facecolors = elev,
              cmap='Greys_r')
    
    
    self.triang.set_mask(depth<0.01)
    plt.tripcolor(self.triang, 
              facecolors = depth,
              cmap='viridis')

    plt.colorbar()
    
    return  
    
  def save_depth_frame(self, frame=-1):

    figsize=(10,6)
    dpi = 160
    name = self.name
    time = self.time[frame] 
    plot_dir = self.plot_dir

    self._depth_frame(figsize,dpi,frame);
    
    if plot_dir is None:
        plt.savefig(name+'_depth_{0:0>10}.png'.format(int(time)))
    else:
        plt.savefig(os.path.join(plot_dir, name+'_depth_{0:0>10}.png'.format(int(time))))
    plt.close()
    
    return    

  def plot_depth_frame(self, frame=-1):
  
    figsize=(5,3)
    dpi = 80
    
    self._depth_frame(figsize,dpi,frame)
    
    plt.show()
    
    return

  def _stage_frame(self, figsize, dpi, frame):
 
    name = self.name
    time = self.time[frame] 
    stage = self.stage[frame,:]
    depth = self.depth[frame,:]
    try:
      elev  = self.elev[frame,:]
    except:
      elev  = self.elev

    ims = []
    
    fig = plt.figure(figsize=figsize, dpi=dpi)

    plt.title('Stage: Time {0:0>4}'.format(time))
    
    self.triang.set_mask(depth>0.01)
    plt.tripcolor(self.triang, 
              facecolors = elev,
              cmap='Greys_r')
    
    
    self.triang.set_mask(depth<0.01)
    plt.tripcolor(self.triang, 
              facecolors = stage,
              cmap='viridis')

    plt.colorbar()
    
    return  
    
  def save_stage_frame(self, frame=-1):

    figsize=(10,6)
    dpi = 160
    name = self.name
    time = self.time[frame] 
    plot_dir = self.plot_dir

    self._stage_frame(figsize,dpi,frame);
    
    if plot_dir is None:
        plt.savefig(name+'_stage_{0:0>10}.png'.format(int(time)))
    else:
        plt.savefig(os.path.join(plot_dir, name+'_stage_{0:0>10}.png'.format(int(time))))
    plt.close()
    
    return    

  def plot_stage_frame(self, frame=-1):
  
    figsize=(5,3)
    dpi = 80
    
    self._stage_frame(figsize,dpi,frame)
    
    plt.show()
    
    return

  def _speed_frame(self, figsize, dpi, frame):
 
    name = self.name
    time = self.time[frame] 
    depth = self.depth[frame,:]
    try:
      elev  = self.elev[frame,:]
    except:
      elev  = self.elev    
    speed = self.speed[frame,:]

    ims = []
    
    fig = plt.figure(figsize=figsize, dpi=dpi)

    plt.title('Speed: Time {0:0>4}'.format(time))
    
    self.triang.set_mask(depth>0.01)
    plt.tripcolor(self.triang, 
              facecolors = elev,
              cmap='Greys_r')
    
    
    self.triang.set_mask(depth<0.01)
    plt.tripcolor(self.triang, 
              facecolors = speed,
              cmap='viridis')

    plt.colorbar()
    
    return  
    
  def save_speed_frame(self, frame=-1):

    figsize=(10,6)
    dpi = 160
    name = self.name
    time = self.time[frame] 
    plot_dir = self.plot_dir

    self._speed_frame(figsize,dpi,frame);
    
    if plot_dir is None:
        plt.savefig(name+'_speed_{0:0>10}.png'.format(int(time)))
    else:
        plt.savefig(os.path.join(plot_dir, name+'_speed_{0:0>10}.png'.format(int(time))))
    plt.close()
    
    return    

  def plot_speed_frame(self, frame=-1):
  
    figsize=(5,3)
    dpi = 80
    
    self._speed_frame(figsize,dpi,frame)
    
    plt.show()
    
    return

  def make_depth_animation(self):
    
    return self._make_quantity_animation(quantity='depth')
  
  def make_speed_animation(self):
    
    return self._make_quantity_animation(quantity='speed')
  
  def make_stage_animation(self):
    
    return self._make_quantity_animation(quantity='stage')
  
  def _make_quantity_animation(self, quantity='depth'):
    import numpy as np
    import glob
    from matplotlib import image, animation
    from matplotlib import pyplot as plt

    plot_dir = self.plot_dir
    name = self.name
    
    if plot_dir is None:
        expression = name+'_'+quantity+'_*.png'
    else:
        expression = os.path.join(plot_dir, name+'_'+quantity+'_*.png')
    img_files = sorted(glob.glob(expression))

    figsize=(10,6)

    fig = plt.figure(figsize=figsize, dpi=80)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis('off')  # so there's not a second set of axes
    im = plt.imshow(image.imread(img_files[0]))

    def init():
      im.set_data(image.imread(img_files[0]))
      return im,

    def animate(i):
      image_i=image.imread(img_files[i])
      im.set_data(image_i)
      return im,

    anim = animation.FuncAnimation(fig, animate, init_func=init,
                            frames=len(img_files), interval=200, blit=True)

    plt.close()
  
    return anim
  
 
  def make_plot_dir(self, clobber=True):
    """
    Utility function to create a directory for storing a sequence of plot
    files, or if the directory already exists, clear out any old plots.  
    If clobber==False then it will abort instead of deleting existing files.
    """

    plot_dir = self.plot_dir
    if plot_dir is None:
      return
    else:
      import os
      if os.path.isdir(plot_dir):
          if clobber:
              os.system("rm %s/*" % plot_dir)
          else:
              raise IOError('*** Cannot clobber existing directory %s' % plot_dir)
      else:
          os.system("mkdir %s" % plot_dir)
      print "Figure files for each frame will be stored in " + plot_dir

  




