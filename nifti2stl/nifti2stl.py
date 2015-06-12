"""
PREREQUISITES
svxconv and mindthegap must be available on your system.
For example, via an alias that defines the necessary classpaths.
alias svxconv="env CLASSPATH=/my/github/AbFab3D/apps/svxconv/classes:/my/github/AbFab3D/classes:/my/github/AbFab3D/lib/commons-io-2.4.jar java svxconv.SVXConv"
"""

"""
TODO
svxconv step does not work, configuration issue.
convert to svg to create thick region outlines, then rasterize again
using rsvglib. Aliasing no problem here.
To print the hippocampus:
- find start and end slices
- take the mean
- from start to mean: print mean + outline in gold, hippocampus in red
- from mean to end: print mean + outline in gold, hippocampus in red
"""

from fancypipe import *
import os,os.path as op, shutil
import nibabel,numpy
import scipy.misc
import scipy.ndimage as ndimage
import math

"""
In the code below, the transpose is needed because images are written with the 
first dimension as rows, second as columns. This must be flipped to align with
the common definition of x- and y axes .
The ::-1 part is a mirroring operation on the y-axis, which is needed because
images are written top to bottom, reversed from the common y-axis direction.
def get_slice(img,dim,i):
    if dim == 0:
        slc = img[i,:,::-1].squeeze();
    elif dim == 1:
        slc = img[:,i,::-1].squeeze();
    elif dim == 2:
        slc = img[:,::-1,i].squeeze();
    else:
        raise Exception('Cannot return a slice for dimension "{}"'.format(dim))
    slc = slc.swapaxes(0,1)
    return slc
"""

def applyMask(A,labelrange,labelinclude,labelexclude):
  if labelrange:
    maskIncl = numpy.logical_and(A>=labelrange[0],A<=labelrange[1])
  else:
    if labelinclude == '*':
      maskIncl = numpy.ones(A.shape)
    else:
      maskIncl = numpy.in1d(A, labelinclude).reshape(A.shape)
  maskExcl = numpy.in1d(A, labelexclude).reshape(A.shape)
  mask = numpy.logical_and(maskIncl,numpy.logical_not(maskExcl))
  return mask

def ball3d(radius):
  r = int(math.ceil(radius+0.5))
  ball = numpy.zeros([2*r-1,2*r-1,2*r-1],float)
  subsample0 = [0.05,0.15,0.25,0.35,0.45]
  subsample = [-v for v in subsample0]
  subsample.extend(subsample0)
  c = r-1
  for i in range(0,r):
    for j in range(0,r):
      for k in range(0,r):
        subi = subsample if i>0 else subsample0
        for di in subi:
          subj = subsample if j>0 else subsample0
          for dj in subj:
             subk = subsample if k>0 else subsample0
             for dk in subk:
               if (i+di)**2+(j+dj)**2+(k+dk)**2<=radius*radius:
                 ball[c-i,c-j,c-k] += 1
                 ball[c-i,c-j,c+k] += 1
                 ball[c-i,c+j,c-k] += 1
                 ball[c-i,c+j,c+k] += 1
                 ball[c+i,c-j,c-k] += 1
                 ball[c+i,c-j,c+k] += 1
                 ball[c+i,c+j,c-k] += 1
                 ball[c+i,c+j,c+k] += 1
  ball /= len(subsample)**3
  return ball
  
def binaryContour(A):
  x3d = numpy.zeros([3,3,3])
  x3d[1,1,:] = 1
  x3d[:,1,1] = 1
  x3d[1,:,1] = 1
  return A-ndimage.minimum_filter(A,footprint=x3d).astype(numpy.uint8)
  
class CombineMasks(FancyTask):
  def main(self,mask1,slicerange1,mask2,slicerange2):
    mn1 = slicerange1[0]
    mx1 = slicerange1[1]
    mn2 = slicerange2[0]
    mx2 = slicerange2[1]
    mn = min(mn1,mn2)
    mx = max(mx1,mx2)
    mask = numpy.zeros([mask1.shape[0],mx-mn+1,mask1.shape[2]])
    mask[:,mn1-mn:mx1-mn,:] = mask1
    mask[:,mn2-mn:mx2-mn,:] = numpy.maximum(mask[:,mn2-mn:mx2-mn,:],mask2)
    return FancyOutput(
      mask = mask,
      slicerange = [mn,mx]
    )

class Svx2Stl(FancyExec):
  inputs = odict((
    ('svxfile', {'type': assertFile}),
    ('stlfile', {'default': None}),
  ))
  AbFab3dRoot = '/my/github/AbFab3D'
  SVXENV = dict(os.environ, CLASSPATH='{0}/apps/svxconv/classes:{0}/classes:{0}/lib/commons-io-2.4.jar:{0}/lib/vecmath.jar'.format(AbFab3dRoot))

  def beforeMain(self,svxfile,stlfile):
    self.setEnv(self.SVXENV)
    if not stlfile: 
      stlfile = self.tempfile('result.stl')
      self.setInput(stlfile = stlfile)
    return FancyInput(
      *['/usr/bin/java','svxconv.SVXConv'],
      **{
        '-input': svxfile,
        '-output': stlfile,
        '-voxelSize':'0.0001'
      }
    )
#endclass

class Mask2Stl(FancyTask):
  inputs = odict((
    ('svxfile', {'default': None}),
    ('stlfile', {'default': None}),
  ))
  def main(self,mask,slicerange,svxfile,stlfile):
    svxdir = self.tempdir('svx')
    densdir = self.tempdir('svx','dens')
    for i in range(0,mask.shape[1]):
      M = mask[:,i,::-1].squeeze()
      densfile = op.join(densdir,'slice{:05d}.png'.format(i+1))
      scipy.misc.imsave(densfile,M)
   
    sz = [mask.shape[0],mask.shape[2]]
    bg = numpy.zeros(sz).astype(numpy.uint8)
    scipy.misc.imsave(op.join(densdir,'slice{:05d}.png'.format(0)),bg)
    scipy.misc.imsave(op.join(densdir,'slice{:05d}.png'.format(mask.shape[1]+1)),bg)
    
    with open(op.join(svxdir,'manifest.xml'),'w') as fp:
      fp.write('<grid gridSizeX="{}" gridSizeY="{}" gridSizeZ="{}" voxelSize="{}" subvoxelBits="8" originX="0" originY="{}" originZ="0" slicesOrientation="Y">'.format(sz[0],int(mask.shape[1]+2),sz[1],0.0001,slicerange[0]*0.0001))
      fp.write('<channels><channel type="DENSITY" slices="dens/slice%05d.png" bits="8"/></channels>')
      fp.write('</grid>')
    
    svxfile = self.tempfile('nifti2stl.svx')
    dozip = FancyExec.fromParent(self).setCwd(svxdir).setInput(
      *['zip','-rD',svxfile,'.']
    )
    
    svx2stl = Svx2Stl.fromParent(self).setInput(
      svxfile = dozip.requestOutput(2),
      stlfile = stlfile
    )

    return FancyOutput(
      svxfile = dozip.requestOutput(2),
      stlfile = svx2stl.requestOutput('stlfile')
    )
#endclass

class Nifti2Mask(FancyTask):
  inputs = odict((
    ('niifile', {'type': assertFile}),
    ('slicerange', {'type': assertList,'default':[]}),
    ('labelinclude', {'type': assertList,'default':'*'}),
    ('labelexclude', {'type': assertList,'default':[0]}),
    ('labelrange', {'type': assertList,'default':None}),
    ('bottompx', {'default':None}),
    ('strokepx', {'default':None})
  ))
  def main(self,niifile,slicerange,labelinclude,labelexclude,labelrange,bottompx,strokepx):
    nii = nibabel.load(niifile)
    nii = nibabel.as_closest_canonical(nii)
    hdr = nii.get_header()
    img = nii.get_data()
    pngdir = self.tempdir('png')
    svgdir = self.tempdir('svg')
    svxdir = self.tempdir('svx')
    densdir = self.tempdir('svx/dens')

    # autodetect slicerange if not specified
    if slicerange == []:
      mask = applyMask(img,labelrange,labelinclude,labelexclude)
      select = numpy.flatnonzero(mask.max(axis=0).max(axis=1))
      slicerange = [select[0],select[-1]]

    srange = range(slicerange[0],slicerange[1])
    mask = img[:,srange,:]
    mask = applyMask(mask,labelrange,labelinclude,labelexclude).astype(float)
    if strokepx:
      if bottompx:
        bottompx = range(0,bottompx) if bottompx>0 else range(mask.shape[1]+bottompx,mask.shape[1])
        solidmask = mask[:,bottompx,:]
      mask = binaryContour(mask)
      if bottompx:
        mask[:,bottompx,:] = solidmask
      radius = strokepx*math.pi/3
      ball = ball3d(radius)-1
      mask = ndimage.grey_dilation(mask,structure=ball)
    
    mask = (mask*255.999).astype(numpy.uint8)
    return FancyOutput(
      mask = mask,
      slicerange = slicerange
    )
#endclass

class Nifti2Stl(FancyTask):
  inputs = odict(Nifti2Mask.inputs)
  inputs['stlfile'] = {'type': str}
  
  def main(self,niifile,slicerange,labelinclude,labelexclude,labelrange,bottompx,strokepx,stlfile):
    nifti2mask = Nifti2Mask.fromParent(self).setInput(
      niifile,slicerange,labelinclude,labelexclude,labelrange,bottompx,strokepx
    )
    mask2stl = Mask2Stl.fromParent(self).setInput(
      mask = nifti2mask.requestOutput('mask'),
      slicerange = nifti2mask.requestOutput('slicerange'),
      stlfile = stlfile
    )
    result = FancyOutput(
      slicerange = slicerange
    )
    if 'mask' in self.requests: 
      result['mask'] = nifti2mask.requestOutput('mask')
    if 'svxfile' in self.requests:
      result['svxfile'] = mask2stl.requestOutput('svxfile')
    if 'stlfile' in self.requests:
      result['stlfile'] = mask2stl.requestOutput('stlfile')
    return result
#endclass

if __name__ == '__main__':
  Nifti2Stl.fromCommandLine().run()
