"""
PREREQUISITE
svxconv must be available on your system.
For example, via an alias that defined the necessary classpaths.
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
import os.path as op, shutil
import nibabel,numpy
import scipy.misc
import math

"""
In the code below, the transpose is needed because images are written with the 
first dimension as rows, second as columns. This must be flipped to align with
the common definition of x- and y axes .
The ::-1 part is a mirroring operation on the y-axis, which is needed because
images are written top to bottom, reversed from the common y-axis direction.
"""
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

class Nifti2Stl(FancyModule):
  inputs = odict((
    ('niifile', {'type': assertFile}),
    ('stlfile', {'type': str}),
    ('slicerange', {'type': assertList}),
    ('labelinclude', {'type': assertList,'default':'*'}),
    ('labelexclude', {'type': assertList,'default':[0]}),
    ('labelrange', {'type': assertList,'default':None})
  ))
  def main(self,niifile,stlfile,slicerange,labelinclude,labelexclude,labelrange):
    nii = nibabel.load(niifile)
    nii = nibabel.as_closest_canonical(nii)
    hdr = nii.get_header()
    img = nii.get_data()
    step = slicerange[2] if 2 in slicerange else 1
    pngdir = self.tempdir('png')
    svgdir = self.tempdir('svg')
    svxdir = self.tempdir('svx')
    densdir = self.tempdir('svx/dens')
    i0 = slicerange[0]
    for i in range(i0,slicerange[1],step):
      slc = get_slice(img,1,i)
      if labelrange:
        maskIncl = numpy.logical_and(slc>labelrange[0],slc<labelrange[1]).astype(numpy.uint8)*255
      else:
        if labelinclude == '*':
          maskIncl = numpy.ones(slc.shape).astype(numpy.uint8)*255
        else:
          maskIncl = numpy.in1d(slc, labelinclude).reshape(slc.shape).astype(numpy.uint8)*255
      maskExcl = numpy.in1d(slc, labelexclude).reshape(slc.shape).astype(numpy.uint8)*255
      mask = numpy.logical_and(maskIncl,numpy.logical_not(maskExcl))
      pngfile = op.join(pngdir,'slice{:04d}.png'.format(i-i0+1))
      scipy.misc.imsave(pngfile,mask)
      svgfile = op.join(svgdir,'slice{:04d}.svg'.format(i-i0+1))
      FancyExec.fromParent(self).setInput(
        *['/my/github/Vectorization-of-brain-atlases/mindthegap/bin/mindthegap',
          '-i',pngfile,
          '-o',svgfile,
          '-t','1.5'
        ]
      ).run()
      densfile = op.join(densdir,'slice{:04d}.png'.format(i-i0+1))      
      if i>i0:
        with open(svgfile,'r') as fp:
          svg = fp.read().replace('fill="#FFFFFF" stroke-width="0"','stroke-width="4px"')
        with open(svgfile,'w') as fp:
          fp.write(svg)
        FancyExec.fromParent(self).setInput(
          *['rsvg',
            svgfile,
            densfile
          ]
        ).run()
      else:
        shutil.copyfile(pngfile,densfile)
      
    slc = get_slice(img,1,1)
    mask = numpy.zeros(slc.shape).astype(numpy.uint8)
    scipy.misc.imsave(op.join(densdir,'slice{:04d}.png'.format(0)),mask)
    scipy.misc.imsave(op.join(densdir,'slice{:04d}.png'.format(i-i0+2)),mask)
    
    with open(op.join(svxdir,'manifest.xml'),'w') as fp:
      numSlices = math.floor((slicerange[1]-slicerange[0])/step)+2
      fp.write('<grid gridSizeX="{}" gridSizeY="{}" gridSizeZ="{}" voxelSize="{}" subvoxelBits="8" originX="0" originY="0" originZ="0" slicesOrientation="Y">'.format(int(img.shape[0]),int(numSlices),int(img.shape[2]),0.0001))
      fp.write('<channels><channel type="DENSITY" slices="dens/slice%04d.png" bits="8"/></channels>')
      fp.write('</grid>')
    
    svxfile = self.tempfile('nifti2stl.svx')
    FancyExec.fromParent(self).setInput(
      *['zip','-rD',svxfile,'.'],
      __workdir__ = svxdir
    ).run()

    stlfile = self.tempfile('nifti2stl.stl')
    FancyExec.fromParent(self).setInput(
      *['svxconv','-input',svxfile,'-output',stlfile,'-voxelSize','0.0002']
    ).run()

    return FancyOutput(
      svxfile = svxfile
    )
#endclass

if __name__ == '__main__':
  Nifti2Stl.fromCommandLine().run()


