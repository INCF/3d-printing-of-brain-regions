from fancypipe import *
import os,os.path as op, shutil
import nibabel,numpy
import scipy.misc
import math
import nifti2stl

class Hippo3d(FancyTask):
  inputs = odict((
    ('mnipath', {'type': assertDir}),
    ('resultspath', {'type': assertDir})
  ))
  def main(self,mnipath,resultspath):
    nii = nibabel.load(op.join(resultspath,'VU','probmap0_hippoR.nii.gz'))
    img = nii.get_data()
    nii = nibabel.load(op.join(resultspath,'VU','probmap0_hippoL.nii.gz'))
    img += nii.get_data()
    mask = nifti2stl.applyMask(img,[0.5,1.0],None,None)
    select = numpy.flatnonzero(mask.max(axis=0).max(axis=1))
    sliceStart = select[0]
    sliceEnd = select[-1]
    midSlice = int((sliceStart+sliceEnd)/2) 

    stlfiles = {}
    for D in [0,1]:
      fullrange = [sliceStart,midSlice] if D==0 else [midSlice+1,sliceEnd]
      # context: contour + grey matter
      contour = nifti2stl.Nifti2Mask.fromParent(self).setInput(
        niifile = op.join(mnipath,'mni_icbm152_t1_tal_nlin_asym_09c_mask.nii.gz'),
        slicerange = fullrange,
        strokepx = 1.1,
        bottompx = -1 if D==0 else 1
      )
      grey = nifti2stl.Nifti2Mask.fromParent(self).setInput(
        niifile = op.join(mnipath,'mni_icbm152_gm_tal_nlin_asym_09c.nii.gz'),
        labelrange = [0.5,1.0],
        slicerange = [midSlice-8,midSlice-1] if D==0 else [midSlice+2,midSlice+9],
        strokepx = 1.1,
        bottompx = None
      )
      context = nifti2stl.CombineMasks.fromParent(self).setInput(
        contour.requestOutput('mask'),
        contour.requestOutput('slicerange'),
        grey.requestOutput('mask'),
        grey.requestOutput('slicerange')
      )
      context_stl = nifti2stl.Mask2Stl.fromParent(self).setInput(
        mask = context.requestOutput('mask'),
        slicerange = context.requestOutput('slicerange'),
        stlfile = self.tempfile('underlay{}.stl'.format(D))
      )

      # hippocampus: left + right
      hippoL = nifti2stl.Nifti2Mask.fromParent(self).setInput(
        niifile = op.join(resultspath,'VU','probmap0_hippoL.nii.gz'),
        labelrange = [0.5,1.0],
        slicerange = fullrange
      )
      hippoR = nifti2stl.Nifti2Mask.fromParent(self).setInput(
        niifile = op.join(resultspath,'VU','probmap0_hippoR.nii.gz'),
        labelrange = [0.5,1.0],
        slicerange = fullrange
      )
      hippo = nifti2stl.CombineMasks.fromParent(self).setInput(
        hippoL.requestOutput('mask'),
        hippoL.requestOutput('slicerange'),
        hippoR.requestOutput('mask'),
        hippoR.requestOutput('slicerange')
      )      
      hippo_stl = nifti2stl.Mask2Stl.fromParent(self).setInput(
        mask = hippo.requestOutput('mask'),
        slicerange = hippo.requestOutput('slicerange'),
        stlfile = self.tempfile('hippo{}.stl'.format(D))
      )

      stlfiles['context_stl_{}'.format(D)] = context_stl.requestOutput('stlfile')
      stlfiles['hippo_stl_{}'.format(D)] = hippo_stl.requestOutput('stlfile')

    return FancyOutput(
      **stlfiles
    )
#endclass

if __name__ == '__main__':
  Hippo3d.fromCommandLine().run()
