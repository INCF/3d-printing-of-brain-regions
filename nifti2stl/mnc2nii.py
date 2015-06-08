from fancypipe import *
import nibabel,numpy

class Mnc2Nii(FancyModule):
  inputs = odict((
    ('mncfile', {'type': assertFile}),
    ('niifile', {'type': assertMatch('\.nii(\.gz)?$',fromStart=False,decompose=False)})
  ))
  def main(self,mncfile,niifile):
    mnc = nibabel.load(mncfile)
    hdr = mnc.get_header()
    q = hdr.get_best_affine()
    print('Affine transformation matrix: {}'.format(q))
    img = mnc.get_data()
    print('Some data... {}'.format(img[:10,:10,:10]))
    print('Niifile {}'.format(niifile))
    nii = nibabel.Nifti1Image(img,q)
    nii = nibabel.as_closest_canonical(nii)
    nibabel.save(nii,niifile)
    return FancyOutput( niifile )
#endclass

if __name__ == '__main__':
  Mnc2Nii.fromCommandLine().run()
