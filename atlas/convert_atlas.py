import nrrd
from scipy.io import savemat

annotation, header = nrrd.read("atlas/annotation_25.nrrd")

print("Atlas shape:", annotation.shape)
print("Atlas datatype:", annotation.dtype)

savemat(
    "atlas/annotationVolume.mat",
    {"annotationVolume": annotation}
)