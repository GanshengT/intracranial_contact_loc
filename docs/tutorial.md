# intracranial_contact_loc
[![DOI](https://zenodo.org/badge/891756132.svg)](https://doi.org/10.5281/zenodo.14217838)

This repository includes algorithms in MATLAB for localizing intracranial contact. CT images, MRI images, and freesurfer software are required.

Updates:
2025 Aug 19, I developed an automated pipeline for estimating intracranial electrode location based on CT images and manufacture information of the electrode lead. The main script has been tested on DIXI electrodes.

## Installation
```bash
git clone https://github.com/GanshengT/intracranial_contact_loc.git
cd intracranial_contact_loc
```

## Tutorial

### Normalize contacts to MNI or another reference brain
After localizing contacts, users can map the contact table to a reference anatomical image with the local plugin in `automated_pipeline/mni_normalization_plugin`.

```matlab
addpath(genpath('automated_pipeline/mni_normalization_plugin'));
reference_mri = fullfile('automated_pipeline', 'mni_normalization_plugin', ...
    'example_data', 'reference', 'mni_icbm152_t1_tal_nlin_asym_09c.nii.gz');
[contact_tbl, paths] = gt_normalize_contacts_to_reference_brain(subject_mri, contact_table_mat, reference_mri, output_dir);
```

The packaged example reference is FreeSurfer's MNI152 NIfTI file, copied from the FreeSurfer installation. You may use another reference brain image if it is the desired registration target. The plugin writes a warped subject MRI, an affine transform, a reference-space contact CSV, and FreeView `.dat` contact files. MATLAB, ANTs, a subject MRI, a contact table, and a reference brain image are required.

The example uses `automated_pipeline/mni_normalization_plugin/example_data/ants/bin` as the packaged ANTs binary folder. These binaries are copied from a local ANTs installation and can be replaced by setting `ants_bin_path` to another compatible ANTs `bin` directory.
