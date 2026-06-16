# intracranial_contact_loc
[![DOI](https://zenodo.org/badge/891756132.svg)](https://doi.org/10.5281/zenodo.14217838)

Welcome to the documentation. Use the navigation for tutorials and examples. This repository includes algorithms in MATLAB for localizing intracranial contact. CT images, elecrtode information (OR sheet) and freesurfer software are required. MRI images are optional and for later coregistration and brain segmentation, and elecrtode labeling.

- [Tutorials](tutorial.md)
- [Examples](examples.md)
- [API Reference](api.md)
API is under reorganization. Examples are currently included in the tutorials.

## Reference-brain normalization plugin
`automated_pipeline/mni_normalization_plugin` provides a local MATLAB plugin for mapping localized contacts to MNI or another reference brain. The minimal example uses FreeSurfer's MNI152 NIfTI reference:

```matlab
addpath(genpath('automated_pipeline/mni_normalization_plugin'));
reference_mri = fullfile('automated_pipeline', 'mni_normalization_plugin', ...
    'example_data', 'reference', 'mni_icbm152_t1_tal_nlin_asym_09c.nii.gz');
[contact_tbl, paths] = gt_normalize_contacts_to_reference_brain(subject_mri, contact_table_mat, reference_mri, output_dir);
```

The packaged reference is copied from FreeSurfer. Users can replace it with any compatible MNI template, atlas-space brain, or study-specific reference image.

The example also includes a packaged ANTs runtime folder at `automated_pipeline/mni_normalization_plugin/example_data/ants/bin`. Users may replace this path with their own ANTs installation.

Developers and maintainers:
Gansheng Tan (g.tan@wustl.edu)

Updates:
2025 Aug 19, I developed an automated pipeline for estimating intracranial electrode location based on CT images and manufacture information of the electrode lead. The main script has been tested on DIXI electrodes.

## License
Distributed under the [MIT License](../LICENSE).
