# intracranial_contact_loc
[![DOI](https://zenodo.org/badge/891756132.svg)](https://doi.org/10.5281/zenodo.14217838)
[![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-blue.svg)](https://ganshengt.github.io/intracranial_contact_loc/)

This repository includes algorithms in MATLAB for localizing intracranial contact. CT images, MRI images, and freesurfer software are required.

Updates:
2025 Aug 19, I developed an automated pipeline for estimating intracranial electrode location based on CT images and manufacture information of the electrode lead. The main script has been tested on DIXI electrodes.

## Installation
```bash
git clone https://github.com/GanshengT/intracranial_contact_loc.git
cd intracranial_contact_loc
```

## Motivation and objective
**Electrode localization** refers to identifying the precise anatomical positions of implanted intracranial electrodes.  
This step is essential in both **clinical practice** and **neuroscience research**:

- **Clinical applications**:  
  - Mapping eloquent cortex for surgical planning.  
  - Delineating seizure onset zones in patients with epilepsy.  
  - Guiding resection or stimulation strategies with minimal risk.  

- **Research applications**:  
  - Linking electrophysiological signals to anatomical structures.  
  - Understanding network-level dynamics of cognition and disease.  
  - Enabling reproducible and anatomically grounded scientific discovery.  

Precise and **objective localization** is thus critical. However, most workflows still rely heavily on **manual localization**, which presents two major limitations:  

1. **Subjectivity** – human raters must distinguish contacts on CT scans that often suffer from streak artifacts and noisy intensity. For example, an isosurface in FreeSurfer may appear as a continuous lead, making individual contacts ambiguous.  
2. **Time-consuming and error-prone** – manual tracing and labeling are labor-intensive, and mistakes can propagate through analysis pipelines.  

To address these issues, **v1.0** of this pipeline introduced an automated approach that leverages CT intensity and a connected-graph searching algorithm to identify electrode contacts.  

Building on this foundation, **v2.0** incorporates **electrode manufacturing specifications** to further improve robustness. By exploiting the known geometry of linear electrode leads, the algorithm can fit multiple contacts simultaneously, mitigating the impact of noisy voxels and imaging artifacts. This allows faster, more accurate, and reproducible electrode localization for both clinical and research use cases.

## Methods
### step 1: metal identification
Our automated electrode localization approach begins by analyzing the CT intensity distribution.  
We identify peaks corresponding to metal contacts, which allows us to select candidate voxels associated with the electrode lead.

<p align="center">
  <img src="figs/BJH079_CT_intensity_distribution_electrode_threshold.png" alt="CT intensity distribution showing electrode threshold" width="600">
</p>

*Figure: Example CT intensity distribution. We find peaks with prominance higher than the noise. The right peak correspond to metal voxels*  

### step 2: lead identification
We visualize the **metal voxel isosurface** to sanity-check the threshold and the spatial distribution of metal (contacts + lead body). We also overlay the **planned trajectory** (exported from ROSA) and restrict analysis to voxels that lie within a **cylindrical neighborhood** of this line segment. This (i) suppresses unrelated metal (e.g., screws/plates) and (ii) disambiguates adjacent leads.

**Note on ROSA coordinates:** Ensure the ROSA trajectory is in **the same RAS frame** as the CT used for localization. Co registration is needed

<p align="center">
  <video src="figs/BJH079_ct_rotation_traj_R_and_metal_voxels.mp4" controls width="700">
    Your browser does not support the video tag. <a href="figs/BJH079_ct_rotation_traj_R_and_metal_voxels.mp4">Download the MP4</a>.
  </video>
</p>

<p align="center">
  <img src="figs/BJH079_cylinder_around_planned_traj_lead_R.png" alt="cylinder around the planned trajectory for localizing CT voxels corresponding to this lead" width="600">
</p>

*Figure: Extracted metal voxels and the planned trajectory for one lead*  

### step 3


## Tutorial
Please use this link (https://ganshengt.github.io/intracranial_contact_loc/) for detailed documentation.


