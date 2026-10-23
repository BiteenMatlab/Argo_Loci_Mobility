# Agrobacterium Chromosome Locus Mobility Analysis

Code for processing microscopy locus tracking movies, filtering trajectories, calculating mean squared displacement (MSD), and comparing locus mobility across *Agrobacterium* strains using Bayesian hierarchical analysis.

The tracking workflow uses MATLAB Live Scripts (`.mlx`) together with the Python segmentation scripts in `Segmentation/`. The separate script `bayes_rank.py` performs statistical analysis on the tracking results.

## Repository contents

```text
.
├── ImageProcess_v5.mlx
├── Filter_and_TrajOrg.mlx
├── MSD_v4.mlx
├── Confine_measure.mlx
├── bayes_rank.py
└── Segmentation/
    ├── Omni_batch.py
    └── Omni_batch.ipynb
```

| File | Purpose |
| --- | --- |
| `ImageProcess_v5.mlx` | Extract phase-contrast reference images, estimate movie drift, run segmentation, divide movies and masks into subregions, run `SMALL-LABS`, combine subregion tracks, and correct track coordinates for drift. |
| `Filter_and_TrajOrg.mlx` | Measure cell morphology, identify spatial locus clusters, filter cells and trajectories, and organize accepted trajectories for MSD analysis. |
| `MSD_v4.mlx` | Calculate individual and ensemble MSD curves using `msdanalyzer`, fit anomalous diffusion parameters, and save plots and per-track measurements. |
| `Confine_measure.mlx` | Calculate confinement area of each trajectory. |
| `Segmentation/Omni_batch.py` | Batch segmentation of phase-contrast reference images using the Omnipose `bact_phase_omni` model. Called from `ImageProcess_v5.mlx`. |
| `Segmentation/Omni_batch.ipynb` | Notebook version of the segmentation workflow, with optional image and mask visualization. |
| `bayes_rank.py` | Fit a Bayesian model with replicate effects nested within strains, summarize posterior strain ranks, and calculate pairwise comparison probabilities. |

## Requirements

### Tracking and MSD analysis

- MATLAB with Live Editor and Python integration (`pyenv` and `pyrunfile`).
- MATLAB functionality used by the scripts includes image processing, `lsqcurvefit`, `dbscan`, and `parfor`; install the corresponding Image Processing, Optimization, Statistics and Machine Learning, and Parallel Computing toolboxes.
- [SMALL-LABS](https://github.com/BiteenMatlab/SMALL-LABS), including `SMALLLABS_main`, on the MATLAB path.
- [`msdanalyzer`](https://tinevez.github.io/msdanalyzer/) for MSD calculation and fitting. Add the parent directory of `@msdanalyzer` to the MATLAB path.
- A Python environment compatible with the segmentation scripts, including `cellpose_omni`, `omnipose`, NumPy, and Matplotlib. See the [Omnipose project](https://github.com/kevinjohncutler/omnipose) for installation resources. The supplied scripts use the `cellpose_omni` API and the `bact_phase_omni` model; verify compatibility with the installed version.


## Data organization

The image-processing script expects a hierarchy of strain or condition, replicate, and movie folders:

```text
data/
└── AtWX278/
    ├── Rep1/
    │   ├── rep1_movie1/
    │   ├── rep1_movie2/
    │   └── ...
    └── Rep2/
        └── ...
```

For processing from raw data, each movie folder should contain a phase-contrast TIFF movie and a fluorescence ND2 movie. The image-processing script searches filenames using substrings including `Phase.tif`, `PhaseRef.tif`, `FP.nd2`, and `FP.mat`. Adjust these searches to match your filenames.

After processing, movie folders contain segmentation masks in `masks/`, subregion data in `Subregion/`, and drift and tracking results in `Results/`. The filtering script specifically loads:

```text
rep1_movie1/
├── masks/
│   └── <segmentation-mask file>
└── Results/
    └── Tracks_all.mat
```

## Bayesian statistical analysis

`bayes_rank.py` performs Bayesian hierarchical analysis of per-track mobility measurements, accounting for replicates nested within strains. It reports posterior strain rankings, pairwise comparison probabilities, and sampling diagnostics.

### Input

- CSV or Excel file with one column per strain.
- Replicate labels in the last column.
- Strictly positive measurements on their original scale; the script applies a `log10` transformation internally.
- Blank cells are allowed for missing measurements.

### Usage

```bash
python bayes_rank.py measurements.xlsx
```

### Model and interpretation

The script transforms measurements with `log10`, centers them by the overall mean log value, and models each observation using a strain mean, a replicate effect nested within strain, and a Student's t residual. Replicate labels are combined with strain names, so `Rep1` in two strains represents two separate replicate groups. 

For each posterior draw, strains are ranked by their modeled means, with rank 1 assigned to the largest mean. The text report contains:

- The most frequent posterior rank for each strain and the probability of that rank.
- Median ranks and 2.5th–97.5th percentile rank intervals.
- Pairwise probabilities `P(μ_row > μ_column)`.

Pairwise entries are posterior probabilities, rather than p-values. Check the reported sampling diagnostics when interpreting ranks and pairwise probabilities. 

## Associated data

The associated [Zenodo dataset](https://doi.org/10.5281/zenodo.23544930) contains processed tracking data for 16 *Agrobacterium* strains, with two replicates per strain and typically 6 movies per replicate (range: 5 to 8). It also includes an `Example_raw&processed/` folder with the complete raw data and processed results for one movie.

## Software references

- **SMALL-LABS:** Isaacoff, B. P., Li, Y., Lee, S. A., and Biteen, J. S. (2019). SMALL-LABS: Measuring Single-Molecule Intensity and Position in Obscuring Backgrounds. *Biophysical Journal*, 116(6), 975–982. [Article](https://doi.org/10.1016/j.bpj.2019.02.006); [software](https://github.com/BiteenMatlab/SMALL-LABS).
- **Omnipose:** [Project and installation resources](https://github.com/kevinjohncutler/omnipose).
- **msdanalyzer:** [Documentation and citation information](https://tinevez.github.io/msdanalyzer/).
