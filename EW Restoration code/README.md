# Joint fit of common and decorrelating directions

This folder contains the two figures in `output/pdf/joint_fit_paper_note.pdf`, their plotting code, the saved fit results, and the numerical inputs used to compute those results.

- **Figure 1:** S-only, D-only, and joint S,D constraints at an upper invariant-mass cut of 3 TeV.
- **Figure 2:** Joint S,D contours for upper invariant-mass cuts of 0.5, 1, 2, and 3 TeV.

## Reproduce the plots

Use Python 3.12 or later. Run these commands from this folder:

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python reproduce_plots.py
```

On Windows, activate the environment with `.venv\Scripts\activate` instead.

The script writes PDF, SVG, and PNG versions of both figures to `output/pdf/`. It also checks the contour equations and writes `plot_validation.json`. Plotting uses only NumPy and Matplotlib; ReportLab and pypdf are needed only to rebuild the note. The pinned versions are the versions used to verify this package.

To rebuild the one-page note:

```sh
python output/joint-sd-fit/make_paper_style_note.py
```

## Recompute the fits

The saved integration results are included, so Mathematica is unnecessary for this step:

```sh
python output/joint-sd-fit/fit_joint.py
python reproduce_plots.py
```

This regenerates `fit_results.json`, `results.txt`, `sensitivity.csv`, and `plot_matrices.wl`. The fitter checks quadrature convergence, equivalence of the channel-space and S,D fits, and the profiled confidence intervals.

## Recompute the integrations (optional)

This step requires a working WolframScript/Mathematica installation. Run from the repository root:

```sh
SCAN_ORDER=96 wolframscript -file output/joint-sd-fit/integrate_joint.wls
SCAN_ORDER=128 wolframscript -file output/joint-sd-fit/integrate_joint.wls
python output/joint-sd-fit/fit_joint.py
python reproduce_plots.py
```

The exported amplitude expressions and tabulated luminosities needed by this integrator are included. The integration starts from those saved inputs; regenerating the amplitudes from the original research notebooks or regenerating luminosities from external PDF grids is outside this package. Neither the original workspace nor a download of those grids is needed to run the included integration script.

## Files and provenance

| Path | Contents |
| --- | --- |
| `reproduce_plots.py` | Exact plotting code extracted from the paper-note builder, with portable paths and standalone PDF/SVG/PNG output |
| `output/joint-sd-fit/make_paper_style_note.py` | One-page note builder, adjusted to use Matplotlib's installed fonts |
| `output/joint-sd-fit/fit_joint.py` | Fit calculation from saved response integrals |
| `output/joint-sd-fit/fit_results.json` | Information matrices and fit results used by both plots |
| `output/joint-sd-fit/response_integrals_order*.json` | Order-96 and order-128 integration results used by the fitter |
| `output/joint-sd-fit/response_integrals_order*.wl` | The same saved integrations in Wolfram Language format |
| `output/joint-sd-fit/integrate_joint.wls` | Integration code for the joint linear response |
| `output/n7-masscut/*_amplitudes.wl` | Exported HEFT channel amplitudes |
| `output/chq3-masscut/*_amplitudes.wl` | Exported SMEFT channel amplitudes |
| `output/n7-masscut/luminosity_1.dat` | CTEQ6L1 luminosity table |
| `output/n7-masscut/luminosity_2.dat` | nCTEQ15 proton luminosity table |
| `output/pdf/` | Both figures, the short note, and the contour validation result |

All scientific input files were copied from the calculation used for the note. Folder names are retained so the integration script resolves its inputs. No machine-specific Python installation or absolute workspace path is required.

## Plot conventions

The horizontal axis is C_Hq^(3) in TeV^-2; the vertical axis is n7/10^-4. Internal fit coordinates are `(n7/10^-4, C_Hq^(3)/0.01)`. Both coefficients float. All displayed contours use Delta chi^2 = 2.30, the nominal 68.3% joint region for two parameters. The black point is the Standard Model; the dotted horizontal line is the selected n7 = 0 SMEFT direction.

The benchmark uses SM Asimov rates, sqrt(s) = 14 TeV, nCTEQ15, interference-order templates, kappa7 = 1, and pT bins of 100-200 and 200-1500 GeV. Fractional channel errors are WZ = (0.088, 0.091) and Wh = (0.15, 0.15). Channel and bin measurements are uncorrelated, while the induced S,D covariance is retained. These fractional errors are held fixed across mass cuts. Source W+ flavor/charge conventions are preserved. The figures describe expected sensitivity in this restricted two-parameter benchmark.

## Verification

The packaged plot script, fitter, and note builder were run successfully. All regenerated fit entries match the original saved results within a relative tolerance of 1e-12. The maximum plotted contour residual is 1.78e-15. Both figures and the rebuilt note were visually checked. The optional Wolfram integration was not rerun during packaging; its saved order-96 and order-128 outputs are included.
