# Self-contacting hydrogen bonds in protein structures

This repository contains a Jupyter-based structural bioinformatics workflow for identifying and analysing self-contacting hydrogen bonds involving cysteine (Cys), serine (Ser), and threonine (Thr) residues in high-resolution protein crystal structures.

The workflow covers dataset preparation, structure annotation, geometric hydrogen-bond detection, post-processing, visualization, and selected downstream quantum-chemical analyses. It also includes a historical implementation of the geometric criteria used by Eswar and Ramakrishnan (2000).

## Repository contents

| File | Purpose |
| --- | --- |
| `Self contacting Hbonds.ipynb` | Main end-to-end research notebook and analysis record |
| `H-bonding parameters.json` | Atom definitions and hydrogen-bond geometry parameters |
| `h_bonds_cys_ser_thr_as_donor.py` | Detects interactions in which Cys, Ser, or Thr side chains act as donors |
| `h_bonds_cys_ser_thr_met_as_acceptor.py` | Detects interactions in which Cys, Ser, Thr, or Met side chains act as acceptors |
| `h_bonds_cys_ser_thr_as_donor_script_used_for_Eswar_and_Ramakrishnan_2000.py` | Historical-criteria donor analysis |
| `ssec_strc_phi_psi_chi1_sasa_value_returner.py` | Retrieves secondary structure, torsion-angle, and solvent-accessibility annotations |
| `post_processing.py` | Filtering, classification, summary-table, and workbook generation utilities |
| `post_processing_script_used_for_Eswar_and_Ramakrishnan_2000.py` | Post-processing for the historical-criteria analysis |
| `range_to_table.py` | Converts Excel worksheet ranges to formatted tables |
| `chimerax_phi_psi_chi1_files_generator.py` | ChimeraX helper for generating phi, psi, and chi1 attribute files |
| `chimera_hadding_script_for_self_contacts.py` | UCSF Chimera helper for adding hydrogens to self-contact model systems |
| `chimera_hadding_script_for_chalcogen_bonds.py` | UCSF Chimera helper for chalcogen-bond model systems |
| `chimera_SELF_AND_NONSELF_CONTACTING_N_systems_creator_for_interaction_energy.py` | UCSF Chimera helper for preparing donor, acceptor, and complex model systems |
| `chimera_script_create_images_from_session_files.py` | UCSF Chimera helper for rendering saved sessions |

## Requirements

The notebook metadata records Python 3.12. The core Python dependencies are listed in `requirements.txt`.

Several stages depend on separately installed scientific software that is not distributed through PyPI:

- [UCSF Chimera](https://www.cgl.ucsf.edu/chimera/) for scripts importing `chimera`
- [UCSF ChimeraX](https://www.rbvi.ucsf.edu/chimerax/) for torsion-angle attribute generation
- MolProbity for hydrogen addition in the structure-preparation workflow
- STRIDE for secondary-structure assignment
- NACCESS for solvent-accessible surface area calculations
- PISCES for selecting non-redundant protein chains
- Gaussian for the optional quantum-chemical calculations

The classic Chimera scripts must be run inside Chimera's Python environment, and the ChimeraX script must be run inside ChimeraX. They are not expected to import successfully in an ordinary Python interpreter.

## Installation

Create a virtual environment and install the Python dependencies:

```bash
python -m venv .venv
```

On Windows PowerShell:

```powershell
.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
jupyter lab
```

On Linux or macOS:

```bash
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
jupyter lab
```

Open `Self contacting Hbonds.ipynb` from JupyterLab.

## Running the workflow

This repository preserves the original research notebook, including its cell outputs and Windows/WSL-oriented paths. It is an analysis record rather than a single-command pipeline.

Before running a section:

1. Read its Markdown heading and notes in the notebook.
2. Change its input and output paths to locations on your computer. Many cells contain absolute paths such as `D:\\...`, `G:\\...`, or `/mnt/d/...`.
3. Place `H-bonding parameters.json` in the location specified by that section, or change the corresponding variable to point to the repository copy.
4. Prepare the required PDB, MolProbity, STRIDE, NACCESS, ChimeraX, spreadsheet, or Gaussian files for that stage.
5. Run only the relevant cells in order. Some dataset-wide stages are documented in the notebook as taking several days.

The principal hydrogen-bond analysis modules are imported directly by the notebook, so keep the notebook and Python files in the same directory unless you also update the import paths.

## Data and reproducibility

Raw protein structures, generated annotations, spreadsheets, Chimera sessions, and Gaussian input/output files are not included. The `.gitignore` excludes common generated and potentially large research outputs so they are not committed accidentally.

For a reproducible analysis, record the following alongside each run:

- PDB identifiers and structure release/deposition dates
- Resolution and non-redundancy filters
- Versions of MolProbity, STRIDE, NACCESS, Chimera/ChimeraX, and Gaussian
- The exact hydrogen-bond parameter JSON used
- Any notebook path or parameter changes

## Notes

- The notebook contains both current and historical/alternative analysis routes; not every section is required for the core self-contacting hydrogen-bond search.
- Some automation cells use Windows-only packages (`pywin32`, PyGetWindow, and PyAutoGUI) and assume a particular desktop layout.
- Review generated structures and spreadsheets before using them in quantitative analysis.

## License

No license has been selected. Unless a license is added, the source remains under the copyright holder's default rights.
