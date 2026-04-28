---
editor_options: 
  markdown: 
    wrap: 72
---

# Replication Materials: Stochastic Modelling of Environmental Risk in Megacities

This repository contains the `R` code required to reproduce the algorithms, Monte Carlo simulations, empirical results, tables, and figures presented in the manuscript:

> **Title:** Stochastic Modelling of Environmental Risk in Megacities: Accuracy and Scalability of Multiscale Geographically and Temporally Weighted Regression

## Repository Structure

The project follows a reproducible, project-oriented workflow. The files are organised as follows:

* **`data/`**: This folder is a placeholder for the analytical datasets. Due to file size constraints, the primary data is hosted on Figshare.
    * **Dataset DOI:** [https://doi.org/10.6084/m9.figshare.32114857](https://doi.org/10.6084/m9.figshare.32114857)
    * **Required File:** `sao_paulo_mgtwr_panel_dual_geom.rds`

* **`script/`**: Contains the R scripts required for the analysis.
    * `Functions_code.R`: Contains all underlying algorithms (2SCALL, TDS-MGTWR, Back-fitting), spatial data utilities, and plotting functions.
    * `Run_functions_code.R`: The primary execution script that runs the Monte Carlo simulation and the empirical application in São Paulo.

* **`outputs/`**: Directory where generated plots, model evaluation metrics, tables, and cached model results (`.rds` files) are saved. Subdirectories include `/simulated`, `/monte_carlo`, and `/real`.

* **`Vanessa.Rproj`**: RStudio Project file. **This is the entry point for the analysis** to ensure relative paths resolve correctly.

## Software and Dependencies

The analysis was performed using `R`. To ensure reproducibility across different environments, this project uses the `here` package for relative file paths, eliminating the need for absolute paths (e.g., `setwd()`).

**Required Packages:**
The script utilises `pacman` to handle package installation and loading. Key dependencies include:

* **Spatial & Census Data:** `sf`, `spdep`, `geobr`, `censobr`, `classInt`
* **Modelling & Math:** `MASS`, `RANN`
* **Data Manipulation:** `dplyr`, `tidyr`, `purrr`, `stringr`, `readr`
* **Visualization & Tables:** `ggplot2`, `ggspatial`, `ggridges`, `patchwork`, `ggdist`, `ggrepel`, `corrplot`, `gtsummary`, `gt`, `knitr`, `kableExtra`

## Instructions for Reproduction

To reproduce the analysis, please follow these steps strictly to ensure file paths and dependencies resolve correctly:

1.  **Download the Repository:** Download this GitHub repository (Click "Code" -> "Download ZIP") and unzip it to your local machine.
2.  **Download the Data:** Visit the [Figshare link](https://doi.org/10.6084/m9.figshare.32114857) and download the processed dataset. Place the `.rds` file inside the `data/` folder of this project.
3.  **Open the Project:** Open the file `Vanessa.Rproj` in RStudio.
    * *Critical:* Opening the project file sets the working directory to the project root automatically. **Do not use `setwd()`.**
4.  **Execute the Analysis:** Open the script `script/Run_functions_code.R` and run the code.
    * *Note on Computation Time:* The Monte Carlo simulation section (`N_sim = 300`) is computationally intensive and may take several hours to complete. For a rapid test, consider lowering the `N_sim` value or skipping directly to the São Paulo empirical application.

## Data Source

The datasets constructed for this empirical analysis integrate information from the following primary sources:

* **Census Data:** Extracted via the `censobr` and `geobr` packages (IBGE Census data for 2000, 2010, and 2022).
* **Urban Infrastructure & Informality:** Geospatial data representing transport hubs (train/bus) and informal settlements (slums, tenements, and irregular lots) within São Paulo (GeoSampa).

## Contact

For questions regarding the code or data, please contact the corresponding author via the journal submission system.
