# GHG Data Useful Scripts

A collection of R scripts and functions for analyzing greenhouse gas flux data, with a focus on ecosystem respiration (Reco) and gross primary productivity (GPP) modeling.

## Getting Started

### Download the Repository

**Option 1: Download as ZIP**
1. Click the green "Code" button at the top of this page
2. Select "Download ZIP"
3. Extract the ZIP file to your desired location

**Option 2: Clone with Git**
```bash
git clone https://github.com/bravemaster3/GHG-data-useful-scripts.git
```

### Usage

1. Open the R project file `GHG-data-useful-scripts.Rproj` in RStudio
2. Open and run the example script: `Reco_Temp_GPP_PAR_fitting_examples_fluxnet.R`
   - A sample dataset is provided in the `data/` folder
   - To use your own data: place your file in the `data/` folder and update the file path in the example script
   - install required packages if any is missing
3. Modify the script for your own data as needed

### Exploring Functions

- **View function code**: Press `F2` while cursor is on a function name in the script
- **Browse all functions**: Open scripts in the `functions/` folder

## Current Features

- Lloyd & Taylor (1994) temperature response fitting for Reco
- Hyperbolic light response (Michaelis-Menten) fitting for GPP


---

More scripts and functions will be added over time.
