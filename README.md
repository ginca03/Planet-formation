# Planet Formation Simulation

## Overview
This project investigates planet formation within protoplanetary disks using the **streaming instability** model. It explores how dust and gas interactions, turbulence, and particle growth contribute to planetesimal formation. The simulation implements numerical methods to ensure stability and minimize numerical diffusion, producing results consistent with astrophysical observations.

## Features
- **Simulates dust settling and turbulent diffusion** within a protoplanetary disk.
- **Uses numerical methods** to solve advection-diffusion equations.
- **Includes particle growth mechanisms** to study the formation of planetesimals.
- **Generates output files** that store the dust density evolution for further analysis.

## Installation and Dependencies
This simulation requires **CERN ROOT** and **Python** with Matplotlib installed.

### Install CERN ROOT:
```
sudo apt update
sudo apt install root-system
```

### Install Python dependencies:
```
pip install matplotlib
```

## Running the Simulation
1. Navigate to the project folder:
   ```sh
   cd Planet-formation
   ```
2. Execute the main simulation script using CERN ROOT:
   ```sh
   root main.cpp
   ```
3. After execution, two output files will be generated:
   - `R1.root` (dust density evolution at 1 AU)
   - `R100_diffusion.root` (dust density evolution with turbulence at 100 AU)

## Plotting the Results
To visualize the results, use the provided Python script:
```sh
python graph_dust.py R1.root
```
or
```sh
python graph_dust.py R100_diffusion.root
```

## Adjusting Log Scale
If you want to set the log scale on the x-axis of the "disk-height vs alpha" graph:
- Open the ROOT canvas
- Navigate to `View >> Editor`
- Check the box for **Log scale** in the x-axis section

## Computational Considerations
- **Execution Time:** The full simulation may take **~20 minutes**.
- **Numerical Stability:** The CFL condition is respected to ensure accuracy and minimize numerical diffusion.
- **Turbulence Effects:** The diffusion coefficient `α ≈ 10⁻⁴` was determined through observational comparisons.
- **Dust Growth:** Larger dust particles can overcome turbulence and settle, leading to planetesimal formation.

## Results and Findings
- **Dust settles faster at larger radial distances**, in agreement with theoretical expectations.
- **Turbulent diffusion prevents small particles from accumulating**, requiring sufficient particle growth for planet formation.
- **The simulation reproduces observed patterns of planetary formation in the Solar System**, confirming the validity of the model.

## Author
**Giancarlo Venturato**

## References
For a detailed explanation of the methodology and findings, refer to the [project report](https://github.com/ginca03/Planet-formation/blob/main/docs/06033934_Project3_Report.pdf).

