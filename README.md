# Object Conformal

R Code for Paper: *Conformal Inference for Random Objects*

This repository contains the R code used in the paper *Conformal Inference for Random Objects*. The structure of the repository is as follows:

- **`Real_Data/`**: Includes the code required to reproduce the results presented in the real data analysis section.
- **`Simulation/`**: Contains the code necessary to implement the proposed methods in the simulation studies.


---

## Simulation Studies:

1. **Proposed Method**:  
   The `Proposed` folder contains R code for implementing the proposed method. Different experimental settings are organized into subfolders. Inside each subfolder, run `Parallel.R` to generate the results.

2. **HPD Method**:  
   The `HPD` folder contains R code for running the HPD-split method. The original code is available at https://github.com/rizbicki/predictionBands.

3. **DCP Method**:  
   The `DCP` folder includes R code for running the DCP methods. The original code is available at https://github.com/kwuthrich/Replication_DCP.

4. **CQR Method**:  
   The `CQR` folder provides Python code for running the CQR methods. The original code is available at https://github.com/msesia/cqr-comparison.

---

## Real Data Analysis:

### **New York Taxi Data**
1. Download the data from [NYC TLC Trip Record Data](https://www.nyc.gov/site/tlc/about/tlc-trip-record-data.page).
2. Run `data_prepare.R` to preprocess the data.
3. Execute `conf.R` as the main script to generate the results presented in Section 7.1 of the paper.

### **U.S. Energy Data**
1. Download the data from [U.S. Energy Information Administration](https://www.eia.gov/electricity/data/state/).
2. Run `data_process.R` to preprocess the data.
3. Execute `split_conf.R` as the main script to produce the results in Section 7.2 of the paper.
4. View `Illustration.html` for an animation showing the dynamic changes in U.S. energy data.

---


