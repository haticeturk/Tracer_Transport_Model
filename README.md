# Tracer Transport Model
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)  
## Introduction
This repository contains Python code for a *Tracer Transport Model*.  The Model compute catchment-scale hydrologic transport using StorAge Selection (SAS) functions to estimate time variable transit time distributions (TTDs). The formulation is based a process-based hydrological and transport model, which builds on the DYNAMITE modeling framework. The model employs StorAge Selection (SAS) functions as function of age rank storage to describe how water of different ages is released from storage compartments into various outflows. Water fluxes and tracer dynamics are simulated simultaneously by integrating SAS functions within the hydrological model. In this conceptualization, the catchment is divided into five interconnected reservoirs: snow, canopy interception, the unsaturated root zone,  fast-response storage, and groundwater, which includes both active and passive components aimed at computing catchment-scale hydrologic transport. The model parameter and flux equations are described (@turk2024soil)  

Türk, H., Stumpp, C., Hrachowitz, M., Schulz, K., Strauss, P., Blöschl, G., & Stockinger, M. (2024). Soil moisture and precipitation intensity control the transit time distribution of quick flow in a flashy headwater catchment. Hydrology and Earth System Sciences Discussions, 2024, 1-33. https://doi.org/10.5194/hess-2024-359

## Repository Structure
- Data folder includes the Data_test.csv. # Example dataset for test runs.
- Model_Run.py is Script demonstrating model calibration & tracer tracking 
- Tracer_Mod_Wettness.py # Model config 1 - SAS paramameter depends on soil moisture  for preferential flow   
- Track_Tracer_Mod_Wetness.py # Variant config 1 allowing manual parameter changes for root zone & GW when running the model 


## Running the Model
### Example Test Run
An example of model calibration and tracer tracking is provided in **`Model_Run.py`**, which uses the synthetic dataset from **`Data/Data_test/`**. To run a test:

### Models 

1. **Tracer_Mod_Wettness.py**  
   - Uses uniform SAS functions for all compartments, except the **root zone’s preferential flow**.  
   - The root zone’s shape parameter (\(\alpha\)) is parameterized by soil moisture.

3. **Track_Tracer_Mod_Wetness.py** 
   - Allows users to modify shape parameters for both the **root zone** and **groundwater** compartments when calling the model.  
   - Ideal if you want to experiment with different SAS  parameters. 

```
@article{turk2024soil,
  title={Soil moisture and precipitation intensity control the transit time distribution of quick flow in a flashy headwater catchment},
  author={T{\"u}rk, Hatice and Stumpp, Christine and Hrachowitz, Markus and Schulz, Karsten and Strauss, Peter and Bl{\"o}schl, G{\"u}nter and Stockinger, Michael},
  journal={Hydrology and Earth System Sciences Discussions},
  volume={2024},
  pages={1--33},
  year={2024},
  publisher={G{\"o}ttingen, Germany}
}
``` 

