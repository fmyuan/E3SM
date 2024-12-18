[![E3SM Logo](https://e3sm.org/wp-content/themes/e3sm/assets/images/e3sm-logo.png)](https://e3sm.org)

<img src="https://gist.githubusercontent.com/ypwong22/2613576d9c70217e7a55ae9f5ef61445/raw/4d9dff783e95857b46ff66a0c99814e074b11078/ERWlogo.svg" width="500" alt="ELM-ERW Logo">

E3SM Land Model with Enhanced Rock Weathering
================================================================================

This is a branched version of the Energy Exascale Earth System Model (E3SM), a
state-of-the-art fully coupled model of the Earth's climate including important 
biogeochemical and cryospheric processes. E3SM is intended to address the most 
challenging and demanding climate-change research problems and Department of 
Energy mission needs while efficiently using DOE Leadership Computing Facilities.
The main branch of E3SM is here: [https://github.com/E3SM-Project/E3SM](https://github.com/E3SM-Project/E3SM)

The ELM-ERW branch is developed for modeling the large-scale carbon dioxide removal
potential of spreading silicate rock powder on soil. Specifically, it considers
the following processes:

- primary mineral dissolution
- CO<sub>2</sub> dissolution and cation exchange equilibria
- advection-diffusion and runoff loss of major ions
- secondary mineral formation
- soil pH dynamics and feedback to soil N<sub>2</sub>O and NO emissions
- phosphorus release during weathering and feedback to ecosystem productivity

The model has been tested in the land-only mode and runs at a speed commensurate
with the original E3SM Land Model. 

Funding Source
--------------------------------------------------------------------------------
This work is supported by the project "Coupled Ecosystem-and-engineering 
Decision-making Framework for Enhanced Weathering", under Oak Ridge National 
Laboratory (ORNL)'s Laboratory Directed Research and Development (LDRD) program. 

License
--------------------------------------------------------------------------------
The E3SM model is available under a BSD 3-clause license.
Please see [LICENSE](LICENSE) for details.
