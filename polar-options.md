# New options in the polar version:

This file describes the new options or new option values in the WRF-Chem-Polar version. Unless specified otherwise, these options are not specific to polar regions and should improve model results in other regions.

## In polar/main

Emission options (&chem namelist):
- seas_opt (5, 6, 7, 8): This option controls the sea spray emissions source function.
  - 0: Disables aerosol emissions from sea-spray.
  - 5: Sea-spray emissions from Monahan et al. (1986). No emissions below 100 nm diameters
  - 6 (recommended): Sea-spray emissions from Ioannidis et al. (2023). Uses Gong et al. (1997) above 200 nm and O'Dowd et al. (1997) below 200 nm; uses the whitecap fraction from Salisbury et al. (2013) as a function of windspeed. Uses the emission dependence on SST from Jaegle et al (2011).
  - 7: Sea-spray emissions from Gong et al. (1997) and O'Dowd et al. (1997). Uses Gong et al. (1997) above 200 nm and O'Dowd et al. (1997) below 200 nm; uses the whitecap fraction from Monahan et al. (1986) and does not use a correction factor for SST.
  - 8: Sea-spray emissions from Salter et al. (2015). Emissions are sensitive to SST, with fine emissions decreasing at higher SST and coarse emissions increasing at higher SST.
- seas_so4_opt (0, 1): This option controls the primary sulfate emission in sea-spray. Only compatible with seas_opt 5,6,7,8
  - 0 (default): No sea-spray sulfate emissions
  - 1: Sea-spray sulfate emissions assuming a sulfate/sodium mass ratio of 0.252 (Calhoun et al., 1991)
- seas_oa_opt (0, 1): This option controls the primary marine organic emissions in sea-spray. Only compatible with seas_opt 5,6,7,8
  - 0 (default): No marine organic emissions from sea-spray
  - 1 (recommended): Marine organic emissions from Vignati et al. (2010). Requires chlorophyll-a input in the wrflowinp files in variable CHLOROA
- seas_leads_opt (0, 1): This option controls the sea-spray emissions from sea ice leads. Only compatible with seas_opt 5,6,7,8
  - 0 (default): Sea-spray emissions from leads use the seas_opt source function weighted by the open ocean fraction (1-seaice_concentration)
  - 1: Sea-spray emissions from leads use the seas_opt source function weighted by the leads fraction, also applying a correction factor to reduce the emissions to account for the reduced wind fetch over leads (Lapere et al., 2024). Requires lead fraction input in the wrflowinp files in variable LEADFRAC
- dms_opt (0, 1, previoulsy dmsemis_opt): This option controls the treatment of dimethyl sulfide (DMS) emissions from the surface ocean in the model. This replaces the old option dmsemis_opt (deprecated but can still be used if needed).
  - 0 (default): No DMS emissions from the ocean surface.
  - 1 (recommended): DMS emissions from the surface ocean are activated for the open ocean grid cell fraction, using the Nightingale et al. (2000) air–sea flux parameterization. In sea ice, emissions are scaled by the open ocean fraction to the power of 0.4 (Loose et al., 2009). Requires oceanic DMS concentration input in the wrflowinp files in variable DMS_OCEAN. DMS_OCEAN can be taken from the climatologies of [Lana 2011](https://doi.org/10.1029/2010GB003850), [Hulswar 2021](https://doi.org/10.5194/essd-14-2963-2022) or the CSIB model [(Hayashida et al., 2019)](https://doi.org/10.5194/gmd-12-1965-2019).
- blowing_snow_opt (0, 1): This option controls the emissions from blowing snow.
  - 0 (default): Disables emissions associated with blowing snow
  - 1: Includes sea salt aerosol emissions from blowing snow (on the main, halogen and mercury branches) and bromine emissions from blowing snow (halogens and mercury branches)
- biomass_burn_opt (6, 7): This option controls the biomass burning emissions.
  - 0 (default): Disables biomass burning emissions.
  - 6: Biomass burning emissions for SAPRC gas-phase chemistry mechanisms. Requires wrffirechemi input files for SAPRC, created by the fire_emis preprocessor.
  - 7: Biomass burning emissions for CBMZ gas-phase chemistry mechanisms. Requires wrffirechemi input files for CBMZ, created by the fire_emis preprocessor.

Aerosol-cloud interaction options:
- aci_wrfchem_opt (0, 1, 2)
- aci_wrf_opt (0, 1, 2, 3)
- mp_morr_icenuc_option (0, 1)

Deposition options:
- mosaic_aer_settling_opt (0, 1)

## In the halogen and mercury branches only
- surface_snow_opt (0, 1)

## New developments that might introduce new polar options:
- The 17-reaction DMS mechanism (von Glasow and Crutzen 2004, Marelle et al 2025).
- Nucleation from MSA
