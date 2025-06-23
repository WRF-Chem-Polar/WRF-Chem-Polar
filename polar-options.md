# New polar options:

## In polar/main
Emission options:
- seas_opt (5, 6, 7, 8): This option controls the sea spray emissions source function.
  - 5: Sea-spray emissions from Monahan et al. (1986). No emissions below 100 nm diameters
  - 6 (recommended): Sea-spray emissions from Ioannidis et al. (2023). Uses Gong et al. (1997) above 200 nm and O'Dowd et al. (1997) below 200 nm; uses the whitecap fraction from Salisbury et al. (2013) as a function of windspeed. Uses the emission dependence on SST from Jaegle et al (2011).
  - 7: Sea-spray emissions from Gong et al. (1997) and O'Dowd et al. (1997). Uses Gong et al. (1997) above 200 nm and O'Dowd et al. (1997) below 200 nm; uses the whitecap fraction from Monahan et al. (1986) and does not use a correction factor for SST.
  - 8: Sea-spray emissions from Salter et al. (2015). Emissions are sensitive to SST, with fine emissions decreasing at higher SST and coarse emissions increasing at higher SST.
- seas_so4_opt (0, 1): This option controls the primary sulfate emission in sea-spray. Only compatible with seas_opt 5,6,7,8
  - 0 (default): No sea-spray sulfate emissions
  - 1: Sea-spray sulfate emissions assuming a sulfate/sodium mass ratio of 0.252 (Calhoun et al., 1991)
- seas_oa_opt (0, 1)
- seas_leads_opt (0, 1)
- dms_opt (0, 1, previoulsy dmsemis_opt)
- blowing_snow_opt (0, 1)
- biomass_burn_opt (6, 7)

New values for native WRF options:
- biomass_burn_opt (6, 7)
- seas_opt (5, 6, 7, 8)

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
