# New polar options:

## In polar/main
Emission options:
- dms_opt (0, 1, previoulsy dmsemis_opt)
- blowing_snow_opt (0, 1)
- seas_so4_opt (0, 1)
- seas_oa_opt (0, 1)
- seas_leads_opt (0, 1)
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
