# New namelist.input options in the polar version

This file describes the new namelist.input options in the WRF-Chem-Polar version.

*Many of these options are not specific to polar regions and can be used for modeling outside the poles.  The main polar specific developments are related to sea ice, which shoudl not negatively impact the results outside the polar regions.*

## Chemical mechanism options (&chem namelist)

### chem_opt (0, 112, 114, 202, 212, 300)
This option controls the mechanism for gas-phase chemistry and aerosols.
  - 0 (default): Disables gas-phase chemistry and aerosols. Used to run WRF (met-only) simulations, but it's better to disable the whole &chem namelist rather than each individual option.
  - 112: MOZCART_KPP option, MOZART chemistry with simple bulk GOCART aerosols. Useful for quick tests or studies focusing on gas-phase chemistry.
  - 114 (not tested): T1_MOZCART_KPP option, MOZART chemistry with updated T1 mechanism, with simple bulk GOCART aerosols. Should predict gas-phase species more accurately than MOZCART_KPP, but we never tested it.
  - 202: MOZART_MOSAIC_4BIN_AQ_KPP: MOZART chemistry + MOSAIC aerosols with 4 size bins, including secondary organic aerosols, aqueous chemistry and aerosol-cloud interactions. Includes simplified DMS chemistry. Useful for simulations in continental regions or for studies not focused on sulfate aerosols or aerosol-cloud interactions.
  - 212 (recommended): MOZART_MOSAIC_4BIN_DMS_AQ_KPP: MOZART chemistry + MOSAIC aerosols with 4 size bins, including secondary organic aerosols, aqueous chemistry and aerosol-cloud interactions, and complex DMS and MSA chemistry. Recommended option for chemistry+aerosol runs in the polar regions.
  - 300: GOCART_SIMPLE: GOCART simple bulk aerosols, no gas-phase chemistry. Recommended for very fast tests of sea-spray or dust emissions and transport.

## Emission options (&chem namelist)

### dust_opt (0, 1, 3, 4, 5, 6, 7, 13)
This option controls mineral dust aerosol emissions.
  - 0 (default): Disables dust emissions
  - 1 (GOCART-only, not recommended): WRF-GOCART dust emissions following Ginoux (2004) and the GOCART model, calculated in gocart_dust_driver
  - 3 (GOCART-only, untested): AFWA dust emissions following Legrand et al. (2019), calculated in gocart_dust_afwa_driver
  - 4 (untested): UoC dust emissions following Shao (2001), calculated in uoc_dust_driver
  - 5 (GOCART and MOSAIC, not recommended): WRF-GOCART dust emissions following Ginoux (2004) and the GOCART model, calculated in mineraldust_emis. Should be identical to options 1 and 13.
  - 6 (recommended, GOCART and MOSAIC): WRF-Chem-Polar (WCP) dust emissions based on Ginoux (2004), with updates from Legrand et al. (2019) and additional corrections to allow for emissions from high-latitude regions.
  - 7 (GOCART and MOSAIC, not recommended): FLEXDUST dust emissions following Groot Zwaaftink et al. (2016) and the FLEXDUST model, not including specific treatment for emissions of Icelandic dust.
  - 13 (MOSAIC-only, not recommended): WRF-GOCART dust emissions following Ginoux (2004) and the GOCART model, calculated in mosaic_dust_gocartemis in mosaic_addemiss.
### dust_erod_opt (0, 1)
This option controls the source of the soil erodibility fraction used for dust emissions.
  - 0 (default): Uses the default EROD from WPS, read from wrfinput
  - 1 (recommended for high-latitude dust): Uses EROD_HL read from wrfinput. EROD_HL needs to be included in the wrfinput file with an additional preprocessor.
### seas_opt (5, 6, 7, 8)
This option controls the sea spray emissions source function.
  - 0 (default): Disables aerosol emissions from sea-spray.
  - 5: Sea-spray emissions from Monahan et al. (1986). No emissions below 100 nm diameters
  - 6 (recommended): Sea-spray emissions from Ioannidis et al. (2023). Uses Gong et al. (1997) above 200 nm and O'Dowd et al. (1997) below 200 nm; uses the whitecap fraction from Salisbury et al. (2013) as a function of windspeed. Uses the emission dependence on SST from Jaegle et al (2011).
  - 7: Sea-spray emissions from Gong et al. (1997) and O'Dowd et al. (1997). Uses Gong et al. (1997) above 200 nm and O'Dowd et al. (1997) below 200 nm; uses the whitecap fraction from Monahan et al. (1986) and does not use a correction factor for SST.
  - 8: Sea-spray emissions from Salter et al. (2015). Emissions are sensitive to SST, with fine emissions decreasing at higher SST and coarse emissions increasing at higher SST.
### seas_so4_opt (0, 1)
This option controls the primary sulfate emission in sea-spray. Only compatible with seas_opt 5,6,7,8
  - 0 (default): No sea-spray sulfate emissions
  - 1: Sea-spray sulfate emissions assuming a sulfate/sodium mass ratio of 0.252 (Calhoun et al., 1991)
### seas_oa_opt (0, 1)
This option controls the primary marine organic emissions in sea-spray. Only compatible with seas_opt 5,6,7,8
  - 0 (default): No marine organic emissions from sea-spray
  - 1 (recommended): Marine organic emissions from Vignati et al. (2010). Requires chlorophyll-a input in the wrflowinp file in variable CHLOROA
### seas_leads_opt (0, 1, 2)
This option controls the sea-spray emissions from sea ice leads. Only compatible with seas_opt 5,6,7,8
  - 0 (default): Sea-spray emissions from leads use the seas_opt source function weighted by the open ocean fraction (1-seaice_concentration)
  - 1: Sea-spray emissions from leads use the seas_opt source function weighted by the leads fraction, also applying a correction factor to reduce the emissions to account for the reduced wind fetch over leads (Lapere et al., 2024). Requires lead fraction input in the wrflowinp file in variable LEADFRAC
  - 2: Sea-spray emissions from leads use the seas_opt source function weighted by the leads fraction, also applying a correction factor to reduce the emissions to account for the reduced wind fetch over leads (Lapere et al., 2024). Lead fraction is calculated online as 1-seaice where seaice is above 80%, ignoring LEADFRAC even if it is available.
### dms_opt (0, 1, 2, 3; previoulsy dmsemis_opt)
This option controls the treatment of dimethyl sulfide (DMS) emissions from the surface ocean in the model. This replaces the old option dmsemis_opt (deprecated but can still be used if needed). Options 1, 2, and 3 require oceanic DMS concentration input in the wrflowinp file in variable DMS_OCEAN. DMS_OCEAN can be taken from the climatologies of [Lana 2011](https://doi.org/10.1029/2010GB003850), [Hulswar 2021](https://doi.org/10.5194/essd-14-2963-2022) or the CSIB model [(Hayashida et al., 2019)](https://doi.org/10.5194/gmd-12-1965-2019). For options 1, 2, and 3, emissions from sea ice regions are scaled by the open ocean fraction to the power of 0.4 (Loose et al., 2009).
  - 0 (default): No DMS emissions from the ocean surface.
  - 1 (recommended): DMS emissions using the Nightingale et al. (2000) sea-air flux parameterization.
  - 2: DMS emissions using the Liss and Merlivat (1986) sea-air flux parameterization.
  - 3: DMS emissions using the Wanninkhof (2014) sea-air flux parameterization.
### nuc_msa_opt (0, 1)
This option controls the activation of new particle formation from gas MSA
  - 0 (default): No nucleation of MSA aerosols
  - 1: Activation of NPF from MSA following [Riccobono et al. (2014)](https://doi.org/10.1126/science.1243527)
### nuc_msa_fac (real)
Tuning factor applied to the [Riccobono et al. (2014)](https://doi.org/10.1126/science.1243527) formula when nuc_msa_opt is activated. The default value is 1.0, i.e. the original formulation of the nucleation rate is used. For Arctic studies, we have found that a value of 0.067 is more appropriate for reproducing observations of aerosol number concentration.
### nuc_sulf_opt (0, 1, 2, 3, 4)
This option controls which scheme is used for the nucleation of sulfate aerosols from sulfuric acid.
  - 0: no nucleation
  - 1: [Napari et al (2002)](https://doi.org/10.1029/2002JD002132) (invalid, deactivates nucleation)
  - 2: [Wexler et al (1994)](https://doi.org/10.1016/1352-2310(94)90129-5) (default)
  - 3: [Kulmala et al (1998)](https://doi.org/10.1029/97JD03718)
  - 4: [Vehkamaki et al. (2002)](https://doi.org/10.1029/2002JD002184)
### blowing_snow_opt (0, 1): This option controls the emissions from blowing snow.
  - 0 (default): Disables emissions associated with blowing snow
  - 1: Includes sea salt aerosol emissions from blowing snow (on the main, halogen and mercury branches) and bromine emissions from blowing snow (halogens and mercury branches)
### biomass_burn_opt (6, 7): This option controls the biomass burning emissions.
  - 0 (default): Disables biomass burning emissions.
  - 6: Biomass burning emissions for SAPRC gas-phase chemistry mechanisms. Requires wrffirechemi input files for SAPRC, created by the fire_emis preprocessor.
  - 7: Biomass burning emissions for CBMZ gas-phase chemistry mechanisms. Requires wrffirechemi input files for CBMZ, created by the fire_emis preprocessor.
### surface_snow_opt (0, 1): This option controls halogen emissions and recycling from surface snow in the halogens and mercury branches.
  - 0 (default): Disables halogen emissions and recycling from surface snow.
  - 1: Halogen emissions and recycling from surface snow, following Toyota et al. (2011).

### pbap_opt (0, 1): This option controls Primary Biological Aerosol Particle (PBAP) emissions from land surfaces
  - 0 (default): Disables PBAP emissions
  - 1: Emissions of fungal spores and bacteria from land surfaces following Hummel et al. (2015) and Burrows et al. (2009) respectively. The option is only active for GOCART aerosols, and PBAP are emitted as OC1 (hydrophobic organic carbon).

## Aerosol-cloud interaction options
### aci_wrfchem_opt (0, 1, 2): This option controls the aerosol-cloud interactions for liquid droplets in WRF-Chem.
This option only works with Thompson aerosol-aware microphysics, but not with other microphysics scheme.
  - 0 (default): Uses the default aerosol-cloud interaction in Thompson aerosol-aware microphysics. Uses the default water-friendly and ice-friendly aerosol number concentrations supplied by the Thompson aerosol-aware microphysics scheme.
  - 1: In Thompson aerosol-aware microphysics, calculate aerosol-cloud interactions for liquid clouds with the Abdul-Razzak and Ghan 2002 (ARG02) parameterization and for ice clouds with the DeMott2010  (D10) parameterization. ARG02 is used, both for the droplet number in the microphysics code and to set cw aerosols for cloud chemistry and wet scavenging. The aerosol parameters are calculated from WRF-Chem aerosols (MOSAIC-4bin only). Ice friendly aerosols for D10 are also calculated from the WRF-Chem dust aerosol.
  - 2: In Thompson aerosol-aware microphysics, calculate aerosol-cloud interactions for liquid clouds with the Thompson and Eidhammer 2014 (TE14) parameterization and for ice clouds with the DeMott2010 (D10) parameterization. Water friendly and ice friendly aerosols for TE14 and D10 are set up from WRF-Chem aerosols (GOCART or MOSAIC-4bin only). For MOSAIC-4bin, cw aerosols for cloud chemistry and wet scavenging are still calculated using Abdul-Razzak and Ghan 2002.

### aci_wrf_opt (0, 1, 2, 3), &phys namelist option: This option controls the aerosol-cloud interactions for liquid droplets in WRF (without chemistry).
This option only works with Thompson aerosol-aware microphysics, but not with other microphysics scheme.
  - 0 (default): Uses the default aerosol-cloud interaction in Thompson aerosol-aware microphysics. Uses the default water-friendly and ice-friendly aerosol number concentrations supplied by the Thompson aerosol-aware microphysics scheme.

For option 1,2,3, an external 4D aerosol climatology is read in Thompson aerosol-aware microphysics in the auxiliary input file auxinput18. The auxinput file needs to contain QNWFA_EXT and QNIFA_EXT fields in #/kg (option 1,2,3), WRF_AER_SO4_EXT in µg/kg (option 2) and WRF_AER_SOLUBLE_EXT in µg/kg (option 3). These namelist options also need to be included in namelist &time_control for a 3-hourly climatology:

auxinput18_inname                    = 'wrf_wfa_ifa_d\<domain\>_\<date\>'

auxinput18_interval_h                = 3,

frames_per_auxinput18                = 1,

io_form_auxinput18                   = 2

  - 1: Calculates aerosol-cloud interactions in Thompson aerosol-aware microphysics with the Thompson and Eidhammer (2014) activation scheme, using aerosols from an external NWFA and NIFA aerosol climatology. The climatology is read from auxiliary input file auxinput18. The aux file needs to contain QNWFA_EXT and QNIFA_EXT fields, the water friendly and ice friendly aerosol numbers in #/kg.
  - 2: Calculates aerosol-cloud interactions in Thompson aerosol-aware microphysics with the Boucher and Lohmann (1995) scheme, using a sulfate aerosol concentration climatology. The option requires additional 3D input fields QNWFA_EXT, QNIFA_EXT and WRF_AER_SO4_EXT (accumulation-mode sulfate mixing ratio in ug/kg). In the radiation driver, cloud droplet number concentrations predicted by the microphysics are overwritten by values calculated from Boucher and Lohmann (1995).
  - 3: Calculates aerosol-cloud interactions in Thompson aerosol-aware microphysics with the LMDZ6 (Madeleine et al., 2020) scheme, using a soluble aerosol concentration climatology. The option requires additional 3D input fields QNWFA_EXT, QNIFA_EXT and WRF_AER_SOLUBLE_EXT (accumulation-mode sulfate+seasalt+ammonium aerosol mixing ratio in ug/kg). In the radiation driver, cloud droplet number concentrations predicted by the microphysics are overwritten by values calculated from Madeleine et al. (2020).

### mp_morr_icenuc_option (0, 1): This option controls the aerosol-cloud interactions for ice crystals in Morrison microphysics in WRF-Chem.
- 0 (default): Use the default Morrison ice nucleation scheme depending on temperature only
- 1 (in development, not recommended, for testing only): Use a simplified aerosol-aware ice nucleation by increasing the heterogeneous freezing temperature in Morrison by one degree for each order of magnitude increase in dust aerosol number.
- 2 (in development, not recommended): Use classical nucleation theory (CNT) for deposition from Keita et al. (2020) and immersion from Hoose et al. (2010) to calculate heterogeneous freezing in Morrison microphysics. CNT aerosol parameters are calculated from WRF-Chem dust aerosols (MOSAIC-4bin only).

## Deposition options
###  mosaic_aer_settling_opt (0, 1)
This option controls aerosol sedimentation above the first vertical level for MOSAIC aerosols.
   - 0 (default): No aerosol sedimentation above the first level in MOSAIC. Settling velocities are calculated in the first level and taken into account in the dry deposition velocity.
   - 1: Includes aerosol sedimentation of MOSAIC aerosols at all vertical levels.
###  wetscav_onoff (0,1)
This option controls wet scavenging of aerosols by precipitation.
  - 0 (default): No wet scavenging, unchanged from upstream WRF-Chem
  - 1 (recommended): For GOCART mechanisms, includes wet scavenging following Luo et al. (2019, 2020), whereas upstream WRF had no wet scavenging for GOCART aerosols. For other aerosol mechanisms, the behavior of this option is unchanged from upstream WRF.
### aer_drydep_opt (0,1,301,311)
This option controls aerosol dry deposition
  - 0 (default): No dry deposition, unchanged from upstream WRF-Chem
  - 1 (MOSAIC aerosols, recommended): Dry deposition is enabled for aerosols. Calculates dry deposition velocities with Binkowski and Shankar (1995). Unchanged from upstream WRF-Chem.
  - 1 (GOCART aerosols, recommended): Dry deposition is enabled for aerosols. Calculates dry deposition velocities with Emmerson et al. (2020).
  - 301 (MOSAIC only): Dry deposition is enabled for aerosols. Calculates dry deposition velocities with Zhang (2001). Unchanged from upstream WRF except we added a mapping to allow using this option with MODIS land use.
  - 311 (MOSAIC only): Dry deposition is enabled for aerosols. Calculates dry deposition velocities with Zhang (2001), with updated parameters. Unchanged from upstream WRF except we added a mapping to allow using this option with MODIS land use.
