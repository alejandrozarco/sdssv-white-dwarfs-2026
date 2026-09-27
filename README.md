# White dwarfs 2026: measurements

Measurements of white dwarfs from public SDSS-V DR20 spectra (Astra 0.8.1), with DESI DR1, SDSS/BOSS, ESO X-shooter, TESS, HST/COS, GALEX, ATLAS, ZTF and Gaia DR3 epoch photometry. Each table gives the measured quantities with existing SIMBAD, MWDD (snapshot 2026-08-05) and SDSS-V SnowWhite classifications. The scripts in `scripts/` download the public data and recompute every table and figure. Data were retrieved 2026-09-23 to 2026-09-27. Methods are in [METHODS.md](METHODS.md).

## Topics

| topic | page | tables | objects |
|---|---|---|---|
| Ca II triplet emission (gaseous discs) | [docs/gas_discs.md](docs/gas_discs.md) | `gas_disc_white_dwarfs.csv`, `gas_disc_epochs_*.csv`, `gas_disc_screen.csv`, `desi_gas_disc_screen.csv` | 6 |
| Carbon lines | [docs/carbon.md](docs/carbon.md) | `carbon_white_dwarfs.csv`, `carbon_screen.csv`, `carbon_*_features.csv` | 9 |
| Zeeman splitting | [docs/zeeman.md](docs/zeeman.md) | `magnetic_zeeman.csv` | 30 |
| Photometric periods | [docs/periodic.md](docs/periodic.md) | `periodic_6021870154194477312.csv`, `periodic_white_dwarfs.csv`, `reflection_3107374277060584064*.csv`, `tess_ffi_*.csv` | 18 |
| TESS amplitude spectra | [docs/zz_ceti.md](docs/zz_ceti.md) | `zz_ceti_objects.csv`, `zz_ceti_tess_sectors.csv` | 5 |
| Eclipse and Balmer emission | [docs/eclipse_and_emission.md](docs/eclipse_and_emission.md) | `eclipsing_4731701084150029824.csv`, `balmer_emission.csv` | 3 |
| Hot white dwarfs with He II lines | [docs/hot_white_dwarfs.md](docs/hot_white_dwarfs.md) | `hot_white_dwarfs.csv` | 6 |
| Short-period white dwarfs with irradiated companions | [docs/irradiated_companions.md](docs/irradiated_companions.md) | `irradiated_companions.csv`, `irradiated_companions_periods.csv`, `desi_halpha_6914922055508553984.csv` | 6 |
| Day-scale periods of hot white dwarfs | [docs/hot_wd_periods.md](docs/hot_wd_periods.md) | `hot_wd_periods.csv` | 9 |
| Periods of DA white dwarfs with emission lines | [docs/dae_wd_periods.md](docs/dae_wd_periods.md) | `dae_wd_periods.csv` | 7 |

## Objects

One page per object with the measurement and its supporting data (tables, input photometry and spectra identifiers, figures, scripts).

<!-- object-index:start -->
### Ca II triplet emission (gaseous discs) (6)

| object | description | measurement |
|---|---|---|
| [WD 0856+048](docs/objects/578709631539357440.md) ★★★★ | White dwarf, G = 18.341; catalogued: SnowWhite DA; SIMBAD WD*/DA; MWDD DA; DESI DR1 DA. | Ca II emission, EW 17.55 A |
| [WD J1959+2208](docs/objects/1827014701883095680.md) ★★★ | White dwarf, G = 16.78; catalogued: SnowWhite DBA/DB; SIMBAD WD*/DB; MWDD DB. | Ca II emission, EW 22.8 A |
| [GALEX J0039-0356](docs/objects/2527617665632689024.md) ★★★ | White dwarf, G = 18.907; catalogued: SnowWhite DA; SIMBAD WD*/DA; MWDD DA. | Ca II emission, EW 34.42 A |
| [SDSS J2054+1610](docs/objects/1764314497240770176.md) ★★★ | White dwarf, G = 18.411; catalogued: SnowWhite DA; SIMBAD WD?. | Ca II emission, EW 24.15 A |
| [WDJ1448+3225](docs/objects/1283510882895711872.md) ★★★ | White dwarf, G = 19.383; catalogued: SIMBAD WD*/DBA; MWDD DBA; DESI DR1 DBA. | Ca II emission (double-peaked) |
| [WDJ1611+4017 (tentative)](docs/objects/1379988076130545536.md) ★★ | White dwarf, G = 17.753; catalogued: SnowWhite DA; SIMBAD WD*/DA; MWDD DA; DESI DR1 DA. | Ca II emission (single-peaked or blended) |

### Carbon lines (9)

| object | description | measurement |
|---|---|---|
| [GALEX J073504.2-794409](docs/objects/5208047381438507520.md) ★★★ | White dwarf, G = 16.564; catalogued: SIMBAD WD*/DA; MWDD DA. | C I/C II lines (CCF 16.7/3.4) |
| [GALEX J212643.1-494857](docs/objects/6466745168812781568.md) ★★ | White dwarf, G = 19.438; catalogued: SIMBAD WD?. | C I/C II lines (CCF 15.0/2.6) |
| [Gaia DR3 5836110898905253760](docs/objects/5836110898905253760.md) ★★★ | White dwarf, G = 17.441; catalogued: SIMBAD PM*. | C I/C II lines (CCF 7.5/8.6) |
| [GALEX J205119.1-161749](docs/objects/6886051830805052288.md) ★★★ | White dwarf, G = 17.548; catalogued: SIMBAD WD*/DA; MWDD DA. | C I/C II lines (CCF 18.2/7.7) |
| [Gaia DR3 883885440381808000](docs/objects/883885440381808000.md) ★★★ | White dwarf, G = 18.654; catalogued: SIMBAD WD?. | C I/C II lines (CCF 7.8/7.3) |
| [GALEX J031529.6-443716](docs/objects/4847399905305694080.md) ★★★ | White dwarf, G = 19.7; catalogued: SIMBAD WD?. | C I/C II lines (CCF 4.5/11.2) |
| [Gaia DR3 2076678981825545088](docs/objects/2076678981825545088.md) ★★★ | White dwarf, G = 18.698; catalogued: SIMBAD WD*/DA; MWDD DA. | C I/C II lines (CCF 3.8/4.6) |
| [GALEX J213644.9-515758](docs/objects/6465542891501713408.md) ★★ | White dwarf, G = 18.644; catalogued: SIMBAD WD*/DC:; MWDD DC:. | C I/C II lines (CCF 1.6/10.5) |
| [GALEX J014648.4+400114](docs/objects/343958710690034944.md) ★★ | White dwarf, G = 18.645; catalogued: SIMBAD WD*/DB; MWDD DB. | C I/C II lines (CCF 3.1/5.9) |

### Zeeman splitting (30)

| object | description | measurement |
|---|---|---|
| [GALEX J095130.1-245723](docs/objects/5660016586818424832.md) ★★ | DA white dwarf, G = 17.165; catalogued: SnowWhite DAH/DA; SIMBAD WD*/DA; MWDD DA. | B = 8.46 MG |
| [GALEX J231613.4-552927](docs/objects/6499095244738784128.md) ★★ | DA white dwarf, G = 16.701; catalogued: SnowWhite DAH; SIMBAD WD*/DA; MWDD DA. Common proper motion with HD 219458 at 164 arcsec. | B = 8.0 MG |
| [Gaia DR3 1980205739970324224](docs/objects/1980205739970324224.md) ★★★★ | DA white dwarf, G = 17.057; catalogued: SnowWhite DA; SIMBAD WD*/DA; MWDD DA:. Listed as ZZ Ceti (P = 1286 s) in Vincent et al. 2020, AJ 160, 252. | B = 5.64 MG |
| [Gaia DR3 4307667617377160704](docs/objects/4307667617377160704.md) ★★ | DA white dwarf, G = 17.895; catalogued: SnowWhite DA/DAH; SIMBAD WD*/DA; MWDD DA:. | B = 8.45 MG |
| [GALEX J091732.4-185850](docs/objects/5680077137810906624.md) ★★ | DA white dwarf, G = 17.615; catalogued: SnowWhite DAH; SIMBAD WD*/DA; MWDD DA:. | B = 6.18 MG |
| [GALEX J190543.6-114356](docs/objects/4198738558061020928.md) ★★ | DA white dwarf, G = 17.858; catalogued: SnowWhite DA/DAH; SIMBAD WD*/DA; MWDD DA. | B = 10.64 MG |
| [Gaia DR3 5807585134758743040](docs/objects/5807585134758743040.md) ★★ | DA white dwarf, G = 17.693; catalogued: SnowWhite DA; SIMBAD WD*/DA:; MWDD DA:. | B = 8.77 MG |
| [GALEX J060732.0-390139](docs/objects/2883030508640778752.md) ★★ | DA white dwarf, G = 17.582; catalogued: SnowWhite DAH; SIMBAD WD*/DA:; MWDD DA:. | B = 8.83 MG |
| [GALEX J215736.7-574324](docs/objects/6412133010376770560.md) ★★ | DA white dwarf, G = 18.645; catalogued: SnowWhite DAH; SIMBAD WD*/DA; MWDD DA. | B = 7.96 MG |
| [GALEX J203016.2-620507](docs/objects/6430762242043644032.md) ★★ | DA white dwarf, G = 18.771; catalogued: SnowWhite DAH/DA; SIMBAD WD*/DA; MWDD DA. | B = 6.96 MG |
| [Gaia DR3 5848754492268362624](docs/objects/5848754492268362624.md) ★★ | DA white dwarf, G = 17.888; catalogued: SnowWhite DAH/DA; SIMBAD WD*/DA; MWDD DA:. | B = 7.11 MG |
| [Gaia DR3 428300220431345536](docs/objects/428300220431345536.md) ★★ | DA white dwarf, G = 17.572; catalogued: SnowWhite DA/DAH; SIMBAD WD*/DA; MWDD DA:. | B = 6.51 MG |
| [GALEX J195006.0+015739](docs/objects/4241409569220727424.md) ★★ | DA white dwarf, G = 18.226; catalogued: SnowWhite DA/DAH; SIMBAD WD*/DA; MWDD DA. | B = 9.74 MG |
| [Cl* Melotte   25  REIDA     503](docs/objects/50526755482496000.md) ★★ | DA white dwarf, G = 18.333; catalogued: SnowWhite DAH; SIMBAD WD*/DA:; MWDD DA:. | B = 10.02 MG |
| [GALEX J155847.0+165740](docs/objects/1199447137276626048.md) ★★ | DA white dwarf, G = 18.283; catalogued: SnowWhite DAH; SIMBAD WD*/DA; MWDD DA. | B = 5.48 MG |
| [Gaia DR3 4236646794083461120](docs/objects/4236646794083461120.md) ★★ | DA white dwarf, G = 18.069; catalogued: SnowWhite DA/DAH; SIMBAD WD*/DA; MWDD DA. | B = 5.44 MG |
| [GALEX J224743.1+202608](docs/objects/2833867392391927936.md) ★★ | DA white dwarf, G = 17.688; catalogued: SnowWhite DAH/DA; SIMBAD WD*/DA:; MWDD DA:. | B = 10.34 MG |
| [Gaia DR3 3115382600062991872](docs/objects/3115382600062991872.md) ★★ | DA white dwarf, G = 18.412; catalogued: SnowWhite DA; SIMBAD WD*/DA; MWDD DA. | B = 4.67 MG |
| [Gaia DR3 2069622487994113408](docs/objects/2069622487994113408.md) ★★ | DA white dwarf, G = 17.58; catalogued: SnowWhite DAH/DA; SIMBAD WD*/DA:; MWDD DA:. | B = 6.81 MG |
| [GALEX J083837.0+161140](docs/objects/657989024107287168.md) ★★ | DA white dwarf, G = 19.111; catalogued: SnowWhite DA/DAH; SIMBAD WD?. | B = 8.62 MG |
| [GALEX J035849.5-605233](docs/objects/4680104332056706304.md) ★★ | DA white dwarf, G = 18.697; catalogued: SnowWhite DAH/DA; SIMBAD WD*/DA; MWDD DA. | B = 6.04 MG |
| [GALEX J052137.6-425430](docs/objects/4800536829945050752.md) ★★ | DA white dwarf, G = 19.525; catalogued: SnowWhite DA/DAH; SIMBAD WD?. | B = 6.93 MG |
| [Gaia DR3 4187365308538865792](docs/objects/4187365308538865792.md) ★★ | DA white dwarf, G = 19.421; catalogued: SnowWhite DA/DAH; SIMBAD WD*/DA; MWDD DA. | B = 5.47 MG |
| [GALEX J004807.2+433940](docs/objects/375892788968158848.md) ★★ | DA white dwarf, G = 19.121; catalogued: SnowWhite DAH; SIMBAD WD*/DA:; MWDD DA:. | B = 8.76 MG |
| [GALEX J033920.8-475633](docs/objects/4833309388918586624.md) ★★ | DA white dwarf, G = 19.696; catalogued: SnowWhite DA/DAH; SIMBAD WD?. | B = 7.27 MG |
| [GALEX J213436.2-521245](docs/objects/6465559379883395328.md) ★★ | DA white dwarf, G = 19.088; catalogued: SnowWhite DA; SIMBAD WD*/DA; MWDD DA. | B = 4.16 MG |
| [Gaia DR3 5665371272869153152](docs/objects/5665371272869153152.md) ★★ | DA white dwarf, G = 19.297; catalogued: SnowWhite DA; SIMBAD WD?. | B = 8.5 MG |
| [Gaia DR3 3327361677328480256](docs/objects/3327361677328480256.md) ★★ | DA white dwarf, G = 17.835; catalogued: SnowWhite DA; SIMBAD WD*/DA; MWDD DA. | B = 3.57 MG |
| [GALEX J040038.6-615458](docs/objects/4679463733391272448.md) ★★ | DA white dwarf, G = 18.514; catalogued: SnowWhite DA; SIMBAD WD*/DA; MWDD DA. Common proper motion with Gaia DR3 4679466653969618816 (G 9.83) at 86 arcsec. | B = 4.28 MG |
| [GALEX J050006.8+080244](docs/objects/3290180587821828480.md) ★★ | DA white dwarf, G = 18.906; catalogued: SnowWhite DA; SIMBAD WD*/DA; MWDD DA. | B = 5.3 MG |

### Photometric periods (18)

| object | description | measurement |
|---|---|---|
| [GALEX J060343.7-380911](docs/objects/2883364038621038208.md) ★★ | White dwarf selected by its Gaia DR3 GLS frequency (GALEX J060343.7-380911). | P = 10.80225 h, 5.0% (ATLAS) |
| [GALEX J043613.3+383720](docs/objects/178685757799822080.md) ★★ | White dwarf selected by its Gaia DR3 GLS frequency (GALEX J043613.3+383720). | P = 175.15604 h, 5.0% (ZTF) |
| [GALEX J211204.8-571801](docs/objects/6456720612064924928.md) ★★ | White dwarf selected by its Gaia DR3 GLS frequency (GALEX J211204.8-571801). | P = 1.02224 h, 3.5% (ATLAS) |
| [GALEX J054140.8-362248](docs/objects/2888030331609338240.md) ★★ | White dwarf selected by its Gaia DR3 GLS frequency (GALEX J054140.8-362248). | P = 16.35274 h, 2.4% (ATLAS) |
| [GALEX J124819.8-261413](docs/objects/3496637913394359680.md) ★★ | White dwarf selected by its Gaia DR3 GLS frequency (GALEX J124819.8-261413). | P = 141.25536 h, 4.5% (ATLAS) |
| [WDJ025503.24+475833.96](docs/objects/437628614520520320.md) ★★ | White dwarf selected by its Gaia DR3 GLS frequency (WDJ025503.24+475833.96). | P = 121.00599 h, 4.5% (ZTF) |
| [GALEX J191430.4-572023](docs/objects/6639666736903611136.md) ★★ | White dwarf selected by its Gaia DR3 GLS frequency (GALEX J191430.4-572023). | P = 89.07682 h, 1.8% (ATLAS) |
| [SDSS J102251.62+161151.6](docs/objects/3890059941364406144.md) ★★★ | DA with narrow Balmer lines; Gentile Fusillo et al. (2021) H-atmosphere fit 22,100 K, 0.32 Msun. | P = 1.45554 h, 2.9% (ZTF) |
| [WDJ072009.19+464840.48](docs/objects/974895286283420160.md) ★ | DA at 51 pc; the TESS period is listed by Oliveira da Rosa et al. (2024). | P = 19.14313 h, 1.8% (ZTF) |
| [GALEX J003750.4+190136](docs/objects/2795150147707769728.md) ★★★ | DA; Gentile Fusillo et al. (2021) H-atmosphere fit 27,000 K, 1.08 Msun. | P = 18.24064 h, 3.0% (ZTF) |
| [Gaia DR3 6170660401283991680](docs/objects/6170660401283991680.md) ★★★ | Gaia XP class DO (Vincent et al. 2024). | P = 27.00552 h, 3.2% (ATLAS) |
| [GALEX J132200.5-422412](docs/objects/6136817910121524096.md) ★★ | Hot-subdwarf candidate in Geier et al. (2019); GALEX J132200.5-422412. | P = 18.45593 h, 3.5% (ATLAS) |
| [Gaia DR3 3123625093275668736](docs/objects/3123625093275668736.md) ★★ | White dwarf and M dwarf in Rebassa-Mansergas et al. (2025); SDSS-V SnowWhite DA_MS. | P = 10.39275 h, 2.9% (ZTF) |
| [Gaia DR3 3354819845628139904](docs/objects/3354819845628139904.md) ★★ | Hot-subdwarf candidate in Geier et al. (2019). | P = 12.58907 h, 4.0% (ZTF) |
| [WDJ043832.74+003117.01](docs/objects/3230486971974872192.md) ★★★ | DO white dwarf, Teff 105.6 kK, log g 8.0 (Kilic et al. 2026a, via the MWDD). | P = 26.09254 h, 2.6% (ZTF) |
| [WDJ080026.64+633414.85](docs/objects/1094376947131876352.md) ★★★ | DESI DR1 class DOA, Teff 78.2 kK (Swan et al. 2026). | P = 21.41901 h, 2.1% (ZTF) |
| [Gaia DR3 6021870154194477312](docs/objects/6021870154194477312.md) ★★★★ | White dwarf with a 103.4-min period; three Gaia sources within 13 arcsec are fitted separately. | 103.4-min period (ATLAS, Gaia, TESS) |
| [WDJ064438.09-004550.51](docs/objects/3107374277060584064.md) ★★★★ | Hot white dwarf with He II 4686 absorption (SIMBAD WD* DO:, MWDD DO:, SDSS-V SnowWhite DA:; VSX type WD without a period); H-alpha, H-beta and Ca II emission whose velocity follows the photometric phase. | P = 14.229 h; emission follows the phase |

### TESS amplitude spectra (pulsation candidates) (5)

| object | description | measurement |
|---|---|---|
| [GALEX J054243.4-261011](docs/objects/2908195134345338496.md) ★★ | White dwarf, G = 17.634; SnowWhite DA (Teff 12173 K, log g 7.93); catalogued: SIMBAD WD*/DA; MWDD DA. | pulsations, 692.1 s |
| [GALEX J033619.1-564435](docs/objects/4729763229265811328.md) ★★ | White dwarf, G = 17.493; SnowWhite DA (Teff 11626 K, log g 8.07); catalogued: SIMBAD WD*/DA; MWDD DA. | pulsations, 605.1 s |
| [L  210-25](docs/objects/6472153670805656832.md) ★★ | White dwarf, G = 16.402; SnowWhite DA (Teff 11191 K, log g 8.16); catalogued: SIMBAD WD*/DA; MWDD DA. | pulsations, 817.3 s |
| [[OHD2001] WD J2324-595](docs/objects/6492083311194727168.md) ★★★ | White dwarf, G = 16.809; SnowWhite DA (Teff 11600 K, log g 7.89); catalogued: SIMBAD WD*/DA; MWDD DA. | pulsations, 170.4 s |
| [GALEX J214927.5-515827](docs/objects/6558472750993181568.md) ★★★ | White dwarf, G = 17.084; SnowWhite DA (Teff 11817 K, log g 8.12); catalogued: SIMBAD WD*/DA; MWDD DA. | pulsations, 1067.1 s |

### Eclipse and Balmer emission (3)

| object | description | measurement |
|---|---|---|
| [Gaia DR3 4731701084150029824](docs/objects/4731701084150029824.md) ★★★ | White dwarf, G = 18.373; catalogued: none found. | eclipses, P = 0.14786971 d |
| [Gaia DR3 2002597083798483200](docs/objects/2002597083798483200.md) ★★ | G = 18.816; catalogued: SIMBAD WD*/DQ:; MWDD DQ:. | Balmer emission, EW 318.3 A |
| [Gaia DR3 1977447164064222976](docs/objects/1977447164064222976.md) ★★ | G = 19.454; catalogued: SIMBAD WD*/DC:; MWDD DC:. | Balmer emission, EW 45.9 A |

### Hot white dwarfs with He II lines (6)

| object | description | measurement |
|---|---|---|
| [GALEX J055029.7-155446](docs/objects/2995107164834343680.md) ★★★ | Hot white dwarf; no earlier spectrum found. | DAO, Teff ~60-110 kK (TMAP) |
| [GALEX J062928.9-415857](docs/objects/5570041179495992704.md) ★★★ | Hot white dwarf; no earlier spectrum found. | DAO, Teff ~60-110 kK (TMAP) |
| [SDSS J081413.43+022524.7](docs/objects/3090786872841030016.md) ★★★ | Hot white dwarf; no earlier spectrum found. | DAO, Teff ~60-110 kK (TMAP) |
| [Gaia DR3 4036084504408126976](docs/objects/4036084504408126976.md) ★★★ | Hot white dwarf; no earlier spectrum found. | DAO, Teff ~60-110 kK (TMAP) |
| [GALEX J190659.9-755815](docs/objects/6365804611201098368.md) ★★★ | Hot white dwarf; no earlier spectrum found. | DAO, Teff ~60-110 kK (TMAP) |
| [WDJ095852.35-175833.41](docs/objects/5671975077144346112.md) ★★★ | Hot white dwarf; no earlier spectrum found. | DAO, Teff ~60-110 kK (TMAP) |

### Short-period white dwarfs with irradiated companions (6)

| object | description | measurement |
|---|---|---|
| [WDJ205249.27-032419.53](docs/objects/6914922055508553984.md) ★★★★ | Low-mass white dwarf: GF21 H-atmosphere Teff 16231 K, 0.301 Msun; G = 17.466, 334 pc. | P = 97.703 min; W1 3.74x model |
| [WDJ212738.67+593755.72](docs/objects/2191618770599895296.md) ★★★★ | Low-mass white dwarf: GF21 H-atmosphere Teff 13777 K, 0.263 Msun; G = 16.871, 241 pc. | P = 130.182 min; W1 3.61x model |
| [WDJ070106.16-534811.37](docs/objects/5503429908930455808.md) ★★★★ | Low-mass white dwarf: GF21 H-atmosphere Teff 16240 K, 0.256 Msun; G = 18.426, 602 pc. | P = 81.493 min; W1 2.21x model |
| [WDJ040444.35-395043.1](docs/objects/4844023064578952320.md) ★★★★ | Low-mass white dwarf: GF21 H-atmosphere Teff 20562 K, 0.243 Msun; G = 17.663, 632 pc. | P = 117.319 min; W1 2.02x model |
| [WDJ194901.41+673005.59](docs/objects/2249098833310553728.md) ★★★★★ | Low-mass white dwarf: GF21 H-atmosphere Teff 16817 K, 0.143 Msun; G = 18.027, 739 pc. | P = 63.68 min; W1 x model |
| [WDJ005615.18-661731.98](docs/objects/4705562733524591232.md) ★★★ | Low-mass white dwarf: GF21 H-atmosphere Teff 15358 K, 0.255 Msun; G = 18.801, 714 pc. | P = 73.522 min; W1 1.69x model |

### Day-scale periods of hot white dwarfs (9)

| object | description | measurement |
|---|---|---|
| [WDJ091433.60+581238.12](docs/objects/1038176780370360576.md) ★★★ | SBSS 0910+584; spectral type DO (MWDD). G = 17.73, 810 pc. | P = 27.23048 h, 3.8% (ZTF) |
| [WDJ151215.73+065156.43](docs/objects/1157401396015448960.md) ★ | GALEX J151215.7+065156; UHE white dwarf, DOZ (Reindl et al. 2021, where the period is published). G = 17.22, 990 pc. | P = 5.4245 h, 1.8% (ZTF) |
| [WDJ221519.86+253059.05](docs/objects/1879989790567353344.md) ★★★ | GALEX J221519.8+253059; spectral type DOZ (MWDD). G = 17.05, 1190 pc. | P = 18.54001 h, 2.2% (ZTF) |
| [WDJ075540.94+400917.91](docs/objects/920621124593362816.md) ★ | KUV 07523+4017; DOZ / PG 1159 (Reindl et al. 2021, where the period is published). G = 17.80, 1052 pc. | P = 20.78523 h, 2.5% (ZTF) |
| [WDJ065819.86+441438.40](docs/objects/953685015492787456.md) ★★★ | Gaia XP class DO (Vincent et al. 2024). G = 17.53, 483 pc. | P = 14.74481 h, 2.8% (ZTF) |
| [WDJ025657.85-145029.92](docs/objects/5157333438398813824.md) ★★★ | Spectral type DA with a photometric temperature above 100 kK (MWDD). G = 17.37, 1072 pc. | P = 28.51831 h, 3.4% (ZTF) |
| [WDJ171743.52+515840.07](docs/objects/1415911839725510528.md) ★★ | Spectral type DA, 68.0 kK (Kilic et al. 2026, via the MWDD). G = 17.59, 833 pc. | P = 30.1191 h, 1.9% (ZTF) |
| [WDJ234931.84-353916.52](docs/objects/2311285729210966144.md) ★★ | GALEX J234931.8-353916; SDSS-V DR20 SnowWhite class DA. G = 17.88, 882 pc. | P = 26.4793 h, 2.6% (ATLAS) |
| [WDJ065134.01+185201.09](docs/objects/3365371721281530880.md) ★★ | No spectrum found. G = 18.10, 1100 pc. | P = 20.17182 h, 3.6% (ATLAS) |

### Periods of DA white dwarfs with emission lines (7)

| object | description | measurement |
|---|---|---|
| [WDJ172406.13+562003.08](docs/objects/1420761029600606592.md) ★★★ | DESI DR1 class DAe; GF21 15.0 kK, 0.11 Msun; W1, W2 excess 1.8, 1.9 mag. G = 16.26, 403 pc. | P = 5.99166 h, 5.8% (ZTF) |
| [WDJ170125.28+343530.59](docs/objects/1337970174853051392.md) ★★★ | DESI DR1 class DAE; GF21 12.1 kK, 0.17 Msun; W1, W2 excess 2.0, 1.7 mag. G = 18.94, 741 pc. | P = 2.73751 h, 3.9% (ZTF) |
| [WDJ161752.97+015840.92](docs/objects/4409006786607484672.md) ★★★ | DESI DR1 class DAE; GF21 22.4 kK, 0.24 Msun; W1, W2 excess 2.1, 2.6 mag. G = 18.95, 1233 pc. | P = 2.04723 h, 5.7% (ZTF) |
| [WDJ132308.63+055900.97](docs/objects/3717349170269867520.md) ★★★ | DESI DR1 class DAe; GF21 22.8 kK, 0.34 Msun; W1, W2 excess 1.9, 2.0 mag. G = 18.99, 955 pc. | P = 2.68284 h, 7.4% (ZTF) |
| [WDJ214656.86+143125.17](docs/objects/1769157090045264128.md) ★★★ | DESI DR1 class DAE; GF21 16.0 kK, 0.24 Msun; W1, W2 excess 1.9, 1.6 mag. G = 19.45, 1010 pc. | P = 6.83088 h, 3.9% (ZTF) |
| [WDJ083531.69+315503.31](docs/objects/709815329316284928.md) ★★★ | DESI DR1 class DAE; GF21 13.2 kK, 0.25 Msun; W1, W2 excess 2.1, 2.4 mag. G = 19.26, 742 pc. | P = 4.42399 h, 4.9% (ZTF) |
| [WDJ075449.34+442357.52](docs/objects/926161868627454976.md) ★★★ | DESI DR1 class DAe; GF21 18.9 kK, 0.21 Msun; W1, W2 excess 2.2, 2.9 mag. G = 19.40, 1449 pc. | P = 3.82105 h, 3.9% (ZTF) |
<!-- object-index:end -->

<img src="figures/gas_discs/578709631539357440_epochs.png" width="720">

<img src="figures/carbon/hot_dq_comparison_sdssv.png" width="430"> <img src="figures/zeeman/overview.png" width="430">

## Layout
- `tables/`: measurement tables (CSV).
- `figures/<topic>/`: one figure per object or per comparison.
- `docs/`: one page per topic; `docs/objects/`: one page per object.
- `data/`: inputs (ATLAS and ZTF photometry, Gaia epoch photometry, object lists); `data/cache/` holds downloads and is not tracked.
- `scripts/`: measurement and figure scripts.

## Reproduction
```
pip install -r requirements.txt
cd scripts
python zeeman_split.py ../data/zeeman_input.csv
python carbon_lines.py 95077848 5208047381438507520
python carbon_lines.py 110600288 6466745168812781568
python carbon_lines.py 102600838 5836110898905253760
python carbon_lines.py 114554634 6886051830805052288
python carbon_lines.py 57623143 883885440381808000
python carbon_lines.py 93091478 4847399905305694080
python carbon_lines.py 67111869 2076678981825545088
python carbon_lines.py 110590717 6465542891501713408
python carbon_lines.py 116892932 343958710690034944
python carbon_screen.py --table
python carbon_screen.py --sample ../tables/carbon_screen_sample.csv   # full 3,480-spectrum screen (about 1 GB of downloads)
python carbon_screen.py --sample-nonda ../tables/carbon_screen_nonda.csv   # 3,286 non-DA spectra
python lamost_compare.py
python cos_lines.py
python galex_colours.py 5208047381438507520
python galex_colours.py 6886051830805052288
python galex_colours.py 4847399905305694080 --window -0.25 0.05 --target-galex 22.858 0.203 19.731 0.019
python tess_periodogram.py 2055170284 102 120
python tess_pixel_test.py 6492083311194727168 2055170284 102 82.34 120
python j0353_eclipse.py
python cv_balmer.py 65701864
python periodic_6021870154194477312.py
python tess_periodogram.py 1251484163 65 120 0.5 50
python periodic_white_dwarfs.py
python hot_dae_wd_periods.py
python reflection_3107374277060584064.py
python tess_ffi_photometry.py 3890059941364406144 155.7148994 16.1977169 16.488695 18.034 45 46 72
python hot_white_dwarfs.py   # downloads 280 TheoSSA model spectra (about 1.4 GB)
python irradiated_companions.py
python desi_halpha_6914922055508553984.py   # after irradiated_companions.py (uses its ephemeris)
python gas_disc_screen.py --sample   # about 51,000 visit files, downloaded and deleted one by one; keeps about 1.5 GB
python gas_disc_screen.py --pass2 ../data/cache/gas_disc_pass2.csv
python gas_disc_screen.py --table ../data/cache/gas_disc_pass2.csv
python gas_disc_epochs.py 578709631539357440 55774610 134.841202 4.636784
python gas_disc_epochs.py 1827014701883095680 63867520 299.804316 22.147851
python gas_disc_epochs.py 2527617665632689024 70254122 9.892213 -3.946573
python gas_disc_epochs.py 1764314497240770176 63203321 313.742432 16.179092
python desi_gas_disc_screen.py --sample   # 44,417 DESI DR1 spectra via SPARCL; keeps about 0.5 GB
python desi_gas_disc_screen.py --pass2 ../data/cache/desi_gas_pass2.csv
python desi_gas_disc_screen.py --table ../data/cache/desi_gas_pass2.csv
python gas_disc_epochs.py 1283510882895711872 - 222.081216 32.416845
python gas_disc_epochs.py 1379988076130545536 60943632 242.822993 40.284241
python wise_excess.py 2527617665632689024 9.892213 -3.946573
python wise_excess.py 578709631539357440 134.841202 4.636784
python ztf_lightcurve.py 578709631539357440 134.841202 4.636784
python ztf_lightcurve.py 1827014701883095680 299.804316 22.147851
python ztf_lightcurve.py 2527617665632689024 9.892213 -3.946573
python ztf_lightcurve.py 1764314497240770176 313.742432 16.179092
python figures.py            # all figures; or name one, e.g. python figures.py gas_discs
python object_pages.py       # object pages and the README object index
```
Arguments for the other TESS light curves and pixel tests are in the table columns (TIC, sector, cadence, frequency).

## Data sources
SDSS-V DR20 and SDSS DR17 (including BOSS); DESI DR1 via SPARCL (NOIRLab Astro Data Lab) and the DESI DR1 white-dwarf catalogues (Swan et al. 2026; Amorim et al. 2026); ESO X-shooter phase 3 spectra (programme 115.28GM.001); Legacy Surveys DR10 (Astro Data Lab); Gaia DR3 (ESA/Gaia/DPAC), including epoch photometry (VizieR I/355/epphot), the Gaia Synthetic Photometry Catalogue (VizieR J/A+A/674/A33) and Gentile Fusillo et al. (2021, VizieR J/MNRAS/508/3877); TESS SPOC light curves and full-frame images (TESScut), and HST/COS program 17420 (MAST); CoRoT faint-star light curves (CDS, B/corot); GALEX GUVcat AIS (Bianchi et al. 2017), GALEX GR6/7 (MAST) and the GALEX CAUSE Kepler catalogue (Olmedo et al. 2015); LAMOST DR10; ATLAS forced photometry (Tonry et al. 2018; Shingles et al. 2021); ZTF public data releases (IRSA); CatWISE2020 (Marocco et al. 2021) and VISTA Hemisphere Survey DR5 (VizieR II/365, II/367); Montreal synthetic photometry of pure-H white dwarfs (Holberg & Bergeron 2006; Bédard et al. 2020); DESI DR1 fibre-assignment and redshift tables (Astro Data Lab); NIST Atomic Spectra Database; TMAP H+He model spectra from TheoSSA (GAVO Data Center); Montreal White Dwarf Database; SIMBAD.
