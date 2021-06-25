extensions[gis]

globals
[
  year
  num                              ; total number of fish (1000s)
  ssb                              ; spawning stock biomass (tonnes)
  tsb                              ; total stock biomass (tonnes)
  rec
  annual_num
  annual_ssb
  annual_tsb
  num_overone
  icesnum                          ; assessment numbers
  icesssb                          ; assessment spawning stock biomass
  icestsb                          ; assessment total stock biomass
  sst-data                         ; temperature data
  land
  directory
  load-data
  t0                               ; von Bertalanffy parameter, hypothetical age at length 0 (years)
  linf                             ; von Bertalanffy parameter, asymptotic length (cm)
  K                                ; von Bertalanffy parameter, growth rate coefficient
  A                                ; aspect ratio of the caudal fin
  y                                ; recruitment
  Rest                             ; assessment recruitment (number)
  SR_SSB                           ; begin year SSB (used for stock recruit relationship)
  boltz                            ; Boltzmann constant (J/K)
  EaS                              ; activation energy for speed arrhenius equation (J)
  refTS                            ; reference temperature for speed arrhenius equation (K)
  Nfeb                             ; connectivity parameter
  Nmar                             ; connectivity parameter
  Napr                             ; connectivity parameter
  Nmay                             ; connectivity parameter
  nspawn                           ; connectivity parameter
  n                                ; connectivity parameter
  Locfeb                           ; connectivity parameter
  Locmar                           ; connectivity parameter
  Locapr                           ; connectivity parameter
  Locmay                           ; connectivity parameter
  loc                              ; connectivity parameter
  raster
  rasterSST
  phyto-data
  coast-data
  ncohort
  recruitment
  Ea                             ; Activation energy (eV)
  Tref                           ; Reference temperature for the energy budget
  Cmax                           ; maximum ingestion rate g/day
  ep                             ; energy content of phytoplankton kJ/g
  ae                             ; Assimilation efficiency
  A0                             ; normalizing constant for AMR
  Ef                             ; energy content of flesh kJ/g
  El                             ; energy content of lipid kJ/g
  Fs                             ; energy costs of synthesising flesh kJ/g
  Ls                             ; energy costs of synthesising lipid kJ/g
  egg_mass
  larval-production ; the number of eggs that survive to larval stage each year
  Lm                            ; threshold length (cm) for maturity
  num-matured                   ; number of individuals that have reached sexual maturity (reset each year)
  rolling_tick
  a_g;Length-weight parameter
  b_g; Length-weight parameter
  eggs_per_bass
  boltzS
  year_spin
  H
  Am
  Pm
  I
  GL


  recruits_ABC
  L0
  L1
  L2
  L3
  L4
  L5
  L6
  L7
  L8
  L9
  L10
  L11
  L12
  L13
  L14
  L15
  L16
  L17
  L18
  L19
  L20
  L21
  L22
  L23
  L24
  L25
  L26
  L27
  L28
  L29
  L30

  N0
  N1
  N2
  N3
  N4
  N5
  N6
  N7
  N8
  N9
  N10
  N11
  N12
  N13
  N14
  N15
  N16
  N17
  N18
  N19
  N20
  N21
  N22
  N23
  N24
  N25
  N26
  N27
  N28
  N29
  N30



  M0
  M1
  M2
  M3
  M4
  M5
  M6
  M7
  M8
  M9
  M10
  M11
  M12
  M13
  M14
  M15
  M16
  M17
  M18
  M19
  M20
  M21
  M22
  M23
  M24
  M25
  M26
  M27
  M28
  M29
  M30


  MT0
  MT1
  MT2
  MT3
  MT4
  MT5
  MT6
  MT7
  MT8
  MT9
  MT10
  MT11
  MT12
  MT13
  MT14
  MT15
  MT16
  MT17
  MT18
  MT19
  MT20
  MT21
  MT22
  MT23
  MT24
  MT25
  MT26
  MT27
  MT28
  MT29
  MT30



  annual_rec
  run-year
  SumSSB



ABC_Rec
ABC_SSB
ABC_Num

ABC_M0
ABC_M1
ABC_M2
ABC_M3
ABC_M4
ABC_M5
ABC_M6
ABC_M7
ABC_M8
ABC_M9
ABC_M10
ABC_M11
ABC_M12
ABC_M13
ABC_M14
ABC_M15
ABC_M16
ABC_M17
ABC_M18
ABC_M19
ABC_M20
ABC_M21
ABC_M22
ABC_M23
ABC_M24
ABC_M25
ABC_M26
ABC_M27
ABC_M28
ABC_M29
ABC_M30


ABC_N0
ABC_N1
ABC_N2
ABC_N3
ABC_N4
ABC_N5
ABC_N6
ABC_N7
ABC_N8
ABC_N9
ABC_N10
ABC_N11
ABC_N12
ABC_N13
ABC_N14
ABC_N15
ABC_N16
ABC_N17
ABC_N18
ABC_N19
ABC_N20
ABC_N21
ABC_N22
ABC_N23
ABC_N24
ABC_N25
ABC_N26
ABC_N27
ABC_N28
ABC_N29
ABC_N30



]

patches-own
[
  sst
  phyto
  coast-patches
  nursery-patches
  offshore-patches
  spawn-patches
  sea
  IVb                              ; label for coastal patches within a particular ICES sea area
  IVc                              ; see above
  VIId                             ; see above
  VIIe                             ; see above
  VIIfg                            ; see above
  VIIa                             ; see above
  CS                               ; label for potential offshore spawning patches in the Celtic Sea
  C                                ; label for potential offshore spawning patches in the Channel
  NS                               ; label for potential offshore spawning patches in the North Sea
  IS                               ; label for potential offshore spawning patches in the Irish Sea

  catch
  chl
  row
  column
  rectangle
  IVb_dist
  IVb_target
  IVc_dist
  IVc_target
  VIId_dist
  VIId_target
  VIIe_dist
  VIIe_target
  VIIfg_dist
  VIIfg_target
  VIIa_dist
  VIIa_target
  processed
  sst-output

]

turtles-own
[
  age                              ; actual age in days
  cohort                           ; cohort age (age 1st January)
  number                           ; number of individuals within the superindividual
  L                                ; length (cm)
  IVbt                             ; affinity to an ICES division for feeding
  IVct                             ; see above
  VIIdt                            ; see above
  VIIet                            ; see above
  VIIfgt                           ; see above
  VIIat                            ; see above
  spawn-trigger                    ; a switch for spawning and feeding migrations (0 = feeding; 1 = spawning)
  spawn-count                      ; counter for the number of days in a spawning ground
  speed                            ; calculated from length (units are "patch lengths")
  r                                ; number of movement repeats in a time step
  M                                ; natural mortality
  FCi                              ; commercial inshore fishing mortality
  FCo                              ; commercial offshore fishing mortality
  FRi                              ; recreational inshore fishing mortality
  ingestion-rate                   ; ingestion of a whiole super-individual (needed to calc. predation mortality) g/day
  func-response                    ; Holling type II functional response adjusts ingestion rate based on chl availability
  assimilation-energy
  MR                               ; either active (AMR) and standard (SMR) metabolic rate depending on the conditions (kJ/day)
  max-growth-rate
  std-mass                          ; standard mass as calculated from length using known relationship (g)
  structural-mass                   ; structural mass, i.e. total mass - gonad mass - lipid mass (g) (energy reserve)
  energy-reserve-max                ; kJ
  energy-reserve                    ; kJ
  growth-costs                      ; energy cost of adding new length and mass
  growth-rate                       ; realised amount of length added in a day if insufficient energy to grow maximally (cm)
  potential-fecundity               ; the potential number of eggs that can be produced
  max-R                             ; energy costs of synthesising a full season's eggs (kJ)
  maintenance-energy                ; 10% of the energy stored at the beginning of spawning is saved for maintenance costs while fasting on the spawning grounds
  gonad-mass                        ; mass of the gonads (eggs, g)
  ER
  FK                                ; condition factor = 100*(W/L^3)
  total-mass                        ; structural mass + lipid mass (energy reserve) + gonad mass (g)
  migrating                         ; boolean. Prevents individuals "moving locally" if true
  development                       ; the number of days remaining before an egg hatches
  embryo-duration
  realised-fecundity
  Wrep                             ; weight as seen by the assessment
  energy-reserve-max-mass
  reach_target
  num_near_me
  mass_near_me
  spawn_chosen


]

breed [juv-bass]                     ; bass with L < 42
breed [mat-bass]                     ; bass with L > 42
breed [eggs egg]                    ; Eggs
breed [ys-larvae]
breed [larvae  ]


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; INITIALISATION ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;




to setup

  clear-all
  reset-ticks


   set year 1985

set ABC_SSB []
set ABC_Rec []
set ABC_Num []

set ABC_M0  []
set ABC_M1  []
set ABC_M2  []
set ABC_M3  []
set ABC_M4  []
set ABC_M5  []
set ABC_M6  []
set ABC_M7  []
set ABC_M8  []
set ABC_M9  []
set ABC_M10 []
set ABC_M11 []
set ABC_M12 []
set ABC_M13 []
set ABC_M14 []
set ABC_M15 []
set ABC_M16 []
set ABC_M17 []
set ABC_M18 []
set ABC_M19 []
set ABC_M20 []
set ABC_M21 []
set ABC_M22 []
set ABC_M23 []
set ABC_M24 []
set ABC_M25 []
set ABC_M26 []
set ABC_M27 []
set ABC_M28 []
set ABC_M29 []

set ABC_N0  []
set ABC_N1  []
set ABC_N2  []
set ABC_N3  []
set ABC_N4  []
set ABC_N5  []
set ABC_N6  []
set ABC_N7  []
set ABC_N8  []
set ABC_N9  []
set ABC_N10 []
set ABC_N11 []
set ABC_N12 []
set ABC_N13 []
set ABC_N14 []
set ABC_N15 []
set ABC_N16 []
set ABC_N17 []
set ABC_N18 []
set ABC_N19 []
set ABC_N20 []
set ABC_N21 []
set ABC_N22 []
set ABC_N23 []
set ABC_N24 []
set ABC_N25 []
set ABC_N26 []
set ABC_N27 []
set ABC_N28 []
set ABC_N29 []



;data

;C:/Users/dm820646/Dropbox/From andrews PC


  gis:load-coordinate-system ("data/esriwkt.txt")
  gis:set-world-envelope-ds [ -9 9 48 57.5 ]     ; sets the coordinate boundaries of the NetLogo world
  set land gis:load-dataset "data/europe_countries.shp"
   foreach gis:feature-list-of land
   [
     gis:set-drawing-color grey
     gis:draw ? 1.0
     gis:fill ? 2.0
   ]


set rolling_tick   1



  ; Assign sst and CHL to patches
set SST-data gis:load-dataset "data/rs_data_bass/sst_1.asc"            ;See load-chl procedure for full details on import of satellite data RB
    gis:apply-raster SST-data SST

set phyto-data gis:load-dataset "data/rs_data_bass/chl_1.asc"  ;RB
    gis:apply-raster phyto-data chl

ask patches with [sst >= 0]
    [set pcolor scale-color blue sst 5 18]


    file-open "data/inputs.txt" ; where are coast patches. TRUE or false with coordinates
    while [not file-at-end?]
    [
      let coast_X file-read
     let coast_Y file-read
    let coast? file-read
      ask patch coast_X coast_Y [set coast-patches coast?]
    ]
    file-close

   ask patches with [(coast-patches = true)]
   [
   ifelse sst > 0
    []

     [set sst [sst] of min-one-of patches with [sst > 0] [distance myself]]

   ifelse (chl > 0)
 []

  [if (any? neighbors with [chl > 0])  [set chl mean [chl] of neighbors with [chl > 0 ]]]
  ]


  ask patches [                                 ; define sea and non-sea patches
        ifelse chl > 0
    [ set sea true]
    [ set sea false]
    if (pxcor <= 9 and pycor >= 28) or (pxcor <= 1 and pycor >= 18)   ; patches in VI and VIIb
    [ set pcolor black
      set sea false]
    if (pxcor >= 2 and pxcor <= 3 and pycor >= 26 and pycor <= 27) or (pxcor >= 34 and pxcor <= 35 and pycor >= 36 and pycor <= 37)
    [ set pcolor black
      set sea false]
                                ; initialise the Arrhenius factors

  ]



  ask patches with [(coast-patches != 0)]
  [
    set coast-patches true     ; extend coastal patches to ICES rectangles intersecting land
        set pcolor green
  ]

  ; Manual edits to give fish enough room to move around the coast
  ask patches with [(pxcor >= 14 and pxcor <= 15 and pycor >= 34 and pycor <= 35) or
    (pxcor >= 18 and pxcor <= 19 and pycor >= 24 and pycor <= 25) or
    (pxcor >= 20 and pxcor <= 21 and pycor >= 20 and pycor <= 21) or
    (pxcor >= 16 and pxcor <= 17 and pycor >= 6 and pycor <= 7)]
  [
    set coast-patches true
    set pcolor green
  ]


  file-open "data/land.txt"     ; patches ~ >50% land set not sea, preventing fish travelling on them
  while [ not file-at-end? ]
  [
    let landx file-read
    let landy file-read
    ask patch landx landy
    [
      set sea false
      set coast-patches 0
      set nursery-patches 0
      set pcolor black
    ]
  ]
  file-close


  ask patches with [ sea = true and coast-patches = 0 ]   ; sea patches that are not coastal are set to offshore
  [
    set offshore-patches true
  ]

  ask patches with [ offshore-patches = true ]     ; assign potential spawning patches the correct sea area
  [
    if (pxcor >= 20 and pycor < 24)
    [ set NS true]

    if (pxcor < 8 and pycor < 16)
    [ set CS true]

    if (pxcor >= 8 and pxcor < 20 and pycor < 16)
    [ set C true]

    if (pxcor < 13 and pycor >= 16 and pycor < 24)
    [ set IS true]
  ]

  ask patches with [coast-patches = true]
  [
    if ( pxcor >= 10 and pycor >= 30 ) or (pxcor >= 14 and pycor >= 26) or (pxcor >= 16 and pycor >= 22)         ; assigns patches the relevant labels depending on which ICES area they are in
    [ set IVb true]

    if ( pycor >= 12 and pycor <= 21 and pxcor > 14)
    [ set IVc true]

    if ( pxcor >= 14 and pycor >= 6 and pycor <= 11 ) or ( pxcor >= 16 and pycor >= 4 and pycor <= 11)
    [ set VIId true]

    if ( pxcor <= 13 and pycor <= 7 ) or ( pxcor >= 8 and pxcor <= 13 and pycor <= 9 ) or ( pxcor >= 10 and pxcor <= 13 and pycor <= 11 ) or (pxcor >= 14 and pxcor <= 15 and pycor <= 5)
    [ set VIIe true]

    if ( pycor <= 15 and pycor >= 12 and pxcor <= 13 ) or ( pycor >= 10 and pycor <= 11 and pxcor <= 9 ) or ( pycor >= 8 and pycor <= 9 and pxcor <= 7 )
      [set VIIfg true]

    if ( pycor >= 16 and pycor <= 27 and pxcor <= 13 )
      [ set VIIa true]
  ]


  ask patches gis:intersecting land             ; patches intersecting land set to coastal (& nursery if south of 54N) patches except odd patch on Danish coast
  [
    if (sea = true ) and (pxcor != 34 or pycor != 35)
    [ set coast-patches true
      set pcolor green


      if (pycor < 24)
      [
        set nursery-patches true
        set pcolor turquoise
      ]
    ]
  ]








;;;;;;;;;;;;;;;;;;;;;;IVb


    ask patch 18 22 ;
   [
      set IVb_dist 1
      set processed true
      set IVb_target true
      set pcolor red
                 ]

  while [10 < count patches with [(sea = true) and (processed != true)]]                    ; The while condition is set at 20 <.... because there are 20 patches that satisfy the conditions to be "shelf-edge", but are separated from the shelf by deep water
  [
  ask patches with [sea = true]           ; i.e. patch 37 26, then more on next step
  [
    if processed != true [stop]
    ask neighbors with [(processed != true) and (sea = true)]      ; ask neighbours on shelf edge that have not yet been processed ;was orginally neighbours4 takes far longer not sure why??
    [
      set IVb_dist ([IVb_dist] of myself) + 1
      set processed true

       ]
  ]
  ]
  ask patches [set processed 0]


;;;;;;;;;;;;;;;IVc

  ask patch 21 16
   [
      set IVc_dist 1
      set processed true
       set IVc_target true
      set pcolor red

   ]

  while [10 < count patches with [(sea = true) and (processed != true)]]
  [
  ask patches with [sea = true]
  [
    if processed != true [stop]
    ask neighbors with [(processed != true) and (sea = true)]
    [
      set IVc_dist ([IVc_dist] of myself) + 1
      set processed true
    ]
  ]
  ]
  ask patches [set processed 0]


;;;;;;;;;;;;;;;VIId

  ask patch 18 10
   [
      set VIId_dist 1
      set processed true
      set VIId_target true
      set pcolor red
       ]

  while [10 < count patches with [(sea = true) and (processed != true)]]
  [
  ask patches with [sea = true]
  [
    if processed != true [stop]
    ask neighbors with [(processed != true) and (sea = true)]
    [
      set VIId_dist ([VIId_dist] of myself) + 1
      set processed true

        ]
  ]
  ]
  ask patches [set processed 0]



;;;;;;;;;;;;VIIe
  ask patch 10 8
   [
      set VIIe_dist 1
      set processed true
      set VIIe_target true
      set pcolor red
       ]

  while [10 < count patches with [(sea = true) and (processed != true)]]
  [
  ask patches with [sea = true]
  [
    if processed != true [stop]
    ask neighbors with [(processed != true) and (sea = true)]
    [
      set VIIe_dist ([VIIe_dist] of myself) + 1
      set processed true

       ]
  ]
  ]

  ask patches [set processed 0]



;;;;;;;;;;;;VIIfg
 ask patch 11 13
   [
      set VIIfg_dist 1
      set processed true
      set VIIfg_target true
      set pcolor red
       ]

  while [10 < count patches with [(sea = true) and (processed != true)]]
  [
  ask patches with [sea = true]
  [
    if processed != true [stop]
    ask neighbors with [(processed != true) and (sea = true)]
    [
      set VIIfg_dist ([VIIfg_dist] of myself) + 1
      set processed true]
  ]
  ]

  ask patches [set processed 0]



;;;;;;;;;;;;VIIa
  ask patch 11 21
   [  set VIIa_dist 1
      set processed true
      set VIIa_target true
      set pcolor red
       ]

  while [10 < count patches with [(sea = true) and (processed != true)]]
  [
  ask patches with [sea = true]
  [ if processed != true [stop]
    ask neighbors with [(processed != true) and (sea = true)]
    [ set VIIa_dist ([VIIa_dist] of myself) + 1
      set processed true]
  ]
  ]

  ask patches [set processed 0]



     ; set parameter values
  set linf 84.55                      ; von bertelanfy asymptotic length, value as used in stock assessment (cm)
  set k  0.09699                       ; von bertelanfy growth parameter, value used in stock assessment
  set t0 -0.73                        ; von bertelanfy, to, hypothetical age at length 0, value taken from stock assessment (years)
  set refTS 279                       ; reference temperature 6 degrees, the temperature at which our swimming speed at the length used in the paper coincides with that from the model in "Effect of temperature on maximum swimming speed and cost of transport in juvenile European sea bass (Dicentrarchus labrax" (K)
  set A 1.76                          ; aspect ratio of the caudal fin, taken from FishBase
  set Ea 0.5                          ; Activation energy (eV)
  set EaS 0.1903656  ;activation enenrgy in ev same working from claurmx paper see og trace
  set Tref 285.15                     ; Reference temperature for the energy budget (12oC)
  set boltz (8.62 * (10 ^ -5))
  set Cmax 0.54                       ; bass (Lanari, D’Agaro and Ballestrazzi, 2012)
  set ep 6.02                         ; energy content of phytoplankton kJ/g (Annis et al., 2011)
  set A0 0.1227808                    ; bass normalizing constant calc in most recent C+R script
  set Ef 7         ;Kj/g                  ; (Peters 1983)
  set El 39.3  ;kj/g                       ; (Schmidt-nielsen 1997)
  set Ls 14.7                        ; Pullar and Webster (1977), see TRACE
  set Fs 3.6                            ; (Sibly et al., 2013; Sibly and Calow, 1986)
  set egg_mass 0.96 * 10 ^ -3   ;value here in grams    ;quoted in paper in miligrams 0.96 mg         ; bass (Cerdá et al., 1994)
  set ncohort 10
  set a_g (0.00001296 *  0.95) ;1.296 * (10 ^ -5) ;Length-weight parameter (ICES, 2012a) Report of the Inter-Benchmark Protocol on New Species (Turbot and Sea bass). ICES C. 2012/ACOM45.
  set b_g 2.969; Length-weight parameter
  set eggs_per_bass 375000 ;bass book
  set GL 0.02485 ;young lifwe stage mort


  set H  8.17E-01 ;0.5
  set Am 5.61E-04 ; 6.60E-04
  set AE 2.18E-03; 0.01
  set PM 6.15E-02 ;0.16686103 ;0.03686103 ;; 0.03712
  set I  6.70E+13; 5E+13


  setup_turtles



 end




to setup-ICES                                                                   ; loads ICES rectangles shapefile
  ask patches
  [
    set row ceiling (pycor / 2) + 24                                            ; the ICES row in which each patch is situated is set. This begins at 26 (lowest lat in model domain)
                                                                                ; and each ICES row spans 6 patches
    ifelse pxcor < 18
    [set column word "E" (ceiling ((pxcor + 1) / 2)) ]                          ; the same is done for each column, but they contain a letter prefix, which is added using the word
    [set column word "F" (ceiling (((pxcor + 1) - 18) / 2))]

    set rectangle word row column
    set plabel rectangle                                                        ; the rectangle name is given by concatenating the row and the column, eg 38D8
  ]
end

to setup_turtles ;so should we be using 1985 stuff here?

  file-open "data/2020/icesagedistribution.txt"

  while [not file-at-end?]
  [
    let ayear file-read
    if ayear = year
    [
      let acounter 0
      while [acounter <= 30]
      [
        let j file-read
        ifelse acounter < 4
        [
          ask n-of ncohort patches with [nursery-patches = true]
          [
            sprout 1[
              set breed juv-bass
              set size 0.4
              set shape "fish"
              set age acounter + (222 + random 93) / 365
              set cohort acounter
              set color yellow
              set number j / 10
              set L linf * (1 - exp(- k  * (age - t0)))
            ]
          ]
        ]
        [
          ifelse acounter < 6
          [
            ask n-of ncohort patches with [coast-patches = true]
            [
              sprout 1[
                set breed juv-bass
                set size 0.6
                set shape "fish"
                set age acounter + (222 + random 93) / 365
                set cohort acounter
                set color yellow
                set number j / 10
                set L linf * (1 - exp(- k  * (age - t0)))
              ]
            ]
          ]
          [
            ask n-of ncohort patches with [viie = true or viifg = true]
            [
              sprout 1[
                set breed mat-bass
                set size 0.8
                set shape "fish"
                set age acounter + (222 + random 93) / 365
                set cohort acounter
                set color yellow
                set number j / 10
                set L linf * (1 - exp(- k  * (age - t0)))
                set spawn-trigger 1
              ]
            ]
          ]
        ]
        set acounter acounter + 1
      ]
    ]
  ]

  file-close


 ask turtles [  ;JW RB



  set L linf * (1 - exp(- k  * (age - t0)))
  set std-mass a_g * (L ^ (b_g))
  set structural-mass std-mass
  set energy-reserve-max ((structural-mass * 0.01) * El)
  set energy-reserve-max-mass (energy-reserve-max / El)
  set energy-reserve (energy-reserve-max * 0.5)
  set total-mass (structural-mass + (energy-reserve / El))
  set maintenance-energy (energy-reserve * 0.1    )
  set potential-fecundity (std-mass ) * (eggs_per_bass)                      ; maximum number of eggs that can be produced, potential fecundity, depends on body length , changed for bass JW pickett and pawson fecundity estimate.
  set max-R (potential-fecundity * egg_mass * (Ef + Fs))
  calculate-speed
  calculate-r



 ]
   ; site fidelity
   ask mat-bass
   [
     let decider random 6
     if decider = 0 [ set ivbt true  ]
     if decider = 1 [ set ivct true ]
     if decider = 2 [ set viidt true ]
     if decider = 3 [ set viiet true ]
     if decider = 4 [ set viifgt true ]
     if decider = 5 [ set viiat true ]
   ]

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; RUN ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



to go

ask turtles [
  if ( number < 0 ) ;should stop error
    [set number  0.01]
 ]


 ; if  (year = 2004 and ticks < 1) [ reset-timer ]
  if not any? turtles [stop]

  tick


  ; PATCH UPDATES
 ifelse rolling_tick < 6935

  [
     if ticks = 1
  [set directory (word "data/temperature/fourpatch/2004")]    ; change directory at the beginning of the year
  ]


  [
  if ticks = 1
  [ set directory (word "data/temperature/fourpatch/" year) ]   ; change directory at the beginning of the year
  ]


ask turtles

[ sea_check ]



update-patches

  if ( ticks = 60 ) ;arbitary date 17th march
  [calc_sst_mean]



;;;;;;;;;;;;;;;;;;;;;;;ssb tsb nums;;;;;;;;;;;

 set num (sum [number] of turtles with [(breed = juv-bass or breed = mat-bass)])  ;to macth with ICES 1000s of in output graphs
 set ssb (sum [total-mass * number] of mat-bass) / 1000    ; ssb (tonnes)
 set tsb (sum [total-mass * number] of turtles) / 1000     ; tsb (tonnes)
 set rec sum [number] of turtles with [(age < 1) and (breed = juv-bass)]


 if ticks = 1 [set annual_num num]   ;change these when correspond with ICES
 if ticks = 1 [set annual_tsb tsb]
 if ticks = 1 [set annual_rec rec]
 if ticks = 1 [set annual_ssb ssb]




if year >= 2004

[

 if ticks = 1 [

 set ABC_SSB  lput annual_ssb ABC_SSB
 set ABC_Rec  lput annual_rec ABC_Rec
 set ABC_Num  lput annual_num ABC_Num

 ]

]



if year >= 2004

[

 if ticks = 1 [

calc-length-distribution
calc-mass-distribution
calc-age-distribution
calc-total-mass-distribution

 ]

]



if year >= 2004

[

 if ticks = 1 [

set ABC_M0  lput M0 ABC_M0
set ABC_M1  lput M1 ABC_M1
set ABC_M2  lput M2 ABC_M2
set ABC_M3  lput M3 ABC_M3
set ABC_M4  lput M4 ABC_M4
set ABC_M5  lput M5 ABC_M5
set ABC_M6  lput M6 ABC_M6
set ABC_M7  lput M7 ABC_M7
set ABC_M8  lput M8 ABC_M8
set ABC_M9  lput M9 ABC_M9
set ABC_M10 lput M10 ABC_M10
set ABC_M11 lput M11 ABC_M11
set ABC_M12 lput M12 ABC_M12
set ABC_M13 lput M13 ABC_M13
set ABC_M14 lput M14 ABC_M14
set ABC_M15 lput M15 ABC_M15
set ABC_M16 lput M16 ABC_M16
set ABC_M17 lput M17 ABC_M17
set ABC_M18 lput M18 ABC_M18
set ABC_M19 lput M19 ABC_M19
set ABC_M20 lput M20 ABC_M20
set ABC_M21 lput M21 ABC_M21
set ABC_M22 lput M22 ABC_M22
set ABC_M23 lput M23 ABC_M23
set ABC_M24 lput M24 ABC_M24
set ABC_M25 lput M25 ABC_M25
set ABC_M26 lput M26 ABC_M26
set ABC_M27 lput M27 ABC_M27
set ABC_M28 lput M28 ABC_M28
set ABC_M29 lput M29 ABC_M29
 ]

]


if year >= 2004

[

if ticks = 1 [

set ABC_N0  lput N0 ABC_N0
set ABC_N1  lput N1 ABC_N1
set ABC_N2  lput N2 ABC_N2
set ABC_N3  lput N3 ABC_N3
set ABC_N4  lput N4 ABC_N4
set ABC_N5  lput N5 ABC_N5
set ABC_N6  lput N6 ABC_N6
set ABC_N7  lput N7 ABC_N7
set ABC_N8  lput N8 ABC_N8
set ABC_N9  lput N9 ABC_N9
set ABC_N10 lput N10 ABC_N10
set ABC_N11 lput N11 ABC_N11
set ABC_N12 lput N12 ABC_N12
set ABC_N13 lput N13 ABC_N13
set ABC_N14 lput N14 ABC_N14
set ABC_N15 lput N15 ABC_N15
set ABC_N16 lput N16 ABC_N16
set ABC_N17 lput N17 ABC_N17
set ABC_N18 lput N18 ABC_N18
set ABC_N19 lput N19 ABC_N19
set ABC_N20 lput N20 ABC_N20
set ABC_N21 lput N21 ABC_N21
set ABC_N22 lput N22 ABC_N22
set ABC_N23 lput N23 ABC_N23
set ABC_N24 lput N24 ABC_N24
set ABC_N25 lput N25 ABC_N25
set ABC_N26 lput N26 ABC_N26
set ABC_N27 lput N27 ABC_N27
set ABC_N28 lput N28 ABC_N28
set ABC_N29 lput N29 ABC_N29


 ]

]



;;;;;;;;;;;;;;;GROWTH/LENGTH-BASED PROCESSES;;;;;;

calculate-catch

 ask turtles

[ natural-mortality ]


 ask turtles [
if (breed != eggs) or (breed != ys-larvae)

[
  calculate-ingestion
  calculate-assimilation
  calc-maintenance
 ]


if (breed != eggs)
[
  calc-growth
  calc-total-mass
  calculate-speed
  calculate-r
]

  transform

  ]


;;;;;;;;;;;;;;;;;;;fishing;;;;;;;;;;;;;;;;;;;;;;;

    if (ticks = 1)
    [ ask patches [ set catch 0 ] ]


  recreational-inshore-F
  commercial-inshore-F
  commercial-offshore-F
  fishing-mortality


icesssbtsb                                       ; load assessment ssb and tsb for comparison
plot-age-distribution                            ; plot histogram of model age distribution against assessment distribution (plots done in 1000s)



;  if (ticks = 1)
;[
;  set num (sum [number] of turtles) / 1000         ; numbers-at-age (000s)
;  set ssb (sum [Wrep * number] of mat-bass) / 1000    ; ssb (tonnes)
;  set tsb (sum [Wrep * number] of turtles) / 1000     ; tsb (tonnes)
;]



;;;;;;;;;;;;;;;;;;;;;movement;;;;;;;;;;;;;;;;;;;;;;;;;;

  ; MOVEMENT

 ask turtles

[ if (breed = juv-bass) or (breed = mat-bass)

  [avoid-corwall]
]


local-movement


if ( ticks >= 274 ) or ( ticks < 152 )


[spawn-migration]

feeding-migration



   drift_eggs

   drift_ys-larvae

   drift_larvae

;;;;;;;;;;;;;;;;;;;;;AGE;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ask turtles [ set age age + (1 / 365) ]
  transform

;;;;;;;;;;;;;;;;;;spawn;;;;;;;;;;;;;;;;;;;


 ifelse rolling_tick < 6935
;ifelse year <= 2003

[

 ; CONNECTIVITY AND SPAWNING
  ; Number and location of spawning super-individuals in spawning months?
  if (ticks = 46) [
    set Nfeb count mat-bass with [spawn-count > 0 and (NS = true or IS = true or C = true or CS = true)]
    set Locfeb (list (count mat-bass with [spawn-count > 0 and NS = true]) (count mat-bass with [spawn-count > 0 and IS = true]) (count mat-bass with [spawn-count > 0 and C = true]) (count mat-bass with [spawn-count > 0 and CS = true]))
  ]
  if (ticks = 75) [
    set Nmar count mat-bass with [spawn-count > 0 and (NS = true or IS = true or C = true or CS = true)]
    set Locmar (list (count mat-bass with [spawn-count > 0 and NS = true]) (count mat-bass with [spawn-count > 0 and IS = true]) (count mat-bass with [spawn-count > 0 and C = true]) (count mat-bass with [spawn-count > 0 and CS = true]))
  ]
  if (ticks = 106) [
    set Napr count mat-bass with [spawn-count > 0 and (NS = true or IS = true or C = true or CS = true)]
    set Locapr (list (count mat-bass with [spawn-count > 0 and NS = true]) (count mat-bass with [spawn-count > 0 and IS = true]) (count mat-bass with [spawn-count > 0 and C = true]) (count mat-bass with [spawn-count > 0 and CS = true]))
  ]
  if (ticks = 136) [
    set Nmay count mat-bass with [spawn-count > 0 and (NS = true or IS = true or C = true or CS = true)]
    set Locmay (list (count mat-bass with [spawn-count > 0 and NS = true]) (count mat-bass with [spawn-count > 0 and IS = true]) (count mat-bass with [spawn-count > 0 and C = true]) (count mat-bass with [spawn-count > 0 and CS = true]))

    set nspawn (list (round (Nfeb / (Nfeb + Nmar + Napr + Nmay) * ncohort)) (round (Nmar / (Nfeb + Nmar + Napr + Nmay) * ncohort)) (round (Napr / (Nfeb + Nmar + Napr + Nmay) * ncohort)) (round (Nmay / (Nfeb + Nmar + Napr + Nmay) * ncohort)) )
    if ncohort - sum nspawn != 0
    [
      let x ncohort - sum nspawn
      let index position max nspawn nspawn
      set nspawn replace-item index nspawn (item index nspawn + x)
    ]
  ]

  if ( ticks = 152 or ticks = 182 or ticks = 213 or ticks = 244)
  [spin_spawn]


]


[

ask mat-bass [

  if ( ticks = 60 ) ;arbitary date 17th march
  [
  spawn
   calc_sst_mean

  ]
]






;;;;;;;;;;;;eggs dvelop;;;;;;;;;;;;;;;;;;

  ask eggs
  [calc-egg-development]


]
;;;;;;;;;;;;;;;;;;next year;;;;;;;;;;;;




set rolling_tick rolling_tick + 1



 if ticks >= 365[

  ask turtles [ set cohort cohort + 1 ]
   reset-ticks
   set year year + 1
    if year = 2015 [
     stop   ]
  ]

;]







end






to spin-up

repeat 6935 [go]

end

to go-ABC

repeat 3655 [go]

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to sea_check

 if (sea != true) [

 let target-patch min-one-of patches in-radius 2 with [sea = true] [distance myself]      ; check on sea patch with SST and CHL present if so carry on if not move turtle back to closest sea patch and carry on. JW
 move-to target-patch


 ]




end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; PATCH UPDATES ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to update-patches



;if rolling_ticks > 20*365 do this only for 2005 so only up to 46 in the list then start at on again. once reached rolling equivalent of 20 years, then can update beyond this and go on with the model?

 ifelse rolling_tick < 6935
[

set year_spin  2004

  if (ticks mod 8 = 0) and (ticks != 1)
  [
    set directory "data/rs_data_bass"
    set raster (word directory "/chl_" ((ticks / 8) + ((year_spin - 2004) * 46)) ".asc")
    set phyto-data gis:load-dataset  raster
    gis:apply-raster phyto-data chl

    set directory "data/rs_data_bass"
    set rasterSST (word directory "/sst_" ((ticks / 8) + ((year_spin - 2004) * 46)) ".asc")
    set sst-data gis:load-dataset  rasterSST
    gis:apply-raster sst-data sst

     ask patches with [(coast-patches = true)]
    [
      ifelse sst > 0
      []
      [set sst [sst] of min-one-of patches with [sst > 0] [distance myself]]

      ifelse chl > 0
      []
      [set chl mean [chl] of neighbors with [chl > 0]]
    ]
  ]


  ask patches with [offshore-patches = true]
  [
    set pcolor scale-color blue SST 5 18   ; recolour offshore patches according to new sst value

    ifelse ( (sst >= 9 and sst <= 15) and (ticks > 31 and ticks < 152) and (pycor < 24))      ; highlight spawning grounds, dependant on sst, location and time of year (Feb-May)
    [
      set spawn-patches true
      if (pxcor mod 2 = 0 and pycor mod 2 = 1)
      [ set plabel "S"
        set plabel-color yellow
      ]
    ]
    [
      set spawn-patches 0                  ; remove spawning label for patches no longer meeting required sst
      set plabel ""
    ]
 ]


]



; once past 2004

[

  if (ticks mod 8 = 0) and (ticks != 1)
  [
    set directory "data/rs_data_bass"
    set raster (word directory "/chl_" ((ticks / 8) + ((year - 2004) * 46)) ".asc")
    set phyto-data gis:load-dataset  raster
    gis:apply-raster phyto-data chl

    set directory "data/rs_data_bass"
    set rasterSST (word directory "/sst_" ((ticks / 8) + ((year - 2004) * 46)) ".asc")
    set sst-data gis:load-dataset  rasterSST
    gis:apply-raster sst-data sst

     ask patches with [(coast-patches = true)]
    [
      ifelse sst > 0
      []
      [set sst [sst] of min-one-of patches with [sst > 0] [distance myself]]

      ifelse chl > 0
      []
      [set chl mean [chl] of neighbors with [chl > 0]]
    ]
  ]




  ask patches with [offshore-patches = true]
  [
    set pcolor scale-color blue SST 5 18   ; recolour offshore patches according to new sst value

    ifelse ( (sst >= 9 and sst <= 15) and (ticks > 31 and ticks < 152) and (pycor < 24))      ; highlight spawning grounds, dependant on sst, location and time of year (Feb-May)
    [
      set spawn-patches true
      if (pxcor mod 2 = 0 and pycor mod 2 = 1)
      [ set plabel "S"
        set plabel-color yellow
      ]
    ]
    [
      set spawn-patches 0                  ; remove spawning label for patches no longer meeting required sst
      set plabel ""
    ]
 ]


]


end














to calculate-ingestion ;JW RB

ifelse (breed != eggs) or  (breed != YS-larvae)

[

 if (any? other turtles-here)
 [
     set num_near_me (number + sum [number] of other turtles-here)
     set mass_near_me (number * total-mass + sum [number * total-mass ] of other turtles-here)
 ]

]

[  ]



ifelse breed != larvae
[
   ifelse SST > 5

  [
    set func-response (chl / (h + chl))
 ifelse (any? other turtles-here)
 [set ingestion-rate (Cmax * exp((- Ea / Boltz) * ((1 / (SST + 273.15)) - (1 / Tref))) * func-response * (total-mass ^ (2 / 3))) * (i * (1 / ((mass_near_me ^ 2 / 3) )))]
 [set ingestion-rate (Cmax * exp((- Ea / Boltz) * ((1 / (SST + 273.15)) - (1 / Tref))) * func-response * (total-mass ^ (2 / 3))) * (i * (1 / ((total-mass ^ 2 / 3) * number )))]
  ]

   [set ingestion-rate 0 ]
]


[
set func-response (chl / (h + chl))
 ifelse (any? other turtles-here)
 [set ingestion-rate (Cmax * exp((- Ea / Boltz) * ((1 / (SST + 273.15)) - (1 / Tref))) * func-response * (total-mass ^ (2 / 3))) * (i * (1 / ((mass_near_me ^ 2 / 3) )))]
 [set ingestion-rate (Cmax * exp((- Ea / Boltz) * ((1 / (SST + 273.15)) - (1 / Tref))) * func-response * (total-mass ^ (2 / 3))) * (i * (1 / ((total-mass ^ 2 / 3) * number )))]
]

end



to calculate-assimilation  ;JW RB

  set assimilation-energy ((ingestion-rate * ep ) * AE) ;absorbed energy is on slider

end


to calc-maintenance ;JW RB

set MR (A0 * (total-mass ^ (0.75)) * 2) * exp((- Ea / boltz) * ((1 / (SST + 273.15)) - (1 / Tref))) ; smr x2 as in peters

ifelse assimilation-energy > MR                                                          ; if enough energy is available in "assimilated" store to pay maintenance costs
 [set assimilation-energy assimilation-energy - MR   ]                                    ; subtract these costs from the store

 [                                                                                       ; if not,
  set energy-reserve energy-reserve + assimilation-energy                          ; add the assimilated energy to the reserves
   set energy-reserve energy-reserve - MR                                                 ; then subtract the costs of maintenance
    set assimilation-energy  0                                                            ; and set assimilation-energy 0 ; chnaged to get rid as reserves will still be there?
  ]
                                                                                          ; the individual will die (in the procedure "calc-starvation") if this results in empty energy reserves
end



to calc-growth                                                                            ; if t < 240 days then grow according to Gompertz, and von bertalanffy if not

ifelse age > 0.20 ; first 70 days

 [if (breed = mat-bass) or  (breed = juv-bass) ;and  (SST > 9)

   [set max-growth-rate (linf - L) * (1 - exp(- k / 365)) * exp((- Ea / boltz) * ((1 / (SST + 273.15)) - (1 / Tref)))]
      ifelse L < Linf
      [grow]
      [stop]
  ]

  [ set max-growth-rate   GL   * exp((- Ea / boltz) * ((1 / (SST + 273.15)) - (1 / Tref)))  ; first 70 days grown as straight line
    grow
  ]
end


to grow  ;JW  RB

  let possible-L (L + max-growth-rate)
  set growth-costs ((((a_g * (possible-L ^ b_g)) - structural-mass)) * (Fs + Ef))


    if (breed = mat-bass) or ( breed = juv-bass ) [



    ifelse (assimilation-energy * 0.5) >= growth-costs                                       ; juveniles and adults allocate energy equally to growth in length and to fat reserves
   [
     set L L + max-growth-rate
     set structural-mass a_g * (L ^ b_g)                                            ; the new structural mass is calculated. Their total mass is calculated later, after they have stored lipid
     set assimilation-energy assimilation-energy - growth-costs
     calc-storage                                                                          ; a procedure "calc-storage" is called that converts remaining energy to lipid stores

    ]
    [

      set growth-rate (max-growth-rate / growth-costs) * (assimilation-energy * 0.5)         ; sub-optimal growth rate if energy is insufficient
      set L L + growth-rate
      set structural-mass a_g * (L ^ b_g)
      set assimilation-energy assimilation-energy * 0.5
      calc-storage
   ]

    ]





 if (breed = larvae)
[
    ifelse assimilation-energy >= growth-costs                                               ; if eneough energy is assimilated to cover max growth individuals grow maximally
    [
      set L L + max-growth-rate                                                             ; length is updated
      set structural-mass a_g * (L ^ b_g)                                           ; new structural mass is calculated
      set total-mass structural-mass                                                        ; because larvae do not store lipid, their structural mass is equal to their total mass
      set assimilation-energy assimilation-energy - growth-costs                              ; subtract the costs of growth from the assimilated energy
    ]
    [
      set growth-rate (max-growth-rate / growth-costs) * assimilation-energy                 ; growth rate is adjusted if there is not enough energy
      set L L + growth-rate
      set structural-mass a_g * (L ^ (b_g))   ;mass in KG
      set assimilation-energy 0
    ]

     ]



  if breed = YS-larvae                                                                      ; yolk sac larvae are nourished by the yolk so grow maximally
  [
    set growth-rate max-growth-rate
    set L L + max-growth-rate
    set structural-mass a_g * (L ^ (b_g))
    set total-mass structural-mass
  ]



end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to calc-storage ; JW RB

  set std-mass a_g * (L ^ (b_g))   ;; changed these to bass fishbase JW

  set structural-mass std-mass
                                                                                           ; juveniles and adults store any remaining energy as lipid
  set energy-reserve-max ((structural-mass * 0.01) * El)                                   ; individual can store no more lipid than the equivalent of 78% of their structural mass (see TRACE for why)


set energy-reserve-max-mass (energy-reserve-max / El)

  if assimilation-energy > 0
  [
    set energy-reserve energy-reserve + (assimilation-energy * (El / (El + Ls)))           ; costs of synthesis are accounted for in storing lipid

  ]
  if energy-reserve > energy-reserve-max
  [
    set energy-reserve energy-reserve-max

  ]


end


;;;; total mass (structural + lipid + gonad) is calculated ;;;;

to calc-total-mass ;JW RB

 ifelse energy-reserve > 0
 [ set total-mass structural-mass + (energy-reserve / El) + gonad-mass ]         ; total mass is the sum of structural, fat and gonad mass
 [set total-mass structural-mass + gonad-mass]
  set FK 100 * (total-mass / (L ^ (3)))                                          ; Fulton's condition factor


end



to calc-reproduction ;JW RB
                                                                                   ; the energy costs of reproduction are calculated changed from robs calculated from length for bass calc from mass
    set potential-fecundity (std-mass ) * (eggs_per_bass)                      ; maximum number of eggs that can be produced, potential fecundity, depends on body length , changed for bass JW pickett and pawson fecundity estimate.
    set max-R (potential-fecundity * (egg_mass * (Ef + Fs)))                         ; the energy cost of producing the maximum number of eggs, equation 8
    set maintenance-energy energy-reserve * 0.1                                   ; at the beginning of spawning 10% of an individuals energy is set aside to cover maintenance costs while spawning; it cannot be allocated to reproduction ; only small amount


end

to spawn                                                                            ; JW Rb, changed from nicolas "old_spawn"

ask mat-bass [

      calc-reproduction
      ifelse (energy-reserve - maintenance-energy) >= max-R     ;if there is enough energy after the maintenance energy (see "calc-reproduction") is set aside, the maximum amount is sent to the reproduction store R

      [
        set energy-reserve energy-reserve - max-R                          ; the energy needed to produce the max number of eggs is taken from the energy reserve
        set gonad-mass  (max-R / (Ef + Fs))    / 1000                            ;gonad mass is calculated from the energy that went in to producing the eggs
        set realised-fecundity potential-fecundity                         ; realised fecundity is set to max potential fecundity

     ]
 [

        set ER energy-reserve - maintenance-energy ;energy reserve is what is left after maintenence is taken out
        set gonad-mass (ER) / (Ef + Fs) / 1000
        set realised-fecundity (ER / (max-R)) * potential-fecundity
      ]
       ]



if who = [who] of max-one-of mat-bass [who]   ; so this needs to be 10 bass split into the 4 dffernt sea areas ;set affitity to chanfge
 [

    ask n-of 1 mat-bass
     [set ivbt  false
      set IVct  false
      set VIIdt false
      set VIIet false
      set VIIfgt false
      set VIIat false

      set ivbt true ;Nsea
      set color red
      deposit-eggs]


ask n-of 1 mat-bass
     [set ivbt  false
      set IVct  false
      set VIIdt false
      set VIIet false
      set VIIfgt false
      set VIIat false

      set ivct true ;Nsea
      set color red
      deposit-eggs]



ask n-of 1 mat-bass
     [set ivbt  false
      set IVct  false
      set VIIdt false
      set VIIet false
      set VIIfgt false
      set VIIat false

      set viidt true ;Channel
      set color red
      deposit-eggs]






ask n-of 1 mat-bass
     [set ivbt  false
      set IVct  false
      set VIIdt false
      set VIIet false
      set VIIfgt false
      set VIIat false

      set viidt true ;Channel
      set color red
      deposit-eggs]


ask n-of 1 mat-bass
     [set ivbt  false
      set IVct  false
      set VIIdt false
      set VIIet false
      set VIIfgt false
      set VIIat false

      set viidt true ;Channel
      set color red
      deposit-eggs]



ask n-of 1 mat-bass
     [set ivbt  false
      set IVct  false
      set VIIdt false
      set VIIet false
      set VIIfgt false
      set VIIat false

      set viiet true ;Channel
      set color red
      deposit-eggs]




ask n-of 1 mat-bass
     [set ivbt  false
      set IVct  false
      set VIIdt false
      set VIIet false
      set VIIfgt false
      set VIIat false

      set viiet true ;Channel
      set color red
      deposit-eggs]



ask n-of 1 mat-bass
     [set ivbt  false
      set IVct  false
      set VIIdt false
      set VIIet false
      set VIIfgt false
      set VIIat false

      set viifgt true ;Celtic and Irish
      set color red
      deposit-eggs]


ask n-of 1 mat-bass
     [set ivbt  false
      set IVct  false
      set VIIdt false
      set VIIet false
      set VIIfgt false
      set VIIat false

      set viifgt true ;Celtic and Irish
      set color red
      deposit-eggs]




ask n-of 1 mat-bass
     [set ivbt  false
      set IVct  false
      set VIIdt false
      set VIIet false
      set VIIfgt false
      set VIIat false

      set viiat true ;Celtic and Irish
      set color red
      deposit-eggs]





 ]





end




to deposit-eggs       ; egg super-individuals are "hatched" and their variables are set. Most must e set to 0 otherwise they are inherited from the parent

  set gonad-mass 0

  hatch-eggs 1

  [

    set color white
    set shape "dot"
    set size 1
    set std-mass 0.001
    set L 0.1
    set energy-reserve 0
    set structural-mass std-mass
    set total-mass std-mass
    set development 1


   if (ivbt =  True  or ivct =  True)
     [
      set number (sum [(realised-fecundity * number)] of mat-bass / ncohort) * 2  ;north sea
      ]

   if (viidt = True  or viiet = True)
    [
    set number (sum [(realised-fecundity * number)] of mat-bass / ncohort) * 5    ;channel
    ]

   if (VIIfgt = True  or viia = True)
    [
    set number (sum [(realised-fecundity * number)] of mat-bass / ncohort) * 3    ;west coast irish and celtic sea
    ]



    set MR 0
    set ingestion-rate 0
    set max-r 0
    set R 0
    set growth-rate 0
    set max-growth-rate 0
    set growth-costs 0
    set func-response 0
    set fk 0
    set age 0
    set embryo-duration 5 ; round ((2.08 * (10 ^ (-11))) * exp(Ea / (Boltz * (SST + 273.15))))  ; embryo duration is set at 5 days for simplicity

      ]
end


;gonad mass is the same as weight of eggs, these need to be depostited and then hang around for egg durations;
;then turn into larvea and be exposed to the mortalirty rate as a function of size before becoming juviniles.


;;;; egg development ;;;;

to calc-egg-development
  set development development + 1                                       ; eggs get 1 day closer to hatching each time-step
end






;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  AGE ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


to transform
  ask eggs
  [
    if development = embryo-duration                    ; when eggs have developed they hatch into larvae
    [
      set breed YS-larvae
      set color black
      set size 0.3
      set shape "fish"
      set L 0.3                                         ; Villamore et al. (2004, cm)
      set std-mass 0.001                                ; g (Sibly et al. 2015)
      set energy-reserve-max 0
      set energy-reserve energy-reserve-max ;
      set larval-production larval-production + number



    ]
  ]

  ask YS-larvae
  [
    ;if L >= 0.61                                        ; size threshold (Sette)

    if L >= 0.35    ;1.5- 5.5mm beraud 2018
    [
      set breed larvae
      set color black
      set size 0.3
      set shape "fish"

      set std-mass a_g * L ^ (b_g) ; same as adults

      set structural-mass std-mass

      set total-mass structural-mass

      set energy-reserve-max ((structural-mass * 0.01) * El)

      set energy-reserve energy-reserve-max

    ]
  ]

  ask larvae
  [
    ;if (L >= 3)
    if (L >= 1.425)

    [
      set breed juv-bass

      set color yellow
      set size 0.3
      set shape "fish"
      set energy-reserve-max ((structural-mass * 0.01) * El)
      set energy-reserve energy-reserve-max * 0.5


    ]
  ]





 ask juv-bass
  [if  L >= 42

   [


ifelse coast-patches = True

[

    set breed mat-bass
    set size 1
    set shape "fish"

    ; newly recruited mat bass allocated current area to detemine migratory behaviour

       if IVb = true [
      set IVbt  true
      set IVct  false
      set VIIdt false
      set VIIet false
      set VIIfgt false
      set VIIat false]



    if IVc = true [
      set IVct true
      set VIIdt false
      set VIIet false
      set VIIfgt false
      set VIIat false
      set IVbt  false]

    if VIId = true [
      set IVct false
      set VIIdt true
      set VIIet false
      set VIIfgt false
      set VIIat false
      set IVbt  false]

    if VIIe = true [
      set IVct false
      set VIIdt false
      set VIIet true
      set VIIfgt false
      set VIIat false
      set IVbt  false]

    if VIIfg = true [
      set IVct false
      set VIIdt false
      set VIIet false
      set VIIfgt true
      set VIIat false
      set IVbt  false]

    if VIIa = true [
      set IVct false
      set VIIdt false
      set VIIet false
      set VIIfgt false
      set VIIat true
      set IVbt  false]

]



[
  ;check in radius for coast patch
  ;move towards it
  ;repeat above

    let coasttarget-patch min-one-of patches in-radius 5 with [coast-patches = true] [distance myself]      ; check on sea patch with SST and CHL present if so carry on if not move turtle back to closest sea patch and carry on. JW
    move-to coasttarget-patch

    set breed mat-bass
    set size 1
    set shape "fish"

    ; newly recruited mat bass allocated current area to detemine migratory behaviour

      if IVb = true [
      set IVbt  true
      set IVct  false
      set VIIdt false
      set VIIet false
      set VIIfgt false
      set VIIat false]

      if IVc = true [
      set IVct true
      set VIIdt false
      set VIIet false
      set VIIfgt false
      set VIIat false
      set IVbt  false]

      if VIId = true [
      set IVct false
      set VIIdt true
      set VIIet false
      set VIIfgt false
      set VIIat false
      set IVbt  false]

    if VIIe = true [
      set IVct false
      set VIIdt false
      set VIIet true
      set VIIfgt false
      set VIIat false
      set IVbt  false]

    if VIIfg = true [
      set IVct false
      set VIIdt false
      set VIIet false
      set VIIfgt true
      set VIIat false
      set IVbt  false]

    if VIIa = true [
      set IVct false
      set VIIdt false
      set VIIet false
      set VIIfgt false
      set VIIat true
      set IVbt  false]

  ]


  ]]

end




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; GROWTH/LENGTH-BASED PROCESSES ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


to calculate-speed ; speed is not kmph is patches per day (12hrs) divide speed by 12 in need kmph which is what rob has used.

    set speed  12 * (10 ^ ( -0.828 + 0.6196 * ( log L 10 ) + 0.3478 * ( log A 10 ))) / 30  *  ((exp ( - (EaS / boltz) * ((1 / (sst + 273)) - (1 / refTS)))))
end


to calculate-r        ; ensures fish never travel more then 0.25 patches in one go, reducing overlap with land

    set r 1
    while [( speed / r ) > 0.25 ][set r r + 1]

end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; MORTALITY ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


to natural-mortality


ifelse (breed = mat-bass) or (breed = juv-bass)                       ; larval mortality is constant
    [
    set number number * exp ( - Am )                ; exponential decay
    if age > 30  [die]                      ; bass above age 30 die
    ]

    [set number number * exp ( - Pm )]


end






;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; MOVEMENT ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to spawn-migration

  ask mat-bass
  [
    if (spawn-trigger = 0)
    [
      if (ticks < 91 or ticks >= 274) and (any? neighbors with [ spawn-patches = true ] or [ sst ] of patch-here < 9) ; No onset after March
      [
        set spawn-trigger 1
        set spawn-count 0
      ]
    ]
    if (spawn-trigger = 1)
    [
      repeat r
      [
        ifelse (spawn-patches = true) ; 1. Am I in a spawning area?
        [
          rt random 360 fd (speed / r)
        ]
        [
          ifelse (any? neighbors with [spawn-patches = true]) ; 2. Am I neighboring a spawning area?
          [
            face one-of neighbors with [ spawn-patches = true ]
            fd (speed / r)
          ]
          [
            ifelse (coast-patches = true) ; 3. Move along coast
            [
              if (ivb = true)
              [
                ifelse any? neighbors with [(coast-patches = true) and (pycor < [ ycor ] of myself)]
                [
                  ifelse ((pxcor = 33 or pxcor = 34) and pycor = 23) ; awkward patches!
                  [
                    face one-of neighbors with [(coast-patches = true) and (pxcor < [ xcor ] of myself)]
                    fd (speed / r)
                  ]
                  [
                    face one-of neighbors with [(coast-patches = true) and (pycor < [ ycor ] of myself)]
                    fd (speed / r)
                  ]
                ]
                [
                  ifelse (pxcor > 29) or (pycor > 33)
                  [
                    face one-of neighbors with [(coast-patches = true) and (pxcor < [ xcor ] of myself)]
                    fd (speed / r)
                  ]
                  [
                    face one-of neighbors with [(coast-patches = true) and (pxcor > [ xcor ] of myself)]
                    fd (speed / r)
                  ]
                ]
              ] ; end of ivb
              if (ivc = true)
              [
                ifelse (any? neighbors with [(coast-patches = true) and (pycor < [ round(ycor) ] of myself)])
                [
                  ifelse (pxcor = 27 or pxcor = 28) and (pycor = 19 or pycor = 20) ; Awkward patches!
                  [
                    ifelse (any? neighbors with [(coast-patches = true) and (pxcor < [ xcor ] of myself)])
                    [
                      face one-of neighbors with [(coast-patches = true) and (pxcor < [ xcor ] of myself)]
                      fd (speed / r)
                    ]
                    [
                      face one-of neighbors with [(coast-patches = true) and (pycor > [ ycor ] of myself)]
                      fd (speed / r)
                    ]
                  ]
                  [
                    ifelse (pxcor = 19 and pycor = 20) ; Awkward patch!
                    [
                      face one-of neighbors with [(coast-patches = true) and (pxcor > [ xcor ] of myself)]
                      fd (speed / r)
                    ]
                    [
                      face one-of neighbors with [(coast-patches = true) and (pycor < [ ycor ] of myself)]
                      fd (speed / r)
                    ]

                  ]
                ]
                [
                  ifelse (pxcor < 22 and pycor > 13)
                  [
                    face one-of neighbors with [(coast-patches = true) and (pxcor > [ xcor ] of myself)]
                    fd (speed / r)
                  ]
                  [
                    ifelse (pxcor = 28 and pycor = 18)
                    [
                      face one-of neighbors with [(coast-patches = true) and (pycor > [ ycor ] of myself)]
                      fd (speed / r)
                    ]
                    [
                      face one-of neighbors with [(coast-patches = true) and (pxcor < [ xcor ] of myself)]
                      fd (speed / r)
                    ]
                  ]
                ]
              ] ; end of ivc
              if (viid = true)
              [
                ifelse (any? neighbors with [(coast-patches = true) and (pxcor < [ xcor ] of myself)])
                [
                  face one-of neighbors with [(coast-patches = true) and (pxcor < [ xcor ] of myself)]
                  fd (speed / r)
                ]
                [
                  ifelse (any? neighbors with [(coast-patches = true) and (pycor < [ ycor ] of myself)])
                  [
                    face one-of neighbors with [(coast-patches = true) and (pycor < [ ycor ] of myself)]
                    fd (speed / r)
                  ]
                  [
                    face one-of neighbors with [(coast-patches = true) and (pycor > [ ycor ] of myself)]
                    fd (speed / r)
                  ]
                ]
              ] ; end of viid
              if (viie = true) ; Move randomly in 7e or head to Cornwall? Currently random
              [
;                    ifelse (any? neighbors with [(coast-patches = true) and pxcor < [ xcor ] of myself])
;                    [
;                      face one-of neighbors with [(coast-patches = true) and (pxcor < [ xcor ] of myself)]
;                      fd (speed / r)
;                    ]
;                    [
;                      face one-of neighbors with [(coast-patches = true) and (pycor < [ ycor ] of myself)]
;                      fd (speed / r)
;                    ]
;                  ]
;                ]
;                [
                  rt random 360 fd (speed / r)
;                ]
              ] ; end of viie
              if (viifg = true)
              [
                ifelse pxcor < 4
                [
                  face one-of neighbors with [coast-patches = true]
                  fd (speed / r)
                ]
                [
                  ifelse any? neighbors with [(coast-patches = true) and (pycor < [ ycor ] of myself)]
                  [
                    face one-of neighbors with [(coast-patches = true) and (pycor < [ ycor ] of myself)]
                    fd (speed / r)
                  ]
                  [
                    ifelse pycor > 13
                    [
                      face one-of neighbors with [(coast-patches = true) and (pxcor > [ xcor ] of myself)]
                      fd (speed / r)
                    ]
                    [
                      face one-of neighbors with [(coast-patches = true) and (pxcor < [ xcor ] of myself)]
                      fd (speed / r)
                    ]
                  ]
                ]
              ]
              if (viia = true)
              [
                ifelse any? neighbors with [(coast-patches = true) and (pycor < [ round(ycor) ] of myself)]
                  [
                    face one-of neighbors with [(coast-patches = true) and (pycor < [ round(ycor) ] of myself)]
                    fd (speed / r)
                  ]
                  [
                    ifelse ((xcor > 7 and ycor < 22) or (xcor <= 7 and ycor >= 22 ))
                    [
                      face one-of neighbors with [(coast-patches = true) and (pxcor < [ xcor ] of myself)]
                      fd (speed / r)
                    ]
                    [
                      face one-of neighbors with [(coast-patches = true) and (pxcor > [ xcor ] of myself)]
                      fd (speed / r)
                    ]
                  ]
              ]
            ]
            [
              ifelse offshore-patches = true ; not neighboring a spawning patch
              [
;                ifelse (CS = true or C = true) and (any? patches with [ (CS = true or C = true) and (spawn-patches = true) ])
;                or (NS = true and any? patches with [ NS = true and spawn-patches = true ])
;                or (IS = true and any? patches with [ IS = true and spawn-patches = true ])
                ifelse (spawn-count > 0)
                [
                  face max-one-of neighbors with [offshore-patches = true ] [ sst ]
                  fd (speed / r)
                ]
                [
                  face min-one-of patches in-radius 5 with [coast-patches = true] [distance myself] ; changed as kept getting get stuck with code below

                  ;one-of neighbors with [ coast-patches = true ]



                  fd (speed / r)
                ]
              ]
              ; Non-sea patch - test!
              [
                face one-of neighbors with [ coast-patches = true ]
                fd (speed / r)
              ]
            ]
          ]
        ]
      ] ; end repeat

      if ([ spawn-patches ] of patch-here = true or spawn-count > 0)
      [
        set spawn-count spawn-count + 1
        if (spawn-count >= 60)
        [
          set spawn-trigger 0
          set spawn-count 0
          reallocate-site
        ]
      ]

      if ticks = 151
      [
        set spawn-trigger 0
        set spawn-count 0
        reallocate-site
      ]

    ] ; end if spawn-trigger
  ] ; end ask

end


to reallocate-site
  if random 101 > site-fidelity?         ; probability of a bass reallocating depends on the site-fidelity? slider
  [
    set ivbt 0
    set ivct 0
    set viiat 0
    set viidt 0
    set viiet 0
    set viifgt 0

    let decider random 6
    if decider = 0 [ set ivbt true ]
    if decider = 1 [ set ivct true ]
    if decider = 2 [ set viidt true ]
    if decider = 3 [ set viiet true ]
    if decider = 4 [ set viifgt true ]
    if decider = 5 [ set viiat true ]
 ]
end


to feeding-migration

  ask mat-bass
  [
    if (spawn-trigger = 0)
    [
      repeat r
      [
       ; Migrating
        if (IVbt = true or IVct = true)
        [
          set color cyan
          ifelse (ivb = true or ivc = true)   ; 1. Am I where I want to be (but may have to move between b and c)?
          [
           if (IVbt = true)
           [
             ifelse (ivb = true)
             [
               rt random 360 fd (speed / r)
             ]
             [
               ifelse (any? neighbors with [(coast-patches = true) and (pycor > [ round(ycor) ] of myself)])
               [
                 face one-of neighbors with [(coast-patches = true) and (pycor > [ round(ycor) ] of myself)]
                 fd (speed / r)
               ]
               [
                 ifelse (xcor < 23)
                 [
                   face one-of neighbors with [(coast-patches = true) and (pxcor < [ xcor ] of myself)]
                   fd (speed / r)
                 ]
                 [
                   face one-of neighbors with [(coast-patches = true) and (pxcor > [ xcor ] of myself)]
                   fd (speed / r)
                 ]
               ]
             ]
           ]
           if (IVct = true)
           [
             ifelse (ivc = true)
             [
               rt random 360 fd (speed / r)
             ]
             [  ; MAY NOT BE RELEVENT IF BUILDING IN A BREAK
               ifelse (any? neighbors with [(coast-patches = true) and (pycor < [ ycor ] of myself)])
                 [
                   ifelse ((pxcor = 33 or pxcor = 34) and pycor = 23) ; awkward patch!
                   [
                     face one-of neighbors with [(coast-patches = true) and (pxcor < [ xcor ] of myself)]
                     fd (speed / r)
                   ]
                   [
                     face one-of neighbors with [(coast-patches = true) and (pycor < [ ycor ] of myself)]
                     fd (speed / r)
                   ]
                 ]
                 [
                   ifelse (pxcor > 29) or (pycor > 33)
                   [
                     face one-of neighbors with [(coast-patches = true) and (pxcor < [ xcor ] of myself)]
                     fd (speed / r)
                   ]
                   [
                     face one-of neighbors with [(coast-patches = true) and (pxcor > [ xcor ] of myself)]
                     fd (speed / r)
                   ]
                 ]
             ]
           ]
          ]
          [
            ifelse any? neighbors with [(ivb = true) or (ivc = true)]   ; 2. Am I neighboring where I want to be?
            [
              face one-of neighbors with [(ivb = true) or (ivc = true)]
              fd (speed / r)
            ]
            [
              ifelse (offshore-patches = true)   ; 3. Am I offshore? Move to coast
              [
                if (NS = true or ycor > 23 ) ; as fish may travel to a random offshore patch north of NS spawning ground when heading to the coast
                [
                  ifelse (xcor > 23)
                  [
                    face one-of neighbors with [ pxcor > [ xcor ] of myself ]
                    fd (speed / r)
                  ]
                  [
                    face one-of neighbors with [ pxcor < [ xcor ] of myself]
                    fd (speed / r)
                  ]
                ]
                if (CS = true or C = true)
                [
                  ifelse (xcor > 11)
                  [
                    face one-of neighbors with [coast-patches = true]
                    fd (speed / r)
                  ]
                  [
                    ifelse (ycor > 7)
                    [
                      ifelse (any? neighbors with [(offshore-patches = true) and (pycor < [ ycor ] of myself)])
                      [
                        face one-of neighbors with [(offshore-patches = true) and (pycor < [ ycor ] of myself)]
                        fd (speed / r)
                      ]
                      [
                        face one-of neighbors with [(offshore-patches = true) and (pxcor < [ xcor ] of myself)]
                        fd (speed / r)
                      ]
                    ]
                    [
                      face one-of neighbors with [ pxcor > [ xcor ] of myself]
                      fd (speed / r)
                    ]
                  ]
                ]
                if (IS = true)
                [
                  face one-of neighbors with [sea = true and pycor < [ ycor ] of myself]
                  fd (speed / r)
                ]
              ]
              [
                ifelse (coast-patches = true)   ; Am I on a coastal patch but in an area I don't want to be?
                [
                  if (viid = true or viie = true)
                  [
                    ifelse (any? neighbors with [(coast-patches = true) and (pxcor > [ round xcor ] of myself)])
                    [
                      ifelse (pxcor = 13 and (pycor = 4 or pycor = 5))     ; One awkward patch that messes things up!
                      [
                        face one-of neighbors with [(coast-patches = true) and (pycor > [ ycor ] of myself)]
                        fd (speed / r)
                      ]
                      [
                        face one-of neighbors with [(coast-patches = true) and (pxcor > [ round xcor ] of myself)]
                        fd (speed / r)
                      ]
                    ]
                    [
                      face one-of neighbors with [(coast-patches = true) and (pycor > [ round ycor ] of myself)]
                      fd (speed / r)
                    ]
                  ]
                  if (viifg = true)   ; May end up here when returning from spawning ground (but shouldn't end up in as far as the Bristol Channel)
                  [
                    ifelse (any? neighbors with [(coast-patches = true) and (pycor < [ round ycor ] of myself)])
                    [
                      ifelse (pxcor = 7 and pycor = 9) ; stop fish jumping Cornwall!
                      [
                        face one-of neighbors with [(coast-patches = true) and (pxcor < [ round xcor ] of myself)]
                        fd (speed / r)
                      ]
                      [
                        face one-of neighbors with [(coast-patches = true) and (pycor < [ round ycor ] of myself)]
                        fd (speed / r)
                      ]
                    ]
                    [
                      ifelse (pxcor = 11 and pycor = 13)
                      [
                        face one-of neighbors with [(coast-patches = true) and (pxcor < [ xcor ] of myself)]
                        fd (speed / r)
                      ]
                      [
                        face one-of neighbors with [(sea = true) and (pycor < [ ycor ] of myself)]
                        fd (speed / r)
                      ]
                    ]
                  ]
                  if (viia = true)
                  [
                    ifelse any? neighbors with [(coast-patches = true) and (pycor < [ round(ycor) ] of myself)]
                      [
                        face one-of neighbors with [(coast-patches = true) and (pycor < [ round(ycor) ] of myself)]
                        fd (speed / r)
                      ]
                      [
                        ifelse ((xcor > 7 and ycor < 22) or (xcor <= 7 and ycor >= 22 ))
                        [
                          face one-of neighbors with [(coast-patches = true) and (pxcor < [ xcor ] of myself)]
                          fd (speed / r)
                        ]
                        [
                          face one-of neighbors with [(coast-patches = true) and (pxcor > [ xcor ] of myself)]
                          fd (speed / r)
                        ]
                      ]
                  ]
                ]
                [; Moved into non-sea patch
                  ifelse (pxcor = 7 and pycor = 8) ; Awkward patch
                  [
                    face one-of neighbors with [(coast-patches = true) and (pxcor > [ xcor ] of myself)]
                    fd (speed / r)
                  ]
                  [
                    ifelse (pxcor = 8 and pycor = 1) ; Awkward patch
                    [
                       face one-of neighbors with [(coast-patches = true) and (pycor > [ ycor ] of myself)]
                       fd (speed / r)
                    ]
                    [
                      ifelse (pxcor = 7 and pycor = 15) ; awkward patch
                      [
                        face one-of neighbors with [ (coast-patches = true) and pycor < [ ycor ] of myself ]
                        fd (speed / r)
                      ]
                      [
                        face one-of neighbors with [ coast-patches = true ]
                        fd (speed / r)
                      ]
                    ]
                  ]
                ]
              ]
            ]
          ]
        ] ; end of IVbc

        if (VIIdt = true or VIIet = true)
        [
          set color violet
          ifelse (viid = true or viie = true)   ; Am I where I want to be (but may have to move between d and e)?
          [
            if (VIIdt = true)
            [
              ifelse (viid = true)
              [
                rt random 360 fd (speed / r)
              ]
              [
                ifelse (any? neighbors with [(coast-patches = true) and (pxcor > [ round xcor ] of myself)])
                [
                  ifelse (pxcor = 13 and (pycor = 4 or pycor = 5))     ; One awkward patch that messes things up!
                  [
                    face one-of neighbors with [(coast-patches = true) and (pycor > [ ycor ] of myself)]
                    fd (speed / r)
                  ]
                  [
                    face one-of neighbors with [(coast-patches = true) and (pxcor > [ round xcor ] of myself)]
                    fd (speed / r)
                  ]
                ]
                [
                  face one-of neighbors with [(coast-patches = true) and (pycor > [ round ycor ] of myself)]
                  fd (speed / r)
                ]
              ]
            ]
            if (VIIet = true)
            [
              ifelse (viie = true)
              [
                rt random 360 fd (speed / r)
              ]
              [ ; MAY NOT BE RELEVANT IF BUILDING IN A BREAK
                ifelse (any? neighbors with [(coast-patches = true) and (pxcor < [ xcor ] of myself)])
                  [
                    face one-of neighbors with [(coast-patches = true) and (pxcor < [ xcor ] of myself)]
                  ]
                  [
                    ifelse (any? neighbors with [(coast-patches = true) and (pycor < [ ycor ] of myself)])
                    [
                      face one-of neighbors with [(coast-patches = true) and (pycor < [ ycor ] of myself)]
                    ]
                    [
                      face one-of neighbors with [(coast-patches = true and pycor > [ ycor ] of myself)] ; to deal with patch 16 5
                    ]
                  ]
                fd (speed / r)
              ]
            ]
          ]
          [
            ifelse any? neighbors with [(viid = true) or (viie = true)]   ; 2. Am I neighboring where I want to be?
            [
              ifelse (pxcor = 8 and pycor = 1) ; Awkward non-sea patch
                [
                   face one-of neighbors with [(coast-patches = true) and (pycor > [ ycor ] of myself)]
                   fd (speed / r)
                ]
                [
                  ifelse (pxcor = 7 and pycor = 9) ; stop fish jumping Cornwall!
                    [
                      face one-of neighbors with [(coast-patches = true) and (pxcor < [ round xcor ] of myself)]
                      fd (speed / r)
                    ]
                    [
                      face one-of neighbors with [(viid = true) or (viie = true)]
                      fd (speed / r)
                    ]
                ]

            ]
            [
              ifelse (offshore-patches = true)   ; Am I offshore? Move to coast
              [
               if (CS = true or C = true)
               [
                 ifelse (ycor > 7) ; Should only have to deal with CSC.
                 [
                   ifelse (any? neighbors with [(offshore-patches = true) and (pycor < [ ycor ] of myself)])
                     [
                       face one-of neighbors with [(offshore-patches = true) and (pycor < [ ycor ] of myself)]
                       fd (speed / r)
                     ]
                     [
                       face one-of neighbors with [(offshore-patches = true) and (pxcor < [ xcor ] of myself)]
                       fd (speed / r)
                     ]
                 ]
                 [
                   face one-of neighbors with [ pxcor > [ xcor ] of myself]
                   fd (speed / r)
                 ]
               ]
               if (NS = true or ycor > 23)
               [
                 face one-of neighbors with [(sea = true) and (pycor < [ ycor ] of myself)]
                 fd (speed / r)
               ]
               if (IS = true)
                [
                  face one-of neighbors with [sea = true and pycor < [ ycor ] of myself]
                  fd (speed / r)
                ]
              ] ; end offshre patches
              [
                ifelse (coast-patches = true)   ; Am I on a coastal patch but in an area I don't want to be?
                [
                  if (ivb = true)
                  [
                    ifelse any? neighbors with [(coast-patches = true) and (pycor < [ ycor ] of myself)]
                    [
                      ifelse ((pxcor = 33 or pxcor = 34) and pycor = 23) ; awkward patches!
                      [
                        face one-of neighbors with [(coast-patches = true) and (pxcor < [ xcor ] of myself)]
                        fd (speed / r)
                      ]
                      [
                        face one-of neighbors with [(coast-patches = true) and (pycor < [ ycor ] of myself)]
                        fd (speed / r)
                      ]
                    ]
                    [
                      ifelse (pxcor > 29) or (pycor > 33)
                      [
                        face one-of neighbors with [(coast-patches = true) and (pxcor < [ xcor ] of myself)]
                        fd (speed / r)
                      ]
                      [
                        face one-of neighbors with [(coast-patches = true) and (pxcor > [ xcor ] of myself)]
                        fd (speed / r)
                      ]
                    ]
                  ] ; end of ivb
                  if (ivc = true)
                  [
                    ifelse (any? neighbors with [(coast-patches = true) and (pycor < [ round(ycor) ] of myself)])
                    [
                      ifelse (pxcor = 27 or pxcor = 28) and (pycor = 19 or pycor = 20) ; Awkward patches!
                      [
                        ifelse (any? neighbors with [(coast-patches = true) and (pxcor < [ xcor ] of myself)])
                        [
                          face one-of neighbors with [(coast-patches = true) and (pxcor < [ xcor ] of myself)]
                          fd (speed / r)
                        ]
                        [
                          face one-of neighbors with [(coast-patches = true) and (pycor > [ ycor ] of myself)]
                          fd (speed / r)
                        ]
                      ]
                      [
                        ifelse (pxcor = 19 and pycor = 20) ; Awkward patch!
                        [
                          face one-of neighbors with [(coast-patches = true) and (pxcor > [ xcor ] of myself)]
                          fd (speed / r)
                        ]
                        [
                          face one-of neighbors with [(coast-patches = true) and (pycor < [ ycor ] of myself)]
                          fd (speed / r)
                        ]

                      ]
                    ]
                    [
                      ifelse (pxcor < 22 and pycor > 13)
                      [
                        face one-of neighbors with [(coast-patches = true) and (pxcor > [ xcor ] of myself)]
                        fd (speed / r)
                      ]
                      [
                        ifelse (pxcor = 28 and pycor = 18)
                        [
                          face one-of neighbors with [(coast-patches = true) and (pycor > [ ycor ] of myself)]
                          fd (speed / r)
                        ]
                        [
                          face one-of neighbors with [(coast-patches = true) and (pxcor < [ xcor ] of myself)]
                          fd (speed / r)
                        ]
                      ]
                    ]
                  ]

                  if (viifg = true)   ; May end up here when returning from spawning ground (but shouldn't end up in as far as the Bristol Channel)
                  [
                    ifelse (any? neighbors with [(coast-patches = true) and (pycor < [ round ycor ] of myself)])
                    [
                      face one-of neighbors with [(coast-patches = true) and (pycor < [ round ycor ] of myself)]
                      fd (speed / r)
                    ]
                    [
                      ifelse (pxcor = 11 and pycor = 13)
                      [
                        face one-of neighbors with [(coast-patches = true) and (pxcor < [ xcor ] of myself)]
                        fd (speed / r)
                      ]
                      [
                        face one-of neighbors with [(sea = true) and (pycor < [ ycor ] of myself)]
                        fd (speed / r)
                      ]
                    ]
                  ]
                 if (viia = true)
                   [
                     ifelse any? neighbors with [(coast-patches = true) and (pycor < [ round(ycor) ] of myself)]
                       [
                         face one-of neighbors with [(coast-patches = true) and (pycor < [ round(ycor) ] of myself)]
                         fd (speed / r)
                       ]
                       [
                         ifelse ((xcor > 7 and ycor < 22) or (xcor <= 7 and ycor >= 22 ))
                         [
                           face one-of neighbors with [(coast-patches = true) and (pxcor < [ xcor ] of myself)]
                           fd (speed / r)
                         ]
                         [
                           face one-of neighbors with [(coast-patches = true) and (pxcor > [ xcor ] of myself)]
                           fd (speed / r)
                         ]
                       ]
                   ]
                ]
                [
                  ; Moved into non-sea patch
                  ifelse (pxcor = 7 and pycor = 8) ; Awkward patch
                  [
                    face one-of neighbors with [(coast-patches = true) and (pxcor > [ xcor ] of myself)]
                    fd (speed / r)
                  ]
                  [
                    ifelse (pxcor = 7 and pycor = 15) ; awkward patch
                    [
                      face one-of neighbors with [ (coast-patches = true) and pycor < [ ycor ] of myself ]
                      fd (speed / r)
                    ]
                    [
                      face one-of neighbors with [ coast-patches = true ]
                      fd (speed / r)
                    ]
                  ]
                ]
              ]
            ]
          ]
        ] ; end of VIIde

        if (VIIfgt = true or VIIat = true)
        [
          set color orange
          ifelse (viia = true or viifg = true)   ; 1. Am I where I want to be (but may have to move between a and fg)?
          [
            if (VIIat = true)
            [
              ifelse (viia = true)
              [
                rt random 360 fd (speed / r)
              ]
              [
                ifelse (any? neighbors with [(coast-patches = true) and (pycor > [ round ycor ] of myself)])
                [
                  face one-of neighbors with [(coast-patches = true) and (pycor > [ round ycor ] of myself)]
                  fd (speed / r)
                ]
                [
                  ifelse (ycor > 12)
                  [
                    face one-of neighbors with [(coast-patches = true) and (pxcor < [ round xcor ] of myself)]
                    fd (speed / r)
                  ]
                  [
                    face one-of neighbors with [(coast-patches = true) and (pxcor > [ xcor ] of myself)]
                    fd (speed / r)
                  ]
                ]
              ]
            ]
            if (VIIfgt = true)
            [
              ifelse (viifg = true)
              [
                rt random 360
                fd (speed / r)
              ]
              [
                ifelse any? neighbors with [(coast-patches = true) and (pycor < [ round(ycor) ] of myself)]
                  [
                    face one-of neighbors with [(coast-patches = true) and (pycor < [ round(ycor) ] of myself)]
                    fd (speed / r)
                  ]
                  [
                    ifelse ((xcor > 7 and ycor < 22) or (xcor <= 7 and ycor >= 22 ))
                    [
                      face one-of neighbors with [(coast-patches = true) and (pxcor < [ xcor ] of myself)]
                      fd (speed / r)
                    ]
                    [
                      face one-of neighbors with [(coast-patches = true) and (pxcor > [ xcor ] of myself)]
                      fd (speed / r)
                    ]
                  ]
              ]
            ]
          ]
          [
            ifelse any? neighbors with [(viia = true) or (viifg = true)]   ; 2. Am I neighboring where I want to be?
              [
                ifelse (pxcor = 8 and pycor = 8) ; stop fish jumping Cornwall!
                [
                  face one-of neighbors with [(coast-patches = true) and pycor < [ ycor ] of myself]
                  fd (speed / r)
                ]
                [
                  face one-of neighbors with [(viia = true) or (viifg = true)]
                  fd (speed / r)
                ]
              ]
              [
                ifelse (offshore-patches = true)   ; Am I offshore? Move to coast
                [
                  if (IS = true)
                  [
                    ; shouldn't be an issue as all IS spawning patches neighbour a coast-patch
                  ]
                  if (CS = true or C = true)
                  [
                    ifelse (pxcor > 13)
                    [
                      face one-of neighbors with [(offshore-patches = true) and (pxcor < [ xcor ] of myself)]
                      fd (speed / r)
                    ]
                    [
                      ifelse (pxcor > 5)
                      [
                        face one-of neighbors with [(pxcor < [ xcor ] of myself and pycor >= [ ycor ] of myself)]
                        fd (speed / r)
                      ]
                      [
                        ifelse (pycor < 8)
                        [
                          face one-of neighbors with [pycor > [ ycor ] of myself]
                          fd (speed / r)
                        ]
                        [
                          face one-of neighbors with [(sea = true) and (pycor > [ ycor ] of myself or pxcor > [ xcor ] of myself)]
                          fd (speed / r)
                        ]
                      ]
                    ]
                  ]
                  if (NS = true or ycor > 23)
                  [
                    face one-of neighbors with [(sea = true) and (pycor < [ ycor ] of myself)]
                    fd (speed / r)
                  ]
                ]
                [
                  ifelse (coast-patches = true)
                  [
                    if (viie = true)
                    [
                      ifelse (pxcor > 5 and pxcor < 8 and pycor > 5)
                      [
                        face one-of neighbors with [(coast-patches = true) and (pycor > [ ycor ] of myself)]
                        fd (speed / r)
                      ]
                      [
                        ifelse pycor > 7
                        [
                          ifelse (pxcor = 11 and (pycor = 9 or pycor = 10)) ; bottleneck
                          [
                            face one-of neighbors with [(coast-patches = true) and (pycor < [ ycor ] of myself)]
                            fd (speed / r)
                          ]
                          [
                            ifelse (any? neighbors with [(coast-patches = true) and (pxcor < [ xcor ] of myself)])
                              [
                                face one-of neighbors with [(coast-patches = true) and (pxcor < [ xcor ] of myself)]
                                fd (speed / r)
                              ]
                              [
                                face one-of neighbors with [(coast-patches = true) and (pycor < [ ycor ] of myself)]
                                fd (speed / r)
                              ]
                          ]

                        ]
                        [
                          face one-of neighbors with [(pycor > [ ycor ] of myself) and (pxcor <= [ xcor ] of myself)]
                          fd (speed / r)
                        ]
                      ]
                    ]
                    if (viid = true)
                    [
                      ifelse (any? neighbors with [(coast-patches = true) and (pxcor < [ xcor ] of myself)])
                      [
                        face one-of neighbors with [(coast-patches = true) and (pxcor < [ xcor ] of myself)]
                      ]
                      [
                        ifelse (any? neighbors with [(coast-patches = true) and (pycor < [ ycor ] of myself)])
                        [
                          face one-of neighbors with [(coast-patches = true) and (pycor < [ ycor ] of myself)]
                        ]
                        [
                          face one-of neighbors with [(coast-patches = true and pycor > [ ycor ] of myself)] ; to deal with patch 16 5
                        ]
                      ]
                    fd (speed / r)
                    ]
                    if (ivb = true)
                    [
                      ifelse any? neighbors with [(coast-patches = true) and (pycor < [ ycor ] of myself)]
                      [
                        ifelse ((pxcor = 33 or pxcor = 34) and pycor = 23) ; awkward patches!
                        [
                          face one-of neighbors with [(coast-patches = true) and (pxcor < [ xcor ] of myself)]
                          fd (speed / r)
                        ]
                        [
                          face one-of neighbors with [(coast-patches = true) and (pycor < [ ycor ] of myself)]
                          fd (speed / r)
                        ]
                      ]
                      [
                        ifelse (pxcor > 29) or (pycor > 33)
                        [
                          face one-of neighbors with [(coast-patches = true) and (pxcor < [ xcor ] of myself)]
                          fd (speed / r)
                        ]
                        [
                          face one-of neighbors with [(coast-patches = true) and (pxcor > [ xcor ] of myself)]
                          fd (speed / r)
                        ]
                      ]
                    ] ; end of ivb
                    if (ivc = true)
                    [
                      ifelse (any? neighbors with [(coast-patches = true) and (pycor < [ round(ycor) ] of myself)])
                      [
                        ifelse (pxcor = 27 or pxcor = 28) and (pycor = 19 or pycor = 20) ; Awkward patches!
                        [
                          ifelse (any? neighbors with [(coast-patches = true) and (pxcor < [ xcor ] of myself)])
                          [
                            face one-of neighbors with [(coast-patches = true) and (pxcor < [ xcor ] of myself)]
                            fd (speed / r)
                          ]
                          [
                            face one-of neighbors with [(coast-patches = true) and (pycor > [ ycor ] of myself)]
                            fd (speed / r)
                          ]
                        ]
                        [
                          ifelse (pxcor = 19 and pycor = 20) ; Awkward patch!
                          [
                            face one-of neighbors with [(coast-patches = true) and (pxcor > [ xcor ] of myself)]
                            fd (speed / r)
                          ]
                          [
                            face one-of neighbors with [(coast-patches = true) and (pycor < [ ycor ] of myself)]
                            fd (speed / r)
                          ]

                        ]
                      ]
                      [
                        ifelse (pxcor < 22 and pycor > 13)
                        [
                          face one-of neighbors with [(coast-patches = true) and (pxcor > [ xcor ] of myself)]
                          fd (speed / r)
                        ]
                        [
                          ifelse (pxcor = 28 and pycor = 18)
                          [
                            face one-of neighbors with [(coast-patches = true) and (pycor > [ ycor ] of myself)]
                            fd (speed / r)
                          ]
                          [
                            face one-of neighbors with [(coast-patches = true) and (pxcor < [ xcor ] of myself)]
                            fd (speed / r)
                          ]
                        ]
                      ]
                    ]
                  ]
                  [
                    face one-of neighbors with [ coast-patches = true ]    ; Moved into non-sea patch
                    fd (speed / r)
                  ]
                ]
              ]
          ]
        ] ; end of VIIafg
      ] ; end repeat
    ] ; end if spawn-trigger
  ] ; end ask

end



to avoid-corwall
 if (pxcor = 8 and pycor = 8) ; stop fish jumping Cornwall!
                [
                  face one-of neighbors with [(coast-patches = true) and pycor < [ ycor ] of myself]
                  fd (speed / r)
                ]

  if (pxcor = 7 and pycor = 9) ; stop fish jumping Cornwall!
                      [
                        face one-of neighbors with [(coast-patches = true) and (pxcor < [ xcor ] of myself)]
                        fd (speed / r)
                      ]

 if (pxcor = 7 and pycor = 8) ; stop fish jumping Cornwall!
                      [
                        face one-of neighbors with [(coast-patches = true) and pycor < [ ycor ] of myself]
                        fd (speed / r)
                      ]

end





to local-movement ; need to add in exclusion to stop jumping cornwall.
  ask juv-bass

 [
    ifelse L < 32
    [
      repeat r
      [ifelse nursery-patches = 0
        [ rt  180 fd (speed / r) ]
        [rt random 360 fd (speed / r)]
        ]]

    [
      repeat r
      [ ifelse coast-patches = 0
        [ rt  180 fd (speed / r) ]
        [rt random 360 fd (speed / r)]
      ]]

      ]
end





to drift_eggs

ask eggs [

 if ((ivbt = true) and (IVb_target = true ))
        or ((ivct = true) and (IVc_target = true ))
        or ((VIIdt = true) and (VIId_target = true ))
        or ((VIIet = true) and (VIIe_target = true ))
        or ((VIIfgt = true) and (VIIfg_target = true ))
        or ((VIIat = true) and (VIIa_target = true ))

[set reach_target  True ]

ifelse reach_target = True

[   repeat r
      [
        ifelse nursery-patches = 0
        [ rt  180 fd (speed / r) ]
        [rt random 360 fd (speed / r)]
      ]]
[


if ivbt = true ;and ivct = false and VIIdt = false and VIIet = false and VIIfgt = false and VIIat = false
[
 face one-of patches with [IVb_target = true ]
 move-to one-of neighbors with [IVb_dist = ([IVb_dist] of myself) - 1]
]


if  ivct = true ;and VIIdt = false and VIIet = false and VIIfgt = false and VIIat = false and ivbt = false
[
 face one-of patches with [IVc_target = true ]
 move-to  one-of neighbors with [IVc_dist = ([IVc_dist] of myself) - 1]
]


if VIIdt = true ;and ivbt = false and ivct = false  and VIIet = false and VIIfgt = false and VIIat = false
[
 face one-of patches with [VIId_target = true ]
 move-to one-of neighbors with [VIId_dist = ([VIId_dist] of myself) - 1]
]


if VIIet = true ;and ivbt = false and ivct = false and VIIdt = false and VIIfgt = false and VIIat = false
[
 face one-of patches with [VIIe_target = true ]
 move-to one-of neighbors with [VIIe_dist = ([VIIe_dist] of myself) - 1]
]


if VIIfgt = true ;and ivbt = false and ivct = false and VIIdt = false and VIIet = false and VIIat = false
[
 face one-of patches with [VIIfg_target = true ]
 move-to one-of neighbors with [VIIfg_dist = ([VIIfg_dist] of myself) - 1]
]


if VIIat = true ;and ivbt = false and ivct = false and VIIdt = false and VIIet = false and VIIfgt = false
[
 face one-of patches with [VIIa_target = true ]
 move-to one-of neighbors with [VIIa_dist = ([VIIa_dist] of myself) - 1]
]

  ] ]

end




to drift_ys-larvae

ask ys-larvae [

 if ((ivbt = true) and (IVb_target = true ))
        or ((ivct = true) and (IVc_target = true ))
        or ((VIIdt = true) and (VIId_target = true ))
        or ((VIIet = true) and (VIIe_target = true ))
        or ((VIIfgt = true) and (VIIfg_target = true ))
        or ((VIIat = true) and (VIIa_target = true ))

[set reach_target  True ]

ifelse reach_target = True
[
      repeat r
      [
        rt random 360 fd (speed / r)
        ifelse nursery-patches = 0
        [ rt  180 fd (speed / r) ]
        [rt random 360 fd (speed / r)]

      ]]

[

if ivbt = true
[
 face one-of patches with [IVb_target = true ]
 move-to one-of neighbors with [IVb_dist = ([IVb_dist] of myself) - 1]
]


if  ivct = true
[
 face one-of patches with [IVc_target = true ]
 move-to  one-of neighbors with [IVc_dist = ([IVc_dist] of myself) - 1]
]


if VIIdt = true
[
 face one-of patches with [VIId_target = true ]
 move-to one-of neighbors with [VIId_dist = ([VIId_dist] of myself) - 1]
]


if VIIet = true
[
 face one-of patches with [VIIe_target = true ]
 move-to one-of neighbors with [VIIe_dist = ([VIIe_dist] of myself) - 1]
]


if VIIfgt = true
[
 face one-of patches with [VIIfg_target = true ]
 move-to one-of neighbors with [VIIfg_dist = ([VIIfg_dist] of myself) - 1]
]


if VIIat = true
[
 face one-of patches with [VIIa_target = true ]
 move-to one-of neighbors with [VIIa_dist = ([VIIa_dist] of myself) - 1]
]

  ] ]

end





to drift_larvae

ask larvae[

 if ((ivbt = true) and (IVb_target = true ))
        or ((ivct = true) and (IVc_target = true ))
        or ((VIIdt = true) and (VIId_target = true ))
        or ((VIIet = true) and (VIIe_target = true ))
        or ((VIIfgt = true) and (VIIfg_target = true ))
        or ((VIIat = true) and (VIIa_target = true ))

[set reach_target  True ]

ifelse reach_target = True

[
    repeat r
      [ ifelse nursery-patches = 0
        [ rt  180 fd (speed / r) ]
        [rt random 360 fd (speed / r)]
      ]
    ]

[

if ivbt = true ;and ivct = false and VIIdt = false and VIIet = false and VIIfgt = false and VIIat = false
[
 face one-of patches with [IVb_target = true ]
 move-to one-of neighbors with [IVb_dist = ([IVb_dist] of myself) - 1]
]


if  ivct = true ;and VIIdt = false and VIIet = false and VIIfgt = false and VIIat = false and ivbt = false
[
 face one-of patches with [IVc_target = true ]
 move-to  one-of neighbors with [IVc_dist = ([IVc_dist] of myself) - 1]
]


if VIIdt = true ;and ivbt = false and ivct = false  and VIIet = false and VIIfgt = false and VIIat = false
[
 face one-of patches with [VIId_target = true ]
 move-to one-of neighbors with [VIId_dist = ([VIId_dist] of myself) - 1]
]


if VIIet = true ;and ivbt = false and ivct = false and VIIdt = false and VIIfgt = false and VIIat = false
[
 face one-of patches with [VIIe_target = true ]
 move-to one-of neighbors with [VIIe_dist = ([VIIe_dist] of myself) - 1]
]


if VIIfgt = true ;and ivbt = false and ivct = false and VIIdt = false and VIIet = false and VIIat = false
[
 face one-of patches with [VIIfg_target = true ]
 move-to one-of neighbors with [VIIfg_dist = ([VIIfg_dist] of myself) - 1]
]


if VIIat = true ;and ivbt = false and ivct = false and VIIdt = false and VIIet = false and VIIfgt = false
[
 face one-of patches with [VIIa_target = true ]
 move-to one-of neighbors with [VIIa_dist = ([VIIa_dist] of myself) - 1]
]

  ] ]

end









;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


to calc_sst_mean
ask patches with [spawn-patches = true]
 [set sst-output  mean [sst] of patches with [sst >= 0]]
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



to recreational-inshore-F
    ;file-open "C:/Users/dm820646/Dropbox/abc_shared/IBM/data/RecreationalinshoreF.txt"
     file-open "data/2020/RecreationalinshoreF.txt"
    while [not file-at-end?]
    [ let fyear file-read
      if (fyear = year)
      [ let counter 0
        while [counter <= 30]
        [ let FFRI file-read
          ask turtles with [cohort = counter and l > angler-min-size]
          [set FRi FFRI]
          set counter counter + 1
        ]]]
    file-close
end






to commercial-inshore-F
;  file-open "C:/Users/dm820646/Dropbox/abc_shared/IBM/data/CommercialinshoreF.txt"
       file-open "data/2020/CommercialinshoreF.txt"
    while [not file-at-end?]
    [ let fyear file-read
      if (fyear = year)
      [ let counter 0
        while [counter <= 30]
        [ let FFCI file-read
          ask turtles with [cohort = counter and l > com-net-hole-size]
          [set FCi FFCI]
          set counter counter + 1
        ]]]
    file-close
end




to commercial-offshore-F
;file-open "C:/Users/dm820646/Dropbox/abc_shared/IBM/data/CommercialoffshoreF.txt"
     file-open "data/2020/CommercialoffshoreF.txt"
    while [not file-at-end?]
    [ let fyear file-read
      if (fyear = year)
      [ let counter 0
        while [counter <= 30]
        [ let FFCO file-read
          ask turtles with [cohort = counter and l > com-net-hole-size]
          [set FCo FFCO]
          set counter counter + 1
        ]]]
    file-close
end





to fishing-mortality
;  recreational-inshore-F
;  commercial-inshore-F
;  commercial-offshore-F
  ask turtles with [(coast-patches = true or sea = false) and (l > angler-min-size)]
  [ set number number * exp ( - FRi / 365) ]

  ;ifelse (ticks > 90 and ticks < 305)
  ;[
    ask turtles with [(coast-patches = true or sea = false) and (l > com-net-hole-size)]
    [set number number * exp ( - (FCi * 365 / 214) / 365) ]
  ;]
  ;[
    ask turtles with [(offshore-patches = true) and ( l > com-net-hole-size )]
    [ set number number * exp ( - (FCo * 365 / 151) / 365)]
  ;]
end




to calculate-catch
 ask patches with [offshore-patches = true]
      [set catch catch + sum [ (((FCo * 365 / 151 ) / 365) / ((FCo * 365 / 151 ) / 365 + Am )) * number * (1 - exp (- (Am + (FCo * 365 / 151 ) / 365))) * total-mass / 1000 ] of turtles-here with [ l > com-net-hole-size ]] ; Baranov catch equation

 ask patches with [coast-patches = true]
     [set catch catch + sum [ (((FCi * 365 / 214 ) / 365 + FRi / 365 ) / ((FCi * 365 / 214 ) / 365 + FRi / 365 + Am )) * number * (1 - exp (- (Am  + (FCi * 365 / 214 ) / 365 + FRi / 365 ))) * total-mass / 1000 ] of turtles-here] ; Baranov catch equation. MLS restrictions now specified when reading in the data.

end







to icesssbtsb                                                                                 ; reads in assessment abundance, ssb and tsb
  file-open "data/2020/icesnumssbtsb check.txt"
  while [not file-at-end?]
  [if year = file-read
  [set icesnum file-read
   set icesssb file-read
   set icestsb file-read]]
  file-close
end







to plot-age-distribution                                                                      ; age distribution plot
  set-current-plot "age distribution"
  set-current-plot-pen "model"
  clear-plot
  let countz 0
  while [countz < 30]
  [ plotxy countz ((sum [number] of (turtles with [cohort = countz])) / 1000)
    set countz (countz + 1) ]
  set-current-plot-pen "ices"
  file-open "data/2020/icesagedistribution.txt"
  while [not file-at-end?]
  [ let yearz file-read
    if yearz = year
    [ let countr 0
      while [countr < 30]
      [ let point file-read
        plotxy countr (point / 1000)
        set countr (countr + 1)]]]
  file-close
end








to recruits-read-in                    ; reads in the number of number of recruits estimated by ICES
  file-open "data/2020/icesagedistribution.txt"
  while [not file-at-end?]
  [
    if year + 1 = file-read
    [ set Rest file-read ]
  ]
  file-close
end






to spin_spawn


ask mat-bass [

      calc-reproduction
      ifelse (energy-reserve - maintenance-energy) >= max-R     ;if there is enough energy after the maintenance energy (see "calc-reproduction") is set aside, the maximum amount is sent to the reproduction store R

      [
        set energy-reserve energy-reserve - max-R                          ; the energy needed to produce the max number of eggs is taken from the energy reserve
        set gonad-mass  (max-R / (Ef + Fs))      / 1000                         ;gonad mass is calculated from the energy that went in to producing the eggs
        set realised-fecundity potential-fecundity                         ; realised fecundity is set to max potential fecundity

     ]
 [

        set ER energy-reserve - maintenance-energy ;energy reserve is what is left after maintenence is taken out
        set gonad-mass (ER) / (Ef + Fs) / 1000
        set realised-fecundity (ER / (max-R)) * potential-fecundity
      ]
       ]


  if (ticks = 152) [ set n item 0 nspawn ]
  if (ticks = 182) [ set n item 1 nspawn ]
  if (ticks = 213) [ set n item 2 nspawn ]
  if (ticks = 244) [ set n item 3 nspawn ]


  if n > 0

 [


   ; Recruits taken from the stock assessment
  recruits-read-in
    set y Rest

  if (ticks = 152) [set loc map round map [? * n ] map [? / sum locfeb] locfeb]
  if (ticks = 182) [set loc map round map [? * n ] map [? / sum locmar] locmar]
  if (ticks = 213) [set loc map round map [? * n ] map [? / sum locapr] locapr]
  if (ticks = 244) [set loc map round map [? * n ] map [? / sum locmay] locmay]

  if n - sum loc != 0
  [
    let x n - sum loc
    let index position max loc loc
    set loc replace-item index loc (item index loc + x)
  ]

  let loopcount 1

  foreach loc
  [
    if ? > 0
    [
      if (loopcount = 1) ; North Sea
      [
        let counter 0
        while [counter < ?]
        [
          let where random 100
          ifelse where < 3
          [
            ask one-of patches with [nursery-patches = true and ivb = true]
            [
              sprout 1 [
                set breed juv-bass
                set size 0.4
                set shape "fish"
                set age 100 / 365
                set cohort -1
                set color yellow
                set number y / ncohort
                set L linf * (1 - exp(- K * (age - t0)))
                set std-mass a_g * (L ^ (b_g))   ;; this really makes total mass...so is too mauch to then add resrves to
                set structural-mass std-mass
                set energy-reserve-max ((structural-mass * 0.01) * El)
                set energy-reserve-max-mass (energy-reserve-max / El)
                set energy-reserve (energy-reserve-max * 0.5)
                set total-mass (structural-mass + (energy-reserve / El))
                set maintenance-energy (energy-reserve * 0.1    )
                set potential-fecundity (std-mass ) * (eggs_per_bass)                      ; maximum number of eggs that can be produced, potential fecundity, depends on body length , changed for bass JW pickett and pawson fecundity estimate.
                set max-R (potential-fecundity * egg_mass * (Ef + Fs))
                calculate-speed
                calculate-r

              ]
            ]
          ]
          [
            ifelse where < 94
            [
              ask one-of patches with [nursery-patches = true and ivc = true]
              [
                sprout 1 [
                  set breed juv-bass
                  set size 0.4
                  set shape "fish"
                  set age 100 / 365
                  set cohort -1
                  set color yellow
                  set number y / ncohort
                  set L linf * (1 - exp(- K * (age - t0)))
                  set std-mass a_g * (L ^ (b_g))   ;; this really makes total mass...so is too mauch to then add resrves to
                  set structural-mass std-mass
                  set energy-reserve-max ((structural-mass * 0.01) * El)
                  set energy-reserve-max-mass (energy-reserve-max / El)
                  set energy-reserve (energy-reserve-max * 0.5)
                  set total-mass (structural-mass + (energy-reserve / El))
                  set maintenance-energy (energy-reserve * 0.1    )
                  set potential-fecundity (std-mass ) * (eggs_per_bass)                      ; maximum number of eggs that can be produced, potential fecundity, depends on body length , changed for bass JW pickett and pawson fecundity estimate.
                  set max-R (potential-fecundity * egg_mass * (Ef + Fs))
                  calculate-speed
                  calculate-r
                ]
              ]
            ]
            [
              ask one-of patches with [nursery-patches = true and viia = true]
              [
                sprout 1 [
                  set breed juv-bass
                  set size 0.4
                  set shape "fish"
                  set age 100 / 365
                  set cohort -1
                  set color yellow
                  set number y / ncohort
                  set L linf * (1 - exp(- K * (age - t0)))
                  set std-mass a_g * (L ^ (b_g))   ;; this really makes total mass...so is too mauch to then add resrves to
                  set structural-mass std-mass
                  set energy-reserve-max ((structural-mass * 0.01) * El)
                  set energy-reserve-max-mass (energy-reserve-max / El)
                  set energy-reserve (energy-reserve-max * 0.5)
                  set total-mass (structural-mass + (energy-reserve / El))
                  set maintenance-energy (energy-reserve * 0.1    )
                  set potential-fecundity (std-mass ) * (eggs_per_bass)                      ; maximum number of eggs that can be produced, potential fecundity, depends on body length , changed for bass JW pickett and pawson fecundity estimate.
                  set max-R (potential-fecundity * egg_mass * (Ef + Fs))
                  calculate-speed
                  calculate-r
                ]
              ]
            ]
          ]
          set counter counter + 1
        ]
      ]
      if (loopcount = 2) ; Irish Sea
      [
        let counter 0
        while [counter < ?]
        [
          let where random 100
          ifelse where < 5
          [
            ask one-of patches with [nursery-patches = true and ivc = true]
            [
              sprout 1 [
                set breed juv-bass
                set size 0.4
                set shape "fish"
                set age 100 / 365
                set cohort -1
                set color yellow
                set number y / ncohort
                set L linf * (1 - exp(- K * (age - t0)))
                set std-mass a_g * (L ^ (b_g))   ;; this really makes total mass...so is too mauch to then add resrves to
                set structural-mass std-mass
                set energy-reserve-max ((structural-mass * 0.01) * El)
                set energy-reserve-max-mass (energy-reserve-max / El)
                set energy-reserve (energy-reserve-max * 0.5)
                set total-mass (structural-mass + (energy-reserve / El))
                set maintenance-energy (energy-reserve * 0.1    )
                set potential-fecundity (std-mass ) * (eggs_per_bass)                      ; maximum number of eggs that can be produced, potential fecundity, depends on body length , changed for bass JW pickett and pawson fecundity estimate.
                set max-R (potential-fecundity * egg_mass * (Ef + Fs))
                calculate-speed
                calculate-r
              ]
            ]
          ]
          [
            ifelse where < 98
            [
              ask one-of patches with [nursery-patches = true and viia = true]
              [
                sprout 1 [
                  set breed juv-bass
                  set size 0.4
                  set shape "fish"
                  set age 100 / 365
                  set cohort -1
                  set color yellow
                  set number y / ncohort
                  set L linf * (1 - exp(- K * (age - t0)))
                  set std-mass a_g * (L ^ (b_g))   ;; this really makes total mass...so is too mauch to then add resrves to
                  set structural-mass std-mass
                  set energy-reserve-max ((structural-mass * 0.01) * El)
                  set energy-reserve-max-mass (energy-reserve-max / El)
                  set energy-reserve (energy-reserve-max * 0.5)
                  set total-mass (structural-mass + (energy-reserve / El))
                  set maintenance-energy (energy-reserve * 0.1    )
                  set potential-fecundity (std-mass ) * (eggs_per_bass)                      ; maximum number of eggs that can be produced, potential fecundity, depends on body length , changed for bass JW pickett and pawson fecundity estimate.
                  set max-R (potential-fecundity * egg_mass * (Ef + Fs))
                  calculate-speed
                  calculate-r
                ]
              ]
            ]
            [
              ifelse where < 99
              [
                ask one-of patches with [nursery-patches = true and viid = true]
                [
                  sprout 1 [
                    set breed juv-bass
                    set size 0.4
                    set shape "fish"
                    set age 100 / 365
                    set cohort -1
                    set color yellow
                    set number y / ncohort
                    set L linf * (1 - exp(- K * (age - t0)))
                    set std-mass a_g * (L ^ (b_g))   ;; this really makes total mass...so is too mauch to then add resrves to
                    set structural-mass std-mass
                    set energy-reserve-max ((structural-mass * 0.01) * El)
                    set energy-reserve-max-mass (energy-reserve-max / El)
                    set energy-reserve (energy-reserve-max * 0.5)
                    set total-mass (structural-mass + (energy-reserve / El))
                    set maintenance-energy (energy-reserve * 0.1    )
                    set potential-fecundity (std-mass ) * (eggs_per_bass)                      ; maximum number of eggs that can be produced, potential fecundity, depends on body length , changed for bass JW pickett and pawson fecundity estimate.
                    set max-R (potential-fecundity * egg_mass * (Ef + Fs))
                    calculate-speed
                    calculate-r
                  ]
                ]
              ]
              [
                ask one-of patches with [nursery-patches = true and viifg = true]
                [
                  sprout 1 [
                    set breed juv-bass
                    set size 0.4
                    set shape "fish"
                    set age 100 / 365
                    set cohort -1
                    set color yellow
                    set number y / ncohort
                    set L linf * (1 - exp(- K * (age - t0)))
                    set std-mass a_g * (L ^ (b_g))   ;; this really makes total mass...so is too mauch to then add resrves to
                    set structural-mass std-mass
                    set energy-reserve-max ((structural-mass * 0.01) * El)
                    set energy-reserve-max-mass (energy-reserve-max / El)
                    set energy-reserve (energy-reserve-max * 0.5)
                    set total-mass (structural-mass + (energy-reserve / El))
                    set maintenance-energy (energy-reserve * 0.1    )
                    set potential-fecundity (std-mass ) * (eggs_per_bass)                      ; maximum number of eggs that can be produced, potential fecundity, depends on body length , changed for bass JW pickett and pawson fecundity estimate.
                    set max-R (potential-fecundity * egg_mass * (Ef + Fs))
                    calculate-speed
                    calculate-r
                  ]
                ]
              ]
            ]
          ]
          set counter counter + 1
        ]
      ]
      if (loopcount = 3) ; Channel
      [
        let counter 0
        while [counter < ?]
        [
          let where random 100
          ifelse where < 4
          [
            ask one-of patches with [nursery-patches = true and ivc = true]
            [
              sprout 1 [
                set breed juv-bass
                set size 0.4
                set shape "fish"
                set age 100 / 365
                set cohort -1
                set color yellow
                set number y / ncohort
                set L linf * (1 - exp(- K * (age - t0)))
                  set std-mass a_g * (L ^ (b_g))   ;; this really makes total mass...so is too mauch to then add resrves to
                set structural-mass std-mass
                set energy-reserve-max ((structural-mass * 0.01) * El)
                set energy-reserve-max-mass (energy-reserve-max / El)
                set energy-reserve (energy-reserve-max * 0.5)
                set total-mass (structural-mass + (energy-reserve / El))
                set maintenance-energy (energy-reserve * 0.1    )
                set potential-fecundity (std-mass ) * (eggs_per_bass)                      ; maximum number of eggs that can be produced, potential fecundity, depends on body length , changed for bass JW pickett and pawson fecundity estimate.
                set max-R (potential-fecundity * egg_mass * (Ef + Fs))
                calculate-speed
                calculate-r
              ]
            ]
          ]
          [
            ifelse where < 5
            [
              ask one-of patches with [nursery-patches = true and viia = true]
              [
                sprout 1 [
                  set breed juv-bass
                  set size 0.4
                  set shape "fish"
                  set age 100 / 365
                  set cohort -1
                  set color yellow
                  set number y / ncohort
                  set L linf * (1 - exp(- K * (age - t0)))
                  set std-mass a_g * (L ^ (b_g))   ;; this really makes total mass...so is too mauch to then add resrves to
                  set structural-mass std-mass
                  set energy-reserve-max ((structural-mass * 0.01) * El)
                  set energy-reserve-max-mass (energy-reserve-max / El)
                  set energy-reserve (energy-reserve-max * 0.5)
                  set total-mass (structural-mass + (energy-reserve / El))
                  set maintenance-energy (energy-reserve * 0.1    )
                  set potential-fecundity (std-mass ) * (eggs_per_bass)                      ; maximum number of eggs that can be produced, potential fecundity, depends on body length , changed for bass JW pickett and pawson fecundity estimate.
                  set max-R (potential-fecundity * egg_mass * (Ef + Fs))
                  calculate-speed
                  calculate-r
                ]
              ]
            ]
            [
              ifelse where < 27
              [
                ask one-of patches with [nursery-patches = true and viid = true]
                [
                  sprout 1 [
                    set breed juv-bass
                    set size 0.4
                    set shape "fish"
                    set age 100 / 365
                    set cohort -1
                    set color yellow
                    set number y / ncohort
                    set L linf * (1 - exp(- K * (age - t0)))
                    set std-mass a_g * (L ^ (b_g))   ;; this really makes total mass...so is too mauch to then add resrves to
                    set structural-mass std-mass
                    set energy-reserve-max ((structural-mass * 0.01) * El)
                    set energy-reserve-max-mass (energy-reserve-max / El)
                    set energy-reserve (energy-reserve-max * 0.5)
                    set total-mass (structural-mass + (energy-reserve / El))
                    set maintenance-energy (energy-reserve * 0.1    )
                    set potential-fecundity (std-mass ) * (eggs_per_bass)                      ; maximum number of eggs that can be produced, potential fecundity, depends on body length , changed for bass JW pickett and pawson fecundity estimate.
                    set max-R (potential-fecundity * egg_mass * (Ef + Fs))
                    calculate-speed
                    calculate-r
                  ]
                ]
              ]
              [
                ifelse where < 72
                [
                  ask one-of patches with [nursery-patches = true and viie = true]
                  [
                    sprout 1 [
                      set breed juv-bass
                      set size 0.4
                      set shape "fish"
                      set age 100 / 365
                      set cohort -1
                      set color yellow
                      set number y / ncohort
                      set L linf * (1 - exp(- K * (age - t0)))
                      set std-mass a_g * (L ^ (b_g))   ;; this really makes total mass...so is too mauch to then add resrves to
                      set structural-mass std-mass
                      set energy-reserve-max ((structural-mass * 0.01) * El)
                      set energy-reserve-max-mass (energy-reserve-max / El)
                      set energy-reserve (energy-reserve-max * 0.5)
                      set total-mass (structural-mass + (energy-reserve / El))
                      set maintenance-energy (energy-reserve * 0.1    )
                      set potential-fecundity (std-mass ) * (eggs_per_bass)                      ; maximum number of eggs that can be produced, potential fecundity, depends on body length , changed for bass JW pickett and pawson fecundity estimate.
                      set max-R (potential-fecundity * egg_mass * (Ef + Fs))
                      calculate-speed
                      calculate-r
                    ]
                  ]
                ]
                [
                  ask one-of patches with [nursery-patches = true and viifg = true]
                  [
                    sprout 1 [
                      set breed juv-bass
                      set size 0.4
                      set shape "fish"
                      set age 60 / 365
                      set cohort -1
                      set color yellow
                      set number y / ncohort
                      set L linf * (1 - exp(- K * (age - t0)))
                      set std-mass a_g * (L ^ (b_g))   ;; this really makes total mass...so is too mauch to then add resrves to
                      set structural-mass std-mass
                      set energy-reserve-max ((structural-mass * 0.01) * El)
                      set energy-reserve-max-mass (energy-reserve-max / El)
                      set energy-reserve (energy-reserve-max * 0.5)
                      set total-mass (structural-mass + (energy-reserve / El))
                      set maintenance-energy (energy-reserve * 0.1    )
                      set potential-fecundity (std-mass ) * (eggs_per_bass)                      ; maximum number of eggs that can be produced, potential fecundity, depends on body length , changed for bass JW pickett and pawson fecundity estimate.
                      set max-R (potential-fecundity * egg_mass * (Ef + Fs))
                      calculate-speed
                      calculate-r
                    ]
                  ]
                ]
              ]
            ]
          ]
          set counter counter + 1
        ]
      ]
      if (loopcount = 4) ; Celtic Sea
      [
        let counter 0
        while [counter < ?]
        [
          let where random 100
          ifelse where < 38
          [
            ask one-of patches with [nursery-patches = true and viia = true]
            [
              sprout 1 [
                set breed juv-bass
                set size 0.4
                set shape "fish"
                set age 100 / 365
                set cohort -1
                set color yellow
                set number y / ncohort
                set L linf * (1 - exp(- K * (age - t0)))
                set std-mass a_g * (L ^ (b_g))   ;; this really makes total mass...so is too mauch to then add resrves to
                set structural-mass std-mass
                set energy-reserve-max ((structural-mass * 0.01) * El)
                set energy-reserve-max-mass (energy-reserve-max / El)
                set energy-reserve (energy-reserve-max * 0.5)
                set total-mass (structural-mass + (energy-reserve / El))
                set maintenance-energy (energy-reserve * 0.1    )
                set potential-fecundity (std-mass ) * (eggs_per_bass)                      ; maximum number of eggs that can be produced, potential fecundity, depends on body length , changed for bass JW pickett and pawson fecundity estimate.
                set max-R (potential-fecundity * egg_mass * (Ef + Fs))
                calculate-speed
                calculate-r
              ]
            ]
          ]
          [
          ifelse where < 40
          [
            ask one-of patches with [nursery-patches = true and viie = true]
            [
              sprout 1 [
                set breed juv-bass
                set size 0.4
                set shape "fish"
                set age 100 / 365
                set cohort -1
                set color yellow
                set number y / ncohort
                set L linf * (1 - exp(- K * (age - t0)))
                set std-mass a_g * (L ^ (b_g))   ;; this really makes total mass...so is too mauch to then add resrves to
                set structural-mass std-mass
                set energy-reserve-max ((structural-mass * 0.01) * El)
                set energy-reserve-max-mass (energy-reserve-max / El)
                set energy-reserve (energy-reserve-max * 0.5)
                set total-mass (structural-mass + (energy-reserve / El))
                set maintenance-energy (energy-reserve * 0.1    )
                set potential-fecundity (std-mass ) * (eggs_per_bass)                      ; maximum number of eggs that can be produced, potential fecundity, depends on body length , changed for bass JW pickett and pawson fecundity estimate.
                set max-R (potential-fecundity * egg_mass * (Ef + Fs))
                calculate-speed
                calculate-r
              ]
            ]
          ]
          [
            ask one-of patches with [nursery-patches = true and viifg = true]
            [
              sprout 1 [
                set breed juv-bass
                set size 0.4
                set shape "fish"
                set age 100 / 365
                set cohort -1
                set color yellow
                set number y / ncohort
                set L linf * (1 - exp(- K * (age - t0)))
                set std-mass a_g * (L ^ (b_g))   ;; this really makes total mass...so is too mauch to then add resrves to
                set structural-mass std-mass
                set energy-reserve-max ((structural-mass * 0.01) * El)
                set energy-reserve-max-mass (energy-reserve-max / El)
                set energy-reserve (energy-reserve-max * 0.5)
                set total-mass (structural-mass + (energy-reserve / El))
                set maintenance-energy (energy-reserve * 0.1    )
                set potential-fecundity (std-mass ) * (eggs_per_bass)                      ; maximum number of eggs that can be produced, potential fecundity, depends on body length , changed for bass JW pickett and pawson fecundity estimate.
                set max-R (potential-fecundity * egg_mass * (Ef + Fs))
                calculate-speed
                calculate-r
              ]
            ]
          ]
          ]
          set counter counter + 1
        ]
      ]
    ]
    set loopcount loopcount + 1
  ]

]

end







;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to calc-length-distribution


set L0 mean [l] of turtles with [(age > 0) and (age < 1) and (breed = juv-bass)]
set L1 mean [l] of turtles with [(age > 1) and (age < 2)]
set L2 mean [l] of turtles with [(age > 2) and (age < 3)]
set L3 mean [l] of turtles with [(age > 3) and (age < 4)]
set L4 mean [l] of turtles with [(age > 4) and (age < 5)]
set L5 mean [l] of turtles with [(age > 5) and (age < 6)]
set L6 mean [l] of turtles with [(age > 6) and (age < 7)]
set L7 mean [l] of turtles with [(age > 7) and (age < 8)]
set L8 mean [l] of turtles with [(age > 8) and (age < 9)]
set L9 mean [l] of turtles with [(age > 9) and (age < 10)]
set L10 mean [l] of turtles with [(age > 10) and (age < 11)]
set L11 mean [l] of turtles with [(age > 11) and (age < 12)]
set L12 mean [l] of turtles with [(age > 12) and (age < 13)]
set L13 mean [l] of turtles with [(age > 13) and (age < 14)]
set L14 mean [l] of turtles with [(age > 14) and (age < 15)]
set L15 mean [l] of turtles with [(age > 15) and (age < 16)]
set L16 mean [l] of turtles with [(age > 16) and (age < 17)]
set L17 mean [l] of turtles with [(age > 17) and (age < 18)]
set L18 mean [l] of turtles with [(age > 18) and (age < 19)]
set L19 mean [l] of turtles with [(age > 19) and (age < 20)]
set L20 mean [l] of turtles with [(age > 20) and (age < 21)]
set L21 mean [l] of turtles with [(age > 21) and (age < 22)]
set L22 mean [l] of turtles with [(age > 22) and (age < 23)]
set L23 mean [l] of turtles with [(age > 23) and (age < 24)]
set L24 mean [l] of turtles with [(age > 24) and (age < 25)]
set L25 mean [l] of turtles with [(age > 25) and (age < 26)]
set L26 mean [l] of turtles with [(age > 26) and (age < 27)]
set L27 mean [l] of turtles with [(age > 27) and (age < 28)]
set L28 mean [l] of turtles with [(age > 28) and (age < 29)]
set L29 mean [l] of turtles with [(age > 29) and (age < 30)]


end



to calc-mass-distribution


set M0 mean [total-mass] of turtles with [(age > 0) and (age < 1)and (breed = juv-bass)]
set M1 mean [total-mass] of turtles with [(age > 1) and (age < 2)]
set M2 mean [total-mass] of turtles with [(age > 2) and (age < 3)]
set M3 mean [total-mass] of turtles with [(age > 3) and (age < 4)]
set M4 mean [total-mass] of turtles with [(age > 4) and (age < 5)]
set M5 mean [total-mass] of turtles with [(age > 5) and (age < 6)]
set M6 mean [total-mass] of turtles with [(age > 6) and (age < 7)]
set M7 mean [total-mass] of turtles with [(age > 7) and (age < 8)]
set M8 mean [total-mass] of turtles with [(age > 8) and (age < 9)]
set M9 mean [total-mass] of turtles with [(age > 9) and (age < 10)]
set M10 mean [total-mass] of turtles with [(age > 10) and (age < 11)]
set M11 mean [total-mass] of turtles with [(age > 11) and (age < 12)]
set M12 mean [total-mass] of turtles with [(age > 12) and (age < 13)]
set M13 mean [total-mass] of turtles with [(age > 13) and (age < 14)]
set M14 mean [total-mass] of turtles with [(age > 14) and (age < 15)]
set M15 mean [total-mass] of turtles with [(age > 15) and (age < 16)]
set M16 mean [total-mass] of turtles with [(age > 16) and (age < 17)]
set M17 mean [total-mass] of turtles with [(age > 17) and (age < 18)]
set M18 mean [total-mass] of turtles with [(age > 18) and (age < 19)]
set M19 mean [total-mass] of turtles with [(age > 19) and (age < 20)]
set M20 mean [total-mass] of turtles with [(age > 20) and (age < 21)]
set M21 mean [total-mass] of turtles with [(age > 21) and (age < 22)]
set M22 mean [total-mass] of turtles with [(age > 22) and (age < 23)]
set M23 mean [total-mass] of turtles with [(age > 23) and (age < 24)]
set M24 mean [total-mass] of turtles with [(age > 24) and (age < 25)]
set M25 mean [total-mass] of turtles with [(age > 25) and (age < 26)]
set M26 mean [total-mass] of turtles with [(age > 26) and (age < 27)]
set M27 mean [total-mass] of turtles with [(age > 27) and (age < 28)]
set M28 mean [total-mass] of turtles with [(age > 28) and (age < 29)]
set M29 mean [total-mass] of turtles with [(age > 29) and (age < 30)]

end


to calc-total-mass-distribution


set MT0 sum [total-mass] of turtles with [(age > 0) and (age < 1)and (breed = juv-bass)]
set MT1 sum [total-mass] of turtles with [(age > 1) and (age < 2)]
set MT2 sum [total-mass] of turtles with [(age > 2) and (age < 3)]
set MT3 sum [total-mass] of turtles with [(age > 3) and (age < 4)]
set MT4 sum [total-mass] of turtles with [(age > 4) and (age < 5)]
set MT5 sum [total-mass] of turtles with [(age > 5) and (age < 6)]
set MT6 sum [total-mass] of turtles with [(age > 6) and (age < 7)]
set MT7 sum [total-mass] of turtles with [(age > 7) and (age < 8)]
set MT8 sum [total-mass] of turtles with [(age > 8) and (age < 9)]
set MT9 sum [total-mass] of turtles with [(age > 9) and (age < 10)]
set MT10 sum [total-mass] of turtles with [(age > 10) and (age < 11)]
set MT11 sum [total-mass] of turtles with [(age > 11) and (age < 12)]
set MT12 sum [total-mass] of turtles with [(age > 12) and (age < 13)]
set MT13 sum [total-mass] of turtles with [(age > 13) and (age < 14)]
set MT14 sum [total-mass] of turtles with [(age > 14) and (age < 15)]
set MT15 sum [total-mass] of turtles with [(age > 15) and (age < 16)]
set MT16 sum [total-mass] of turtles with [(age > 16) and (age < 17)]
set MT17 sum [total-mass] of turtles with [(age > 17) and (age < 18)]
set MT18 sum [total-mass] of turtles with [(age > 18) and (age < 19)]
set MT19 sum [total-mass] of turtles with [(age > 19) and (age < 20)]
set MT20 sum [total-mass] of turtles with [(age > 20) and (age < 21)]
set MT21 sum [total-mass] of turtles with [(age > 21) and (age < 22)]
set MT22 sum [total-mass] of turtles with [(age > 22) and (age < 23)]
set MT23 sum [total-mass] of turtles with [(age > 23) and (age < 24)]
set MT24 sum [total-mass] of turtles with [(age > 24) and (age < 25)]
set MT25 sum [total-mass] of turtles with [(age > 25) and (age < 26)]
set MT26 sum [total-mass] of turtles with [(age > 26) and (age < 27)]
set MT27 sum [total-mass] of turtles with [(age > 27) and (age < 28)]
set MT28 sum [total-mass] of turtles with [(age > 28) and (age < 29)]
set MT29 sum [total-mass] of turtles with [(age > 29) and (age < 30)]

end




to calc-age-distribution


set N0 sum [number] of turtles with [(age > 0) and (age < 1)and (breed = juv-bass)]
set N1 sum [number] of turtles with  [(age > 1) and (age < 2)]
set N2 sum [number] of turtles with  [(age > 2) and (age < 3)]
set N3 sum [number] of turtles with  [(age > 3) and (age < 4)]
set N4 sum [number] of turtles with [(age > 4) and (age < 5)]
set N5 sum [number] of turtles with [(age > 5) and (age < 6)]
set N6 sum [number] of turtles with [(age > 6) and (age < 7)]
set N7 sum [number] of turtles with  [(age > 7) and (age < 8)]
set N8 sum [number] of turtles with  [(age > 8) and (age < 9)]
set N9 sum [number] of turtles with [(age > 9) and (age < 10)]
set N10 sum [number] of turtles with  [(age > 10) and (age < 11)]
set N11 sum [number] of turtles with  [(age > 11) and (age < 12)]
set N12 sum [number] of turtles with [(age > 12) and (age < 13)]
set N13 sum [number] of turtles with [(age > 13) and (age < 14)]
set N14 sum [number] of turtles with [(age > 14) and (age < 15)]
set N15 sum [number] of turtles with  [(age > 15) and (age < 16)]
set N16 sum [number] of turtles with  [(age > 16) and (age < 17)]
set N17 sum [number] of turtles with [(age > 17) and (age < 18)]
set N18 sum [number] of turtles with [(age > 18) and (age < 19)]
set N19 sum [number] of turtles with  [(age > 19) and (age < 20)]
set N20 sum [number] of turtles with [(age > 20) and (age < 21)]
set N21 sum [number] of turtles with [(age > 21) and (age < 22)]
set N22 sum [number] of turtles with  [(age > 22) and (age < 23)]
set N23 sum [number] of turtles with  [(age > 23) and (age < 24)]
set N24 sum [number] of turtles with  [(age > 24) and (age < 25)]
set N25 sum [number] of turtles with  [(age > 25) and (age < 26)]
set N26 sum [number] of turtles with  [(age > 26) and (age < 27)]
set N27 sum [number] of turtles with  [(age > 27) and (age < 28)]
set N28 sum [number] of turtles with  [(age > 28) and (age < 29)]
set N29 sum [number] of turtles with  [(age > 29) and (age < 30)]

end












;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
@#$#@#$#@
GRAPHICS-WINDOW
12
86
495
616
-1
-1
13.14
1
10
1
1
1
0
0
0
1
0
35
0
37
0
0
1
ticks
30.0

BUTTON
802
281
865
314
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
800
242
863
275
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

MONITOR
727
241
787
286
year
year
17
1
11

SLIDER
877
298
1049
331
site-fidelity?
site-fidelity?
0
100
100
1
1
NIL
HORIZONTAL

PLOT
514
372
1032
519
plot 1
NIL
NIL
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"eggs" 1.0 0 -16777216 true "" "plot sum [number] of eggs"
"yslarvae" 1.0 0 -11033397 true "" "plot sum [number] of ys-larvae"
"larvae" 1.0 0 -2674135 true "" "plot sum [number] of larvae"
"Juv bass" 1.0 0 -13840069 true "" "plot sum [number] of juv-bass"
"Mat bass" 1.0 0 -4699768 true "" "plot sum [number] of mat-bass"

MONITOR
727
290
799
335
NIL
rolling_tick
17
1
11

PLOT
1072
331
1546
491
Rec
NIL
NIL
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"rec" 1.0 0 -16777216 true "" "plot sum [number] of turtles with [(age < 1) and (breed = juv-bass)]"
"mat" 1.0 0 -2674135 true "" ""

SLIDER
878
223
1051
256
com-net-hole-size
com-net-hole-size
0
60
0
1
1
NIL
HORIZONTAL

SLIDER
877
262
1050
295
angler-min-size
angler-min-size
0
80
0
1
1
NIL
HORIZONTAL

PLOT
1071
178
1543
328
num
NIL
NIL
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"model" 1.0 0 -2674135 true "" "plot (sum [number] of turtles with [(breed = juv-bass) or (breed = mat-bass)] )"
"ICES" 1.0 0 -16777216 true "" "plot icesnum"

PLOT
1074
23
1540
173
age distribution
NIL
NIL
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"ICES" 1.0 0 -16777216 true "" ""
"model" 1.0 0 -2674135 true "" ""

SWITCH
517
95
621
128
spin_up
spin_up
0
1
-1000

BUTTON
796
201
868
234
NIL
spin-up
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
798
160
871
193
NIL
go-ABC
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

@#$#@#$#@
## WHAT IS IT?



## HOW IT WORKS



## HOW TO USE IT




## THINGS TO NOTICE




## THINGS TO TRY


## EXTENDING THE MODEL



## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270

@#$#@#$#@
NetLogo 5.3.1
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="Validation outputs" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>mean [l] of turtles with [(age &gt; 0) and (age &lt; 1)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 1) and (age &lt; 2)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 2) and (age &lt; 3)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 3) and (age &lt; 4)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 4) and (age &lt; 5)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 5) and (age &lt; 6)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 6) and (age &lt; 7)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 7) and (age &lt; 8)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 8) and (age &lt; 9)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 9) and (age &lt; 10)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 10) and (age &lt; 11)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 11) and (age &lt; 12)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 12) and (age &lt; 13)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 13) and (age &lt; 14)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 14) and (age &lt; 15)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 15) and (age &lt; 16)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 16) and (age &lt; 17)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 17) and (age &lt; 18)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 18) and (age &lt; 19)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 19) and (age &lt; 20)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 20) and (age &lt; 21)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 21) and (age &lt; 22)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 22) and (age &lt; 23)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 23) and (age &lt; 24)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 24) and (age &lt; 25)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 25) and (age &lt; 26)]</metric>
    <metric>sum [number] of turtles with [(age &gt; 0) and (age &lt; 1)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 1) and (age &lt; 2)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 2) and (age &lt; 3)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 3) and (age &lt; 4)]</metric>
    <metric>sum [number] of turtles with [(age &gt; 4) and (age &lt; 5)]</metric>
    <metric>sum [number] of turtles with [(age &gt; 5) and (age &lt; 6)]</metric>
    <metric>sum [number] of turtles with [(age &gt; 6) and (age &lt; 7)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 7) and (age &lt; 8)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 8) and (age &lt; 9)]</metric>
    <metric>sum [number] of turtles with [(age &gt; 9) and (age &lt; 10)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 10) and (age &lt; 11)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 11) and (age &lt; 12)]</metric>
    <metric>sum [number] of turtles with [(age &gt; 12) and (age &lt; 13)]</metric>
    <metric>sum [number] of turtles with [(age &gt; 13) and (age &lt; 14)]</metric>
    <metric>sum [number] of turtles with [(age &gt; 14) and (age &lt; 15)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 15) and (age &lt; 16)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 16) and (age &lt; 17)]</metric>
    <metric>sum [number] of turtles with [(age &gt; 17) and (age &lt; 18)]</metric>
    <metric>sum [number] of turtles with [(age &gt; 18) and (age &lt; 19)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 19) and (age &lt; 20)]</metric>
    <metric>sum [number] of turtles with [(age &gt; 20) and (age &lt; 21)]</metric>
    <metric>sum [number] of turtles with [(age &gt; 21) and (age &lt; 22)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 22) and (age &lt; 23)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 23) and (age &lt; 24)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 24) and (age &lt; 25)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 25) and (age &lt; 26)]</metric>
    <metric>sum [number] of turtles with [(age &lt; 1) and (breed = juv-bass)]</metric>
    <metric>annual_ssb</metric>
    <metric>year</metric>
    <enumeratedValueSet variable="Mc">
      <value value="9.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="absorbed-energy">
      <value value="0.9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Me">
      <value value="0.03712"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="i">
      <value value="55000000000000"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="sensitivity" repetitions="5" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <exitCondition>year = 2005</exitCondition>
    <metric>mean [l] of turtles with [(age &gt; 0) and (age &lt; 1)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 1) and (age &lt; 2)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 2) and (age &lt; 3)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 3) and (age &lt; 4)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 4) and (age &lt; 5)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 5) and (age &lt; 6)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 6) and (age &lt; 7)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 7) and (age &lt; 8)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 8) and (age &lt; 9)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 9) and (age &lt; 10)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 10) and (age &lt; 11)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 11) and (age &lt; 12)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 12) and (age &lt; 13)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 13) and (age &lt; 14)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 14) and (age &lt; 15)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 15) and (age &lt; 16)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 16) and (age &lt; 17)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 17) and (age &lt; 18)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 18) and (age &lt; 19)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 19) and (age &lt; 20)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 20) and (age &lt; 21)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 21) and (age &lt; 22)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 22) and (age &lt; 23)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 23) and (age &lt; 24)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 24) and (age &lt; 25)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 25) and (age &lt; 26)]</metric>
    <metric>sum [number] of turtles with [(age &gt; 0) and (age &lt; 1)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 1) and (age &lt; 2)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 2) and (age &lt; 3)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 3) and (age &lt; 4)]</metric>
    <metric>sum [number] of turtles with [(age &gt; 4) and (age &lt; 5)]</metric>
    <metric>sum [number] of turtles with [(age &gt; 5) and (age &lt; 6)]</metric>
    <metric>sum [number] of turtles with [(age &gt; 6) and (age &lt; 7)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 7) and (age &lt; 8)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 8) and (age &lt; 9)]</metric>
    <metric>sum [number] of turtles with [(age &gt; 9) and (age &lt; 10)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 10) and (age &lt; 11)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 11) and (age &lt; 12)]</metric>
    <metric>sum [number] of turtles with [(age &gt; 12) and (age &lt; 13)]</metric>
    <metric>sum [number] of turtles with [(age &gt; 13) and (age &lt; 14)]</metric>
    <metric>sum [number] of turtles with [(age &gt; 14) and (age &lt; 15)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 15) and (age &lt; 16)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 16) and (age &lt; 17)]</metric>
    <metric>sum [number] of turtles with [(age &gt; 17) and (age &lt; 18)]</metric>
    <metric>sum [number] of turtles with [(age &gt; 18) and (age &lt; 19)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 19) and (age &lt; 20)]</metric>
    <metric>sum [number] of turtles with [(age &gt; 20) and (age &lt; 21)]</metric>
    <metric>sum [number] of turtles with [(age &gt; 21) and (age &lt; 22)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 22) and (age &lt; 23)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 23) and (age &lt; 24)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 24) and (age &lt; 25)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 25) and (age &lt; 26)]</metric>
    <metric>sum [number] of turtles with [(age &lt; 1) and (breed = juv-bass)]</metric>
    <metric>annual_ssb</metric>
    <metric>year</metric>
    <enumeratedValueSet variable="linf">
      <value value="84.55"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k">
      <value value="0.09699"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t0">
      <value value="-0.73"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="refTS">
      <value value="279"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="A">
      <value value="1.76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Ea">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EaS">
      <value value="0.1903656"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Tref">
      <value value="285.15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boltz">
      <value value="9.482E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Cmax">
      <value value="0.54"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ep">
      <value value="6.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="A0">
      <value value="0.1227808"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Ef">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="El">
      <value value="39.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Ls">
      <value value="14.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fs">
      <value value="3.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="egg-mass">
      <value value="0.0096"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ncohort">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="a_g">
      <value value="1.2312E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="b_g">
      <value value="2.969"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="eggs_per_bass">
      <value value="375000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Gl">
      <value value="0.02485"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AM">
      <value value="6.6E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AE">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="PM">
      <value value="0.03712"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="H">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="I">
      <value value="50000000000000"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Outputs" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>mean [l] of turtles with [(age &gt; 0) and (age &lt; 1)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 1) and (age &lt; 2)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 2) and (age &lt; 3)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 3) and (age &lt; 4)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 4) and (age &lt; 5)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 5) and (age &lt; 6)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 6) and (age &lt; 7)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 7) and (age &lt; 8)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 8) and (age &lt; 9)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 9) and (age &lt; 10)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 10) and (age &lt; 11)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 11) and (age &lt; 12)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 12) and (age &lt; 13)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 13) and (age &lt; 14)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 14) and (age &lt; 15)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 15) and (age &lt; 16)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 16) and (age &lt; 17)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 17) and (age &lt; 18)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 18) and (age &lt; 19)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 19) and (age &lt; 20)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 20) and (age &lt; 21)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 21) and (age &lt; 22)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 22) and (age &lt; 23)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 23) and (age &lt; 24)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 24) and (age &lt; 25)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 25) and (age &lt; 26)]</metric>
    <metric>sum [number] of turtles with [(age &gt; 0) and (age &lt; 1)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 1) and (age &lt; 2)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 2) and (age &lt; 3)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 3) and (age &lt; 4)]</metric>
    <metric>sum [number] of turtles with [(age &gt; 4) and (age &lt; 5)]</metric>
    <metric>sum [number] of turtles with [(age &gt; 5) and (age &lt; 6)]</metric>
    <metric>sum [number] of turtles with [(age &gt; 6) and (age &lt; 7)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 7) and (age &lt; 8)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 8) and (age &lt; 9)]</metric>
    <metric>sum [number] of turtles with [(age &gt; 9) and (age &lt; 10)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 10) and (age &lt; 11)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 11) and (age &lt; 12)]</metric>
    <metric>sum [number] of turtles with [(age &gt; 12) and (age &lt; 13)]</metric>
    <metric>sum [number] of turtles with [(age &gt; 13) and (age &lt; 14)]</metric>
    <metric>sum [number] of turtles with [(age &gt; 14) and (age &lt; 15)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 15) and (age &lt; 16)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 16) and (age &lt; 17)]</metric>
    <metric>sum [number] of turtles with [(age &gt; 17) and (age &lt; 18)]</metric>
    <metric>sum [number] of turtles with [(age &gt; 18) and (age &lt; 19)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 19) and (age &lt; 20)]</metric>
    <metric>sum [number] of turtles with [(age &gt; 20) and (age &lt; 21)]</metric>
    <metric>sum [number] of turtles with [(age &gt; 21) and (age &lt; 22)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 22) and (age &lt; 23)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 23) and (age &lt; 24)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 24) and (age &lt; 25)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 25) and (age &lt; 26)]</metric>
    <metric>sum [number] of turtles with [(age &lt; 1) and (breed = juv-bass)]</metric>
    <metric>annual_ssb</metric>
    <metric>year</metric>
    <enumeratedValueSet variable="linf">
      <value value="84.55"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k">
      <value value="0.09699"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t0">
      <value value="-0.73"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="refTS">
      <value value="279"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="A">
      <value value="1.76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Ea">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EaS">
      <value value="0.1903656"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Tref">
      <value value="285.15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boltz">
      <value value="8.62E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Cmax">
      <value value="0.54"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ep">
      <value value="6.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="A0">
      <value value="0.1227808"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Ef">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="El">
      <value value="39.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Ls">
      <value value="14.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fs">
      <value value="3.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="egg-mass">
      <value value="0.0096"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ncohort">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="a_g">
      <value value="1.2312E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="b_g">
      <value value="2.969"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="eggs_per_bass">
      <value value="375000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Gl">
      <value value="0.02485"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AM">
      <value value="6.6E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AE">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="PM">
      <value value="0.03712"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="H">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="I">
      <value value="50000000000000"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="check en outputs" repetitions="5" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <exitCondition>year = 2005</exitCondition>
    <metric>mean [l] of turtles with [(age &gt; 0) and (age &lt; 1)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 1) and (age &lt; 2)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 2) and (age &lt; 3)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 3) and (age &lt; 4)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 4) and (age &lt; 5)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 5) and (age &lt; 6)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 6) and (age &lt; 7)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 7) and (age &lt; 8)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 8) and (age &lt; 9)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 9) and (age &lt; 10)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 10) and (age &lt; 11)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 11) and (age &lt; 12)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 12) and (age &lt; 13)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 13) and (age &lt; 14)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 14) and (age &lt; 15)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 15) and (age &lt; 16)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 16) and (age &lt; 17)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 17) and (age &lt; 18)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 18) and (age &lt; 19)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 19) and (age &lt; 20)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 20) and (age &lt; 21)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 21) and (age &lt; 22)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 22) and (age &lt; 23)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 23) and (age &lt; 24)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 24) and (age &lt; 25)]</metric>
    <metric>mean [l] of turtles with [(age &gt; 25) and (age &lt; 26)]</metric>
    <metric>sum [number] of turtles with [(age &gt; 0) and (age &lt; 1)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 1) and (age &lt; 2)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 2) and (age &lt; 3)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 3) and (age &lt; 4)]</metric>
    <metric>sum [number] of turtles with [(age &gt; 4) and (age &lt; 5)]</metric>
    <metric>sum [number] of turtles with [(age &gt; 5) and (age &lt; 6)]</metric>
    <metric>sum [number] of turtles with [(age &gt; 6) and (age &lt; 7)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 7) and (age &lt; 8)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 8) and (age &lt; 9)]</metric>
    <metric>sum [number] of turtles with [(age &gt; 9) and (age &lt; 10)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 10) and (age &lt; 11)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 11) and (age &lt; 12)]</metric>
    <metric>sum [number] of turtles with [(age &gt; 12) and (age &lt; 13)]</metric>
    <metric>sum [number] of turtles with [(age &gt; 13) and (age &lt; 14)]</metric>
    <metric>sum [number] of turtles with [(age &gt; 14) and (age &lt; 15)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 15) and (age &lt; 16)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 16) and (age &lt; 17)]</metric>
    <metric>sum [number] of turtles with [(age &gt; 17) and (age &lt; 18)]</metric>
    <metric>sum [number] of turtles with [(age &gt; 18) and (age &lt; 19)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 19) and (age &lt; 20)]</metric>
    <metric>sum [number] of turtles with [(age &gt; 20) and (age &lt; 21)]</metric>
    <metric>sum [number] of turtles with [(age &gt; 21) and (age &lt; 22)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 22) and (age &lt; 23)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 23) and (age &lt; 24)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 24) and (age &lt; 25)]</metric>
    <metric>sum [number] of turtles with  [(age &gt; 25) and (age &lt; 26)]</metric>
    <metric>sum [number] of turtles with [(age &lt; 1) and (breed = juv-bass)]</metric>
    <metric>annual_ssb</metric>
    <metric>year</metric>
    <enumeratedValueSet variable="linf">
      <value value="84.55"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k">
      <value value="0.09699"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="t0">
      <value value="-0.73"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="refTS">
      <value value="279"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="A">
      <value value="1.76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Ea">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EaS">
      <value value="0.1903656"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Tref">
      <value value="285.15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boltz">
      <value value="8.62E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Cmax">
      <value value="0.54"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ep">
      <value value="6.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="A0">
      <value value="0.1227808"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Ef">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="El">
      <value value="39.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Ls">
      <value value="14.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fs">
      <value value="3.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="egg-mass">
      <value value="0.0096"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ncohort">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="a_g">
      <value value="1.2312E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="b_g">
      <value value="2.969"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="eggs_per_bass">
      <value value="375000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Gl">
      <value value="0.02485"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AM">
      <value value="6.6E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AE">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="PM">
      <value value="0.03712"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="H">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="I">
      <value value="55000000000000"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180

@#$#@#$#@
0
@#$#@#$#@
