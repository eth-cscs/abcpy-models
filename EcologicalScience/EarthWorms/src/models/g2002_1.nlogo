;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;                                   Definition of Model Parameters                                   ;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

globals [ empirical_data food_per_patch half_per_patch speed_in_patches reference_T Boltz T
  Arrhenius adult_size juvenile_size cocoon_size ]

patches-own 
[ change_in_food                    ; change in food per timestep
  current_food                      ; food (g)
  func_response                     ; relative functional response
  ]

turtles-own 
[ mass                              ; individual mass (g)
  BMR                               ; energy cost of maintenance (kJ)
  ingestion_rate                    ; ingestion rate (g/day)
  
  energy_left                       ; intermediate daily energy reserve (kJ)
  energy_reserve_max                ; maximum energy reserve (kJ) (dependent on mass)
  energy_reserve                    ; energy reserve (kJ): amount of energy stored as tissue (7 kJ/g)
  hatchlings                        ; number of hatchlings produced
  R                                 ; energy accumulated towards cocoon production (kJ)
  ]

breed [ adults adult ]
breed [ juveniles juvenile ]
breed [ cocoons cocoon ]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;                                          Setup Functions                                           ;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup-interface
  set initial_number_juveniles 5    ; number of juveniles to start the simulation with
  set initial_number_adults 0       ; number of adults to start the simulation with
  set total_food 150                ; total amount of food initially available (g)
  set scape_size 0.0144             ; size of the scape (m)
  set temperature 20                ; temperature (C)
end

to setup-parameters
  set B_0 967                       ; kJ/day
  set activation_energy 0.25        ; eV
  set energy_tissue 7               ; kJ/g
  set energy_food 10.6              ; kJ/g
  
  set energy_synthesis 3.6          ; kJ/g  
  set half_saturation_coef 3.5      ; g
  set max_ingestion_rate 0.7        ; g/day 
  set mass_birth 0.011              ; g
  set mass_cocoon 0.015             ; g
  
  set mass_maximum  0.5             ; g
  set mass_sexual_maturity 0.25     ; g
  set growth_constant  0.177        ; g
  set max_reproduction_rate 0.182   ; kJ/g/day
  set speed 0.004                   ; m/day
end

to setup
  clear-all
  
  ;; from Gunadi et al. (2002)
  set empirical_data [ 0.011 0.032 0.079 0.131 0.256 0.368 0.414 0.442 0.462 0.448 
    0.45 0.459 0.494 0.495  0.49 0.482 0.447 0.439 0.422 0.412 0.4 0.402 0.387 0.402 0.393 0.341 0.333 ]

  set-default-shape adults "worm"
  set-default-shape juveniles "worm"
  set-default-shape cocoons "dot"
  
  set adult_size sqrt(count patches) * 0.3
  set juvenile_size sqrt(count patches) * 0.2
  set cocoon_size sqrt(count patches) * 0.1
  
  set food_per_patch total_food / count patches
  set half_per_patch ((half_saturation_coef * scape_size) / 0.01) / count patches
  set speed_in_patches speed / (sqrt(scape_size) / sqrt(count patches))
  
  set reference_T 298.15            ; Kelvins
  set Boltz (8.62 * (10 ^ -5))      ; eV K-1
  
  set T 273.15 + temperature
  set Arrhenius (e ^ ((- activation_energy / Boltz ) * ((1 /  T ) - (1 / reference_T))))
  
  setup-patches
  setup-turtles
  
  reset-ticks
end

to setup-patches
  ask patches 
  [ set current_food food_per_patch
    update-patch ]
end

to setup-turtles
  create-adults initial_number_adults
  [ set color red
    set size adult_size
    setxy random-xcor random-ycor
    set mass mass_sexual_maturity
    set energy_reserve_max (mass / 2) * energy_tissue
    set energy_reserve energy_reserve_max ]
  
  create-juveniles initial_number_juveniles
  [ set color pink
    set size juvenile_size
    setxy random-xcor random-ycor
    set mass mass_birth
    set energy_reserve_max (mass / 2) * energy_tissue
    set energy_reserve energy_reserve_max ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;                                            Run Functions                                           ;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to go
  if not any? turtles with [ breed != cocoons ] or ticks > 182 [ stop ]
  update-experiment
  calc-ingestion
  
  ask adults
  [ calc-assimilation
    calc-maintenance
    calc-reproduction
    calc-growth
    update-reserves
    move ]
  
  ask juveniles
  [ calc-assimilation
    calc-maintenance
    calc-growth
    transform-juvenile
    update-reserves
    move ]
  
  ask patches
  [ update-patch ]
  
  tick
end

to go-7 
  repeat 7 [ go ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;                                      Ingestion & Assimilation                                      ;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to calc-ingestion
  ask turtles with [ breed != cocoons ] [ set ingestion_rate (max_ingestion_rate * Arrhenius) *
    func_response * (mass ^ (2 / 3)) ]
  ask patches [ correct-ingestion-rate ]
end

to calc-assimilation
  set energy_left (ingestion_rate * energy_food)
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;                                      Maintenance & Starvation                                      ;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to calc-maintenance
  set BMR B_0 * (mass ^ (3 / 4)) * e ^ (- activation_energy / (Boltz * T))
  
  ifelse energy_left > BMR
  [ set energy_left energy_left - BMR ]
  [ set energy_reserve energy_reserve - (BMR - energy_left)
    set energy_left 0 ]
  
  if energy_reserve < (energy_reserve_max * 0.5) and breed != cocoons
  [ onset-starvation-strategy ]
end

to onset-starvation-strategy
  set mass mass - (BMR / (energy_tissue + energy_synthesis))
  set energy_reserve energy_reserve + BMR
  
  if mass < mass_sexual_maturity 
  [ set breed juveniles
    set color pink
    set size juvenile_size ]
  
  if mass < mass_birth 
  [ die ]  
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;                                            Reproduction                                            ;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to calc-reproduction
  let max_R (max_reproduction_rate * Arrhenius) * mass
  
  let energy_for_R min list energy_left max_R
  set energy_left energy_left - energy_for_R
  
  if energy_for_R < max_R and energy_reserve > 0.5 * energy_reserve_max
  [ set energy_reserve energy_reserve - (max_R - energy_for_R)
    set energy_for_R max_R ] 
  
  set R R + energy_for_R
  
  if R >= (mass_cocoon * (energy_tissue + energy_synthesis))
  [ reproduce ]
end

to reproduce
  hatch-cocoons 1
  [ set color white
    set size cocoon_size
    set mass mass_cocoon
    set ingestion_rate 0
    set energy_left 0
    set energy_reserve_max 0
    set energy_reserve 0
    set hatchlings 0
    set R 0 ]
  
  set hatchlings hatchlings + 1
  set R (R - (mass_cocoon * (energy_tissue + energy_synthesis))) 
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;                                                Growth                                              ;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to calc-growth
  let max_G (growth_constant * Arrhenius) * (mass_maximum ^ (1 / 3) * mass ^ (2 / 3) - mass)
  let energy_for_G min list (max_G * (energy_tissue + energy_synthesis)) energy_left
  let to_grow 0
  
  if max_G > 0 and (energy_tissue + energy_synthesis) > 0
  [ set to_grow (energy_for_G / (max_G * (energy_tissue + energy_synthesis))) * max_G ]
  
  if (mass + to_grow) < mass_maximum
  [ set mass mass + to_grow
    set energy_left energy_left - energy_for_G ]
end

to transform-juvenile
  if mass >= mass_sexual_maturity 
  [ set breed adults
    set color red
    set size adult_size ]
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;                                           Energy Reserves                                          ;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to update-reserves
  set energy_reserve_max (mass / 2) * energy_tissue
  
  if energy_left > 0
  [ ifelse (energy_tissue + energy_synthesis > 0)
    [ set energy_reserve energy_reserve + (energy_left *
        (energy_tissue / (energy_tissue + energy_synthesis))) ]
    [ die ] ]
  
  if energy_reserve > energy_reserve_max and breed != cocoons 
  [ set energy_reserve energy_reserve_max ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;                                              Movement                                              ;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to move   
    rt random 90
    lt random 90
    fd speed_in_patches
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;                                           Patch Functions                                          ;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to update-patch
  if current_food > 0 
  [ set current_food max list (current_food - sum [ ingestion_rate ] of
      turtles-here with [ breed != cocoons ]) 0
    set func_response current_food / (current_food + half_per_patch) ]

  ifelse current_food > 0 
  [ set pcolor scale-color green current_food (food_per_patch * 2) (0 - food_per_patch * 0.5) ]
  [ set pcolor brown ]
end

to correct-ingestion-rate
  if current_food = 0
  [ ask turtles-here with [ breed != cocoons ] [ set ingestion_rate 0 ] ]
  
  if sum [ ingestion_rate ] of turtles-here with [ breed != cocoons ] > current_food
  [ ask turtles-here with [ breed != cocoons ] [ set ingestion_rate current_food /
      count turtles-here with [ breed != cocoons ] ]
    set current_food 0 ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;                                        Experiment Simulation                                       ;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to update-experiment
  if remainder ticks 7 = 0
  [ ask cocoons [ die ]
    ask turtles [ setxy random-xcor random-ycor ] ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;                                             Reporters                                              ;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report mean-mass
  ifelse any? turtles with [ breed != cocoons ]
  [ report mean [ mass ] of turtles with [ breed != cocoons ] ]
  [ report 0 ]
end

to-report mean-hatchlings
  ifelse any? turtles with [ breed != cocoons ]
  [ report mean [ hatchlings ] of turtles with [ breed != cocoons ] ]
  [ report 0 ]
end
@#$#@#$#@
GRAPHICS-WINDOW
267
21
577
352
-1
-1
150.0
1
10
1
1
1
0
1
1
1
0
1
0
1
1
1
1
ticks
30.0

BUTTON
603
220
667
253
NIL
Setup\n
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
813
221
876
254
NIL
Go\n
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
33
22
245
55
initial_number_juveniles
initial_number_juveniles
0
10
5
1
1
NIL
HORIZONTAL

SLIDER
33
170
245
203
temperature
temperature
0
30
20
1
1
NIL
HORIZONTAL

PLOT
601
21
852
201
Food Available
time (ticks)
total food (g)
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot sum [ current_food ] of patches"

INPUTBOX
33
102
131
162
total_food
150
1
0
Number

SLIDER
33
61
245
94
initial_number_adults
initial_number_adults
0
10
0
1
1
NIL
HORIZONTAL

PLOT
866
21
1213
203
Average Body Mass
time (ticks)
mass (g)
0.0
2.0
0.0
0.55
true
true
"" "if ((remainder ticks 7) = 0)\n[ set-current-plot-pen \"model  \"\n  plot mean-mass\n  set-current-plot-pen \"data  \"\n  plot item (ticks / 7) empirical_data ]"
PENS
"data  " 1.0 0 -16777216 true "" ""
"model  " 1.0 0 -7500403 true "" ""

BUTTON
744
221
807
254
Go 7
go-7
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
674
221
737
254
Go
go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
600
268
662
328
B_0
967
1
0
Number

INPUTBOX
669
268
779
328
activation_energy
0.25
1
0
Number

INPUTBOX
786
268
870
328
energy_tissue
7
1
0
Number

INPUTBOX
996
270
1098
330
energy_synthesis
3.6
1
0
Number

INPUTBOX
726
336
845
396
max_ingestion_rate
0.7
1
0
Number

INPUTBOX
600
336
719
396
half_saturation_coef
3.5
1
0
Number

INPUTBOX
852
337
924
397
mass_birth
0.011
1
0
Number

INPUTBOX
703
405
833
465
mass_sexual_maturity
0.25
1
0
Number

INPUTBOX
599
404
696
464
mass_maximum
0.5
1
0
Number

INPUTBOX
930
338
1012
398
mass_cocoon
0.015
1
0
Number

INPUTBOX
842
406
941
466
growth_constant
0.177
1
0
Number

INPUTBOX
950
405
1084
465
max_reproduction_rate
0.182
1
0
Number

BUTTON
880
221
1037
254
Set Basic Parameters
setup-parameters
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
1041
221
1183
254
Set Basic Interface
setup-interface
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
879
269
988
329
energy_food
10.6
1
0
Number

INPUTBOX
1092
406
1186
466
speed
0.0040
1
0
Number

INPUTBOX
139
102
246
162
scape_size
0.0144
1
0
Number

@#$#@#$#@
## CREDITS AND REFERENCES

Copyright (C) 2015 Elske van der Vaart, Mark Beaumont, Alice Johnston, Richard Sibly <elskevdv at gmail.com>

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

This program accompanies the following paper:
van der Vaart, E., Johnston, A.S.A., Beaumont, M.A., & Sibly, R.M. "Calibration and evaluation of individual-based models using Approximate Bayesian Computation" (2015) Ecological Modelling.

The earthworm IBM it runs was originally described here:
Johnston, A.S.A., Hodson, M.E., Thorbek, P., Alvarez, T. & Sibly, R.M. "An energy budget agent-based model of earthworm populations and its application to study the effects of pesticides" (2014) Ecological Modelling, 280, 5 - 17.

The empirical data it plots was taken from:
Gunadi, B., Blount, C. & Edwards, C.A. "The growth and fecundity of Eisenia fetida (Savigny) in cattle solids pre-composted for different periods" (2002)
Pedobiologia, 46, 15 - 23.
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
0
Rectangle -7500403 true true 151 225 180 285
Rectangle -7500403 true true 47 225 75 285
Rectangle -7500403 true true 15 75 210 225
Circle -7500403 true true 135 75 150
Circle -16777216 true false 165 76 116

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

worm
true
0
Polygon -7500403 true true 165 210 165 225 135 255 105 270 90 270 75 255 75 240 90 210 120 195 135 165 165 135 165 105 150 75 150 60 135 60 120 45 120 30 135 15 150 15 180 30 180 45 195 45 210 60 225 105 225 135 210 150 210 165 195 195 180 210

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270

@#$#@#$#@
NetLogo 5.1.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="Temp16.5" repetitions="50" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <exitCondition>ticks = 365</exitCondition>
    <metric>count juveniles + count adults + count cocoons</metric>
    <metric>count adults</metric>
    <metric>sum [hatchlings] of turtles with [breed = adults]</metric>
    <metric>count juveniles</metric>
    <metric>sum [mass] of turtles with [breed = juveniles]</metric>
  </experiment>
  <experiment name="getal:random" repetitions="100" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="189"/>
    <exitCondition>ticks = 189</exitCondition>
    <metric>sum [mass] of turtles with [breed = juveniles]</metric>
    <metric>sum [mass] of turtles with [breed = adults]</metric>
  </experiment>
  <experiment name="experiment" repetitions="100" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="180"/>
    <exitCondition>ticks = 180</exitCondition>
    <metric>sum [hatchlings] of turtles with [breed = adults]</metric>
  </experiment>
  <experiment name="experiment" repetitions="25" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <exitCondition>ticks = 189</exitCondition>
    <metric>sum [mass] of turtles</metric>
    <enumeratedValueSet variable="Initial_number_cocoons">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Initial_number_adults">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food_dynamics">
      <value value="&quot;depleting&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="soil_moisture">
      <value value="80"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food_density_patch">
      <value value="150"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Initial_number_juveniles">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Temperature">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="temp?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="moisture?">
      <value value="false"/>
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
