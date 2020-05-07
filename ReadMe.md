# Mark-recapture modeling for mussel populations in Texas.

Completed dissertation research of V. A. Sotola, PhD, Texas State University.
Submitted to <em>Freshwater Biology</em>.


## Data files
Surveys included three secondary sampling dates within each of 5 primary sampling events at **upper** and **lower** sites in the river.


### Upper Site
<b>`data_SITE_tags.csv`</b> Floy and PIT tag IDs associated with individual mussels at upper site.
<b>`data_SITE_floy.csv`</b> Floy tag detections at upper site. 
<b>`data_SITE_pit.csv`</b>  PIT tag detections at upper site. 


### Lower Site
<b>`data_alt_tags.csv`</b> Floy and PIT tag IDs associated with individual mussels at lower site.

<b>`data_alt_floy.csv`</b> Floy tag detections at upper site keyed to `floy_id` in `data_alt_tags.csv`. 

<b>`data_alt_pit.csv`</b>  PIT tag detections at upper site keyed to `floy_id` in `data_alt_tags.csv`. 


## Model files
Robust-design mark-recapture models used in this study, modified from  Riecke et al. (2018) 

<b>`FWB-P-Feb-20-0083-R1.R`</b> Contains data manipulation, JAGS model, and plotting code.
