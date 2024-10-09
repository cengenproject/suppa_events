# Scripts

First call `src/suppa2_generate_events.sh` to create the ioe files with events descriptions. Then, as indicated in a comment in `src/suppa2_psiPerEvent.sh`, group these ioe for each event type in a single one.

Separately, bring in quantifications from stringtie (see `stringtie_quantifs` repo). If needed, prefilter by removing outliers with `R/str_q_remove_outliers.R`

Run `src/suppa2_psiPerEvent.sh`.

Run `src/suppa2_dpsi.sh`, which will call `src/split_events.R` to group PSI by neuron type.


# Analysis

Main analysis in `R/dpsi.R`. This script uses accessory functions defined in `R/sequence_properties.R`, `R/sequence_conservation.R`, and `R/extract_event_coordinates.R`.

Additional analyses directly on the PSIs (as opposed to deltaPSIs) are in `R/psi.R`, not included in paper.


# Older versions

I initially thought about using Salmon for quantifications: the fact that we need to deduplicate UMIs makes it non-trivial, deleted (can be found in old commits).



