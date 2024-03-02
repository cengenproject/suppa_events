# Scripts

First call `suppa2_generate_events.sh` to create the ioe files with events descriptions. Then, as indicated in a comment in `suppa2_psiPerEvent.sh`, group these ioe for each event type in a single one.

Separately, bring in quantifications from stringtie (see `stringtie_quantifs`).

Run `suppa2_psiPerEvent.sh`.

Run `suppa2_dpsi.sh`, which will call `split_events.R` to group PSI by neuron type.


# Analysis

See in `R/`


# Older versions

I initially thought about using Salmon for quantifications: the fact that we need to deduplicate UMIs makes it non-trivial, deleted (can be found in old commits).



