# Thresholding for Predicted Synergistic Pairs

This project thresholds highly predicted (>0.7) ENCODE-rE2G Extended enhancer-gene pairs in K562 based on inter-enhancer-contact, EP300 values, and distance, and can be easily applied to other cell types.

Resources (K562 Hi-C, K562 ENCODE Predictions) are linked in the config.yaml file along with the thresholds used for each parameter.

All thresholds were determined from a screen which identified synergistic EE pairs. These thresholds were used to nominate ENCODE predicted enhancers to further test for synergistic interactions.

The current workflow is as follows:
1. Filter for all EG pairs with a Score > 0.7

2. Filter for all EG pairs where the distance between the center of the enhancer its predicted gene's TSS is < 2Mb

3. For every gene G, gather all predicted enhancers. Separately for downstream and upstream enhancers, create all possible enhancer-enhancer EE pairs.

E.g. If the following enhancers (after filters 1 and 2) predicted for gene G are:

upstream |---E1-------E2----E3--G----------E4-----E5--| downstream

Then step three will result in the following EE pairs: E1 E2; E1 E3; E2 E3; E5 E4. The first element in each EE pair will always be farther from the gene TSS than the second.

4. For each EE pair, calculate the contact using Hi-C with the "observed" counts using "SCALE" normalization. Resolution can be found in the config.yaml. This is done separately for each chromosome to speed up the workflow. 

5. Pairs are then filtered such that the element farther from the TSS (element 1) has a EP300 overlap value larger than the value indicated by 'p300-threshold' in the config.yaml, and the closer element (element 2) has a value below that threshold.

6. EE pairs are nominated as possibly synergistic if they are less than the EE-distance-threshold OR they are greater than the EE-distance-threshold and the contact between the EE pair is greater than the EE-contact-threshold.
