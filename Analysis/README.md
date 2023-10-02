# Analysis 

* ### exampleFVAsimilarity.m

The code is used to calculate the similarity between 2 FVA calculated on the same model for different conditions (e.g. healthy vs pathological, test vs reference treatment). 

Similarity is calculated using the `FVAsimilarity` function. The function output contains: 
- `overallSim` a variable between 0 and 1, where 1 indicates identical FVA between the two models
- `rxnSim` a arrat containing values between 0 and 1 for each reactions, indicating the similarity for each reaction between the two conditions 

* ### manualGapfilling.m 

We are looking to improve the Arabidopsis model () to enable it to produce Luteine and S_epsilon_Carotene. To acheive this, we manually incoporating and eliminating reactions within the model and assissing their impact on the production of the targeted molecules. 

* ### Automated gap filling:

See COBRA toolbox tutorial:
https://cobrapy.readthedocs.io/en/latest/gapfilling.html
