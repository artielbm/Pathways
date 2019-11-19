# Type 2 Diabetes Pathways
## Matlab and xpp files for Pathways to Diabetes Model

This is supplementary material for:

J. Ha and A. Sherman. 2019. Type 2 Diabetes: One Disease, Many Pathways
doi: https://doi.org/10.1101/648816.

### Matlab equation definition file for all figures:

[pathway.m](./pathway_matlab/pathway.m)

### Matlab run files for individual figures:

#### Figures 1 - 4. Matlab will generate Fig1 - Fig6, corresponding to Panels A - F (a counter will tick off 365 iterations):

Figure 1: [FIG1.m](./pathway_matlab/FIG1.m)

Figure 2: [FIG2.m](./pathway_matlab/FIG2.m)

Figure 3: [FIG3.m](./pathway_matlab/FIG3.m)

Figure 4: [FIG4.m](./pathway_matlab/FIG4.m)

#### Figures 6, 7. Matlab will generate Panels A, B (IGT first case):

Figure 6: [FIG6.m](./pathway_matlab/FIG6.m)

Figure 7: [FIG7.m](./pathway_matlab/FIG7.m)

### XPP source file for all figures:

pathway.ode: [pathway.ode](./pathway_xpp/pathway.ode)

### XPP set files for individual figures:

#### Figures 1 - 4. xpp will plot the longitudinal response of glucose to meals each day. The envelope of the plot will show the fasting and peak postprandial values. Use the xpp menus to display the other variables:

FIG1.set (IGT first pathway): [FIG1.set](./pathway_xpp/FIG1.set)

FIG2.set (Sustained IGT): [FIG2.set](./pathway_xpp/FIG2.set)

FIG3.set (IFG first pathway IGT): [FIG3.set](./pathway_xpp/FIG3.set)

FIG4.set (IGT without CGI): [FIG4.set](./pathway_xpp/FIG4.set)
