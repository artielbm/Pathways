# Type 2 Diabetes Pathways
## Matlab and xpp files for Pathways to Diabetes Model

This is supplementary material for:

J. Ha and A. Sherman. 2019. Type 2 Diabetes: One Disease, Many Pathways
doi: https://doi.org/10.1101/648816.

Matlab equation definition file for all figures:
[pathway.m](./pathway_matlab/pathway.m)

Matlab run files for individual figures:
Figures 1 - 4. Matlab will generate Fig1 - Fig6, corresponding to Panels A - F (a counter will tick off 365 iterations):

Figure 1: [FIG1.m](./pathway_matlab/FIG1.m)

Figure 2

Figure 3

Figure 4

Figures 6, 7. Matlab will generate Panels A, B (IGT first case):

Figure 6
Figure 7
XPP source file for all figures:
pathway.ode

XPP set files for individual figures:
Figures 1 - 4. xpp will plot the longitudinal response of glucose to meals each day. The envelope of the plot will show the fasting and peak postprandial values. Use the xpp menus to display the other variables:

FIG1.set (IGT first pathway)
FIG2.set (Sustained IGT )
FIG3.set (IFG first pathway IGT )
FIG4.set (IGT without CGI)
