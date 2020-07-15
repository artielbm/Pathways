# Type 2 Diabetes Pathways
## Matlab and xpp files for Pathways to Diabetes Model

This is supplementary material for:

J. Ha and A. Sherman. 2019. Type 2 Diabetes: One Disease, Many Pathways [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/32663101/).

### Matlab equation definition file for all figures except Figs. 9 and S15:

[pathway.m](./pathway_matlab/pathway.m)

### Matlab equation definition files for Figs. 9 and S15:

[pathway_FIG9.m](./pathway_matlab/pathway_FIG9.m)

[pathway_FIGS15.m](./pathway_matlab/pathway_FIGS15.m)

### Matlab run files for individual figures:

#### Figures 1 - 4. Matlab will generate Fig1 - Fig6, corresponding to Panels A - F (a counter will tick off 365 iterations):

Figure 1: [FIG1.m](./pathway_matlab/FIG1.m)

Figure 2: [FIG2.m](./pathway_matlab/FIG2.m)

Figure 3: [FIG3.m](./pathway_matlab/FIG3.m)

Figure 4: [FIG4.m](./pathway_matlab/FIG4.m)

#### Figures 6, 7. Matlab will generate Panels A, B (IGT first case):

Figure 6: [FIG6.m](./pathway_matlab/FIG6.m)

Figure 7: [FIG7.m](./pathway_matlab/FIG7.m)

##### Supplemental Figures

FigS8: [FIGS8.m](./pathway_matlab/FIGS8.m)

FigS9: [FIGS9.m](./pathway_matlab/FIGS9.m)

FigS10: [FIGS10.m](./pathway_matlab/FIGS10.m)

FigS11: [FIGS11.m](./pathway_matlab/FIGS11.m)

FigS12: [FIGS12.m](./pathway_matlab/FIGS12.m)

FigS13: [FIGS13.m](./pathway_matlab/FIGS13.m)

FigS14: [FIGS14.m](./pathway_matlab/FIGS14.m)

FigS14: [FIGS15.m](./pathway_matlab/FIGS15.m)

### XPP source file for all figures except Figs. 9, S10 and S11:

pathway.ode: [pathway.ode](./pathway_xpp/pathway.ode)

### XPP source files Figs. 9, S10 and S11:

pathway_FIG9.ode (insulin reistance induced by insulin): [pathway_FIG9.ode](./pathway_xpp/pathway_FIG9.ode)

pathway_FIGS10_11.ode (non-insulin dependent glucose uptake): [pathway_FIGS10_11.ode](./pathway_xpp/pathway_FIGS10_11.ode)

### XPP set files for individual figures:

#### Figures 1 - 4. xpp will plot the longitudinal response of glucose to meals each day. The envelope of the plot will show the fasting and peak postprandial values. Use the xpp menus to display the other variables:

FIG1.set (IGT-first pathway): [FIG1.set](./pathway_xpp/FIG1.set)

FIG2.set (Sustained IGT): [FIG2.set](./pathway_xpp/FIG2.set)

FIG3.set (IFG-first pathway IGT): [FIG3.set](./pathway_xpp/FIG3.set)

FIG4.set (IGT without CGI): [FIG4.set](./pathway_xpp/FIG4.set)

FIG9_no_cl_no_FE.set (default IGT-first pathway; use pathway_FIG9.ode): [FIG9_no_cl_no_FE.set](./pathway_xpp/FIG9_no_cl_no_FE.set)

FIG9_modest_cl_FE.set (reduced clearance delays T2D; use pathway_FIG9.ode): [FIG9_modest_cl_FE.set](./pathway_xpp/FIG9_modest_cl_FE.set)

FIG9_balanced_cl_FE.set (reduced clearance and induced resistance cancel out; use pathway_FIG9.ode): [FIG9_balanced_cl_FE.set](./pathway_xpp/FIG9_balanced_cl_FE.set)

FIG9_extreme_cl_FE.set (reduced clearance and severe induced resistance accelerates T2D; use pathway_FIG9.ode): [FIG9_extreme_cl_FE.set](./pathway_xpp/FIG9_extreme_cl_FE.set)


##### Supplemental Figures

FIGS8.set (mild GOF KATP defect): [FIGS8.set](./pathway_xpp/FIGS8.set)

FIGS9.set: (severe GOF KATP defect): [FIGS9.set](./pathway_xpp/FIGS9.set)

FIGS10.set (impaired glucose effectiveness; use pathway_FIGS10_11.ode): [FIGS10.set](./pathway_xpp/FIGS10.set)

FIGS11.set (enhanced glucose effectiveness; use pathway_FIGS10_11.ode): [FIGS11.set](./pathway_xpp/FIGS11.set)

FIGS12_fast_si.set  (IGT-first pathway, same as Fig. 1): [FIGS12_fast_si.set](./pathway_xpp/FIGS12_fast_si.set)

FIGS12_slow_si.set (protected from T2D): [FIGS12_slow_si.set](./pathway_xpp/FIGS12_slow_si.set)

FIGS13.set (IFG-first pathway driven by GOF KATP defect): [FIGS13.set](./pathway_xpp/FIGS13.set)

FIGS14.set: (IFG-first pathway driven by impaired vesicle trafficking - small RRP): [FIGS14.set](./pathway_xpp/FIGS14.set)

FIGS15.set (impaired incretin signaling accelerates T2D): [FIGS15.set](./pathway_xpp/FIGS15.set)
