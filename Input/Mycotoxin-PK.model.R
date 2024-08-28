# Generic mycotoxin model

States  = { 
  Q_GI,       # Quantity of mycotoxin in the GI compartment (nmol)
  Q_fec,  # Quantity of mycotoxin in feces (nmol)
  Qcpt,       # Quantity in central compartment (nmol)
  Qu,     # Quantity of mycotoxin in urine (nmol)
  Qmet,   # Quantity of metabolite central compartment (nmol)
  Qu_met,     # Quantity of metabolite in urine (nmol)
  AUC,         # AUC of central compartment (nmol-hr/L)
};  

Outputs = {
  Ccpt, # Central compartment concentration
  Cmet, # Central compartment metabolite concentration
  Ccpt_out, # Central compartment concentration (truncated at 1e-15)
  Cmet_out, # Central compartment metabolite concentration (truncated at 1e-15)
  Qu_out, # Amount excreted in urine (truncated at 1e-15)
  Q_fec_out, # Amount excreted in feces (truncated at 1e-15)
  Qu_met_out, # Amount of metabolite excreted in urine (truncated at 1e-15)
};

#population mean
M_lnFgutabs = 0.2;
M_lnkgutelim = -1.05; # Gut elimination
M_lnkufrac = 0.32; # Fraction of elimination that is urine
M_lnVdist = 2.1; # Volume of distribution (L/kg)
M_lnVdistmet = 0; # Volume of distribution for metabolite (L/kg)
M_lnCltot = 0;
M_lnClmet = 0;

#population variance
SD_lnFgutabs = 0.2;
SD_lnkgutelim = 0.2;
SD_lnkufrac = 0.2;
SD_lnVdist = 0.2; 
SD_lnVdistmet = 0.2;
SD_lnCltot = 0.2;
SD_lnClmet = 0.2;

#individual log transformed, z-score
lnFgutabs = 0;
lnkgutelim = 0;
lnkufrac = 0;
lnVdist = 0; 
lnVdistmet = 0;
lnCltot = 0;
lnClmet =0;

#individual parameters
# Oral input modeling
InitDose    = 50; # ingested input at t=0 (nmol)
ConstDoseRate = 0; # Constant dose rate per hour (nmol/h) 
Fgutabs = 3133890;
Fgutabs_tmp = 0.5; 
kgutabs    = 0.35; # Intestinal absorption rate (/h); kgutelim * Fgutabs/(1-Fgutabs)
Cltot = 0.5;
Clmet = 0.5;

# Distribution volumes (L/kg)
Vdist = 0.74;
Vdistmet = 1.0;

# Body weight (kg)
BW = 70;

# Elimination rate constants (/h)
ku       = 0.32;     # Urinary excretion rate constant 
ku_tmp   = 0.1;     # Set it to avoid the value geq ktot
kmet       = 0.2;   # Metabolism rate constant 
kgutelim   = 0.35;     #gut elimination rate

#GSD - Mesidual errors
GSD_Ccpt = 1.1; # Central compartment concentration
GSD_Cmet = 1.1; # Central compartment metabolite concentration
GSD_Qu = 1.1;  # Amount excreted in urine
GSD_Q_fec = 1.1;  # Amount excreted in feces
GSD_Qu_met = 1.1; # Amount of metabolite excreted in urine

Initialize {
  Q_GI = InitDose;
  Fgutabs = exp(M_lnFgutabs +SD_lnFgutabs * lnFgutabs);
  Fgutabs = Fgutabs_tmp; #if gut abs tmp is bigger than 1 write 1, otherwise (:) write gut abs temp
  kgutelim = exp(M_lnkgutelim + SD_lnkgutelim * lnkgutelim);
  kgutabs = kgutelim * Fgutabs/(1-Fgutabs);
  Cltot = exp(M_lnCltot + SD_lnCltot * lnCltot);
  Clmet = exp(M_lnClmet + SD_lnClmet * lnClmet);
  Vdist = exp(M_lnVdist + SD_lnVdist * lnVdist);
  ku_tmp = (Cltot/Vdist) * exp(M_lnkufrac + SD_lnkufrac * lnkufrac);
  ku = ku_tmp;
  kmet = (Cltot/Vdist) - ku;
  Vdistmet = exp(M_lnVdistmet + SD_lnVdistmet * lnVdistmet);
}

Dynamics { 
  dt (Q_GI)  = ConstDoseRate - Q_GI * (kgutabs + kgutelim);
  dt (Q_fec) = Q_GI * kgutelim;
  dt (Qcpt) = Q_GI * kgutabs - Qcpt * (ku + kmet);
  dt (Qu) = Qcpt * ku;
  dt (Qmet) = Qcpt * kmet - Qmet * (Clmet/Vdistmet);
  dt (Qu_met) = Qmet * (Clmet/Vdistmet);
  dt (AUC) = Qcpt/(Vdist*BW);
}

CalcOutputs { # truncate at Q values of 1e-15
  Ccpt = Qcpt / (1000*Vdist*BW);
  Cmet = Qmet / (1000*Vdistmet*BW);
  Ccpt_out = (Ccpt < 1e-15 ? 1e-15 : Ccpt);
  Cmet_out = (Cmet < 1e-15 ? 1e-15 : Cmet);
  Qu_out = (Qu < 1e-15 ? 1e-15 : Qu);
  Q_fec_out = (Q_fec < 1e-15 ? 1e-15 : Q_fec);
  Qu_met_out = (Qu_met < 1e-15 ? 1e-15 : Qu_met);
}

End.
