#Example run with default parameters
Integrate (Lsodes, 1e-4, 1e-6, 1);
OutputFile("Mycotoxin-default-run.out");

M_lnFgutabs = -1.76;
M_lnkgutelim = -0.7; # Gut elimination
M_lnkufrac = -0.7; # Fraction of elimination that is urine
M_lnCltot = -3.33;
M_lnClmet = -3.33;
M_lnVdist = -0.10; # Volume of distribution (L/kg)
M_lnVdistmet = -0.10; # Volume of distribution for metabolite (L/kg)

  Simulation {
    InitDose    = 13172;    # ingested dose (ng)
    BW         = 66;      # body weight (kg)

    PrintStep(Ccpt_out, Cmet_out, Qu_out, Qu_met_out, Q_fec_out, 0, 48, 0.1);
  }

End.
