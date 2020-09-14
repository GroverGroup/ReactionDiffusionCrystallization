# ReactionDiffusionCrystallization
Code used in Salami et al. 2020 Reaction Chemistry &amp; Engineering

The code included in this repository was used for all modeling in the following study: Salami, H.; Lagerman, C.E.; Harris, P.R.; McDonald, M.A.; Bommarius, A.S.; Rousseau, R.W.; Grover, M.A. Model development for enzymatic reactive crystallization of β-lactam antibiotics: a reaction–diffusion-crystallization approach. Reaction Chemistry & Engineering 2020, doi:10.1039/D0RE00276C

Free use of this code is allowed provided the above study it cited.

Main.m is the main simulation of the reactive crystallization of beta-lactam antibiotics with an enzyme immobilized on a porous support
pH_help.m calculates the pH of a solution of the relevant species using a charge balance approach. The pKa values are corrected for ionic strength
rate_calc.m calculates the rate of change of concentration within porous enzyme carrier beads and the rate of crystal growth outside the carrier beads
RxnDiffusion.m supports rate_calc.m by providing reaction rates within the beads
titr.m calculates the amount of sodium hydroxide needed to titrate the solution to the desired pH value
Crystal.m supports rate_calc.m by providing crystal growth rate and nucleation rates outside the beads
Csat.m supports Crystal.m by by providing the saturation concentration of each species at the given pH value.
