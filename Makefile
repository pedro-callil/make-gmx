################################################################################
# GMX MAKEFILE - AVOID RERUNNING SIMULATIONS WHILE STILL SCRIPTING THEM EASILY #
#                                                                              #
# Sometimes one needs to rerun a single part of a simulation, while keeping    #
# others equal; for instance, we might make a mistake in our production .mdp   #
# file, and while the production resulting from it is useless, rerunning the   #
# equilibration is a waste of time; however, avoiding that might either mean   #
# checking for file existence for each script action, or editing the script    #
# to make ad-hoc exclusions.                                                   #
#                                                                              #
# In the interest of saving time and reducing complexity, a template to script #
# our simulations using GNU Make, an utility that already performs that check, #
# is exhibited below. By executing it (i.e. executing "make" in its directory) #
# one can rerun only the sections of the simulation that are affected by file  #
# deletions and similar actions.                                               #
################################################################################

PKSRC	= to-packmol.inp
MINSRC	= Minimization.mdp
EQUISRC	= Equilibration.mdp
PRODSRC	= Production.mdp
TOPSRC	= topol.top

BOXTYPE	= triclinic
BOXX	= 13.5
BOXY	= 13.5
BOXZ	= 13.5

ITPS	= $(shell cat $(TOPSRC) |grep ^\#include |sed 's/^\#include//g;s/\"//g')
PDBS	= $(shell cat $(PKSRC) |grep ^structure |sed 's/^structure\ //g')

SRC	= $(PKSRC) $(MINSRC) $(EQUISRC) $(PRODSRC) $(TOPSRC)

PKOUT	= from-packmol
EDITOUT	= from-editconf
MINOUT	= from-minimization
EQUIOUT	= from-equilibration
PRODOUT	= from-production

.PHONY: all clean

all: $(PRODOUT).gro

clean:
	rm from-*

$(PRODOUT).gro: $(SRC) $(ITPS) $(EQUIOUT).gro
	gmx grompp -f $(PRODSRC) \
		-c $(EQUIOUT).gro \
		-p $(TOPSRC) \
		-o $(PRODOUT).tpr
	gmx mdrun -deffnm $(PRODOUT) -v

$(EQUIOUT).gro: $(SRC) $(ITPS) $(MINOUT).gro
	gmx grompp -f $(EQUISRC) \
		-c $(MINOUT).gro \
		-p $(TOPSRC) \
		-o $(EQUIOUT).tpr \
		-maxwarn 1
	gmx mdrun -deffnm $(EQUIOUT) -v

$(MINOUT).gro: $(SRC) $(ITPS) $(EDITOUT).gro
	gmx grompp -f $(MINSRC) \
		-c $(EDITOUT).gro \
		-p $(TOPSRC) \
		-o $(MINOUT).tpr
	gmx mdrun -deffnm $(MINOUT) -v

$(EDITOUT).gro: $(PKOUT).pdb
	gmx editconf -f $(PKOUT).pdb \
		-bt $(BOXTYPE) \
		-box $(BOXX) $(BOXY) $(BOXZ) \
		-o $(EDITOUT).gro

$(PKOUT).pdb: $(SRC) $(PDBS)
	packmol < $(PKSRC)



