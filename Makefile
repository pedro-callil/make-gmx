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

################################################################################
#                                                                              #
#                     EDITIONS USUALLY SHOULD BE MADE HERE                     #
#                                                                              #
################################################################################

# Input files
PKSRC	= to-packmol.inp
MINSRC	= Minimization.mdp
EQUISRC	= Equilibration.mdp
PRODSRC	= Production.mdp
TOPSRC	= topol.top

# Output files
PKOUT	= from-packmol
EDITOUT	= from-editconf
MINOUT	= from-minimization
NDXOUT  = from-ndx
EQUIOUT	= from-equilibration
PRODOUT	= from-production

# Runtime variables
export ROC_ENABLE_PRE_VEGA=1

################################################################################
#                                                                              #
#        INTERNAL VARIABLES, SHELL ONE-LINERS AND MAKE TARGETS/COMMANDS        #
#                                                                              #
################################################################################

ITPS	= $(shell cat $(TOPSRC) |grep ^\#include |sed 's/^\#include//g;s/\"//g')
PDBS	= $(shell cat $(PKSRC) |grep ^structure |sed 's/^structure\ //g')

SRC	= $(PKSRC) $(MINSRC) $(EQUISRC) $(PRODSRC) $(TOPSRC)

NDXSTR  = $$(cat $(EQUISRC) |grep ^tc-grps |tr '[:blank:]' '\n' \
		|grep '_' \
		|sed 's/^/\"/g;s/$$/\"\\nq/g;s/_/\"\ |\ \"/g')

BOXX	= $$(cat $(PKSRC) |grep box |tr -d '[A-Za-z]' \
	  	|tr -s '[:blank:]' '-' \
		|cut -d '-' -f 1,2,5 \
		|sed 's/^-/-(/g;s/$$/)\/9.5/' \
		|bc -l |cut -c -8 \
		| awk '{if($$1>max) {max=$$1}} END {print max}')
BOXY	= $$(cat $(PKSRC) |grep box |tr -d '[A-Za-z]' \
	  	|tr -s '[:blank:]' '-' \
		|cut -d '-' -f 1,3,6 \
		|sed 's/^-/-(/g;s/$$/)\/9.5/' \
		|bc -l |cut -c -8 \
		|awk '{if($$1>max) {max=$$1}} END {print max}')
BOXZ	= $$(cat $(PKSRC) |grep box |tr -d '[A-Za-z]' \
	  	|tr -s '[:blank:]' '-' \
		|cut -d '-' -f 1,4,7 \
		|sed 's/^-/-(/g;s/$$/)\/9.5/' \
		|bc -l |cut -c -8 \
		|awk '{if($$1>max) {max=$$1}} END {print max}')
BOXTYPE	= triclinic

.PHONY: all clean

all: $(PRODOUT).gro

clean:
	rm *from-* mdout.mdp

$(PRODOUT).gro: $(PRODSRC) $(ITPS) $(EQUIOUT).gro $(NDXOUT).ndx
	gmx grompp -f $(PRODSRC) \
		-c $(EQUIOUT).gro \
		-p $(TOPSRC) \
		-o $(PRODOUT).tpr \
		-n $(NDXOUT).ndx
	gmx mdrun -deffnm $(PRODOUT) -v
	mv $(PRODSRC) $(PRODSRC).bak
	mv $(PRODOUT).xtc $(PRODOUT).old.xtc
	mv $(PRODOUT).edr $(PRODOUT).old.edr
	mv $(PRODOUT).log $(PRODOUT).old.log
	sed 's/;energygrps/energygrps/g' $(PRODSRC).bak > $(PRODSRC)
	gmx grompp -f $(PRODSRC) \
		-c $(MINOUT).gro \
		-p $(TOPSRC) \
		-o $(PRODOUT).tpr \
		-n $(NDXOUT).ndx
	gmx mdrun -deffnm $(PRODOUT) -v \
		-s $(PRODOUT).tpr \
		-rerun $(PRODOUT).old.xtc
	mv $(PRODSRC).bak $(PRODSRC)

$(EQUIOUT).gro: $(EQUISRC) $(ITPS) $(MINOUT).gro $(NDXOUT).ndx
	gmx grompp -f $(EQUISRC) \
		-c $(MINOUT).gro \
		-p $(TOPSRC) \
		-o $(EQUIOUT).tpr \
		-n $(NDXOUT).ndx \
		-maxwarn 1
	gmx mdrun -deffnm $(EQUIOUT) -v
	mv $(EQUISRC) $(EQUISRC).bak
	mv $(EQUIOUT).xtc $(EQUIOUT).old.xtc
	mv $(EQUIOUT).edr $(EQUIOUT).old.edr
	mv $(EQUIOUT).log $(EQUIOUT).old.log
	sed 's/;energygrps/energygrps/g' $(EQUISRC).bak > $(EQUISRC)
	gmx grompp -f $(EQUISRC) \
		-c $(MINOUT).gro \
		-p $(TOPSRC) \
		-o $(EQUIOUT).tpr \
		-n $(NDXOUT).ndx \
		-maxwarn 1
	gmx mdrun -deffnm $(EQUIOUT) -v \
		-s $(EQUIOUT).tpr \
		-rerun $(EQUIOUT).old.xtc
	mv $(EQUISRC).bak $(EQUISRC)

$(MINOUT).gro: $(MINSRC) $(ITPS) $(EDITOUT).gro
	gmx grompp -f $(MINSRC) \
		-c $(EDITOUT).gro \
		-p $(TOPSRC) \
		-o $(MINOUT).tpr
	gmx mdrun -deffnm $(MINOUT) -v

$(NDXOUT).ndx: $(EDITOUT).gro
	echo -e "$(NDXSTR)" | \
		gmx make_ndx -f $(EDITOUT).gro -o $(NDXOUT).ndx

$(EDITOUT).gro: $(PKOUT).pdb
	gmx editconf -f $(PKOUT).pdb \
		-bt $(BOXTYPE) \
		-box $(BOXX) $(BOXY) $(BOXZ) \
		-o $(EDITOUT).gro

$(PKOUT).pdb: $(PKSRC) $(PDBS)
	packmol < $(PKSRC)



