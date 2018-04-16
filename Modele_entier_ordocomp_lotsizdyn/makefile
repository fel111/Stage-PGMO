# Makefile for the Branch-and-Price for the VRPTF designed by SUN

#pour mac osX
# FRAMEWORKS = -framework IOKit -framework CoreFoundation -framework Ruby

#pour les include

INCLUDES_PATH= \
	-I/usr/include \
	-I/home/fgoupil/Téléchargements/scipoptsuite-5.0.1/scip/src \
	-I/home/fgoupil/CPLEX_Studio1271/cplex/include \
	-I/home/fgoupil/CPLEX_Studio1271/concert/include \
	-I/home/fgoupil/Documents/Stage_PGMO	

# Les librairies que l'on va utiliser
LIBS= \
	-lscip-5.0.1.linux.x86_64.gnu.opt \
	-lnlpi.cppad-5.0.1.linux.x86_64.gnu.opt \
	-llpispx2-5.0.1.linux.x86_64.gnu.opt \
	-lsoplex-3.1.1.linux.x86_64.gnu.opt \
	-ltpinone-5.0.1.linux.x86_64.gnu.opt \
	-lpthread \
	-lm \
	-lilocplex \
	-lconcert \
	-lcplex \
	-lstdc++

# Les endroits o� on peut trouver des librairies
LIBS_PATH= \
	-L/usr/lib \
	-L/usr/local/lib \
	-L/opt/local/lib \
	-L/Developer/usr/lib/gcc/i686-apple-darwin10/4.2.1 \
	-L/home/fgoupil/Téléchargements/scipoptsuite-5.0.1/scip/lib/static \
	-L/home/fgoupil/Téléchargements/scipoptsuite-5.0.1/soplex/lib \
	-L/home/fgoupil/CPLEX_Studio1271/cplex/lib/x86-64_linux/static_pic \
	-L/home/fgoupil/CPLEX_Studio1271/concert/lib/x86-64_linux/static_pic

# Le compilateur que nous utiliserons
COMPILO=g++
#/usr/bin/g++-4.2
#g++
#gcc

# Les options de compilation des sources
#    -c    : Ne fait pas un ex�cutable, mais seulement un .o
#    -g    : Fait un objet d�buggable
#    -W    : Affiche un maximum de warnings si le code C++ est douteux
#    -Wall : Affiche un maximum de warnings si le code C++ est douteux
#    -ansi : Assure que le code C++ suit la norme ANSI, et est donc portable.
#    -O2   : Lorsque votre programme est au point, l'option -O2 permet d'optimiser votre programme (le rendre plus rapide) 
#FLAGS= -O3 -ansi
FLAGS= -W -Wall -DIL_STD -std=c++11 -g -O3
#-m64 
#-DNDEBUG deactivate assertions
#FLAGS= -g -W -Wall -DIL_STD 


OBJ=	modele_entier.o \
		../Ordo_compact/ordo_comp.o \
		../Ordo_compact/ordo_cplex.o \
		../Lot_Sizing/Lot_Sizing_Dynamique/lotsizdyn2.o \
		../Lot_Sizing/Lot_Sizing_Compact/lotsizentCPX.o \
		../Lot_Sizing/Lot_Sizing_Compact/lotsizcontcom.o \
		../Lot_Sizing/Lot_Sizing_Compact/lotsizcontCPX.o \
		../data_struct.o \
		../Modele_compact/compact.o \
		../Modele_compact/modele_entier_cplex.o \
		../Lecteur_Fichiers/lecteur_taches.o \
		../Lecteur_Fichiers/lecteur_pwd.o



# Liste des objets que l'on va obtenir, et qu'il faudra linker.
# pour exe
EXEC_NAME1 = modele_entier
#EXEC_NAME1 = example
OBJ1= $(OBJ)

# pour exe
#EXEC_NAME2 = eval
#OBJ2= main_simple_eval.o $(OBJ)

# Apr�s avoir d�fini les macro, on d�finit des cibles. Une cible est une expression
# du style :
#   <nom> : <dep1> <dep2> ... <depn>
#   <tab> <commande>
#
# Pour que make ex�cute une cible <nom>, il faut taper � l'endroit o� se trouve le makefile
#   make <nom>
# Si on ne donne pas de <nom>, la premi�re cible est ex�cut�e.
# L'ex�cution de la cible <nom> se d�roule comme suit :
#   Pour chacun des <depi>, regarder si <depi> n'est pas un nom de fichier.
#     Si <depi> n'est pas un nom de fichier
#        Alors ex�cuter la cible <depi> (qui doit �tre d�finie dans ce fichier makefile)
#              puis ex�cuter <commande>
#     Si <depi> est un fichier 
#        Si il est a jour (sa date est ant�rieur � la date de la derni�re commande make)
#           Ne rien faire
#        Si il n'est pas a jour et qu'il existe une cible <depi>, l'ex�cuter
#              puis ex�cuter <commande>
#        Si il n'est pas a jour et qu'il n'existe pas de cible <depi>
#               ex�cuter <commande>.
#
# Ce m�canisme permet, entre autre chose, de ne recompiler que les objets dont le source
# C++ correspondant � �t� modifi�.

######################
#                    #
# Cibles Principales #
#                    #
######################

all : $(EXEC_NAME1) $(EXEC_NAME2)

# Pour compiler algo_evo, il faut v�rifier que tous les .o
# sont � jour (les recompiler dans le cas contraire) puis lancer
# la commande de linkage.
$(EXEC_NAME1) : $(OBJ1)
	$(COMPILO) -o $(EXEC_NAME1) $(OBJ1) $(LIBS_PATH) $(LIBS) 
#	$(COMPILO) -g -o $(EXEC_NAME1) $(OBJ1) $(LIBS_PATH) $(LIBS) $(FRAMEWORKS)

# Taper 'make clean' provoque l'effacement de tous les .o, et des *~ laiss�s 
# par emacs. Il n'y a pas de d�pendence pour cette cible.
clean : 
	rm -f *.o *~ 
 
# Taper 'make clear' fait un 'make clean' puis efface en plus l'ex�cutable.
clear : 
	make clean;
	rm *.exe;
	rm -f $(EXEC_NAME1) $(EXEC_NAME2);
	rm -f -r html


###########################
#                         #
# Compilation des sources #
#                         #
###########################
# $@  	Le nom de la cible
# $< 	Le nom de la premi�re d�pendance
# $^ 	La liste des d�pendances
# $? 	La liste des d�pendances plus r�centes que la cible
# $* 	Le nom du fichier sans suffixe
# Source bidon, bien utile pour faire des copier-coller et obtenir les suivantes
# par en remplacement de xx.
# Le sens est que xx.o s'obtent � partir de xx.cc (la seule d�pendance ici) quand
# le fichier xx.cc vient d'�tre modifi�
#xx.o : xx.cc
#	$(COMPILO) $(FLAGS) $(INCLUDE_PATH) -o xx.o xx.cc
#%.o : %.cc
#	$(COMPILO) $(FLAGS) -o $< $@
#%.o: %.cc
#	$(COMPILO) $(FLAGS) $(INCLUDES_PATH) -o $@ -c $<
#%.o: %.c
#	gcc $(FLAGS) $(INCLUDES_PATH) -o $@ -c $<

%.o: %.CPP
	$(COMPILO) $(FLAGS) $(INCLUDES_PATH) -o $@ -c $<

%.o: %.cpp
	$(COMPILO) $(FLAGS) $(INCLUDES_PATH) -o $@ -c $<


# Cree une documentation avec Doxygen et le ficheir de configuration 'doxygen_config'
#doc :
#	doxygen DoxygenOutput/doxygen_config

#clearDoc :
#	rm -f -r html
