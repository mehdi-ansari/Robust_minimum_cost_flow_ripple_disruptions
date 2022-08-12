#===========================MakeFile==============================
# For Cplex + Concert Projects
# by Hassan Baalbaki www.hassanbaalbaki.com
#================================================================


SYSTEM = x86-64_linux
LIBFORMAT = static_pic

EXDIR = /home/mehdi/RevisionApril2022/MaxModel_BR

#------------------------------------------------------------
#
# When you adapt this makefile to compile your CPLEX programs
# please copy this makefile and set CPLEXDIR and CONCERTDIR to
# the directories where CPLEX and CONCERT are installed.
#
#------------------------------------------------------------
ILOGSTUDIODIR =/opt/cplex/12.10.0
CPLEXDIR = $(ILOGSTUDIODIR)/cplex
CONCERTDIR = $(ILOGSTUDIODIR)/concert




# ---------------------------------------------------------------------
# Compiler selection
# ---------------------------------------------------------------------

CCC = g++
CC = gcc
JAVAC = javac

# ---------------------------------------------------------------------
# Compiler options
# ---------------------------------------------------------------------

CCOPT = -m64 -O -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD -std=c++11


# ---------------------------------------------------------------------
# Link options and libraries
# ---------------------------------------------------------------------

CPLEXBINDIR = $(CPLEXDIR)/bin/$(BINDIST)
CPLEXLIBDIR = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CCLNFLAGS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lilocplex -lcplexdistmip -lcplex -lm -lpthread -ldl



all:${DIRS}
	make all_cpp

execute: all
	make execute_cpp


CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR = $(CPLEXDIR)/include

#EXINC = $(EXDIR)/include
#EXINC = $(EXDIR)/src/cpp
#EXDATA = $(EXDIR)/data
EXSRCCPP = $(EXDIR)


CFLAGS = $(COPT) -I$(CPLEXINCDIR)
CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR)


#------------------------------------------------------------
#
#Build and Execute directory
#
#------------------------------------------------------------



INC_PATH = $(EXSRCCPP)
SRC_PATH = $(EXSRCCPP)
OBJ_PATH = ${EXDIR}/obj



SOURCES = ${wildcard ${SRC_PATH}/*.cpp}
OBJECTS = ${patsubst ${SRC_PATH}/%.cpp, ${OBJ_PATH}/%.o, ${SOURCES}}


INCLUDES = -I${INC_PATH}

DIRS = ${OBJ_PATH}


#------------------------------------------------------------
# make all : to compile all.
# make execute : to compile and execute some instances
#------------------------------------------------------------



CPP_EX = globalplanner




all_cpp: $(CPP_EX)







# ------------------------------------------------------------

clean :
	/bin/rm -rf *.o *~ $(OBJ_PATH)/*.o
	/bin/rm -rf ${SRC_PATH}/*~
	/bin/rm -rf $(CPP_EX)
	/bin/rm -rf *.mps *.ord *.sos *.lp *.sav *.net *.msg *.log *.clp

# ------------------------------------------------------------
#
# The Cpp and *.o Files
#


all: ${DIRS} globalplanner

${OBJ_PATH}:
	mkdir -p $@


${OBJ_PATH}/%.o : ${SRC_PATH}/%.cpp
	$(CCC) ${INCLUDES} -c $(CCFLAGS) $< -o $@

globalplanner: ${OBJECTS}
	$(CCC) $(CCFLAGS) $^ -o globalplanner $(CCLNFLAGS)


# Local Variables:
# mode: makefile
# End:


#======================== END OF MakeFile =========================