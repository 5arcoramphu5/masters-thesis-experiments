# a list of all the programs in your package
PROGS = experiments/comparing-trajectories-test/main experiments/calculate-normal-form/main experiments/periodic-orbit/main

# a list of all your units to be linked with your programs
OTHERS = shared/test_cases/test_cases_collection

# path to normal-forms repository
NFDIR = ../normal-forms/.obj/source

# files to link from normal-forms repository
NFOBJFILES = ${NFDIR}/NormalFormFinder/helperFunctions.o

# path to directory, where script capd-config is located
CAPDBINDIR = ~/libraries/CAPD/build/bin/

# setting compiler and linker flags
CAPDFLAGS = `${CAPDBINDIR}capd-config --cflags`
CAPDLIBS = `${CAPDBINDIR}capd-config --libs`
CXXFLAGS += ${CAPDFLAGS} 

# directory where object and dependancy files will be created
OBJDIR = .obj/

OTHERS_OBJ = ${OTHERS:%=${OBJDIR}%.o}
OBJ_FILES = ${OTHERS_OBJ} ${PROGS:%=${OBJDIR}%.o}

.PHONY: all
all: ${PROGS}

# rule to link executables
${PROGS}: % : ${OBJDIR}%.o ${OTHERS_OBJ}
	${CXX} -o $@.exe $< ${OTHERS_OBJ} ${NFOBJFILES} ${CAPDLIBS}

# include files with dependencies
-include ${OBJ_FILES:%=%.d}

#rule to compile .cpp files and generate corresponding files with dependencies
${OBJ_FILES}: ${OBJDIR}%.o : %.cpp
	@mkdir -p ${dir $@}
	$(CXX) -g ${CXXFLAGS} -MT $@ -MD -MP -MF ${@:%=%.d} -c -o $@ $< -std=c++20

# rule to clean all object files, dependencies and executables
.PHONY: clean
clean:
	rm -r ${OBJDIR}*


	

