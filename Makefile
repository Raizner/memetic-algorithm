ARCH_MACH = $(shell uname -m)
OS_TYPE   = $(shell uname -s)

############### name of target
#
TARGET		= ./dolphin

########################################################################
##########        Linux Cluster Environment                   ##########
########################################################################

CPPFLAGS := $(CPPFLAGS) -D__LINUX__ -DEBUG

############### C++ compiler commands
#
CCC		= /usr/bin/g++ 
CC		= /usr/bin/gcc

############### Other Packages
#

############### libraries needed to generate program
#
LDLIBS		:=  -lm

############### objects needed to generate target
#
OBJS = 	Dolphin.o \
ObjectiveFunctions/ObjectiveFunction.o \
ObjectiveFunctions/FFMS.o ObjectiveFunctions/FWeierstrass.o ObjectiveFunctions/FSphere.o ObjectiveFunctions/FRosenbrock.o ObjectiveFunctions/FRastrigin.o ObjectiveFunctions/FGriewank.o ObjectiveFunctions/FAckley.o ObjectiveFunctions/FBump.o ObjectiveFunctions/FScaffer.o ObjectiveFunctions/FGriewankRosenbrock.o ObjectiveFunctions/FSchwefel102.o ObjectiveFunctions/FElliptic.o ObjectiveFunctions/FStep.o ObjectiveFunctions/FSchwefel102Noisy.o \
Rng/Bernoulli.o Rng/Binomial.o Rng/Cauchy.o Rng/DiffGeometric.o Rng/DiscreteUniform.o Rng/Erlang.o Rng/Geometric.o Rng/GlobalRng.o Rng/HyperGeometric.o Rng/LogNormal.o Rng/NegExponential.o Rng/Normal.o Rng/Poisson.o Rng/Rng.o Rng/Uniform.o Rng/Weibull.o Rng/ihs.o \
Elements/Chromosome_Binary.o Elements/Chromosome_Real.o \
Operators/Mutations/Mutation_Gaussian.o Operators/Mutations/Mutation_BitFlip.o Operators/Scalings/Scaling_Linear.o \
LocalSearches/LocalSearch.o LocalSearches/LocalSearch_DSCG.o LocalSearches/LocalSearch_DFP.o LocalSearches/LocalSearch_ES.o \
Utilities/Statistics.o \
MemeticAlgorithms/MALamarcBaldwin.o \


PURIFY	= purify $(PFLAGS)

CCFLAGS		=  $(CPPFLAGS)
COMPILE.cc	=  $(CCC) $(CCFLAGS) $(CPPFLAGS)  -c -w
LINK.cc		=  $(CCC) $(CCFLAGS) $(CPPFLAGS)  $(LDFLAGS)



############### options for preprocessor (include directories)
#
CPPFLAGS	:=\
		$(CPPFLAGS)  			\
		-Wall 				\
		-I.


############### options for compiler
#

CCFLAGS:= -g -Wall -O3 

CXXFLAGS:=$(CXXFLAGS) $(CFLAGS)


############### options for linker
#

LDFLAGS		:= $(LDFLAGS) 			\

###############	make targets
#
all:		$(TARGET) $(INCLUDE)


.c:
		$(LINK.c) -o $@ $< $(LDLIBS)
.c.o:
		$(COMPILE.c) -o $@ $<
.cpp:
		$(LINK.cc) -o $@ $< $(LDLIBS)
.cpp.o:
		$(COMPILE.cc) -o $@ $<
.f:
		$(LINK.f) -o $@ $< $(LDLIBS)
.f.o:
		$(COMPILE.f) -o $@ $<
		
############### suffixes list
#
.SUFFIXES:	.o .so .a .c .cpp .h

.PHONY:		all clean depend

############### make shared object (dynamic link library)
#
$(TARGET):	$(OBJS)
		$(LINK.cc) -o $@ $(OBJS) $(LDLIBS)

###############	clean up
#
clean:
		-$(RM) $(OBJS)


###############	update dependencies
#
depend:
		makedepend -Y -- $(CPPFLAGS) *.c *.cc *.cpp 2> /dev/null

cleandepend:
		makedepend

# DO NOT DELETE
