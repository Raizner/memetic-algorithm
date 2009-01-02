/*!
 * <File comment goes here!!>
 * 
 * Copyright (c) 2008 by <Quang Huy / NTU>
 */

// Mutations
#include "Operators/Mutations/Mutation.h"
#include "Operators/Mutations/Mutation_BitFlip.h"
#include "Operators/Mutations/Mutation_Gaussian.h"

// Crossovers
#include "Operators/Crossovers/Crossover.h"
#include "Operators/Crossovers/Crossover_nPoint.h"
#include "Operators/Crossovers/Crossover_Uniform.h"

// Scaling & Selections 
#include "Operators/Scalings/Scaling.h"
#include "Operators/Scalings/Scaling_Linear.h"

#include "Operators/Selections/Selection.h"
#include "Operators/Selections/Selection_RouletteWheel.h"

// Recombinations
#include "Operators/Recombinations/Recombination.h"
#include "Operators/Recombinations/Recombination_KeepBest.h"
#include "Operators/Recombinations/Recombination_MuPlusLambda.h"
#include "Operators/Recombinations/Recombination_MuCommaLambda.h"
