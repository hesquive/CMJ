% Hugo Esquivel, 2021
% -

clearvars; close all; clc;

% ----------------------------------------------------------------------------------------------------------------------
% Input:
% ----------------------------------------------------------------------------------------------------------------------
numRealizations=5e4;
numRandomBasisVectorsPerDimension=4;

x=getRealizationsFromRandomVariables(numRealizations); % this should be understood as the simulation input

fact=getStructuralResponse();
f=getEvaluationOfStructuralResponse(fact,x); % this should be understood as the simulation output
% ----------------------------------------------------------------------------------------------------------------------


% ----------------------------------------------------------------------------------------------------------------------
% Approximating simulation output (aka structural response):
% ----------------------------------------------------------------------------------------------------------------------
Psi=getRandomBasisVectors(numRandomBasisVectorsPerDimension);
Psinum=getNumericRandomBasisVectors(Psi);

coeff=getCoefficientsOfApproximateStructuralResponse(x,f,Psinum);
fapprox=getApproximateStructuralResponse(coeff,Psi);

fprintf('Mean-Squared Error: %.4f\n\n',computeMeanSquaredError(x,f,fapprox))
% ----------------------------------------------------------------------------------------------------------------------


% ----------------------------------------------------------------------------------------------------------------------
% Plots:
% ----------------------------------------------------------------------------------------------------------------------
%plotActualStructuralResponse(fact)
%plotApproximateStructuralResponse(fapprox)
%plotRandomBasisVectors(Psinum)
plotApproximationError(fact,fapprox)
% ----------------------------------------------------------------------------------------------------------------------
