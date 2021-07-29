% Hugo Esquivel, 2021
% -

clearvars; close all; clc;

% ----------------------------------------------------------------------------------------------------------------------
% Input:
% ----------------------------------------------------------------------------------------------------------------------
numRandomVariables=2; % number of random variables

probabilityInfo.name=cell(numRandomVariables,1); % distribution's name
probabilityInfo.pars=cell(numRandomVariables,1); % distribution's parameters (if any)

% 1st random variable (this can be understood as fc after normalization):
rv=1;
probabilityInfo.name{rv}='Beta'; % options: Uniform, Beta, Normal

if strcmpi(probabilityInfo.name{rv},'Beta')
    probabilityInfo.pars{rv}=[2,5]; % parameters: alpha and beta
end

% 2nd random variable (this can be understood as fy after normalization):
rv=2;
probabilityInfo.name{rv}='Beta'; % options: Uniform, Beta, Normal

if strcmpi(probabilityInfo.name{rv},'Beta')
    probabilityInfo.pars{rv}=[2,5]; % parameters: alpha and beta
end

numRealizations=5e4; % number of realizations
numRandomBasisVectorsPerDimension=4; % cardinality of random basis over each random dimension

x=getRealizationsFromRandomVariables(numRealizations,probabilityInfo); % this can be understood as the simulation input

fact=getStructuralResponse(probabilityInfo); % actual structural response (hidden by the simulation; i.e. unknown)
f=getEvaluationOfStructuralResponse(fact,x,probabilityInfo); % this can be understood as the simulation output
% ----------------------------------------------------------------------------------------------------------------------


% ----------------------------------------------------------------------------------------------------------------------
% NISP: Approximating simulation output (aka structural response):
% ----------------------------------------------------------------------------------------------------------------------
probabilityInfo=getProbabilityDistribution(probabilityInfo);

Psi=getRandomBasisVectors(numRandomBasisVectorsPerDimension,probabilityInfo);
Psinum=getNumericRandomBasisVectors(Psi,probabilityInfo);

coeff=getCoefficientsOfApproximateStructuralResponse(x,f,Psinum);
g=getApproximateStructuralResponse(coeff,Psi);

fprintf('Mean-Squared Error: %.4f\n\n',computeMeanSquaredError(x,f,g,probabilityInfo))
% ----------------------------------------------------------------------------------------------------------------------


% ----------------------------------------------------------------------------------------------------------------------
% SOBOL: Preprocessing (simplifying notation of pdf and supp for simplicity):
% ----------------------------------------------------------------------------------------------------------------------
syms x [1,numRandomVariables]

pdf=cell(numRandomVariables,1);
supp=cell(numRandomVariables,1);

for i=1:numRandomVariables
    pdf{i}=probabilityInfo.pdf{i}(x(i));
    supp{i}=probabilityInfo.supp{i};
end
% ----------------------------------------------------------------------------------------------------------------------


% ----------------------------------------------------------------------------------------------------------------------
% SOBOL: Decomposition of structural response:
% ----------------------------------------------------------------------------------------------------------------------
g0=int(int(g*pdf{1}*pdf{2},x1,supp{1}),x2,supp{2}); % total mean

g1=int(g*pdf{2},x2,supp{2})-g0;
g2=int(g*pdf{1},x1,supp{1})-g0;

g12=g-(g0+g1(x1)+g2(x2));
% ----------------------------------------------------------------------------------------------------------------------


% ----------------------------------------------------------------------------------------------------------------------
% SOBOL: Checking orthogonality of decomposition:
% ----------------------------------------------------------------------------------------------------------------------
tol=1e-8;

check1=abs(int(g1*pdf{1},x1,supp{1}))<tol;
check2=abs(int(g2*pdf{2},x2,supp{2}))<tol;

check12=abs(int(int(g12*pdf{1}*pdf{2},x1,supp{1}),x2,supp{2}))<tol;

if ~(check1 && check2 && check12)
    error('Orthogonality check did not pass.')
end
% ----------------------------------------------------------------------------------------------------------------------


% ----------------------------------------------------------------------------------------------------------------------
% SOBOL: Decomposition of total variance:
% ----------------------------------------------------------------------------------------------------------------------
Dg=int(int(g^2*pdf{1}*pdf{2},x1,supp{1}),x2,supp{2})-g0^2; % total variance

D1=int(g1^2*pdf{1},x1,supp{1});
D2=int(g2^2*pdf{2},x2,supp{2});

D12=int(int(g12^2*pdf{1}*pdf{2},x1,supp{1}),x2,supp{2});
% ----------------------------------------------------------------------------------------------------------------------


% ----------------------------------------------------------------------------------------------------------------------
% SOBOL: Contributions to total variance:
% ----------------------------------------------------------------------------------------------------------------------
S1=D1/Dg; % contribution to total variance due to x1
S2=D2/Dg; % contribution to total variance due to x2

S12=D12/Dg; % contribution to total variance due to the interaction between x1 and x2
% ----------------------------------------------------------------------------------------------------------------------


% ----------------------------------------------------------------------------------------------------------------------
% SOBOL: Summary:
% ----------------------------------------------------------------------------------------------------------------------
fprintf('Summary:\n\n')
fprintf('Total mean: %.6f.\n',g0)
fprintf('Total variance: %.4f.\n',Dg)
fprintf('\n')
fprintf('Random variable x1 contributes to the total variance response by: %.2f%%.\n',S1*100)
fprintf('Random variable x2 contributes to the total variance response by: %.2f%%.\n',S2*100)
fprintf('Random variables x1 and x2 contribute to the total variance response by: %.2f%%.\n',S12*100)
% ----------------------------------------------------------------------------------------------------------------------


% ----------------------------------------------------------------------------------------------------------------------
% NISP: Plots:
% ----------------------------------------------------------------------------------------------------------------------
plotJointProbabilityDistribution(probabilityInfo)
%plotActualStructuralResponse(fact,probabilityInfo)
%plotApproximateStructuralResponse(g,probabilityInfo)
%plotRandomBasisVectors(Psinum,probabilityInfo)
plotApproximationError(fact,g,probabilityInfo)
% ----------------------------------------------------------------------------------------------------------------------
