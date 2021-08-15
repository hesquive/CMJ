% Hugo Esquivel, 2021
% -

clearvars; close all; clc;


% ----------------------------------------------------------------------------------------------------------------------
% INPUT:
% ----------------------------------------------------------------------------------------------------------------------
syms x1 x2 x3 % random variables x1, x2 and x3

g=getStructuralResponse(x1,x2,x3); % structural response g

probabilityDistribution='Uniform'; % options: Uniform or Beta

switch probabilityDistribution
    case 'Uniform'
        f1=1/2; % probability density of x1 (assuming that x1 is a uniformly distributed random variable in [-1,1])
        f2=1/2; % probability density of x2 (assuming that x2 is a uniformly distributed random variable in [-1,1])
        f3=1/2; % probability density of x3 (assuming that x3 is a uniformly distributed random variable in [-1,1])
        
    case 'Beta'
        % probability density of x1 (assuming that x1 is a beta distributed random variable in [-1,1])
        alpha1=2;
        beta1=5;
        f1=(x1+1)^(alpha1-1)*(1-x1)^(beta1-1)/(2^(alpha1+beta1-1)*beta(alpha1,beta1));
        
        % probability density of x2 (assuming that x2 is a beta distributed random variable in [-1,1])
        alpha2=2;
        beta2=5;
        f2=(x2+1)^(alpha2-1)*(1-x2)^(beta2-1)/(2^(alpha2+beta2-1)*beta(alpha2,beta2));
        
        % probability density of x3 (assuming that x3 is a beta distributed random variable in [-1,1])
        alpha3=2;
        beta3=5;
        f3=(x3+1)^(alpha3-1)*(1-x3)^(beta3-1)/(2^(alpha3+beta3-1)*beta(alpha3,beta3));
end

f12=f1*f2; % joint probability density of x1 and x2
f13=f1*f3; % joint probability density of x1 and x3
f23=f2*f3; % joint probability density of x2 and x3

f123=f1*f2*f3; % joint probability density of x1, x2 and x3
% ----------------------------------------------------------------------------------------------------------------------



% ----------------------------------------------------------------------------------------------------------------------
% DECOMPOSITION OF STRUCTURAL RESPONSE:
% ----------------------------------------------------------------------------------------------------------------------
g0=int(int(int(g*f123,x1,-1,1),x2,-1,1),x3,-1,1); % total mean

g1=int(int(g*f23,x2,-1,1),x3,-1,1)-g0; 
g2=int(int(g*f13,x1,-1,1),x3,-1,1)-g0;
g3=int(int(g*f12,x1,-1,1),x2,-1,1)-g0;

g12=int(g*f3,x3,-1,1)-(g0+g1+g2);
g13=int(g*f2,x2,-1,1)-(g0+g1+g3);
g23=int(g*f1,x1,-1,1)-(g0+g2+g3);

g123=g-(g0+g1+g2+g3+g12+g13+g23);
% ----------------------------------------------------------------------------------------------------------------------



% ----------------------------------------------------------------------------------------------------------------------
% CHECKING ORTHOGONALITY OF DECOMPOSITION:
% ----------------------------------------------------------------------------------------------------------------------
check1=int(g1*f1,x1,-1,1)==0;
check2=int(g2*f2,x2,-1,1)==0;
check3=int(g3*f3,x3,-1,1)==0;

check12=int(int(g12*f12,x1,-1,1),x2,-1,1)==0;
check13=int(int(g13*f13,x1,-1,1),x3,-1,1)==0;
check23=int(int(g23*f23,x2,-1,1),x3,-1,1)==0;

check123=int(int(int(g123*f123,x1,-1,1),x2,-1,1),x3,-1,1)==0;

if ~(check1 && check2 && check3 && check12 && check13 && check23 && check123)
    error('Orthogonality check did not pass.')
end
% ----------------------------------------------------------------------------------------------------------------------



% ----------------------------------------------------------------------------------------------------------------------
% DECOMPOSITION OF TOTAL VARIANCE:
% ----------------------------------------------------------------------------------------------------------------------
Dg=int(int(int(g^2*f123,x1,-1,1),x2,-1,1),x3,-1,1)-g0^2; % total variance

D1=int(g1^2*f1,x1,-1,1);
D2=int(g2^2*f2,x2,-1,1);
D3=int(g3^2*f3,x3,-1,1);

D12=int(int(g12^2*f12,x1,-1,1),x2,-1,1);
D13=int(int(g13^2*f13,x1,-1,1),x3,-1,1);
D23=int(int(g23^2*f23,x2,-1,1),x3,-1,1);

D123=int(int(int(g123^2*f123,x1,-1,1),x2,-1,1),x3,-1,1);
% ----------------------------------------------------------------------------------------------------------------------



% ----------------------------------------------------------------------------------------------------------------------
% CONTRIBUTIONS TO TOTAL VARIANCE:
% ----------------------------------------------------------------------------------------------------------------------
S1=D1/Dg; % contribution to total variance due to x1
S2=D2/Dg; % contribution to total variance due to x2
S3=D3/Dg; % contribution to total variance due to x3

S12=D12/Dg; % contribution to total variance due to the interaction between x1 and x2
S13=D13/Dg; % contribution to total variance due to the interaction between x1 and x3
S23=D23/Dg; % contribution to total variance due to the interaction between x2 and x3

S123=D123/Dg; % contribution to total variance due to the interaction between x1, x2 and x3
% ----------------------------------------------------------------------------------------------------------------------



% ----------------------------------------------------------------------------------------------------------------------
% SUMMARY:
% ----------------------------------------------------------------------------------------------------------------------
fprintf('Summary:\n\n')
fprintf('Total mean: %.6f.\n',g0)
fprintf('Total variance: %.4f.\n',Dg)
fprintf('\n')
fprintf('Random variable x1 contributes to the total variance response by: %.2f%%.\n',S1*100)
fprintf('Random variable x2 contributes to the total variance response by: %.2f%%.\n',S2*100)
fprintf('Random variable x3 contributes to the total variance response by: %.2f%%.\n',S3*100)
fprintf('Random variables x1 and x2 contribute to the total variance response by: %.2f%%.\n',S12*100)
fprintf('Random variables x1 and x3 contribute to the total variance response by: %.2f%%.\n',S13*100)
fprintf('Random variables x2 and x3 contribute to the total variance response by: %.2f%%.\n',S23*100)
fprintf('Random variables x1, x2 and x3 contribute to the total variance response by: %.2f%%.\n',S123*100)
% ----------------------------------------------------------------------------------------------------------------------



% ----------------------------------------------------------------------------------------------------------------------
% FUNCTIONS:
% ----------------------------------------------------------------------------------------------------------------------
function g=getStructuralResponse(x1,x2,x3)
p1=legendreP(2,x1)+jacobiP(3,3,5,x1); % a univariate function defined on x1 = [-1,1]
p2=legendreP(5,x2)+jacobiP(2,2,7,x2); % a univariate function defined on x2 = [-1,1]
p3=legendreP(3,x3)+jacobiP(2,4,1,x3); % a univariate function defined on x3 = [-1,1]

% Structural response of interest (notice that this g is a function of x1, x2 and x3):
g=p1*p2*p3/1e3; % a trivariate function defined on x1 x x2 x x3 = [-1,1]^3
end
% ----------------------------------------------------------------------------------------------------------------------
