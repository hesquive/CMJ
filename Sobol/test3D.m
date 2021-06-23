% Hugo Esquivel, 2021
% -

clearvars; close all; clc;

% ----------------------------------------------------------------------------------------------------------------------
% Input:
% ----------------------------------------------------------------------------------------------------------------------
syms x1 x2 x3 % random variables x1, x2 and x3 (these can be understood as fc, fy and M, respectively)

f=getStructuralResponse(x1,x2,x3); % structural response f (this can be understood as the roof drift ratio)

probabilityDistribution='Uniform'; % options: Uniform or Beta

switch probabilityDistribution
    case 'Uniform'
        pdf1=1/2; % probability density of x1 (assuming that x1 is a uniformly distributed random variable in [-1,1])
        pdf2=1/2; % probability density of x2 (assuming that x2 is a uniformly distributed random variable in [-1,1])
        pdf3=1/2; % probability density of x3 (assuming that x3 is a uniformly distributed random variable in [-1,1])
        
    case 'Beta'
        % probability density of x1 (assuming that x1 is a beta distributed random variable in [-1,1])
        alphax1=2;
        betax1=5;
        pdf1=(x1+1)^(alphax1-1)*(1-x1)^(betax1-1)/(2^(alphax1+betax1-1)*beta(alphax1,betax1));
        
        % probability density of x2 (assuming that x2 is a beta distributed random variable in [-1,1])
        alphax2=5;
        betax2=1;
        pdf2=(x2+1)^(alphax2-1)*(1-x2)^(betax2-1)/(2^(alphax2+betax2-1)*beta(alphax2,betax2));
        
        % probability density of x3 (assuming that x3 is a beta distributed random variable in [-1,1])
        alphax3=2;
        betax3=2;
        pdf3=(x3+1)^(alphax3-1)*(1-x3)^(betax3-1)/(2^(alphax3+betax3-1)*beta(alphax3,betax3));
end

pdf12=pdf1*pdf2; % joint probability density of x1 and x2
pdf13=pdf1*pdf3; % joint probability density of x1 and x3
pdf23=pdf2*pdf3; % joint probability density of x2 and x3

pdf=pdf1*pdf2*pdf3; % joint probability density of x1, x2 and x3
% ----------------------------------------------------------------------------------------------------------------------


% ----------------------------------------------------------------------------------------------------------------------
% Decomposition of structural response:
% ----------------------------------------------------------------------------------------------------------------------
f0=int(int(int(f*pdf,x1,-1,1),x2,-1,1),x3,-1,1); % total mean

f1=int(int(f*pdf23,x2,-1,1),x3,-1,1)-f0; 
f2=int(int(f*pdf13,x1,-1,1),x3,-1,1)-f0;
f3=int(int(f*pdf12,x1,-1,1),x2,-1,1)-f0;

f12=int(f*pdf3,x3,-1,1)-(f0+f1+f2);
f13=int(f*pdf2,x2,-1,1)-(f0+f1+f3);
f23=int(f*pdf1,x1,-1,1)-(f0+f2+f3);

f123=f-(f0+f1+f2+f3+f12+f13+f23);
% ----------------------------------------------------------------------------------------------------------------------


% ----------------------------------------------------------------------------------------------------------------------
% Checking orthogonality of decomposition:
% ----------------------------------------------------------------------------------------------------------------------
check1=int(f1*pdf1,x1,-1,1)==0;
check2=int(f2*pdf2,x2,-1,1)==0;
check3=int(f3*pdf3,x3,-1,1)==0;

check12=int(int(f12*pdf12,x1,-1,1),x2,-1,1)==0;
check13=int(int(f13*pdf13,x1,-1,1),x3,-1,1)==0;
check23=int(int(f23*pdf23,x2,-1,1),x3,-1,1)==0;

check123=int(int(int(f123*pdf,x1,-1,1),x2,-1,1),x3,-1,1)==0;

if ~(check1 && check2 && check3 && check12 && check13 && check23 && check123)
    error('Orthogonality check did not pass.')
end
% ----------------------------------------------------------------------------------------------------------------------


% ----------------------------------------------------------------------------------------------------------------------
% Decomposition of total variance:
% ----------------------------------------------------------------------------------------------------------------------
Df=int(int(int(f^2*pdf,x1,-1,1),x2,-1,1),x3,-1,1)-f0^2; % total variance

D1=int(f1^2*pdf1,x1,-1,1);
D2=int(f2^2*pdf2,x2,-1,1);
D3=int(f3^2*pdf3,x3,-1,1);

D12=int(int(f12^2*pdf12,x1,-1,1),x2,-1,1);
D13=int(int(f13^2*pdf13,x1,-1,1),x3,-1,1);
D23=int(int(f23^2*pdf23,x2,-1,1),x3,-1,1);

D123=int(int(int(f123^2*pdf,x1,-1,1),x2,-1,1),x3,-1,1);
% ----------------------------------------------------------------------------------------------------------------------


% ----------------------------------------------------------------------------------------------------------------------
% Contributions to total variance:
% ----------------------------------------------------------------------------------------------------------------------
S1=D1/Df; % contribution to total variance due to x1
S2=D2/Df; % contribution to total variance due to x2
S3=D3/Df; % contribution to total variance due to x3

S12=D12/Df; % contribution to total variance due to the interaction between x1 and x2
S13=D13/Df; % contribution to total variance due to the interaction between x1 and x3
S23=D23/Df; % contribution to total variance due to the interaction between x2 and x3

S123=D123/Df; % contribution to total variance due to the interaction between x1, x2 and x3

disp(double([S1;S2;S3;S12;S13;S23;S123]))
% ----------------------------------------------------------------------------------------------------------------------


% ----------------------------------------------------------------------------------------------------------------------
% Summary:
% ----------------------------------------------------------------------------------------------------------------------
fprintf('Summary:\n\n')
fprintf('Total mean: %.6f.\n',f0)
fprintf('Total variance: %.4f.\n',Df)
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
% Functions:
% ----------------------------------------------------------------------------------------------------------------------
function f=getStructuralResponse(x1,x2,x3)
p1=legendreP(2,x1)+jacobiP(3,3,5,x1); % a univariate function defined on x1 = [-1,1]
p2=legendreP(5,x2)+jacobiP(2,2,7,x2); % a univariate function defined on x2 = [-1,1]
p3=legendreP(3,x3)+jacobiP(2,4,1,x3); % a univariate function defined on x3 = [-1,1]

% Structural response of interest (notice that this f is a function of x1, x2 and x3):
f=p1*p2*p3/1e3; % a trivariate function defined on x1 x x2 x x3 = [-1,1]^3
end
% ----------------------------------------------------------------------------------------------------------------------
