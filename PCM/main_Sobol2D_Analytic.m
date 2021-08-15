% Hugo Esquivel, 2021
% -

clearvars; close all; clc;


% ----------------------------------------------------------------------------------------------------------------------
% INPUT:
% ----------------------------------------------------------------------------------------------------------------------
syms x1 x2 % random variables x1 and x2

g=getStructuralResponse(x1,x2); % structural response g
plotStructuralResponse(g) % plots structural response g

probabilityDistribution='Uniform'; % options: Uniform or Beta

switch probabilityDistribution
    case 'Uniform'
        f1=1/2; % probability density of x1 (assuming that x1 is a uniformly distributed random variable in [-1,1])
        f2=1/2; % probability density of x2 (assuming that x2 is a uniformly distributed random variable in [-1,1])
        
    case 'Beta'
        % probability density of x1 (assuming that x1 is a beta distributed random variable in [-1,1])
        alpha1=2;
        beta1=5;
        f1=(x1+1)^(alpha1-1)*(1-x1)^(beta1-1)/(2^(alpha1+beta1-1)*beta(alpha1,beta1));
        
        % probability density of x2 (assuming that x2 is a beta distributed random variable in [-1,1])
        alpha2=2;
        beta2=5;
        f2=(x2+1)^(alpha2-1)*(1-x2)^(beta2-1)/(2^(alpha2+beta2-1)*beta(alpha2,beta2));
end

f12=f1*f2; % joint probability density of x1 and x2
% ----------------------------------------------------------------------------------------------------------------------



% ----------------------------------------------------------------------------------------------------------------------
% DECOMPOSITION OF STRUCTURAL RESPONSE:
% ----------------------------------------------------------------------------------------------------------------------
g0=int(int(g*f12,x1,-1,1),x2,-1,1); % total mean

g1=int(g*f2,x2,-1,1)-g0;
g2=int(g*f1,x1,-1,1)-g0;

g12=g-(g0+g1+g2);
% ----------------------------------------------------------------------------------------------------------------------



% ----------------------------------------------------------------------------------------------------------------------
% CHECKING ORTHOGONALITY OF DECOMPOSITION:
% ----------------------------------------------------------------------------------------------------------------------
check1=int(g1*f1,x1,-1,1)==0;
check2=int(g2*f2,x2,-1,1)==0;

check12=int(int(g12*f12,x1,-1,1),x2,-1,1)==0;

if ~(check1 && check2 && check12)
    error('Orthogonality check did not pass.')
end
% ----------------------------------------------------------------------------------------------------------------------



% ----------------------------------------------------------------------------------------------------------------------
% DECOMPOSITION OF TOTAL VARIANCE:
% ----------------------------------------------------------------------------------------------------------------------
Dg=int(int(g^2*f12,x1,-1,1),x2,-1,1)-g0^2; % total variance

D1=int(g1^2*f1,x1,-1,1);
D2=int(g2^2*f2,x2,-1,1);

D12=int(int(g12^2*f12,x1,-1,1),x2,-1,1);
% ----------------------------------------------------------------------------------------------------------------------



% ----------------------------------------------------------------------------------------------------------------------
% CONTRIBUTIONS TO TOTAL VARIANCE:
% ----------------------------------------------------------------------------------------------------------------------
S1=D1/Dg; % contribution to total variance due to x1
S2=D2/Dg; % contribution to total variance due to x2

S12=D12/Dg; % contribution to total variance due to the interaction between x1 and x2
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
fprintf('Random variables x1 and x2 contribute to the total variance response by: %.2f%%.\n',S12*100)
% ----------------------------------------------------------------------------------------------------------------------



% ----------------------------------------------------------------------------------------------------------------------
% FUNCTIONS:
% ----------------------------------------------------------------------------------------------------------------------
function plotStructuralResponse(g)
g=matlabFunction(g);

numPoints=50;

x1=linspace(-1,1,numPoints); % x1 = [-1,1]
x2=linspace(-1,1,numPoints); % x2 = [-1,1]

[X1,X2]=meshgrid(x1,x2);
G=g(X1,X2);

surf(X1,X2,G)
colormap(jet)
title('\bfseries{Structural Response, $g$}','Interpreter','LaTeX')
xlabel('$x_1$','Interpreter','LaTeX')
ylabel('$x_2$','Interpreter','LaTeX')
zlabel('$g(x_1,x_2)$','Interpreter','LaTeX')
end

function g=getStructuralResponse(x1,x2)
p1=legendreP(2,x1)+jacobiP(3,3,5,x1); % a univariate function defined on x1 = [-1,1]
p2=legendreP(5,x2)+jacobiP(2,2,7,x2); % a univariate function defined on x2 = [-1,1]

% Structural response of interest (notice that this g is a function of x1 and x2):
g=p1*p2/1e2; % a bivariate function defined on x1 x x2 = [-1,1]^2
end
% ----------------------------------------------------------------------------------------------------------------------
