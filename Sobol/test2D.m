% Hugo Esquivel, 2021
% -

clearvars; close all; clc;

% ----------------------------------------------------------------------------------------------------------------------
% Input:
% ----------------------------------------------------------------------------------------------------------------------
syms x1 x2 % random variables x1 and x2 (these can be understood as fc and fy, respectively)

f=getStructuralResponse(x1,x2); % structural response f (this can be understood as the roof drift ratio)
plotStructuralResponse(f) % plots structural response f

probabilityDistribution='Uniform'; % options: Uniform or Beta

switch probabilityDistribution
    case 'Uniform'
        pdf1=1/2; % probability density of x1 (assuming that x1 is a uniformly distributed random variable in [-1,1])
        pdf2=1/2; % probability density of x2 (assuming that x2 is a uniformly distributed random variable in [-1,1])
        
    case 'Beta'
        % probability density of x1 (assuming that x1 is a beta distributed random variable in [-1,1])
        alphax1=2;
        betax1=5;
        pdf1=(x1+1)^(alphax1-1)*(1-x1)^(betax1-1)/(2^(alphax1+betax1-1)*beta(alphax1,betax1));
        
        % probability density of x2 (assuming that x2 is a beta distributed random variable in [-1,1])
        alphax2=5;
        betax2=1;
        pdf2=(x2+1)^(alphax2-1)*(1-x2)^(betax2-1)/(2^(alphax2+betax2-1)*beta(alphax2,betax2));
end

pdf12=pdf1*pdf2; % joint probability density of x1 and x2
% ----------------------------------------------------------------------------------------------------------------------


% ----------------------------------------------------------------------------------------------------------------------
% Decomposition of structural response:
% ----------------------------------------------------------------------------------------------------------------------
f0=int(int(f*pdf12,x1,-1,1),x2,-1,1); % total mean

f1=int(f*pdf2,x2,-1,1)-f0;
f2=int(f*pdf1,x1,-1,1)-f0;

f12=f-(f0+f1+f2);
% ----------------------------------------------------------------------------------------------------------------------


% ----------------------------------------------------------------------------------------------------------------------
% Checking orthogonality of decomposition:
% ----------------------------------------------------------------------------------------------------------------------
check1=int(f1*pdf1,x1,-1,1)==0;
check2=int(f2*pdf2,x2,-1,1)==0;

check12=int(int(f12*pdf12,x1,-1,1),x2,-1,1)==0;

if ~(check1 && check2 && check12)
    error('Orthogonality check did not pass.')
end
% ----------------------------------------------------------------------------------------------------------------------


% ----------------------------------------------------------------------------------------------------------------------
% Decomposition of total variance:
% ----------------------------------------------------------------------------------------------------------------------
Df=int(int(f^2*pdf12,x1,-1,1),x2,-1,1)-f0^2; % total variance

D1=int(f1^2*pdf1,x1,-1,1);
D2=int(f2^2*pdf2,x2,-1,1);

D12=int(int(f12^2*pdf12,x1,-1,1),x2,-1,1);
% ----------------------------------------------------------------------------------------------------------------------


% ----------------------------------------------------------------------------------------------------------------------
% Contributions to total variance:
% ----------------------------------------------------------------------------------------------------------------------
S1=D1/Df; % contribution to total variance due to x1
S2=D2/Df; % contribution to total variance due to x2

S12=D12/Df; % contribution to total variance due to the interaction between x1 and x2

disp(double([S1;S2;S12]))
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
fprintf('Random variables x1 and x2 contribute to the total variance response by: %.2f%%.\n',S12*100)
% ----------------------------------------------------------------------------------------------------------------------



% ----------------------------------------------------------------------------------------------------------------------
% Functions:
% ----------------------------------------------------------------------------------------------------------------------
function plotStructuralResponse(f)
f=matlabFunction(f);

numPoints=50;

x1=linspace(-1,1,numPoints); % x1 = [-1,1]
x2=linspace(-1,1,numPoints); % x2 = [-1,1]

[x1grid,x2grid]=meshgrid(x1,x2);
fgrid=f(x1grid,x2grid);

surf(x1grid,x2grid,fgrid)
colormap(jet)
title('\bfseries{Structural Response, $f$}','Interpreter','LaTeX')
xlabel('$x_1$','Interpreter','LaTeX')
ylabel('$x_2$','Interpreter','LaTeX')
zlabel('$f(x_1,x_2)$','Interpreter','LaTeX')
end

function f=getStructuralResponse(x1,x2)
p1=legendreP(2,x1)+jacobiP(3,3,5,x1); % a univariate function defined on x1 = [-1,1]
p2=legendreP(5,x2)+jacobiP(2,2,7,x2); % a univariate function defined on x2 = [-1,1]

% Structural response of interest (notice that this f is a function of x1 and x2):
f=p1*p2/1e2; % a bivariate function defined on x1 x x2 = [-1,1]^2
end
% ----------------------------------------------------------------------------------------------------------------------
