function f=getStructuralResponse()
% Hugo Esquivel, 2021
% -

syms x1 x2

p1=legendreP(2,x1)+jacobiP(3,3,5,x1); % a univariate function defined on x1 = [-1,1]
p2=legendreP(5,x2)+jacobiP(2,2,7,x2); % a univariate function defined on x2 = [-1,1]

% Structural response of interest (notice that this f is a function of x1 and x2):
f=p1*p2/1e2; % a bivariate function defined on x1 x x2 = [-1,1]^2
end
