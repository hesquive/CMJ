function x=getRealizationsFromRandomVariables(numRealizations)
% Hugo Esquivel, 2021
% -
% Remark: Assuming that the random variables are mutually independent and uniformly distributed in [-1,1]:

numRandomVariables=2;

x=rand(numRealizations,numRandomVariables);
x=2*x-1; % coordinate transformation to go from [0,1] to [-1,1]
end
