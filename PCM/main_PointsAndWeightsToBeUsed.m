% Hugo Esquivel, 2021
% -

clearvars; close all; clc;

% ----------------------------------------------------------------------------------------------------------------------
% INPUT:
% ----------------------------------------------------------------------------------------------------------------------
numRandomVariables=5; % assumed to be normally distributed (all of them)

numPointsPerDimension=11; % number of quadrature points to be used on each random dimension

fileName=sprintf('Hermite-%d.txt',numPointsPerDimension); % file containing Hermite's quadrature information
% ----------------------------------------------------------------------------------------------------------------------


% ----------------------------------------------------------------------------------------------------------------------
% BODY:
% ----------------------------------------------------------------------------------------------------------------------
d=load(fileName,'r');

z1=d(:,1)*sqrt(2); % probabilists' version of Hermite quadrature points
w1=d(:,2)/sqrt(pi); % probabilists' version of Hermite quadrature weights

numPoints=numPointsPerDimension^numRandomVariables; % because a full tensor product will be used below

z=zeros(numPoints,numRandomVariables); % quadrature points
w=zeros(numPoints,1); % quadrature weights

k=0;

for i5=1:numPointsPerDimension
    for i4=1:numPointsPerDimension
        for i3=1:numPointsPerDimension
            for i2=1:numPointsPerDimension
                for i1=1:numPointsPerDimension
                    k=k+1;
                    
                    z(k,1)=z1(i1);
                    z(k,2)=z1(i2);
                    z(k,3)=z1(i3);
                    z(k,4)=z1(i4);
                    z(k,5)=z1(i5);
                    
                    w(k)=w1(i1)*w1(i2)*w1(i3)*w1(i4)*w1(i5);
                end
            end
        end
    end
end

% now let's ignore those tiny weights that do not contribute significantly to the computation of the mean and variance:
tol=1e-10; % any weight less than this tolerance will be ignored
idx=w>tol; % points to be saved are denoted with 1

% thus, realizations to be used for seismic simulations:
p=z(idx,:); % realizations (note that they are normalized wrt the standard normal distribution^)
v=w(idx); % corresponding weights

fprintf('Number of points to be used: %d (instead of %d points).\n',length(v),length(w))
fprintf('Thus, computational cost is reduced by %.1f%%.\n',(1-length(v)/length(w))*100)
fprintf('Precision of resulting weights is %.14f%%.\n',sum(v)*100)
fprintf('\n')

% Footnote:
% ^ to denormalize, take each column of p, say the ith column p(:,i), and perform the following operation:
%        ith_sigma * p(:,i) + ith_mu,
%   where ith_sigma and ith_mu are the standard deviation and the mean of the ith random variable, respectively.
% ----------------------------------------------------------------------------------------------------------------------


% ----------------------------------------------------------------------------------------------------------------------
% A CHECK:
% ----------------------------------------------------------------------------------------------------------------------
fprintf('Checking this approach with a polynomial function:\n')

syms x1 x2 x3 x4 x5

% analytical evaluation of mean and variance:
fun=(x2+1)*(x3-2)*(x5-4)+(x1-1)^2*(x3+2)*(x4+1); % polynomial function
pdf=(1/sqrt(2*pi))^5*exp(-1/2*(x1^2+x2^2+x3^2+x4^2+x5^2)); % joint pdf

mean1=double(int(int(int(int(int(fun*pdf,x5,-Inf,Inf),x4,-Inf,Inf),x3,-Inf,Inf),x2,-Inf,Inf),x1,-Inf,Inf)); % exact mean
var1=double(int(int(int(int(int((fun-mean1)^2*pdf,x5,-Inf,Inf),x4,-Inf,Inf),x3,-Inf,Inf),x2,-Inf,Inf),x1,-Inf,Inf)); % exact variance

% numerical evaluation of mean and variance:
funnum=matlabFunction(fun,'Vars',[x1,x2,x3,x4,x5]); % polynomial function
pdfnum=matlabFunction(pdf,'Vars',[x1,x2,x3,x4,x5]); % joint pdf

mean2=sum(funnum(p(:,1),p(:,2),p(:,3),p(:,4),p(:,5)).*v); % approximate mean
var2=sum((funnum(p(:,1),p(:,2),p(:,3),p(:,4),p(:,5))-mean2).^2.*v); % approximate variance

fprintf('- Exact solution:       [%f, %f].\n',mean1,var1)
fprintf('- Approximate solution: [%f, %f].\n',mean2,var2)
fprintf('- Error:                [%f%%, %f%%].\n',abs(1-mean2/mean1)*100,abs(1-var2/var1)*100)
fprintf('\n')
% ----------------------------------------------------------------------------------------------------------------------
