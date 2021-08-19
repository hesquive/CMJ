% Hugo Esquivel, 2021
% -

clearvars; close all; clc;


% ----------------------------------------------------------------------------------------------------------------------
% INPUT:
% ----------------------------------------------------------------------------------------------------------------------
numRandomVariables=5; % do not modify... they are assumed here to be identically distributed (all of them)
numPoints1=11; % number of quadrature points per random dimension

probabilityDistribution='Normal'; % options: Uniform, Beta or Normal
% ----------------------------------------------------------------------------------------------------------------------



% ----------------------------------------------------------------------------------------------------------------------
% BODY:
% ----------------------------------------------------------------------------------------------------------------------
if strcmpi(probabilityDistribution,'Uniform')
    fileName=sprintf('Legendre-%d.txt',numPoints1); % file containing Legendre's quadrature information
    
    d=load(fileName,'r');
    
    z1=d(:,1); % Legendre quadrature points per random dimension
    y1=d(:,2)/2; % Legendre quadrature weights per random dimension (probabilists' version)
    
elseif strcmpi(probabilityDistribution,'Beta')
    betafun=@(x,y) beta(x,y);
    
    alpha=2; % do not modify
    beta=5; % do not modify
    
    fileName=sprintf('Jacobi_%d_%d-%d.txt',beta-1,alpha-1,numPoints1); % file containing Jabobi(4,1)'s quadrature information
    
    d=load(fileName,'r');
    
    z1=d(:,1); % Jacobi(4,1) quadrature points per random dimension
    y1=d(:,2)/(2^(alpha+beta-1)*betafun(alpha,beta)); % Jacobi(4,1) quadrature weights per random dimension (probabilists' version)
    
elseif strcmpi(probabilityDistribution,'Normal')
    fileName=sprintf('Hermite-%d.txt',numPoints1); % file containing Hermite's quadrature information
    
    d=load(fileName,'r');
    
    z1=d(:,1)*sqrt(2); %  Hermite quadrature points per random dimension (probabilists' version)
    y1=d(:,2)/sqrt(pi); % Hermite quadrature weights per random dimension (probabilists' version)
end

numPoints=numPoints1^numRandomVariables; % because a full tensor product will be used below

z=zeros(numPoints,numRandomVariables); % quadrature points
y=zeros(numPoints,1); % quadrature weights

k=0;

for i5=1:numPoints1
    for i4=1:numPoints1
        for i3=1:numPoints1
            for i2=1:numPoints1
                for i1=1:numPoints1
                    k=k+1;
                    
                    z(k,1)=z1(i1);
                    z(k,2)=z1(i2);
                    z(k,3)=z1(i3);
                    z(k,4)=z1(i4);
                    z(k,5)=z1(i5);
                    
                    y(k)=y1(i1)*y1(i2)*y1(i3)*y1(i4)*y1(i5);
                end
            end
        end
    end
end

% now let's ignore those tiny weights that do not contribute significantly to the computation of the mean and variance:
if strcmpi(probabilityDistribution,'Uniform')
    tol=1e-10; % any weight less than this tolerance will be ignored
    
elseif strcmpi(probabilityDistribution,'Beta')
    tol=1e-10; % any weight less than this tolerance will be ignored
    
elseif strcmpi(probabilityDistribution,'Normal')
    tol=1e-10; % any weight less than this tolerance will be ignored
end

idx=y>tol; % points to be saved are denoted with 1

% thus, realizations to be used for seismic simulations:
p=z(idx,:); % realizations (note that they are normalized)
v=y(idx); % corresponding weights

fprintf('Number of points to be used: %d (instead of %d points).\n',length(v),length(y))
fprintf('Thus, computational cost is reduced by %.1f%%.\n',(1-length(v)/length(y))*100)
fprintf('Precision of resulting weights is %.14f%%.\n',sum(v)*100)
fprintf('\n')
% ----------------------------------------------------------------------------------------------------------------------



% ----------------------------------------------------------------------------------------------------------------------
% GENERATE FILE WITH POINTS AND WEIGHTS TO BE USED:
% ----------------------------------------------------------------------------------------------------------------------
fileName=sprintf('pointsAndWeightsToBeUsed_%s.txt',probabilityDistribution);

fid=fopen(fileName,'w');

for i=1:numPoints
    if idx(i)==1
        fprintf(fid,'%d\t%20.14e\t%21.14e\t%21.14e\t%21.14e\t%21.14e\t%21.14e\n',...
            i,y(i),z(i,1),z(i,2),z(i,3),z(i,4),z(i,5)); % id, weight, point_x1, point_x2, point_x3, point_x4, point_x5
    end
end

fclose(fid);
% ----------------------------------------------------------------------------------------------------------------------



% ----------------------------------------------------------------------------------------------------------------------
% A CHECK:
% ----------------------------------------------------------------------------------------------------------------------
if strcmpi(probabilityDistribution,'Normal')
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
end
% ----------------------------------------------------------------------------------------------------------------------
