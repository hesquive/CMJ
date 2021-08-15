% Hugo Esquivel, 2021
% -

clearvars; close all; clc;


% ----------------------------------------------------------------------------------------------------------------------
% INPUT:
% ----------------------------------------------------------------------------------------------------------------------
numRandomVariables=3; % do not modify
numPoints1=11; % number of quadrature points per random dimension

probabilityDistribution='Normal'; % options: Uniform, Beta or Normal
% ----------------------------------------------------------------------------------------------------------------------



% ----------------------------------------------------------------------------------------------------------------------
% QUADRATURE RULE:
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

X1=zeros(numPoints1,numPoints1,numPoints1);
X2=zeros(numPoints1,numPoints1,numPoints1);
X3=zeros(numPoints1,numPoints1,numPoints1);

for i3=1:numPoints1
    for i2=1:numPoints1
        for i1=1:numPoints1
            X1(i1,i2,i3)=z1(i1);
            X2(i1,i2,i3)=z1(i2);
            X3(i1,i2,i3)=z1(i3);
        end
    end
end

W1=zeros(numPoints1,1,1);
W2=zeros(1,numPoints1,1);
W3=zeros(1,1,numPoints1);

for i1=1:numPoints1
    W1(i1,:,:)=y1(i1);
end

for i2=1:numPoints1
    W2(:,i2,:)=y1(i2);
end

for i3=1:numPoints1
    W3(:,:,i3)=y1(i3);
end

W12=zeros(numPoints1,numPoints1,1);
W13=zeros(numPoints1,1,numPoints1);
W23=zeros(1,numPoints1,numPoints1);

for i2=1:numPoints1
    for i1=1:numPoints1
        W12(i1,i2,:)=y1(i1)*y1(i2);
    end
end

for i3=1:numPoints1
    for i1=1:numPoints1
        W13(i1,:,i3)=y1(i1)*y1(i3);
    end
end

for i3=1:numPoints1
    for i2=1:numPoints1
        W23(:,i2,i3)=y1(i2)*y1(i3);
    end
end

W123=zeros(numPoints1,numPoints1,numPoints1);

for i3=1:numPoints1
    for i2=1:numPoints1
        for i1=1:numPoints1
            W123(i1,i2,i3)=y1(i1)*y1(i2)*y1(i3);
        end
    end
end

G=getStructuralResponse(X1,X2,X3); % structural response G
% ----------------------------------------------------------------------------------------------------------------------



% ----------------------------------------------------------------------------------------------------------------------
% DECOMPOSITION OF STRUCTURAL RESPONSE:
% ----------------------------------------------------------------------------------------------------------------------
G0=zeros(numPoints1,numPoints1,numPoints1);

G0(:,:,:)=sum(sum(sum(G.*W123))); % total mean

G1=zeros(numPoints1,numPoints1,numPoints1);
G2=zeros(numPoints1,numPoints1,numPoints1);
G3=zeros(numPoints1,numPoints1,numPoints1);

for i1=1:numPoints1
    G1(i1,:,:)=sum(sum(G(i1,:,:).*W23))-G0(i1,:,:);
end

for i2=1:numPoints1
    G2(:,i2,:)=sum(sum(G(:,i2,:).*W13))-G0(:,i2,:);
end

for i3=1:numPoints1
    G3(:,:,i3)=sum(sum(G(:,:,i3).*W12))-G0(:,:,i3);
end

G12=zeros(numPoints1,numPoints1,numPoints1);
G13=zeros(numPoints1,numPoints1,numPoints1);
G23=zeros(numPoints1,numPoints1,numPoints1);

for i2=1:numPoints1
    for i1=1:numPoints1
        G12(i1,i2,:)=sum(G(i1,i2,:).*W3)-(G0(i1,i2,:)+G1(i1,i2,:)+G2(i1,i2,:));
    end
end

for i3=1:numPoints1
    for i1=1:numPoints1
        G13(i1,:,i3)=sum(G(i1,:,i3).*W2)-(G0(i1,:,i3)+G1(i1,:,i3)+G3(i1,:,i3));
    end
end

for i3=1:numPoints1
    for i2=1:numPoints1
        G23(:,i2,i3)=sum(G(:,i2,i3).*W1)-(G0(:,i2,i3)+G2(:,i2,i3)+G3(:,i2,i3));
    end
end

G123=zeros(numPoints1,numPoints1,numPoints1);

G123(:,:,:)=G-(G0+G1+G2+G3+G12+G13+G23);
% ----------------------------------------------------------------------------------------------------------------------



% ----------------------------------------------------------------------------------------------------------------------
% DECOMPOSITION OF TOTAL VARIANCE:
% ----------------------------------------------------------------------------------------------------------------------
DG=sum(sum(sum(G.^2.*W123)))-G0(1,1,1)^2; % total variance

D1=sum(G1(:,1,1).^2.*W1);
D2=sum(G2(1,:,1).^2.*W2);
D3=sum(G3(1,1,:).^2.*W3);

D12=sum(sum(G12(:,:,1).^2.*W12));
D13=sum(sum(G13(:,1,:).^2.*W13));
D23=sum(sum(G23(1,:,:).^2.*W23));

D123=sum(sum(sum(G123.^2.*W123)));
% ----------------------------------------------------------------------------------------------------------------------



% ----------------------------------------------------------------------------------------------------------------------
% CONTRIBUTIONS TO TOTAL VARIANCE:
% ----------------------------------------------------------------------------------------------------------------------
S1=D1/DG; % contribution to total variance due to x1
S2=D2/DG; % contribution to total variance due to x2
S3=D3/DG; % contribution to total variance due to x3

S12=D12/DG; % contribution to total variance due to the interaction between x1 and x2
S13=D13/DG; % contribution to total variance due to the interaction between x1 and x3
S23=D23/DG; % contribution to total variance due to the interaction between x2 and x3

S123=D123/DG; % contribution to total variance due to the interaction between x1, x2 and x3
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
function G=getStructuralResponse(X1,X2,X3)
syms x1 x2 x3

p1=legendreP(2,x1)+jacobiP(3,3,5,x1); % a univariate function defined on x1 = [-1,1]
p2=legendreP(5,x2)+jacobiP(2,2,7,x2); % a univariate function defined on x2 = [-1,1]
p3=legendreP(3,x3)+jacobiP(2,4,1,x3); % a univariate function defined on x3 = [-1,1]

% Structural response of interest (notice that this g is a function of x1, x2 and x3):
g=p1*p2*p3/1e3; % a trivariate function defined on x1 x x2 x x3 = [-1,1]^3

g=matlabFunction(g,'Vars',[x1,x2,x3]);

G=g(X1,X2,X3);
end
% ----------------------------------------------------------------------------------------------------------------------
