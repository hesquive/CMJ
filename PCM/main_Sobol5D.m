% Hugo Esquivel, 2021
% -

clearvars; close all; clc;


% ----------------------------------------------------------------------------------------------------------------------
% INPUT:
% ----------------------------------------------------------------------------------------------------------------------
numRandomVariables=5; % do not modify
numPoints1=11; % number of quadrature points per random dimension

probabilityDistribution='Uniform'; % options: Uniform, Beta or Normal
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

X1=zeros(numPoints1,numPoints1,numPoints1,numPoints1,numPoints1);
X2=zeros(numPoints1,numPoints1,numPoints1,numPoints1,numPoints1);
X3=zeros(numPoints1,numPoints1,numPoints1,numPoints1,numPoints1);
X4=zeros(numPoints1,numPoints1,numPoints1,numPoints1,numPoints1);
X5=zeros(numPoints1,numPoints1,numPoints1,numPoints1,numPoints1);

for i5=1:numPoints1
    for i4=1:numPoints1
        for i3=1:numPoints1
            for i2=1:numPoints1
                for i1=1:numPoints1
                    X1(i1,i2,i3,i4,i5)=z1(i1);
                    X2(i1,i2,i3,i4,i5)=z1(i2);
                    X3(i1,i2,i3,i4,i5)=z1(i3);
                    X4(i1,i2,i3,i4,i5)=z1(i4);
                    X5(i1,i2,i3,i4,i5)=z1(i5);
                end
            end
        end
    end
end

W1=zeros(numPoints1,1,1,1,1);
W2=zeros(1,numPoints1,1,1,1);
W3=zeros(1,1,numPoints1,1,1);
W4=zeros(1,1,1,numPoints1,1);
W5=zeros(1,1,1,1,numPoints1);

for i1=1:numPoints1
    W1(i1,:,:,:,:)=y1(i1);
end

for i2=1:numPoints1
    W2(:,i2,:,:,:)=y1(i2);
end

for i3=1:numPoints1
    W3(:,:,i3,:,:)=y1(i3);
end

for i4=1:numPoints1
    W4(:,:,:,i4,:)=y1(i4);
end

for i5=1:numPoints1
    W5(:,:,:,:,i5)=y1(i5);
end

W12=zeros(numPoints1,numPoints1,1,1,1);
W13=zeros(numPoints1,1,numPoints1,1,1);
W14=zeros(numPoints1,1,1,numPoints1,1);
W15=zeros(numPoints1,1,1,1,numPoints1);
W23=zeros(1,numPoints1,numPoints1,1,1);
W24=zeros(1,numPoints1,1,numPoints1,1);
W25=zeros(1,numPoints1,1,1,numPoints1);
W34=zeros(1,1,numPoints1,numPoints1,1);
W35=zeros(1,1,numPoints1,1,numPoints1);
W45=zeros(1,1,1,numPoints1,numPoints1);

for i2=1:numPoints1
    for i1=1:numPoints1
        W12(i1,i2,:,:,:)=y1(i1)*y1(i2);
    end
end

for i3=1:numPoints1
    for i1=1:numPoints1
        W13(i1,:,i3,:,:)=y1(i1)*y1(i3);
    end
end

for i4=1:numPoints1
    for i1=1:numPoints1
        W14(i1,:,:,i4,:)=y1(i1)*y1(i4);
    end
end

for i5=1:numPoints1
    for i1=1:numPoints1
        W15(i1,:,:,:,i5)=y1(i1)*y1(i5);
    end
end

for i3=1:numPoints1
    for i2=1:numPoints1
        W23(:,i2,i3,:,:)=y1(i2)*y1(i3);
    end
end

for i4=1:numPoints1
    for i2=1:numPoints1
        W24(:,i2,:,i4,:)=y1(i2)*y1(i4);
    end
end

for i5=1:numPoints1
    for i2=1:numPoints1
        W25(:,i2,:,:,i5)=y1(i2)*y1(i5);
    end
end

for i4=1:numPoints1
    for i3=1:numPoints1
        W34(:,:,i3,i4,:)=y1(i3)*y1(i4);
    end
end

for i5=1:numPoints1
    for i3=1:numPoints1
        W35(:,:,i3,:,i5)=y1(i3)*y1(i5);
    end
end

for i5=1:numPoints1
    for i4=1:numPoints1
        W45(:,:,:,i4,i5)=y1(i4)*y1(i5);
    end
end

W123=zeros(numPoints1,numPoints1,numPoints1,1,1);
W124=zeros(numPoints1,numPoints1,1,numPoints1,1);
W125=zeros(numPoints1,numPoints1,1,1,numPoints1);
W134=zeros(numPoints1,1,numPoints1,numPoints1,1);
W135=zeros(numPoints1,1,numPoints1,1,numPoints1);
W145=zeros(numPoints1,1,1,numPoints1,numPoints1);
W234=zeros(1,numPoints1,numPoints1,numPoints1,1);
W235=zeros(1,numPoints1,numPoints1,1,numPoints1);
W245=zeros(1,numPoints1,1,numPoints1,numPoints1);
W345=zeros(1,1,numPoints1,numPoints1,numPoints1);

for i3=1:numPoints1
    for i2=1:numPoints1
        for i1=1:numPoints1
            W123(i1,i2,i3,:,:)=y1(i1)*y1(i2)*y1(i3);
        end
    end
end

for i4=1:numPoints1
    for i2=1:numPoints1
        for i1=1:numPoints1
            W124(i1,i2,:,i4,:)=y1(i1)*y1(i2)*y1(i4);
        end
    end
end

for i5=1:numPoints1
    for i2=1:numPoints1
        for i1=1:numPoints1
            W125(i1,i2,:,:,i5)=y1(i1)*y1(i2)*y1(i5);
        end
    end
end

for i4=1:numPoints1
    for i3=1:numPoints1
        for i1=1:numPoints1
            W134(i1,:,i3,i4,:)=y1(i1)*y1(i3)*y1(i4);
        end
    end
end

for i5=1:numPoints1
    for i3=1:numPoints1
        for i1=1:numPoints1
            W135(i1,:,i3,:,i5)=y1(i1)*y1(i3)*y1(i5);
        end
    end
end

for i5=1:numPoints1
    for i4=1:numPoints1
        for i1=1:numPoints1
            W145(i1,:,:,i4,i5)=y1(i1)*y1(i4)*y1(i5);
        end
    end
end

for i4=1:numPoints1
    for i3=1:numPoints1
        for i2=1:numPoints1
            W234(:,i2,i3,i4,:)=y1(i2)*y1(i3)*y1(i4);
        end
    end
end

for i5=1:numPoints1
    for i3=1:numPoints1
        for i2=1:numPoints1
            W235(:,i2,i3,:,i5)=y1(i2)*y1(i3)*y1(i5);
        end
    end
end

for i5=1:numPoints1
    for i4=1:numPoints1
        for i2=1:numPoints1
            W245(:,i2,:,i4,i5)=y1(i2)*y1(i4)*y1(i5);
        end
    end
end

for i5=1:numPoints1
    for i4=1:numPoints1
        for i3=1:numPoints1
            W345(:,:,i3,i4,i5)=y1(i3)*y1(i4)*y1(i5);
        end
    end
end

W1234=zeros(numPoints1,numPoints1,numPoints1,numPoints1,1);
W1235=zeros(numPoints1,numPoints1,numPoints1,1,numPoints1);
W1245=zeros(numPoints1,numPoints1,1,numPoints1,numPoints1);
W1345=zeros(numPoints1,1,numPoints1,numPoints1,numPoints1);
W2345=zeros(1,numPoints1,numPoints1,numPoints1,numPoints1);

for i4=1:numPoints1
    for i3=1:numPoints1
        for i2=1:numPoints1
            for i1=1:numPoints1
                W1234(i1,i2,i3,i4,:)=y1(i1)*y1(i2)*y1(i3)*y1(i4);
            end
        end
    end
end

for i5=1:numPoints1
    for i3=1:numPoints1
        for i2=1:numPoints1
            for i1=1:numPoints1
                W1235(i1,i2,i3,:,i5)=y1(i1)*y1(i2)*y1(i3)*y1(i5);
            end
        end
    end
end

for i5=1:numPoints1
    for i4=1:numPoints1
        for i2=1:numPoints1
            for i1=1:numPoints1
                W1245(i1,i2,:,i4,i5)=y1(i1)*y1(i2)*y1(i4)*y1(i5);
            end
        end
    end
end

for i5=1:numPoints1
    for i4=1:numPoints1
        for i3=1:numPoints1
            for i1=1:numPoints1
                W1345(i1,:,i3,i4,i5)=y1(i1)*y1(i3)*y1(i4)*y1(i5);
            end
        end
    end
end

for i5=1:numPoints1
    for i4=1:numPoints1
        for i3=1:numPoints1
            for i2=1:numPoints1
                W2345(:,i2,i3,i4,i5)=y1(i2)*y1(i3)*y1(i4)*y1(i5);
            end
        end
    end
end

W12345=zeros(numPoints1,numPoints1,numPoints1,numPoints1,numPoints1);

for i5=1:numPoints1
    for i4=1:numPoints1
        for i3=1:numPoints1
            for i2=1:numPoints1
                for i1=1:numPoints1
                    W12345(i1,i2,i3,i4,i5)=y1(i1)*y1(i2)*y1(i3)*y1(i4)*y1(i5);
                end
            end
        end
    end
end

G=getStructuralResponse(X1,X2,X3,X4,X5); % structural response G
% ----------------------------------------------------------------------------------------------------------------------



% ----------------------------------------------------------------------------------------------------------------------
% DECOMPOSITION OF STRUCTURAL RESPONSE:
% ----------------------------------------------------------------------------------------------------------------------
G0=zeros(numPoints1,numPoints1,numPoints1,numPoints1,numPoints1);

G0(:,:,:,:,:)=sum(sum(sum(sum(sum(G.*W12345))))); % total mean

G1=zeros(numPoints1,numPoints1,numPoints1,numPoints1,numPoints1);
G2=zeros(numPoints1,numPoints1,numPoints1,numPoints1,numPoints1);
G3=zeros(numPoints1,numPoints1,numPoints1,numPoints1,numPoints1);
G4=zeros(numPoints1,numPoints1,numPoints1,numPoints1,numPoints1);
G5=zeros(numPoints1,numPoints1,numPoints1,numPoints1,numPoints1);

for i1=1:numPoints1
    G1(i1,:,:,:,:)=sum(sum(sum(sum(G(i1,:,:,:,:).*W2345))))-G0(i1,:,:,:,:);
end

for i2=1:numPoints1
    G2(:,i2,:,:,:)=sum(sum(sum(sum(G(:,i2,:,:,:).*W1345))))-G0(:,i2,:,:,:);
end

for i3=1:numPoints1
    G3(:,:,i3,:,:)=sum(sum(sum(sum(G(:,:,i3,:,:).*W1245))))-G0(:,:,i3,:,:);
end

for i4=1:numPoints1
    G4(:,:,:,i4,:)=sum(sum(sum(sum(G(:,:,:,i4,:).*W1235))))-G0(:,:,:,i4,:);
end

for i5=1:numPoints1
    G5(:,:,:,:,i5)=sum(sum(sum(sum(G(:,:,:,:,i5).*W1234))))-G0(:,:,:,:,i5);
end

G12=zeros(numPoints1,numPoints1,numPoints1,numPoints1,numPoints1);
G13=zeros(numPoints1,numPoints1,numPoints1,numPoints1,numPoints1);
G14=zeros(numPoints1,numPoints1,numPoints1,numPoints1,numPoints1);
G15=zeros(numPoints1,numPoints1,numPoints1,numPoints1,numPoints1);
G23=zeros(numPoints1,numPoints1,numPoints1,numPoints1,numPoints1);
G24=zeros(numPoints1,numPoints1,numPoints1,numPoints1,numPoints1);
G25=zeros(numPoints1,numPoints1,numPoints1,numPoints1,numPoints1);
G34=zeros(numPoints1,numPoints1,numPoints1,numPoints1,numPoints1);
G35=zeros(numPoints1,numPoints1,numPoints1,numPoints1,numPoints1);
G45=zeros(numPoints1,numPoints1,numPoints1,numPoints1,numPoints1);

for i2=1:numPoints1
    for i1=1:numPoints1
        G12(i1,i2,:,:,:)=sum(sum(sum(G(i1,i2,:,:,:).*W345)))-(G0(i1,i2,:,:,:)+G1(i1,i2,:,:,:)+G2(i1,i2,:,:,:));
    end
end

for i3=1:numPoints1
    for i1=1:numPoints1
        G13(i1,:,i3,:,:)=sum(sum(sum(G(i1,:,i3,:,:).*W245)))-(G0(i1,:,i3,:,:)+G1(i1,:,i3,:,:)+G3(i1,:,i3,:,:));
    end
end

for i4=1:numPoints1
    for i1=1:numPoints1
        G14(i1,:,:,i4,:)=sum(sum(sum(G(i1,:,:,i4,:).*W235)))-(G0(i1,:,:,i4,:)+G1(i1,:,:,i4,:)+G4(i1,:,:,i4,:));
    end
end

for i5=1:numPoints1
    for i1=1:numPoints1
        G15(i1,:,:,:,i5)=sum(sum(sum(G(i1,:,:,:,i5).*W234)))-(G0(i1,:,:,:,i5)+G1(i1,:,:,:,i5)+G5(i1,:,:,:,i5));
    end
end

for i3=1:numPoints1
    for i2=1:numPoints1
        G23(:,i2,i3,:,:)=sum(sum(sum(G(:,i2,i3,:,:).*W145)))-(G0(:,i2,i3,:,:)+G2(:,i2,i3,:,:)+G3(:,i2,i3,:,:));
    end
end

for i4=1:numPoints1
    for i2=1:numPoints1
        G24(:,i2,:,i4,:)=sum(sum(sum(G(:,i2,:,i4,:).*W135)))-(G0(:,i2,:,i4,:)+G2(:,i2,:,i4,:)+G4(:,i2,:,i4,:));
    end
end

for i5=1:numPoints1
    for i2=1:numPoints1
        G25(:,i2,:,:,i5)=sum(sum(sum(G(:,i2,:,:,i5).*W134)))-(G0(:,i2,:,:,i5)+G2(:,i2,:,:,i5)+G5(:,i2,:,:,i5));
    end
end

for i4=1:numPoints1
    for i3=1:numPoints1
        G34(:,:,i3,i4,:)=sum(sum(sum(G(:,:,i3,i4,:).*W125)))-(G0(:,:,i3,i4,:)+G3(:,:,i3,i4,:)+G4(:,:,i3,i4,:));
    end
end

for i5=1:numPoints1
    for i3=1:numPoints1
        G35(:,:,i3,:,i5)=sum(sum(sum(G(:,:,i3,:,i5).*W124)))-(G0(:,:,i3,:,i5)+G3(:,:,i3,:,i5)+G5(:,:,i3,:,i5));
    end
end

for i5=1:numPoints1
    for i4=1:numPoints1
        G45(:,:,:,i4,i5)=sum(sum(sum(G(:,:,:,i4,i5).*W123)))-(G0(:,:,:,i4,i5)+G4(:,:,:,i4,i5)+G5(:,:,:,i4,i5));
    end
end

G123=zeros(numPoints1,numPoints1,numPoints1,numPoints1,numPoints1);
G124=zeros(numPoints1,numPoints1,numPoints1,numPoints1,numPoints1);
G125=zeros(numPoints1,numPoints1,numPoints1,numPoints1,numPoints1);
G134=zeros(numPoints1,numPoints1,numPoints1,numPoints1,numPoints1);
G135=zeros(numPoints1,numPoints1,numPoints1,numPoints1,numPoints1);
G145=zeros(numPoints1,numPoints1,numPoints1,numPoints1,numPoints1);
G234=zeros(numPoints1,numPoints1,numPoints1,numPoints1,numPoints1);
G235=zeros(numPoints1,numPoints1,numPoints1,numPoints1,numPoints1);
G245=zeros(numPoints1,numPoints1,numPoints1,numPoints1,numPoints1);
G345=zeros(numPoints1,numPoints1,numPoints1,numPoints1,numPoints1);

for i3=1:numPoints1
    for i2=1:numPoints1
        for i1=1:numPoints1
            G123(i1,i2,i3,:,:)=sum(sum(G(i1,i2,i3,:,:).*W45))-(G0(i1,i2,i3,:,:)+G1(i1,i2,i3,:,:)+G2(i1,i2,i3,:,:)...
                +G3(i1,i2,i3,:,:)+G12(i1,i2,i3,:,:)+G13(i1,i2,i3,:,:)+G23(i1,i2,i3,:,:));
        end
    end
end

for i4=1:numPoints1
    for i2=1:numPoints1
        for i1=1:numPoints1
            G124(i1,i2,:,i4,:)=sum(sum(G(i1,i2,:,i4,:).*W35))-(G0(i1,i2,:,i4,:)+G1(i1,i2,:,i4,:)+G2(i1,i2,:,i4,:)...
                +G4(i1,i2,:,i4,:)+G12(i1,i2,:,i4,:)+G14(i1,i2,:,i4,:)+G24(i1,i2,:,i4,:));
        end
    end
end

for i5=1:numPoints1
    for i2=1:numPoints1
        for i1=1:numPoints1
            G125(i1,i2,:,:,i5)=sum(sum(G(i1,i2,:,:,i5).*W34))-(G0(i1,i2,:,:,i5)+G1(i1,i2,:,:,i5)+G2(i1,i2,:,:,i5)...
                +G5(i1,i2,:,:,i5)+G12(i1,i2,:,:,i5)+G15(i1,i2,:,:,i5)+G25(i1,i2,:,:,i5));
        end
    end
end

for i4=1:numPoints1
    for i3=1:numPoints1
        for i1=1:numPoints1
            G134(i1,:,i3,i4,:)=sum(sum(G(i1,:,i3,i4,:).*W25))-(G0(i1,:,i3,i4,:)+G1(i1,:,i3,i4,:)+G3(i1,:,i3,i4,:)...
                +G4(i1,:,i3,i4,:)+G13(i1,:,i3,i4,:)+G14(i1,:,i3,i4,:)+G34(i1,:,i3,i4,:));
        end
    end
end

for i5=1:numPoints1
    for i3=1:numPoints1
        for i1=1:numPoints1
            G135(i1,:,i3,:,i5)=sum(sum(G(i1,:,i3,:,i5).*W24))-(G0(i1,:,i3,:,i5)+G1(i1,:,i3,:,i5)+G3(i1,:,i3,:,i5)...
                +G5(i1,:,i3,:,i5)+G13(i1,:,i3,:,i5)+G15(i1,:,i3,:,i5)+G35(i1,:,i3,:,i5));
        end
    end
end

for i5=1:numPoints1
    for i4=1:numPoints1
        for i1=1:numPoints1
            G145(i1,:,:,i4,i5)=sum(sum(G(i1,:,:,i4,i5).*W23))-(G0(i1,:,:,i4,i5)+G1(i1,:,:,i4,i5)+G4(i1,:,:,i4,i5)...
                +G5(i1,:,:,i4,i5)+G14(i1,:,:,i4,i5)+G15(i1,:,:,i4,i5)+G45(i1,:,:,i4,i5));
        end
    end
end

for i4=1:numPoints1
    for i3=1:numPoints1
        for i2=1:numPoints1
            G234(:,i2,i3,i4,:)=sum(sum(G(:,i2,i3,i4,:).*W15))-(G0(:,i2,i3,i4,:)+G2(:,i2,i3,i4,:)+G3(:,i2,i3,i4,:)...
                +G4(:,i2,i3,i4,:)+G23(:,i2,i3,i4,:)+G24(:,i2,i3,i4,:)+G34(:,i2,i3,i4,:));
        end
    end
end

for i5=1:numPoints1
    for i3=1:numPoints1
        for i2=1:numPoints1
            G235(:,i2,i3,:,i5)=sum(sum(G(:,i2,i3,:,i5).*W14))-(G0(:,i2,i3,:,i5)+G2(:,i2,i3,:,i5)+G3(:,i2,i3,:,i5)...
                +G5(:,i2,i3,:,i5)+G23(:,i2,i3,:,i5)+G25(:,i2,i3,:,i5)+G35(:,i2,i3,:,i5));
        end
    end
end

for i5=1:numPoints1
    for i4=1:numPoints1
        for i2=1:numPoints1
            G245(:,i2,:,i4,i5)=sum(sum(G(:,i2,:,i4,i5).*W13))-(G0(:,i2,:,i4,i5)+G2(:,i2,:,i4,i5)+G4(:,i2,:,i4,i5)...
                +G5(:,i2,:,i4,i5)+G24(:,i2,:,i4,i5)+G25(:,i2,:,i4,i5)+G45(:,i2,:,i4,i5));
        end
    end
end

for i5=1:numPoints1
    for i4=1:numPoints1
        for i3=1:numPoints1
            G345(:,:,i3,i4,i5)=sum(sum(G(:,:,i3,i4,i5).*W12))-(G0(:,:,i3,i4,i5)+G3(:,:,i3,i4,i5)+G4(:,:,i3,i4,i5)...
                +G5(:,:,i3,i4,i5)+G34(:,:,i3,i4,i5)+G35(:,:,i3,i4,i5)+G45(:,:,i3,i4,i5));
        end
    end
end

G1234=zeros(numPoints1,numPoints1,numPoints1,numPoints1,numPoints1);
G1235=zeros(numPoints1,numPoints1,numPoints1,numPoints1,numPoints1);
G1245=zeros(numPoints1,numPoints1,numPoints1,numPoints1,numPoints1);
G1345=zeros(numPoints1,numPoints1,numPoints1,numPoints1,numPoints1);
G2345=zeros(numPoints1,numPoints1,numPoints1,numPoints1,numPoints1);

for i4=1:numPoints1
    for i3=1:numPoints1
        for i2=1:numPoints1
            for i1=1:numPoints1
                G1234(i1,i2,i3,i4,:)=sum(G(i1,i2,i3,i4,:).*W5)-(G0(i1,i2,i3,i4,:)+G1(i1,i2,i3,i4,:)+G2(i1,i2,i3,i4,:)...
                    +G3(i1,i2,i3,i4,:)+G4(i1,i2,i3,i4,:)+G12(i1,i2,i3,i4,:)+G13(i1,i2,i3,i4,:)+G14(i1,i2,i3,i4,:)...
                    +G23(i1,i2,i3,i4,:)+G24(i1,i2,i3,i4,:)+G34(i1,i2,i3,i4,:)+G123(i1,i2,i3,i4,:)+G124(i1,i2,i3,i4,:)...
                    +G134(i1,i2,i3,i4,:)+G234(i1,i2,i3,i4,:));
            end
        end
    end
end

for i5=1:numPoints1
    for i3=1:numPoints1
        for i2=1:numPoints1
            for i1=1:numPoints1
                G1235(i1,i2,i3,:,i5)=sum(G(i1,i2,i3,:,i5).*W4)-(G0(i1,i2,i3,:,i5)+G1(i1,i2,i3,:,i5)+G2(i1,i2,i3,:,i5)...
                    +G3(i1,i2,i3,:,i5)+G5(i1,i2,i3,:,i5)+G12(i1,i2,i3,:,i5)+G13(i1,i2,i3,:,i5)+G15(i1,i2,i3,:,i5)...
                    +G23(i1,i2,i3,:,i5)+G25(i1,i2,i3,:,i5)+G35(i1,i2,i3,:,i5)+G123(i1,i2,i3,:,i5)+G125(i1,i2,i3,:,i5)...
                    +G135(i1,i2,i3,:,i5)+G235(i1,i2,i3,:,i5));
            end
        end
    end
end

for i5=1:numPoints1
    for i4=1:numPoints1
        for i2=1:numPoints1
            for i1=1:numPoints1
                G1245(i1,i2,:,i4,i5)=sum(G(i1,i2,:,i4,i5).*W3)-(G0(i1,i2,:,i4,i5)+G1(i1,i2,:,i4,i5)+G2(i1,i2,:,i4,i5)...
                    +G4(i1,i2,:,i4,i5)+G5(i1,i2,:,i4,i5)+G12(i1,i2,:,i4,i5)+G14(i1,i2,:,i4,i5)+G15(i1,i2,:,i4,i5)...
                    +G24(i1,i2,:,i4,i5)+G25(i1,i2,:,i4,i5)+G45(i1,i2,:,i4,i5)+G124(i1,i2,:,i4,i5)+G125(i1,i2,:,i4,i5)...
                    +G145(i1,i2,:,i4,i5)+G245(i1,i2,:,i4,i5));
            end
        end
    end
end

for i5=1:numPoints1
    for i4=1:numPoints1
        for i3=1:numPoints1
            for i1=1:numPoints1
                G1345(i1,:,i3,i4,i5)=sum(G(i1,:,i3,i4,i5).*W2)-(G0(i1,:,i3,i4,i5)+G1(i1,:,i3,i4,i5)+G3(i1,:,i3,i4,i5)...
                    +G4(i1,:,i3,i4,i5)+G5(i1,:,i3,i4,i5)+G13(i1,:,i3,i4,i5)+G14(i1,:,i3,i4,i5)+G15(i1,:,i3,i4,i5)...
                    +G34(i1,:,i3,i4,i5)+G35(i1,:,i3,i4,i5)+G45(i1,:,i3,i4,i5)+G134(i1,:,i3,i4,i5)+G135(i1,:,i3,i4,i5)...
                    +G145(i1,:,i3,i4,i5)+G345(i1,:,i3,i4,i5));
            end
        end
    end
end

for i5=1:numPoints1
    for i4=1:numPoints1
        for i3=1:numPoints1
            for i2=1:numPoints1
                G2345(:,i2,i3,i4,i5)=sum(G(:,i2,i3,i4,i5).*W1)-(G0(:,i2,i3,i4,i5)+G2(:,i2,i3,i4,i5)+G3(:,i2,i3,i4,i5)...
                    +G4(:,i2,i3,i4,i5)+G5(:,i2,i3,i4,i5)+G23(:,i2,i3,i4,i5)+G24(:,i2,i3,i4,i5)+G25(:,i2,i3,i4,i5)...
                    +G34(:,i2,i3,i4,i5)+G35(:,i2,i3,i4,i5)+G45(:,i2,i3,i4,i5)+G234(:,i2,i3,i4,i5)+G235(:,i2,i3,i4,i5)...
                    +G245(:,i2,i3,i4,i5)+G345(:,i2,i3,i4,i5));
            end
        end
    end
end

G12345=zeros(numPoints1,numPoints1,numPoints1,numPoints1,numPoints1);

G12345(:,:,:,:,:)=G-(G0+G1+G2+G3+G4+G5+G12+G13+G14+G15+G23+G24+G25+G34+G35+G45+G123+G124+G125+G134+G135+G145+G234...
    +G235+G245+G345+G1234+G1235+G1245+G1345+G2345);
% ----------------------------------------------------------------------------------------------------------------------



% ----------------------------------------------------------------------------------------------------------------------
% DECOMPOSITION OF TOTAL VARIANCE:
% ----------------------------------------------------------------------------------------------------------------------
DG=sum(sum(sum(sum(sum(G.^2.*W12345)))))-G0(1,1,1,1,1)^2; % total variance

D1=sum(G1(:,1,1,1,1).^2.*W1);
D2=sum(G2(1,:,1,1,1).^2.*W2);
D3=sum(G3(1,1,:,1,1).^2.*W3);
D4=sum(G4(1,1,1,:,1).^2.*W4);
D5=sum(G5(1,1,1,1,:).^2.*W5);

D12=sum(sum(G12(:,:,1,1,1).^2.*W12));
D13=sum(sum(G13(:,1,:,1,1).^2.*W13));
D14=sum(sum(G14(:,1,1,:,1).^2.*W14));
D15=sum(sum(G15(:,1,1,1,:).^2.*W15));
D23=sum(sum(G23(1,:,:,1,1).^2.*W23));
D24=sum(sum(G24(1,:,1,:,1).^2.*W24));
D25=sum(sum(G25(1,:,1,1,:).^2.*W25));
D34=sum(sum(G34(1,1,:,:,1).^2.*W34));
D35=sum(sum(G35(1,1,:,1,:).^2.*W35));
D45=sum(sum(G45(1,1,1,:,:).^2.*W45));

D123=sum(sum(sum(G123(:,:,:,1,1).^2.*W123)));
D124=sum(sum(sum(G124(:,:,1,:,1).^2.*W124)));
D125=sum(sum(sum(G125(:,:,1,1,:).^2.*W125)));
D134=sum(sum(sum(G134(:,1,:,:,1).^2.*W134)));
D135=sum(sum(sum(G135(:,1,:,1,:).^2.*W135)));
D145=sum(sum(sum(G145(:,1,1,:,:).^2.*W145)));
D234=sum(sum(sum(G234(1,:,:,:,1).^2.*W234)));
D235=sum(sum(sum(G235(1,:,:,1,:).^2.*W235)));
D245=sum(sum(sum(G245(1,:,1,:,:).^2.*W245)));
D345=sum(sum(sum(G345(1,1,:,:,:).^2.*W345)));

D1234=sum(sum(sum(sum(G1234(:,:,:,:,1).^2.*W1234))));
D1235=sum(sum(sum(sum(G1235(:,:,:,1,:).^2.*W1235))));
D1245=sum(sum(sum(sum(G1245(:,:,1,:,:).^2.*W1245))));
D1345=sum(sum(sum(sum(G1345(:,1,:,:,:).^2.*W1345))));
D2345=sum(sum(sum(sum(G2345(1,:,:,:,:).^2.*W2345))));

D12345=sum(sum(sum(sum(sum(G12345.^2.*W12345)))));
% ----------------------------------------------------------------------------------------------------------------------



% ----------------------------------------------------------------------------------------------------------------------
% CONTRIBUTIONS TO TOTAL VARIANCE:
% ----------------------------------------------------------------------------------------------------------------------
S1=D1/DG; % contribution to total variance due to x1
S2=D2/DG; % contribution to total variance due to x2
S3=D3/DG; % contribution to total variance due to x3
S4=D4/DG; % contribution to total variance due to x4
S5=D5/DG; % contribution to total variance due to x5

S12=D12/DG; % contribution to total variance due to the interaction between x1 and x2
S13=D13/DG; % contribution to total variance due to the interaction between x1 and x3
S14=D14/DG; % contribution to total variance due to the interaction between x1 and x4
S15=D15/DG; % contribution to total variance due to the interaction between x1 and x5
S23=D23/DG; % contribution to total variance due to the interaction between x2 and x3
S24=D24/DG; % contribution to total variance due to the interaction between x2 and x4
S25=D25/DG; % contribution to total variance due to the interaction between x2 and x5
S34=D34/DG; % contribution to total variance due to the interaction between x3 and x4
S35=D35/DG; % contribution to total variance due to the interaction between x3 and x5
S45=D45/DG; % contribution to total variance due to the interaction between x4 and x5

S123=D123/DG; % contribution to total variance due to the interaction between x1, x2 and x3
S124=D124/DG; % contribution to total variance due to the interaction between x1, x2 and x4
S125=D125/DG; % contribution to total variance due to the interaction between x1, x2 and x5
S134=D134/DG; % contribution to total variance due to the interaction between x1, x3 and x4
S135=D135/DG; % contribution to total variance due to the interaction between x1, x3 and x5
S145=D145/DG; % contribution to total variance due to the interaction between x1, x4 and x5
S234=D234/DG; % contribution to total variance due to the interaction between x2, x3 and x4
S235=D235/DG; % contribution to total variance due to the interaction between x2, x3 and x5
S245=D245/DG; % contribution to total variance due to the interaction between x2, x4 and x5
S345=D345/DG; % contribution to total variance due to the interaction between x3, x4 and x5

S1234=D1234/DG; % contribution to total variance due to the interaction between x1, x2, x3 and x4
S1235=D1235/DG; % contribution to total variance due to the interaction between x1, x2, x3 and x5
S1245=D1245/DG; % contribution to total variance due to the interaction between x1, x2, x4 and x5
S1345=D1345/DG; % contribution to total variance due to the interaction between x1, x3, x4 and x5
S2345=D2345/DG; % contribution to total variance due to the interaction between x2, x3, x4 and x5

S12345=D12345/DG; % contribution to total variance due to the interaction between x1, x2, x3, x4 and x5
% ----------------------------------------------------------------------------------------------------------------------



% ----------------------------------------------------------------------------------------------------------------------
% SUMMARY:
% ----------------------------------------------------------------------------------------------------------------------
fprintf('Summary:\n\n')
fprintf('Total mean: %.6f.\n',G0(1,1,1,1))
fprintf('Total variance: %.4f.\n',DG)
fprintf('\n')
fprintf('Random variable x1 contributes to the total variance response by: %.2f%%.\n',S1*100)
fprintf('Random variable x2 contributes to the total variance response by: %.2f%%.\n',S2*100)
fprintf('Random variable x3 contributes to the total variance response by: %.2f%%.\n',S3*100)
fprintf('Random variable x4 contributes to the total variance response by: %.2f%%.\n',S4*100)
fprintf('Random variable x5 contributes to the total variance response by: %.2f%%.\n',S5*100)
fprintf('\n')
fprintf('Random variables x1 and x2 contribute to the total variance response by: %.2f%%.\n',S12*100)
fprintf('Random variables x1 and x3 contribute to the total variance response by: %.2f%%.\n',S13*100)
fprintf('Random variables x1 and x4 contribute to the total variance response by: %.2f%%.\n',S14*100)
fprintf('Random variables x1 and x5 contribute to the total variance response by: %.2f%%.\n',S15*100)
fprintf('Random variables x2 and x3 contribute to the total variance response by: %.2f%%.\n',S23*100)
fprintf('Random variables x2 and x4 contribute to the total variance response by: %.2f%%.\n',S24*100)
fprintf('Random variables x2 and x5 contribute to the total variance response by: %.2f%%.\n',S25*100)
fprintf('Random variables x3 and x4 contribute to the total variance response by: %.2f%%.\n',S34*100)
fprintf('Random variables x3 and x5 contribute to the total variance response by: %.2f%%.\n',S35*100)
fprintf('Random variables x4 and x5 contribute to the total variance response by: %.2f%%.\n',S45*100)
fprintf('\n')
fprintf('Random variables x1, x2 and x3 contribute to the total variance response by: %.2f%%.\n',S123*100)
fprintf('Random variables x1, x2 and x4 contribute to the total variance response by: %.2f%%.\n',S124*100)
fprintf('Random variables x1, x2 and x5 contribute to the total variance response by: %.2f%%.\n',S125*100)
fprintf('Random variables x1, x3 and x4 contribute to the total variance response by: %.2f%%.\n',S134*100)
fprintf('Random variables x1, x3 and x5 contribute to the total variance response by: %.2f%%.\n',S135*100)
fprintf('Random variables x1, x4 and x5 contribute to the total variance response by: %.2f%%.\n',S145*100)
fprintf('Random variables x2, x3 and x4 contribute to the total variance response by: %.2f%%.\n',S234*100)
fprintf('Random variables x2, x3 and x5 contribute to the total variance response by: %.2f%%.\n',S235*100)
fprintf('Random variables x2, x4 and x5 contribute to the total variance response by: %.2f%%.\n',S245*100)
fprintf('Random variables x3, x4 and x5 contribute to the total variance response by: %.2f%%.\n',S345*100)
fprintf('\n')
fprintf('Random variables x1, x2, x3 and x4 contribute to the total variance response by: %.2f%%.\n',S1234*100)
fprintf('Random variables x1, x2, x3 and x5 contribute to the total variance response by: %.2f%%.\n',S1235*100)
fprintf('Random variables x1, x2, x4 and x5 contribute to the total variance response by: %.2f%%.\n',S1245*100)
fprintf('Random variables x1, x3, x4 and x5 contribute to the total variance response by: %.2f%%.\n',S1345*100)
fprintf('Random variables x2, x3, x4 and x5 contribute to the total variance response by: %.2f%%.\n',S2345*100)
fprintf('\n')
fprintf('Random variables x1, x2, x3, x4 and x5 contribute to the total variance response by: %.2f%%.\n',S12345*100)
% ----------------------------------------------------------------------------------------------------------------------



% ----------------------------------------------------------------------------------------------------------------------
% FUNCTIONS:
% ----------------------------------------------------------------------------------------------------------------------
function G=getStructuralResponse(X1,X2,X3,X4,X5)
syms x1 x2 x3 x4 x5

p1=legendreP(2,x1)+jacobiP(3,3,5,x1); % a univariate function defined on x1 = [-1,1]
p2=legendreP(5,x2)+jacobiP(2,2,7,x2); % a univariate function defined on x2 = [-1,1]
p3=legendreP(3,x3)+jacobiP(2,4,1,x3); % a univariate function defined on x3 = [-1,1]
p4=legendreP(4,x4)+jacobiP(5,3,2,x4); % a univariate function defined on x4 = [-1,1]
p5=legendreP(3,x5)+jacobiP(4,5,4,x5); % a univariate function defined on x5 = [-1,1]

% Structural response of interest (notice that this g is a function of x1, x2, x3, x4 and x5):
g=p1*p2*p3*p4*p5/1e5; % a pentavariate function defined on x1 x x2 x x3 x x4 x x5 = [-1,1]^5

g=matlabFunction(g,'Vars',[x1,x2,x3,x4,x5]);

G=g(X1,X2,X3,X4,X5);
end
% ----------------------------------------------------------------------------------------------------------------------
