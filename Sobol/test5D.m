% Hugo Esquivel, 2021
% -

clearvars; close all; clc;

% ----------------------------------------------------------------------------------------------------------------------
% Input:
% ----------------------------------------------------------------------------------------------------------------------
syms x1 x2 x3 x4 x5 % random variables x1, x2, x3, x4 and x5

f=getStructuralResponse(x1,x2,x3,x4,x5); % structural response f

probabilityDistribution='Beta'; % options: Uniform or Beta

switch probabilityDistribution
    case 'Uniform'
        pdf1=1/2; % probability density of x1 (assuming that x1 is a uniformly distributed random variable in [-1,1])
        pdf2=1/2; % probability density of x2 (assuming that x2 is a uniformly distributed random variable in [-1,1])
        pdf3=1/2; % probability density of x3 (assuming that x3 is a uniformly distributed random variable in [-1,1])
        pdf4=1/2; % probability density of x4 (assuming that x4 is a uniformly distributed random variable in [-1,1])
        pdf5=1/2; % probability density of x5 (assuming that x5 is a uniformly distributed random variable in [-1,1])
        
    case 'Beta'
        % probability density of x1 (assuming that x1 is a beta distributed random variable in [-1,1])
        alpha1=2;
        beta1=5;
        pdf1=(x1+1)^(alpha1-1)*(1-x1)^(beta1-1)/(2^(alpha1+beta1-1)*beta(alpha1,beta1));
        
        % probability density of x2 (assuming that x2 is a beta distributed random variable in [-1,1])
        alpha2=5;
        beta2=1;
        pdf2=(x2+1)^(alpha2-1)*(1-x2)^(beta2-1)/(2^(alpha2+beta2-1)*beta(alpha2,beta2));
        
        % probability density of x3 (assuming that x3 is a beta distributed random variable in [-1,1])
        alpha3=2;
        beta3=2;
        pdf3=(x3+1)^(alpha3-1)*(1-x3)^(beta3-1)/(2^(alpha3+beta3-1)*beta(alpha3,beta3));
        
        % probability density of x4 (assuming that x4 is a beta distributed random variable in [-1,1])
        alpha4=2;
        beta4=4;
        pdf4=(x4+1)^(alpha4-1)*(1-x4)^(beta4-1)/(2^(alpha4+beta4-1)*beta(alpha4,beta4));
        
        % probability density of x5 (assuming that x5 is a beta distributed random variable in [-1,1])
        alpha5=1;
        beta5=3;
        pdf5=(x5+1)^(alpha5-1)*(1-x5)^(beta5-1)/(2^(alpha5+beta5-1)*beta(alpha5,beta5));
end

pdf12=pdf1*pdf2; % joint probability density of x1 and x2
pdf13=pdf1*pdf3; % joint probability density of x1 and x3
pdf14=pdf1*pdf4; % joint probability density of x1 and x4
pdf15=pdf1*pdf5; % joint probability density of x1 and x5
pdf23=pdf2*pdf3; % joint probability density of x2 and x3
pdf24=pdf2*pdf4; % joint probability density of x2 and x4
pdf25=pdf2*pdf5; % joint probability density of x2 and x5
pdf34=pdf3*pdf4; % joint probability density of x3 and x4
pdf35=pdf3*pdf5; % joint probability density of x3 and x5
pdf45=pdf4*pdf5; % joint probability density of x4 and x5

pdf123=pdf1*pdf2*pdf3; % joint probability density of x1, x2 and x3
pdf124=pdf1*pdf2*pdf4; % joint probability density of x1, x2 and x4
pdf125=pdf1*pdf2*pdf5; % joint probability density of x1, x2 and x5
pdf134=pdf1*pdf3*pdf4; % joint probability density of x1, x3 and x4
pdf135=pdf1*pdf3*pdf5; % joint probability density of x1, x3 and x5
pdf145=pdf1*pdf4*pdf5; % joint probability density of x1, x4 and x5
pdf234=pdf2*pdf3*pdf4; % joint probability density of x2, x3 and x4
pdf235=pdf2*pdf3*pdf5; % joint probability density of x2, x3 and x5
pdf245=pdf2*pdf4*pdf5; % joint probability density of x2, x4 and x5
pdf345=pdf3*pdf4*pdf5; % joint probability density of x3, x4 and x5

pdf1234=pdf1*pdf2*pdf3*pdf4; % joint probability density of x1, x2, x3 and x4
pdf1235=pdf1*pdf2*pdf3*pdf5; % joint probability density of x1, x2, x3 and x5
pdf1245=pdf1*pdf2*pdf4*pdf5; % joint probability density of x1, x2, x4 and x5
pdf1345=pdf1*pdf3*pdf4*pdf5; % joint probability density of x1, x3, x4 and x5
pdf2345=pdf2*pdf3*pdf4*pdf5; % joint probability density of x2, x3, x4 and x5

pdf12345=pdf1*pdf2*pdf3*pdf4*pdf5; % joint probability density of x1, x2, x3, x4 and x5
% ----------------------------------------------------------------------------------------------------------------------


% ----------------------------------------------------------------------------------------------------------------------
% Decomposition of structural response:
% ----------------------------------------------------------------------------------------------------------------------
f0=int(int(int(int(int(f*pdf12345,x1,-1,1),x2,-1,1),x3,-1,1),x4,-1,1),x5,-1,1); % total mean

f1=int(int(int(int(f*pdf2345,x2,-1,1),x3,-1,1),x4,-1,1),x5,-1,1)-f0;
f2=int(int(int(int(f*pdf1345,x1,-1,1),x3,-1,1),x4,-1,1),x5,-1,1)-f0;
f3=int(int(int(int(f*pdf1245,x1,-1,1),x2,-1,1),x4,-1,1),x5,-1,1)-f0;
f4=int(int(int(int(f*pdf1235,x1,-1,1),x2,-1,1),x3,-1,1),x5,-1,1)-f0;
f5=int(int(int(int(f*pdf1234,x1,-1,1),x2,-1,1),x3,-1,1),x4,-1,1)-f0;

f12=int(int(int(f*pdf345,x3,-1,1),x4,-1,1),x5,-1,1)-(f0+f1+f2);
f13=int(int(int(f*pdf245,x2,-1,1),x4,-1,1),x5,-1,1)-(f0+f1+f3);
f14=int(int(int(f*pdf235,x2,-1,1),x3,-1,1),x5,-1,1)-(f0+f1+f4);
f15=int(int(int(f*pdf234,x2,-1,1),x3,-1,1),x4,-1,1)-(f0+f1+f5);
f23=int(int(int(f*pdf145,x1,-1,1),x4,-1,1),x5,-1,1)-(f0+f2+f3);
f24=int(int(int(f*pdf135,x1,-1,1),x3,-1,1),x5,-1,1)-(f0+f2+f4);
f25=int(int(int(f*pdf134,x1,-1,1),x3,-1,1),x4,-1,1)-(f0+f2+f5);
f34=int(int(int(f*pdf125,x1,-1,1),x2,-1,1),x5,-1,1)-(f0+f3+f4);
f35=int(int(int(f*pdf124,x1,-1,1),x2,-1,1),x4,-1,1)-(f0+f3+f5);
f45=int(int(int(f*pdf123,x1,-1,1),x2,-1,1),x3,-1,1)-(f0+f4+f5);

f123=int(int(f*pdf45,x4,-1,1),x5,-1,1)-(f0+f1+f2+f3+f12+f13+f23);
f124=int(int(f*pdf35,x3,-1,1),x5,-1,1)-(f0+f1+f2+f4+f12+f14+f24);
f125=int(int(f*pdf34,x3,-1,1),x4,-1,1)-(f0+f1+f2+f5+f12+f15+f25);
f134=int(int(f*pdf25,x2,-1,1),x5,-1,1)-(f0+f1+f3+f4+f13+f14+f34);
f135=int(int(f*pdf24,x2,-1,1),x4,-1,1)-(f0+f1+f3+f5+f13+f15+f35);
f145=int(int(f*pdf23,x2,-1,1),x3,-1,1)-(f0+f1+f4+f5+f14+f15+f45);
f234=int(int(f*pdf15,x1,-1,1),x5,-1,1)-(f0+f2+f3+f4+f23+f24+f34);
f235=int(int(f*pdf14,x1,-1,1),x4,-1,1)-(f0+f2+f3+f5+f23+f25+f35);
f245=int(int(f*pdf13,x1,-1,1),x3,-1,1)-(f0+f2+f4+f5+f24+f25+f45);
f345=int(int(f*pdf12,x1,-1,1),x2,-1,1)-(f0+f3+f4+f5+f34+f35+f45);

f1234=int(f*pdf5,x5,-1,1)-(f0+f1+f2+f3+f4+f12+f13+f14+f23+f24+f34+f123+f124+f134+f234);
f1235=int(f*pdf4,x4,-1,1)-(f0+f1+f2+f3+f5+f12+f13+f15+f23+f25+f35+f123+f125+f135+f235);
f1245=int(f*pdf3,x3,-1,1)-(f0+f1+f2+f4+f5+f12+f14+f15+f24+f25+f45+f124+f125+f145+f245);
f1345=int(f*pdf2,x2,-1,1)-(f0+f1+f3+f4+f5+f13+f14+f15+f34+f35+f45+f134+f135+f145+f345);
f2345=int(f*pdf1,x1,-1,1)-(f0+f2+f3+f4+f5+f23+f24+f25+f34+f35+f45+f234+f235+f245+f345);

f12345=f-(f0+f1+f2+f3+f4+f5+f12+f13+f14+f15+f23+f24+f25+f34+f35+f45+f123+f124+f125+f134+f135+f145+f234+f235+f245+f345...
    +f1234+f1235+f1245+f1345+f2345);
% ----------------------------------------------------------------------------------------------------------------------


% ----------------------------------------------------------------------------------------------------------------------
% Checking orthogonality of decomposition:
% ----------------------------------------------------------------------------------------------------------------------
check1=int(f1*pdf1,x1,-1,1)==0;
check2=int(f2*pdf2,x2,-1,1)==0;
check3=int(f3*pdf3,x3,-1,1)==0;
check4=int(f4*pdf4,x4,-1,1)==0;
check5=int(f5*pdf5,x5,-1,1)==0;

check12=int(int(f12*pdf12,x1,-1,1),x2,-1,1)==0;
check13=int(int(f13*pdf13,x1,-1,1),x3,-1,1)==0;
check14=int(int(f14*pdf14,x1,-1,1),x4,-1,1)==0;
check15=int(int(f15*pdf15,x1,-1,1),x5,-1,1)==0;
check23=int(int(f23*pdf23,x2,-1,1),x3,-1,1)==0;
check24=int(int(f24*pdf24,x2,-1,1),x4,-1,1)==0;
check25=int(int(f25*pdf25,x2,-1,1),x5,-1,1)==0;
check34=int(int(f34*pdf34,x3,-1,1),x4,-1,1)==0;
check35=int(int(f35*pdf35,x3,-1,1),x5,-1,1)==0;
check45=int(int(f45*pdf45,x4,-1,1),x5,-1,1)==0;

check123=int(int(int(f123*pdf123,x1,-1,1),x2,-1,1),x3,-1,1)==0;
check124=int(int(int(f124*pdf124,x1,-1,1),x2,-1,1),x4,-1,1)==0;
check125=int(int(int(f125*pdf125,x1,-1,1),x2,-1,1),x5,-1,1)==0;
check134=int(int(int(f134*pdf134,x1,-1,1),x3,-1,1),x4,-1,1)==0;
check135=int(int(int(f135*pdf135,x1,-1,1),x3,-1,1),x5,-1,1)==0;
check145=int(int(int(f145*pdf145,x1,-1,1),x4,-1,1),x5,-1,1)==0;
check234=int(int(int(f234*pdf234,x2,-1,1),x3,-1,1),x4,-1,1)==0;
check235=int(int(int(f235*pdf235,x2,-1,1),x3,-1,1),x5,-1,1)==0;
check245=int(int(int(f245*pdf245,x2,-1,1),x4,-1,1),x5,-1,1)==0;
check345=int(int(int(f345*pdf345,x3,-1,1),x4,-1,1),x5,-1,1)==0;

check1234=int(int(int(int(f1234*pdf1234,x1,-1,1),x2,-1,1),x3,-1,1),x4,-1,1)==0;
check1235=int(int(int(int(f1235*pdf1235,x1,-1,1),x2,-1,1),x3,-1,1),x5,-1,1)==0;
check1245=int(int(int(int(f1245*pdf1245,x1,-1,1),x2,-1,1),x4,-1,1),x5,-1,1)==0;
check1345=int(int(int(int(f1345*pdf1345,x1,-1,1),x3,-1,1),x4,-1,1),x5,-1,1)==0;
check2345=int(int(int(int(f2345*pdf2345,x2,-1,1),x3,-1,1),x4,-1,1),x5,-1,1)==0;

check12345=int(int(int(int(int(f12345*pdf12345,x1,-1,1),x2,-1,1),x3,-1,1),x4,-1,1),x5,-1,1)==0;

if ~(check1 && check2 && check3 && check4 && check5...
        && check12 && check13 && check14 && check15 && check23 && check24 && check25 && check34 && check35 && check45...
        && check123 && check124 && check125 && check134 && check135 && check145 && check234 && check235 && check245...
        && check345 && check1234 && check1235 && check1245 && check1345 && check2345 && check12345)
    
    error('Orthogonality check did not pass.')
end
% ----------------------------------------------------------------------------------------------------------------------


% ----------------------------------------------------------------------------------------------------------------------
% Decomposition of total variance:
% ----------------------------------------------------------------------------------------------------------------------
Df=int(int(int(int(int(f^2*pdf12345,x1,-1,1),x2,-1,1),x3,-1,1),x4,-1,1),x5,-1,1)-f0^2; % total variance

D1=int(f1^2*pdf1,x1,-1,1);
D2=int(f2^2*pdf2,x2,-1,1);
D3=int(f3^2*pdf3,x3,-1,1);
D4=int(f4^2*pdf4,x4,-1,1);
D5=int(f5^2*pdf5,x5,-1,1);

D12=int(int(f12^2*pdf12,x1,-1,1),x2,-1,1);
D13=int(int(f13^2*pdf13,x1,-1,1),x3,-1,1);
D14=int(int(f14^2*pdf14,x1,-1,1),x4,-1,1);
D15=int(int(f15^2*pdf15,x1,-1,1),x5,-1,1);
D23=int(int(f23^2*pdf23,x2,-1,1),x3,-1,1);
D24=int(int(f24^2*pdf24,x2,-1,1),x4,-1,1);
D25=int(int(f25^2*pdf25,x2,-1,1),x5,-1,1);
D34=int(int(f34^2*pdf34,x3,-1,1),x4,-1,1);
D35=int(int(f35^2*pdf35,x3,-1,1),x5,-1,1);
D45=int(int(f45^2*pdf45,x4,-1,1),x5,-1,1);

D123=int(int(int(f123^2*pdf123,x1,-1,1),x2,-1,1),x3,-1,1);
D124=int(int(int(f124^2*pdf124,x1,-1,1),x2,-1,1),x4,-1,1);
D125=int(int(int(f125^2*pdf125,x1,-1,1),x2,-1,1),x5,-1,1);
D134=int(int(int(f134^2*pdf134,x1,-1,1),x3,-1,1),x4,-1,1);
D135=int(int(int(f135^2*pdf135,x1,-1,1),x3,-1,1),x5,-1,1);
D145=int(int(int(f145^2*pdf145,x1,-1,1),x4,-1,1),x5,-1,1);
D234=int(int(int(f234^2*pdf234,x2,-1,1),x3,-1,1),x4,-1,1);
D235=int(int(int(f235^2*pdf235,x2,-1,1),x3,-1,1),x5,-1,1);
D245=int(int(int(f245^2*pdf245,x2,-1,1),x4,-1,1),x5,-1,1);
D345=int(int(int(f345^2*pdf345,x3,-1,1),x4,-1,1),x5,-1,1);

D1234=int(int(int(int(f1234^2*pdf1234,x1,-1,1),x2,-1,1),x3,-1,1),x4,-1,1);
D1235=int(int(int(int(f1235^2*pdf1235,x1,-1,1),x2,-1,1),x3,-1,1),x5,-1,1);
D1245=int(int(int(int(f1245^2*pdf1245,x1,-1,1),x2,-1,1),x4,-1,1),x5,-1,1);
D1345=int(int(int(int(f1345^2*pdf1345,x1,-1,1),x3,-1,1),x4,-1,1),x5,-1,1);
D2345=int(int(int(int(f2345^2*pdf2345,x2,-1,1),x3,-1,1),x4,-1,1),x5,-1,1);

D12345=int(int(int(int(int(f12345^2*pdf12345,x1,-1,1),x2,-1,1),x3,-1,1),x4,-1,1),x5,-1,1);
% ----------------------------------------------------------------------------------------------------------------------


% ----------------------------------------------------------------------------------------------------------------------
% Contributions to total variance:
% ----------------------------------------------------------------------------------------------------------------------
S1=D1/Df; % contribution to total variance due to x1
S2=D2/Df; % contribution to total variance due to x2
S3=D3/Df; % contribution to total variance due to x3
S4=D4/Df; % contribution to total variance due to x4
S5=D5/Df; % contribution to total variance due to x5

S12=D12/Df; % contribution to total variance due to the interaction between x1 and x2
S13=D13/Df; % contribution to total variance due to the interaction between x1 and x3
S14=D14/Df; % contribution to total variance due to the interaction between x1 and x4
S15=D15/Df; % contribution to total variance due to the interaction between x1 and x5
S23=D23/Df; % contribution to total variance due to the interaction between x2 and x3
S24=D24/Df; % contribution to total variance due to the interaction between x2 and x4
S25=D25/Df; % contribution to total variance due to the interaction between x2 and x5
S34=D34/Df; % contribution to total variance due to the interaction between x3 and x4
S35=D35/Df; % contribution to total variance due to the interaction between x3 and x5
S45=D45/Df; % contribution to total variance due to the interaction between x4 and x5

S123=D123/Df; % contribution to total variance due to the interaction between x1, x2 and x3
S124=D124/Df; % contribution to total variance due to the interaction between x1, x2 and x4
S125=D125/Df; % contribution to total variance due to the interaction between x1, x2 and x5
S134=D134/Df; % contribution to total variance due to the interaction between x1, x3 and x4
S135=D135/Df; % contribution to total variance due to the interaction between x1, x3 and x5
S145=D145/Df; % contribution to total variance due to the interaction between x1, x4 and x5
S234=D234/Df; % contribution to total variance due to the interaction between x2, x3 and x4
S235=D235/Df; % contribution to total variance due to the interaction between x2, x3 and x5
S245=D245/Df; % contribution to total variance due to the interaction between x2, x4 and x5
S345=D345/Df; % contribution to total variance due to the interaction between x3, x4 and x5

S1234=D1234/Df; % contribution to total variance due to the interaction between x1, x2, x3 and x4
S1235=D1235/Df; % contribution to total variance due to the interaction between x1, x2, x3 and x5
S1245=D1245/Df; % contribution to total variance due to the interaction between x1, x2, x4 and x5
S1345=D1345/Df; % contribution to total variance due to the interaction between x1, x3, x4 and x5
S2345=D2345/Df; % contribution to total variance due to the interaction between x2, x3, x4 and x5

S12345=D12345/Df; % contribution to total variance due to the interaction between x1, x2, x3, x4 and x5

disp(double([S1;S2;S3;S4;S5;S12;S13;S14;S15;S23;S24;S25;S34;S35;S45;S123;S124;S125;S134;S135;S145;S234;S235;S245;...
    S345;S1234;S1235;S1245;S1345;S2345;S12345]))
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
% Functions:
% ----------------------------------------------------------------------------------------------------------------------
function f=getStructuralResponse(x1,x2,x3,x4,x5)
p1=legendreP(2,x1)+jacobiP(3,3,5,x1); % a univariate function defined on x1 = [-1,1]
p2=legendreP(5,x2)+jacobiP(2,2,7,x2); % a univariate function defined on x2 = [-1,1]
p3=legendreP(3,x3)+jacobiP(2,4,1,x3); % a univariate function defined on x3 = [-1,1]
p4=legendreP(4,x4)+jacobiP(5,3,2,x4); % a univariate function defined on x4 = [-1,1]
p5=legendreP(3,x5)+jacobiP(4,5,4,x5); % a univariate function defined on x5 = [-1,1]

% Structural response of interest (notice that this f is a function of x1, x2, x3, x4 and x5):
f=p1*p2*p3*p4*p5/1e5; % a pentavariate function defined on x1 x x2 x x3 x x4 x x5 = [-1,1]^5
end
% ----------------------------------------------------------------------------------------------------------------------
