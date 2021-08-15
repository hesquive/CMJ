% Hugo Esquivel, 2021
% -

clearvars; close all; clc;


% ----------------------------------------------------------------------------------------------------------------------
% INPUT:
% ----------------------------------------------------------------------------------------------------------------------
syms x1 x2 x3 x4 x5 % random variables x1, x2, x3, x4 and x5

g=getStructuralResponse(x1,x2,x3,x4,x5); % structural response g

probabilityDistribution='Beta'; % options: Uniform or Beta

switch probabilityDistribution
    case 'Uniform'
        f1=1/2; % probability density of x1 (assuming that x1 is a uniformly distributed random variable in [-1,1])
        f2=1/2; % probability density of x2 (assuming that x2 is a uniformly distributed random variable in [-1,1])
        f3=1/2; % probability density of x3 (assuming that x3 is a uniformly distributed random variable in [-1,1])
        f4=1/2; % probability density of x4 (assuming that x4 is a uniformly distributed random variable in [-1,1])
        f5=1/2; % probability density of x5 (assuming that x5 is a uniformly distributed random variable in [-1,1])
        
    case 'Beta'
        % probability density of x1 (assuming that x1 is a beta distributed random variable in [-1,1])
        alpha1=2;
        beta1=5;
        f1=(x1+1)^(alpha1-1)*(1-x1)^(beta1-1)/(2^(alpha1+beta1-1)*beta(alpha1,beta1));
        
        % probability density of x2 (assuming that x2 is a beta distributed random variable in [-1,1])
        alpha2=2;
        beta2=5;
        f2=(x2+1)^(alpha2-1)*(1-x2)^(beta2-1)/(2^(alpha2+beta2-1)*beta(alpha2,beta2));
        
        % probability density of x3 (assuming that x3 is a beta distributed random variable in [-1,1])
        alpha3=2;
        beta3=5;
        f3=(x3+1)^(alpha3-1)*(1-x3)^(beta3-1)/(2^(alpha3+beta3-1)*beta(alpha3,beta3));
        
        % probability density of x4 (assuming that x4 is a beta distributed random variable in [-1,1])
        alpha4=2;
        beta4=5;
        f4=(x4+1)^(alpha4-1)*(1-x4)^(beta4-1)/(2^(alpha4+beta4-1)*beta(alpha4,beta4));
        
        % probability density of x5 (assuming that x5 is a beta distributed random variable in [-1,1])
        alpha5=2;
        beta5=5;
        f5=(x5+1)^(alpha5-1)*(1-x5)^(beta5-1)/(2^(alpha5+beta5-1)*beta(alpha5,beta5));
end

f12=f1*f2; % joint probability density of x1 and x2
f13=f1*f3; % joint probability density of x1 and x3
f14=f1*f4; % joint probability density of x1 and x4
f15=f1*f5; % joint probability density of x1 and x5
f23=f2*f3; % joint probability density of x2 and x3
f24=f2*f4; % joint probability density of x2 and x4
f25=f2*f5; % joint probability density of x2 and x5
f34=f3*f4; % joint probability density of x3 and x4
f35=f3*f5; % joint probability density of x3 and x5
f45=f4*f5; % joint probability density of x4 and x5

f123=f1*f2*f3; % joint probability density of x1, x2 and x3
f124=f1*f2*f4; % joint probability density of x1, x2 and x4
f125=f1*f2*f5; % joint probability density of x1, x2 and x5
f134=f1*f3*f4; % joint probability density of x1, x3 and x4
f135=f1*f3*f5; % joint probability density of x1, x3 and x5
f145=f1*f4*f5; % joint probability density of x1, x4 and x5
f234=f2*f3*f4; % joint probability density of x2, x3 and x4
f235=f2*f3*f5; % joint probability density of x2, x3 and x5
f245=f2*f4*f5; % joint probability density of x2, x4 and x5
f345=f3*f4*f5; % joint probability density of x3, x4 and x5

f1234=f1*f2*f3*f4; % joint probability density of x1, x2, x3 and x4
f1235=f1*f2*f3*f5; % joint probability density of x1, x2, x3 and x5
f1245=f1*f2*f4*f5; % joint probability density of x1, x2, x4 and x5
f1345=f1*f3*f4*f5; % joint probability density of x1, x3, x4 and x5
f2345=f2*f3*f4*f5; % joint probability density of x2, x3, x4 and x5

f12345=f1*f2*f3*f4*f5; % joint probability density of x1, x2, x3, x4 and x5
% ----------------------------------------------------------------------------------------------------------------------



% ----------------------------------------------------------------------------------------------------------------------
% DECOMPOSITION OF STRUCTURAL RESPONSE:
% ----------------------------------------------------------------------------------------------------------------------
g0=int(int(int(int(int(g*f12345,x1,-1,1),x2,-1,1),x3,-1,1),x4,-1,1),x5,-1,1); % total mean

g1=int(int(int(int(g*f2345,x2,-1,1),x3,-1,1),x4,-1,1),x5,-1,1)-g0;
g2=int(int(int(int(g*f1345,x1,-1,1),x3,-1,1),x4,-1,1),x5,-1,1)-g0;
g3=int(int(int(int(g*f1245,x1,-1,1),x2,-1,1),x4,-1,1),x5,-1,1)-g0;
g4=int(int(int(int(g*f1235,x1,-1,1),x2,-1,1),x3,-1,1),x5,-1,1)-g0;
g5=int(int(int(int(g*f1234,x1,-1,1),x2,-1,1),x3,-1,1),x4,-1,1)-g0;

g12=int(int(int(g*f345,x3,-1,1),x4,-1,1),x5,-1,1)-(g0+g1+g2);
g13=int(int(int(g*f245,x2,-1,1),x4,-1,1),x5,-1,1)-(g0+g1+g3);
g14=int(int(int(g*f235,x2,-1,1),x3,-1,1),x5,-1,1)-(g0+g1+g4);
g15=int(int(int(g*f234,x2,-1,1),x3,-1,1),x4,-1,1)-(g0+g1+g5);
g23=int(int(int(g*f145,x1,-1,1),x4,-1,1),x5,-1,1)-(g0+g2+g3);
g24=int(int(int(g*f135,x1,-1,1),x3,-1,1),x5,-1,1)-(g0+g2+g4);
g25=int(int(int(g*f134,x1,-1,1),x3,-1,1),x4,-1,1)-(g0+g2+g5);
g34=int(int(int(g*f125,x1,-1,1),x2,-1,1),x5,-1,1)-(g0+g3+g4);
g35=int(int(int(g*f124,x1,-1,1),x2,-1,1),x4,-1,1)-(g0+g3+g5);
g45=int(int(int(g*f123,x1,-1,1),x2,-1,1),x3,-1,1)-(g0+g4+g5);

g123=int(int(g*f45,x4,-1,1),x5,-1,1)-(g0+g1+g2+g3+g12+g13+g23);
g124=int(int(g*f35,x3,-1,1),x5,-1,1)-(g0+g1+g2+g4+g12+g14+g24);
g125=int(int(g*f34,x3,-1,1),x4,-1,1)-(g0+g1+g2+g5+g12+g15+g25);
g134=int(int(g*f25,x2,-1,1),x5,-1,1)-(g0+g1+g3+g4+g13+g14+g34);
g135=int(int(g*f24,x2,-1,1),x4,-1,1)-(g0+g1+g3+g5+g13+g15+g35);
g145=int(int(g*f23,x2,-1,1),x3,-1,1)-(g0+g1+g4+g5+g14+g15+g45);
g234=int(int(g*f15,x1,-1,1),x5,-1,1)-(g0+g2+g3+g4+g23+g24+g34);
g235=int(int(g*f14,x1,-1,1),x4,-1,1)-(g0+g2+g3+g5+g23+g25+g35);
g245=int(int(g*f13,x1,-1,1),x3,-1,1)-(g0+g2+g4+g5+g24+g25+g45);
g345=int(int(g*f12,x1,-1,1),x2,-1,1)-(g0+g3+g4+g5+g34+g35+g45);

g1234=int(g*f5,x5,-1,1)-(g0+g1+g2+g3+g4+g12+g13+g14+g23+g24+g34+g123+g124+g134+g234);
g1235=int(g*f4,x4,-1,1)-(g0+g1+g2+g3+g5+g12+g13+g15+g23+g25+g35+g123+g125+g135+g235);
g1245=int(g*f3,x3,-1,1)-(g0+g1+g2+g4+g5+g12+g14+g15+g24+g25+g45+g124+g125+g145+g245);
g1345=int(g*f2,x2,-1,1)-(g0+g1+g3+g4+g5+g13+g14+g15+g34+g35+g45+g134+g135+g145+g345);
g2345=int(g*f1,x1,-1,1)-(g0+g2+g3+g4+g5+g23+g24+g25+g34+g35+g45+g234+g235+g245+g345);

g12345=g-(g0+g1+g2+g3+g4+g5+g12+g13+g14+g15+g23+g24+g25+g34+g35+g45+g123+g124+g125+g134+g135+g145+g234+g235+g245+g345...
    +g1234+g1235+g1245+g1345+g2345);
% ----------------------------------------------------------------------------------------------------------------------



% ----------------------------------------------------------------------------------------------------------------------
% CHECKING ORTHOGONALITY OF DECOMPOSITION:
% ----------------------------------------------------------------------------------------------------------------------
check1=int(g1*f1,x1,-1,1)==0;
check2=int(g2*f2,x2,-1,1)==0;
check3=int(g3*f3,x3,-1,1)==0;
check4=int(g4*f4,x4,-1,1)==0;
check5=int(g5*f5,x5,-1,1)==0;

check12=int(int(g12*f12,x1,-1,1),x2,-1,1)==0;
check13=int(int(g13*f13,x1,-1,1),x3,-1,1)==0;
check14=int(int(g14*f14,x1,-1,1),x4,-1,1)==0;
check15=int(int(g15*f15,x1,-1,1),x5,-1,1)==0;
check23=int(int(g23*f23,x2,-1,1),x3,-1,1)==0;
check24=int(int(g24*f24,x2,-1,1),x4,-1,1)==0;
check25=int(int(g25*f25,x2,-1,1),x5,-1,1)==0;
check34=int(int(g34*f34,x3,-1,1),x4,-1,1)==0;
check35=int(int(g35*f35,x3,-1,1),x5,-1,1)==0;
check45=int(int(g45*f45,x4,-1,1),x5,-1,1)==0;

check123=int(int(int(g123*f123,x1,-1,1),x2,-1,1),x3,-1,1)==0;
check124=int(int(int(g124*f124,x1,-1,1),x2,-1,1),x4,-1,1)==0;
check125=int(int(int(g125*f125,x1,-1,1),x2,-1,1),x5,-1,1)==0;
check134=int(int(int(g134*f134,x1,-1,1),x3,-1,1),x4,-1,1)==0;
check135=int(int(int(g135*f135,x1,-1,1),x3,-1,1),x5,-1,1)==0;
check145=int(int(int(g145*f145,x1,-1,1),x4,-1,1),x5,-1,1)==0;
check234=int(int(int(g234*f234,x2,-1,1),x3,-1,1),x4,-1,1)==0;
check235=int(int(int(g235*f235,x2,-1,1),x3,-1,1),x5,-1,1)==0;
check245=int(int(int(g245*f245,x2,-1,1),x4,-1,1),x5,-1,1)==0;
check345=int(int(int(g345*f345,x3,-1,1),x4,-1,1),x5,-1,1)==0;

check1234=int(int(int(int(g1234*f1234,x1,-1,1),x2,-1,1),x3,-1,1),x4,-1,1)==0;
check1235=int(int(int(int(g1235*f1235,x1,-1,1),x2,-1,1),x3,-1,1),x5,-1,1)==0;
check1245=int(int(int(int(g1245*f1245,x1,-1,1),x2,-1,1),x4,-1,1),x5,-1,1)==0;
check1345=int(int(int(int(g1345*f1345,x1,-1,1),x3,-1,1),x4,-1,1),x5,-1,1)==0;
check2345=int(int(int(int(g2345*f2345,x2,-1,1),x3,-1,1),x4,-1,1),x5,-1,1)==0;

check12345=int(int(int(int(int(g12345*f12345,x1,-1,1),x2,-1,1),x3,-1,1),x4,-1,1),x5,-1,1)==0;

if ~(check1 && check2 && check3 && check4 && check5...
        && check12 && check13 && check14 && check15 && check23 && check24 && check25 && check34 && check35 && check45...
        && check123 && check124 && check125 && check134 && check135 && check145 && check234 && check235 && check245...
        && check345 && check1234 && check1235 && check1245 && check1345 && check2345 && check12345)
    
    error('Orthogonality check did not pass.')
end
% ----------------------------------------------------------------------------------------------------------------------



% ----------------------------------------------------------------------------------------------------------------------
% DECOMPOSITION OF TOTAL VARIANCE:
% ----------------------------------------------------------------------------------------------------------------------
Dg=int(int(int(int(int(g^2*f12345,x1,-1,1),x2,-1,1),x3,-1,1),x4,-1,1),x5,-1,1)-g0^2; % total variance

D1=int(g1^2*f1,x1,-1,1);
D2=int(g2^2*f2,x2,-1,1);
D3=int(g3^2*f3,x3,-1,1);
D4=int(g4^2*f4,x4,-1,1);
D5=int(g5^2*f5,x5,-1,1);

D12=int(int(g12^2*f12,x1,-1,1),x2,-1,1);
D13=int(int(g13^2*f13,x1,-1,1),x3,-1,1);
D14=int(int(g14^2*f14,x1,-1,1),x4,-1,1);
D15=int(int(g15^2*f15,x1,-1,1),x5,-1,1);
D23=int(int(g23^2*f23,x2,-1,1),x3,-1,1);
D24=int(int(g24^2*f24,x2,-1,1),x4,-1,1);
D25=int(int(g25^2*f25,x2,-1,1),x5,-1,1);
D34=int(int(g34^2*f34,x3,-1,1),x4,-1,1);
D35=int(int(g35^2*f35,x3,-1,1),x5,-1,1);
D45=int(int(g45^2*f45,x4,-1,1),x5,-1,1);

D123=int(int(int(g123^2*f123,x1,-1,1),x2,-1,1),x3,-1,1);
D124=int(int(int(g124^2*f124,x1,-1,1),x2,-1,1),x4,-1,1);
D125=int(int(int(g125^2*f125,x1,-1,1),x2,-1,1),x5,-1,1);
D134=int(int(int(g134^2*f134,x1,-1,1),x3,-1,1),x4,-1,1);
D135=int(int(int(g135^2*f135,x1,-1,1),x3,-1,1),x5,-1,1);
D145=int(int(int(g145^2*f145,x1,-1,1),x4,-1,1),x5,-1,1);
D234=int(int(int(g234^2*f234,x2,-1,1),x3,-1,1),x4,-1,1);
D235=int(int(int(g235^2*f235,x2,-1,1),x3,-1,1),x5,-1,1);
D245=int(int(int(g245^2*f245,x2,-1,1),x4,-1,1),x5,-1,1);
D345=int(int(int(g345^2*f345,x3,-1,1),x4,-1,1),x5,-1,1);

D1234=int(int(int(int(g1234^2*f1234,x1,-1,1),x2,-1,1),x3,-1,1),x4,-1,1);
D1235=int(int(int(int(g1235^2*f1235,x1,-1,1),x2,-1,1),x3,-1,1),x5,-1,1);
D1245=int(int(int(int(g1245^2*f1245,x1,-1,1),x2,-1,1),x4,-1,1),x5,-1,1);
D1345=int(int(int(int(g1345^2*f1345,x1,-1,1),x3,-1,1),x4,-1,1),x5,-1,1);
D2345=int(int(int(int(g2345^2*f2345,x2,-1,1),x3,-1,1),x4,-1,1),x5,-1,1);

D12345=int(int(int(int(int(g12345^2*f12345,x1,-1,1),x2,-1,1),x3,-1,1),x4,-1,1),x5,-1,1);
% ----------------------------------------------------------------------------------------------------------------------



% ----------------------------------------------------------------------------------------------------------------------
% CONTRIBUTIONS TO TOTAL VARIANCE:
% ----------------------------------------------------------------------------------------------------------------------
S1=D1/Dg; % contribution to total variance due to x1
S2=D2/Dg; % contribution to total variance due to x2
S3=D3/Dg; % contribution to total variance due to x3
S4=D4/Dg; % contribution to total variance due to x4
S5=D5/Dg; % contribution to total variance due to x5

S12=D12/Dg; % contribution to total variance due to the interaction between x1 and x2
S13=D13/Dg; % contribution to total variance due to the interaction between x1 and x3
S14=D14/Dg; % contribution to total variance due to the interaction between x1 and x4
S15=D15/Dg; % contribution to total variance due to the interaction between x1 and x5
S23=D23/Dg; % contribution to total variance due to the interaction between x2 and x3
S24=D24/Dg; % contribution to total variance due to the interaction between x2 and x4
S25=D25/Dg; % contribution to total variance due to the interaction between x2 and x5
S34=D34/Dg; % contribution to total variance due to the interaction between x3 and x4
S35=D35/Dg; % contribution to total variance due to the interaction between x3 and x5
S45=D45/Dg; % contribution to total variance due to the interaction between x4 and x5

S123=D123/Dg; % contribution to total variance due to the interaction between x1, x2 and x3
S124=D124/Dg; % contribution to total variance due to the interaction between x1, x2 and x4
S125=D125/Dg; % contribution to total variance due to the interaction between x1, x2 and x5
S134=D134/Dg; % contribution to total variance due to the interaction between x1, x3 and x4
S135=D135/Dg; % contribution to total variance due to the interaction between x1, x3 and x5
S145=D145/Dg; % contribution to total variance due to the interaction between x1, x4 and x5
S234=D234/Dg; % contribution to total variance due to the interaction between x2, x3 and x4
S235=D235/Dg; % contribution to total variance due to the interaction between x2, x3 and x5
S245=D245/Dg; % contribution to total variance due to the interaction between x2, x4 and x5
S345=D345/Dg; % contribution to total variance due to the interaction between x3, x4 and x5

S1234=D1234/Dg; % contribution to total variance due to the interaction between x1, x2, x3 and x4
S1235=D1235/Dg; % contribution to total variance due to the interaction between x1, x2, x3 and x5
S1245=D1245/Dg; % contribution to total variance due to the interaction between x1, x2, x4 and x5
S1345=D1345/Dg; % contribution to total variance due to the interaction between x1, x3, x4 and x5
S2345=D2345/Dg; % contribution to total variance due to the interaction between x2, x3, x4 and x5

S12345=D12345/Dg; % contribution to total variance due to the interaction between x1, x2, x3, x4 and x5
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
function g=getStructuralResponse(x1,x2,x3,x4,x5)
p1=legendreP(2,x1)+jacobiP(3,3,5,x1); % a univariate function defined on x1 = [-1,1]
p2=legendreP(5,x2)+jacobiP(2,2,7,x2); % a univariate function defined on x2 = [-1,1]
p3=legendreP(3,x3)+jacobiP(2,4,1,x3); % a univariate function defined on x3 = [-1,1]
p4=legendreP(4,x4)+jacobiP(5,3,2,x4); % a univariate function defined on x4 = [-1,1]
p5=legendreP(3,x5)+jacobiP(4,5,4,x5); % a univariate function defined on x5 = [-1,1]

% Structural response of interest (notice that this g is a function of x1, x2, x3, x4 and x5):
g=p1*p2*p3*p4*p5/1e5; % a pentavariate function defined on x1 x x2 x x3 x x4 x x5 = [-1,1]^5
end
% ----------------------------------------------------------------------------------------------------------------------
