function plotJointProbabilityDistribution(probabilityInfo)
% Hugo Esquivel, 2021
% -

numRandomVariables=length(probabilityInfo.name);

if numRandomVariables==2
    numPoints=50;
    
    s=getProbabilitySupports(probabilityInfo);
    
    x1=linspace(s(1,1),s(1,2),numPoints);
    x2=linspace(s(2,1),s(2,2),numPoints);
    
    [x1grid,x2grid]=meshgrid(x1,x2);
    jpdfgrid=getJointPDF(x1grid,x2grid,probabilityInfo);
    
    figure
    surf(x1grid,x2grid,jpdfgrid)
    colormap(turbo)
    title('\bfseries{Joint Probability Distribution, $f$}','Interpreter','LaTeX')
    xlabel('$x_1$','Interpreter','LaTeX')
    ylabel('$x_2$','Interpreter','LaTeX')
    zlabel('$f(x_1,x_2)$','Interpreter','LaTeX')
end
end

function z=getJointPDF(x,y,probabilityInfo)
syms x1 x2

jointPDF=matlabFunction(probabilityInfo.pdf{1}(x1)*probabilityInfo.pdf{2}(x2),'Vars',[x1,x2]);

z=jointPDF(x,y)  +0*x+0*y;
end
