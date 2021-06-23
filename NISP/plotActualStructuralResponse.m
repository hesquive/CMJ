function plotActualStructuralResponse(fact)
% Hugo Esquivel, 2021
% -

fact=matlabFunction(fact);

numPoints=50;

x1=linspace(-1,1,numPoints); % x1 = [-1,1]
x2=linspace(-1,1,numPoints); % x2 = [-1,1]

[x1grid,x2grid]=meshgrid(x1,x2);
factgrid=fact(x1grid,x2grid);

figure
surf(x1grid,x2grid,factgrid)
colormap(jet)
title('\bfseries{Actual Structural Response, $f_\mathrm{act}$}','Interpreter','LaTeX')
xlabel('$x_1$','Interpreter','LaTeX')
ylabel('$x_2$','Interpreter','LaTeX')
zlabel('$f_\mathrm{act}(x_1,x_2)$','Interpreter','LaTeX')
end
