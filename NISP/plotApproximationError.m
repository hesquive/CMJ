function plotApproximationError(fact,fapprox)
% Hugo Esquivel, 2021
% -

fact=matlabFunction(fact);
fapprox=matlabFunction(fapprox);

numPoints=50;

x1=linspace(-1,1,numPoints); % x1 = [-1,1]
x2=linspace(-1,1,numPoints); % x2 = [-1,1]

[x1grid,x2grid]=meshgrid(x1,x2);
factgrid=fact(x1grid,x2grid);
fapproxgrid=fapprox(x1grid,x2grid);

epsgrid=factgrid-fapproxgrid;

figure
surf(x1grid,x2grid,epsgrid)
colormap(jet)
title('\bfseries{Approximation Error, $\epsilon$}','Interpreter','LaTeX')
xlabel('$x_1$','Interpreter','LaTeX')
ylabel('$x_2$','Interpreter','LaTeX')
zlabel('$\epsilon(x_1,x_2)$','Interpreter','LaTeX')

figure
surf(x1grid,x2grid,factgrid)
hold on
surf(x1grid,x2grid,fapproxgrid)
colormap(jet)
title('\bfseries{Actual and Approximate Structural Response, $f_\mathrm{act},\,f_\mathrm{approx}$}','Interpreter','LaTeX')
xlabel('$x_1$','Interpreter','LaTeX')
ylabel('$x_2$','Interpreter','LaTeX')
zlabel('$f_\mathrm{act}(x_1,x_2),\,f_\mathrm{approx}(x_1,x_2)$','Interpreter','LaTeX')
end
