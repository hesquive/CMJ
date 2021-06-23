function plotApproximateStructuralResponse(fapprox)
% Hugo Esquivel, 2021
% -

fapprox=matlabFunction(fapprox);

numPoints=50;

x1=linspace(-1,1,numPoints); % x1 = [-1,1]
x2=linspace(-1,1,numPoints); % x2 = [-1,1]

[x1grid,x2grid]=meshgrid(x1,x2);
fapproxgrid=fapprox(x1grid,x2grid);

figure
surf(x1grid,x2grid,fapproxgrid)
colormap(jet)
title('\bfseries{Approximate Structural Response, $f_\mathrm{approx}$}','Interpreter','LaTeX')
xlabel('$x_1$','Interpreter','LaTeX')
ylabel('$x_2$','Interpreter','LaTeX')
zlabel('$f_\mathrm{approx}(x_1,x_2)$','Interpreter','LaTeX')
end
