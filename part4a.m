figure;
hold on;

x0 = -20; % Punto de origen
y0 = -15;
v1 = 1; % Direccion de inicio
v2 = 5;
d = 5.2; % Distancia sobre el plano al punto de destino
N = 100; % Cant. Iteraciones

h = d / N;
y = zeros(N, 4); % u, v, p, q
y(1, :) = [x0, y0, v1, v2];

u = y(:, 1);
v = y(:, 2);
p = y(:, 3);
q = y(:, 4);

for i = 1:N-1 % va de 2 a N
  u(i+1) = u(i) + h * p(i);
  v(i+1) = v(i) + h * q(i);
  p(i+1) = p(i) + h * (-2*v(i)*p(i)*q(i)) /...
  	(u(i)*u(i) + v(i)*v(i) + 1);
  q(i+1) = q(i) + h * (-2*u(i)*p(i)*q(i)) /...
  	(u(i)*u(i) + v(i)*v(i) + 1);
endfor

plot3(u, v, u .* v, '-k', 'lineWidth', 1);

colormap(jet(256));
[X,Y] = meshgrid(-25:.5:25);
surf(X, Y, X .* Y);
shading interp;
