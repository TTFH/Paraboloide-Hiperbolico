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

for i = 1:N-1 % Iterar de 2 a N
  u(i+1) = u(i) + h * p(i);
  v(i+1) = v(i) + h * q(i);
  p(i+1) = p(i) + h * (-2*v(i)*p(i)*q(i)) /...
  	(u(i)*u(i) + v(i)*v(i) + 1);
  q(i+1) = q(i) + h * (-2*u(i)*p(i)*q(i)) /...
  	(u(i)*u(i) + v(i)*v(i) + 1);
endfor
% ".*" es el producto coordenada a coordenada
plot3(u, v, u .* v, '-b', 'lineWidth', 1);

% ------------------------
q0 = [-20; -15; 1; 5];
f = @(t, q) [
q(3);
q(4);
(-2 * q(2) * q(3) * q(4)) / (q(1) * q(1) + q(2) * q(2) + 1);
(-2 * q(1) * q(3) * q(4)) / (q(1) * q(1) + q(2) * q(2) + 1);
];

t = [0: 0.01: 5.2];
[t, qa] = ode45(f, t, q0);
u = qa(:, 1);
v = qa(:, 2);

plot3(u, v, u .* v, '-r', 'lineWidth', 1);
% ------------------------

colormap(jet(256));
[X,Y] = meshgrid(-25:.5:25);
surf(X, Y, X .* Y);
shading interp;


