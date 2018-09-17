figure;
hold on;
% Punto de origen:
u0 = -20; v0 = -15;
% Direccion de inicio:
p0 = 1; q0 = 5;
% Cant. Iteraciones:
N = 100;
% Distancia sobre el plano al punto de destino:
d = 5.2;

h = d / N;
y = zeros(N, 4); % y(u, v, p, q)
y(1, :) = [u0, v0, p0, q0];

u = y(:, 1);
v = y(:, 2);
p = y(:, 3);
q = y(:, 4);
% Euler hacia adelante
for i = 1:N-1
  u(i+1) = u(i) + h * p(i);
  v(i+1) = v(i) + h * q(i);
  p(i+1) = p(i) + h * (-2*v(i)*p(i)*q(i)) /...
  	(u(i)*u(i) + v(i)*v(i) + 1);
  q(i+1) = q(i) + h * (-2*u(i)*p(i)*q(i)) /...
  	(u(i)*u(i) + v(i)*v(i) + 1);
endfor
% Geodesica (Euler hacia adelante)
plot3(u, v, u .* v, '-b', 'lineWidth', 1);

% Hallar solucion de octave
y0 = [u0; v0; p0; q0];
f = @(s, y) [ % y'(u', v', p', q')
  y(3);
  y(4);
  (-2 * y(2) * y(3) * y(4)) / (y(1) * y(1) + y(2) * y(2) + 1);
  (-2 * y(1) * y(3) * y(4)) / (y(1) * y(1) + y(2) * y(2) + 1);
];
s = [h: h: d];
[s, qa] = ode45(f, s, y0);
u_octave = qa(:, 1);
v_octave = qa(:, 2);
% Geodesica (Solucion de octave)
plot3(u_octave, v_octave, u_octave .* v_octave, '-r', 'lineWidth', 1);

% Superficie (Paraboloide Hiperbolico)
colormap(jet(256));
[X,Y] = meshgrid(-25: 0.5: 25);
surf(X, Y, X .* Y);
shading interp;

% valor_real -  estimado
ErAbsU = abs(u_octave - u);
Error_Absoluto_U = sum(ErAbsU) / N
ErAbsV = abs(v_octave - v);
Error_Absoluto_V = sum(ErAbsV) / N

% (valor_real -  estimado) / valor_real x100to
ErRelU = abs(u_octave - u) ./ abs(u_octave);
Error_Relativo_U = 100 * sum(ErRelU) / N
ErRelV = abs(v_octave - v) ./ abs(v_octave);
Error_Relativo_V = 100 * sum(ErRelV) / N
