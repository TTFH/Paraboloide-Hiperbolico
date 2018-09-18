figure;
hold on;
% Punto de origen:
u0 = 1; v0 = 0;
% Direccion de inicio:
p0 = 1; q0 = 0;
% Cant. Iteraciones:
N = 1000;
% Distancia sobre el plano al punto de destino:
d = 5;
h = d / N;

% Hallar solucion de octave
y0 = [u0; v0; p0; q0];
f = @(s, y) [ % y'(u', v', p', q')
  y(3);
  y(4);
  y(3) - 2 * tan(y(2)) * y(3) * y(4);
  y(4) + sin(y(2)) * cos(y(2)) * y(3) * y(3);
];
s = [h: h: d];
[s, qa] = ode45(f, s, y0);
u = qa(:, 1);
v = qa(:, 2);

% Geodesica (Solucion de octave)
plot3(cos(u) .* cos(v), sin(u) .* cos(v), sin(v), '-k', 'lineWidth', 2);

% Superficie (Esfera)
colormap(jet(256));
[X,Y] = meshgrid(-25: 0.5: 25);
surf(cos(X) .* cos(Y), sin(X) .* cos(Y), sin(Y));
shading interp;
