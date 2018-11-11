figure;
hold on;
u0 = 0; v0 = 0;
p0 = 0.7; q0 = 0.71;

N = 100;
h = 0.1;
d = h * N;

%h = 1;
%d = 10;
%N = d / h;

y = zeros(N, 4); % y(u, v, p, q)
y(1, :) = [u0, v0, p0, q0];

u = y(:, 1);
v = y(:, 2);
p = y(:, 3);
q = y(:, 4);

% Metodo del Trapecio
for i = 1:N-1
  if (i <= 7)
    salida1(i,1) = u(i);
    salida1(i,2) = v(i);
  endif
  % Predictor
  u(i+1) = u(i) + h * p(i);
  v(i+1) = v(i) + h * q(i);
  p(i+1) = p(i) + h * (-2*v(i)*p(i)*q(i)) /...
    (u(i)*u(i) + v(i)*v(i) + 1);
  q(i+1) = q(i) + h * (-2*u(i)*p(i)*q(i)) /...
    (u(i)*u(i) + v(i)*v(i) + 1);
  % Corrector
  j = 1;
  eps = 0.01;
  while (((h/2)*abs(p(i)) > eps) && ((h/2)*abs(p(i)) > eps) && (j < i))
    u(i+1) = u(i) + (h/2) * (p(i) + p(i+1));
    v(i+1) = v(i) + (h/2) * (q(i) + q(i+1));
    p(i+1) = p(i) + (h/2) * ((-2*v(i)*p(i)*q(i)) /...
      (u(i)*u(i) + v(i)*v(i) + 1) + (-2*v(i+1)*p(i+1)*q(i+1)) /...
      (u(i+1)*u(i+1) + v(i+1)*v(i+1) + 1));
    q(i+1) = q(i) + (h/2) * ((-2*u(i)*p(i)*q(i)) /...
      (u(i)*u(i) + v(i)*v(i) + 1) + (-2*u(i+1)*p(i+1)*q(i+1)) /...
      (u(i+1)*u(i+1) + v(i+1)*v(i+1) + 1));
    j = j + 1;
  endwhile
  iteraciones(i) = j;
endfor

% Geodesica (Metodo del Trapecio)
plot3(u, v, u .* v, '-b', 'lineWidth', 1);

% Hallar solucion de octave
y0 = [u0; v0; p0; q0];
f = @(s, y) [ % y'(u', v', p', q')
  y(3);
  y(4);
  (-2 * y(2) * y(3) * y(4)) / (y(1) * y(1) + y(2) * y(2) + 1);
  (-2 * y(1) * y(3) * y(4)) / (y(1) * y(1) + y(2) * y(2) + 1);
];
s = [0: h: h*(N-1)];
[s, qa] = ode45(f, s, y0);
u_octave = qa(:, 1);
v_octave = qa(:, 2);

% Geodesica (Solucion de octave)
plot3(u_octave, v_octave, u_octave .* v_octave, '-r', 'lineWidth', 1);

% Superficie (Paraboloide Hiperbolico)
colormap(jet(256));
[X,Y] = meshgrid(-d: d / 64: d);
surf(X, Y, X .* Y);
shading interp;

abs(u(N) - u_octave(N))
abs(v(N) - v_octave(N))
