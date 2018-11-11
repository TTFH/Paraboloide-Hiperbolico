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

du = @(s, u, v, p, q) p;
dv = @(s, u, v, p, q) q;
dp = @(s, u, v, p, q) (-2 * v * p * q) / (u * u + v * v + 1);
dq = @(s, u, v, p, q) (-2 * u * p * q) / (u * u + v * v + 1);

% Runge-Kutta
for i = 1:N-1
  if (i <= 7)
    salida1(i,1) = u(i);
    salida1(i,2) = v(i);
  endif

  uk1 = du(i*h, u(i), v(i), p(i), q(i));
  vk1 = dv(i*h, u(i), v(i), p(i), q(i));
  pk1 = dp(i*h, u(i), v(i), p(i), q(i));
  qk1 = dq(i*h, u(i), v(i), p(i), q(i));

  uk2 = du(i*h + (h/2), u(i) + (h/2)*uk1, v(i) + (h/2)*vk1, p(i) + (h/2)*pk1, q(i) + (h/2)*qk1);
  vk2 = dv(i*h + (h/2), u(i) + (h/2)*uk1, v(i) + (h/2)*vk1, p(i) + (h/2)*pk1, q(i) + (h/2)*qk1);
  pk2 = dp(i*h + (h/2), u(i) + (h/2)*uk1, v(i) + (h/2)*vk1, p(i) + (h/2)*pk1, q(i) + (h/2)*qk1);
  qk2 = dq(i*h + (h/2), u(i) + (h/2)*uk1, v(i) + (h/2)*vk1, p(i) + (h/2)*pk1, q(i) + (h/2)*qk1);

  uk3 = du(i*h + (h/2), u(i) + (h/2)*uk2, v(i) + (h/2)*vk2, p(i) + (h/2)*pk2, q(i) + (h/2)*qk2);
  vk3 = dv(i*h + (h/2), u(i) + (h/2)*uk2, v(i) + (h/2)*vk2, p(i) + (h/2)*pk2, q(i) + (h/2)*qk2);
  pk3 = dp(i*h + (h/2), u(i) + (h/2)*uk2, v(i) + (h/2)*vk2, p(i) + (h/2)*pk2, q(i) + (h/2)*qk2);
  qk3 = dq(i*h + (h/2), u(i) + (h/2)*uk2, v(i) + (h/2)*vk2, p(i) + (h/2)*pk2, q(i) + (h/2)*qk2);

  uk4 = du(i*h + (h/2), u(i) + h*uk3, v(i) + h*vk3, p(i) + h*pk3, q(i) + h*qk3);
  vk4 = dv(i*h + (h/2), u(i) + h*uk3, v(i) + h*vk3, p(i) + h*pk3, q(i) + h*qk3);
  pk4 = dp(i*h + (h/2), u(i) + h*uk3, v(i) + h*vk3, p(i) + h*pk3, q(i) + h*qk3);
  qk4 = dq(i*h + (h/2), u(i) + h*uk3, v(i) + h*vk3, p(i) + h*pk3, q(i) + h*qk3);

  u(i+1) = u(i) + (h/6) * (uk1 + 2*uk2 + 2*uk3 + uk4);
  v(i+1) = v(i) + (h/6) * (vk1 + 2*vk2 + 2*vk3 + vk4);
  p(i+1) = p(i) + (h/6) * (pk1 + 2*pk2 + 2*pk3 + pk4);
  q(i+1) = q(i) + (h/6) * (qk1 + 2*qk2 + 2*qk3 + qk4);
endfor

% Geodesica (Runge-Kutta)
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
