
n=512
ak = 0.35
k = 2*pi
a = ak/k
g = 9.81
h = 0.5

x = linspace(-0.5,0.5,n);
y = linspace(-0.5,0.5,n);
[X,Y] = meshgrid(x,y);

f = StokesWaveSurface(X, Y, a, k, h);
[u, v] = StokesWaveVelocities(X, Y, a, k, g, h);

l = ones(n,n)*9;
p = zeros(n,n);

output_matrix_double(l,n,x,y,'test_input_matrix_l.bin')
output_matrix_double(p,n,x,y,'test_input_matrix_p.bin')
output_matrix_double((f > Y),n,x,y,'test_input_matrix_f.bin')
output_matrix_double(u.*(f > Y),n,x,y,'test_input_matrix_u.bin')
output_matrix_double(v.*(f > Y),n,x,y,'test_input_matrix_v.bin')

subplot(1,3,1)
contourf(x,y,(f > Y))
hold on 
streamslice(x,y,v.*(f > Y),u.*(f > Y))
hold off
subplot(1,3,2)
contourf(x,y,u.*(f > Y))
subplot(1,3,3)
contourf(x,y,v.*(f > Y))

function [u, v] = StokesWaveVelocities(x, y, a, k, g, h)
  % StokesWaveVelocities - Calculate velocities for a Stokes wave
  %
  % Syntax:
  %   [u, v] = StokesWaveVelocities(x, y, a, k, g, h)
  %
  % Inputs:
  %   x - Horizontal position (array)
  %   y - Vertical position (array)
  %   a - Wave amplitude
  %   k - Wavenumber
  %   g - Acceleration due to gravity
  %   h - Water depth
  %
  % Outputs:
  %   u - Horizontal velocity
  %   v - Vertical velocity

  % Calculate alpha and sigma
  alpha = 1 / tanh(k * h);
  sigma = sqrt(g * k * tanh(k * h) * (1 + k^2 * a^2 * (9/8 * (alpha^2 - 1)^2 + alpha^2)));
  A = a * g / sigma;

  % Horizontal velocity
  u1 = A * cosh(k * (y + h)) / cosh(k * h) * k * cos(k * x);
  u2 = a * k * 3/8 * A / alpha * (alpha^2 - 1)^2 * cosh(2 * k * (y + h)) * 2 * k * cos(2 * k * x) / cosh(2 * k * h);
  u3 = (a * k)^2 * 1/64 * (alpha^2 - 1) * (alpha^2 + 3) * (9 * alpha^2 - 13) * ...
    cosh(3 * k * (y + h)) / cosh(3 * k * h) * a^2 * k^2 * A * 3 * k * cos(3 * k * x);
  u = u1 + u2 + u3;

  % Vertical velocity
  v1 = A * k * sinh(k * (y + h)) / cosh(k * h) * sin(k * x);
  v2 = a * k * 3/8 * A / alpha * (alpha^2 - 1)^2 * 2 * k * sinh(2 * k * (y + h)) * sin(2 * k * x) / cosh(2 * k * h);
  v3 = (a * k)^2 * 1/64 * (alpha^2 - 1) * (alpha^2 + 3) * (9 * alpha^2 - 13) * ...
    3 * k * sinh(3 * k * (y + h)) / cosh(3 * k * h) * a^2 * k^2 * A * sin(3 * k * x);
  v = v1 + v2 + v3;
end



function eta = StokesWaveSurface(x, y, a, k, h)
  % StokesWaveSurface - Calculate the surface elevation for a Stokes wave
  %
  % Syntax:
  %   eta = Wave(x, y, a, k, g, h)
  %
  % Inputs:
  %   x - Horizontal position (array)
  %   y - Vertical position (array)
  %   a - Wave amplitude
  %   k - Wavenumber
  %   h - Water depth
  %
  % Outputs:
  %   eta - Surface elevation

  % First-order term
  eta1 = a * cos(k * x);

  % Calculate alpha
  alpha = 1 / tanh(k * h);

  % Second-order term
  eta2 = 1/4 * alpha * (3 * alpha^2 - 1) * a^2 * k * cos(2 * k * x);

  % Third-order term
  eta3 = -3/8 * (alpha^4 - 3 * alpha^2 + 3) * a^3 * k^2 * cos(k * x) + ...
       3/64 * (8 * alpha^6 + (alpha^2 - 1)^2) * a^3 * k^2 * cos(3 * k * x);

  % Combine terms to get the surface elevation
  eta = eta1 + a * k * eta2 + (a * k)^2 * eta3 - y;
end
