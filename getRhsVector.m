% This function calculates the right hand side vector R
function R = getRhsVector(A, U0, step)
global n ne xL xR dt tau sigma

% For the first step with the exact initial 
% condition, there is no need to multiply to A
if step == 0
  A = eye(n*ne);
end

%%%% calculate g(v) and K(z) %%%%
g1 = zeros(n, 1);
 
for i = 1:n
  g1(i) = guxL(step*dt) * phi(i, 1, xL);
end

G = zeros(n*ne, 1);
K = zeros(n*ne, 1);

G(1:n) = tau * g1;

K(1:n) = -g1;

%%%% calculate h(w) %%%%
h1 = zeros(n, 1);

for i = 1:n
  h1(i) = (-1-sigma)*gqxL(step*dt)*phi(i, 1, xL);
end

H = zeros(n*ne, 1);
H(1:n) = h1;

%%%% calculate s(nu) %%%%
  S = zeros(ne,1);
  S(ne,1) = gpxR(step*dt);

%%%% calculate the right hand side vector %%%% 
m=3*n*ne+2*ne;
R = zeros(m, 1);

R(1:n*ne) = -G + (1/dt)*A*U0;
R(n*ne+1:2*n*ne) = H;
R(2*n*ne+1:3*n*ne) = K;
R(3*n*ne+ne+1:m) = S;