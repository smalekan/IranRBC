//********** preamble block
var A c k y i r w h;
varexo e;
parameters psi beta alpha delta gamma rho sigmae;
psi =2.3;
beta = 0.94;
alpha = 0.66;
delta = 0.0425;
rho= 0.9;
sigmae = 0.01;
//********** model block
model;
1/exp(c)= beta*(exp(r) + 1 - delta)/exp(c(+1));
exp(w) = psi*exp(c);
exp(k) = (1-delta)*exp(k(-1))+ exp(i);
exp(y) = A*(exp(k(-1))^alpha)*(exp(h)^(1-alpha));
exp(y) = exp(c) + exp(i);
exp(r) = A*alpha*(exp(k(-1))^(alpha-1))*(exp(h)^(1-alpha));
exp(w) = A*(1-alpha)*(exp(k(-1))^alpha)*(exp(h)^(-alpha));
A = rho*A(-1) + e;
end;
//******** steady state or initial value block
initval;
A=exp(1);
r =log( 1/beta - 1 + delta);
w = log((1 - alpha)*A*(alpha*A/(1/beta - 1 + delta))^(alpha/(1-alpha)));
h=log(0.5);
y=log(A^(1/(1-alpha))*(alpha/exp(r))^(alpha/(1-alpha))*exp(h));
k = log(alpha*exp(y)/exp(r));
i =log(delta*exp(k));
c =log(exp(y) - exp(i));
end;
//*************** shocks block
shocks;
var e = sigmae^2;
end;
//************ computation block
steady;
stoch_simul(hp_filter=1600,periods = 1000,drop=100,irf=75,order=1);