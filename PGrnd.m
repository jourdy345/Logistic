% Polyá-Gamma(1,c) random variable sampler
%
% DESCRIPTION
%
%    Generates PG(1,c) random variables
%    The algorithm is credited to Devroye
%
% REFERENCE
%
%    "Non-Uniform Random Variate Generation, Luc Devroye, Springer-Verlag, 1986,
%    the University of California, 16 Dec 2010, ISBN:0387963057, 9780387963051"
function [X] = PGrnd(num,n,z)
  z = z*ones(1,num);
  n = n*ones(1,num);
  X = zeros(1,num);
  total_trials = 0;
  for i=1:num
    for j=1:n(i)
      [tmp,ntrials] = PG1(z(i));
      X(i) = X(i)+tmp;
      total_trials = total_trials+ntrials;
    end
  end
end

function [X,num_trials] = PG1(Z)
  PI = 3.141592653589793238462643383280;
  Z = abs(Z)*0.5;
  fz = 0.125*PI*PI+Z*Z*0.5;
  num_trials = 0;
  total_iter = 0;
  while 1
    num_trials = num_trials+1;
    if rand() < mass_texpon(Z)
      X = 0.64+exprnd(1,1,1)/fz;
    else
      X = rtigauss(Z);
    end

    S = a_coef(0,X);
    Y = rand()*S;
    n = 0;

    while 1
      n = n+1;
      total_iter = total_iter+1;
      if mod(n,2)==1
        S = S-a_coef(n,X);
        if Y <= S
          break
        end
      else
        S = S+a_coef(n,X);
        if Y > S
          break
        end
      end
    end
    if Y <= S
      break
    end
  end
  X = 0.25*X;
end

function [out] = mass_texpon(Z)
  TRUNC = 0.64;
  x = TRUNC;
  PI = 3.141592653589793238462643383280;
  fz = 0.125*PI*PI+0.5*Z*Z;
  b = sqrt(1/x)*(x*Z-1);
  a = -sqrt(1/x)*(x*Z+1);
  x0 = log(fz)+fz*TRUNC;
  xb = x0-Z+pnorm(b,0,1,1,1);
  xa = x0+Z+pnorm(a,0,1,1,1);
  qdivp = 4/PI*(exp(xb)+exp(xa));
  out = 1/(1+qdivp);
end

function [X] = rtigauss(Z,R)
  if nargin < 2
    R = 0.64;
  end
  Z = abs(Z);
  mu = 1/Z;
  X = R+1;
  if mu > R
    alph = 0;
    while rand() > alph
      E = exprnd(1,1,2);
      while E(1)^2 > 2*E(2)/R
        E = exprnd(1,1,2);
      end
      X = R/(1+R*E(1))^2;
      alph = exp(-0.5*Z*Z*X);
    end
  else
    while X > R
      lambda = 1;
      Y = randn()^2;
      X = mu+0.5*mu*mu/lambda*Y-0.5*mu/lambda*sqrt(4*mu*lambda*Y+(mu*Y)^2);
      if rand() > mu/(mu+X)
        X = mu*mu/X;
      end
    end
  end
end

function [out] = a_coef(n,x)
  PI = 3.141592653589793238462643383280;
  if x > 0.64
    out = PI*(n+0.5)*exp(-(n+0.5)^2*PI^2*x/2);
  else
    out = (2/PI/x)^1.5*PI*(n+0.5)*exp(-2*(n+0.5)^2/x);
  end
end