% Kolmogorov-Smirnov random variable sampler
%
% DESCRIPTION
%
%    Generates Kolmogorov-Smirnov random variables
%    The algorithm is credited to Devroye
%
% REFERENCE
%
%    "Non-Uniform Random Variate Generation, Luc Devroye, Springer-Verlag, 1986,
%    the University of California, 16 Dec 2010, ISBN:0387963057, 9780387963051"

function [out] = KSrnd(n)
  t = 0.75;
  PI = 3.141592653589793238462643383280;
  tprime = PI*PI/(8*t*t);
  p = 0.373;
  q = 1-p;
  out = zeros(n,1);
  for i=1:n
    X = 0;
    if rand() < exp(log(p)-log(p+q))
      while 1
        accept = 0;
        G = 0;
        while ~accept
          E = exprnd(1,1,2);
          E(1) = E(1)/(1-0.5/tprime);
          E(2) = E(2)*2;
          G = tprime+E(1);
          accept = (E(1)*E(1) <= tprime*E(2)*(G+tprime));
          if ~accept
            accept = (G/tprime-1-log(G/tprime) <= E(2));
          end
        end
        X = PI/sqrt(8*G);
        W = 0;
        Z = 0.5/G;
        P = exp(-G);
        n = 1;
        Q = 1;
        U = rand();
        go = 1;
        success = 0;
        while go
          W = W+Z*Q;
          if U >= W
            out(i) = X;
            success = 1;
            break
          end
          n = n+2;
          Q = P^(n*n-1);
          W = W-n*n*Q;
          go = U >= W;
        end
        if success
          break
        end
      end
    else
      while 1
        E = exprnd(1);
        U = rand();
        X = sqrt(t*t+0.5*E);
        W = 0;
        n = 1;
        Z = exp(-2*X*X);
        go = 1;
        success = 0;
        while go
          n = n+1;
          W = W+n*n*Z^(n*n-1);
          if U >= W
            out(i) = X;
            success = 1;
            break
          end
          n = n+1;
          W = W-n*n*Z^(n*n-1);
          go = U >= W;
        end
        if success
          break
        end
      end
    end
  end
end
