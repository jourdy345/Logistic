function [out] = tnormrnd(mu,sig,low,high)
  if low > high
    error('The lower bound should be smaller than the upper bound.')
  end
  draw = 0;
  typ = 0;
  c_mu = mu;
  c_sig = sig;
  c_low = low;
  c_high = high;
  c_stdlow = (c_low-c_mu)/c_sig;
  c_stdhigh = (c_high-c_mu)/c_sig;

  if 0 <= c_stdhigh && 0 > c_stdlow
    typ = 1;
  end
  if 0 < c_stdlow && c_stdhigh == Inf
    typ = 2;
  end
  if 0 > c_stdhigh && c_stdlow == -Inf
    typ = 3;
  end
  if (0 > c_stdhigh || 0 < c_stdlow) && ~(isinf(c_stdhigh) || isinf(c_stdlow))
    typ = 4;
  end

  if typ == 1
    z = 0;
    valid = 0;
    while ~valid
      z = randn();
      if c_stdlow <= z && z <= c_stdhigh
        valid = 1;
      end
    end
    draw = z;
  end
  if typ == 3
    c_stdlow = -c_stdhigh;
    c_stdhigh = Inf;
    c_sig = -c_sig;
    typ = 2;
  end
  if typ == 2
    alphastar = (c_stdlow+sqrt(c_stdlow*c_stdlow+4))*0.5;
    alph = alphastar;
    e_ = 0;
    z = 0;
    rho = 0;
    u = 0;
    valid = 0;
    while ~valid
      e_ = exprnd(1,1);
      z = c_stdlow+e_/alph;
      rho = exp((-alph-z)^2*0.5);
      u = rand();
      if u <= rho
        valid = 1;
      end
    end
    draw = z;
  end
  if typ == 4
    val1 = (2*sqrt(exp(1)))/(c_stdlow+sqrt(c_stdlow*c_stdlow+4));
    val2 = exp((c_stdlow*c_stdlow-c_stdlow*sqrt(c_stdlow*c_stdlow+4))*0.25);
    if c_stdhigh > c_stdlow+val1*val2
      valid = 0;
      while ~valid
        alphastar = (c_stdlow+sqrt(c_stdlow*c_stdlow+4))*0.5;
        alph = alphastar;
        e_ = 0;
        rho = 0;
        u = 0;
        valid1 = 0;
        while ~valid1
          e_ = exprnd(1,1);
          draw = c_stdlow+e_/alph;
          rho = exp((-alph-draw)^2*0.5);
          u = rand();
          if u <= rho
            valid1 = 1;
          end
        end
        if draw <= c_stdhigh
          valid = 1;
        end
      end
    else
      valid = 0;
      rho = 0;
      z = 0;
      u = 0;
      while ~valid
        z = (c_stdhigh-c_stdlow)*rand()+c_stdlow;
        if 0 < c_stdlow
          rho = exp((c_stdlow*c_stdlow-z*z)*0.5);
        elseif c_stdhigh < 0
          rho = exp((c_stdhigh*c_stdhigh-z*z)*0.5);
        elseif 0 < c_stdhigh && c_stdlow < 0
          rho = exp(-z*z*0.5);
        end
        u = rand();
        if u <= rho
          valid = 1;
        end
      end
      draw = z;
    end
  end
  out = c_mu+c_sig*draw;
end