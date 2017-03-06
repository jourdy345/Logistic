% TRANSLATION OF R CODE
%
% AUTHOR  
%   
%     Daeyoung Lim,
%     Department of Statistics,
%     Korea University
% 
% 
% 
% DESCRIPTION
% 
%     This code evaluates the CDF of normal distribution at a given point.
%     The code is translated from the R core development team's pnorm function.
%     The main computation evaluates near-minimax approximations derived
%     from those in "Rational Chebyshev approximations for the error
%     function" by W. J. Cody, Math. Comp., 1969, 631-637.  This
%     transportable program uses rational functions that theoretically
%     approximate the normal distribution function to at least 18
%     significant decimal digits.  The accuracy achieved depends on the
%     arithmetic system, the compiler, the intrinsic functions, and
%     proper selection of the machine-dependent constants.
%  
% 
%  
% REFERENCE
%  
%    Cody, W. D. (1993).
%    ALGORITHM 715: SPECFUN - A Portable FORTRAN Package of
%    Special Function Routines and Test Drivers".
%    ACM Transactions on Mathematical Software. 19, 22-32.




function [out] = pnorm(x,mu,sig,lower_tail,log_p)
  if lower_tail ~= 0 && lower_tail ~= 1
    error('lower_tail is not Boolean')
  end
  if log_p ~= 0 && log_p ~= 1
    error('log_p is not Boolean')
  end
  if isinf(x) && (mu == x)
    out = NaN;
    return
  end
  if log_p
    R_D_0 = -Inf;
    R_D_1 = 0;
  else
    R_D_0 = 0;
    R_D_1 = 1;
  end
  if lower_tail
    R_DT_0 = R_D_0;
    R_DT_1 = R_D_1;
  else
    R_DT_0 = R_D_1;
    R_DT_1 = R_D_0;
  end

  if sig <= 0
    if sig <= 0
      if x < mu
        out = R_DT_0;
        return
      else
        out = R_DT_1;
        return
      end
    end
  end
  p_ = (x-mu)/sig;
  if isinf(p_)
    if x<mu
      out = R_DT_0;
      return
    else
      out = R_DT_1;
      return
    end
  end
  x = p_;
  a = [2.2352520354606839287,161.02823106855587881,1067.6894854603709582,18154.981253343561249,0.065682337918207449113];
  b = [47.20258190468824187,976.09855173777669322,10260.932208618978205,45507.789335026729956];
  c = [0.39894151208813466764,8.8831497943883759412,93.506656132177855979,597.27027639480026226,2494.5375852903726711,6848.1904505362823326,11602.651437647350124,9842.7148383839780218,1.0765576773720192317e-8];
  d = [22.266688044328115691,235.38790178262499861,1519.377599407554805,6485.558298266760755,18615.571640885098091,34900.952721145977266,38912.003286093271411,19685.429676859990727];
  p = [0.21589853405795699,0.1274011611602473639,0.022235277870649807,0.001421619193227893466,2.9112874951168792e-5,0.02307344176494017303];
  q = [1.28426009614491121,0.468238212480865118,0.0659881378689285515,0.00378239633202758244,7.29751555083966205e-5];

  epsilon = 2.2204460492503131e-16*0.5;
  y = abs(x);
  M_SQRT_32 = 5.656854249492380195206754896838;
  M_1_SQRT_2PI = 0.398942280401432677939946059934;
  xnum = 0;
  xden = 0;
  cum = 0;
  ccum = 0;
  xsq = 0;
  if y <= 0.67448975
    if y > epsilon
      xsq = x*x;
      xnum = a(5)*xsq;
      xden = xsq;
      for i=1:3
        xnum = (xnum+a(i))*xsq;
        xden = (xden+b(i))*xsq;
      end
    else
      xnum = 0;
      xden = 0;
    end
    tmp = x*(xnum+a(4))/(xden+b(4));
    if lower_tail
      cum = 0.5+tmp;
    else
      ccum = 0.5-tmp;
    end
    if log_p
      if lower_tail
        cum = log(cum);
      else
        ccum = log(ccum);
      end
    end
  elseif y <= M_SQRT_32
    xnum = c(9)*y;
    xden = y;
    for i=1:7
      xnum = (xnum+c(i))*y;
      xden = (xden+d(i))*y;
    end
    tmp = (xnum+c(8))/(xden+d(8));
    xsq = fix(y*16)/16;
    del = (y-xsq)*(y+xsq);
    if log_p
      cum = (-xsq*xsq*0.5)+(-del*0.5)+log(tmp);
      if (lower_tail && x > 0) || (~lower_tail && x <= 0)
        ccum = log1p(-exp(-xsq*xsq*0.5)*exp(-del*0.5)*tmp);
      end
    else
      cum = exp(-xsq*xsq*0.5)*exp(-del*0.5)*tmp;
      ccum = 1-cum;
    end
    if x > 0
      tmp = cum;
      if lower_tail
        cum = ccum;
        ccum = tmp;
      end
    end
  elseif (log_p && y < 1.0e170) || (lower_tail && -37.5193 < x && x < 8.2924) || (~lower_tail && -8.2924  < x && x< 37.5193)
    xsq = 1/(x*x);
    xnum = p(6)*xsq;
    xden = xsq;
    for i=1:4
      xnum = (xnum+p(i))*xsq;
      xden = (xden+q(i))*xsq;
    end
    tmp = xsq*(xnum+p(5))/(xden+q(5));
    tmp = (M_1_SQRT_2PI-tmp)/y;
    xsq = fix(x*16)/16;
    del = (x-xsq)*(x+xsq);
    if log_p
      cum = (-xsq*xsq*0.5)+(-del*0.5)+log(tmp);
      if (lower_tail && x > 0) || (~lower_tail && x <= 0)
        ccum = log1p(-exp(-xsq*xsq*0.5)*exp(-del*0.5)*tmp);
      end
    else
      cum = exp(-xsq*xsq*0.5)*exp(-del*0.5)*tmp;
      ccum = 1-cum;
    end
    if x > 0
      tmp = cum;
      if lower_tail
        cum = ccum;
        ccum = tmp;
      end
    end
  else
    if x > 0
      cum = R_D_1;
      ccum = R_D_0;
    else
      cum = R_D_0;
      ccum = R_D_1;
    end
  end
  if log_p
    if cum > -1.0e-37
      cum = -0;
    end
    if ccum > -1.0e-37
      ccum = -0;
    end
  else
    if cum < 1.0e-37
      cum = 0;
    end
    if ccum < 1.0e-37
      ccum = 0;
    end
  end
  if lower_tail
    out = cum;
  else
    out = ccum;
  end
end