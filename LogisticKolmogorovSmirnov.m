function [beta_save,ystar_save,nu_save] = LogisticKolmogorovSmirnov(Y,X,thinin,burnin,nmc)
  [n,p] = size(X);
  XtX = X.'*X;
  betas = zeros(p,1);
  nu = 1;
  ystar = zeros(n,1);
  beta_save = zeros(p,nmc);
  ystar_save = zeros(n,nmc);
  nu_save = zeros(1,nmc);
  disp('Burning in ... ')
  for burn=1:burnin
    if mod(burn,1000) == 0
      disp(['(Burnin) ',num2str(burn),' out of ',num2str(burnin)])
    end
    % ystar
    Xbeta = X*betas;
    for i=1:n
      if Y(i) == 1
        ystar(i) = tnormrnd(Xbeta(i),2*nu,0,Inf);
      else
        ystar(i) = tnormrnd(Xbeta(i),2*nu,-Inf,0);
      end
    end

    % betas
    Sigmainv = 0.25*XtX/(nu*nu)+eye(p);
    Sigmainv = 0.5*(Sigmainv+Sigmainv.');
    SigmainvChol = chol(Sigmainv);
    mu = SigmainvChol\(SigmainvChol.'\((X.'*ystar)*0.25/(nu*nu)));
    betas = mu+SigmainvChol\randn(p,1);

    % nu
    nustar = KSrnd(1);
    LR = nu/nustar*exp(-0.125*(1/(nustar*nustar)-1/(nu*nu))*sum((ystar-X*betas).^2));
    if rand() < LR
      nu = nustar;
    end
  end
  disp('Starting main iterations ... ')
  for iter=1:nmc
    if mod(iter,1000) == 0
      disp(['iter = ',num2str(iter),' out of ',num2str(nmc)])
    end
    for thin=1:thinin
      % ystar
      Xbeta = X*betas;
      for i=1:n
        if Y(i) == 1
          ystar(i) = tnormrnd(Xbeta(i),2*nu,0,Inf);
        else
          ystar(i) = tnormrnd(Xbeta(i),2*nu,-Inf,0);
        end
      end

      % betas
      Sigmainv = 0.25*XtX/(nu*nu)+eye(p);
      Sigmainv = 0.5*(Sigmainv+Sigmainv.');
      SigmainvChol = chol(Sigmainv);
      mu = SigmainvChol\(SigmainvChol.'\((X.'*ystar)*0.25/(nu*nu)));
      betas = mu+SigmainvChol\randn(p,1);

      % nu
      nustar = KSrnd(1);
      LR = nu/nustar*exp(-0.125*(1/(nustar*nustar)-1/(nu*nu))*sum((ystar-X*betas).^2));
      if rand() < LR
        nu = nustar;
      end
    end
    beta_save(:,iter) = betas;
    ystar_save(:,iter) = ystar;
    nu_save(iter) = nu;
  end
end