function [beta_save,omega_save] = LogisticPolyaGamma(Y,X,thinin,burnin,nmc)
  [n,p] = size(X);
  betas = zeros(p,1);
  omega = PGrnd(n,1,2);

  beta_save = zeros(p,nmc);
  omega_save = zeros(n,nmc);
  disp('Burning in ... ')
  for burn=1:burnin
    if mod(burn,1000) == 0
      disp(['(Burnin) ',num2str(burn),' out of ',num2str(burnin)])
    end
    % betas
    Sigmainv = X.'*diag(omega)*X+eye(p);
    Sigmainv = 0.5*(Sigmainv+Sigmainv.');
    SigmainvChol = chol(Sigmainv);
    mu = SigmainvChol\(SigmainvChol.'\(X.'*(Y-0.5*ones(n,1))));
    betas = mu+SigmainvChol\randn(p,1);

    % omegas
    for j=1:n
      omega(j) = PGrnd(1,1,X(j,:)*betas);
    end
  end
  disp('Starting main iterations ... ')
  for iter=1:nmc
    if mod(iter,1000) == 0
      disp(['iter = ',num2str(iter),' out of ',num2str(nmc)])
    end
    for thin=1:thinin
      % betas
      Sigmainv = X.'*diag(omega)*X+eye(p);
      Sigmainv = 0.5*(Sigmainv+Sigmainv.');
      SigmainvChol = chol(Sigmainv);
      mu = SigmainvChol\(SigmainvChol.'\(X.'*(Y-0.5*ones(n,1))));
      betas = mu+SigmainvChol\randn(p,1);

      % omegas
      for j=1:n
        omega(j) = PGrnd(1,1,X(j,:)*betas);
      end
    end
    beta_save(:,iter) = betas;
    omega_save(:,iter) = omega;
  end
end