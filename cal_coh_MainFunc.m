function coh = cal_coh_MainFunc(D,T,L1_flat0,R1_flat0,L_SigmaPlus,L_SigmaMinus,R_SigmaPlus,R_SigmaMinus)

  tol=1e-10;

  %%% constructing the identity matrix on the D^2-dim Hilbert space
  I_D2 = sparse(1:D^2,1:D^2,ones(1,D^2),D^2,D^2);

  Q_proj = I_D2 - R1_flat0*L1_flat0;
  T_proj = Q_proj*T*Q_proj;
  
  %X = ((I_D2-T_proj)+1e-16*I_D2)\I_D2;                                             % FOR LATER: we must use T_reduced here and     
  %X = linsolve((I_D2-T_proj)+SmallC*I_D2,I_D2);                                   % everywhere else relevant, according to the           
  %X = I_D2 * pinv(I_D2-T_proj,SmallC);                                            % desired symmetry sector.  
  [X1,~] = gmres((I_D2-T_proj)+1e-16*I_D2,R_SigmaMinus,[],tol,min(D^2,50));
  [X2,~] = gmres((I_D2-T_proj)+1e-16*I_D2,R_SigmaPlus,[],tol,min(D^2,50));

  %coh = 2*( L_SigmaPlus * X * R_SigmaMinus + L_SigmaMinus * X * R_SigmaPlus);   % Note that we set this to be always a purely real scalar and the 'trace' is hidden in this expression.
  %coh = 2*2*real( L_SigmaPlus * X * R_SigmaMinus );                             % Are there more efficient ways to calculate the inverse of a sparse matrix?
  coh = 2*real( L_SigmaPlus * X1 + L_SigmaMinus * X2);                           % The 'abs' function is because we acknowledge that the coherence can be sensibly both large
  %coh = 2*2*real( L_SigmaPlus * X1);                                            % and negative or positive.
                                                                                   
  %disp(abs(L_SigmaPlus * X1 - L_SigmaMinus * X2)); %DEBUG                                                                                 
  %disp(full(abs(L_SigmaPlus * X * R_SigmaMinus - L_SigmaMinus * X * R_SigmaPlus))); %DEBUG 
  
end
