function coh = cal_coh_MainFunc(D,T,L1_flat0,R1_flat0,L_SigmaPlus,R_SigmaMinus)

  %%% constructing the identity matrix on the D^2-dim Hilbert space
  I_D2 = sparse(1:D^2,1:D^2,ones(1,D^2),D^2,D^2);

  Q_proj = I_D2 - R1_flat0*L1_flat0;
  T_proj = Q_proj*T*Q_proj;
  
  X = ((I_D2-T_proj)+1e-15*I_D2)\I_D2;  
  %X = linsolve((I_D2-T_proj)+1e-15*I_D2,I_D2);
  %X = I_D2 * pinv(I_D2-T_proj,1e-14);

  coh = 2*abs( L_SigmaPlus * X * R_SigmaMinus );   % Note that we set this to be always a purely real scalar and the 'trace' is hidden in this expression.
                                                   % Are there more efficient ways to calculate the inverse of a sparse matrix?
                                                   % The 'abs' function is because we acknowledge that the coherence can be sensibly both large
                                                   % and negative or positive.
                                                                                                                                                                        
end
