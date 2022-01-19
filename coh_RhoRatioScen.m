function coh = coh_RhoRatioScen_AmatrixMostGeneral(D,x)

    %%% constructing the identity matrix on the D-dim Hilbert space (this is the left
    %%% eigenmarix of the T-matrix in the MPS language; i.e. a row vector, which is also the
    %%% left eigenvector in the flattened notation):
    L1 = sparse(1:D,1:D,ones(1,D),D,D);
    L1_flat = reshape(L1,[1,D^2]);	% the reshaped or flattened version to produce the vector (1|.
  
    [Am,A0,Ap] = Amatrices_RhoRatioScen_AmatrixMostGeneral(D,x);
    
    [T,R1,L_SigmaPlus,L_SigmaMinus,R_SigmaPlus,R_SigmaMinus] = T_func(D,Am,A0,Ap);
    
    R1_flat = reshape(R1,[D^2,1]);
    N0 = L1_flat*R1_flat;
    %N0=1;
    L1_flat0 = (1/sqrt(N0))*L1_flat;
    R1_flat0 = (1/sqrt(N0))*R1_flat;

 
    coh = cal_coh_MainFunc(D,T,L1_flat0,R1_flat0,L_SigmaPlus,L_SigmaMinus,R_SigmaPlus,R_SigmaMinus);

end
