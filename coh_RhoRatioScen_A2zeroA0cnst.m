function coh = coh_RhoRatioScen_A2zeroA0cnst(D,x)

    %%% constructing the identity matrix on the D-dim Hilbert space (this is the left
    %%% eigenmarix of the T-matrix in the MPS language; i.e. a row vector, which is also the
    %%% left eigenvector in the flattened notation):
    L1 = sparse(1:D,1:D,ones(1,D),D,D);
    L1_flat = reshape(L1,[1,D^2]); % the reshaped or flattened version to produce the vector (1|.
  
    [A0,A1,A2,A3] = Amatrices_RhoRatioScen_A2zeroA0cnst(D,x);
    
    [T,R1,L_SigmaPlus,R_SigmaMinus] = T_func(D,A0,A1,A2,A3);
    
    R1_flat = reshape(R1,[D^2,1]);
    N0 = L1_flat*R1_flat;  
    L1_flat0 = (1/sqrt(N0))*L1_flat;
    R1_flat0 = (1/sqrt(N0))*R1_flat;
   
    coh = cal_coh_MainFunc(D,T,L1_flat0,R1_flat0,L_SigmaPlus,R_SigmaMinus);

end
