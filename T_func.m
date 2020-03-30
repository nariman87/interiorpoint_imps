function [T,R1,L_SigmaPlus,L_SigmaMinus,R_SigmaPlus,R_SigmaMinus] = T_func(D,Am,A0,Ap)

  %%% constructing the identity matrix on the D-dim Hilbert space
  I_D = sparse(1:D,1:D,ones(1,D),D,D);

  %%% constructing the identity matrix on the D-dim Hilbert space (this is 
  %%% the left eigen-marix of the T-matrix in the MPS language; i.e. a row 
  %%% vector, which is also the left eigenvector in the flattened notation):
  L1 = I_D;
  L1_flat = reshape(L1,[1,D^2]);	% the reshaped or flattened version to produce the vector (1|.
  
  
  %%% constructing the rho_ss on the D-dim Hilbert space (this is the right
  %%% eigenmarix of the T-matrix in the MPS language; i.e. a column vector, 
  %%% which is also the right eigenvector in the flattened notation):
  R1 = eye(D);
  %R1(1,1) = 1.0;
  for ii = 2:D
     R1(ii,ii) = R1(ii-1,ii-1) * ( Am(ii,ii-1)/Ap(ii-1,ii) )^2;   
  end
  %disp(R1); %DEBUG
  N0 = trace(R1);
  R1 = R1/N0;
  R1 = sparse(R1);
  R1_flat = reshape(R1,[D^2,1]);	% the reshaped or flattened version to produce the vector |1).


  T = kron(conj(Am),Am)+kron(conj(A0),A0)+kron(conj(Ap),Ap);   % the D^2xD^2-size transfer matrix, which has a block
                                                               % tridiagonal form and is reducible.            
  %N0 = 1;  
  N0 = L1_flat*R1_flat;
  L1_flat0 = (1/sqrt(N0))*L1_flat;
  R1_flat0 = (1/sqrt(N0))*R1_flat;

  %%% Calculating the terms that were previously appearing in the a2-coefficient 
  %%% of the coherence, and now are the last elements envloping the expectation 
  %%% value (If faster speed is desired, it should be more efficient to calculate 
  %%% these inside the loop of the A-matrices function).  
  L_SigmaPlus =  sqrt(2) * L1_flat0 * ( kron(conj(A0),Am) + kron(conj(Ap),A0) );
  L_SigmaMinus =  sqrt(2) * L1_flat0 * ( kron(conj(Am),A0) + kron(conj(A0),Ap) );
  R_SigmaPlus = ( kron(conj(A0),Am) + kron(conj(Ap),A0) ) * R1_flat0 * sqrt(2);
  R_SigmaMinus = ( kron(conj(Am),A0) + kron(conj(A0),Ap) ) * R1_flat0 * sqrt(2); 

end