function [T,R1,L_SigmaPlus,R_SigmaMinus] = T_func(D,A0,A1,A2,A3)

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
  R1 = zeros(D);
  R1(1,1) = 1.0;
  for ii = 2:D
     R1(ii,ii) = R1(ii-1,ii-1) * ( A0(ii,ii-1)/A3(ii-1,ii) )^2;   
  end 
  R1 = sparse(R1);
  R1_flat = reshape(R1,[D^2,1]);	% the reshaped or flattened version to produce the vector |1).

  T = kron(conj(A0),A0)+kron(conj(A1),A1)+kron(conj(A2),A2)+kron(conj(A3),A3);   % the D^2xD^2-size transfer matrix, which has a block
                                                           			             % tridiagonal form and is reducible.            
  
  N0 = L1_flat*R1_flat;  
  L1_flat0 = (1/sqrt(N0))*L1_flat;
  R1_flat0 = (1/sqrt(N0))*R1_flat;

  %%% Calculating the terms that were previously appearing in the a2-coefficient 
  %%% of the coherence, and now are the last elements envloping the expectation 
  %%% value (If faster speed is desired, it should be more efficient to calculate 
  %%% these inside the loop of the A-matrices function).  
  R_SigmaMinus = ( kron(conj(A0),A2) + kron(conj(A1),A3) ) * R1_flat0;   % this is an \rho_ss-operator-based vector.
  L_SigmaPlus =  L1_flat0 * ( kron(conj(A2),A0) + kron(conj(A3),A1) );   % this is an I-operator-based vector.  

end