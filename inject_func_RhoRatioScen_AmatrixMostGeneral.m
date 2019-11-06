function [C,Ceq,GradC,GradCeq] = inject_func_RhoRatioScen_AmatrixMostGeneral(a_array)

  C = [];     % there are no nonlinear-inequality constraints for a_array.
  GradC = [];
  GradCeq = [];
  
  D = (size(a_array,2)+2)/3;

  %%% constructing the identity matrix on the D-dim Hilbert space
  I_D = sparse(1:D,1:D,ones(1,D),D,D);

  %%% constructing the identity matrix on the D-dim Hilbert space (this is 
  %%% the left eigen-marix of the T-matrix in the MPS language; i.e. a row 
  %%% vector, which is also the left eigenvector in the flattened notation):
  L1 = I_D;
  L1_flat = reshape(L1,[1,D^2]);	% the reshaped or flattened version to produce the vector (1|.
  L1_flat0 = L1_flat;

  [Am,A0,Ap] = Amatrices_RhoRatioScen_AmatrixMostGeneral(D,a_array);
  T = kron(conj(Am),Am)+kron(conj(A0),A0)+kron(conj(Ap),Ap);  
  
%   T_r = T;
%   for ii = 1:D
%       index = 2*ii-1;
%       T_r(ii:index-1,:) = [];
%       T_r(:,ii:index-1) = [];
%       T_r(index+1:D+ii-1,:) = [];
%       T_r(:,index+1:D+ii-1) = [];
%   end
  
  %%% therefore we use a LU-decomposition method:
%  [SpLeft, ~] = spspaces(I_D-T_r,1,tol);
%  null_size_cond = size(SpLeft{3},1);

  %%% constructing the rho_ss on the D-dim Hilbert space (this is the right
  %%% eigenmarix of the T-matrix in the MPS language; i.e. a column vector, 
  %%% which is also the right eigenvector in the flattened notation):
  R1 = eye(D);
  %R1(1,1) = 1.0;
  for ii = 2:D
     R1(ii,ii) = R1(ii-1,ii-1) * ( Am(ii,ii-1)/Ap(ii-1,ii) )^2;   
  end
  %disp(R1); %DEBUG
  R1 = sparse(R1);
  R1_flat = reshape(R1,[D^2,1]);	% the reshaped or flattened version to produce the vector |1).
  N0 = L1_flat*R1_flat;
  R1 = R1/N0;
  R1_flat0 = reshape(R1,[D^2,1]);

%  FixedPoint_cond = trace(Am*R1*Am' + A0*R1*A0' + Ap*R1*Ap' - R1);  
%  FixedPoint_cond = sum(T*R1_flat0 - R1_flat0);
%  ortho_cond = trace(Am'*Am + A0'*A0 + Ap'*Ap)/D; 
%  ortho_cond = sum(L1_flat0*T - L1_flat0);
%  Sz_null_cond = trace(Am'*Am - Ap'*Ap);
  Sz_null_cond = L1_flat0*(-kron(conj(Am),Am) + kron(conj(Ap),Ap))*R1_flat0;
  IdentityExpectation = L1_flat0 * T * R1_flat0;
  
  %Ceq = (null_size_cond-1)^2 + (ortho_cond-1)^2 + Sz_null_cond^2;
  %Ceq = (ortho_cond-1)^2 + Sz_null_cond^2 + FixedPoint_cond^2;
  %Ceq = (ortho_cond-1)^2 + Sz_null_cond^2;
  Ceq = abs(IdentityExpectation-1) + abs(Sz_null_cond);
  %Ceq = abs(IdentityExpectation-1);

end