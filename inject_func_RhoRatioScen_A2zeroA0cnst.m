function [C,Ceq,GradC,GradCeq] = inject_func_RhoRatioScen_A2zeroA0cnst(D,a_array)

  C = [];     % there are no nonlinear-inequality constraints for a_array.
  GradC = [];

  %%% constructing the identity matrix on the D-dim Hilbert space
  I_D = sparse(1:D,1:D,ones(1,D),D,D);

  [A0,A1,A2,A3] = Amatrices_RhoRatioScen_A2zeroA0cnst(D,a_array);
  T = kron(conj(A0),A0)+kron(conj(A1),A1)+kron(conj(A2),A2)+kron(conj(A3),A3);  
  
  T_r = T;
  for ii = 1:D
      index = 2*ii-1;
      T_r(ii:index-1,:) = [];
      T_r(:,ii:index-1) = [];
      T_r(index+1:D+ii-1,:) = [];
      T_r(:,index+1:D+ii-1) = [];
  end
  
  %%% Matlab default null-space finder, which is not optimized for sparse matrices:
  %null_basis = null(I_D-T_r');  
  %null_size = size(null_basis,2);
  
  %%% therefore we use an LU-decomposition approach:
  [SpLeft, ~] = spspaces(I_D-T_r,1,1e-8);
  null_size = size(SpLeft{3},1);
  
  Ceq = null_size-1;
  GradCeq = [];

end