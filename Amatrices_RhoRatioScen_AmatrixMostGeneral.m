function [Am,A0,Ap] = Amatrices_RhoRatioScen_AmatrixMostGeneral(D,a_array)

  %%% initializing A-matrices: 
  %%% Obviousy, it is more efficient to work with A-matrices as row-vectors 
  %%% due to their sparse form. But the following method will be more effective 
  %%% for the structure we want to work with.
  Am = zeros(D);
  A0 = zeros(D);
  Ap = zeros(D);

  for ii = 1:D 

    if ( ii == D)
      
      A0(ii,ii) = a_array(2*D-2+ii);
      
    else
   
      Am(ii+1,ii) = a_array(ii);
      Ap(ii,ii+1) = a_array(D-1+ii);
      A0(ii,ii) = a_array(2*D-2+ii);
        
    end

  end
  
  Am = sparse(Am);
  A0 = sparse(A0);
  Ap = sparse(Ap);

end