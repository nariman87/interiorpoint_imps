function [Am,A0,Ap] = Amatrices_RhoRatioScen_AmatrixMostGeneral(D,a_array)

  %%% initializing A-matrices: 
  %%% Obviousy, it is more efficient to work with A-matrices as row-vectors 
  %%% due to their sparse form. But the following method will be more effective 
  %%% for the structure we want to work with.
  Am = zeros(D);
  A0 = zeros(D);
  Ap = zeros(D);

  for ii = 1:D 

    if ( ii == 1)
        
      Am(ii+1,ii) = a_array(ii);
      %A0(ii,ii) = sign(a_array(2*D-2+ii))*sqrt(1-Am(ii+1,ii)^2);
      A0(ii,ii) = sqrt(1-Am(ii+1,ii)^2);
         
    elseif ( ii == D)
      
      Ap(ii-1,ii) = a_array(D-1+ii-1);  
      %A0(ii,ii) = sign(a_array(2*D-2+ii))*sqrt(1-Ap(ii-1,ii)^2);
      A0(ii,ii) = sqrt(1-Ap(ii-1,ii)^2);
      
    else
   
      Am(ii+1,ii) = a_array(ii);
      Ap(ii-1,ii) = a_array(D-1+ii-1);
      %A0(ii,ii) = sign(a_array(2*D-2+ii))*sqrt(1-Am(ii+1,ii)^2-Ap(ii-1,ii)^2);
      A0(ii,ii) = sqrt(1-Am(ii+1,ii)^2-Ap(ii-1,ii)^2);
        
    end

  end
  
  Am = sparse(Am);
  A0 = sparse(A0);
  Ap = sparse(Ap);

end