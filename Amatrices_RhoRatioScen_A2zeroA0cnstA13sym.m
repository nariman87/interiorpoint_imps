function [A0,A1,A2,A3] = Amatrices_RhoRatioScen_A2zeroA0cnstA13sym(D,a_array)

  %%% initializing A-matrices: 
  %%% Obviousy, it is more efficient to work with A-matrices as row-vectors 
  %%% due to their sparse form. But the following method will be more 
  %%% effective for the structure we want to work with.
  A0 = zeros(D);
  A1 = zeros(D);
  A2 = zeros(D);
  A3 = zeros(D);

  for ii = 1:D 

    if ( ii == 1)
      
      A0(ii+1,ii) = sqrt(a_array(1));      
      A3(ii,ii+1) = sqrt(a_array(ii+1));

      A1(ii,ii) = sqrt( 1 - A0(ii+1,ii)^2 );
      
    elseif (ii == D)
         
      A1(ii,ii) = sqrt( 1 - A3(ii-1,ii)^2 );
      
    else
          
      A0(ii+1,ii) = sqrt(a_array(1));
      
      %%% filling ALL the remaining elements of A3:  
      if (ii > floor(D/2)) && (ii<D-1)  
        A3(ii,ii+1) = sqrt(a_array(1)) * sqrt(a_array(D-ii+1)/a_array(D-ii));
      elseif ii == D-1
        A3(ii,ii+1) = a_array(1)/sqrt(a_array(D-ii+1));  
      elseif ii == floor(D/2)
        A3(ii,ii+1) = sqrt(a_array(1)) * sqrt(a_array(ii)/a_array(D-ii));
      else
        A3(ii,ii+1) = sqrt(a_array(1)) * sqrt(a_array(ii)/a_array(ii+1));
      end
        
      A1(ii,ii) = sqrt( 1 - A0(ii+1,ii)^2 - A3(ii-1,ii)^2 );
        
    end

  end
  
  A0 = sparse(A0);
  A1 = sparse(A1);
  A2 = sparse(A2);
  A3 = sparse(A3);

end