function FourTimeDev = FourTimeCorr_generic(SmallC,tol,coh,D,R_SigmaMinus,L_SigmaPos,R_SigmaUp,L_SigmaUp,B3,L,d_array)

          I_D = sparse(1:D,1:D,ones(1,D),D,D);
          I_D2 = sparse(1:D^2,1:D^2,ones(1,D^2),D^2,D^2);

          first_g4 = kron(I_D,B3) * expv( 1.0, d_array(3)*L+SmallC*I_D2, R_SigmaMinus, tol );
          first_g4 = kron(conj(B3),I_D) * expv( 1.0, d_array(2)*L+SmallC*I_D2, first_g4, tol );
          first_g4 = L_SigmaPos * expv( 1.0, d_array(1)*L+SmallC*I_D2, first_g4, tol );
          f =  ( d_array(1) + 4*d_array(2) + d_array(3) ) / 2;
          FourTimeDev_first = full(first_g4/exp(-f*4/full(coh))) - 1;
          %FourTimeDev_first = 0;

          %second_g4 = kron(conj(B3),I_D) * expv( 1.0, d_array(3)*L+SmallC*I_D2, R_SigmaMinus, tol );
          %second_g4 = kron(I_D,B3) * expv( 1.0, d_array(2)*L+SmallC*I_D2, second_g4, tol );
          %second_g4 = L_SigmaPos * expv( 1.0, d_array(1)*L+SmallC*I_D2, second_g4, tol );
          %f = ( d_array(1) + d_array(3) ) / 2;
          %FourTimeDev_second = full(second_g4/exp(-f*abs(lambda2))) - 1;
          FourTimeDev_second = 0;

          %OnePairMiddle_g4 = kron(conj(B3),B3) * expv( 1.0, d_array(3)*L+SmallC*I_D2, R_SigmaMinus, tol );
          %OnePairMiddle_g4 = L_SigmaPos * expv( 1.0, d_array(1)*L+SmallC*I_D2, OnePairMiddle_g4, tol );
          %f = ( d_array(1) + d_array(3) ) / 2;
          %FourTimeDev_OnePairMiddle = full(OnePairMiddle_g4/exp(-f*abs(lambda2))) - 1;
          FourTimeDev_OnePairMiddle = 0;

          %OnePairLeft_g4 = kron(conj(B3),I_D) * expv( 1.0, d_array(3)*L+SmallC*I_D2, R_SigmaMinus, tol );
          %OnePairLeft_g4 = L_SigmaUp * expv( 1.0, d_array(2)*L+SmallC*I_D2, OnePairLeft_g4, tol );
          %f = d_array(3) / 2;
          %FourTimeDev_OnePairLeft = full(OnePairLeft_g4/exp(-f*abs(lambda2))) - 1;
          FourTimeDev_OnePairLeft = 0;

          %OnePairRight_g4 = kron(I_D,B3) * expv( 1.0, d_array(2)*L+SmallC*I_D2, R_SigmaUp, tol );
          %OnePairRight_g4 = L_SigmaPos * expv( 1.0, d_array(1)*L+SmallC*I_D2, OnePairRight_g4, tol );
          %f = d_array(1) / 2;
          %FourTimeDev_OnePairRight = full(OnePairRight_g4/exp(-f*abs(lambda2))) - 1;
          FourTimeDev_OnePairRight = 0;

          %TwoPair_g4 = L_SigmaUp * expv( 1.0, d_array(2)*L+SmallC*I_D2, R_SigmaUp, tol );
          %FourTimeDev_TwoPair = full(TwoPair_g4) - 1;
          FourTimeDev_TwoPair = 0;

          %FourTimeDev = max([abs(FourTimeDev_first),abs(FourTimeDev_second),abs(FourTimeDev_OnePairMiddle),abs(FourTimeDev_OnePairLeft),abs(FourTimeDev_OnePairRight),abs(FourTimeDev_TwoPair)]);   
          FourTimeDev = abs(FourTimeDev_first);
              
end
