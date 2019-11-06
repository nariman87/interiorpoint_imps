%%% File name: MakePoissonianUltLaser_MPS_ExactSols_RhoAnsatz_FixedA0.m (a Matlab function) | Author: Seyed Nariman Saadatmand | Contact:  n.saadatmand@griffith.edu.au | Created in:  10/May/2018

%%% Notes: 
% - This function was written in an all-contained form and avoids calling
%   external functions as much as possible.


function [coh] = PumpPoissUltLaser_iMPS_ExactNumerics_RhoPoissAnsatz_FixedA0(p, D_min, DeltaD, D_max, calc_LambdaPrimes, search4corrs, calc_corrsVtime, calc_final_g4)

    %%% Output details:
    % [...]		...
  

    %%% printing the welcome message to STDOUT:
    VER="clean-9.26.01";   %(only for this tool)
    fprintf('\nMain Matlab function | ''ultimate-coherence laser'' project | A^{0}-matrix constant and Poissonian-distribution for rho^{ss} | version: %s\n', VER);
    fprintf('Copyright (C) Seyed Nariman Saadatmand 2018\n');
    fprintf('Contact: n.saadatmand@griffith.edu.au\n');
    fprintf('EXAMPLE USAGE: run ./PumpPoissUltLaser_iMPS_ExactNumerics_RhoPoissAnsatz_FixedA0(4,50,1,50,''yes'',''yes'',''yes'',''yes'') in a Matlab command-line environment.\n');
    fprintf('OPTIONS DESCRIPTION: ...\n');
    fprintf('PREREQUISITES: this program requires the presence of functions ... in CWD\n\n');
    
    %%% some initial values/settings for global usage:
    %DIR='.';
    DIR=strcat(getenv('HOME'),'/Dropbox/AcademiaJobs-eDesktop/MyPapers/FourtimeCoherence-FirstPublication/MatlabCollection');
    setenv('EDITOR','vim');
    %beta = 0.01; % when necessary, this is our choice for the fixed arbirtrary flux parameter (which also specifies time discretization).
    tol = 1e-14;
    SmallC = 1e-18;
        
    
    %%% setting up the 'non-reccuring' output files:
    filename_coh = strcat('Coherence_D-p',num2str(p),'-iMPS_ExactNumerics_RhoPoissAnsatz_FixedA0-PumpPoissUltLaser.out');
    if exist(fullfile(DIR,filename_coh),'file')
     fprintf('NOTE: file %s already exist; new data will be attached to its end ...\n', fullfile(DIR,filename_coh));
     FileID_final = fopen( fullfile(DIR,filename_coh) , 'at');
     if FileID_final==-1
       error('ERROR: cannot open the following file for writing: %s', fullfile(DIR,filename_coh));
     end
    else
     edit(fullfile(DIR,filename_coh));
     FileID_final = fopen( fullfile(DIR,filename_coh) , 'at');
     if FileID_final==-1
       error('ERROR: cannot open the following file for writing: %s', fullfile(DIR,filename_coh));
     end
     fprintf(FileID_final,'#D\t#L1_flat0*T*R1_flat0\t#coh_max\t#lambda2_prime\t#lambda3_prime\t#lambda4_prime\t#1-g2(tau)/N^2\t#LastError\t#1-<b''(0)b''(eps*tau)b((2-eps)*tau)b(2*tau)>/(N^2*exp(-(4-eps)*tau*ell))\t#LastError\t#N*t1_max\t#N*t2_max\t#N*t3_max\t#FourTimeCorr_abs\n');  
    end
    

    for D = D_min : DeltaD : D_max 
        
        %%% constructing the identity matrix on the D-dim and D^2-dim Hilbert spaces:
        I_D = sparse(1:D,1:D,ones(1,D),D,D);
        I_D2 = sparse(1:D^2,1:D^2,ones(1,D^2),D^2,D^2);
    
        %%% constructing the identity matrix on the D-dim Hilbert space (this is the left
        %%% eigenmarix of the T-matrix in the MPS language; i.e. a row vector, which is also the
        %%% left eigenvector in the flattened notation):
        L1 = I_D;
        L1_flat = reshape(L1,[1,D^2]);	% the reshaped or flattened version to produce the vector (1|.
        
        
        %%% constructing the rho_ss on the D-dim Hilbert space (this is the right
        %%% eigenmarix of the T-matrix in the MPS language; i.e. a column vector, which is also the
        %%% right eigenvector in the flattened notation):
        R1 = zeros(D);
        %R1(1,1) = 1.0;
        for ii = 1:D
           R1(ii,ii) = sin(pi*ii/(D+1))^p;
           %mu = (D+1)/2;
           %sigma = xi*mu;
           %R1(ii,ii) = exp(-( (ii-mu)/(sqrt(2)*sigma) )^2 ) + 1e-15;
           %mu = (p*(D+1)/2)^2;
           %k = ii - 1 + mu - (D-1)/2;
           %R1(ii,ii) = (mu/k)*R1(ii-1,ii-1);
        end
        R1 = sparse(R1);
        R1_flat = reshape(R1,[D^2,1]);	% the reshaped or flattened version to produce the vector |1).

        N0 = L1_flat*R1_flat;  
        L1_flat0 = L1_flat;
        R1_flat0 = (1/N0)*R1_flat;
       
        
        %%% finding B-matrices for different scenarios confiremd by 'most-general'
        %%% numerical calculations and physical arguments:        
        B0 = zeros(D);        
        %B0(2,1) = 1.0;
        for ii = 2:D
           B0(ii,ii-1) = 1.0;   
        end
        B0 = sparse(B0);
        
        B3 = zeros(D);        
        %B3(D-1,D) = 1.0;
        for ii = 2:D
           B3(ii-1,ii) = sqrt(R1(ii-1,ii-1)/R1(ii,ii));
        end
        B3 = sparse(B3);
        
        L0 = zeros(D);        
        L0(1,1) = 1.0;
        L0(D,D) = R1(D-1,D-1)/R1(D,D);
        for ii = 2:D-1
           L0(ii,ii) = 1 + (R1(ii-1,ii-1)/R1(ii,ii));
        end
        L0 = sparse(L0);
            
        fprintf('NOTES: calculating coh-related quantities ...\n');
        %%% calculating coh-related quantities:
        L = kron(conj(B0),B0)+kron(conj(B3),B3) - 0.5*(kron(conj(L0),I_D)+kron(I_D,L0));
        Q_proj = I_D2 - R1_flat0*L1_flat0;
        L_proj = Q_proj*L*Q_proj;
        IdentityExpectation = L1_flat0 * (I_D2+L) * R1_flat0;
        
        R_SigmaMinus = kron(I_D,B3) * R1_flat0;      % this is an \rho_ss-operator-based vector.
        L_SigmaPos = L1_flat0 * kron(conj(B3),I_D);  % this is an I-operator-based vector.
        R_SigmaUp = kron(conj(B3),B3) * R1_flat0;    % this is an \rho_ss-operator-based vector.
        L_SigmaUp = L1_flat0 * kron(conj(B3),B3);    % this is an I-operator-based vector.
        
        
        X = -(L_proj+SmallC*I_D2)\I_D2;
        %X = -inv(L_proj+SmallC*I_D2);
        coh = 2*abs( L_SigmaPos * X * R_SigmaMinus ); 
        
        
        lambda2_prime=-9999;       
        lambda3_prime=-9999;
        lambda4_prime=-9999;
        if strcmp(calc_LambdaPrimes,'yes')
          fprintf('NOTES: now calculating three smallest lambda^prime eignevalues ...\n');
          OPTS.tol = tol;
          OPTS.maxit = max(500,D);
          OPTS.p = min(30,D);
          [~,lambda_s] = eigs(L+SmallC*I_D2,4,'sm',OPTS);
          %lambda1_prime=lambda_s(1,1); % must be exactly zero, should be only calculated for testing purposes.      
          lambda2_prime=lambda_s(2,2);
          lambda3_prime=lambda_s(3,3);
          lambda4_prime=lambda_s(4,4);
        end
        
        
        Nt1 = -9999;
        Nt2 = -9999;
        Nt3 = -9999;
        FourTimeCorr_abs = -9999;
        CurlyNtau = sqrt((3.0/8)*full(coh));
        if strcmp(search4corrs,'yes')
          
          fprintf('NOTE: now finding the maximum of generic four-time correlator deviation in [0,2*tau]^3 window ...\n');  
          
          A = ones(1,3);
          b = 2*CurlyNtau;
          Lb = SmallC*ones(3,1);
          Ub = b*ones(3,1);
          d0_array = [0.01; 100; 0.01]; % a well-educated, yet acceptable, initial guess.
          %d0_array = Ub/3;
          options = optimoptions('fmincon','Algorithm','interior-point','StepTolerance',SmallC,'OptimalityTolerance',tol,'FunctionTolerance',tol,'ConstraintTolerance',tol,'Display','iter','FiniteDifferenceType','forward','FunValCheck','on','MaxFunctionEvaluations',5000,'MaxIterations',500,'UseParallel',true);
          
          F = @(x) -FourTimeCorr_generic(SmallC,tol,coh,D,R_SigmaMinus,L_SigmaPos,R_SigmaUp,L_SigmaUp,B3,L,x); % generally speaking, six separate optimisations are necessary           
          %F = @(x) -FourTimeCorr_generic(SmallC,tol,coh,D,R_SigmaMinus,L_SigmaPos,R_SigmaUp,L_SigmaUp,B3,L,x);
          %F = @(x) -FourTimeCorr_generic(SmallC,tol,coh,D,R_SigmaMinus,L_SigmaPos,R_SigmaUp,L_SigmaUp,B3,L,x);
          %F = @(x) -FourTimeCorr_generic(SmallC,tol,coh,D,R_SigmaMinus,L_SigmaPos,R_SigmaUp,L_SigmaUp,B3,L,x);
          %F = @(x) -FourTimeCorr_generic(SmallC,tol,coh,D,R_SigmaMinus,L_SigmaPos,R_SigmaUp,L_SigmaUp,B3,L,x);
          %F = @(x) -FourTimeCorr_generic(SmallC,tol,coh,D,R_SigmaMinus,L_SigmaPos,R_SigmaUp,L_SigmaUp,B3,L,x);
          
          [d1_array,FourTimeCorr_abs] = fmincon(F,d0_array,[],[],A,b,Lb,Ub,[],options);
          %[d1_array,FourTimeCorr_abs] = fmincon(F,d0_array,A,b,[],[],Lb,Ub,[],options);

          %d1_array = d0_array; 
          %FourTimeCorr_abs = FourTimeCorr_generic(SmallC,tol,coh,D,R_SigmaMinus,L_SigmaPos,R_SigmaUp,L_SigmaUp,B3,L,d1_array);
          
          Nt1 = d1_array(1);
          Nt2 = d1_array(2);
          Nt3 = d1_array(3);
          
        end    
        
        
        %%% if requested, calculating the desired correlators vs distance:
        if strcmp(calc_corrsVtime,'yes')
         
         %%% setting up the file for printing the CorrUp vs the distance:  
         filename_corrs = strcat('corrs_vs_tau-p',num2str(p),'-D',num2str(D),'-iMPS_ExactNumerics-PumpPoissUltLaser.out');          
         if exist(fullfile(DIR,filename_corrs),'file')
          %fprintf('NOTE: file %s already exist; new data will be attached to its end ...\n', fullfile(DIR,filename_corrs)); 
          fprintf('NOTE: file %s already exist; if any current data are available, will be used for calculations ...\n', fullfile(DIR,filename_corrs));
          FileID_corrs = fopen( fullfile(DIR,filename_corrs) , 'at');
          if FileID_corrs==-1
            error('ERROR: cannot open the following file for writing: %s', fullfile(DIR,filename_corrs));
          end
         else  
          edit(fullfile(DIR,filename_corrs));
          FileID_corrs = fopen( fullfile(DIR,filename_corrs) , 'at');
          if FileID_corrs==-1
            error('ERROR: cannot open the following file for writing: %s', fullfile(DIR,filename_corrs));
          end
          fprintf(FileID_corrs,'#CurlyN*tau\t#<b''(0)b''(eps*t)b((2-eps)*t)b(2*t)>/N^2\t#<b''(0)b''(eps*t)b((2-eps)*t)b(2*t)>/(N^2*exp(-(4-eps)*t*ell))\t#<b''(0)b''(4*t/3)b(2*t/3)b(2*t)>/N^2\t#<b''(0)b''(4*t/3)b(2*t/3)b(2*t)>/(N^2*exp(-2*t*ell/3))\t#g1_normalised\t#exp(-(4-eps)*ell*t)\t#exp(-2*ell*t/3)\n');   % printing the file header
         end   
            

         fprintf('NOTE: now calculating the correlations of interest vs distance ...\n');
          
         CurlyNtau_max = CurlyNtau;
         
         tt_max = 50; 
         for tt = 0:tt_max

          CurlyNtime = tt*((CurlyNtau_max-1)/tt_max) + 1;

          g1 = L_SigmaPos * expv( 1.0, CurlyNtime*L+SmallC*I_D2,  R_SigmaMinus, tol, D);
          %TwoPoint = L_SigmaPos * expm(CurlyNtime*L) * R_SigmaMinus;
          %TwoPoint = L_SigmaPos * padm(CurlyNtime*L) * R_SigmaMinus;
          %TwoPoint = L_SigmaPos * exp_eigs( CurlyNtime*L+SmallC*I_D2, tol ) * R_SigmaMinus;

          first_g4 = kron(I_D,B3) * expv( 1.0, SmallC*CurlyNtime*L+SmallC*I_D2, R_SigmaMinus, tol );
          first_g4 = kron(conj(B3),I_D) * expv( 1.0, 2*(1-SmallC)*CurlyNtime*L+SmallC*I_D2, first_g4, tol );
          first_g4 = L_SigmaPos * expv( 1.0, SmallC*CurlyNtime*L+SmallC*I_D2, first_g4, tol );  

          %second_g4 = kron(I_D,B3) * expv( 1.0, (2.0/3)*CurlyNtau*L+SmallC*I_D2, R_SigmaMinus, tol );
          %second_g4 = kron(conj(B3),I_D) * expv( 1.0, (2.0/3)*CurlyNtau*L+SmallC*I_D2, second_g4, tol );
          %second_g4 = L_SigmaPos * expv( 1.0, (2.0/3)*CurlyNtau*L+SmallC*I_D2, second_g4, tol );  
          second_g4 = kron(conj(B3),I_D) * expv( 1.0, (2.0/3)*CurlyNtau*L+SmallC*I_D2, R_SigmaMinus, tol );
          second_g4 = kron(I_D,B3) * expv( 1.0, (2.0/3)*CurlyNtau*L+SmallC*I_D2, second_g4, tol );
          second_g4 = L_SigmaPos * expv( 1.0, (2.0/3)*CurlyNtau*L+SmallC*I_D2, second_g4, tol );

          fprintf('%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\n', CurlyNtime, full(first_g4), full(first_g4)/exp(-(4-SmallC)*4*CurlyNtime/full(coh)), full(second_g4), full(second_g4)/exp(-(8.0/3)*CurlyNtau/full(coh)), full(g1), exp(-(4-SmallC)*4*CurlyNtime/full(coh)), exp(-(8.0/3)*CurlyNtau/full(coh)) );
          fprintf(FileID_corrs, '%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\n', CurlyNtime, full(first_g4), full(first_g4)/exp(-(4-SmallC)*4*CurlyNtime/full(coh)), full(second_g4), full(second_g4)/exp(-(8.0/3)*CurlyNtau/full(coh)), full(g1), exp(-(4-SmallC)*4*CurlyNtime/full(coh)), exp(-(8.0/3)*CurlyNtau/full(coh)) );
         
         end

         fclose(FileID_corrs);
                  
        end
        

        g2_tau0scale = -9999;
        g2_err = -9999;
        first_g4_tau0scale = -9999;
        first_err = -9999;
        %second_g4_tau0scale = -9999;
        if strcmp(calc_final_g4,'yes')
 
          fprintf('NOTES: now calculating end-points g4-correlations ...\n');
          
          %g2_tau0scale = L_SigmaUp * expm(2*sqrt(0.3*full(coh))*L+SmallC*I_D2) * R_SigmaUp;
          [g2_w, g2_err] = expv( 1.0, 2*CurlyNtau*L+SmallC*I_D2, R_SigmaUp, tol, D );
          g2_tau0scale = L_SigmaUp * g2_w;
          %disp(L_SigmaUp*R_SigmaUp); %DEBUG
          %disp(hump); %DEBUG

          [first_w1, first_err1] = expv( 1.0, SmallC*CurlyNtau*L+SmallC*I_D2, R_SigmaMinus, tol );
          first_g4_tau0scale = kron(I_D,B3) * first_w1;
          [first_w2, first_err2] = expv( 1.0, 2*(1-SmallC)*CurlyNtau*L+SmallC*I_D2, first_g4_tau0scale, tol );
          first_g4_tau0scale = kron(conj(B3),I_D) * first_w2;
          [first_w3, first_err3] = expv( 1.0, SmallC*CurlyNtau*L+SmallC*I_D2, first_g4_tau0scale, tol ); 
          first_g4_tau0scale = L_SigmaPos * first_w3;
          
          %disp(first_err1*norm(first_w2)*norm(first_w3)); %DEBUG
          %disp(norm(first_w1)*first_err2*norm(first_w3)); %DEBUG
          %disp(norm(first_w1)*norm(first_w2)*first_err3); %DEBUG
          first_err = max( [first_err1*norm(first_w2)*norm(first_w3), norm(first_w1)*first_err2*norm(first_w3), norm(first_w1)*norm(first_w2)*first_err3] );
 
          %second_g4_tau0scale = kron(conj(B3),I_D) * expv( 1.0, SmallC*CurlyNtau*L+SmallC*I_D2, R_SigmaMinus, tol );
          %second_g4_tau0scale = kron(I_D,B3) * expv( 1.0, 2*(1-SmallC)*CurlyNtau*L+SmallC*I_D2, second_g4_tau0scale, tol );
          %second_g4_tau0scale = L_SigmaPos * expv( 1.0, SmallC*CurlyNtau*L+SmallC*I_D2, second_g4_tau0scale, tol );

          first_g4_tau0scale = first_g4_tau0scale/exp(-(4-SmallC)*sqrt(6/full(coh)));
          first_err = first_err/exp(-(4-SmallC)*sqrt(6/full(coh)));
          %second_g4_tau0scale = second_g4_tau0scale/exp(-SmallC*sqrt(6/full(coh)));
 
        end      
        
        fprintf('NOTES: now printing the final results ...\n');
        fprintf('%i\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\n', D, full(IdentityExpectation), full(coh), lambda2_prime, lambda3_prime, lambda4_prime, 1-full(g2_tau0scale), g2_err, full(first_g4_tau0scale)-1, first_err, Nt1, Nt2, Nt3, FourTimeCorr_abs);
        fprintf(FileID_final, '%i\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\n', D, full(IdentityExpectation), full(coh), lambda2_prime, lambda3_prime, lambda4_prime, 1-full(g2_tau0scale), g2_err, full(first_g4_tau0scale)-1, first_err, Nt1, Nt2, Nt3, FourTimeCorr_abs);
        
    end   

    fclose(FileID_final);    
    
end    % the main function ends here. 
