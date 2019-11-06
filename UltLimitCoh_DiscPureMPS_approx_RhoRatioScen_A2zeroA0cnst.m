%%% File name:   UltLimitCoh_DiscPureMPS_approx_RhoRatioScen_A2zeroA0cnst.m (a Matlab function)
%%% Description: This is the main Matlab tool for the "ultimate quantum limit to the coherence" project,
%%%              with rho_ss having an optimum distribution (Heisenberg limit for the phase uncertainty). This will 
%%%              use discrete pure-state MPS to simulate a 'laser' as a time lattice and maximize an expression for 
%%%              the coherence to find the optimize A-matrices. There are
%%%              no more arbitrary scenarios used for A-matrices.  
%%% Author: Seyed Nariman Saadatmand | Contact: n.saadatmand@griffith.edu.au | Created in: 30/Nov/2017


%%% Notes: 
% - We assume that all elements on 'non-zero' diag of A-matrices are
%   strictly non-zero, bounded, and definitive. Note that the zero blocks of the T-matrix 
%   have no contribution to the coherence. For the function expansion (in
%   the polynomial basis), we drop the 'definitive' assumptions,
%   considering that all elements are just real and bounded, which makes the forms of bounds simpler.
% - The A-matrices are size DxD and the transfer operator, T-matrix, is therefore size D^2xD^2. 
% - A-matrices are satisfying the fixed-point equation for the left \rho_ss.
% - When performing an optimization to find A-matrices, we always use 'functions' to define 
%   all quantities prior to the optimization/search algorithm to avoid symbolic programming. 
%   However, due to the simple form of A-matrices, a 'direct search' algorithm
%   is always more reliable and efficient enough for small to moderate D.
% - The Heisenberg limit for the phase uncertainty (equivalent to a uniform distribution
%   for \rho_ss) and our other assumptions here simply set a relation between A0 and A3 elements.
% - For simplicity, we set A1 = A2 (note that, these matrices are diagonal and 
%   therefore only equally contribute to any expression only involving the T-matrix); 
% - Now, due to two orthogonality/fixed-point relations, only elements of
%   A0 (or A3) are not uniquely determined. We numerically approximate non-diag elements of these A-matrices, 
%   e.g. A0_(i+1,i), by expanding in a polynomial basis, i.e. creating a degree-p polynomial as 
%   \sqrt(2)a0 + \sum_{n=1 to p=Nbasis-1} a_n * ( (i-1)/(D-2) )^n. Note that, one can always suppose 
%   that A0 is a continuos real-valued bounded function of (i-1)/(D-2). For example, in the 
%   simple case of n=4, this can produce reliable results only with FIVE free 
%   parameters, which is more than enough to find a stable global maximum. 
%   Otherwise, one must increase n, but still sticking relatively 
%   small p's; no significant change in the magnitude of coh_max is expected with the 
%   further increase in n after reaching a certain p_saturated.
% - Free parameters are arranged in a single array as follows: a_array=['A0(Nbasis-1 elements)'];
%   i.e. optimizing over Nbasis-1 free parameters and Lanczos diagonalizing a DxD sparse matrix.
%   Here, Nbasis = D/2 as we now also exploiting the fact that at optimal
%   point A1=A2 are constant matrices.


function [coh] = UltLimitCoh_DiscPureMPS_approx_RhoRatioScen_A2zeroA0cnst(D_min, DeltaD, D_max, coefficient0, D_rescale, CalculateCohCorrs)

    %%% Function inputs description:
    % [...]		...

    %%% File output details:
    % [...]		...
  

    %%% printing the welcome message to STDOUT:
    DIR=strcat(getenv('HOME'),'/Dropbox/AcademiaJobs-eDesktop/MyPapers/FourtimeCoherence_FirstPaper/MatlabCollection');
    VER="clean-8.21.91";   % (only for this tool)
    fprintf('\nMain Matlab function | ''ultimate quantum limit to the coherence'' project | no scenarios for any invovled MPS matrices | version: %s\n', VER);
    fprintf('Copyright (C) Seyed Nariman Saadatmand 2017\n');
    fprintf('Contact: n.saadatmand@griffith.edu.au\n');
    fprintf('EXAMPLE USAGE: run ./UltLimitCoh_DiscPureMPS_approx_RhoRatioScen_A2zeroA0cnst(50,1,50,0.01,9999,''no'') in a Matlab command-line environment.\n');
    fprintf('OPTIONS DESCRIPTION: ...\n');
    fprintf('PREREQUISITES: this program requires the presence of functions ''cal_coh_MainFunc'', ''T_func'', and ... in CWD\n\n');
   
    
    %%% initial values/settings for global usage:
    setenv('EDITOR','vim');
        
        
    %%% setting up the 'non-reccuring' output files:
    FilenameFinal = strcat('Coherence_D-RhoRatioScen_A2zeroA0cnst-UltQuantumLimitProj.out');
    if exist(fullfile(DIR,FilenameFinal),'file')
     fprintf('NOTE: file %s already exist; new data will be attached to its end ...\n', fullfile(DIR,FilenameFinal));
     FileID_final = fopen( fullfile(DIR,FilenameFinal) , 'at');
     if FileID_final==-1
       error('ERROR: cannot open the following file for writing: %s', fullfile(DIR,FilenameFinal));
     end
    else
     %error('ERROR: this file does not exist: %s', fullfile(DIR,filename4));   
     edit(fullfile(DIR,FilenameFinal));
     FileID_final = fopen( fullfile(DIR,FilenameFinal) , 'at');
     if FileID_final==-1
       error('ERROR: cannot open the following file for writing: %s', fullfile(DIR,FilenameFinal));
     end
     fprintf(FileID_final,'#D\t#null_space_size\t#(L1_flat0*T)*R1_flat0\t#<SigmaZ_a>\t#<SigmaZ_b>\t#coh_max\t#finding_algorithm\n');   % printing the file header  
    end
    
    
    for D = D_min : DeltaD : D_max
        
        Nbasis=D;
        
        %%% constructing the identity matrix on the D-dim Hilbert space
        I_D = sparse(1:D,1:D,ones(1,D),D,D);
        
        %%% setting up the file for printing the optimized outputs:
        filename_opt = strcat('OptimizedMatrices-D',num2str(D),'-RhoRatioScen_A2zeroA0cnst-UltQuantumLimitProj.out');          
        if exist(fullfile(DIR,filename_opt),'file')
          fprintf('NOTE: file %s already exist; new data will be attached to its end ...\n', fullfile(DIR,filename_opt)); 
          FileID_opt = fopen( fullfile(DIR,filename_opt) , 'at');
          if FileID_opt==-1
            error('ERROR: cannot open the following file for writing: %s', fullfile(DIR,filename_opt));
          end
        else  
          edit(fullfile(DIR,filename_opt));
          FileID_opt = fopen( fullfile(DIR,filename_opt) , 'at');
          if FileID_opt==-1
            error('ERROR: cannot open the following file for writing: %s', fullfile(DIR,filename_opt));
          end
          fprintf(FileID_opt,'#matrix_index(rows_of_nonzero_diag)\t#A0\t#A1\t#A2\t#A3\t#rho_ss\n');   % printing the file header
        end   
        
        %%% if necessary, setting up a file for printing out the rescaling of optimized A-matrices results:
        if D == D_rescale
            filename_rescale = strcat('CoherenceEps0-AmatricesRescaling-D',num2str(D),'-RhoRatioScen_A2zeroA0cnst-UltQuantumLimitProj.out');          
            if exist(fullfile(DIR,filename_rescale),'file')
              fprintf('NOTE: file %s already exist; new data will be attached to its end ...\n', fullfile(DIR,filename_rescale)); 
              FileID_rescale = fopen( fullfile(DIR,filename_rescale) , 'at');
              if FileID_rescale==-1
                error('ERROR: cannot open the following file for writing: %s', fullfile(DIR,filename_rescale));
              end
            else  
              edit(fullfile(DIR,filename_rescale));
              FileID_rescale = fopen( fullfile(DIR,filename_rescale) , 'at');
              if FileID_rescale==-1
                error('ERROR: cannot open the following file for writing: %s', fullfile(DIR,filename_rescale));
              end
              fprintf(FileID_rescale,'#eps0\t#null_space_size\t#(L1_flat0*T)*R1_flat0\t#<SigmaZ_a>\t#<SigmaZ_b>\t#coh_max\n');   % printing the file header
            end
        end
        
        %%% setting up the file for printing the correlations vs the distance:
        if strcmp(CalculateCohCorrs,'yes')
            filename_corrs = strcat('CohCorrs_vs_r-D',num2str(D),'-RhoRatioScen_A2zeroA0cnst-UltQuantumLimitProj.out');          
            if exist(fullfile(DIR,filename_corrs),'file')
              fprintf('NOTE: file %s already exist; new data will be attached to its end ...\n', fullfile(DIR,filename_corrs)); 
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
              fprintf(FileID_corrs,'#r\t#<sigma+(r)sigma-(0)>\n');   % printing the file header
            end
        end 
        
    
        %%% constructing the identity matrix on the D-dim Hilbert space (this
        %%% is the left eigenmarix of the T-matrix in the MPS language; i.e. 
        %%% a row vector, which is also the left eigenvector in the flattened 
        %%% notation):
        L1 = sparse(1:D,1:D,ones(1,D),D,D);
        L1_flat = reshape(L1,[1,D^2]); % the reshaped or flattened version to produce the vector (1|.
               
        %%% setting up initial values, restricting bounds, and optimization options:

        %a_array0 = coefficient0*ones(1,Nbasis);
        %a_array0 = [0.05, 0.1, -0.1, 0.2, 0.3]
        a_array0 = coefficient0*ones(1,Nbasis);
        %a_array0 = coefficient0*rand(1,Nbasis);

        a_array0(1) = 2*a_array0(1);  % it is best to start from an initial vector that 
                                      % intuitively resembles all orthogonality relations.

        Lb = zeros(1,Nbasis);
        %Lb = -1*ones(1,Nbasis);
        Ub = ones(1,Nbasis);
        
        %%% building the NbasisxNbasis-size A_cond and b_cond arrays;
        %%% arbitrarily we decided to set the D necessary conditions in
        %%% the first D-rows of A_cond.
        A_cond = zeros(Nbasis);
        b_cond = zeros(1,Nbasis);

        for ii = 2:D-1 

              A_cond(ii,1)=1;
              A_cond(ii,ii)=1;
              b_cond(ii)=1;   

        end
           
        options = optimoptions('fmincon','Algorithm','interior-point','FunctionTolerance',1e-5,'ConstraintTolerance',1e-14,'Display','iter','FiniteDifferenceType','forward','FunValCheck','on','MaxFunctionEvaluations',1000000,'MaxIterations',5000,'StepTolerance',1e-8,'UseParallel',true);
        %options = optimoptions('patternsearch','Cache','on','CacheTol',1e-14,'ConstraintTolerance',1e-14,'Display','iter','FunctionTolerance',1e-14,'MaxFunctionEvaluations',1e6,'MaxIterations',5000,'MeshTolerance',1e-14,'StepTolerance',1e-14,'UseParallel',false);

        F = @(x) -coh_RhoRatioScen_A2zeroA0cnst(D,x);

        [x,Fval] = fmincon(F,a_array0,A_cond,b_cond,[],[],Lb,Ub,[],options);
        %[x,Fval] = patternsearch(F,a_array0,[],[],[],[],Lb,Ub,[],options);

        coh = -Fval;
            
        fprintf('NOTES: re-calculating numerical values of some quantities for sanity-checking purposes ...\n');
        %%% re-calculating some more quantities for sanity-checking purposes:
        [A0,A1,A2,A3] = Amatrices_RhoRatioScen_A2zeroA0cnst(D,x);
        [~,Ceq,~,~] = inject_func_RhoRatioScen_A2zeroA0cnst(D,x);
        [T,R1,~,~] = T_func(D,A0,A1,A2,A3);
        R1_flat = reshape(R1,[D^2,1]);
        N0 = L1_flat*R1_flat;  
        L1_flat0 = (1/sqrt(N0))*L1_flat;
        R1_flat0 = (1/sqrt(N0))*R1_flat;
        IdentityExpectation = ( L1_flat0 * T ) * R1_flat0;
        SigmaZ_a = L1_flat0 * ( -kron(conj(A0),A0)+kron(conj(A1),A1)-kron(conj(A2),A2)+kron(conj(A3),A3) ) * R1_flat0;   
        SigmaZ_b = L1_flat0 * ( -kron(conj(A0),A0)-kron(conj(A1),A1)+kron(conj(A2),A2)+kron(conj(A3),A3) ) * R1_flat0; 

            
        % printing the optimized A-matrices and rho_ss:
        fprintf('NOTES: now printing the optimized A-matrices and rho_ss to a file ...\n');
        for ii = 1:D
          if ii == 1  
            fprintf(FileID_opt, '%i\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n', ii, -9999, full(A1(ii,ii)), full(A2(ii,ii)), full(A3(ii,ii+1)), full(R1(ii,ii))/N0);
          elseif ii == D
            fprintf(FileID_opt, '%i\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n', ii, full(A0(ii,ii-1)), full(A1(ii,ii)), full(A2(ii,ii)), -9999, full(R1(ii,ii))/N0);
          else
            fprintf(FileID_opt, '%i\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n', ii, full(A0(ii,ii-1)), full(A1(ii,ii)), full(A2(ii,ii)), full(A3(ii,ii+1)), full(R1(ii,ii))/N0);
          end    
        end
            
        fclose(FileID_opt);
             
            
        fprintf('%i\t%i\t%.8f\t%.8f\t%.8f\t%.8f\t%s\n', D, Ceq+1, full(IdentityExpectation), full(SigmaZ_a), full(SigmaZ_b), full(coh), strcat('DirectOpt_coefficient0_',num2str(coefficient0,10)));
        fprintf(FileID_final, '%i\t%i\t%.8f\t%.8f\t%.8f\t%.8f\t%s\n', D, Ceq+1, full(IdentityExpectation), full(SigmaZ_a), full(SigmaZ_b), full(coh), strcat('DirectOpt_coefficient0_',num2str(coefficient0,10)));
           
        
        %%% Performing the rescaling of the optimized A0 and A3 matrices 
        %%% in the case that an appropriate D_rescale value is specified -- only for
        %%% plotting purposes.
        if D == D_rescale

          fprintf('NOTES: now performing the rescaling of the optimized A0 and A3 matrices and printing to the output ...\n');

          %for eps0 = 0.001 : 0.001 : 5
          %for ii = 0:211 
          
              %eps0 = 1.041^ii * 0.005;
              eps0 = sqrt(0.000089)/(0.25*(A0(D/2,(D/2)-1)+A0((D/2)+1,D/2)+A3(D/2,(D/2)+1)+A3((D/2)+1,(D/2)+2))); % a single suitable choice for eps0.  

              A0_new = eps0 * A0;
              A3_new = eps0 * A3;

              A1_new = sqrt( I_D - A0_new'*A0_new - A3_new'*A3_new );
              A2_new = A2;

              a_array_new = cat(2,diag(A0_new,-1),diag(A3_new,+1));
              a_array_new = a_array_new.^2;

              [~,Ceq_new,~,~] = inject_func_RhoRatioScen_A2zeroA0cnst(D,a_array_new);
              [T_new,R1_new,L_SigmaPlus_new,R_SigmaMinus_new] = T_func(D,A0_new,A1_new,A2_new,A3_new);
              R1_flat_new = reshape(R1_new,[D^2,1]);
              N0_new = L1_flat*R1_flat_new;  
              L1_flat0_new = (1/sqrt(N0_new))*L1_flat;
              R1_flat0_new = (1/sqrt(N0_new))*R1_flat_new;
              IdentityExpectation_new = ( L1_flat0_new * T_new ) * R1_flat0_new;
              SigmaZ_a_new = L1_flat0_new * ( -kron(conj(A0_new),A0_new)+kron(conj(A1_new),A1_new)-kron(conj(A2_new),A2_new)+kron(conj(A3_new),A3_new) ) * R1_flat0_new;   
              SigmaZ_b_new = L1_flat0_new * ( -kron(conj(A0_new),A0_new)-kron(conj(A1_new),A1_new)+kron(conj(A2_new),A2_new)+kron(conj(A3_new),A3_new) ) * R1_flat0_new; 

              coh_new = cal_coh_FinalFunc(D,T_new,L1_flat0_new,R1_flat0_new,L_SigmaPlus_new,R_SigmaMinus_new);

              fprintf(FileID_rescale, '%.8f\t%i\t%.8f\t%.8f\t%.8f\t%.8f\n', eps0, Ceq_new+1, full(IdentityExpectation_new), full(SigmaZ_a_new), full(SigmaZ_b_new), full(coh_new));

           %end
              
           fclose(FileID_rescale);
              
        end 
        
        % if requested, calculating the correlations vs distance:
        if strcmp(CalculateCohCorrs,'yes')

         fprintf('NOTES: now calculating the correlations vs distance ...\n');   

         for r = 0:1e4:1e7               
           %corr = L_SigmaPlus * T^r * R_SigmaMinus;
           corr = L_SigmaPlus_new * T_new^r * R_SigmaMinus_new;
           fprintf(FileID_corrs, '%i\t%.14f\n', r, full(corr));
         end 

         fclose(FileID_corrs);

        end
                      
    end
        
    fclose(FileID_final);    
    
end    % the main function ends here. 