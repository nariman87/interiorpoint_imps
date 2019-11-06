% Nariman Saadatmand (all rights reserved @ 2017) (n.saadatmand@griffith.edu.au)
% UltLimitCoh_DiscPureMPS_approx_RhoRatioScen_AmatrixMostGeneral.m (a Matlab function)
% Part of ``Interior-point iMPS for infinite-range Hamiltonians'' project. 


function [coh] = UltLimitCoh_DiscPureMPS_approx_RhoRatioScen_AmatrixMostGeneral(D_min, DeltaD, D_max, coefficient0)
  

    %%% printing the welcome message to STDOUT:
    %DIR='.';
    DIR=strcat(getenv('HOME'),'/Dropbox/AcademiaJobs--eDesktop/MyPapers/InteriorPoint--InfiniteRangeHams/MatlabCollection');
    %DIR='/media/nariman/1TB-internal-HDD/Dropbox/AcademiaJobs--eDesktop/MyPapers/InteriorPoint--InfiniteRangeHams/MatlabCollection';
    VER="clean-8.21.91";   % (only for this tool)
    fprintf('\nMain Matlab function | ''ultimate quantum limit to the coherence'' project | no scenarios for any invovled MPS matrices | version: %s\n', VER);
    fprintf('Copyright (C) Seyed Nariman Saadatmand 2017\n');
    fprintf('Contact: n.saadatmand@griffith.edu.au\n');
    fprintf('EXAMPLE USAGE: run ./UltLimitCoh_DiscPureMPS_approx_RhoRatioScen_AmatrixMostGeneral(10,10,50,0.01,40) in a Matlab command-line environment.\n');
    fprintf('OPTIONS DESCRIPTION: ...\n');
    fprintf('PREREQUISITES: this program requires the presence of functions, ''cal_coh_FinalFunc'', ''T_func'', and ... in CWD\n\n');
   
    
    %%% initial values/settings for global usage:
    setenv('EDITOR','vim');
    tol = 1e-16;
    SmallC = 1e-32;
        
        
    %%% setting up the 'non-reccuring' output files:
    FilenameFinal = strcat('Coherence_D-RhoRatioScen_AmatrixMostGeneral-UltQuantumLimitProj.out');
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
     fprintf(FileID_final,'#D\t#null_space_size\t#(L1_flat0*T)*R1_flat0\t#SigmaZ1\t#coh_min\t#finding_algorithm\n');   % printing the file header  
    end
    
    
    for D = D_min : DeltaD : D_max
        
        Nbasis=3*D-2;
        
        %%% constructing the identity matrix on the D-dim Hilbert space
        I_D = sparse(1:D,1:D,ones(1,D),D,D);
        
        %%% setting up the file for printing the optimized outputs:
        filename_opt = strcat('OptimizedMatrices-D',num2str(D),'-RhoRatioScen_AmatricesMostGeneral-UltQuantumLimitProj.out');          
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
          fprintf(FileID_opt,'#matrix_index(rows_of_nonzero_diag)\t#Am\t#A0\t#Ap\t#rho_ss\n');   % printing the file header
        end   

    
        %%% constructing the identity matrix on the D-dim Hilbert space (this
        %%% is the left eigenmarix of the T-matrix in the MPS language; i.e. 
        %%% a row vector, which is also the left eigenvector in the flattened 
        %%% notation):
        L1 = sparse(1:D,1:D,ones(1,D),D,D);
        L1_flat = reshape(L1,[1,D^2]);	% the reshaped or flattened version to produce the vector (1|.
               
        %%% setting up initial values, restricting bounds, and optimization options:
        %a_array0 = [0.1, 0.1, -0.1, 0.2, 0.3];
        %a_array0 = coefficient0*ones(1,Nbasis);
        a_array0 = coefficient0*(2*rand(1,Nbasis)-ones(1,Nbasis));
        
%         [Am,A0,Ap] = Amatrices_RhoRatioScen_AmatrixMostGeneral(D,a_array0); %DEBUG
%         disp(abs(Am)); %DEBUG
%         disp(abs(A0)); %DEBUG
%         disp(abs(Ap)); %DEBUG

        %Lb = zeros(1,Nbasis);
        Lb = -1*ones(1,Nbasis);
        Ub = ones(1,Nbasis);
        
        %%% building the (4*D-3)x(4*D-4)-size A_cond and b_cond arrays;
        %%% arbitrarily, we decided to set the D necessary conditions in
        %%% the first D-rows of A_cond.
%         A_cond = zeros(Nbasis);
%         b_cond = zeros(1,Nbasis);
%         for ii = 2:D-1 
%               A_cond(ii,ii)=1;
%               %A_cond(ii,D-1+ii-1)=1;
%               b_cond(ii)=0.5;   
%         end
       
           
        opts = optimoptions(@fmincon,'Algorithm','interior-point','StepTolerance',tol,'FunctionTolerance',tol,'OptimalityTolerance',tol,'ConstraintTolerance',SmallC,'Display','iter','FiniteDifferenceType','forward','FunValCheck','on','MaxFunctionEvaluations',1e5,'MaxIterations',500,'UseParallel',true);
        %opts = optimoptions('patternsearch','Cache','on','CacheTol',tol,'ConstraintTolerance',tol,'Display','iter','FunctionTolerance',tol,'MaxFunctionEvaluations',1e5,'MaxIterations',1e5,'MeshTolerance',tol,'StepTolerance',SmallC,'UseParallel',true);

        F = @(x) coh_RhoRatioScen_AmatrixMostGeneral(SmallC,D,x);
        nonlcon = @inject_func_RhoRatioScen_AmatrixMostGeneral;
        
        [x,coh] = fmincon(F,a_array0,[],[],[],[],Lb,Ub,nonlcon,opts);        
        %[x,coh] = patternsearch(F,a_array0,A_cond,b_cond,[],[],Lb,Ub,[],opts);
        
        %problem = createOptimProblem('fmincon','objective',F,'x0',a_array0,'Aineq',A_cond,'bineq',b_cond,'lb',Lb,'ub',Ub,'options',opts);
        %gs = GlobalSearch;
        %[x,coh] = run(gs,problem);

            
        fprintf('NOTES: re-calculating numerical values of some quantities for sanity-checking purposes ...\n');
        %%% re-calculating some more quantities for sanity-checking purposes:
        [Am,A0,Ap] = Amatrices_RhoRatioScen_AmatrixMostGeneral(D,x);
        [~,Ceq,~,~] = inject_func_RhoRatioScen_AmatrixMostGeneral(x);
        [T,R1,~,~] = T_func(D,Am,A0,Ap);
        R1_flat = reshape(R1,[D^2,1]);
        %N0 = L1_flat*R1_flat;
        N0=1;
        L1_flat0 = (1/sqrt(N0))*L1_flat;
        R1_flat0 = (1/sqrt(N0))*R1_flat;
        %ortho_cond = trace(Am'*Am + A0'*A0 + Ap'*Ap - I_D);
        %ortho_cond = sum(L1_flat0*T - L1_flat0);
        %FixedPoint_cond = trace(Am*R1*Am' + A0*R1*A0' + Ap*R1*Ap' - R1);
        %FixedPoint_cond = sum(T*R1_flat0 - R1_flat0);
        %disp(T*R1_flat0 - R1_flat0); %DEBUG
        IdentityExpectation = L1_flat0 * T * R1_flat0;
        SigmaZ1 = L1_flat0 * ( -kron(conj(Am),Am)+kron(conj(Ap),Ap) ) * R1_flat0; 
        %SigmaP_ident = sqrt(2)*L1_flat0 * ( kron(conj(A0),Am)+kron(conj(Ap),A0) ) * R1_flat0;
        %SigmaM_ident = sqrt(2)*L1_flat0 * ( kron(conj(Am),A0)+kron(conj(A0),Ap) ) * R1_flat0;
        T_r = T;
        for ii = 1:D
           index = 2*ii-1;
           T_r(ii:index-1,:) = [];
           T_r(:,ii:index-1) = [];
           T_r(index+1:D+ii-1,:) = [];
           T_r(:,index+1:D+ii-1) = [];
        end
        [SpLeft,~] = spspaces(I_D-T_r,1,tol);
        null_size_cond = size(SpLeft{3},1);

%        disp(full(SigmaP_ident)); %DEBUG
%        disp(full(SigmaM_ident)); %DEBUG      
%        disp(abs(Am)); %DEBUG
%        disp(abs(A0)); %DEBUG
%        disp(abs(Ap)); %DEBUG
%        disp(R1); %DEBUG
         
        % printing the optimized A-matrices and rho_ss:
        fprintf('NOTES: now printing the optimized A-matrices and rho_ss to a file ...\n');
        for ii = 1:D
          if ii == 1  
            fprintf(FileID_opt, '%i\t%.32f\t%.32f\t%.32f\t%.32f\n', ii, 9999, full(A0(ii,ii)), full(Ap(ii,ii+1)), full(R1(ii,ii)));
          elseif ii == D
            fprintf(FileID_opt, '%i\t%.32f\t%.32f\t%.32f\t%.32f\n', ii, full(Am(ii,ii-1)), full(A0(ii,ii)), 9999, full(R1(ii,ii)));
          else
            fprintf(FileID_opt, '%i\t%.32f\t%.32f\t%.32f\t%.32f\n', ii, full(Am(ii,ii-1)), full(A0(ii,ii)), full(Ap(ii,ii+1)), full(R1(ii,ii)));
          end    
        end
            
        fclose(FileID_opt);
             
            
        fprintf('%i\t%i\t%.16f\t%.16f\t%.16f\t%.16f\t%s\n', D, null_size_cond, Ceq, full(IdentityExpectation), full(SigmaZ1), full(coh), strcat('DirectOpt_coefficient0_',num2str(coefficient0,10)));
        fprintf(FileID_final, '%i\t%i\t%.16f\t%.16f\t%.16f\t%.16f\t%s\n', D, null_size_cond, Ceq, full(IdentityExpectation), full(SigmaZ1), full(coh), strcat('DirectOpt_coefficient0_',num2str(coefficient0,10)));
              
               
    end
        

    fclose(FileID_final);    
    
end    % the main function ends here. 
