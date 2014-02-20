function [ Coefs ] = l1( D,Y,error )
addpath('L1s')
l1method = 'OMPerr';

if strcmp(l1method,'OMPerr')
    Coefs= OMPerr(D,Y,error);
else
	numOfAtoms = size(D,2);
    numOfSignals = size(Y,2);
    mydim=size(D,1);
    Coefs = zeros([numOfAtoms,numOfSignals]);
    for iii = 1:numOfSignals
        %[s, err_mse, iter_time]=greed_omp_qr(Data(:,iii),myDic,numOfAtoms);
        switch l1method
            case 'GPSR_BB' 
                [x,x_debias,objective,times,debias_start,mses,taus]= GPSR_BB(Y(:,iii),D,0.35,'Initialization',0,'Monotone',0,'StopCriterion',1,'ToleranceA',2);
            case 'greed_omp_qr'
                [x, err_mse, iter_time]=greed_omp_qr(Y(:,iii),D,numOfAtoms);
            case 'SolveOMP' %Orthogonal Matching Pursuit
                [x, k] = SolveOMP(D, Y(:,iii),  'lambda',0.5 ,'isNonnegative',1,'tolerance',1e-5);
            case 'SolvePDIPA' %Primal-Dual Interior-Point Method
                [x, pditer, timeSteps, errorSteps] = SolvePDIPA(D, Y(:,iii),  'lambda',0.5 ,'isNonnegative',1,'tolerance',1e-5,'groundtruth',0);
            case 'SolveL1LS' %Gradient Projection
                [x, pditer, timeSteps, errorSteps] = SolveL1LS(D, Y(:,iii),  'lambda',1e-3 ,'isNonnegative',0,'tolerance',1e-5,'groundtruth',0);
            case 'SolveHomotopy' %Homotopy
                [x, total_iter, timeSteps, errorSteps, epsSteps]  = SolveHomotopy(D, Y(:,iii),  'lambda',1e-3 ,'isNonnegative',0,'tolerance',1e-5,'groundtruth',0);
            case 'SolveSpaRSA' %Iterative Thresholding
                [x,iter]= SolveSpaRSA(D, Y(:,iii)); %'LAMBDA',1e-3 ,'TOLERANCE',1e-5,'GROUNDTRUTH',0
            case 'SolveTFOCS' % TFOCS  https://github.com/cvxr/TFOCS
                %[x, nIter, timeSteps, errorSteps]= SolveTFOCS(D, Y(:,iii),'groundtruth',0,'maxtime',5); %'LAMBDA',1e-3 ,'TOLERANCE',1e-5,'GROUNDTRUTH',0
            case 'SolveFISTA' %Proximal Gradient %Fast IST Algorithm
                [x,nIter, timeSteps, errorSteps] = SolveFISTA(D, Y(:,iii),'groundtruth',0);
            case 'SolveSesopPCD' %SesopPCD
                %[x, nIter, timeSteps, errorSteps] = SolveSesopPCD(D,Y(:,iii));%,'groundtruth',0
            case 'SolvePALM' %Primal Augmented Lagrange Multiplier 
                [x, nIter, timeSteps, errorSteps] = SolvePALM(D, Y(:,iii),'groundtruth',0)%
            case 'SolveDALM' %Dual Augmented Lagrange Multiplier
                [x, nIter, timeSteps, errorSteps] = SolveDALM(D, Y(:,iii),'groundtruth',0);
            case 'l1magic'
                x = l1eq_pd(zeros([numOfAtoms,1]), D, 0, Y(:,iii), 1e-3); %pdtol, pdmaxiter, cgtol, cgmaxiter)
            case 'NESTA' %Nesterov's Algorithm
            case 'YALL1' %Alternating Direction Method
            case 'FPC_AS' %Fixed-Point Continuation and Active Set: FPC_AS.m %Bregman Iterative Regularizationcase 
            case 'AMP' %Approximate Message Passing
            otherwise
                disp('otherwise');
                
        end
        Coefs(:,iii)=x';
    end
end

%http://www.eecs.berkeley.edu/~yang/software/l1benchmark/
%     l1-Homotopy
%     SpaRSA
%     YALL1:
%     NESTA
%     FISTA 
%     SpaRSA
%     TFOCS

end

