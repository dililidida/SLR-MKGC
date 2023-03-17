warning off;
clear all;
addpath('./libs');
addpath('./core');
addpath('./KerDataset');
addpath('./ClusteringMeasure');
addpath('./LibADMM/proximal_operators');
fprintf('°Ô°Ô°Ô°Ô°Ô°Ô°Ô°Ô°Ô°ÔGraph Reconstruction°Ô°Ô°Ô°Ô°Ô°Ô°Ô°Ô°Ô°Ô\n');
ds = {'bbcsport2view_Kmatrix', 'YALE_Kmatrix', 'proteinFold_Kmatrix', 'flower17_Kmatrix','caltech101_mit_Kmatrix','UCI_DIGIT_Kmatrix','mfeat_Kmatrix','CCV_Kmatrix','flower102_Kmatrix'};

for di=1:length(ds)
    results=[];
    load(ds{di});
    if(size(Y,2)~=1)
        Y = Y';
    end;
    if ~isempty(find(Y==0,1))
        Y = Y + 1;
    end
    num_views = size(KH,3);
    num_cluster = size(unique(Y),1);
    nmu_samples = length(Y);
    
    arr_rate=[0.05,0.1,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1];
    arr_alpha=[1e-3,1e-2,1e-1,1,1e1,1e2,1e3];
    arr_beta=[1e-3,1e-2,1e-1,1,1e1,1e2,1e3];
    for rate = arr_rate
        clc;fprintf('%s',ds{di});
        neighbor = round( rate * nmu_samples / num_cluster);
        for alpha = arr_alpha
            fprintf('\n');
            for beta = arr_beta
                fprintf('.')
                clear U;clear R;   clear obj;            
                KHL = V10_LocalKernelCalculation(KH,neighbor,num_cluster);
                SumS = zeros(nmu_samples);
                for p=1:num_views
                    DSp = diag(1./sqrt(sum(KHL(:,:,p))+eps));
                    LapS = speye(nmu_samples) + DSp * KHL(:,:,p) * DSp;
                    LapS = (LapS + LapS')/2;
                    [U{p}, ~] = eigs(LapS, num_cluster, 'LA');
                    SumS = SumS + KHL(:,:,p);
                    
                    R{p} = eye(num_cluster);
                end
                S = SumS/num_views;
                
                %Initialization
                gamma = ones(num_views, 1)/num_views;
                maxIter = 20;
                iter = 0;
                tic
                while iter < maxIter
                    iter = iter + 1;
                    
                    %For F
                    SumH = zeros(nmu_samples, num_cluster);
                    for p=1:num_views
                        SumH = SumH + gamma(p)*(U{p}*R{p});
                    end
                    L = SumH*SumH';
                    M = L+2*alpha*S - alpha*eye(nmu_samples);
                    [F, ~, ~]=eig1(M, num_cluster, 1);
                    
                    %For S
                    TS = F*F';
                    N = (alpha/(alpha+beta))*(TS);
                    S = project_fantope(N,num_cluster); %beta\|S\|_F^2
                    
                    %For R
                    for p=1:num_views
                        if gamma(p)>1e-4
                            A =  gamma(p)*U{p}'*F;
                            [Uh,~,Vh] = svd(A,'econ');
                            R{p} = Uh*Vh';
                        end
                    end
                    
                    
                    coef = zeros(1,num_views);
                    for p=1:num_views
                        coef(1,p) = trace(F'*U{p}* R{p});
                    end
                    gamma = coef/norm(coef,2);
                    
                    
                    obj(iter) = -trace(TS*L)+alpha*norm(S-TS,'fro')^2+beta*norm(S,'fro')^2;
                    if iter > 5
                        if (abs(obj(iter)-obj(iter-1))/obj(iter) < 10^-5)
                            break
                        end
                    end
                end
                t=toc;
                S = S-diag(diag(S));
                S = abs(S)+abs(S');
                loop = 5;
                for i=1:loop
                    result=clustering8(S,num_cluster, Y);
%                     arr_ACC(i) = result(7)*100;
%                     arr_MIhat(i) = result(4)*100;
%                     arr_Pur(i) = result(8)*100;
%                     arr_Precision(i) = result(2)*100;
                    results = [results; i, rate,neighbor, alpha, beta, t, result];
                end
            end
        end
        
    end
    
    str=[ds{di} '_ker_results' '.mat'];
    save(str,'results');
    clc;
end
