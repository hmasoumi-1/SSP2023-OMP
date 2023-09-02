close all
clear all
clc
rng(17)
% rng(23592) % d_max/d_min ~ 2.81
% rng(14595849) % d_max/d_min ~ 2.4251

% rng(13680112) % d_max/d_min ~ 2.8303


N = 32;
K = 2;
M = 20;

UnitDFTmtx = dftmtx(N)/sqrt(N);
P_col_nrmlzd = ((myZC(N)).'); % 1 x N and its l2_norm is equal to N
x1 = load('P_ncol_nrmlzd243');
 P_ncol_nrmlzd1 = x1.P_ncol_nrmlzd243;
x2 = load('P_ncol_nrmlzd283');
P_ncol_nrmlzd2 = x2.P_ncol_nrmlzd283;
d_minZC  = min(abs(P_col_nrmlzd*UnitDFTmtx));
d_maxZC  = max(abs(P_col_nrmlzd*UnitDFTmtx));
d_minRND1  = min(abs(P_ncol_nrmlzd1*UnitDFTmtx));
d_maxRND1  = max(abs(P_ncol_nrmlzd1*UnitDFTmtx));
d_minRND2  = min(abs(P_ncol_nrmlzd2*UnitDFTmtx));
d_maxRND2  = max(abs(P_ncol_nrmlzd2*UnitDFTmtx));

% eta = [0.01,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65, 0.7]; % Use this for Probability plot
eta = [6 8 10,12 14 16 18 20 22 24];      % Use this for MSE plot
counter1 = 0;
counter2 = 0;
itr_mntCarlo = 1;  
alpha = 1e-4;
while itr_mntCarlo <= 1000
    CircShifts = randperm(N,M);
    A_rnd1 = csmat(P_ncol_nrmlzd1,UnitDFTmtx,CircShifts);
    d_min_rnd1 = min(sqrt(sum((abs(A_rnd1)).^2)));
    d_max_rnd1 = max(sqrt(sum((abs(A_rnd1)).^2)));
    mu_rnd1 = mutual_coherence(A_rnd1);
    flag1 = (d_min_rnd1 - (K-1)*mu_rnd1*d_max_rnd1);
    if flag1<=0
        counter1 = counter1 + 1;
    end
    A_rnd2 = csmat(P_ncol_nrmlzd2,UnitDFTmtx,CircShifts);
    d_min_rnd2 = min(sqrt(sum((abs(A_rnd2)).^2)));
    d_max_rnd2 = max(sqrt(sum((abs(A_rnd2)).^2)));
    mu_rnd2 = mutual_coherence(A_rnd2);
    flag2 = (d_min_rnd2 - (K-1)*mu_rnd2*d_max_rnd2);
    if flag2<=0
        counter2 = counter2 + 1;
    end
    if ((flag1>0) && (flag2>0))
        A_zc = csmat(P_col_nrmlzd,UnitDFTmtx,CircShifts); % columns of A has equal l_2-norm
        mu_zc = mutual_coherence(A_zc);
        d_min_zc = min(sqrt(sum((abs(A_zc)).^2)));
        d_max_zc = max(sqrt(sum((abs(A_zc)).^2)));
        rho = zeros(length(eta),1);
        supp_detect_indicator1 = zeros(length(eta),1);
        NMSE1_1 = zeros(length(eta),1);
        UpperBound1_1 = zeros(length(eta),1);
        supp_detect_indicator2 = zeros(length(eta),1);
        NMSE2_1 = zeros(length(eta),1);
        UpperBound2_1 = zeros(length(eta),1);
        supp_detect_indicator_eqnorm = zeros(length(eta),1);
        NMSE_eqnorm_1 = zeros(length(eta),1);
        UpperBound_eqnorm_1 = zeros(length(eta),1);
        
        % Sparse signal construction
        support = randsample(N,K);
        x = zeros(N,1);
        x(support) = randn(K,1) + 1j*randn(K,1);
        x=x./norm(x);
        x_abs_min = min(abs(nonzeros(x)));
            
            % noise
        v_tmp = (randn(M,1) + 1j*randn(M,1))/sqrt(2);
            
        for itr_eta = 1:length(eta)            
            
            % noise
            sigma = x_abs_min/((eta(itr_eta))*sqrt(2*log(N)));
            v = sigma * v_tmp;
            
            % measurements
            y1 = A_rnd1*x + v;
             y2 = A_rnd2*x + v;
            y_eqnorm = A_zc*x + v;
            % OMP reconstruction
            [supphat1,xhat1] = omp(y1,A_rnd1,K);
            [supphat2,xhat2] = omp(y2,A_rnd2,K);
            [supphat_eqnorm,xhat_eqnorm] = omp(y_eqnorm,A_zc,K);
            
            % Examine supprt recovery
            supp = sort(support,'descend');
            if numel(supphat1) == K
                ss1 = sum(abs(supphat1-supp));
                if ss1 < 1 % exact support recovery
                    supp_detect_indicator1(itr_eta,1) = 1;
                end
            end
            
            if numel(supphat2) == K
                ss2 = sum(abs(supphat2-supp));
                if ss2 < 1 % exact support recovery
                    supp_detect_indicator2(itr_eta,1) = 1;
                end
            end
            
            if numel(supphat_eqnorm) == K
                ss_eq = sum(abs(supphat_eqnorm-supp));
                if ss_eq < 1 % exact support recovery
                    supp_detect_indicator_eqnorm(itr_eta,1) = 1;
                end
            end
            
            % NMSE
            NMSE1_1(itr_eta,1) = ((norm(x-xhat1))^2);
            NMSE2_1(itr_eta,1) = ((norm(x-xhat2))^2);
            NMSE_eqnorm_1(itr_eta,1) = ((norm(x-xhat_eqnorm))^2);
            
            rho(itr_eta,1) = sigma*sqrt(2*(1+alpha)*log(N));
            UpperBound1_1(itr_eta,1) = ((d_max_rnd1/d_min_rnd1)^2) * ((K*((rho(itr_eta,1))^2)) / ((d_min_rnd1-(K-1)*mu_rnd1*d_max_rnd1)^2));
            UpperBound2_1(itr_eta,1) = ((d_max_rnd2/d_min_rnd2)^2) * ((K*((rho(itr_eta,1))^2)) / ((d_min_rnd2-(K-1)*mu_rnd2*d_max_rnd2)^2));
            UpperBound_eqnorm_1(itr_eta,1) = ((d_max_zc/d_min_zc)^2) * ((K*((rho(itr_eta,1))^2)) / ((d_min_zc-(K-1)*mu_zc*d_max_zc)^2));
            dkfjdkj = 0;
        end
 
        SUPP1(itr_mntCarlo,:) = supp_detect_indicator1.';
        NMSE1_2(itr_mntCarlo,:) = NMSE1_1.';
        UpperBound1_2(itr_mntCarlo,:) = UpperBound1_1.';
        
        SUPP2(itr_mntCarlo,:) = supp_detect_indicator2.';
        NMSE2_2(itr_mntCarlo,:) = NMSE2_1.';
        UpperBound2_2(itr_mntCarlo,:) = UpperBound2_1.';
        
        SUPP_eqnorm(itr_mntCarlo,:) = supp_detect_indicator_eqnorm.';
        NMSE_eqnorm_2(itr_mntCarlo,:) = NMSE_eqnorm_1.';
        UpperBound_eqnorm_2(itr_mntCarlo,:) = UpperBound_eqnorm_1.';
        
        itr_mntCarlo = itr_mntCarlo + 1;
        
    else
        
    end
        
end
Prob_exact_supp1 = (mean(SUPP1)).';
NMSE1 = (mean(NMSE1_2)).';
UpperBound1 = (mean(UpperBound1_2)).';
Prob_exact_supp2 = (mean(SUPP2)).';
NMSE2 = (mean(NMSE2_2)).';
UpperBound2 = (mean(UpperBound2_2)).';
Prob_exact_supp_eqnorm = (mean(SUPP_eqnorm)).';
NMSE_eqnorm = (mean(NMSE_eqnorm_2)).';
UpperBound_eqnorm = (mean(UpperBound_eqnorm_2)).';
clear SUPP1 NMSE1_2 SUPP2 NMSE2_2 SUPP_eqnorm NMSE_eqnorm_2

arnd1 = sqrt(sum((abs(A_rnd1)).^2));
d_maxmin1 = max(arnd1)/min(arnd1);
arnd2 = sqrt(sum((abs(A_rnd2)).^2));
d_maxmin2 = max(arnd2)/min(arnd2);

% % Probability of exact support recovery
% prob_plot(eta,Prob_exact_supp_eqnorm,Prob_exact_supp1,Prob_exact_supp2)  % Comment the nmse_plot line and its "eta" to run this line






% Mean square error 
nmse_plot(eta,NMSE_eqnorm,NMSE1,NMSE2,UpperBound_eqnorm,UpperBound1,UpperBound2) % Comment the prob_plot line and its "eta" to run this line

