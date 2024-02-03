% 1D case - Compares the FM algorithm and the Spectral algorithm with the analytical
% bound, as a function of the SNR for different distributions.
tic
sigma_vec_reduced=logspace(0,3,30).';
num_unique_sigma=length(sigma_vec_reduced);
num_rep=50;
sigma_vec=repelem(sigma_vec_reduced,num_rep);
f_vec=[0,0.1,1,2];
err_squared_a_FM=zeros(length(sigma_vec),length(f_vec));
err_squared_a_Spectral_good_true=zeros(length(sigma_vec),length(f_vec));
err_squared_a_Spectral_bad_true=zeros(length(sigma_vec),length(f_vec));
err_squared_a_Spectral_good_emp=zeros(length(sigma_vec),length(f_vec));
err_squared_a_Spectral_bad_emp=zeros(length(sigma_vec),length(f_vec));
dist_from_circ_true=zeros(1,length(f_vec));
bound=zeros(1,length(f_vec));
B=10;
N=1e6;

%a=exp(2*pi*1i*rand(B,1));%.*exp(-(1/2)*sqrt((1:B).'));%+1i*[0;ones(B-2,1)]; % k=0 to B
a=[0.707671873258695 + 0.706541237153593i;0.988238403497549 + 0.152921083741304i;-0.249817202634132 + 0.968293016223941i;-0.411632223458529 + 0.911350049437968i;-0.568871246227187 - 0.822426595639955i;0.963616145197972 - 0.267289963735642i;0.919568728484166 - 0.392929196667815i;-0.965195103437267 + 0.261530901234871i;0.0597903158092206 + 0.998210958733390i;0.0872124228383627 - 0.996189737601559i];
%a=[-0.648691381546365 + 0.761051569545368i;-0.999143727285276 + 0.0413740525750913i;0.0665012586830974 + 0.997786341154039i;0.243882166531877 - 0.969804871532267i;0.136404716517738 + 0.990653195276589i;0.137628083145165 + 0.990483978027806i;0.702785303999646 - 0.711402007645554i;-0.835561661271188 - 0.549396678378867i;0.948991591066542 + 0.315301379770201i;-0.929691461089495 + 0.368339228398607i];
a_symm_1B=[conj(flip(a(1:end))); 1 ;a];
a_half_1B= [1;a];
k_symm_1B=-B:B;
k_half_1B=0:B;
k_symm_2B=-2*B:2*B;
k_half_2B=0:2*B;


weights=cell(1,B); %max number of approximations x max freq-1
for k=2:B
    weights{k}=exp(-abs(ceil(k/2)-(1:k-1).')); %lorentzian of length n centered at ceil(n/2)
    weights{k}=weights{k}./sum(weights{k}); %normalize to sum 1.
end

for f_index=1:length(f_vec)
    temp=[0.00195583374369666 + 0.00530869159003379i;0.00195583374369666 + 0.00474988194897760i;0.00195583374369666 + 0.00419107230792141i;0.00195583374369666 + 0.00363226266686522i;0.00195583374369666 + 0.00307345302580903i;0.00195583374369666 + 0.00251464338475285i;0.00195583374369666 + 0.00195583374369666i;0.00195583374369666 + 0.00139702410264047i;0.00195583374369666 + 0.000838214461584282i;0.00195583374369666 + 0.000279404820528094i;0.00195583374369666 - 0.000279404820528094i;0.00195583374369666 - 0.000838214461584282i;0.00195583374369666 - 0.00139702410264047i;0.00195583374369666 - 0.00195583374369666i;0.00195583374369666 - 0.00251464338475285i;0.00195583374369666 - 0.00307345302580903i;0.00195583374369666 - 0.00363226266686522i;0.00195583374369666 - 0.00419107230792141i;0.00195583374369666 - 0.00474988194897760i;0.00195583374369666 - 0.00530869159003379i]...
        .*exp(1i*(f_vec(f_index)*rand(2*B,1)).^0.2)*3;
    %temp=[0.00586750123108998 + 0.0159260747701014i;0.00586750123108998 + 0.0142496458469328i;0.00586750123108998 + 0.0125732169237642i;0.00586750123108998 + 0.0108967880005957i;0.00586750123108998 + 0.00922035907742709i;0.00586750123108998 + 0.00754393015425855i;0.00586750123108998 + 0.00586750123108998i;0.00586750123108998 + 0.00419107230792141i;0.00586750123108998 + 0.00251464338475285i;0.00586750123108998 + 0.000838214461584282i;0.00586750123108998 - 0.000838214461584282i;0.00586750123108998 - 0.00251464338475285i;0.00586750123108998 - 0.00419107230792141i;0.00586750123108998 - 0.00586750123108998i;0.00586750123108998 - 0.00754393015425855i;0.00586750123108998 - 0.00922035907742709i;0.00586750123108998 - 0.0108967880005957i;0.00586750123108998 - 0.0125732169237642i;0.00586750123108998 - 0.0142496458469328i;0.00586750123108998 - 0.0159260747701014i];
    %temp=[0.00195583374369666 + 0.00530869159003379i;0.00195583374369666 + 0.00474988194897760i;0.00195583374369666 + 0.00419107230792141i;0.00195583374369666 + 0.00363226266686522i;0.00195583374369666 + 0.00307345302580903i;0.00195583374369666 + 0.00251464338475285i;0.00195583374369666 + 0.00195583374369666i;0.00195583374369666 + 0.00139702410264047i;0.00195583374369666 + 0.000838214461584282i;0.00195583374369666 + 0.000279404820528094i;0.00195583374369666 - 0.000279404820528094i;0.00195583374369666 - 0.000838214461584282i;0.00195583374369666 - 0.00139702410264047i;0.00195583374369666 - 0.00195583374369666i;0.00195583374369666 - 0.00251464338475285i;0.00195583374369666 - 0.00307345302580903i;0.00195583374369666 - 0.00363226266686522i;0.00195583374369666 - 0.00419107230792141i;0.00195583374369666 - 0.00474988194897760i;0.00195583374369666 - 0.00530869159003379i]...
    %    +f(index)/(5e4*B)*(rand(2*B,1)+1i*rand(2*B,1));
    rho_hat_symm_2B_non_normalized =[conj(flip(temp)); 1/(2*pi); temp];
    rho_hat_symm_2B_temp=(1/(2*pi)) * (1/abs(rho_hat_symm_2B_non_normalized(k_symm_2B==0)))...
        * rho_hat_symm_2B_non_normalized;
    rho_func_temp = @(x) real(rho_hat_symm_2B_temp.'*exp(1i*x'*k_symm_2B).');
    
    
    t=-pi:1e-4:pi;
    min_rho_func=min(rho_func_temp(t));
    rho_hat_symm_2B =rho_hat_symm_2B_temp;
    rho_hat_symm_2B(k_symm_2B==0)=rho_hat_symm_2B(k_symm_2B==0)-min(min_rho_func,0); %make >=0
    rho_hat_symm_2B=rho_hat_symm_2B.*(1/(2*pi))*(1/rho_hat_symm_2B(k_symm_2B==0)); %renormalize
    rho_func = @(x) real(rho_hat_symm_2B.'*exp(1i*x.'*k_symm_2B).'); %the final rho_func
    
    rho_hat_symm_1B=rho_hat_symm_2B(-B<=k_symm_2B & k_symm_2B<=B);
    rho_hat_half_2B=rho_hat_symm_2B(0<=k_symm_2B);
    rho_hat_half_1B=rho_hat_symm_2B(0<=k_symm_2B & k_symm_2B<=B);
    
    
    for sigma_index=1:length(sigma_vec)
        sigma=sigma_vec(sigma_index); %noise
        
        disp("f_index="+f_index+", sigma_index="+sigma_index)
        done_percent=100.*(f_index*sigma_index)./...
            (length(sigma_vec)*length(f_vec));
        disp(done_percent+"% done");
        %% Moment Calculation
        M_1_true=2*pi*a_symm_1B.*rho_hat_symm_1B;
        T=toeplitz(conj(rho_hat_half_2B));
        M_2_true=2*pi*diag(a_symm_1B)*T*diag(conj(a_symm_1B))+...
            (sigma^2/(2*pi))*(k_symm_1B==k_symm_1B.'|k_symm_1B==-k_symm_1B.');
        M_2_true_wave = 2*pi*diag(a_symm_1B)*T*diag(conj(a_symm_1B));
        rotations=randsample(t,N,true,abs(rho_func(t))).';
        %absolute value to avoid numeric negatives, rho should always be positive.
        a_rot=a_half_1B.*exp(-1i.*k_half_1B.'.*rotations.')... %KEEP -!!!!! though without also works for a
            +(sigma/sqrt(2*pi))*randn(B+1,N);%N colomns of rotated a vectors
        a_rot_symm=[conj(flip(a_rot(2:end,:),1)) ; a_rot];
        
        M_1_emp=mean(a_rot_symm,2);
        M_2_emp=(1/N)*(a_rot_symm*a_rot_symm');
        M_2_emp_wave = M_2_emp-(sigma^2/(2*pi))*(k_symm_1B==k_symm_1B.'|k_symm_1B==-k_symm_1B.');
        moment_type="emp";
        %moment_type="true";
        if strcmp(moment_type,"emp")==1
            M_1=M_1_emp;
            M_2=M_2_emp;
            M_2_wave=M_2_emp_wave;
        end
        if strcmp(moment_type,"true")==1
            M_1=M_1_true;
            M_2=M_2_true;
            M_2_wave=M_2_true_wave;
        end
        
        %% FM
        S=2*pi*diag(1./M_1)*M_2_wave*diag(1./M_1)';
        
        rho_approx_half_1B=zeros(length(k_half_1B),1);
        rho_approx_tilde_half_1B= zeros(length(k_half_1B),1);
        
        rho_abs_half_1B = sqrt(abs(1./(2*pi*diag(S(k_symm_1B>=0,k_symm_1B>=0)))));
        rho_approx_half_1B(k_half_1B==0)=1/(2*pi);
        %arbitrarily pick phase of 1'st coeff to be 1:
        rho_approx_half_1B(k_half_1B==1) = rho_abs_half_1B(k_half_1B==1);
        
        for k=2:B
            rho_approx_tilde_half_1B(k_half_1B==k) = sum(weights{k}.*flip(rho_approx_half_1B(1<=k_half_1B & k_half_1B <= k-1))./...
                conj(rho_approx_half_1B(1<=k_half_1B & k_half_1B <= k-1))./...
                S(k_symm_1B==k, 1<=k_symm_1B & k_symm_1B <= k-1).');
            
            rho_approx_half_1B(k_half_1B==k)=rho_abs_half_1B(k_half_1B==k).*...
                rho_approx_tilde_half_1B(k_half_1B==k)./abs(rho_approx_tilde_half_1B(k_half_1B==k));
        end
        rho_approx_symm_1B= [conj(flip(rho_approx_half_1B(2:end))); rho_approx_half_1B];
        %[err_squared_rho,theta_min]=circ_error_continuous_unrestricted(rho_approx_symm_1B,rho_hat_symm_1B,B)
        
        a_approx_symm_1B=M_1./(2*pi.*rho_approx_symm_1B);
        [err_squared_a_FM(sigma_index,f_index),theta_min]=circ_error_continuous_unrestricted(a_approx_symm_1B,a_symm_1B,B);
        
        
        %% Spectral algorithm
        W=(1/sqrt(2*B+1))*dftmtx(2*B+1);
        h=(1/(2*B+1))*(k_half_2B.'.*[0 ; conj(flip(rho_hat_half_2B(2:end)))]...
            +(2*B+1-k_half_2B.').*rho_hat_half_2B);
        C_h=circul(h);
        T=toeplitz(conj(rho_hat_half_2B));
        
        M_2_T=M_2_true-(sigma^2/(2*pi))*(k_symm_1B==k_symm_1B.'|k_symm_1B==-k_symm_1B.');
        M_2_C=2*pi*diag(a_symm_1B)*C_h*diag(conj(a_symm_1B));
        P_a_symm_1B_true=abs(a_symm_1B).^2;
        P_a_symm_1B_emp = diag(M_2_emp_wave);
        a_symm_1B_wave_true = a_symm_1B./sqrt(P_a_symm_1B_true);
        a_symm_1B_wave_emp = a_symm_1B./sqrt(P_a_symm_1B_emp);
        M_2_T_true_mathcal=(1/(2*pi))*W'*diag(1./sqrt(P_a_symm_1B_true))*...
            (M_2_true_wave)...
            *diag(1./sqrt(P_a_symm_1B_true))*W;
        
        M_2_T_emp_mathcal=(1/(2*pi))*W'*diag(1./sqrt(P_a_symm_1B_emp))*...
            (M_2_emp_wave)...
            *diag(1./sqrt(P_a_symm_1B_emp))*W;
        
        M_2_C_mathcal=(1/(2*pi))*W'*diag(1./sqrt(P_a_symm_1B_true))*M_2_C*diag(1./sqrt(P_a_symm_1B_true))*W;
        
        %norm(M_2_T_mathcal-circul(ifft(a_symm_1B_wave))*W'*T*W*circul(ifft(a_symm_1B_wave))')
        
        M_2_T_true_mathcal=(M_2_T_true_mathcal'+M_2_T_true_mathcal)/2; %force hermitian
        M_2_T_emp_mathcal=(M_2_T_emp_mathcal'+M_2_T_emp_mathcal)/2;
        M_2_C_mathcal=(M_2_C_mathcal'+M_2_C_mathcal)/2;
        
        [V_T_true,D_T_true]=eig(M_2_T_true_mathcal,"vector");
        [V_T_emp,D_T_emp]=eig(M_2_T_emp_mathcal,"vector");
        [V_C,D_C]=eig(M_2_C_mathcal,"vector");
        
        [D_T_true_sorted,I_T_true_sort]=sort(D_T_true,"descend");
        [D_T_emp_sorted,I_T_emp_sort]=sort(D_T_emp,"descend");
        [D_C_sorted,I_C_sort]=sort(D_C,"descend");
        
        V_T_true_sorted=V_T_true(:,I_T_true_sort);
        V_T_emp_sorted=V_T_emp(:,I_T_emp_sort);
        V_C_sorted=V_C(:,I_C_sort);
        
        [~,kappa_true]=max(min(abs(D_T_true_sorted-D_T_true_sorted.')+diag(inf(1,2*B+1)),[],2));
        [~,kappa_emp]=max(min(abs(D_T_emp_sorted-D_T_emp_sorted.')+diag(inf(1,2*B+1)),[],2));
        %only with kappa_true
        delta_kappa=max(min(abs(D_T_true_sorted(kappa_true)-D_C_sorted(1:end~=kappa_true))),...
            min(abs(D_C_sorted(kappa_true)-D_T_true_sorted(1:end~=kappa_true))));
        
        
        
        u_true=fft(V_T_true_sorted(:,kappa_true));
        u_emp=fft(V_T_emp_sorted(:,kappa_emp));
        
        u_true_renormilized_bad=exp(1i* (angle(a_symm_1B(k_symm_1B==0))-angle(u_true(k_symm_1B==0)))) *...
            sqrt(P_a_symm_1B_true).* u_true;
        u_emp_renormilized_bad=exp(1i* (angle(a_symm_1B(k_symm_1B==0))-angle(u_emp(k_symm_1B==0)))) *...
            sqrt(P_a_symm_1B_emp).* u_emp;
        
        [err_squared_a_Spectral_bad_true(sigma_index,f_index),l_bad_true]=circ_error_continuous(u_true_renormilized_bad,a_symm_1B,B);
        [err_squared_a_Spectral_bad_emp(sigma_index,f_index),l_bad_emp]=circ_error_continuous(u_emp_renormilized_bad,a_symm_1B,B);
        
        u_true_renormilized_good=exp(1i* (angle(a_symm_1B(k_symm_1B==0))-angle(u_true(k_symm_1B==0)))) *...
            sqrt(P_a_symm_1B_true).* u_true./abs(u_true);
        u_emp_renormilized_good=exp(1i* (angle(a_symm_1B(k_symm_1B==0))-angle(u_emp(k_symm_1B==0)))) *...
            sqrt(P_a_symm_1B_emp).* u_emp./abs(u_emp);
        
        [err_squared_a_Spectral_good_true(sigma_index,f_index),l_good_true]=circ_error_continuous(u_true_renormilized_good,a_symm_1B,B);
        [err_squared_a_Spectral_good_emp(sigma_index,f_index),l_good_emp]=circ_error_continuous(u_emp_renormilized_good,a_symm_1B,B);
    end
    dist_from_circ_true(f_index)=sum(abs((conj(rho_hat_half_2B(2:end))-flip(rho_hat_half_2B(2:end))).^2)...
        .*(k_half_2B(2:end).'.*(flip(k_half_2B(2:end).')))./(2*B+1));
    
    %abs(norm(T-C_h,"fro").^2-dist_from_circ(index))
    bound_rho=(2*B+1)*2*(1-sqrt(1-dist_from_circ_true(f_index)./(delta_kappa.^2)));
    bound(f_index)=bound_rho*max(abs(a_symm_1B)).^2;
    if imag(bound(f_index))~=0
        bound(f_index) = -1;
    end
end
close all

energy_of_a=sum(P_a_symm_1B_true);
SNR=energy_of_a./((2*B+1)*sigma_vec_reduced.^2);
figure
for f_index=1:length(f_vec)
    subplot(2,2,f_index)
    lg_cells=cell(1,4);
    
    err_FM_for_plot=mean(reshape(err_squared_a_FM(:,f_index),num_rep,num_unique_sigma),1).'./energy_of_a;
    loglog(SNR,err_FM_for_plot,'-blue',"linewidth",1.15);
    lg_cells{min(find(cellfun(@isempty,lg_cells)))}=...
        "$\textnormal{Frequency Marching Algorithm}$";
    
    hold on
    
    err_spectral_for_plot=mean(reshape(err_squared_a_Spectral_bad_emp(:,f_index),num_rep,num_unique_sigma),1).'./energy_of_a;
    loglog(SNR,err_spectral_for_plot,'-red',"linewidth",1.15);
    lg_cells{min(find(cellfun(@isempty,lg_cells)))}=...
        "$\textnormal{Spectral Algorithm}$";
    
    if bound(f_index)~=0
        yline(bound(f_index)./energy_of_a,"--black","linewidth",1.15);
        lg_cells{min(find(cellfun(@isempty,lg_cells)))}=...
            "$\textnormal{Bound on Noiseless Spectral Algorithm} = "+round(bound(f_index)./energy_of_a,2)+"$";
    end
    
    %legend("Frequency Marching Algorithm","Spectral Algorithm", "Bound on Noiseless Spectral Algorithm = "+bound(f_index)./energy_of_a,...
    %    "fontsize",10.5,"location","southwest")
    %xlabel("$\textnormal{SNR}=\frac{N\cdot \textnormal{Total Energy}}{\sigma^2}$","interpreter","latex","fontsize",15)
    %ylabel("$\frac{||\underline{a}_{\textnormal{est}}-\underline{a}||^2}{\textnormal{Total Energy}}$","interpreter","latex","fontsize",20)
    xlabel("$\textnormal{SNR}$","fontsize",15,"interpreter","latex")
    ylabel("$\textnormal{Relative Squared Error}$","fontsize",15,"interpreter","latex")
    title("$||T-C_{\underline{h}}||^2 = "+round(dist_from_circ_true(f_index),4)+"$",...
        "interpreter","latex","fontsize",15)
    set(gca,'XLim',[min(SNR), max(SNR)],"YLim",...
        [min([err_FM_for_plot; err_spectral_for_plot]), max([err_FM_for_plot; err_spectral_for_plot])]);
    
    %SNR_critical_point= min(SNR(SNR>=1e4 & err_spectral_for_plot>=err_FM_for_plot));
    cutoff=1e-3;
    [~,SNR_critical_point_index]=min(abs(err_spectral_for_plot(SNR>=cutoff)-err_FM_for_plot(SNR>=cutoff)));
    SNR_half=SNR(SNR>=cutoff);
    SNR_critical_point=SNR_half(SNR_critical_point_index);
    if ~isempty(SNR_critical_point) & SNR_critical_point ~= max(SNR)
        xline(SNR_critical_point,"--","color",[0.5 0.5 0.5],"linewidth",1.15);
        loglog(SNR_critical_point,err_spectral_for_plot(SNR==SNR_critical_point),...
            "o","color","m","linewidth",1.2);
        
        %lg_cells{min(find(cellfun(@isempty,lg_cells)))} = "Critical SNR = "+sprintf("%10e",SNR_critical_point);
        lg_cells{min(find(cellfun(@isempty,lg_cells)))} ="$\textnormal{Critical SNR} = "...
            +round(SNR_critical_point./cutoff,2) + "\times {10}^{"+log10(cutoff)+"}$";
    end
    lg_cells(find(cellfun(@isempty,lg_cells)))=[];
    lg=legend(lg_cells);
    lg.FontSize=10.5;
    lg.Location="southwest";
    lg.Interpreter="Latex";
    %legend(cat(2,lg_cells,{"fontsize",10.5,"location","southwest"}))
end
sgtitle("$\textnormal{FM Algorithm vs. Spectral Algorithm}$","interpreter","latex","fontsize",25)
toc