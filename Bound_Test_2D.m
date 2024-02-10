% 2D case - Plots a graph of the analitic bound compared to the empirical errors of
% of the spectral algorithm, for exact given moments, for different
% distributions.

B=10;
Q=2;
M=2*B+1;
image_size=51;
W_M=(1/sqrt(M))*dftmtx(M);
W_Q=(1/sqrt(Q))*dftmtx(Q);
W=kron(W_M,W_Q);
%% Relevant frequencies
k_symm_1B=-B:B;
k_half_1B=0:B;
k_symm_2B=-2*B:2*B;
k_half_2B=0:2*B;

k_symm_1B_Q_mat=repmat(k_symm_1B,[Q,1]);
k_half_1B_Q_mat=repmat(k_half_1B,[Q,1]);
k_symm_2B_Q_mat=repmat(k_symm_2B,[Q,1]);
k_half_2B_Q_mat=repmat(k_half_2B,[Q,1]);

k_symm_1B_Q_vec=vec(k_symm_1B_Q_mat);
k_half_1B_Q_vec=vec(k_half_1B_Q_mat);
k_symm_2B_Q_vec=vec(k_symm_2B_Q_mat);
k_half_2B_Q_vec=vec(k_half_2B_Q_mat);

q_full_1B_mat=repmat((0:Q-1).',[1,2*B+1]);
q_full_1B_vec=vec(q_full_1B_mat);

%% Definition of Picture
%{
%a_mat=exp(2*pi*1i*rand(Q,B));
a_mat=[0.783747689358183 + 0.621079350345597i,0.999576110044461 + 0.0291135746411087i,0.410387184340826 - 0.911911376685700i,0.862541008291170 + 0.505987162896502i,-0.0619778226989724 + 0.998077526795136i,-0.908573643505186 + 0.417724711176763i,0.415250137373605 + 0.909707273473837i,0.610227673494066 + 0.792226095569958i,0.681298484913590 - 0.732005720233419i,-0.951327583025404 - 0.308181488371775i;0.971480085211998 - 0.237121158980994i,0.155878799637023 - 0.987776189135839i,0.678545581509888 - 0.734558298444315i,-0.808213527760008 + 0.588889542737620i,0.309426180590958 - 0.950923466302568i,0.846501185527152 - 0.532386835769938i,-0.0866176052946857 + 0.996241632563614i,0.656260788186847 + 0.754534146270650i,-0.877199368559842 - 0.480126304005743i,0.613131441753812 + 0.789980908081260i];
%a_mat=[-0.646641616810157 - 0.762793956064904i;0.877812699563756 - 0.479004033891773i;0.0401582335407562 + 0.999193332783444i;0.816802514387931 + 0.576917370590931i;-0.308579375282474 - 0.951198596061978i;0.375659048440682 + 0.926757939984677i;-0.839511784037794 - 0.543341480527376i;0.805090157386554 - 0.593152458040337i;0.780834669228819 + 0.624737720431799i;-0.400828934416310 - 0.916152915912341i].';
a_k_0=ones(Q,1);
a_symm_1B_mat=[conj(flip(a_mat,2)) a_k_0 a_mat];
a_symm_1B=vec(a_symm_1B_mat);
%}
isUniformPowerSpectrum=[0,1];
images=cell(length(isUniformPowerSpectrum),1);
a_symm_1B_cell=cell(length(isUniformPowerSpectrum),1);

[a_symm_1B_cell{1},images{1}]=Generate_Picture_cut('.\Flower_Images\Oxalis_tetraphylla_flower.jpg',B,Q,image_size,isUniformPowerSpectrum(1));
[a_symm_1B_cell{2},images{2}]=Generate_Picture_cut('.\Flower_Images\Oxalis_tetraphylla_flower.jpg',B,Q,image_size,isUniformPowerSpectrum(2));

%% Spectral Algorithm
%f_vec=logspace(-3,-20,20).';
f_vec=logspace(log10(500*B),log10(10*B),40).';
bound=zeros(length(f_vec),length(isUniformPowerSpectrum));
dist_from_circ=zeros(length(f_vec),length(isUniformPowerSpectrum));
err_squared_T=zeros(length(f_vec),length(isUniformPowerSpectrum));
P_a_cell=cellfun(@(x) abs(x).^2,a_symm_1B_cell,"UniformOutput",false);
for UPS=[1,2]
    a_symm_1B=a_symm_1B_cell{UPS};
    a_symm_1B_mat=mat(a_symm_1B,Q,B);

    for f_index=1:length(f_vec)
        %% Rho Definition
        %temp=[0.00195583374369666 + 0.00530869159003379i;0.00195583374369666 + 0.00474988194897760i;0.00195583374369666 + 0.00419107230792141i;0.00195583374369666 + 0.00363226266686522i;0.00195583374369666 + 0.00307345302580903i;0.00195583374369666 + 0.00251464338475285i;0.00195583374369666 + 0.00195583374369666i;0.00195583374369666 + 0.00139702410264047i;0.00195583374369666 + 0.000838214461584282i;0.00195583374369666 + 0.000279404820528094i;0.00195583374369666 - 0.000279404820528094i;0.00195583374369666 - 0.000838214461584282i;0.00195583374369666 - 0.00139702410264047i;0.00195583374369666 - 0.00195583374369666i;0.00195583374369666 - 0.00251464338475285i;0.00195583374369666 - 0.00307345302580903i;0.00195583374369666 - 0.00363226266686522i;0.00195583374369666 - 0.00419107230792141i;0.00195583374369666 - 0.00474988194897760i;0.00195583374369666 - 0.00530869159003379i]...
        %    .*exp(1i*(f_vec(f_index)*rand(2*B,1)).^0.2)*3;
        temp=[0.00195583374369666 + 0.00530869159003379i;0.00195583374369666 + 0.00474988194897760i;0.00195583374369666 + 0.00419107230792141i;0.00195583374369666 + 0.00363226266686522i;0.00195583374369666 + 0.00307345302580903i;0.00195583374369666 + 0.00251464338475285i;0.00195583374369666 + 0.00195583374369666i;0.00195583374369666 + 0.00139702410264047i;0.00195583374369666 + 0.000838214461584282i;0.00195583374369666 + 0.000279404820528094i;0.00195583374369666 - 0.000279404820528094i;0.00195583374369666 - 0.000838214461584282i;0.00195583374369666 - 0.00139702410264047i;0.00195583374369666 - 0.00195583374369666i;0.00195583374369666 - 0.00251464338475285i;0.00195583374369666 - 0.00307345302580903i;0.00195583374369666 - 0.00363226266686522i;0.00195583374369666 - 0.00419107230792141i;0.00195583374369666 - 0.00474988194897760i;0.00195583374369666 - 0.00530869159003379i]...
        .*exp(1i*(rand(2*B,1)).^0.2./f_vec(f_index))*3;
        %temp=[0.00195583374369666 + 0.00530869159003379i;0.00195583374369666 + 0.00474988194897760i;0.00195583374369666 + 0.00419107230792141i;0.00195583374369666 + 0.00363226266686522i;0.00195583374369666 + 0.00307345302580903i;0.00195583374369666 + 0.00251464338475285i;0.00195583374369666 + 0.00195583374369666i;0.00195583374369666 + 0.00139702410264047i;0.00195583374369666 + 0.000838214461584282i;0.00195583374369666 + 0.000279404820528094i;0.00195583374369666 - 0.000279404820528094i;0.00195583374369666 - 0.000838214461584282i;0.00195583374369666 - 0.00139702410264047i;0.00195583374369666 - 0.00195583374369666i;0.00195583374369666 - 0.00251464338475285i;0.00195583374369666 - 0.00307345302580903i;0.00195583374369666 - 0.00363226266686522i;0.00195583374369666 - 0.00419107230792141i;0.00195583374369666 - 0.00474988194897760i;0.00195583374369666 - 0.00530869159003379i]...
        %.*exp(1i*(0.8*ones(2*B,1)).^0.3./f_vec(f_index).^2)*3;
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

        varrho=repelem(rho_hat_symm_1B,Q);

        h=(1/(2*B+1))*(k_half_2B.'.*[0 ; conj(flip(rho_hat_half_2B(2:end)))]...
        +(2*B+1-k_half_2B.').*rho_hat_half_2B);
        zeta=repmat(h.',[Q,1]);

        %%
        C=BCCB(zeta);
        T=BTTB(repmat(rho_hat_half_2B.',Q,1));

        a_tilde_symm_1B=a_symm_1B./abs(a_symm_1B);
        %A=BCCB(ifft2(mat(a_tilde_symm_1B,Q,B)));

        P_a=P_a_cell{UPS};
        v_p=1./sqrt(P_a);

        M1=2*pi.*a_symm_1B.*varrho;

        M_2_C=2*pi*diag(a_symm_1B)*C*diag(a_symm_1B)';
        M_2_T=2*pi*diag(a_symm_1B)*T*diag(a_symm_1B)';
        M_2_C_mathcal=(1/(2*pi))*W'*diag(v_p)*M_2_C*diag(v_p)*W;
        M_2_T_mathcal=(1/(2*pi))*W'*diag(v_p)*M_2_T*diag(v_p)*W;

        M_2_T_mathcal=(M_2_T_mathcal'+M_2_T_mathcal)/2; %force hermitian
        M_2_C_mathcal=(M_2_C_mathcal'+M_2_C_mathcal)/2;

        [V_T,D_T]=eig(M_2_T_mathcal,"vector");
        [V_C,D_C]=eig(M_2_C_mathcal,"vector");

        [D_T_sorted,I_T_sort]=sort(D_T,"descend");
        [D_C_sorted,I_C_sort]=sort(D_C,"descend");

        V_T_sorted=V_T(:,I_T_sort);
        V_C_sorted=V_C(:,I_C_sort);

        [~,kappa]=max(min(abs(D_T_sorted-D_T_sorted.')+diag(inf(1,Q*(2*B+1))),[],2));
        delta_kappa=max(min(abs(D_T_sorted(kappa)-D_C_sorted(1:end~=kappa))),...
            min(abs(D_C_sorted(kappa)-D_T_sorted(1:end~=kappa))));

        v_C_kappa=V_C_sorted(:,kappa); %unique eigenvector
        v_T_kappa=V_T_sorted(:,kappa);

        a_est_tilde_C=vec(fft2(mat(v_C_kappa,Q,B)));
        a_est_C=a_est_tilde_C...
            .*exp(1i.*(angle(M1(k_symm_1B_Q_vec==0 & q_full_1B_vec==0))...
            -angle(a_est_tilde_C(k_symm_1B_Q_vec==0 & q_full_1B_vec==0))))...
            .*sqrt(P_a);
        [err_squared_C,l_C]=Circ_Error_Continuous_2D(a_est_C,a_symm_1B,Q,B);

        a_est_tilde_T=vec(fft2(mat(v_T_kappa,Q,B)));
        a_est_T=a_est_tilde_T...
            .*exp(1i.*(angle(M1(k_symm_1B_Q_vec==0 & q_full_1B_vec==0))...
            -angle(a_est_tilde_T(k_symm_1B_Q_vec==0 & q_full_1B_vec==0 ))))...
            .*sqrt(P_a);
        [err_squared_T(f_index,UPS),l_T]=Circ_Error_Continuous_2D(a_est_T,a_symm_1B,Q,B);

        dist_from_circ(f_index,UPS)=(Q^2)*sum(abs((conj(rho_hat_half_2B(2:end))-flip(rho_hat_half_2B(2:end))).^2)...
            .*(k_half_2B(2:end).'.*(flip(k_half_2B(2:end).')))./(2*B+1));

        %dist_from_circ(f_index,UPS)=norm(T-C,"fro").^2;
        bound(f_index,UPS)=2*Q*(2*B+1)*max(P_a).*(1-sqrt(1-dist_from_circ(f_index,UPS)./delta_kappa.^2));
        if imag(bound(f_index,UPS))~=0
            bound(f_index,UPS)=-1;
        end
    end
end
%% Plots
close all
figures=cell(length(isUniformPowerSpectrum),1);
image_figures=cell(length(isUniformPowerSpectrum),1);
for UPS=[1,2]
    figures{UPS}=figure;
    [dist_sorted,I_dist]=sort(dist_from_circ(:,UPS),"descend");
    %loglog(dist_sorted,err_squared_T(I_dist),"o")
    err_relative_db=10*log10(err_squared_T(I_dist,UPS)./sum(P_a_cell{UPS}));
    bound_relative_db=10*log10(bound(I_dist,UPS)./sum(P_a_cell{UPS}));
    dist_db=10*log10(dist_sorted./(Q^2)); %div by Q^2
    plot(dist_db,err_relative_db...
        ,".black","MarkerSize",10)
    hold on
    %loglog(dist_sorted,bound(I_dist),"*")
    plot(dist_db,bound_relative_db,"-blue","linewidth",1.2)
    axis tight
    set(gca,'xdir','reverse')
    legend("Relative Error[dB]","Bound[dB]","Interpreter","latex","Fontsize",15)
    xlabel("$||T-C_{\underline{h}}||^2_F[\textnormal{dB}]$","Interpreter","latex","Fontsize",15)
    if isUniformPowerSpectrum(UPS)==1
        title_string='With Uniform Spectrum';
        figure_string='With_Uniform_Spectrum';
    else
        title_string='Without Uniform Spectrum';
        figure_string='Without_Uniform_Spectrum';
    end
    title(title_string,"Interpreter","latex","FontSize",20)
    saveas(figures{UPS},[pwd '/Figures_Thesis/' figure_string '.fig']);
    saveas(figures{UPS},[pwd '/Figures_Thesis/' figure_string '.png']);
    print(figures{UPS},'-depsc',[pwd '/Figures_Thesis/' figure_string '.eps']);
    image_figures{UPS}=figure;
    imagesc(images{UPS})
    set(gca,'ydir','normal')
    axis off
    title(title_string,"Interpreter","latex","FontSize",20)
    colormap gray
    colorbar
    saveas(image_figures{UPS},[pwd '/Figures_Thesis/' figure_string '_image' '.fig']);
    saveas(image_figures{UPS},[pwd '/Figures_Thesis/' figure_string '_image' '.png']);
    print(image_figures{UPS},'-depsc',[pwd '/Figures_Thesis/' figure_string '_image' '.eps']);
end
