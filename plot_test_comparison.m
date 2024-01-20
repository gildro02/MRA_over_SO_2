energy_of_a=sum(P_a_symm_1B_true);
SNR=energy_of_a./((2*B+1)*sigma_vec_reduced.^2);
figure
for f_index=2
    %subplot(2,2,f_index)
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
%sgtitle("$\textnormal{FM Algorithm vs. Spectral Algorithm}$","interpreter","latex","fontsize",25)