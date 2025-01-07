figure
loglog(N_vec_reduced, mean(reshape(M_1_error, [num_rep_N, num_unique_N]), 1)./energy_of_a);
hold on
loglog(N_vec_reduced, mean(reshape(M_2_error, [num_rep_N, num_unique_N]), 1)./energy_of_a);
loglog(N_vec_reduced, N_vec_reduced .^ (-1), "--")
title("Moment Errors for SNR=inf", "FontSize", 20);
ylabel("Relative Moment Errors", "FontSize", 15);
xlabel("N", "FontSize", 15);
legend("Relative Squared M_1 Error", "Relative Squared M_2 Error", "Slope -1", "FontSize", 10)
grid on
