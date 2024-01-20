[a_symm_1B,image]=Generate_Picture_cut([pwd '\Flower_Images\Oxalis_tetraphylla_flower.jpg'],B,Q,0);

a_symm_1B_mat=mat(a_symm_1B,Q,B);
    a_half_1B_mat=a_symm_1B_mat(:,k_symm_1B>=0);
    a_half_1B_vec=vec(a_half_1B_mat)