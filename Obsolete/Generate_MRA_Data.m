function [data] = Generate_MRA_Data(rho,N,sigma,x,possible_rotations)
L=length(x);
%if isempty(possible_rotations)
 %   possible_rotations=1:L;
group_elements=randsample(possible_rotations,N,true,rho); %non-uniform rotations (from distribution rho)
data=zeros(L,N);
%rotate pictures
for i=1:N
    data(:,i)=circshift(x,group_elements(i));
end

%add noise
data=data+sigma*randn(L,N);

end

