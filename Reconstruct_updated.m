close all
picture=imread('flower-1.jpg');
Full_Picture=im2double(rgb2gray(imresize(picture,103/1000)));

%notice: for even L no stupid summation of two lines, for odd L there is
%this bug.

L=size(Full_Picture,1);
radius=(L)/2;
r=@(x,y) sqrt((x-radius).^2+(y-radius).^2);
theta=@(x,y) atan((y-radius)./(x-radius));
%{
figure
imagesc(Full_Picture);
%}
%% Circular Picture
%{
Circle_Picture=Full_Picture;
for i=1:L
    for j=1:L
        if(r(i,j)>radius)
            Circle_Picture(i,j)=0;
        end
    end
end
%now Circle_Picture is a circular pic;

figure
imagesc(Circle_Picture);
%}
%%

base = fb_basis([L L], inf,0);
coeffVec=base.expand(Full_Picture);
Approx_Full_Picture=base.evaluate(coeffVec);


figure
imagesc(Approx_Full_Picture);
colorbar

%{
figure
imagesc(Full_Picture-Approx_Full_Picture);
colorbar
%}


ks=base.indices.ks;
ells=base.indices.ells;
sgns=base.indices.sgns;

alpha=0.5*pi;
new_coeffVec=zeros(base.count,1); %spin by alpha

k_max=base.k_max;
ell_max=base.ell_max;
for ell=0:ell_max
    for k=1:(k_max(ell+1))
        
        a_lk_index=(sgns==1)&(ks==k)&(ells==ell); %cosines
        b_lk_index=(sgns==-1)&(ks==k)&(ells==ell); %sines
        
        a_lk=coeffVec(a_lk_index);
        b_lk=coeffVec(b_lk_index);
        
        if isempty(b_lk)
            b_lk=0;
        end
        %spin
        a_lk_new=a_lk*cos(ell*alpha)-b_lk*sin(ell*alpha);
        b_lk_new=a_lk*sin(ell*alpha)+b_lk*cos(ell*alpha);
        
        new_coeffVec(a_lk_index)=a_lk_new;
        new_coeffVec(b_lk_index)=b_lk_new;
    end
end

Approx_Rotated_Picture=base.evaluate(new_coeffVec);
figure
imagesc(Approx_Rotated_Picture);
colorbar