close all
a=rgb2gray(sunflower);

L=size(a,1);
c=(L-1)/2;
dx=1/L;
r=@(x,y) sqrt((x-c).^2+(y-c).^2);
theta=@(x,y) atan((y-c)./(x-c));

for i=1:L
    for j=1:L
        if(r(i,j)>c)
            a(i,j)=0;
        end
    end
end
%now a is a circular pic;

figure
imagesc(a);

%radPic=zeros(3,L^2);
%{
for i=1:L
    for j=1:L
        [theta0 r0]=cart2pol(i-c,j-c);
        radPic(:,L*(i-1)+j)=[theta0 ;r0; a(i,j)];
        %radPic(2,L*(i-1)+j)=theta(i,j);
        %radPic(3,L*(i-1)+j)=a(i,j);
    end
end
%}
