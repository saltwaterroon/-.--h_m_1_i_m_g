clc;clear all
I1=im2double(imread('/Users/h-xiao16/Desktop/数字图像处理/综合作业1/phone.bmp'));
block_size=16;
kernel_size=8;
repeat=block_size-kernel_size;
I1o=imresize(I1,[512+repeat 512+repeat],'bicubic');
graph=zeros(size(I1o));
figure(1),imshow(I1);
num=(512)/kernel_size;
I2=double(imread('/Users/h-xiao16/Desktop/数字图像处理/综合作业1/phone.bmp'));
I3=double(imread('/Users/h-xiao16/Desktop/数字图像处理/综合作业1/latent.bmp'));
 % h=fspecial('gaussian',3,3);
I1=Normalization(I1,0.1,2,repeat);
I1=histeq(I1o);
%I1=histeq(I1o);
%I1=imfilter(I1,h);
figure(2),imshow(I1,[]);
im_block=cell(num,num);
im_block0=cell(num,num);
I4=zeros(size(I1));
for i=1:num
    for j=1:num
       % im_block0{i,j}=I1(max(1,(-19+8*i)):min(512,(8*(i)+12)),max(1,(-19+8*j)):min((8*(j)+12),512));
               im_block0{i,j}=I1((kernel_size*(i-1)+1):kernel_size*i+repeat,(kernel_size*(j-1)+1):kernel_size*j+repeat);

    end
end
fs=zeros(num,num);
agl=zeros(num,num);
bg=ones(num,num);
im_dft_abs=cell(num,num);
im_dft_angle=cell(num,num);
for i=1:num
    for j=1:num
        im_block{i,j}=im_block0{i,j}-mean2(im_block0{i,j});%减去直流分量
        temp=im_block0{i,j};
        F1=fft2(im_block{i,j});
        V(i,j)=var(im_block{i,j}(:));
        im_dft_abs{i,j}=abs(fftshift(fft2(im_block{i,j})));
        im_dft_angle{i,j}=angle(fftshift(fft2(im_block{i,j})));
        im_dft_abs{i,j}=(im_dft_abs{i,j}-min(im_dft_abs{i,j}(:)))/(max(im_dft_abs{i,j}(:))-min((im_dft_abs{i,j}(:))));
%        im_abs_max=max(im_dft_abs{i,j});
 %       [m,centre]=max(im_abs_max);%centre为幅度谱中心的位置
        [x,y] = sort(im_dft_abs{i,j}(:),'descend');
        [M,N]=size(im_block{i,j});
        im_dft_abs{i,j}(M/2+1,N/2+1) = 0;
           [x1,y1] = ind2sub(size(im_dft_abs{i,j}),y(1));
           [x2,y2] = ind2sub(size(im_dft_abs{i,j}),y(2));
           %[x_centre,y_centre] = ind2sub(size(im_dft_abs{i,j}),m);
%           if((y1+y2==2*centre)&&(m-x(2)<0.975))
 %                  I4((8*(i-1)+1):8*i,(8*(j-1)+1):8*j)=I1((8*(i-1)+1):8*i,(8*(j-1)+1):8*j);
  %      
   %        end
   if(V(i,j)>0.01)
             if((x(1)==x(2))&&(y(1)+y(2)==y(3)+y(4))&&(x(3)==x(4)))%是指纹
                 bg(i,j)=0;
      %   I4((8*(i-1)+1):8*i,(8*(j-1)+1):8*j)=I1((8*(i-1)+1):8*i,(8*(j-1)+1):8*j);
        
      %  temp=I4((8*(i-1)+1):8*i,(8*(j-1)+1):8*j);
        temp=im_block{i,j};
        fs(i,j)=(sqrt((x1-x2)^2+(y1-y2)^2))*0.5*2*pi;%频率图
        agl(i,j)=atan((y1-y2)/(x1-x2));%方向图
          %  if(agl(i,j)>pi/2)
           %     agl(i,j)=pi+agl(i,j);
            %end
           % if(fs(i,j)>4)
            %    fs(i,j)=-1;
           % end
           
        %{
        h=fspecial('sobel');
        Gy=filter2(h,temp);  % Gradient magnitudes along the y direction
        Gx=filter2(transpose(h),temp); % Gradient magnitudes along the x direction
        Gx2=reshape(Gx,1,num);
        Gy2=reshape(Gy,1,num);
        Vx=2*Gx2*((Gy2)');
        Vx2=num*mean2(2*Gx.*Gy);
        Vy=num*mean2((Gx.^2-Gy.^2));
        theta(i,j)=atand(Vy/Vx);
             
        %}
            end
       % [xi,yi]=sort(im_abs_max,'descend');
       % if((im_abs_max(centre-1)==im_abs_max(centre+1))&&(im_abs_max(centre+1)>0))
   end  
    end
        %subplot(num,num,i+num*(j-1)),imshow(im_dft_abs{i,j});
end


        h = fspecial('gaussian',[5,5],1.2);
        newaglcos=cos(2*agl);
        newaglsin=sin(2*agl);
        newfs=cos(2*fs)+1i*sin(2*fs);
            aglcos2=imfilter(newaglcos,h);  
            aglsin2=imfilter(newaglsin,h);  
            newfs2=imfilter(fs,h);
            result_o=0.5*atan2(aglsin2,aglcos2);
            result_fs=0.5*angle(newfs2);
            for i=1:num
                for j=1:num
                    
                      %  if(result_o(i,j)>0.01&&result_o(i,j)<pi/2) result_o(i,j)=pi+result_o(i,j);
                      %  end
                        
                end
            end

figure(2),quiver(1:num,1:num,cos(agl),sin(agl));set(gca,'YDir','reverse');axis image;
figure(3),quiver(1:num,1:num,cos(result_o),sin(result_o));set(gca,'YDir','reverse');axis image;
figure(4),imshow(newfs2,[]);
%figure(5),imshow(freq,[]);
%figure(4),imshow(fre,[]);
agl=result_o;
fs2=newfs2;



for i=1:64
    for j=1:64
        im_block{i,j}=im_block0{i,j}-mean2(im_block0{i,j});%减去直流分量
        temp=im_block0{i,j};
        F1=fft2(im_block{i,j});
        F0=fft2(im_block0{i,j});
        V(i,j)=var(im_block{i,j}(:));
        im_dft_abs{i,j}=abs(fftshift(fft2(im_block{i,j})));
        im_dft_angle{i,j}=angle(fftshift(fft2(im_block{i,j})));
        im_dft_abs{i,j}=(im_dft_abs{i,j}-min(im_dft_abs{i,j}(:)))/(max(im_dft_abs{i,j}(:))-min((im_dft_abs{i,j}(:))));
%        im_abs_max=max(im_dft_abs{i,j});
 %       [m,centre]=max(im_abs_max);%centre为幅度谱中心的位置
        [x,y] = sort(im_dft_abs{i,j}(:),'descend');
        [M,N]=size(im_block{i,j});
        im_dft_abs{i,j}(M/2+1,N/2+1) = 0;
           [x1,y1] = ind2sub(size(im_dft_abs{i,j}),y(1));
           [x2,y2] = ind2sub(size(im_dft_abs{i,j}),y(2));
           %[x_centre,y_centre] = ind2sub(size(im_dft_abs{i,j}),m);
%           if((y1+y2==2*centre)&&(m-x(2)<0.975))
 %                  I4((8*(i-1)+1):8*i,(8*(j-1)+1):8*j)=I1((8*(i-1)+1):8*i,(8*(j-1)+1):8*j);
  %      
   %        end
   if(V(i,j)>0.01)
             if((x(1)==x(2))&&(y(1)+y(2)==y(3)+y(4))&&(x(3)==x(4)))%是指纹
                    agl(i,j)=0.5*atan2(aglsin2(i,j),aglcos2(i,j));
            [mag, phase] = imgaborfilt(im_block0{i,j},  max(fs2(i, j),2), (agl(i,j)+pi/2) * 180 / pi);
                temp = mag.*cos(phase); 
        %{
        h=fspecial('sobel');
        Gy=filter2(h,temp);  % Gradient magnitudes along the y direction
        Gx=filter2(transpose(h),temp); % Gradient magnitudes along the x direction
        Gx2=reshape(Gx,1,num);
        Gy2=reshape(Gy,1,num);
        Vx=2*Gx2*((Gy2)');
        Vx2=num*mean2(2*Gx.*Gy);
        Vy=num*mean2((Gx.^2-Gy.^2));
        theta(i,j)=atand(Vy/Vx);
             
        %}
            end
       % [xi,yi]=sort(im_abs_max,'descend');
       % if((im_abs_max(centre-1)==im_abs_max(centre+1))&&(im_abs_max(centre+1)>0))
   end  
        if(i == 1)
            if(j == 1)
                graph(1:block_size,1:block_size) = temp;
            else
                 graph(1:block_size, ((j - 1) * (block_size - repeat) + 1) : ((j - 1) * (block_size - repeat) + block_size))...
                    =  graph(1:block_size, ((j - 1) * (block_size - repeat) + 1) : ((j - 1) * (block_size - repeat) + block_size)) + temp; %先直接加上
                graph(1:block_size, ((j - 1) * (block_size - repeat) + 1) : ((j - 1) * (block_size - repeat) + repeat))...
                    = graph(1:block_size, ((j - 1) * (block_size - repeat) + 1) : ((j - 1) * (block_size - repeat) + repeat)) / 2; %再把左侧重复部分除2
            end
        else
            if(j ==1)
                 graph(((i - 1) * (block_size - repeat) + 1) : (i * (block_size - repeat) + repeat), 1:block_size)...
                     = graph(((i - 1) * (block_size - repeat) + 1) : (i * (block_size - repeat) + repeat), 1:block_size)...
                     + temp;%先直接加上
                 graph(((i - 1) * (block_size - repeat) + 1) : ((i - 1) * (block_size - repeat) + repeat), 1:block_size)...
                     =graph(((i - 1) * (block_size - repeat) + 1) : ((i - 1) * (block_size - repeat) + repeat), 1:block_size) / 2;%上方
            else
                graph(((i - 1) * (block_size - repeat) + 1) : (i * (block_size - repeat) + repeat), ((j - 1) * (block_size - repeat) + 1) : (j * (block_size - repeat) + repeat))...
                     = graph(((i - 1) * (block_size - repeat) + 1) : (i * (block_size - repeat) + repeat), ((j - 1) * (block_size - repeat) + 1) : (j * (block_size - repeat) + repeat))...
                     + temp;%先直接加上
                 graph(((i - 1) * (block_size - repeat) +1) : ((i - 1) * (block_size - repeat) + repeat), ((j - 1) * (block_size - repeat) + repeat + 1) : ((j - 1) * (block_size - repeat) + block_size))...
                     = graph(((i - 1) * (block_size - repeat) + 1) : ((i - 1) * (block_size - repeat) + repeat), ((j - 1) * (block_size - repeat) + repeat + 1) : ((j - 1) * (block_size - repeat) + block_size))...
                     / 2;%上方中右
                 graph(((i - 1) * (block_size - repeat) + 1) : ((i - 1) * (block_size - repeat) + block_size), ((j - 1) * (block_size - repeat) + 1) : ((j - 1) * (block_size - repeat) + repeat))...
                     = graph(((i - 1) * (block_size - repeat) + 1) : ((i - 1) * (block_size - repeat) + block_size), ((j - 1) * (block_size - repeat) + 1) : ((j - 1) * (block_size - repeat) + repeat))...
                     /2;%左方全
            end
        end
    end
        %subplot(num,num,i+num*(j-1)),imshow(im_dft_abs{i,j});
end
   %graph = imfilter(graph, fspecial('gaussian',ceil(block_size/2),1));
figure(10),imshow(graph);

%{
im=double(I4);
gb=cell(num,num);
sx=3;
sy=3;%dev
u=11;
v=11;






for i=1:num
    for j=1:num
        [x,y] = meshgrid(-5:5);
    xpi=-x*sin(result_o(i,j))+y*cos(result_o(i,j));
    ypi=x*cos(result_o(i,j))+y*sin(result_o(i,j));
        %figure(i+8*(j-1)),imshow(im_dft_abs{32+i,32+j},[]);
        gb{i,j}= exp(-(xpi.^2/sx^2 + ypi.^2/sy^2)/2).*cos(2*pi*newfs2(i,j)*xpi);
                    filter{i,j} = imrotate(gb{i,j},result_o(i,j),'crop'); 
                   

    end
end
enhanced=zeros(522,522);
pad=zeros(522,522);
pad(6:517,6:517)=im;
for i=6:517
    for j=6:517
        subi=ceil((i-5)/8);subj=ceil((j-5)/8);
         if(newfs2(subi,subj)~=0)
                        for m=-5:5
                            for block_size=-5:5
                                
        enhanced(i,j)=enhanced(i,j)+pad(i-m,j-block_size)*gb{subi,subj}(m+6,block_size+6);
                            end
                        end
                        
         end
     %   enhanced((8*(i-1)+1):8*i,(8*(j-1)+1):8*j)= imfilter(I1((8*(i-1)+1):8*i,(8*(j-1)+1):8*j),gb{i,j});

    end
end
 figure(6),imshow(enhanced,[]);
%}