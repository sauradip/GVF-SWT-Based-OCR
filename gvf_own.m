% % <author>Vijeta Khare</author>
% % <email>kharevijeta@gmail.com</email>
% % <date>8-05-2017</date>
% % <summary>Contains a charachter segmentation of Licence Plate Images</summary> 
% %  All Right Reserved to Author

% clear all;
clc;
clear all;







str='torsotext25.jpg';

%Input Directory
Dir_img = dir('C:\Users\SAURADIP\Documents\MATLAB\html\GVF\ip');
%Output Directory
output_dir=('C:\Users\SAURADIP\Documents\MATLAB\html\GVF\op');
[img_num,img_n]=size(Dir_img);

bod=imread(str);
%I=imresize(bod,0.5); 
I=bod;
[BW1,thresh,v,h] = edge(rgb2gray(I),'sobel');
edgedir=atand(v./h);
BW2 = edge(rgb2gray(I),'canny');
%figure;

%imshowpair(BW1,BW2,'montage')
%title('Sobel Filter                                   Canny Filter');
%[mag,dc,fuse]=mgdc(bod);
[mgdir,mag,dc,fuse,gray]=mgdc(bod);
[sparse]=sparsegen(mag,0);
[sparse_dc]=sparsegen(dc,0);
% for mg

[df1,df2,select,centroid] = kmeans2(sparse);
select=round(double(select)/double(max(max(select))),1);
% for dc

[df11,df21,select1,centroid1] = kmeans2(sparse_dc);
select1=round(double(select1)/double(max(max(select1))),1);

% for mg
sparse2=zeros(size(select,1),size(select,2));
for ms=1:size(select,1)
     for ns=1:size(select,2)
         if(select(ms,ns)==1)
             sparse2(ms,ns)=sparse(ms,ns);
         else
              sparse2(ms,ns)=0;
         end
     end
end
sparsemg_nor=round(double(sparse)/double(max(max(sparse))),1);
sparse3=zeros(size(select1,1),size(select1,2));
for mss=1:size(select1,1)
     for nss=1:size(select1,2)
         if(select1(mss,nss)==1)
             sparse3(mss,nss)=dc(mss,nss);
         else
              sparse3(mss,nss)=0;
         end
     end
end

%[final,sparse6 ,sparse_mg,sparse7,sparse_dc1]=sparsegenmod(sparse2,sparse3,0);
[final]=sparsegenmod(sparse2,sparse3,0);

%sparse4=double(sparse4)/double(max(max(sparse4)));







% 
% [sparse3]=sparsegen(sparse2,0);
% G1=round(double(sparse3)/double(max(max(sparse3))),1);
% [df11,df22,select2,centroid2] = kmeans2(G1);
% 
sparse9=zeros(size(final,1),size(final,2));
for ms=1:size(final,1)
     for ns=1:size(final,2)
         if(final(ms,ns))
             sparse9(ms,ns)=sparse2(ms,ns);
         end
     end
end




text_weighted=round(double(sparse9)/double(max(max(sparse9))),1);
count=0;weight_left=0;weight_right=0;
 freq_left=0;
            freq_right=0;
             cont=0;
        siz=0;
det=zeros(size(final,1),size(final,2));
for ii=1:size(final,1)-2
    for jj=1:size(final,2)-2
        freq_left=0;
        freq_right=0;
        cont=0;
        siz=0;
 x_axis_t=ii; %52
 y_axis_t=jj; %24
 left_point=1;
        right_point=1;
 
 if(final(x_axis_t,y_axis_t)==1  ||  final(x_axis_t,y_axis_t)==0)
     
text_MG_normal = sparsemg_nor(x_axis_t:x_axis_t+2,y_axis_t:y_axis_t+2);


 matmg = text_MG_normal;
 %matmg=text_pixel;
 for xa=1:3
     for ya=1:3
         if(matmg(xa,ya)~=0 )
         if(round(matmg(xa,ya),1)<=0.1 && round(matmg(xa,ya),1)>=0.0) 
             matmg(xa,ya)=0.1;
         elseif((round(matmg(xa,ya),1)<=0.2 && round(matmg(xa,ya),1)>0.1 ) )
                matmg(xa,ya)=0.2;
         elseif((round(matmg(xa,ya),1)<=0.3 && round(matmg(xa,ya),1)>0.2 ) )  
                matmg(xa,ya)=0.3;
         elseif((round(matmg(xa,ya),1)<=0.4 && round(matmg(xa,ya),1)>0.3 ))  
                matmg(xa,ya)=0.4;
         elseif((round(matmg(xa,ya),1)<=0.5 && round(matmg(xa,ya),1)>0.4 ) )  
                matmg(xa,ya)=0.5;
         elseif((round(matmg(xa,ya),1)<=0.6 && round(matmg(xa,ya),1)>0.5 ))  
                matmg(xa,ya)=0.6;
          elseif((round(matmg(xa,ya),1)<=0.7 && round(matmg(xa,ya),1)>0.6 ))  
                matmg(xa,ya)=0.7;  
          elseif((round(matmg(xa,ya),1)<=0.8 && round(matmg(xa,ya),1)>0.7 ))  
                matmg(xa,ya)=0.8;
           elseif((round(matmg(xa,ya),1)<=0.9 && round(matmg(xa,ya),1)>0.8 ))  
                matmg(xa,ya)=0.9;
           elseif((round(matmg(xa,ya),1)<=1.0 && round(matmg(xa,ya),1)>0.9 ))  
                matmg(xa,ya)=1.0;
         end     
         end      
         
     end
 end
 
 d = unique(nonzeros(matmg));
     if(size(d)~=0)
MG_Freq_text_nor = [d,histc(nonzeros(matmg(:)),d)];
 
text_mg=MG_Freq_text_nor; 

w_freq=zeros(numel(text_mg(:,1)),2); % a frequency table of neighbours vs weight frequency


       % pos=find(text_mg==matmg(s,t)); %position of the value in frequency table
        w_freq(1:numel(text_mg(:,1)),2)= (text_mg(:,1)) ;
        w_freq(1:numel(text_mg(:,1)),1)= text_mg(:,2);
        
%      bar3( w_freq(:,2),w_freq(:,1)) 
%      set(gca,'YLim',[0 1])
    
        max_freq=max(w_freq(:,1));   
        position=find(w_freq(:,1)==max_freq);  % location of max frequency element in the array
        mid_point=ceil(size(w_freq,1)/2);
        if(size(w_freq,1)>2 & position < size(w_freq,1) & position >1)
            %mid_point=ceil(size(w_freq,1)/2);
            mid_point=position;
        left_point=mid_point-1;
        right_point=mid_point+1;
            freq_left=w_freq(1:left_point,1);
            freq_right=w_freq(right_point:size(w_freq,1),1);
            
            if(size(freq_left,1)==size(freq_right,1))
%                 if all(freq_left == freq_right)
%                 count=1;
%                 end

                for i=1:mid_point-1
                   if(freq_left(i,1)==freq_right(size(freq_left,1)-i+1,1))
                       count=count+1;
                       
                   end
                end
                
                
                weight_left=sum(sum(w_freq(left_point,1).*w_freq(left_point,2)));
                weight_right=sum(sum(w_freq(right_point,1).*w_freq(right_point,2)));
            else
                count=0;
            end
            
                
        elseif(size(w_freq,1)==2)
                
             if(w_freq(1,1)<w_freq(2,1))
                 left_point=1;
                  freq_left=w_freq(1:left_point,1);
             else
                 right_point=2;
                 freq_right=w_freq(right_point:size(w_freq,1),1);
             end
             
                  weight_left=sum(sum(w_freq(left_point,1).*w_freq(left_point,2)));
                weight_right=sum(sum(w_freq(right_point,1).*w_freq(right_point,2)));
            elseif(size(w_freq,1)>2 & position == size(w_freq,1) & position >1)  
                mid_point=position;
                left_point=mid_point-1;
                %right_point=mid_point+1;
                    freq_left=w_freq(1:left_point,1);
                    %freq_right=w_freq(right_point:size(w_freq,1),1);
                         weight_left=sum(sum(w_freq(left_point,1).*w_freq(left_point,2)));
                weight_right=sum(sum(w_freq(right_point,1).*w_freq(right_point,2)));
            elseif(size(w_freq,1)>2 & position < size(w_freq,1) & position ==1)  
                mid_point=position;
                %left_point=mid_point-1;
                right_point=mid_point+1;
                    %freq_left=w_freq(1:left_point,1);
                    freq_right=w_freq(right_point:size(w_freq,1),1);
                         weight_left=sum(sum(w_freq(left_point,1).*w_freq(left_point,2)));
                weight_right=sum(sum(w_freq(right_point,1).*w_freq(right_point,2)));
                
        end
        
        cont=numel(nonzeros(w_freq(:,1)));
        siz=size(w_freq,1);
       % if(size(weight_left,1)==size(weight_right,1))
      %         if( position==mid_point & cont==siz & count==mid_point-1 & size(w_freq,1)>2 & all(diff(freq_left))>0 & sum(diff(freq_right))<0 &  size(freq_left,1)==size(freq_right,1) & (weight_left+weight_right)<max_freq*w_freq(mid_point,1)& numel(unique(freq_left))==numel(unique(freq_right)) )

        if(  position==mid_point & cont==siz & size(w_freq,1)>2 & mod(size(w_freq,1),2)~=0 & ( (all(diff(freq_left))>0) | (sum(diff(freq_right))<0)  )  & numel(unique(freq_left))==numel(unique(freq_right))& (weight_left+weight_right)<max_freq*w_freq(mid_point,1)  &  size(freq_left,1)==size(freq_right,1) )
            det(x_axis_t,y_axis_t)=1;
        else
            det(x_axis_t,y_axis_t)=0;
        end
       
% left_point=1;
%         right_point=1;
        

     end
 end
    end
end

[mag,dir]=imgradient(rgb2gray(I));
mag=round(mag./max(max(mag)),2);
posdir=dir;
s11=size(posdir,1);
s21=size(posdir,2);
for i=1:s11
    for j=1:s21
      if(posdir(i,j)<0)
          posdir(i,j)=180+posdir(i,j);
      end
    end
end


s1=size(det,1);
s2=size(det,2);
x=size(BW1,1);
y=size(BW1,2);
temp=zeros(x,y);
[L,N]=bwlabeln(BW1);
for i=1:x
    for j=1:y
        if(i<=s1 && j<=s2)
        if det(i,j) && BW1(i,j)
            rc=[0 0];
          [r,c]= find(L==L(i,j));
          rc=[r c];
          if(size(rc,1)>0)
          for k=1:size(rc,1)
              xval=rc(k,1);
              yval=rc(k,2);
              if(  round(posdir(xval,yval),0)>=60 && round(posdir(xval,yval),0)<140 | ( round(posdir(xval,yval),0) >158 && round(posdir(xval,yval),0)<163) | ( round(posdir(xval,yval),0) >170 && round(posdir(xval,yval),0)<178))
                  temp(xval,yval)=1;
              else
                  temp(xval,yval)=0;
              end
          end
          end
              
          
        end
        end
            
    end
end
temp5=temp;
for curr=6:6
    %read image from input directory
    img_route=['C:\Users\SAURADIP\Documents\MATLAB\html\GVF\ip\',str];
    temp1=imread(['C:\Users\SAURADIP\Documents\MATLAB\html\GVF\ip\',Dir_img(3).name]);
    I_orig=imread(img_route);
    I_gray=rgb2gray(I_orig);
    [m,n]=size(I_gray);
    im=I_orig;          
    %________double line cut______________

%chek if image is double line

    %For single line image
    name=Dir_img(curr+2).name;
    up=1;
%     figure,imshow(I_orig);
   %this function will perform the basic calculations and calls the major
%functions like stroke width and reconstruction
I_orig_size=I_orig;
%     I_orig=imresize(I_orig,[40 180]);
    I_gray=rgb2gray(I_orig);
    [m,n]=size(I_gray);
    %normalize gray
    if(mod(m,2)==1)
    I_gray=I_gray(1:(m-1),:);
    m=m-1;
    end
    I_gray=255-I_gray;
    I_gray=imresize(I_gray,[m n]);
    I=zeros(m,n);
    for i=1:m
        for j=1:n
            if(I_gray(i,j)<256/2)
                I(i,j)=1;
            end
        end
    end
    
    %_____________________Compute the GVF gradient and laplace____________________
    % Compute its edge map
%     disp(' Compute edge map ...');
    f = 1 - I/0.99; 
canI=double(edge(I_gray,'canny'));
% imwrite(canI,'D:\Dropbox\double line segmentation 2017\outputs\result7\c2.jpg');
 f = 1 - canI/0.99; 
    % Compute the GVF of the edge map f
%     disp(' Compute GVF ...');
    [u,v] = GVF(f, 0.2, 80); 
% Compute angles from GVF
%figure,quiver(u,v);
[r,c]=size(u);
 angle=[];
 for i=1:r
     for j=1:c
         if u(i,j)>0 & v(i,j)>0
             angle(i,j)=atand(v(i,j)/u(i,j));
         end
         if u(i,j)<0 & v(i,j)>0
             angle(i,j)=atand(v(i,j)/u(i,j));
         end
         if u(i,j)<0 & v(i,j)<0
             angle(i,j)=atand(v(i,j)/u(i,j));
         end
         if u(i,j)>0 & v(i,j)<0
             angle(i,j)=atand(v(i,j)/u(i,j));
         end
         if u(i,j)>0 & v(i,j)==0
             angle(i,j)=0;
         end
         if u(i,j)<0 & v(i,j)==0
             angle(i,j)=180;
         end
         if u(i,j)==0 & v(i,j)>0
             angle(i,j)=270;
         end
         if u(i,j)==0 & v(i,j)<0
             angle(i,j)=90;
         end
         if u(i,j)==0 & v(i,j)==0
             angle(i,j)=500;
         end
     end
 end
%     disp(' Nomalizing the GVF external force ...');
    mag = sqrt(u.*u+v.*v);
    px = u./(mag+1e-10); py = v./(mag+1e-10); 
    [fx,fy] = gradient(f); %gradient of image f
 SqrMagf = fx.*fx + fy.*fy;
 %normalize gradient with magnitude
    fx1=fx.*SqrMagf;
    fy1=fy.*SqrMagf;
    px1=px;
    
    px_paper=-px;
    py_paper=-py;
    %find laplace of image
    h=fspecial ('laplacian');
    laplaceV=filter2(h,f);
    % figure,quiver(flipud(fx1),flipud(fy1)); %gradient
ppx=flipud(px_paper);
ppy=flipud(py_paper);
   % fig1=figure,quiver(ppx,ppy,'LineWidth',0.75); %gvf
    % saveas(fig1,'C:\Users\SAURADIP\Documents\MATLAB\html\GVF\op\abc.jpg');



    plotted_image=zeros(size(canI));
% I1=input_image;
  [m,n]=size(plotted_image);
ods=plotted_image;   %ODS
ids=plotted_image;   %IDS
p3=plotted_image;   %zero vectors
p4=plotted_image; 
p5=plotted_image;  %try with py information
odr=plotted_image;   %holes
ppx=px;
ppy=py;
ods_x=px;   %ODS
ods_y=py;
ids_x=px;   %IDS
ids_y=py;
p3x=px;   %background with near 0's value
p3y=py;
p4x=px;   %IDR (Inside direction ring)
p4y=py;
odr_x=px;   %ODR (holes)
odr_y=py;
p5x=px;   %try with py information
p5y=py;

p3(abs(px)<0.004)=1;
p3(px==1)=1;
p3(px==-1)=1;
px(abs(px)<0.004)=Inf;

% % p3(abs(py)<0.006)=1;
% % py(abs(py)<0.006)=Inf;

for i=2:m-2
    for j=2:n-2
        if(temp(i,j)==1)
        %21,26=-.9, 21,27 = anything, 21,28=.9
       if ((py(i,j)>.85 && py(i+2,j)<-.85) )%|| (py(i,j)<-.85 && py(i+1,j)>.85))
%          p5(i+1,j)=1;  
        p5(i+4,j)=1;
       end
        end
    end
end

c=0;
for i=2:m-2
    for j=2:n-2
      if(temp(i,j)==1)  
% For IDR-------------i=25 j=7
a=(px(i,j)>.8 && px(i,j)<1) || (px(i+2,j)>-1 && px(i+2,j)<-.8);%><
b=px(i-1,j)>-.1 && px(i-1,j) <.1 && py(i-1,j)>-1 && py(i-1,j)<-.8;  %down
c=px(i+1,j)>-.1 && px(i+1,j)<.1 && py(i+1,j)>.8 && py(i+1,j)<1; %up
d=px(i-1,j-1)>.5 && px(i-1,j-1)<.7 && py(i-1,j-1)>-.9 && py(i-1,j-1)<-.7;%right down
e=px(i-1,j+1)<-.7 && px(i-1,j+1)>-.95 && py(i-1,j+1)>-.9 && py(i-1,j+1)<-.3;%left down
f=px(i+1,j-1)>.5 && px(i+1,j-1)<.7 && py(i+1,j-1)>.7 && py(i+1,j-1)<.95;%right up
g=px(i+1,j+1)<-.3 && px(i+1,j+1)>-.9 && py(i+1,j+1)>.7 && py(i+1,j+1)<.95;%left up
if (a)
    c=c+1;
end
if (b)
    c=c+1;
end
if (c)
    c=c+1;
end
if (d)
    c=c+1;
end
if (e)
    c=c+1;
end
if (f)
    c=c+1;
end
if (g)
    c=c+1;
end
if(c>3)
            p4(i-1,j)=1;
            p4(i+1,j)=1;
            p4(i,j+1)=1;
            p4(i-1,j+1)=1;
            p4(i+1,j+1)=1;
            p4(i-1,j-1)=1;
            p4(i,j-1)=1;
            p4(i+1,j-1)=1;
            c=0;
 
end
%i,j loop ends
      end
    end
%     c=0;
end
%imshow(p4);
swd=zeros(m,n);
ii=1;
for i=2:m-2
    for j=2:n-2
        
        if(temp(i,j)==1)
        % FOR ODR 
%         For each pixel search 8 neighbour if opposite then
        %8 neighbours of (i,j) are (i-1,j-1),(i-1,j),(i-1,j+1)
        %(i,j-1),(i,j+1)
        %(i+1,j-1),(i+1,j),(i+1,j+1)
if (px(i,j)<-.005 && px(i-1,j)<-.005 && px(i+1,j)<-.005 && px(i,j+1)>.005 &&px(i-1,j+1)>.005 && px(i+1,j+1)>.005)% && px(i+1,j-1)>.005)
    if abs(px(i,j))>.1     
%             p6(i,j)=1;
% if immediatly next pixel is white thn dont take it bn
            odr(i-1,j)=1;
            odr(i+1,j)=1;
            odr(i,j+1)=1;
            odr(i-1,j+1)=1;
            odr(i+1,j+1)=1;
            odr(i-1,j-1)=1;
            odr(i,j-1)=1;
            odr(i+1,j-1)=1;
    end
          
end   


%For ODS (p1)----------

        for a=-1:1:1
            for b=-1:1:1
                if ((a~=-1 && b~=0)||(a~=0 && b~=0)||(a~=1 && b~=0)) 
        if px(i,j)-px(i+a,j+b)<0    %<> ODS %.05 & .01
            if ((abs(px(i,j)+px(i+a,j+b))<0.2) && (abs(py(i,j)-py(i+a,j+b))>0.01))  %<> ODS
                ispositive=px(i,j)>0;
                ispositive1=px(i+a,j+b)>0;
                if ((ispositive==0 && ispositive1==1) || (ispositive==1 && ispositive1==0))
            ods(i,j)=1;
            ods(i+a,j+b)=1;
                end
            end
        end
        %for IDS (p2) which is Inside direction-----------
        if (px(i,j)-px(i+a,j+b)>=1 && px(i,j)~=Inf)  %>< IDS
            if ((abs(px(i,j)+px(i+a,j+b))<0.1))    %>< IDS .05
                ispositive=px(i,j)>0;
                ispositive1=px(i+a,j+b)>0;
                if ((ispositive==0 && ispositive1==1) || (ispositive==1 && ispositive1==0))
                if  b~=-1
                    ids(i,j)=1;
                    ids(i+a,j+b)=1;
                end
                end
            end
        end
         if (px(i,j)-px(i+a+1,j+b+1)>=1 && px(i,j)~=Inf)  %>< IDS
            if ((abs(px(i,j)+px(i+a+1,j+b+1))<0.1))    %>< IDS .05
                ispositive=px(i,j)>0;
                ispositive1=px(i+a+1,j+b+1)>0;
                if ((ispositive==0 && ispositive1==1) || (ispositive==1 && ispositive1==0))
                if  b~=-1
                    ids(i,j)=1;
                    ids(i+a+1,j+b+1)=1;
                end
                end
            end
         end
        
         %ii=ii+1;
        %-----------------------
                end                
            end
        end       
        end 
    end
end

ods_x(ods==0)=0;   ods_y(ods==0)=0;
ids_x(ids==0)=0;   ids_y(ids==0)=0; 
p3x(p3==0)=0;   p3y(p3==0)=0;
odr_x(odr==0)=0;   odr_y(odr==0)=0;
p4x(p4==0)=0;   p4y(p4==0)=0;
p5x(p5==0)=0;   p5y(p5==0)=0;
%fig1=figure,
%quiver((ppx),(ppy),'color',[0 0 1]);
 hold on


hold off    


%fig2=figure,
%imshow(input_image);
hold on

 quiver(flipud(ids_x),flipud(ids_y),'Linewidth',.75,'color',[1 0 0]);    %IDS


hold off

plotted_image=ods+ids+p3;
abc=ids_x+ids_y;

    fd=strcat(output_dir,'\');
    fn=strcat(name,'_o.jpg');
    fileName=strcat(fd,fn);


end


     
% METHOD 1 FROM GITHUB 
%mag=round(mag,4)/max(max(round(mag,4)));
%mag=gray;
%[magn,dir]=imgradient(rgb2gray(I));

%mag = imtranslate(mag,[-3, -3]);
%dir=imtranslate(dir,[-3, -3]);
%mag=imcrop(mag,[0 0 size(mag,2)-3 size(mag,1)-3]);
%dir=imcrop(dir,[0 0 size(mag,2) size(mag,1)]);
 [edgePointRows, edgePointCols] = find((temp));
%edgePointRows(1,1)= 57;
%edgePointCols(1,1)=109;
% Getting size of the representative image
[m,n] = size(I_gray);

% Initializing Stoke Width array with infinity
swtMap = zeros(m,n);
for i=1:m
    for j=1:n
        swtMap(i,j) = inf;
    end
end

% Set the maximum stroke width, this number is variable for now but must be
% made to be more dynamic in the future
maxStrokeWidth = 350;

% Initialize container for all stoke points found
strokePointsX = zeros(size(edgePointCols));
strokePointsY = zeros(size(strokePointsX));

sizeOfStrokePoints = 0;
gvf=zeros(size(angle,1),size(angle,2));
dist=zeros(size(angle,1),size(angle,2));
 ngb=zeros(size(gvf,1),size(gvf,2)); 
% Iterate through all edge points and compute stoke widths
for i=1:size(edgePointRows)
    step = 1;
    new=0;
    initialX = edgePointRows(i);
    initialY = edgePointCols(i);
    isStroke = 0;
    neighX=zeros(9,1);
    neighY=zeros(size(neighX));
    if(initialX<size(dir,1) && initialY<size(dir,2))
    initialTheta = dir(initialX,initialY);
    initialSeed=I_gray(initialX,initialY);
    end
    sizeOfRay = 0;
    pointOfRayX = zeros(maxStrokeWidth,1);
    pointOfRayY = zeros(maxStrokeWidth,1);
    
    % Record first point of the ray
    pointOfRayX(sizeOfRay+1) = initialX;
    pointOfRayY(sizeOfRay+1) = initialY;
    
    % Increase the size of the ray
    sizeOfRay = sizeOfRay + 1;
    
    nextX=initialX;
    nextY=initialY;
    flag=0;
    % Follow the ray
    while step < 5  && flag==0 && BW2(initialX,initialY)==1
       % gvf=zeros(size(angle,1),size(angle,2));
        nextX = (round(initialX + cosd(initialTheta) * 1 * step));
        nextY = (round(initialY + sind(initialTheta) * 1 * step));
        
      
        step = step + 1;
   
       
        
        % Break loop if out of bounds.  For some reason this is really
        % slow.
        if nextX < 1 | nextY < 1 | nextX > m | nextY > n 
            break
        else
%             prev=mag(nextX,nextY);
%          curr=mag(nextXX,nextYY);
        end
          
        % Record next point of the ray
        pointOfRayX(sizeOfRay+1) = nextX;
        pointOfRayY(sizeOfRay+1) = nextY;
        
        % Increase size of the ray
        sizeOfRay = sizeOfRay + 1;
        
        % Another edge pixel has been found value  = seed pixel then i
        % reach a stroke width pair
      %  if abs(curr-prev)<=1.0  | curr > prev
            
          %  oppositeTheta = theta(nextX,nextY);
            
            % Gradient direction roughtly opposite
%             if abs(initialSeed-curr)<=1.0
%                 isStroke = 1;
%                 strokePointsX(sizeOfStrokePoints+1) = initialX;
%                 strokePointsY(sizeOfStrokePoints+1) = initialY;
%                 sizeOfStrokePoints = sizeOfStrokePoints + 1;
%                 break;
%             end
            
         %   break
         
         if((round(angle(nextX,nextY))+round(angle(initialX,initialY)))<=10 )
             % got a pixel pair with opposite and equal angle, now job is
             % to find neighbours
             gvf(nextX,nextY)=1;
             gvf(initialX,initialY)=1;
             swd=sizeOfRay;
%              for iii=1:sizeOfRay
%                  dist(pointOfRayX(iii),pointOfRayY(iii))=swd;
%              end
            % swd=sizeOfRay;
%              dist(initialX,initialY)=swd;
%              dist(nextX,nextY)=swd;
             [gvf]=calc_neigh(gvf,initialX,initialY,initialTheta,angle,px,I_gray,bod,swd);

            else
            gvf(nextX,nextY)=0;
             gvf(initialX,initialY)=0;
         end
     
    end
    
%     
%            mid_gvf=gvf;
%      % if flag = 1 means i found stroke width
%      neigh_x=zeros(9,1);
%      neigh_y=zeros(9,1);
%      neigh_loop=1;
%     while( neigh_loop==1 )
%        % while (neigh_loop==1)
%             x=0;
%             y=0;
%             k=1;
%             if initialX>=2 && initialY>=2 &&initialX<ms && initialY<N
%                 % calculating 8 neighbour pixel of initial pixel
%         for x=initialX-1:initialX+1
%             for y=initialY-1:initialY+1
%                 if((x~=initialX && y ~=initialY )|| (x ~= nextX && y~=nextY))
%               neigh_x(k)=x;  
%               neigh_y(k)=y;
%               k=k+1;
%                 end
%             end
%         end
%             end
%         
%         for m=1:size(nonzeros(neigh_x),1)
%             initialTheta=angle(neigh_x(m),neigh_y(m));
%             initialX=neigh_x(m);
%             initialY=neigh_y(m);
%             %step2=1;
%             %flag2=0;
%             sizeOfRay2=1;
%             %pointOfRayX1 = zeros(maxStrokeWidth,1);
%             %pointOfRayY1 = zeros(maxStrokeWidth,1);
%             [flag2,gvf1,initialX,initialY]=swt_neigh(initialX,initialY,initialTheta,angle,gvf,m,n,I_gray,1);  % find new swd for neighbour
%             gvf=gvf1+gvf;
%            if(flag2==1)
%                neigh_loop=1; % flag set to find 8 neighbours i.e successful swd found
%                break;
%            %else
%               % neigh_loop=0; % flag set not to find 8 neighbours i.e unsuccessful
%                %flag=0;
%               % break;
%            end
%         end
%         
%     end   
%     %end 
% %     mid_gvf=gvf;
% %      % if flag = 1 means i found stroke width
% %      neigh_x=zeros(9,1);
% %      neigh_y=zeros(9,1);
% %      neigh_loop=1;
% %     while(flag==1 && neigh_loop==1 | flag2==1  )
% %        % while (neigh_loop==1)
% %             x=0;
% %             y=0;
% %             k=1;
% %             if initialX>=2 && initialY>=2 &&initialX<ms && initialY<N
% %                 % calculating 8 neighbour pixel of initial pixel
% %         for x=initialX-1:initialX+1
% %             for y=initialY-1:initialY+1
% %                 if((x~=initialX && y ~=initialY )|| (x ~= nextX && y~=nextY))
% %               neigh_x(k)=x;  
% %               neigh_y(k)=y;
% %               k=k+1;
% %                 end
% %             end
% %         end
% %             end
% %         
% %         for m=1:size(nonzeros(neigh_x),1)
% %             initialTheta=angle(neigh_x(m),neigh_y(m));
% %             initialX=neigh_x(m);
% %             initialY=neigh_y(m);
% %             %step2=1;
% %             %flag2=0;
% %             sizeOfRay2=1;
% %             %pointOfRayX1 = zeros(maxStrokeWidth,1);
% %             %pointOfRayY1 = zeros(maxStrokeWidth,1);
% %             [flag2,gvf1]=swt_neigh(initialX,initialY,initialTheta,angle,gvf,m,n,I_gray,1);  % find new swd for neighbour
% %             gvf=gvf1+gvf;
% %            if(flag2==1)
% %                neigh_loop=1; % flag set to find 8 neighbours i.e successful swd found
% %                break;
% %            else
% %                neigh_loop=0; % flag set not to find 8 neighbours i.e unsuccessful
% %                flag=0;
% %                break;
% %            end
% %         end
% %         
% %     end   
% %     %end

end

mat=temp;
size1=size(gvf,1);
size2=size(gvf,2);
size11=size(mat,1);
size12=size(mat,2);
det_temp=zeros(size1*size2,1);
if(size(mat,1)<size(gvf,1))
    % smaller matrice converted to larger matrice
det_crop = padarray(mat,(size(gvf)-size(mat)),'post');
gvf=(det_crop)+gvf;
elseif(size(mat,1)>size(gvf,1))
   % larger matrice cropped to smaller
   det_crop = imcrop(mat,[0 0 size2 size1]);
gvf=(det_crop)+gvf;
else
    gvf=gvf+mat;
end


for xs=1:size(gvf,1)
for ys=1:size(gvf,2)
if(gvf(xs,ys)==2 || gvf(xs,ys)==1 )
gvf(xs,ys)=1;
else
gvf(xs,ys)=0;
end
end
end
%imshow(imfuse(gvf,I_gray))
imshow(gvf);

s1=size(det,1);
s2=size(det,2);
x=size(BW2,1);
y=size(BW2,2);
%temp1=zeros(x,y);
[L1,N]=bwlabeln(BW2);
for i=1:x
    for j=1:y
        if(i<=s1 && j<=s2)
        if gvf(i,j) && BW2(i,j)
            rc=[0 0];
          [r,c]= find(L1==L1(i,j));
          rc=[r c];
          if(size(rc,1)>0)
          for k=1:size(rc,1)
              xval=rc(k,1);
              yval=rc(k,2);
%               if(  round(posdir(xval,yval),0)>=60 && round(posdir(xval,yval),0)<140 | ( round(posdir(xval,yval),0) >158 && round(posdir(xval,yval),0)<163) | ( round(posdir(xval,yval),0) >170 && round(posdir(xval,yval),0)<178))
%                   temp(xval,yval)=1;
%               else
%                   temp(xval,yval)=0;
%               end
                gvf(xval,yval)=1;
          end
          end
              
          
        end
        end
            
    end
end

CC=bwconncomp(gvf,8);
[L2,N2]=bwlabeln(gvf);
stats = regionprops(CC,'BoundingBox');
k1=0;
img=gvf;
% while(k1>=k)
% [img,k1]=restore(img,BW2,L1,s1,s2);
% 
% end
% when restore() does not do well try restore_test()
% [lb_ang3]=restore_test(gvf,BW2,L1,s1,s2);
% [lb_ang4]=restore_test(lb_ang3,BW2,L1,s1,s2);
% [lb_ang5]=restore_test(lb_ang4,BW2,L1,s1,s2);
% [lb_ang6]=restore_test(lb_ang5,BW2,L1,s1,s2);
% 
[lb_ang3]=restore(gvf,BW2,L1,s1,s2);
[lb_ang4]=restore(lb_ang3,BW2,L1,s1,s2);
[lb_ang5]=restore(lb_ang4,BW2,L1,s1,s2);
[lb_ang6]=restore(lb_ang5,BW2,L1,s1,s2);
%imshow(lb_ang4)
imshow(bod);
hold on;
fill=lb_ang4;
 exp=lb_ang5;
% for straight lines or disconnected we can simply pass components as it is
% , for character component combination weuse below modification


 fill=imfill(fill,'holes');
 fill=bwmorph(fill,'thicken');
 fill=bwmorph(fill,'majority');
 fill=imfill(fill,'holes');
 fill=bwmorph(fill,'thicken',2);
 fill=bwmorph(fill,'majority');
 
CC_fn=bwconncomp(exp,8);


stats_fn = regionprops(CC_fn,'BoundingBox');

for m_i=1:size(stats_fn)
thisBB = stats_fn(m_i).BoundingBox;
if (  (thisBB(4)<2*thisBB(3) || thisBB(3)>3*thisBB(4) ) && (thisBB(3)*thisBB(4)>=30)  )
rectangle('Position', [thisBB(1),thisBB(2),thisBB(3),thisBB(4)],...
  'EdgeColor','r','LineWidth',2 )
end
end
