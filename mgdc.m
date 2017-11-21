function [MG_dir,mag,pic,fuse,DC]=mgdc(bod)

%bod=imread('torsotext1.jpg');
resize=bod;           
B=rgb2gray(resize);
[Xmag,Xdir]=imgradient(B);         
Xmag=double(Xmag)/double(max(max(Xmag)));
C=double(B)/double(max(max(B)));
 [cr,cc]=size(C);
%DC=zeros(cr+3,cc+3);
D=double(C);


DC=D;

win=3;  % sobel window size , generaly 3 

 
switch win
    case 3
       Sx=[-1 -2 -1 ; 0 0 0 ; 1 2 1]; %for 3 x 3
       Sy=[-1 0 1 ; -2 0 2 ; -1 0 1];
    case 7 
    Sx = [ 3 2 1   0 -1 -2 -3
           4 3  2  0 -2 -3 -4 
           5 4  3  0 -3 -4 -5 
           6 5  4  0 -4 -5 -6 
           5 4  3  0 -3 -4 -5 
           4 3  2  0 -2 -3 -4
           3 2  1  0 -1 -2 -3 ];
       Sy=transpose(Sx);
    case 5
       Sx=[-1  -4  -6  -4  -1 ; -2  -8  -12  -8  -2 ; 0  0  0  0  0  ; 2  8  12  8  2 ; 1  4  6  4  1 ];
       Sy=[1  2   0   -2  -1 ; 4   8   0   -8  -4 ; 6   12   0   -12  -6 ; 4   8   0   -8  -4 ; 1   2   0   -2  -1];
end 
        [r,c]=size(D);
        
        O=zeros(r,c);
        O1=zeros(r,c);
       
        f=zeros(r-win,c-win);
        m=zeros(r-win,c-win);
         p=zeros(r-win,c-win);
        px=zeros(win,win);
         py=zeros(win,win);
         xyz=zeros(r-win,c-win);
         xy=zeros(r-win,c-win);
          pic=zeros(r-win,c-win);
          fuse=zeros(r-win,c-win);
          MG_dir=zeros(r-win,c-win);
          DC_dir=zeros(r-win,c-win);
          dc=xyz;
          dc_lam1=zeros(r-win,c-win);
          dc_lam2=zeros(r-win,c-win);
          mg_division=DC_dir;
          dc_mag=DC_dir;
          DC_quad=DC_dir;
          MG_quad=MG_dir;
          mag=m;
          %means=mag; % for which i need to find k means cluster
          if(r>0 )
            for k=1:(r-win)
            for l=1:(c-win)
        patch=D(k:(k+(win-1)),l:(l+(win-1)));
        gx=Sx.*patch;
        gy=Sy.*patch;
%         if(max(max(gx))~=0 && max(max(gy))~=0)
%         gx=gx./max(max(gx));
%         gy=gy./max(max(gy));
%         end
       Gx=sum(sum(gx));
       Gy=sum(sum(gy));
%         O(k,l)=sum(sum(gx));
%         O1(k,l)=sum(sum(gy));
% %             for i=k:(k+2)
% %             for j=l:(l+2)
% %       O(i,j)=sum(sum(gx));
%        % O1(i,j)=sum(sum(gy));
% %             end
% %             end
%         sums=0.0;
%        
%        
%   
%        %m(k,l)=log(sum(f(:))/double(9));
%         
%             end
%             end
%  %f=power(sqrt((power(O,2)+power(O1,2))/2),4);
%  if(max(max(O))~=0 && max(max(O1))~=0)
%  Ix_binary=round(double(O)/double(max(max(O))),1);
%  Iy_binary=round(double(O1)/double(max(max(O1))),1);
%  f=power(sqrt((power(Ix_binary,2)+power(Iy_binary,2))/2),4);
%  end

 O(k,l)=sum(sum(gx));
        O1(k,l)=sum(sum(gy));

%             for i=k:(k+2)
%             for j=l:(l+2)
%       O(i,j)=sum(sum(gx));
       % O1(i,j)=sum(sum(gy));
%             end
%             end
        sums=0.0;
       
       
  
       %m(k,l)=log(sum(f(:))/double(9));
        
            end
            end
 %f=power(sqrt((power(O,2)+power(O1,2))/2),4);
%  if(max(max(O))~=0 && max(max(O1))~=0)
%  Ix_binary=abs(round(double(O)/double(max(max(O))),1));
%  Iy_binary=abs(round(double(O1)/double(max(max(O1))),1));
%  f=sqrt((power(Ix_binary,2)+power(Iy_binary,2))/2);
%  end


Ix_binary=(O);
Iy_binary=(O1);
 f=sqrt((power(Ix_binary,2)+power(Iy_binary,2))/2);
 
   for kk=1:(r-win)
       for ll=1:(c-win)
        G=f(kk:(kk+(win-1)),ll:(ll+(win-1)));
        if ( (sum(sum(G))/double(win*win))~=0)
            if(log(sum(sum(G))/double(win*win)) > 0 )
               m(kk,ll)=((log1p(sum(sum(G))/double(win*win))));   % Magnitude Gradient of Image
            else
            m(kk,ll)= abs((log1p(sum(sum(G))/double(win*win))));
            end
        end
       % if(Iy_binary(kk,ll)~=0)
        mg_division(kk,ll)=round((Ix_binary(kk,ll)./Iy_binary(kk,ll)),2); % division of Iy / Ix for each pixel of image
        
       % end
px=(Ix_binary(kk:(kk+(win-1)),ll:(ll+(win-1))));
py=(Iy_binary(kk:(kk+(win-1)),ll:(ll+(win-1))));
      T=cov(px,py);
      Ei=eig(T);
      dc_lam1(kk,ll)=Ei(1,1);
      dc_lam2(kk,ll)=Ei(2,1);
      
       x1=Ei(1,1);
     y1=Ei(2,1);
    %DC_dir(kk,ll) =atan(x1/y1); % direction of each edge point for DC
%         if(x1<0)
%          x1=0;
%         end
%         if(y1<0)
%          y1=0;
%         end
     su=abs(x1-y1);
     ad=abs(x1+y1);
     if( dc_lam2(kk,ll)~=0)
         dc_mag(kk,ll)=round((dc_lam1(kk,ll)./dc_lam2(kk,ll)),1); 
     end
     if(ad>0 && su < ad )
     dc(kk,ll)=power(round(su/ad,2),2);
     end
     %xyz(kk,ll)=round(su.*su,4);
     %xy(kk,ll)=round(ad.*ad,4);
      if(ad > 0  && (dc(kk,ll)>0))
     pic(kk,ll)=abs(dc(kk,ll));  % Directional Coherence for the image
      else
          pic(kk,ll)=0;
       
      end
      
      
       end
   end
   
   mg_division_normalised = round(double(mg_division)/double(max(max(mg_division))),4);  % division of Iy / Ix for each pixel of image and normalised within 0 to 1 
        MG_dir =atand(mg_division);  % direction of each edge point for MG by tan-1(Iy/Ix)
   
   %dc_mag=round(abs(dc_lam1./dc_lam2),1);
   dc_mag_normalised = round(double(dc_mag)/double(max(max(dc_mag))),1);  % normalised within 0 and 1 result of lambda1/lambda2 rounded to 1 decimal place 
   DC_dir = atand(dc_mag);  % direction of DC values in DC matrix tan-1(lambda1/lambda2)
   
   mag=double(m)/double(max(max(m)));        
     %dc=double(pic)/double(max(max(pic)));
     [rr,cc]=size(m);
     
     % FOR 5 X 5 Window
%          for is=1:(rr-5)
%             for js=1:(cc-5)
%         matmg = m(is:(is+4),js:(js+4));
%         matdc = pic(is:(is+4),js:(js+4));
%         maxmg= max(matmg(:));
%         minmg= min(matmg(:));
%          maxdc= max(matdc(:));
%         mindc= min(matdc(:));
%         diffmg=maxmg - minmg;
%         diffdc=maxdc - mindc;
% %         if(diffmg > diffdc)
% %             result=diffmg;
% %         else
% %             result=diffdc;
% %         end
%         %fuse(is,js)=result;
% %                         if ((diffmg+diffdc)<0)
% %              fuse(is,js)=0;
% %          
% %                         else
%              fuse(is,js)=diffmg+diffdc;
%              newfuse=double(fuse)/double(max(max(fuse)));
%                       %  end            
%            end
%          end
         mg1=imresize(imgradient(im2bw(bod)),.25);
      fuse=(mag+pic);
      
      %fuse=double(fuse)/double(max(max(fuse)));         
          end
end