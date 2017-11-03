function[gvf]=calc_neigh(gvf,initialX,initialY,initialTheta,angle,px,I_gray,bod,swd)
step2=1;
[mag,dir]=imgradient(rgb2gray(bod));
nextX = (round(initialX + cosd(initialTheta) * 1 * step2));
        nextY = (round(initialY + sind(initialTheta) * 1 * step2));
   ngb=zeros(size(gvf,1),size(gvf,2)); 
   ms=size(gvf,1);
   N=size(gvf,2);
           mid_gvf=gvf;
           neigh_loop=1;
           gvf1=0;
     % if flag = 1 means i found stroke width
     neigh_x=zeros(9,1);
     neigh_y=zeros(9,1);
%       for xx=1:size(gvf,1)
%           for yy=1:size(gvf,2)
%   %   for xx=86:86
%   %       for yy=26:26
%              if(gvf(xx,yy)==1)
%      neigh_loop=1;
%      initialX=xx;
%      initialY=yy;
%              end
disp(initialX);
disp(initialY);
    while( neigh_loop==1 )
       % while (neigh_loop==1)
            x=0;
            y=0;
            k=1;
            if initialX>=2 && initialY>=2 &&initialX<ms && initialY<N
                % calculating 8 neighbour pixel of initial pixel
        for x=initialX-1:initialX+1
            for y=initialY-1:initialY+1
                if((x~=initialX && y ~=initialY ) && (x~=nextX && y~=nextY))
              neigh_x(k)=x;  
              neigh_y(k)=y;
              k=k+1;
              
                end
            end
        end
            end
        
        for m=1:size(nonzeros(neigh_x),1)
            initialTheta=dir(neigh_x(m),neigh_y(m));
            initialX=neigh_x(m);
            initialY=neigh_y(m);
            %step2=1;
            %flag2=0;
           % sizeOfRay2=1;
            %pointOfRayX1 = zeros(maxStrokeWidth,1);
            %pointOfRayY1 = zeros(maxStrokeWidth,1);
            swd1=swd;
            [flag2,gvf,initialX,initialY]=swt_neigh(initialX,initialY,initialTheta,angle,px,gvf,ms,N,I_gray,swd1,bod,1);  % find new swd for neighbour
            disp(flag2);
           if(flag2==1)
               neigh_loop=1; % flag set to find 8 neighbours i.e successful swd found
               break;
           else
               neigh_loop=0; % flag set not to find 8 neighbours i.e unsuccessful
              flag=0;
               break;
           end
        end
        
    end   
%          end
%      end



end