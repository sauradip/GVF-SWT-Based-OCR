function [flag2,gvf,initialX,initialY] = swt_neigh(initialX,initialY,initialTheta,angle,px,gvf,m,n,I_gray,swd,bod,isnei)

disp('Stroke Width Started');
step2=1;flag2=0;
sizeOfRay2=1;
 pointOfRayX1 = zeros(350,1);
            pointOfRayY1 = zeros(350,1);
[mag,dir]=imgradient(rgb2gray(bod));
while step2<5 && flag2==0   % loop to calculate stroke width of neighbour and so on
               nextX = (round(initialX + cosd(initialTheta) * 1 * step2));
        nextY = (round(initialY + sind(initialTheta) * 1 * step2));
        
      
        step2 = step2 + 1;
   
       
        
        % Break loop if out of bounds.  For some reason this is really
        % slow.
        if nextX < 1 | nextY < 1 | nextX > m | nextY > n 
            break
        else
%             prev=mag(nextX,nextY);
%          curr=mag(nextXX,nextYY);
        end
          
        % Record next point of the ray
        pointOfRayX1(sizeOfRay2+1) = nextX;
        pointOfRayY1(sizeOfRay2+1) = nextY;
%         disp(angle(nextX,nextY));
%         disp(nextY);
        % Increase size of the ray
        sizeOfRay2 = sizeOfRay2 + 1;

%disp(sizeOfRay2);
       % gvf(nextX,nextY)=1;
        %     gvf(initialX,initialY)=1;
        if(isnei==1)
          %  disp(round(px(nextX,nextY))+round(px(initialX,initialY)));
        if((abs(angle(initialX,initialY)+angle(nextX,nextY))<=10)  && sizeOfRay2<=swd && abs(I_gray(nextX,nextY)-I_gray(initialX,initialY))<=10 && I_gray(initialX,initialY)>230 && I_gray(initialX,initialY) < 245 )
             gvf(nextX,nextY)=1;
             gvf(initialX,initialY)=1;
%disp('yes');
            % disp(round(px(nextX,nextY))+round(px(initialX,initialY)));
              
%              for iii=1:sizeOfRay2
%                 if(pointOfRayX1(iii)>2 & pointOfRayY1(iii)>2)
%                  gvf(pointOfRayX1(iii),pointOfRayY1(iii))=1;
%                 end
%              end
              
             flag2=1; % means one of the neighbour of initial seed has got a stroke width , so we must find 8 neighbour
        else
            gvf(nextX,nextY)=0;
             gvf(initialX,initialY)=0;
        end
        else
           if((round(angle(nextX,nextY))+round(angle(initialX,initialY)))==0 )
             gvf(nextX,nextY)=1;
             gvf(initialX,initialY)=1;
%              for iii=1:sizzeOfRay2
%                  gvf(pointOfRayX1(iii),pointOfRayY1(iii))=1;
%              end
            % disp('yes');
             flag2=1; % means one of the neighbour of initial seed has got a stroke width , so we must find 8 neighbour
          end 
        end
end
disp('Stroke Width Finished');
end