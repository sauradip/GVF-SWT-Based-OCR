function[ lb_ang,k ]=restore(gvf,BW2,L1,s1,s2)


hold on;
label=0;
lb_ang=zeros(s1+3,s2+3); % matrix angle component labelled
flags=0;flags1=0;
in_angle=Inf;
% finding connected component on image 
CC=bwconncomp(gvf,8);
[L2,N2]=bwlabeln(gvf);
stats = regionprops(CC,'BoundingBox');
orient= regionprops(CC,'Orientation');
CC3=bwconncomp(BW2,8);
orient_cann=regionprops(CC3,'Orientation');

% bounding box for components
for k = 1 : size(stats)
  thisBB = stats(k).BoundingBox;


if((thisBB(4)<2*thisBB(3) || thisBB(3)>3*thisBB(4) ) && thisBB(3)*thisBB(4)>=10 && thisBB(3)*thisBB(4)<=0.1 * size(gvf,1)* size(gvf,2) && thisBB(3)<=0.2*s1 && thisBB(4)<=0.2*s2)
  rectangle('Position', [thisBB(1),thisBB(2),thisBB(3),thisBB(4)],...
  'EdgeColor','r','LineWidth',2 )

comp_angle=orient(k).Orientation;
% disp(comp_angle) % component orientation
% if(comp_angle<0)
% comp_angle=90+comp_angle;
% end

steps=1;
if comp_angle < 0 
in_angle=180 + (round(comp_angle));
else
in_angle=abs(round(comp_angle));
 end% initial component angle
init_x=0;
init_y=0;
flag_com=0;

% 
% % finding the component's x , y axis within each bounding box 
y_coord=round(thisBB(1));
y_coord1=round(thisBB(1))+thisBB(3);
x_coord=round(thisBB(2));
x_coord1=round(thisBB(2))+thisBB(4);
if  x_coord1<=size(BW2,1) && y_coord1<=size(BW2,2) && x_coord>0 && y_coord > 0
label=mode(nonzeros(L2(x_coord:x_coord1,y_coord:y_coord1)));
end
rc3=[0 0];
          [r3,c3]= find(L2==label);
          rc3=[r3 c3];
          if(size(rc3,1)>0)
          for k3=1:size(rc3,1)
              xval3=rc3(k3,1);
              yval3=rc3(k3,2);
% if( abs(90-in_angle)<=50 )
%                 lb_ang(xval3,yval3)=1;
% 
% else
% lb_ang(xval3,yval3)=0;
% 
% end




          end
          end

% random x, y 

init_x=rc3(1,1);
init_y=rc3(1,2);

nxt_x=1;nxt_y=1;
nxt_x1=1;nxt_y1=1;
flags=0;flags1=0;
while(steps<=8 && nxt_x<s1 && nxt_y<s2 && nxt_x>0 && nxt_y>0 && nxt_x1>0 && nxt_y1>0 && nxt_x1<s1 && nxt_y1 <s2)
% moving on left of component 
nxt_x=round(init_x + cosd(in_angle) * 1 * steps);
nxt_y=round(init_y + sind(in_angle) * 1 * steps);
% moving on right on component
nxt_x1=round(init_x - cosd(in_angle) * 1 * steps);
nxt_y1=round(init_y - sind(in_angle) * 1 * steps);
steps=steps+1;
if ( nxt_x>0 && nxt_y >0 && nxt_x1>0 && nxt_y1>0 && nxt_x<s1 && nxt_y<s2 && nxt_x1<s1 && nxt_y1<s2 )

 if((L2(nxt_x,nxt_y)~=k && BW2(nxt_x,nxt_y)==1 )|| (L2(nxt_x1,nxt_y1)~=k && BW2(nxt_x1,nxt_y1)==1 ) )
% lb_ang(nxt_x,nxt_y)=1;
%lb_ang(nxt_x1,nxt_y1)=1;
% for left component angle
if( BW2(nxt_x,nxt_y)~=0  & L1(nxt_x,nxt_y)~=0)
temp_angle1=abs(round(orient_cann(L1(nxt_x,nxt_y)).Orientation));
if(temp_angle1<0 )
temp_angle1=180+temp_angle1;
else
temp_angle1=abs(temp_angle1);
disp('init')
disp(in_angle)
disp('angle1')
disp(temp_angle1)
end
if( abs(in_angle - temp_angle1)<=110)
%rc1=[0 0];
lb_ang(nxt_x,nxt_y)=1;


          [r1,c1]= find(L1==L1(nxt_x,nxt_y));
          rc1=[r1 c1];
          if(size(rc1,1)>0)
          for k1=1:size(rc1,1)
              xval1=rc1(k1,1);
              yval1=rc1(k1,2);

                lb_ang(xval1,yval1)=1;
          end
          end
flags=1;
% rc3=[0 0];
%           [r5,c5]= find(L2==label);
%           rc5=[r5 c5];
%           if(size(rc5,1)>0)
%           for k5=1:size(rc5,1)
%               xval5=rc5(k5,1);
%               yval5=rc5(k5,2);
% 
%                 lb_ang(xval5,yval5)=1;
% 
% 
% 
% 
%           end
%           end

end

%           [r5,c5]= find(L2==label);
%           rc5=[r5 c5];
%           if(size(rc5,1)>0)
%           for k5=1:size(rc5,1)
%               xval5=rc5(k5,1);
%               yval5=rc5(k5,2);
% 
%                 lb_ang(xval5,yval5)=0;
% 
% 
% 
% 
%           end
%           end
end
% for right component angle
if( BW2(nxt_x1,nxt_y1)~=0 &  L1(nxt_x1,nxt_y1)~=0  )
temp_angle2=abs(round(orient_cann(L1(nxt_x1,nxt_y1)).Orientation));
if(temp_angle2<0 )
temp_angle2=180+temp_angle2;
else
temp_angle2=abs(temp_angle2);
end
disp('init')
disp(in_angle)
disp('angle2')
disp(temp_angle2)
if( abs(in_angle - temp_angle2)<=110  )
rc2=[0 0];


          [r2,c2]= find(L1==L1(nxt_x1,nxt_y1));
          rc2=[r2 c2];
          if(size(rc2,1)>0)
          for k1=1:size(rc2,1)
              xval1=rc2(k1,1);
              yval1=rc2(k1,2);

                lb_ang(xval1,yval1)=1;
          end
          end
flags=1;
% rc3=[0 0];
%           [r4,c4]= find(L2==label);
%           rc4=[r4 c4];
%           if(size(rc4,1)>0)
%           for k4=1:size(rc4,1)
%               xval4=rc4(k4,1);
%               yval4=rc4(k4,2);
% 
%                 lb_ang(xval4,yval4)=1;
% 
% 
% 
% 
%           end
%           end

end
         


end
% if(flags1==1 || flags==1)
%  [r4,c4]= find(L2==label);
%           rc4=[r4 c4];
%           if(size(rc4,1)>0)
%           for k4=1:size(rc4,1)
%               xval4=rc4(k4,1);
%               yval4=rc4(k4,2);
% 
%                 lb_ang(xval4,yval4)=1;
% 
% 
% 
% 
%           end
%           end
% end
rc3=[0 0];
          [r3,c3]= find(L2==label);
          rc3=[r3 c3];
          if(size(rc3,1)>0)
          for k3=1:size(rc3,1)
              xval3=rc3(k3,1);
              yval3=rc3(k3,2);
if( (flags==1 || flags1==1 ) && in_angle <90 && in_angle>70 )
                lb_ang(xval3,yval3)=1;
elseif((flags==1 & flags1==1 ) && in_angle <90 && in_angle>70)
lb_ang(xval3,yval3)=1;
else
lb_ang(xval3,yval3)=0;

end




          end
          end
end
end
% for com_x=x_coord:x_coord1
% for com_y=y_coord:y_coord1
% if  com_x<=size(BW2,1) && com_y<=size(BW2,2)
% if( flag_com==0 && BW2(com_x,com_y)~=0)
% init_x=com_x;
% init_y=com_y;
% flag_com=1;
% end
% end
% 
% end
% 
% end
% nxt_x=1;nxt_y=1;
% nxt_x1=1;nxt_y1=1;
% while(steps<10 && nxt_x<s1 && nxt_y<s2 && nxt_x>0 && nxt_y>0 && nxt_x1>0 && nxt_y1>0 && nxt_x1<s1 && nxt_y1 <s2)
% % 
% % moving on right of component 
% nxt_x=round(init_x + cosd(in_angle+90) * 1 * steps);
% nxt_y=round(init_y + sind(in_angle+90) * 1 * steps);
% % moving on left on component
% nxt_x1=round(init_x + cosd(in_angle-90) * 1 * steps);
% nxt_y1=round(init_y + sind(in_angle-90) * 1 * steps);
% steps=steps+1;
% % lb_ang(nxt_x,nxt_y)=1;
% % lb_ang(nxt_x1,nxt_y1)=1;
% 
% % % i am extending the inital x , y along the line of orientation and if i
% % % found same component number i will skip and i will only accept different
% % % component
% if ( nxt_x>0 && nxt_y >0 && nxt_x1>0 && nxt_y1>0 && nxt_x<s1 && nxt_y<s2 && nxt_x1<s1 && nxt_y1<s2 )
% 
% if((L1(nxt_x,nxt_y)~=k && BW2(nxt_x,nxt_y)==1 )|| (L1(nxt_x1,nxt_y1)~=k && BW2(nxt_x1,nxt_y1)==1 ) )
% disp('ok')
% % i entered this loop means i got a pixel of a character component now
% % trace it 
% 
% % towards right side component
% if( BW2(nxt_x1,nxt_y1)~=0 & BW2(nxt_x1,nxt_y1)~=0 & L1(nxt_x,nxt_y)~=0)
% disp('x')
% if ( in_angle<=abs(orient(L1(nxt_x,nxt_y)).Orientation)+5 || in_angle>=abs(orient(L1(nxt_x,nxt_y)).Orientation)-5 )
%  rc1=[0 0];
%           [r1,c1]= find(L1==L1(nxt_x,nxt_y));
%           rc1=[r1 c1];
%           if(size(rc1,1)>0)
%           for k1=1:size(rc1,1)
%               xval1=rc1(k1,1);
%               yval1=rc1(k1,2);
% 
%                 lb_ang(xval1,yval1)=1;
%           end
%           end
% end
% end
% % towards left of component 
% if( BW2(nxt_x1,nxt_y1)~=0 & BW2(nxt_x1,nxt_y1) ~=0 & L1(nxt_x1,nxt_y1)~=0)
% if( in_angle<=abs(orient(L1(nxt_x1,nxt_y1)).Orientation)+5 || in_angle>=abs(orient(L1(nxt_x1,nxt_y1)).Orientation)-5 )
% disp('x1')
% rc2=[0 0];
%           [r2,c2]= find(L1==L1(nxt_x1,nxt_y1));
%           rc2=[r2 c2];
%           if(size(rc2,1)>0)
%           for k2=1:size(rc2,1)
%               xval2=rc1(k2,1);
%               yval2=rc1(k2,2);
% 
%                 lb_ang(xval2,yval2)=1;
%           end
%           end
% end
% end
% end
% end
% 
% end


 end
 end
if(label==0)
rc3=[0 0];
          [r3,c3]= find(L2==label);
          rc3=[r3 c3];
          if(size(rc3,1)>0)
          for k3=1:size(rc3,1)
              xval3=rc3(k3,1);
              yval3=rc3(k3,2);



lb_ang(xval3,yval3)=0;

           end



          end
end
end


end