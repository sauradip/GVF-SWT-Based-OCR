function [df1,df2,select,centroid] = kmeans2(G)

% % 
 means=G; % for which i need to find k means cluster
nrows = size(means,1);
ncols = size(means,2);
ab = reshape(means,nrows*ncols,1);
[idx3,centroid] = kmeans(ab,2);


% Average / Max  / Min cluster will depend on the centroid value , higher
% centroid means max , lower centroid means min and rest is avg

cluster1=ab(idx3==1);
cluster2=ab(idx3==2);
%cluster3=ab(idx3==3);
[r1,c1]=size(cluster1); 
[r2,c2]=size(cluster2); 
%[r3,c3]=size(cluster3);



[rf,cf]=size(means);
df3=zeros(rf,cf);
df2=zeros(rf,cf);
df1=zeros(rf,cf);
df4=zeros(rf,cf);

%converting the clustered array into real clustered matrix of m x n size
 
 for ccol=1:r1
     % process 1 %
        value=round(cluster1(ccol,1),2);
       df1(round(means,2)==value)=1;        %max cluster
 end
      
 for ccol1=1:r2
       value1=round(cluster2(ccol1,1),2);
       df2(round(means,2)==value1)=1;       %avg cluster
       
 end
 
%  for ccol=1:r3
%        value2=round(cluster3(ccol,1),4);
%        df3(round(fuse,4)==value2)=1;        %min cluster
%       
%         % process 2 %
%        
%  end
 
 for i=1:rf
     for j=1:cf
       if (df1(i,j) ~=1)
           df1(i,j)=0;
       end
       if (df2(i,j) ~=1)
           df2(i,j)=0;
       end
%         if(df2(i,j))
%             df4(i,j)=D(i,j);
%         end
%        if (df3(i,j) ~=1)
%            df3(i,j)=0;
%        end
           
     end
 end
 
 
  if(centroid(1,1)>centroid(2,1))
             select=df1;
         else
             select=df2;
 end
 
%  for ms=1:size(df1,1)
%      for ns=1:size(df1,2)
%          if(select(ms,ns))
%              select(ms,ns)=sparse(ms,ns);
%          end
%      end
%  end
 
end