function [ sparse ] = sparsegen( mag , val )
x=size(mag,1);
y=size(mag,2);
F=zeros(x,y);
%DC=zeros(x+5,y+5);
sparse=zeros(x,y);
weight=zeros(5,5);
% mean_weight=zeros(200,1);
% median_weight=zeros(200,1);

m= 1;
n=1;
%sp=sparse(D);
for i=1:x-4
    for j =1:y-4
       
       temp=round(mag(i:i+4,j:j+4),1);
        %m=i+4;
       %n=j+4;
       freq = unique((temp));
       freq_table = [freq,histc((temp(:)),freq)];
       for p =1:5
           for q=1:5
               index=find(freq_table(:,1)==temp(p,q));
               weight(p,q)=freq_table(index,1).*freq_table(index,2);
               
           end
       end
%        mean_weight(m,1)=mean((round(weight,1)));
%        median_weight(m,1)=median((round(weight,1)));
%        zero_count(m,1)=25-numel(nonzeros(round(weight,1)));
%        nonzero_count(m,1)=numel(nonzeros(round(weight,1)));
       
      % m=m+1;
       sparse(i:i+4,j:j+4)=round(weight,1);
       %mul=freq_table(:,1).*freq_table(:,2);
      
       %F(i+2,j+2)=high;
       if ( val==1)
       high=max(max(round(weight,1)));
       sparse(i+2,j+2)=high;
       end
       
    end
end

end