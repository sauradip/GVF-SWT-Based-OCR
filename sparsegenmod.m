function [tot] = sparsegenmod( mag , dc, val )

% mag=double(mag)/double(max(max(mag)));
% dc=double(dc)/double(max(max(dc)));
sparse6=mag;
sparse_dc_1=dc;
% x=size(mag,1);
% y=size(mag,2);

x=size(mag,1);
y=size(mag,2);
F=zeros(x+5,y+5);
%DC=zeros(x+5,y+5);
sparse=zeros(x,y);
sparse_dc=zeros(x,y);
tot=zeros(x,y);
weight=zeros(5,5);
dim=(x-4)*(y-4);
mean_weight=zeros(dim,1);
median_weight=zeros(dim,1);
nonzero_count=zeros(dim,1);
zero_count=zeros(dim,1);
sd_wrt_mean=zeros(dim,1);
sd_wrt_median=zeros(dim,1);

mean_weight_dc=zeros(dim,1);
median_weight_dc=zeros(dim,1);
nonzero_count_dc=zeros(dim,1);
zero_count_dc=zeros(dim,1);
sd_wrt_mean_dc=zeros(dim,1);
sd_wrt_median_dc=zeros(dim,1);

m= 1;
n=1;
sd_mean=zeros(5,5);
sd_median=zeros(5,5);

sd_mean_dc=zeros(5,5);
sd_median_dc=zeros(5,5);
%sp=sparse(D);

% for text (152,177) to (163,191)

% for i=143:159
%     for j =147:157
for i=1:x-4
    for j =1:y-4
       
       %temp=round(mag(i:i+4,j:j+4),1);
       temp=round(mag(i:i+4,j:j+4),2);
       temp1=round(dc(i:i+4,j:j+4),1);
        %m=i+4;
       %n=j+4;
%        freq = unique((temp));
%        freq_table = [freq,histc((temp(:)),freq)];
%        for p =1:5
%            for q=1:5
%                index=find(freq_table(:,1)==temp(p,q));
%                weight(p,q)=freq_table(index,1).*freq_table(index,2);
%                
%            end
%        end

% for mg 
       mean_weight(m,1)=mean(mean(((temp))));
       median_weight(m,1)=median(median(((temp))));
       zero_count(m,1)=(25-numel(nonzeros(round(temp,1))))/25;
       nonzero_count(m,1)=(numel(nonzeros(round(temp,1))))/25;
       
       for xx=1:5
           for yy =1:5
              sd_mean(xx,yy)= power(temp(xx,yy) -  mean(mean(((temp)))),2);
              sd_median(xx,yy)= power(temp(xx,yy) -  median(median(((temp)))),2);
           end
       end
       sd_wrt_mean(m,1)=sqrt(sum(sum(sd_mean))/25);
       sd_wrt_median(m,1)=sqrt(sum(sum(sd_median))/25);
       
       
       
       %for dc 
       
         mean_weight_dc(m,1)=mean(mean(((temp1))));
       median_weight_dc(m,1)=median(median(((temp1))));
       zero_count_dc(m,1)=(25-numel(nonzeros(round(temp1,1))))/25;
       nonzero_count_dc(m,1)=(numel(nonzeros(round(temp1,1))))/25;
       
       for xx=1:5
           for yy =1:5
              sd_mean_dc(xx,yy)= power(temp1(xx,yy) -  mean(mean(((temp1)))),2);
              sd_median_dc(xx,yy)= power(temp1(xx,yy) -  median(median(((temp1)))),2);
           end
       end
       sd_wrt_mean_dc(m,1)=sqrt(sum(sum(sd_mean_dc))/25);
       sd_wrt_median_dc(m,1)=sqrt(sum(sum(sd_median_dc))/25);
       
       
      % for dc comparing sd(median)>const and median>const 
  % if( round(sqrt(sum(sum(sd_median_dc))/25)*((25-numel(nonzeros(temp1)))/25),3) >= 0.203 && round(sqrt(sum(sum(sd_median_dc))/25)*((25-numel(nonzeros(temp1)))/25),3) <= 0.265   ) %<.105 >0.082
  %   if( round(median(median(((temp1))))*((25-numel(nonzeros(temp1)))/25),3) >= 0.210 && round(median(median(((temp1))))*((25-numel(nonzeros(temp1)))/25),3) <= 0.380  )   %<.158 >0.095
      %for mg comparing sd(median) > const and median>const
   if( round(sqrt(sum(sum(sd_median))/25)*((25-numel(nonzeros(temp)))/25),4) >=   round(sqrt(sum(sum(sd_mean))/25)*((25-numel(nonzeros(temp)))/25),4) + 0.27 )
      if(round(median(median(((temp))))*((25-numel(nonzeros(temp)))/25),3) >=  round(mean(mean(((temp))))*((25-numel(nonzeros(temp)))/25),3)+ 0.65 ) 

%    if( round(sqrt(sum(sum(sd_mean))/25)*((25-numel(nonzeros(temp)))/25),3) < round(sqrt(sum(sum(sd_median))/25)*((25-numel(nonzeros((temp))))/25),3) )
 %        if(round(mean(mean(((temp))))*((25-numel(nonzeros(temp)))/25),3) < round(median(median(((temp))))*((25-numel(nonzeros((temp))))/25),3)) 
               sparse(i,j)=1;
               sparse_dc(i,j)=1;
               tot(i,j)=1;
               
           else
               sparse(i,j)=0;
               sparse_dc(i,j)=0;
               tot(i,j)=0;
%            end
%         end
     end
    end
       
%        if( round(mean(mean(((temp))))*((25-numel(nonzeros(temp)))/25),2) < round(median(median(((temp))))*((25-numel(nonzeros((temp))))/25),2) )
%            sparse(i,j)=1;
%        else
%            sparse(i,j)=0;
%        end
      
       
%        if( ( sqrt(sum(sum(sd_mean))/25)*((25-numel(nonzeros(round(temp,1))))/25) < sqrt(sum(sum(sd_median))/25)*((25-numel(nonzeros(round(temp,1))))/25) ) &&  (mean(mean(((temp))))*((25-numel(nonzeros(round(temp,1))))/25) < median(median(((temp))))*((25-numel(nonzeros(round(temp,1))))/25)))
%            sparse(i,j)=1;
%        end
        m=m+1;
       %mul=freq_table(:,1).*freq_table(:,2);
      
       %F(i+2,j+2)=high;
%        if ( val==1)
%        high=max(max(round(weight,1)));
%        sparse(i+2,j+2)=high;
%        end
       
    end
end

% mg final matrix with text normalised
sparse4=zeros(size(sparse,1),size(sparse,2));
for ms=1:size(sparse,1)
     for ns=1:size(sparse,2)
         if(sparse(ms,ns)==1)
             sparse4(ms,ns)=mag(ms,ns);
         else
              sparse4(ms,ns)=0;
         end
     end
end
% dc with text normalised
sparse7=zeros(size(sparse_dc,1),size(sparse_dc,2));
for mss=1:size(sparse_dc,1)
     for nss=1:size(sparse_dc,2)
         if(sparse_dc(mss,nss)==1)
             sparse7(mss,nss)=dc(mss,nss);
         else
              sparse7(mss,nss)=0;
         end
     end
end
end