function NNCO(n, m, crifile)
% n: number of microclusters; m: number of dimensions; crifile: output_expand.txt

% Max Version
% NNCO(150,2,'quadratic_clu.ascii');
% NNCO(150,2,'cubic_clu.ascii');
% NNCO(450,9,'ovalhelixes_clu.ascii');
% NNCO(150,4,'iris_output_expand.txt');
% NNCO(500,16,'image_output_expand.txt');

% height of orientation plot
kk=1.5;   
        
% linewidth
linewidth=600/n; 

CRI=load(crifile);
Selected=m;
for i=1:n
   if (CRI(i)>=0)
       R_dist(i)=CRI(i);
   else
       R_dist(i)=0;
   end
   
   %Index(i)=CRI(i+n*1);
   
   for j=1:m
      Color_val(i+(j-1)*n)=(CRI(i+(2+(Selected-1)*m+j-1)*n)+255)/2;
   end
end

clear CRI;

maxY=max(R_dist); 
Y=linspace(0,-maxY/kk,m+1); % Y axis for the microclusters

%hold on

for i=1:n   
   %line([Index(i)+1,Index(i)+1],[0,R_dist(i)],'Color',[0,0,0], 'LineWidth',linewidth);
   line([i,i],[0,R_dist(i)],'Color',[0,0,0], 'LineWidth',linewidth);
   
   for j=1:m
      %line([Index(i)+1,Index(i)+1],[Y(j),Y(j+1)],'Color',[Color_val(i+(j-1)*n)/255, Color_val(i+(j-1)*n)/255, Color_val(i+(j-1)*n)/255], 'LineWidth',linewidth);
      line([i,i],[Y(j),Y(j+1)],'Color',[Color_val(i+(j-1)*n)/255, Color_val(i+(j-1)*n)/255, Color_val(i+(j-1)*n)/255], 'LineWidth',linewidth);
   end
end    

% linewidth for dotted lines
linewidth=150/n;  

% Parallel lines
for i=1:m
   line([0,n],[Y(i),Y(i)],'Color',[0,0,0], 'LineWidth',linewidth, 'LineStyle',':');    
end
line([0,n],[-maxY/kk,-maxY/kk],'Color',[0,0,0], 'LineWidth',linewidth, 'LineStyle',':');

% Vertical Lines
% line([n,n],[-maxY/kk,0],'Color',[0,0,0], 'LineWidth',linewidth, 'LineStyle',':');
% line([0,0],[-maxY/kk,0],'Color',[0,0,0], 'LineWidth',linewidth, 'LineStyle',':');

xlabel('Microcluster Order');
ylabel('Orientation && NNC');
set(gca,'YTick',0:max(R_dist)/5:max(R_dist))
title('NNCO Plot');

