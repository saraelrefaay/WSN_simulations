clc;
clear;
%Field Dimensions - x and y maximum (in meters)
xm = 100;
ym = 100;
%x and y Coordinates of the Sink
%sink.x =0.5 * xm;
%sink.y = ym + 50;
sink.x=50;
sink.y=175;
sink.cluster_ids={};
%sink.x=0.5*xm;
%sink.y=0.5*ym;

%Number of Nodes in the field
n = 100;
%Optimal Election Probability of a node to become cluster head
p=0.11;

packetLength =50; 
ctrPacketLength = 20;
%Energy Model (all values in Joules)
%Initial Energy 
Eo = 0.5;
%Eelec=Etx=Erx
ETX=50*0.000000001;
ERX=50*0.000000001;
%Transmit Amplifier types
Efs=10*0.000000000001;
Emp=0.0013*0.000000000001;
%Data Aggregation Energy
EDA=5*0.000000001;

INFINITY = 999999999999999;
%maximum number of rounds
rmax=9999;

max_cluster_nodes = n*p;
num_of_nodes = n;

do=sqrt(Efs/Emp);

max_rounds = 100;
%Creation of the random Sensor Network
%figure(1);
[xc,yc] = meshgrid(1:33:100, 1:33:100);
blocks =[];
for x=1:9
    blocks(x).nodes = {} ;
end

rng(5);
random_nodes = round(rand(1,10)*100)

% for m=1:length(random_nodes)
%     S(random_nodes(m)).integrity = S(random_nodes(m)).integrity + 10;
% end

for i=1:1:n
    S(i).xd=rand(1,1)*xm;
    XR(i)=S(i).xd;
    S(i).yd=rand(1,1)*ym;
    YR(i)=S(i).yd;
    S(i).g=0;
    %initially there are no cluster heads only nodes
    S(i).type='N';
    S(i).E=Eo;
%     S(i).ENERGY=0;
    S(i).id = i;
    S(i).nearest=0;
    S(i).integrity=0;
    
    if(S(i).xd >0 && S(i).xd <xc(1,2) && S(i).yd >= 0 && S(i).yd < yc(2))
        blocks(1).nodes{end+1} = S(i);
    elseif(S(i).xd >=0 && S(i).xd <xc(1,2) && S(i).yd >= yc(2) && S(i).yd < yc(3))
        blocks(2).nodes{end+1} = S(i);
    elseif(S(i).xd >=0 && S(i).xd <xc(1,2) && S(i).yd >= yc(3) && S(i).yd <= yc(4))
        blocks(3).nodes{end+1} = S(i);
    elseif(S(i).xd >=xc(1,2) && S(i).xd <xc(1,3) && S(i).yd >= 0 && S(i).yd < yc(2))
        blocks(4).nodes{end+1} = S(i);
    elseif(S(i).xd >=xc(1,2) && S(i).xd <xc(1,3) && S(i).yd >= yc(2) && S(i).yd < yc(3))
        blocks(5).nodes{end+1} = S(i);
    elseif(S(i).xd >=xc(1,2) && S(i).xd <xc(1,3) && S(i).yd >= yc(3) && S(i).yd <= yc(4))
        blocks(6).nodes{end+1} = S(i);
    elseif(S(i).xd >=xc(1,3) && S(i).xd <=xc(1,4) && S(i).yd >= 0 && S(i).yd < yc(2))
        blocks(7).nodes{end+1}= S(i);
    elseif(S(i).xd >=xc(1,3) && S(i).xd <=xc(1,4) && S(i).yd >= yc(2) && S(i).yd < yc(3))
        blocks(8).nodes{end+1} = S(i);
    elseif(S(i).xd >=xc(1,3) && S(i).xd <=xc(1,4) && S(i).yd >= yc(3) && S(i).yd <= yc(4))
        blocks(9).nodes{end+1} = S(i);
%     else
%         S(i)
    end
        
    % hold on;
end



results=[];

for round=1:max_rounds
     fprintf('round %d\n',round);
    failed_heads =0;
    blocks_count = length(blocks);
    m=1;
    while m<=blocks_count
        fprintf('blocks %d,%d\n',blocks_count, m);
         nodes_count = length(blocks(m).nodes);
         i=1;
         while i <=nodes_count
%          for i=1:nodes_count
            fprintf('working with %d,%d\n',nodes_count, i);
            for j=1:nodes_count
                if(i~=j)
%                     fprintf('working with %d,%d\n',nodes_count, i);
                    distance = ecludian_distance(blocks(m).nodes{i}.xd, blocks(m).nodes{j}.xd, blocks(m).nodes{i}.yd, blocks(m).nodes{j}.yd);
                    if(distance <= 10)
                        blocks(m).nodes{i}.nearest= blocks(m).nodes{i}.nearest+ 1;
                    end
                end
            end
            if(~isempty(find(random_nodes == blocks(m).nodes{i}.id)))
                blocks(m).nodes{i}.integrity = blocks(m).nodes{i}.integrity + 10;
            end
            blocks(m).nodes{i}.E = blocks(m).nodes{i}.E - 0.01;
            if(blocks(m).nodes{i}.E <= 0 )
                fprintf('low energy deleting %d,%d\n',nodes_count, i);
                blocks(m).nodes{i} = [];
                blocks(m).nodes= blocks(m).nodes(~cellfun('isempty',blocks(m).nodes));
                nodes_count= nodes_count -1;
                i = i-1;
%                 if(i <= 0)
%                     i=1;
%                 end
            else
                blocks(m).nodes{i}.g = blocks(m).nodes{i}.E + blocks(m).nodes{i}.nearest - blocks(m).nodes{i}.integrity;
%                 grade(i) = blocks(m).nodes{i}.g;
            end
            i=i+1;
         end
         grade = [];
         if length(blocks(m).nodes) > 0
            for i=1:length(blocks(m).nodes)
                grade(i) = blocks(m).nodes{i}.g;
            end
            grade
            [grade_max, grade_index] = max(grade);
            blocks(m).nodes{grade_index}.type = 'C';
            blocks(m).head = blocks(m).nodes{grade_index};
         else
             fprintf('empty block deleting %d, %d\n',blocks_count, m);
             blocks(m) = [];
%              blocks= blocks(~cellfun('isempty',blocks));
             blocks_count = blocks_count - 1;
             m=m-1;
%              if(m<= 0)
%                  m=1;
%              end             
         end
         m = m+1;
        %     [row,col] = find(dist == min(dist));
        %     blocks(m).head= blocks(m).nodes{col};
        %     blocks(m).head;
        %     blocks(m).head.type = 'C';
        %     blocks(m).nodes{col}.type = 'C';
    end
    fprintf('length of block %d\n',length(blocks));
    
    for m=1:length(blocks)
        fprintf('head %d\n',blocks(m).head.id);
        if(isempty(find(random_nodes == blocks(m).head.id)) == 0)
            blocks(m).head.id
            failed_heads = failed_heads+1;
        end
    end
    results(end+1,:) = [round, failed_heads];
end
 
 results
% 
% outliers = [];
% for m=1:length(blocks)
%     bl= length(blocks(m).nodes);
%     if(length(blocks(m).nodes) > max_cluster_nodes-1)
%         dist=[];
%        for n=1:length(blocks(m).nodes)
%            if(blocks(m).nodes{n}.type == 'N')
%             dist(end+1) = ecludian_distance(blocks(m).head.xd, blocks(m).nodes{n}.xd, blocks(m).head.yd, blocks(m).nodes{n}.yd);
%            end          
%        end
%        [dist, sorted_index] = sort(dist);
%        sub_index = sorted_index(1:max_cluster_nodes-1);       
%        outliers = [outliers sorted_index(max_cluster_nodes:end)];
%        out = length(outliers);
%        for v=1:(max_cluster_nodes-1)
%           blocks(m).head.ids{end+1} = S(sub_index(v));
%        end
%     else
%         for v=1:length(blocks(m).nodes)
%             if(blocks(m).nodes{v}.type == 'N')
%                 blocks(m).head.ids{end+1} = blocks(m).nodes{v};
%             end
%         end
%     end
%     head_ids= length(blocks(m).head.ids)
%     
% end
% length(outliers)
% 
% % for k=2:length(dist)
% %     if (dist(k)<shortest)
% %         shortest= dist(k);
% %     end
% %    
% % end
% 
% 
% 
% S(n+1).xd=sink.x;
% S(n+1).yd=sink.y;
%  
% 
% %random cluster head%
% 
% %cluster_head= ;
% % collected_head_index = n*p;
% 
% % 
% % node_location_x= n;
% % node_location_y = n;
% % head_pos_x = n*p;
% % head_pos_y = n*p;
%  
% % S(1).node_ids =[];
% % S(1).type = 'C';
% % head_pos_x(1) = S(1).xd;
% % head_pos_y(1) = S(1).yd;
% % node_location_x(1)=S(1).xd;
% % node_location_y(1)=S(1).yd;


