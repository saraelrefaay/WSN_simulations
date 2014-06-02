function [ id ] = get_closest_head( S,head, sink, cluster_head )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
min= ecludian_distance(head.xd, sink.x, head.yd, sink.y);
id = head.id;
for h=1:length(cluster_head)
    next_head = S(cluster_head(h));
    if head.id ~=  next_head.id
        v1 = [head.xd-sink.x, head.yd-sink.y];
        v2 = [head.xd - next_head.xd, head.yd - next_head.yd ];
        if acos(dot(v1,v2)/(norm(v1)*norm(v2)))*180/pi < 90
            dist = ecludian_distance(head.xd,next_head.xd, head.yd, next_head.yd);
            if dist < min
                min = dist;
                id = cluster_head(h);
            end
        end
    end
end