function [ sent ] = send_sink(S, sink,head_id, cluster_head)
%Choose shortest path to send data from head to sink
%   Detailed explanation goes here
   % fprint('working with head = %c', head_index)
    %head_id
    head = S(head_id);
    sent= false;
    sink_dist = ecludian_distance(head.xd, sink.x, head.yd, sink.y);
    if(sink_dist < 10 || length(cluster_head) == 1)
        disp('dist sink less than 10')
        sent = true;
        return
    else
         head.next = get_closest_head(S, head, sink,cluster_head);
         if head.next == head.id
            % disp('final node')
             sent = true;
            return
         else
            %fprint('next head = %d', head.next)
             head_left = cluster_head(find(cluster_head~=head_id));
            % head_left
             sent = send_sink(S, sink, head.next, head_left);
            
         end
    end
end

