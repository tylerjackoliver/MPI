function [average_io, average_time, average_per, ...
    average_gather] = time_per_iteration(runtimes)

    tot_io = zeros(6,1);
    tot_time = zeros(6,1);
    tot_per = zeros(6,1);
    tot_gather = zeros(6,1);
    
    for i = 1:5
        
        for j = 1:6
            
            tot_io(j) = tot_io(j) + runtimes((i-1)*6 + j, 3);
            tot_time(j) = tot_time(j) + runtimes((i-1)*6 + j, 3);
            tot_per(j) = tot_per(j) + runtimes((i-1)*6 + j, 3);
            tot_gather(j) = tot_gather(j) + runtimes((i-1)*6 + j, 3);
            
        end
        
    end
    
    average_io = tot_io / 5;
    average_time = tot_time / 5;
    average_per = tot_per / 5;
    average_gather = tot_gather / 5;

end