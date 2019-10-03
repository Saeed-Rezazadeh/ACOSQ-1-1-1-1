function [T , codebook , SDR , Distortion] = COSQ_1(Pr , f , T , codebook , delta)
FileID = fopen ('Results.txt' , 'a') ;
D = [1 2] ;

while  abs ((D(2) - D(1)) / D(2)) >= (0.001 /4)
    D(1) = D(2) ;
    %% Optimal Partitions
    parfor u_index = 1 : length(T)
        summation = 0 ;
        temp = zeros (2 , 1);
        u = T(u_index , 1) ;
        for i = 1 : 2
            for j = 1 : 2
                summation = summation + Pr(j , i) * (u - codebook(j)) ^ 2 ;
            end
            temp (i) = summation ;
            summation = 0 ;
        end
        [~  , partition_index ] = min(temp) ;
        T_u(u_index , 1) = partition_index ;
    end
    T(: , 2) = T_u ;
    
    %% Optimal Codebook
    parfor j = 1 : 2
        numerator = 0 ;
        denominator = 0 ;
        for i = 1 : 2
            u_index = find (T(: , 2) == i) ;
            if (isempty(u_index) ~= 1)
                u = T(u_index , 1) ;
                numerator = numerator + Pr(j , i) * sum (u .* f(u_index)) ;
                denominator = denominator + Pr(j , i) * sum (f(u_index)) ;
            end
        end
        codebook(j) = numerator / denominator ;
    end
    
    %% Distortion
    D(2) = distortion_1 (Pr , T , codebook , f , delta) ;
    fprintf (FileID , 'Overall D_1 = %f\n' ,D(2)) ;
end
SDR = 10 * log10(1 / D (2)) ;
Distortion = D(2);
fprintf (FileID , 'SDR_1 = %4.2f\n' , SDR) ;
fclose (FileID) ;
end