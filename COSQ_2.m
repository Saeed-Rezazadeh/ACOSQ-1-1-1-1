function [T , codebook , SDR , Distortion] = COSQ_2(Pr_z , f , T  , y_1 , codebook , delta)
FileID = fopen ('Results.txt' , 'a') ;
D = [1 2] ;

while  abs ((D(2) - D(1)) / D(2)) >= (0.001 /4)
    D(1) = D(2) ;
    %% Optimal Partitions
    parfor u_index = 1 : length(T)
        summation = 0 ;
        temp = zeros (2 , 1);
        u = T(u_index , 1) ;
        x_1 = T(u_index , 2) ;
        for x_2 = 1 : 2
            for y_2 = 1 : 2
                summation = summation + Pr_z(xor(x_1 - 1 , y_1 - 1) + 1 , xor(x_2 - 1 , y_2 - 1) + 1 ) * (u  - codebook(y_2)) ^ 2 ;
            end
            temp (x_2) = summation ;
            summation = 0 ;
        end
        [~  , partition_index ] = min(temp) ;
        T_u(u_index , 1) = partition_index ;
    end
    T(: , 3) = T_u ;
    
    %% Optimal Codebook
    for y_2 = 1 : 2
        numerator = 0 ;
        denominator = 0 ;
        for x_2 = 1 : 2
            u_index = find (T(: , 3) == x_2) ;
            if(isempty(u_index)~= 1)
                for u_i = 1 : length(u_index)
                    x_1 = T(u_index(u_i) , 2) ; 
                    numerator = numerator + Pr_z(xor(x_1 - 1 , y_1 - 1) + 1 , xor(x_2 - 1 , y_2 - 1) + 1 ) * T(u_index(u_i) , 1).* f(u_index(u_i)) ;
                    denominator = denominator + Pr_z(xor(x_1 - 1 , y_1 - 1) + 1 , xor(x_2 - 1 , y_2 - 1) + 1 ) * f(u_index(u_i)) ;
                end
            end
        end
        codebook(y_2) = numerator / denominator ;
    end
    
    %% Distortion
    D(2) = distortion_2 (Pr_z , T , y_1 , codebook , f , delta) ;
    fprintf (FileID , 'Overall D_2 = %f\n' ,D(2)) ;
end
SDR = 10 * log10(1 / D (2)) ;
Distortion = D(2);
fprintf (FileID , 'SDR_2 = %4.2f\n' , SDR) ;
fclose (FileID) ;
end