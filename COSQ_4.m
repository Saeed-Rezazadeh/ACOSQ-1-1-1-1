function [T , codebook , SDR ,Distortion] = COSQ_4(Pr_z , f , T , y_1 , y_2 , y_3 , codebook , delta)
FileID = fopen ('Results.txt' , 'a') ;
D = [1 2] ;
y_1_2 = (y_1 - 1) * 2 + y_2 ;
while  abs ((D(2) - D(1)) / D(2)) >= (0.001 /4)
    D(1) = D(2) ;
    %% Optimal Partitions
    parfor u_index = 1 : length(T)
        summation = 0 ;
        temp = zeros (2 , 1);
        u = T(u_index , 1) ;
        x_3 = T(u_index , 4 + y_1_2) ;
        for x_4 = 1 : 2
            for y_4 = 1 : 2
                summation = summation + Pr_z(xor(x_3 - 1 , y_3 - 1) + 1 , xor(x_4 - 1 , y_4 - 1) + 1 ) * (u  - codebook(y_4)) ^ 2 ;
            end
            temp (x_4) = summation ;
            summation = 0 ;
            
        end
        [~  , partition_index ] = min(temp) ;
        T_u(u_index , 1) = partition_index ;
    end
    T(: , 9) = T_u ;
    
    %% Optimal Codebook
    parfor y_4 = 1 : 2
        numerator = 0 ;
        denominator = 0 ;
        for x_4 = 1 : 2
            u_index = find (T(: ,  9 ) == x_4) ;
            if (isempty(u_index)~=1)
                for u_i = 1 : length(u_index)
                    x_3 = T(u_index(u_i) , 4 + y_1_2) ;
                    
                    numerator = numerator + Pr_z(xor(x_3 - 1 , y_3 - 1) + 1 , xor(x_4 - 1 , y_4 - 1) + 1 ) *  T(u_index(u_i) , 1) * f(u_index(u_i)) ;
                    denominator = denominator + Pr_z(xor(x_3 - 1 , y_3 - 1) + 1 , xor(x_4 - 1 , y_4 - 1) + 1 ) * f(u_index(u_i)) ;
                end
            end
        end
        codebook(y_4) = numerator / denominator ;
    end
    
    %% Distortion
    D(2) = distortion_4 (Pr_z , T ,  y_1 , y_2 , y_3 , codebook , f , delta) ;
    fprintf (FileID , 'Overall D_4 = %f\n' ,D(2)) ;
end

SDR = 10 * log10(1 / D (2)) ;
Distortion = D(2);
fprintf (FileID , 'SDR_4 = %4.2f\n' , SDR) ;
fclose (FileID) ;
end