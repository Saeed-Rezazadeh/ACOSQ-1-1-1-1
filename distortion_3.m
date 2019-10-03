function D = distortion_3 (Pr_z , T , y_1 , y_2 , codebook , f , delta)
summation = 0 ;
for y_3 = 1 : 2
    for x_3 = 1 : 2
        u_index = find (T(: , 5) == x_3) ;
        if(isempty(u_index) ~= 1)
            for u_i = 1 : length(u_index)
                x_2 = T(u_index(u_i) , 2 + y_1) ; 
                summation = summation + Pr_z(xor(x_2 - 1 , y_2 - 1) + 1 , xor(x_3 - 1 , y_3 - 1) + 1 ) ... 
                    * delta * f(u_index(u_i)) .* (T(u_index(u_i) , 1) - codebook(y_3)) .^ 2 ;
            end
        end
    end
end
D = summation ;
end