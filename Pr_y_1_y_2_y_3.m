function Probability_y_1_2_3 = Pr_y_1_y_2_y_3( y_1, y_2 , y_3 , Pr_1 , Pr_z , f  , T , delta)
y_1_2 = (y_1 - 1) * 2 + y_2 ;
summation = 0 ;
for x_1 = 1 : 2
    for x_2 = 1 : 2
        for x_3 = 1 : 2
            u_index = find(T(: , 2) == x_1 & T(: , 2 + y_1) == x_2 & T(: , 4 + y_1_2) == x_3) ;
            
            summation = summation + Pr_1(x_1 , y_1) * Pr_z(xor(x_1 - 1 , y_1 - 1) + 1 , xor(x_2 - 1 , y_2 - 1) + 1 ) ...
                * Pr_z(xor(x_2 - 1 , y_2 - 1) + 1 , xor(x_3 - 1 , y_3 - 1) + 1 ) ...
                * delta * sum(f(u_index)) ;
        end
    end
end
Probability_y_1_2_3 = summation ;
end