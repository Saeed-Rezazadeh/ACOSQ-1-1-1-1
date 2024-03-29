function [f_u_given_y_1_y_2_y_3] = generate_pdf_rate_4(Pr_z ,  f_u_given_y_1_y_2 , T , y_1 , y_2 , y_3 , delta)
y_1_2 = (y_1 - 1) * 2 + y_2 ;
%% numerator
numerator = zeros (length(T) , 1) ;
for u_index = 1 : length(T)
    x_3 = T(u_index , 4 + y_1_2) ;
    x_2 = T(u_index , 2 + y_1) ;
    numerator(u_index) = Pr_z(xor(x_2 - 1 , y_2 - 1) + 1 , xor(x_3 - 1 , y_3 - 1) + 1 ) * f_u_given_y_1_y_2(u_index) ;
end

%% denominator
denominator = 0 ;
for x_3 = 1 : 2
    u_index_x_3 = find(T(: , 4 + y_1_2) == x_3) ;
    for u_i = 1 : length(u_index_x_3)
        x_2 = T(u_index_x_3(u_i) , 2 + y_1) ;
        denominator = denominator + Pr_z(xor(x_2 - 1 , y_2 - 1) + 1 , xor(x_3 - 1 , y_3 - 1) + 1 ) * delta * f_u_given_y_1_y_2(u_index_x_3(u_i)) ;
    end
end
f_u_given_y_1_y_2_y_3 = numerator ./ denominator ;
f_u_given_y_1_y_2_y_3 = f_u_given_y_1_y_2_y_3 ./ (sum(f_u_given_y_1_y_2_y_3) * delta ) ;
end