function [f_u_given_y_1] = generate_pdf_rate_2(Pr , f , T , y_1 , delta)
denominator = 0 ;
f_u_given_y_1 = zeros (length(T) , 1) ;
for u_index = 1 : length(T)
    x_1 = T(u_index , 2) ;
    numerator = Pr(x_1 , y_1) * f(u_index) ;
    for x_1 = 1 : 2
        u_index_x_1 = find (T(: , 2) == x_1) ;
        denominator = denominator + Pr(x_1 , y_1) * delta * sum(f(u_index_x_1)) ;
    end
    f_u_given_y_1(u_index) = numerator / denominator ;
    denominator = 0 ;
end
f_u_given_y_1 = f_u_given_y_1 ./ (sum(f_u_given_y_1) * delta ) ;
end