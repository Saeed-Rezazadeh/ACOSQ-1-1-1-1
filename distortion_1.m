function D = distortion_1 (Pr , T , codebook , f , delta)
summation = 0 ;
parfor j = 1 : 2
    for i = 1 : 2
        u_index = find (T(: , 2) == i) ;
        if (isempty(u_index) ~= 1)
            u = T(u_index , 1) ;
            summation = summation + Pr(j , i) * delta * sum(f(u_index) .* (u - codebook(j)) .^ 2) ;
        end
    end
end
D = summation ;
end