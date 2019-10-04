function [SNR , T , codebook]  = COVQ (Pr , T , initial_codebook , k , numLevel )

codebook = initial_codebook ;
FileID = fopen ('Results.txt' , 'a') ;

D = [1 2] ;

while  abs ((D(2) - D(1)) / D(2)) >= (0.001 /4)
    D(1) = D(2) ;
    %% Finding the optimal partitions
    
    parfor u = 1 : length(T)
        summation = 0 ;
        temp = [] ;
        for i = 1 : numLevel
            for j = 1 : numLevel
                summation = summation + Pr (j , i) * sum((T(u , 1:k) - codebook (j , : )) .^2) ;
            end
            temp(i) = summation ;
            summation = 0 ;
        end
        [~ , partition_index ] = min (temp) ;
        T_u (u , 3) = partition_index ;
    end
    T(: , k + 1) = T_u(: , 3) ;
    %% Finding the optimal centroids
    parfor j = 1 : numLevel
        numerator = zeros(1 , k) ;
        denominator = 0 ;
        for i = 1 : numLevel
            x_index = find (T (: , k + 1) == i) ;
            numerator = numerator + Pr (j , i) * sum(T(x_index , 1 : k));
            denominator = denominator + Pr (j , i) * length(x_index) ;
        end
        if (denominator == 0)
            codebook (j , : ) = zeros (1 , k) ;
        else
            codebook ( j , : ) = numerator / denominator ;
        end
    end
    %% Finding D_2
    summation = 0 ;
    parfor i = 1 : length(T)
        for j = 1 : numLevel
            
            summation = summation + ...
                Pr(j , T(i , k + 1 ) ) *...
                sum( (T(i , 1:k ) - codebook(j , : )) .^ 2) ;
            
        end
    end
    D(2) = (1 / (k * length(T))) * summation ;
    fprintf (FileID , 'D = %f \n' , D(2)) ;
end

SNR = 10 * log10 (1 / D(2)) ;

end
% End of the COVQ function