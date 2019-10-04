function [T , indexed_T] = initialization (sigma , u , p , alpha , k , numLevel , epsilon )

T = generate_source(u , sigma , alpha , k , p ) ;

indexed_T = cat (2 , T , zeros(length(T) , 1)) ; 

end
% End of function