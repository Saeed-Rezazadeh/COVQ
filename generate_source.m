function T = generate_source (u , sigma , alpha , k , p)

data = zeros (alpha * k , 1) ;
T = zeros (alpha , k) ;
X_n = sigma * randn +  u ;

W_n = sqrt(1 - p ^ 2) * randn (1 , k * alpha) ;
for i = 1 : alpha * k
    X_n  = p * X_n +  W_n(i) ;
    data (i) = X_n ;
end

j = 0 ;
for i = 1 : alpha
    T(i , : ) = data (1 + k * j  : k + k * j)  ;
    j = j + 1 ;
end

end