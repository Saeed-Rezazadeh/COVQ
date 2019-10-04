% The script corresponds to the Algorithm 2 using traning set.
clc ;
clear all;
close all
%% The Results.txt
% This file contains the ultimate SDR values for different channel parameters delta and epsilon.
% Also, since the ACOSQ is an iterative algorithm, the distirtion
% value for every iteration is provided in this file for given a epsilon and delta.
FileID = fopen ('Results.txt' , 'a') ;

sigma = 1 ; % Determine the source's standard devision. 
u = 0 ; % Determine the source's mean value. 
rho = [0 , 0.5 , 0.9]; % Determine the amount of correlation between successive source samples. 
k = 2 ; % Determine the source's blocklength. 
alpha = 200000 ; % Determine the size of traning set. Ideally choose alpha = 50000 * k 
QuantizationLevels = [2 4 8 16 64 256] ; % Determine the number of Quantization levels. Better be chosen as 2^k to insure one-bit per source sample.  

%% Channel's cross-over probability epsilon
epsilon = unique ([0 10^-5 : 10^-5 : 10^ -4 , 10^-4 : 10^ -4 : 10^ -3 , 10^-3 , 0.005 0.01 0.05 0.1]);

% Since the design converges to a locally optimal solution, to avoid
% bad local optimums, we use a increase-decrease method.
SIZE = length(epsilon) ;
noise_index = [1 : SIZE , SIZE : -1 : 1 , 1 : SIZE , SIZE : -1 : 1] ;

SNR = zeros (length (rho) , length(noise_index) , length(QuantizationLevels)) ;
%% Determine the noise correlation. 
delta = 0 ; 
for rho_index = 1 : length(rho)
    for j = 1 : length(QuantizationLevels)
        numLevel = QuantizationLevels(j) ;
        p = rho (rho_index) ;
        [T , indexed_T] = initialization (sigma , u , p , alpha , k , numLevel) ;
        
        % Choose the initial codebook using the splitting algorithm. 
        [I  , codebook] = kmeans (T , numLevel , 'MaxIter' , 1000 ,  'OnlinePhase','on') ;
        for h = 1 : length(noise_index)
            i = noise_index(h) ;
            Pr = Channel_with_Memory(numLevel , epsilon(i) , delta) ; 
            [SNR(rho_index , h , j) , indexed_T , codebook]  = ...
                COVQ (Pr , indexed_T , codebook , k , numLevel ) ;
            fprintf (FileID , 'SDR = %f \n' , SNR(rho_index , h , j)) ;
            
            Data = ['T\Data_p_' num2str(rho_index) , '_i_' , num2str(h) '_numLevel_' num2str(numLevel)];
            save (Data , 'indexed_T' , 'codebook') ;
        end
        fprintf (FileID , '\nRate = %d \n' , log2(numLevel) / k) ;
    end
end
final_SNR = max(SNR , [] , 2) ;
save ('MyData' , 'final_SNR')
