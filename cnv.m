%Digital Signal Processing
%Paul Adkisson
%10/20/20
%Purpose:
%   Performs convolution between two real sequences h and x using the DFT
%   method.

function y = cnv(h, x)
    %Zero Pad h and x
    N = size(x, 1) + length(h) - 1;
    x_bar = [x; zeros([N - size(x,  1), size(x, 2)])];
    h_bar = [h; zeros([N - length(h), 1])];
    
    %Take DFT
    X_bar = fft(x_bar);
    H_bar = fft(h_bar);
    
    %Extract non-redundant parts
    X_bar_nr = X_bar(1:floor(N/2+1), :); %+1 due to Matlab indexing
    H_bar_nr = H_bar(1:floor(N/2+1)); %+1 due to Matlab indexing
    
    %Solve for Y_nr
    Y_nr = X_bar_nr.*H_bar_nr;
    
    y = ifft(Y_nr, N, 'symmetric');
end