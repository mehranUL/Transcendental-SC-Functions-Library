function [sigmoid_vdc, sigmoid_lfsr] = sigmoid_SC_calculator(X, N)
format long
if X > 1 || X < 0
    fprintf("Error");
else
vd(:,1) = vdcorput(N-1,2);
vd(:,2) = vdcorput(N-1,4);
vd(:,3) = vdcorput(N-1,8);
vd(:,4) = vdcorput(N-1,16);
vd(:,5) = vdcorput(N-1,32);
vd(:,6) = vdcorput(N-1,64);
vd(:,7) = vdcorput(N-1,128);
vd(:,8) = vdcorput(N-1,256);
vd(:,9) = vdcorput(N-1,512);
vd(:,10) = vdcorput(N-1,1024);

seed_mat = generateBinaryValues(N);


if N == 1024
    [~,lfval] = LFSR3(seed_mat(randi(N),:),N/2,N);
    [~,lfval2] = LFSR3_2(seed_mat(randi(N),:),N/2,N);
    lfval = lfval/N;
    lfval2 = lfval2/N;
elseif N == 512
    [~,lfval] = LFSR3(seed_mat(randi(N),:),N/2,N);
    [~,lfval2] = LFSR3_2(seed_mat(randi(N),:),N/2,N);
    lfval = lfval/N;
    lfval2 = lfval2/N;
elseif N == 256
    [~,lfval] = LFSR3(seed_mat(randi(N),:),N/2,N);
    [~,lfval2] = LFSR3_2(seed_mat(randi(N),:),N/2,N);
    lfval = lfval/N;
    lfval2 = lfval2/N;
else
    lfval = rand(1,N);
    lfval2 = rand(1,N);
end
    
X2_stream_input = zeros(1, N);

X2_stream_vdc10 = zeros(1, N);
X3_stream_vdc12 = zeros(1, N);
X4_stream_vdc2 = zeros(1, N);
X5_stream_vdc2 = zeros(1, N);
X2_stream_lfsr = zeros(1, N);
X2_stream_lfsr_coeff10 = zeros(1, N);
X2_stream_lfsr_coeff12 = zeros(1, N);
X2_stream_lfsr_coeff2 = zeros(1, N);
X3_stream_lfsr_coeff2 = zeros(1, N);


    for k = 1:N
        
        if X > vd(k,log2(N))
            X2_stream_input(k) = 1;
        end
        if 1/10 > vd(k,1)
        
            X2_stream_vdc10(k) = 1;
        end
        if 1/12 > vd(k,2)
        
            X3_stream_vdc12(k) = 1;
        end
        if 1/2 > vd(k,2)
        
            X4_stream_vdc2(k) = 1;
        end
        if 1/2 > vd(k,5)
        
            X5_stream_vdc2(k) = 1;
        end
        if X > lfval(k)
            X2_stream_lfsr(k) = 1;
        end
        if 1/10 > lfval2(k)
            X2_stream_lfsr_coeff10(k) = 1;
        end
        if 1/12 > lfval2(k)
            X2_stream_lfsr_coeff12(k) = 1;
        end
        if 1/2 > lfval2(k)
            X2_stream_lfsr_coeff2(k) = 1;
        end
        if 1/2 > lfval2(k)
            X3_stream_lfsr_coeff2(k) = 1;
        end
    end

    input = X2_stream_input;
    %input = X2_stream_vdc(i,:);
    %input = X2_stream_lfsr(i,:);
    n1_v = and(input, circshift(input,2));
    n2_v = not(and(n1_v, X2_stream_vdc10));
    n3_v = not(and3(X3_stream_vdc12, n2_v, circshift(n1_v,0)));
    n4_v = not(and3(X4_stream_vdc2, n3_v, circshift(input,0)));
    y_v = not(and(n4_v, X5_stream_vdc2));
    sigmoid_vdc = sum(y_v)/N;

    n1_lf = and(X2_stream_lfsr, circshift(X2_stream_lfsr,2));
    n2_lf = not(and(n1_lf, X2_stream_lfsr_coeff10));
    n3_lf = not(and3(circshift(X2_stream_lfsr_coeff12, 1), n2_lf, circshift(n1_lf,1)));
    n4_lf = not(and3(circshift(X2_stream_lfsr_coeff2, 1), n3_lf, circshift(X2_stream_lfsr,4)));
    y_lf = not(and(n4_lf, circshift(X3_stream_lfsr_coeff2,0)));
    sigmoid_lfsr = sum(y_lf)/N;

end
end