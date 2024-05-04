
clear
format long
N = 1024;

sobol = net(sobolset(256), N);
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

%z1 = sqrt(-2.*(log(vd(1:end,1)))).*sin(2*pi.*vd(1:end,8));
% z1 = sqrt(-2.*(log(vd(1:end,1)))).*tanh(vd(1:end,1));
% z1 = z1 - floor(z1);


seed_mat = generateBinaryValues(N);
%for jj = 1:1024

%LFSR3([true false true true true false false false false false],N/2,N);

if N == 1024
    %[~,lfval] = LFSR3([true false true true true false false true false false],N/2,N); %N=1024
    %[~,lfval2] = LFSR3_2([false true false true false false false true false true],N/2,N); %N=1024
    [~,lfval] = LFSR3(seed_mat(randi(N),:),N/2,N);
    [~,lfval2] = LFSR3_2(seed_mat(randi(N),:),N/2,N);
    lfval = lfval/N;
    lfval2 = lfval2/N;
elseif N == 512
    %[~,lfval] = LFSR3([true false true false true false false false false],N/2,N); %N=512
    %[~,lfval2] = LFSR3_2([false false true true true true false false true],N/2,N); %N=512
    [~,lfval] = LFSR3(seed_mat(randi(N),:),N/2,N);
    [~,lfval2] = LFSR3_2(seed_mat(randi(N),:),N/2,N);
    lfval = lfval/N;
    lfval2 = lfval2/N;
elseif N == 256
    %[~,lfval] = LFSR3([true false true false true false false false],N/2,N); %N=256
    %[~,lfval2] = LFSR3_2([false false false false true true false false],N/2,N); %N=256
    [~,lfval] = LFSR3(seed_mat(randi(N),:),N/2,N);
    [~,lfval2] = LFSR3_2(seed_mat(randi(N),:),N/2,N);
    lfval = lfval/N;
    lfval2 = lfval2/N;
else
    lfval = rand(1,N);
    lfval2 = rand(1,N);
end


sigmoid_sobol = zeros(1,N);
abs_sigmoid_sob = zeros(1,N);
sigmoid_vdc = zeros(1,N);
abs_sigmoid_vd = zeros(1,N);
sigmoid_lfsr = zeros(1,N);
abs_sigmoid_lf = zeros(1,N);

EE = zeros(1,N);

X2_stream_sobol = zeros(N, N);
X3_stream_sobol = zeros(N, N);
X2_stream_vdc = zeros(N, N);
X3_stream_vdc = zeros(N, N);
X4_stream_vdc = zeros(N, N);
X5_stream_vdc = zeros(N, N);
X2_stream_lfsr = zeros(N, N);
X2_stream_lfsr_coeff = zeros(N, N);

for i = 1:N
    for k = 1:N
        %if i/N > sobol(k,1)
        if i/N > vd(k,7) %7
            X2_stream_sobol(i,k) = 1;
        end
        if i/N > sobol(k,220)
        %if i/N > vd(k,10)
            X3_stream_sobol(i,k) = 1;
        end
        if i/N > vd(k,2) %2-->1024  1-->256
        %if i/N > z1(k)
            X2_stream_vdc(i,k) = 1;
        end
        if i/N > vd(k,2)
        %if i/N > z1(k)
            X3_stream_vdc(i,k) = 1;
        end
        if i/N > vd(k,5)
        %if i/N > z1(k)
            X4_stream_vdc(i,k) = 1;
        end
        if i/N > vd(k,9)
        %if i/N > z1(k)
            X5_stream_vdc(i,k) = 1;
        end
        if i/N > lfval(k)
            X2_stream_lfsr(i,k) = 1;
        end
        if i/N > lfval2(k)
            X2_stream_lfsr_coeff(i,k) = 1;
        end
    end
end



for i = 1:N
    %EE(i) = sin((i-1)/N);
    EE(i) = 1/2 + (i/N)/4 - ((i/N)^3)/48 + ((i/N)^5)/480;

    n1_s(i,:) = and(X2_stream_sobol(i,:), circshift(X2_stream_sobol(i,:),2));
    n2_s(i,:) = not(and(n1_s(i,:), X2_stream_sobol(ceil((1/10)*N),:)));
    n3_s(i,:) = not(and3(X2_stream_sobol(ceil((1/12)*N),:), n2_s(i,:), circshift(n1_s(i,:),1)));
    n4_s(i,:) = not(and3(X2_stream_sobol(ceil((1/2)*N),:), n3_s(i,:), circshift(X2_stream_sobol(i,:),4)));
    y_s(i,:) = not(and(n4_s(i,:), X3_stream_sobol(ceil((1/2)*N),:)));
    sigmoid_sobol(i) = sum(y_s(i,:))/N;
    %abs_sigmoid_sob(i) = abs(sigmoid_sobol(i) - EE(i));
    abs_sigmoid_sob(i) = (sigmoid_sobol(i) - EE(i))^2;

    input = X2_stream_sobol(i,:);
    %input = X2_stream_vdc(i,:);
    %input = X2_stream_lfsr(i,:);
    n1_v(i,:) = and(input, circshift(input,1));
    n2_v(i,:) = not(and(n1_v(i,:), X2_stream_vdc(ceil((1/480)*N),:)));
    n3_v(i,:) = and(n2_v(i,:), X2_stream_vdc(ceil((1/48)*N),:));
    n4_v(i,:) = not(and(n3_v(i,:), circshift(n1_v(i,:),0)));
    n5_v(i,:) = and(n4_v(i,:), X2_stream_vdc(ceil((1/4)*N), :));
    n6_v(i,:) = not(and(n5_v(i,:), circshift(input,0)));
    y_v(i,:) = not(and(n6_v(i,:), X2_stream_vdc(ceil((1/2)*N),:)));
    sigmoid_vdc(i) = sum(y_v(i,:))/N;
    %abs_sigmoid_vd(i) = abs(sigmoid_vdc(i) - EE(i));
    abs_sigmoid_vd(i) = (sigmoid_vdc(i) - EE(i))^2;

    n1_lf(i,:) = and(X2_stream_lfsr(i,:), circshift(X2_stream_lfsr(i,:),1));
    n2_lf(i,:) = not(and(n1_lf(i,:), X2_stream_lfsr_coeff(ceil((1/480)*N),:)));
    n3_lf(i,:) = and(n2_lf(i,:), X2_stream_lfsr_coeff(ceil((1/48)*N),:));
    n4_lf(i,:) = not(and(n3_lf(i,:), circshift(n1_lf(i,:),1)));
    n5_lf(i,:) = and(n4_lf(i,:), X2_stream_lfsr_coeff(ceil((1/4)*N), :));
    n6_lf(i,:) = not(and(n5_lf(i,:), circshift(X2_stream_lfsr(i,:),2)));
    y_lf(i,:) = not(and(n6_lf(i,:), X2_stream_lfsr_coeff(ceil((1/2)*N),:)));
    sigmoid_lfsr(i) = sum(y_lf(i,:))/N;
    %abs_sigmoid_lf(i) = abs(sigmoid_lfsr(i) - EE(i));
    abs_sigmoid_lf(i) = (sigmoid_lfsr(i) - EE(i))^2;
end
%MSE_sobol = mean(abs_sigmoid_sob)
MSE_vdc = mean(abs_sigmoid_vd)
MSE_lfsr = mean(abs_sigmoid_lf)
%end
%lfsr_lfsr = mean(MSE_lfsr)