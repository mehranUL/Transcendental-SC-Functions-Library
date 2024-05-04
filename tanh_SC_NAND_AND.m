
clear


format long
N = 1024;
%N = 512;
%N = 256;

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

seed_mat = generateBinaryValues(N);
for jj = 1:1000
%z1 = sqrt(-2.*(log(vd(1:end,1)))).*sin(2*pi.*vd(1:end,8));
% z1 = sqrt(-2.*(log(vd(1:end,1)))).*tanh(vd(1:end,1));
% z1 = z1 - floor(z1);

if N == 1024
    %[~,lfval] = LFSR3([false true true true true false false false true true],N/2,N); %N=1024
    %[~,lfval] = LFSR3([0 1 1 1 1 0 0 0 1 1],N/2,N); %N=1024
    [~,lfval] = LFSR3(seed_mat(randi(N),:),N/2,N); %N=1024
    %[~,lfval2] = LFSR3_2([false true false true false false false true false true],N/2,N); %N=1024
    %[~,lfval2] = LFSR3_2([0 1 0 1 0 0 0 1 0 1],N/2,N); %N=1024
    [~,lfval2] = LFSR3_2(seed_mat(randi(N),:),N/2,N);
elseif N == 512
    [~,lfval] = LFSR3([true true true false false false false false false],N/2,N); %N=512
    [~,lfval2] = LFSR3_2([false false true true true true false false true],N/2,N); %N=512
elseif N == 256
    [~,lfval] = LFSR3([false false false true true false false false],N/2,N); %N=256
    [~,lfval2] = LFSR3_2([false false false false true true false false],N/2,N); %N=256
end
lfval = lfval/N;
lfval2 = lfval2/N;
%lfval = rand(1,N);

sin_sob = zeros(1,N);
abs_sin_sob = zeros(1,N);
sin_vdc = zeros(1,N);
abs_sin_vd = zeros(1,N);
sin_lfsr = zeros(1,N);
abs_sin_lf = zeros(1,N);

EE = zeros(1,N);

X2_stream_sobol = zeros(N, N);
X2_stream_vdc = zeros(N, N);
X3_stream_vdc = zeros(N, N);
X4_stream_vdc = zeros(N, N);
X5_stream_vdc = zeros(N, N);
X2_stream_lfsr_input = zeros(N, N);
X2_stream_lfsr_coeff = zeros(N, N);

for i = 1:N
    for k = 1:N
        %if i/N > sobol(k,6)
        if i/N > vd(k,1)%1
            X2_stream_sobol(i,k) = 1;
        end
        if i/N > vd(k,2)%5 for sobol, 2 for vdc, 1 for LFSR----2
        %if i/N > z1(k)
            X2_stream_vdc(i,k) = 1;
        end
        if i/N > vd(k,6)%6 %5,4
        %if i/N > z1(k)
            X3_stream_vdc(i,k) = 1;
        end
        if i/N > vd(k,6)%9 %6
        %if i/N > z1(k)
            X4_stream_vdc(i,k) = 1;
        end
        if i/N > vd(k,9)
        %if i/N > z1(k)
            X5_stream_vdc(i,k) = 1;
        end
        if i/N > lfval(k)
            X2_stream_lfsr_input(i,k) = 1;
        end
        if i/N > lfval2(k)
            X2_stream_lfsr_coeff(i,k) = 1;
        end
    end
end

%N=1024, 1/42->24 1/20->51 1/6->170
%N=512, 1/42->12 1/20->26 1/6->85
%N=256, 1/42->6 1/20->13 1/6->42


for i = 1:N
    %EE(i) = sin((i-1)/N);
    EE(i) = (i/N) - ((i/N)^3)/3 + 2*((i/N)^5)/15 - 17*((i/N)^7)/315;

    n1_s(i,:) = and(X2_stream_sobol(i,:), circshift(X2_stream_sobol(i,:),1));
    n2_s(i,:) = not(and(n1_s(i,:), X2_stream_sobol(ceil((17/315)*N),:)));
    n3_s(i,:) = and(n2_s(i,:), X2_stream_sobol(ceil((2/15)*N),:));
    n4_s(i,:) = not(and(n3_s(i,:), circshift(n1_s(i,:), 1)));%1
    n5_s(i,:) = and(n4_s(i,:), X2_stream_sobol(ceil((1/3)*N),:));
    n6_s(i,:) = not(and(n5_s(i,:), circshift(n1_s(i,:), 2)));%2
    y_s(i,:) = and(n6_s(i,:), circshift(X2_stream_sobol(i,:),2));%2
    tanh_sobol(i) = sum(y_s(i,:))/N;
    %abs_sin_sob(i) = abs(sin_sobol(i) - EE(i));
    abs_tanh_sob(i) = (tanh_sobol(i) - EE(i))^2;

    input = X2_stream_sobol(i,:);
    %input = X2_stream_lfsr_input(i,:);
    n1_v(i,:) = and(input, circshift(input,1));
    n2_v(i,:) = not(and(n1_v(i,:), X2_stream_vdc(ceil((17/315)*N),:)));
    n3_v(i,:) = and(n2_v(i,:), X3_stream_vdc(ceil((2/15)*N),:));
    n4_v(i,:) = not(and(n3_v(i,:), circshift(n1_v(i,:), 0)));%1
    n5_v(i,:) = and(n4_v(i,:), X4_stream_vdc(ceil((1/3)*N),:));
    n6_v(i,:) = not(and(n5_v(i,:), circshift(n1_v(i,:), 0)));%2
    y_v(i,:) = and(n6_v(i,:), circshift(input,0));%2
    tanh_vdc(i) = sum(y_v(i,:))/N;
    %abs_sin_vd(i) = abs(sin_vdc(i) - EE(i));
    abs_tanh_vd(i) = (tanh_vdc(i) - EE(i))^2;

    n1_lf(i,:) = and(X2_stream_lfsr_input(i,:), circshift(X2_stream_lfsr_input(i,:),1));
    n2_lf(i,:) = not(and(n1_lf(i,:), circshift(X2_stream_lfsr_coeff(ceil((17/315)*N),:),0)));
    n3_lf(i,:) = and(n2_lf(i,:), circshift(X2_stream_lfsr_coeff(ceil((2/15)*N),:),0));
    n4_lf(i,:) = not(and(n3_lf(i,:), circshift(n1_lf(i,:), 1)));%1
    n5_lf(i,:) = and(n4_lf(i,:), circshift(X2_stream_lfsr_coeff(ceil((1/3)*N),:),0));
    n6_lf(i,:) = not(and(n5_lf(i,:), circshift(n1_lf(i,:), 2)));%2
    y_lf(i,:) = and(n6_lf(i,:), circshift(X2_stream_lfsr_input(i,:),2));%2
    tanh_lfsr(i) = sum(y_lf(i,:))/N;
    %abs_sin_lf(i) = abs(sin_lfsr(i) - EE(i));
    abs_tanh_lf(i) = (tanh_lfsr(i) - EE(i))^2;
end
%MSE_sobol = mean(abs_tanh_sob)
%MSE_vdc = mean(abs_tanh_vd)
MSE_lfsr(jj) = mean(abs_tanh_lf);
end
%lfsr_vdc = mean(MSE_vdc)
lfsr_lfsr = mean(MSE_lfsr)