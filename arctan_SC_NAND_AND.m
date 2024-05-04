
clear
%for jj = 1:100

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

%z1 = sqrt(-2.*(log(vd(1:end,1)))).*sin(2*pi.*vd(1:end,8));
% z1 = sqrt(-2.*(log(vd(1:end,1)))).*tanh(vd(1:end,1));
% z1 = z1 - floor(z1);
seed_mat = generateBinaryValues(N);
%for jj = 1:1024

if N == 1024
    %[~,lfval] = LFSR3([false true true true true false false false true true],N/2,N); %N=1024
    %[~,lfval2] = LFSR3_2([false true false true false false false true false true],N/2,N); %N=1024
%     s1 = randi(N);
%     s2 = randi(N);
    [~,lfval] = LFSR3(seed_mat(randi(N),:),N/2,N);
    [~,lfval2] = LFSR3_2(seed_mat(randi(N),:),N/2,N);
elseif N == 512
    %[~,lfval] = LFSR3([true true true false false false false false false],N/2,N); %N=512
    %[~,lfval2] = LFSR3_2([false false true true true true false false true],N/2,N); %N=512
    [~,lfval] = LFSR3(seed_mat(randi(N),:),N/2,N);
    [~,lfval2] = LFSR3_2(seed_mat(randi(N),:),N/2,N);
elseif N == 256
    %[~,lfval] = LFSR3([false false false true true false false false],N/2,N); %N=256
    %[~,lfval2] = LFSR3_2([false false false false true true false false],N/2,N); %N=256
    [~,lfval] = LFSR3(seed_mat(randi(N),:),N/2,N);
    [~,lfval2] = LFSR3_2(seed_mat(randi(N),:),N/2,N);
end
lfval = lfval/N;
lfval2 = lfval2/N;
%lfval = rand(1,N);

arctan_sob = zeros(1,N);
abs_arctan_sob = zeros(1,N);
arctan_vdc = zeros(1,N);
abs_arctan_vd = zeros(1,N);
arctan_lfsr = zeros(1,N);
abs_arctan_lf = zeros(1,N);

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
        %if i/N > sobol(k,9)
        if i/N > vd(k,2)%2
            X2_stream_sobol(i,k) = 1;
        end
        if i/N > vd(k,1)%5 for sobol, 10 for vdc, 1 for LFSR
        %if i/N > z1(k)
            X2_stream_vdc(i,k) = 1;
        end
        if i/N > vd(k,6)%6 %5,4
        %if i/N > z1(k)
            X3_stream_vdc(i,k) = 1;
        end
        if i/N > vd(k,9)%9 %1
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
    %EE(i) = arctan((i-1)/N);
    EE(i) = i/N - ((i/N)^3)/3 + ((i/N)^5)/5 - ((i/N)^7)/7;

    n1_s(i,:) = and(X2_stream_sobol(i,:), circshift(X2_stream_sobol(i,:),1));
    n2_s(i,:) = not(and(n1_s(i,:), X2_stream_sobol(ceil((1/7)*N),:)));
    n3_s(i,:) = and(n2_s(i,:), X2_stream_sobol(ceil((1/5)*N),:));
    n4_s(i,:) = not(and(n3_s(i,:), circshift(n1_s(i,:), 1)));%1
    n5_s(i,:) = and(n4_s(i,:), X2_stream_sobol(ceil((1/3)*N),:));
    n6_s(i,:) = not(and(n5_s(i,:), circshift(n1_s(i,:), 2)));%2
    y_s(i,:) = and(n6_s(i,:), circshift(X2_stream_sobol(i,:),2));%2
    arctan_sobol(i) = sum(y_s(i,:))/N;
    %abs_arctan_sob(i) = abs(arctan_sobol(i) - EE(i));
    abs_arctan_sob(i) = (arctan_sobol(i) - EE(i))^2;

    input = X2_stream_sobol(i,:);
    %input = X2_stream_lfsr_input(i,:);
    n1_v(i,:) = and(input, circshift(input,1));
    n2_v(i,:) = not(and(n1_v(i,:), X2_stream_vdc(ceil((1/7)*N),:)));
    n3_v(i,:) = and(n2_v(i,:), X2_stream_vdc(ceil((1/5)*N),:));
    n4_v(i,:) = not(and(n3_v(i,:), circshift(n1_v(i,:), 0)));%1
    n5_v(i,:) = and(n4_v(i,:), X2_stream_vdc(ceil((1/3)*N),:));
    n6_v(i,:) = not(and(n5_v(i,:), circshift(n1_v(i,:), 0)));%2
    y_v(i,:) = and(n6_v(i,:), circshift(input,0));%2
    arctan_vdc(i) = sum(y_v(i,:))/N;
    %abs_arctan_vd(i) = abs(arctan_vdc(i) - EE(i));
    abs_arctan_vd(i) = (arctan_vdc(i) - EE(i))^2;

    n1_lf(i,:) = and(X2_stream_lfsr_input(i,:), circshift(X2_stream_lfsr_input(i,:),1));
    n2_lf(i,:) = not(and(n1_lf(i,:), circshift(X2_stream_lfsr_coeff(ceil((1/7)*N),:),0)));
    n3_lf(i,:) = and(n2_lf(i,:), circshift(X2_stream_lfsr_coeff(ceil((1/5)*N),:),1));
    n4_lf(i,:) = not(and(n3_lf(i,:), circshift(n1_lf(i,:), 1)));%1
    n5_lf(i,:) = and(n4_lf(i,:), circshift(X2_stream_lfsr_coeff(ceil((1/3)*N),:),1));
    n6_lf(i,:) = not(and(n5_lf(i,:), circshift(n1_lf(i,:), 2)));%2
    y_lf(i,:) = and(n6_lf(i,:), circshift(X2_stream_lfsr_input(i,:),2));%2
    arctan_lfsr(i) = sum(y_lf(i,:))/N;
    %abs_arctan_lf(i) = abs(arctan_lfsr(i) - EE(i));
    abs_arctan_lf(i) = (arctan_lfsr(i) - EE(i))^2;
end
%MSE_sobol = mean(abs_arctan_sob)
MSE_vdc = mean(abs_arctan_vd)
MSE_lfsr = mean(abs_arctan_lf)
%end
%lfsr_vdc = mean(MSE_vdc)
%lfsr_lfsr = mean(MSE_lfsr)