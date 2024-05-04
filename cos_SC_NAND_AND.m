
clear
%for jj = 1:100
format long
N = 256;

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

%LF3(LF3 == -1) = 0;
%lfval = lfval/N;
%lfval2 = lfval2/N;
else
    lfval = rand(1,N);
    lfval2 = rand(1,N);
end

cos_sobol = zeros(1,N);
abs_cos_sob = zeros(1,N);
cos_vdc = zeros(1,N);
abs_cos_vd = zeros(1,N);
cos_lfsr = zeros(1,N);
abs_cos_lf = zeros(1,N);

EE = zeros(1,N);

X2_stream_sobol = zeros(N, N);
X3_stream_sobol = zeros(N, N);
X2_stream_vdc = zeros(N, N);
X3_stream_vdc = zeros(N, N);
X4_stream_vdc = zeros(N, N);
X5_stream_vdc = zeros(N, N);
X2_stream_lfsr_input = zeros(N, N);
X2_stream_lfsr_coeff = zeros(N, N);

for i = 1:N
    for k = 1:N
        %if i/N > sobol(k,12)%12
        if i/N > vd(k,2) %2
            X2_stream_sobol(i,k) = 1;
        end
        if i/N > sobol(k,9)
        %if i/N > vd(k,10)
            X3_stream_sobol(i,k) = 1;
        end
        if i/N > vd(k,4) %3
        %if i/N > z1(k)
            X2_stream_vdc(i,k) = 1;
        end
        if i/N > vd(k,2)%2
        %if i/N > z1(k)
            X3_stream_vdc(i,k) = 1;
        end
        if i/N > vd(k,3)%3
        %if i/N > z1(k)
            X4_stream_vdc(i,k) = 1;
        end
        if i/N > vd(k,5)%5
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

%N=1024, 1/56->18 1/30->34 1/12->85 1/2->512
%N=512, 1/56->9 1/30->17 1/12->43 1/2->256
%N=256, 1/56->4 1/30->8 1/12->21 1/2->128




for i = 1:N
    %EE(i) = cos((i-1)/N);
    EE(i) = 1 - ((i/N)^2)/factorial(2) + ((i/N)^4)/factorial(4) - ((i/N)^6)/factorial(6) + ((i/N)^8)/factorial(8);

    n1_s(i,:) = and(X2_stream_sobol(i,:), circshift(X2_stream_sobol(i,:),1));
    n2_s(i,:) = not(and(n1_s(i,:), X2_stream_sobol(ceil((1/40320)*N),:)));
    n3_s(i,:) = and(n2_s(i,:), X2_stream_sobol(ceil((1/720)*N),:));
    n4_s(i,:) = not(and(n3_s(i,:), circshift(n1_s(i,:),1)));
    n5_s(i,:) = and(n4_s(i,:), X2_stream_sobol(ceil((1/24)*N),:));
    n6_s(i,:) = not(and(n5_s(i,:), circshift(n1_s(i,:),2)));
    n7_s(i,:) = and(n6_s(i,:), X2_stream_sobol(ceil((1/2)*N),:));
    y_s(i,:) = not(and(n7_s(i,:), circshift(n1_s(i,:),3)));
    cos_sobol(i) = sum(y_s(i,:))/N;
    %abs_cos_sob(i) = abs(cos_sobol(i) - EE(i));
    abs_cos_sob(i) = (cos_sobol(i) - EE(i))^2;

    input = X2_stream_sobol(i,:);
    %input = X2_stream_vdc(i,:);
    %input = X2_stream_lfsr_input(i,:);
    n1_v(i,:) = and(input, circshift(input,1));
    n2_v(i,:) = not(and(n1_v(i,:), X2_stream_vdc(ceil((1/40320)*N),:)));
    n3_v(i,:) = and(n2_v(i,:), X3_stream_vdc(ceil((1/720)*N),:));
    n4_v(i,:) = not(and(n3_v(i,:), circshift(n1_v(i,:),0)));
    n5_v(i,:) = and(n4_v(i,:), X4_stream_vdc(ceil((1/24)*N),:));
    n6_v(i,:) = not(and(n5_v(i,:), circshift(n1_v(i,:),0)));
    n7_v(i,:) = and(n6_v(i,:), X5_stream_vdc(ceil((1/2)*N),:));
    y_v(i,:) = not(and(n7_v(i,:), circshift(n1_v(i,:),0)));
    cos_vdc(i) = sum(y_v(i,:))/N;
    %abs_cos_vd(i) = abs(cos_vdc(i) - EE(i));
    abs_cos_vd(i) = (cos_vdc(i) - EE(i))^2;

    n1_lf(i,:) = and(X2_stream_lfsr_input(i,:), circshift(X2_stream_lfsr_input(i,:),1));
    n2_lf(i,:) = not(and(n1_lf(i,:), X2_stream_lfsr_coeff(ceil((1/40320)*N),:)));
    n3_lf(i,:) = and(n2_lf(i,:), X2_stream_lfsr_coeff(ceil((1/720)*N),:));
    n4_lf(i,:) = not(and(n3_lf(i,:), circshift(n1_lf(i,:),1)));
    n5_lf(i,:) = and(n4_lf(i,:), X2_stream_lfsr_coeff(ceil((1/24)*N),:));
    n6_lf(i,:) = not(and(n5_lf(i,:), circshift(n1_lf(i,:),2)));
    n7_lf(i,:) = and(n6_lf(i,:), X2_stream_lfsr_coeff(ceil((1/2)*N),:));
    y_lf(i,:) = not(and(n7_lf(i,:), circshift(n1_lf(i,:),3)));
    cos_lfsr(i) = sum(y_lf(i,:))/N;
    %abs_cos_lf(i) = abs(cos_lfsr(i) - EE(i));
    abs_cos_lf(i) = (cos_lfsr(i) - EE(i))^2;
end
%MSE_sobol = mean(abs_cos_sob)
MSE_vdc = mean(abs_cos_vd)
MSE_lfsr = mean(abs_cos_lf)
%end
%lfsr_vdc = mean(MSE_vdc)
%lfsr_lfsr = mean(MSE_lfsr)