clear
%for jj = 1:100
%clear
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

ln_sob = zeros(1,N);
abs_ln_sob = zeros(1,N);
ln_vdc = zeros(1,N);
abs_ln_vd = zeros(1,N);
ln_lfsr = zeros(1,N);
abs_ln_lf = zeros(1,N);

EE = zeros(1,N);

%X2_stream_unary = zeros(N, N);
X2_stream_sobol = zeros(N, N);
X2_stream_vdc = zeros(N, N);
X3_stream_vdc = zeros(N, N);
X4_stream_vdc = zeros(N, N);
X5_stream_vdc = zeros(N, N);
X2_stream_lfsr_input = zeros(N, N);
X2_stream_lfsr_coeff = zeros(N, N);



for i = 1:N
    for k = 1:N
        %if i/N > sobol(k,4)
        if i/N > vd(k,5)%5
            X2_stream_sobol(i,k) = 1;
        end
        if i/N > vd(k,9)%9
        %if i/N > sobol(k,2)
            X2_stream_vdc(i,k) = 1;
        end
        if i/N > vd(k,8)%6 %5,4---10
        %if i/N > z1(k)
            X3_stream_vdc(i,k) = 1;
        end
        if i/N > vd(k,9)%9 
        %if i/N > z1(k)
            X4_stream_vdc(i,k) = 1;
        end
        if i/N > vd(k,9)%9
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
%     EE(i) = sin((i-1)/N);
%     if i == N
%         EE(i) = sin(i/N);
%     end
    EE(i) = i/N - ((i/N)^2)/2 + ((i/N)^3)/3 - ((i/N)^4)/4 + ((i/N)^5)/5;
    %EE(i) = (i-1)/N - (((i-1)/N)^3)/factorial(3) + (((i-1)/N)^5)/factorial(5) - (((i-1)/N)^7)/factorial(7);
    



    input = X2_stream_sobol(i,:);
    n1_v(i,:) = not(and(input, circshift(X2_stream_vdc(ceil((1/5)*N),:),0)));%best with 15 Delays
    %n1_v(i,:) = and(X4_stream_vdc(i,:), circshift(X4_stream_vdc(i,:),3));
    n2_v(i,:) = and(n1_v(i,:), X2_stream_vdc(ceil((1/4)*N), :));
    n3_v(i,:) = not(and(n2_v(i,:), circshift(input,0)));
    n4_v(i,:) = and(n3_v(i,:), X2_stream_vdc(ceil((1/3)*N), :));
    n5_v(i,:) = not(and(n4_v(i,:), circshift(input,0)));
    n6_v(i,:) = and(n5_v(i,:), X2_stream_vdc(ceil((1/2)*N), :));
    n7_v(i,:) = not(and(n6_v(i,:), circshift(input,0)));
    y_v(i,:) = (and(n7_v(i,:), circshift(input,0)));
    ln_vdc(i) = sum(y_v(i,:))/N;
    %abs_ln_vd(i) = abs(ln_vdc(i) - EE(i)); % Absolute Error
    abs_ln_vd(i) = (ln_vdc(i) - EE(i))^2; % Square Error


    n1_lf(i,:) = not(and(X2_stream_lfsr_input(i,:), circshift(X2_stream_lfsr_coeff(ceil((1/5)*N),:),0)));
    n2_lf(i,:) = and(n1_lf(i,:), circshift(X2_stream_lfsr_coeff(ceil((1/4)*N), :),0));
    n3_lf(i,:) = not(and(n2_lf(i,:), circshift(X2_stream_lfsr_input(i,:),1)));
    n4_lf(i,:) = and(n3_lf(i,:), circshift(X2_stream_lfsr_coeff(ceil((1/3)*N), :),0));
    n5_lf(i,:) = not(and(n4_lf(i,:), circshift(X2_stream_lfsr_input(i,:),2)));
    n6_lf(i,:) = and(n5_lf(i,:), circshift(X2_stream_lfsr_coeff(ceil((1/2)*N), :),0));
    n7_lf(i,:) = not(and(n6_lf(i,:), circshift(X2_stream_lfsr_input(i,:),3)));
    y_lf(i,:) = (and(n7_lf(i,:), circshift(X2_stream_lfsr_input(i,:),4)));
    ln_lfsr(i) = sum(y_lf(i,:))/N;
    %abs_ln_lf(i) = abs(ln_lfsr(i) - EE(i)); % Absolute Error
    abs_ln_lf(i) = (ln_lfsr(i) - EE(i))^2; % Square Error
end
%MSE_sobol = mean(abs_ln_sob)
MSE_vdc = mean(abs_ln_vd)
MSE_lfsr = mean(abs_ln_lf)
%end
%lfsr_vdc = mean(MAE_vdc)
%lfsr_lfsr = mean(MSE_lfsr)