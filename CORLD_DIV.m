%CORDIV Stochastic Divider using CORLD design.
function [Z_sobol, Z_CORDIV_vd] = CORLD_DIV(X,Y,N)
% X is dividend
% Y is divisor
% N is the bit-stream length.


sobol_seq = net(sobolset(100), N);
sobol_seq_new = sobol_seq(:,2);

vd(:,1) = vdcorput(N-1,2);

X1_stream_sobol = zeros(1, N);
X2_stream_sobol = zeros(1, N);

%X1_sobol = zeros(1, N);


X1_stream_vandercorput = zeros(1, N);
X2_stream_vandercorput = zeros(1, N);

Z_CORDIV_sobol = zeros(1, N+1);
Z_CORDIV_vd = zeros(1, N+1);

% X_CORLD = zeros(2^bp,N);
% X_CORLD_vd = zeros(2^bp,N);

%for i = 1:2^bp
 %   for j = 1:2^bp
        if X < Y
            m = X;
            s = X;
            for z=1:N
                if (((Y-1)/N) > sobol_seq_new(z))
                    X2_stream_sobol(z) = 1;                    
                end
                if (((Y-1)/N) > vd(z))
                    X2_stream_vandercorput(z) = 1;
                end
            end
            for d = 1:N
                if m>0
                    X1_stream_sobol(d) = X2_stream_sobol(d);
                    if X2_stream_sobol(d) == 1
                        m = m-1;
                    end
                else
                    X1_stream_sobol(d) = 0;
                end
            end
            for d = 1:N
                if s>0
                    X1_stream_vandercorput(d) = X2_stream_vandercorput(d);
                    if X2_stream_vandercorput(d) == 1
                        s = s-1;
                    end
                else
                    X1_stream_vandercorput(d) = 0;
                end         
            end

            X1_sobol = and(X1_stream_sobol, X2_stream_sobol);
            X2_sobol = or(X1_stream_sobol, X2_stream_sobol);
% 
%             X1_vd = and(X1_stream_vandercorput, X2_stream_vandercorput);
%             X2_vd = or(X1_stream_vandercorput, X2_stream_vandercorput);

            
            for k = 1:N
                Z_CORDIV_sobol(k+1) = bitor(bitand(not(X2_stream_sobol(k)),Z_CORDIV_sobol(k)),and(X2_stream_sobol(k),X1_stream_sobol(k)));
                %Z_CORDIV_sobol(k+1) = bitor(bitand(not(X2_sobol(k)),Z_CORDIV_sobol(k)),and(X2_sobol(k),X1_sobol(k)));

                Z_CORDIV_vd(k+1) = bitor(bitand(not(X2_stream_vandercorput(k)),Z_CORDIV_vd(k)),and(X2_stream_vandercorput(k),X1_stream_vandercorput(k)));
                %Z_CORDIV_vd(k+1) = bitor(bitand(not(X2_vd(k)),Z_CORDIV_vd(k)),and(X2_vd(k),X1_vd(k)));
            end
        else    % X >= Y
            m = Y;
            s = Y;
            for z=1:N
                if (((X-1)/N) > sobol_seq_new(z))
                    X1_stream_sobol(z) = 1;                    
                end
                if (((X-1)/N) > vd(z))
                    X1_stream_vandercorput(z) = 1;
                end
            end
            for d = 1:N
                if m>0
                    X2_stream_sobol(d) = X1_stream_sobol(d);
                    if X1_stream_sobol(d) == 1
                        m = m-1;
                    end
                else
                    X2_stream_sobol(d) = 0;
                end
            end
            for d = 1:N
                if s>0
                    X2_stream_vandercorput(d) = X1_stream_vandercorput(d);
                    if X1_stream_vandercorput(d) == 1
                        s = s-1;
                    end
                else
                    X2_stream_vandercorput(d) = 0;
                end         
            end

            
            for k = 1:N
                Z_CORDIV_sobol(k+1) = bitor(bitand(not(X1_stream_sobol(k)),Z_CORDIV_sobol(k)),and(X2_stream_sobol(k),X1_stream_sobol(k)));
                %Z_CORDIV_sobol(k+1) = bitor(bitand(not(X2_sobol(k)),Z_CORDIV_sobol(k)),and(X2_sobol(k),X1_sobol(k)));

                Z_CORDIV_vd(k+1) = bitor(bitand(not(X1_stream_vandercorput(k)),Z_CORDIV_vd(k)),and(X2_stream_vandercorput(k),X1_stream_vandercorput(k)));
                %Z_CORDIV_vd(k+1) = bitor(bitand(not(X2_vd(k)),Z_CORDIV_vd(k)),and(X2_vd(k),X1_vd(k)));
            end
        end
    %end
    %end
    Z_sobol = sum(Z_CORDIV_sobol(1:end-1))/N;
    %Z_sobol = sum(Z_CORDIV_sobol)/N;
    Z_vd = sum(Z_CORDIV_vd)/N;
end
    


