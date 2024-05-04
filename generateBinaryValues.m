function binaryMatrix = generateBinaryValues(N)
    % Initialize a matrix to store the binary values
    if pow2(log2(N)) ~= N
        msg = 'N should be powers of 2!';
        error(msg)
    end
    binaryMatrix = zeros(N, log2(N));
    
    % Iterate from 0 to 1023
    for i = 0:N-1
        % Convert the decimal number to binary format with 10 bits
        binaryFormat = dec2bin(i, log2(N));
        % Convert the binary format string to an array of integers
        binaryArray = arrayfun(@(x) str2double(x), num2str(binaryFormat')');
        % Store the binary array in the matrix
        binaryMatrix(i+1, :) = binaryArray;
    end
end
