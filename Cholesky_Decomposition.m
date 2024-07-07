M = [6.01 -3; -3 6.01];
N = [3.12; 2.82];

n = height(M);

R = zeros(size(M));
for i = 1:n
    sum = 0;
    for k = 1:(i-1)
        sum = sum + R(k,i)^2;
    end
    R(i,i) = sqrt(M(i,i) - sum);
    for j = (i+1):n
        sum = 0;
        for k = 1:(i-1)
            sum = sum + R(k,i)*R(k,j);
        end
        R(i,j) = (M(i,j) - sum)/R(i,i);
    end
end

Z = zeros(n,1);
for i = 1:n
    sum = 0;
    for j = 1:(i-1)
        sum = sum + R(j,i)*Z(j);
    end
    Z(i) = (N(i) - sum)/R(i,i);
end

xhat = zeros(n,1);
for i = n:-1:1
    sum = 0;
    for j = i+1:n
        sum = sum + R(i,j)*xhat(j);
    end
    xhat(i) = (Z(i) - sum)/R(i,i);
end

S = zeros(size(M));
for i = 1:n
    S(i,i) = 1/R(i,i);
end
for i = 1:n
    for j = i+1:n
        sum = 0;
        for k = i:(j-1)
            sum = sum + R(k,j)*S(i,k);
        end
        S(i,j) = -S(j,j)*sum;
    end
end

P = S*transpose(S);