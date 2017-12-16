function [X] = solveByPartialPivotAndGEM(A,B)
[m,n] = size(A);
[o,p] = size(B);
%strict check for correct input of n*n linear system
if m==n && m>0 && m==o && p==1
    X = performGEWithPP(A,B);
else
    disp('Provide proper inputs!');
end
end

function [X] = performGEWithPP(A, B)
X = zeros(size(B));
n = length(A);
step = 1;
%Change A to upper diagonal matrix
while step < n
    [A,B] = getPivotedMatrix(A, B, step);
    row = step + 1;
    while row <= n
        m = A(row, step)/A(step, step);
        col = step + 1;
        while col <= n
            A(row, col) = A(row, col) - m*A(step, col);
            col = col + 1;
        end
        B(row) = B(row) - B(step)*m;
        row = row + 1;
    end
    step = step + 1;
end
%Backward substitution to find X
X(n) = B(n)/A(n,n);
i = n-1;
while i >=1
    sum = 0;
    j = n;
    while j >= (i+1)
        sum = sum + A(i,j)*X(j);
        j = j - 1;
    end
    X(i) = (B(i) - sum)/A(i, i);
    i = i - 1;
end
end

function [Ap, Bp] = getPivotedMatrix(A, B, diagIndex)
[m,n] = size(A);
[o,p] = size(B);
if m > 0 && n > 0 && diagIndex > 0 && diagIndex <= n && m==o
    Ap = A;
    Bp = B;
    rowIndexWithMaxElemValue = diagIndex;
    rowIndex = diagIndex + 1;
    while rowIndex <= m
        if abs(Ap(rowIndexWithMaxElemValue, diagIndex)) < abs(Ap(rowIndex, diagIndex))
            rowIndexWithMaxElemValue = rowIndex;
        end
        rowIndex = rowIndex + 1;
    end
    %swaps row_diagIndex of matrix with row in which highest element is found
    Ap([rowIndexWithMaxElemValue diagIndex],:) = Ap([diagIndex rowIndexWithMaxElemValue],:);
    Bp([rowIndexWithMaxElemValue diagIndex],:) = Bp([diagIndex rowIndexWithMaxElemValue],:);
else
    Ap = [];
    Bp = [];
end
end