%Implmentation of iterative numerical methods: Jacobi, Gauss-Sidel and SOR
function [jacobiSoln, gsSoln, sorSoln] = iterativeSystemOfEquationsSolver(A, B)
[m1, n1] = size(A);
[m2, n2] = size(B);
if (m1 > 0 && m2 > 0 && m1 == n1 && m1 == m2)
    error = 10^-5;
    Xo = zeros(size(B));
    X_actual = [1;1];
    [jacobiSoln, jerrorVectors] = solveByJacobi(A, B, Xo, error);
    [gsSoln, gserrorVectors] = solveByGS(A, B, Xo, error);
    [sorSoln, sorerrorVectors] = solveBySOR(A, B, Xo, error);
    [em,en] = size(jerrorVectors);
    jerrorVectors = jerrorVectors - repmat(X_actual, [1, en]);
    jerrorVectors = max(jerrorVectors,[],1); 
    
    n = zeros(1,en);
    for id = 1:en 
        n(id) = norm(jerrorVectors(:,id));
    end
    n
    
    [em,en] = size(gserrorVectors);
    gserrorVectors = gserrorVectors - repmat(X_actual, [1, en]);
    gserrorVectors = max(gserrorVectors,[],1);
    [em,en] = size(sorerrorVectors);
    sorerrorVectors = sorerrorVectors - repmat(X_actual, [1, en]);
    sorerrorVectors = max(sorerrorVectors,[],1);
    plot(1:1:length(jerrorVectors), log(jerrorVectors), '-o', ...
         1:1:length(gserrorVectors), log(gserrorVectors), '-x',...
         1:1:length(sorerrorVectors), log(sorerrorVectors), '-.');
     legend('Jacobi errors', 'GS errors', 'SOR errors');
     xlabel('step count');
     ylabel('logarithmic error');
else
    disp('Please provide corect input for A and B');
end
end

%Jacobi Method
function [X, eVector] = solveByJacobi(A, B, Xk_minus1, error)
    eVector = Xk_minus1;
    failSafeLoopLimit = 10^6;
    Xk = zeros(size(Xk_minus1));
    stepCount = 1;
    while stepCount <= failSafeLoopLimit
    for i = 1:length(Xk_minus1)
        sum = 0;
        for j = 1:length(Xk_minus1)
            if j ~= i
                sum = sum + A(i, j)*Xk_minus1(j);
            end
        end
        Xk(i) = (-1/A(i, i))*sum + B(i)/A(i, i);
    end
    eVector = [eVector Xk];
    if norm(Xk - Xk_minus1, Inf) < error
        break;
    end
    Xk_minus1 = Xk;
    stepCount = stepCount + 1;
    end
    X = Xk;
end

%Gauss Seidel
function [X, eVector] = solveByGS(A, B, Xk_minus1, error)
    eVector = Xk_minus1;
    failSafeLoopLimit = 10^6;
    Xk = zeros(size(Xk_minus1));
    stepCount = 1;
    while stepCount <= failSafeLoopLimit
    for i = 1:length(Xk_minus1)
        sumL = 0;
        sumU = 0;
        for j = 1:length(Xk_minus1)
            if j ~= i && j < i
                sumL = sumL + A(i, j)*Xk(j);
            end
            if j ~= i && j > i
                sumU = sumU + A(i, j)*Xk_minus1(j);
            end
        end
        Xk(i) = (1/A(i, i))*(-sumL -sumU + B(i));
    end
    eVector = [eVector Xk];
    if norm(Xk - Xk_minus1, Inf) < error
        break;
    end
    Xk_minus1 = Xk;
    stepCount = stepCount + 1;
    end
    X = Xk;
end

%SOR
function [X, eVector] = solveBySOR(A, B, Xk_minus1, error)
    eVector = Xk_minus1;
    % optimum value of W as found previously
    w = 1.0170;
    failSafeLoopLimit = 10^6;
    Xk = zeros(size(Xk_minus1));
    stepCount = 1;
    while stepCount <= failSafeLoopLimit
    for i = 1:length(Xk_minus1)
        sumL = 0;
        sumU = 0;
        for j = 1:length(Xk_minus1)    
            if j ~= i && j < i
                sumL = sumL + A(i, j)*Xk(j);
            end
            if j ~= i && j > i
                sumU = sumU + A(i, j)*Xk_minus1(j);
            end
        end
        Xk(i) = (1/A(i, i))*(-w*sumL + (1 - w)*A(i, i)*Xk_minus1(i)- w*sumU + w*B(i));
    end
    eVector = [eVector Xk];
    if norm(Xk - Xk_minus1, Inf) < error
        break;
    end
    Xk_minus1 = Xk;
    stepCount = stepCount + 1;
    end
    X = Xk;
end