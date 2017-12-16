function NumericalIntegrationImpl()
dydt = @(t,y) 2*y;
exactFn = @(t) exp(2*t);
y_o = 1;
t_start = 0;
t_end = 1;

h_values = [0.1 0.05 0.001];
for index = 1:length(h_values)
    h = h_values(index);
    exactSoln = exactFn((0:h:1)');
    y_Euler_h = runEuler(dydt, t_end, h, y_o, t_start);
    error_Euler_h = exactSoln - y_Euler_h;
    e_Euler_by_hp1 = error_Euler_h/h;
    e_Euler_by_hp2 = error_Euler_h/h^2;
    e_Euler_by_hp3 = error_Euler_h/h^4;
    
    y_RK2_h = runRKSecondOrder(dydt, t_end, h, y_o, t_start);
    error_RK2_h = exactSoln - y_RK2_h;
    e_RK2_by_hp1 = error_RK2_h/h;
    e_RK2_by_hp2 = error_RK2_h/h^2;
    e_RK2_by_hp3 = error_RK2_h/h^4;
    
    y_AB4_h = ABFourthOrder(dydt, t_end, h, y_o, t_start);
    error_AB4_h = exactSoln - y_AB4_h;
    e_AB4_by_hp1 = error_AB4_h/h;
    e_AB4_by_hp2 = error_AB4_h/h^2;
    e_AB4_by_hp3 = error_AB4_h/h^4;
    
    T = table(exactSoln, y_Euler_h, error_Euler_h, e_Euler_by_hp1, e_Euler_by_hp2, e_Euler_by_hp3,...
        y_RK2_h, error_RK2_h, e_RK2_by_hp1, e_RK2_by_hp2, e_RK2_by_hp3,...
        y_AB4_h, error_AB4_h, e_AB4_by_hp1, e_AB4_by_hp2, e_AB4_by_hp3);
    
    filename = 'P1_table.xlsx';
    writetable(T,filename,'Sheet',index,'Range','D1')
    
end
end

function  y_soln = runEuler(dydt, t_end, h, y_o, t_start)
steps = floor((t_end - t_start)/h);
y_soln = zeros(steps + 1, 1);
y_soln(1) = y_o;
tn = t_start;
for idx = 2:steps+1
    y_soln(idx) = y_soln(idx - 1) + h*dydt(tn, y_soln(idx - 1));
    tn = tn + h;
end
end

function  y_soln = runRKSecondOrder(dydt, t_end, h, y_o, t_start)
%a = 1/2, b = 0, c = 1
steps = floor((t_end - t_start)/h);
y_soln = zeros(steps + 1, 1);
y_soln(1) = y_o;
tn = t_start;
for idx = 2:steps+1
    k1 = dydt(tn, y_soln(idx - 1));
    k2 = dydt(tn + h/2, y_soln(idx - 1) + (h/2)*k1);
    y_soln(idx) = y_soln(idx - 1) + h*k2;
    tn = tn + h;
end
end

function  y_soln = ABFourthOrder(dydt, t_end, h, y_o, t_start)
steps = floor((t_end - t_start)/h);
y_soln = zeros(steps + 1, 1);
y_soln(1) = y_o;
tn = t_start;
%using fourth order RK to start method and find first 4 solutions
for idx = 2:4
    k1 = dydt(tn, y_soln(idx - 1));
    k2 = dydt(tn + h/2, y_soln(idx - 1) + h/2*k1);
    k3 = dydt(tn + h/2, y_soln(idx - 1) + h/2*k2);
    k4 = dydt(tn + h, y_soln(idx - 1) + h*k3);
    y_soln(idx) = y_soln(idx - 1) + h/6*(k1 + 2*k2 + 2*k3 + k4);
    tn = tn + h;
end
%AB 4th order
for idx = 5:steps+1
    y_soln(idx) = y_soln(idx - 1) + h/24*( 55*dydt(tn, y_soln(idx - 1)) ...
        - 59*dydt(tn, y_soln(idx - 2)) + 37*dydt(tn, y_soln(idx - 3))...
        - 9*dydt(tn, y_soln(idx - 4)));
    tn = tn + h;
end
end