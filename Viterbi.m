function [x_opt] = Viterbi(N, s1, s_1, r_n)

%-------------------------------------------------------------------------%
%                                 FORWARD                                 %
%-------------------------------------------------------------------------%

pointer_pi(1) = -1;  
pointer_0(1) = -1;

for n=1:1:N
    if n==1
        w_3pi2(n) = real(r_n(:,n)'*s_1*exp(1i*0)); % This is for the 1st step.  
        w_pi2(n) = real(r_n(:,n)'*s1*exp(1i*0));   
    elseif n~=1
        if mod(n,2)==0 % Even bits
            % Even symbols can end up ONLY with phase pi or 0
            
            % From 3pi/2 to pi with symbol -1
            % From 3pi/2 to 0 with symbol +1
            w3pi2_pi(n) = real(r_n(:,n)'*s_1*exp(1i*3*pi/2));
            w3pi2_0(n) = real(r_n(:,n)'*s1*exp(1i*3*pi/2));
            
            % From pi/2 to pi with symbol +1
            % From pi/2 to 0 with symbol -1
            wpi2_pi(n) = real(r_n(:,n)'*s1*exp(1i*pi/2));
            wpi2_0(n) = real(r_n(:,n)'*s_1*exp(1i*pi/2));

            % The cost may be the weight(3pi/2, pi) + the weight of the
            % last symbol 3pi/2 due to the memory property of the phase.
            % The cost may be the weight(pi/2, pi) + the weight of the
            % last symbol pi/2 due to the memory property of the phase.
            total_cost_1 = w3pi2_pi(n) + w_3pi2(n-1);
            total_cost_2 = wpi2_pi(n) + w_pi2(n-1);
            [w_pi(n),pointer_pi(n)] = max([total_cost_1 0 total_cost_2 0]);
            
            % The cost may be the weight(3pi/2, 0) + the weight of the
            % last symbol 0 due to the memory property of the phase.
            % The cost may be the weight(pi/2, 0) + the weight of the
            % last symbol 0 due to the memory property of the phase.
            total_cost_1 = w3pi2_0(n) + w_3pi2(n-1);
            total_cost_2 = wpi2_0(n) + w_pi2(n-1);
            [w_0(n),pointer_0(n)] = max([total_cost_1 0 total_cost_2 0]);
            
        elseif mod(n,2)~=0 % Odd bits
            % Odd symbols can end up ONLY with phase 3pi/2 or pi/2
            
            % From 0 to 3pi/2 with symbol -1
            % From 0 to pi/2 with symbol +1
            w0_3pi2(n) = real(r_n(:,n)'*s_1);
            w0_pi2(n) = real(r_n(:,n)'*s1);
            
            % From pi to 3pi/2 with symbol +1
            % From pi to pi/2 with symbol -1
            wpi_3pi2(n) = real(r_n(:,n)'*s1*exp(1j*pi));
            wpi_pi2(n) = real(r_n(:,n)'*s_1*exp(1j*pi));
            
            % The cost may be the weight(pi, 3pi/2) + the weight of the
            % last symbol pi due to the memory property of the phase.
            % The cost may be the weight(0, 3pi/2) + the weight of the
            % last symbol 0 due to the memory property of the phase.
            total_cost_1 = wpi_3pi2(n) + w_pi(n-1);
            total_cost_2 = w0_3pi2(n) + w_0(n-1);
            [w_3pi2(n),pointer_3pi2(n)] = max([0 total_cost_1 0 total_cost_2]);

            % The cost may be the weight(pi, pi/2) + the weight of the
            % last symbol pi due to the memory property of the phase.
            % The cost may be the weight(0, pi/2) + the weight of the
            % last symbol 0 due to the memory property of the phase.
            total_cost_1 = wpi_pi2(n) + w_pi(n-1);
            total_cost_2 = w0_pi2(n) + w_0(n-1);
            [w_pi2(n),pointer_pi2(n)] = max([0 total_cost_1 0 total_cost_2]);
        end
    end
end

%-------------------------------------------------------------------------%
%                                 BACKWARD                                %
%-------------------------------------------------------------------------%

% Keeping ONLY the route that gives the maximum weight sum
if mod(n,2)~=0
    [~,route(N+1)] = max([w_3pi2(n) 0 w_pi2(n) 0]);
elseif mod(n,2)==0
    [~,route(N+1)] = max([0 w_pi(n) 0 w_0(n)]);
end

route(1) = 4;
for n=N:-1:1
    if n~=1
        if mod(n,2)==0 % Even symbols
            [~,p] = max([0 w_pi(n) 0 w_0(n)]);
            enter = [0 pointer_pi(n) 0 pointer_0(n)];
            route(n) = enter(p);
        elseif mod(n,2)~=0 % Odd symbols
            [~,p] = max([w_3pi2(n) 0 w_pi2(n) 0]);
            enter = [pointer_3pi2(n) 0 pointer_pi2(n) 0];
            route(n) = enter(p);
        end
    end
    route;
    
    % Restoring the symbols 
    if route(n)-route(n+1)==-1
        x_opt(n) = -1;
    elseif route(n)-route(n+1)==1
        x_opt(n) = 1;
    elseif route(n)-route(n+1)==-3
        x_opt(n) = 1;
    elseif route(n)-route(n+1)==3
        x_opt(n) = -1;
    end  
    
end

end