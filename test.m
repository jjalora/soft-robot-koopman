classdef test
    %data: This class stores function that can be used to mess with
    %experimental and simulated data files with fields t , y , u , (x) ,...
    %   Detailed explanation goes here
    
    properties
        Property1
    end
    
    methods
        function obj = data(varargin)
            %data: Construct an instance of this class
            %   Detailed explanation goes here
%             obj.Property1 = inputArg1 + inputArg2;
        end
    end
   
    methods(Static)
        function G_opt = solve_G_QP(U, X)
            % Function to solve the QP problem for minimizing ||U - GX||^2
        
            % Number of rows and columns
            [m, ~] = size(U); 
            [n, ~] = size(X);
            % Helper function to vectorize a matrix
            vec = @(M) M(:);
            
            % Formulate the QP problem
            
            % Hessian
            H = kron(X * X', eye(m)); 
            
            % Linear term
            f = -2 * vec(U * X');  
            
            % No equality or inequality constraints
            A = [];
            b = [];
            Aeq = [];
            beq = [];
            
            % Solve the QP problem
            options = optimoptions('quadprog', 'Display', 'off');
            G_vec = quadprog(H, f, A, b, Aeq, beq, [], [], [], options);
            
            % Reshape the solution back into matrix form
            G_opt = reshape(G_vec, m, n);
        end
    end
end