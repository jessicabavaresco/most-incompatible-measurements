%   See-saw algorithm to find optimal set of general N k-outcome measurements
%   to steer  a given quantum state rho_AB. The optimal set is the one that,
%   when locally applied to state rho_AB, will generate the assemblage that
%   is most robust to generalized noise. Calculates a lower bound for the
%   generalized robustness of steering of the input quantum state when
%   subjected to arbitrary sets of N k-outcome measurements.
%   
%   authors:     Jessica Bavaresco, Marco Tulio Quintino, Leonardo Guerini,
%                Thiago O. Maciel, Daniel Cavalcanti, Marcelo Terra Cunha
%
%   requires:    Yalmip (https://yalmip.github.io) and QETLAB (http://www.qetlab.com)
%
%   last update: May, 2017

function [opt_M, gr] = seesaw_gr(rho_AB,d,N,k)
%seesaw_gr Calculates a lower bound for the generalized robustness of
%   steering gr of a quantum state rho_AB subjected to N k-outcome
%   measurements using a see-saw algorithm. This see-saw consists in the
%   iteration of two SDPs (GR_sdp1 and GR_sdp2). GR_sdp1 optimizes over
%   steering inequalitites for fixes state and measurements. GR_sdp2
%   optimizes over measurements for fixed state and steering inequality.
%   The final set of measurements is a candidate for the optimal set of N
%   k-outcome measurements to steer the state rho_AB. Next, SDP gr_ste
%   calculates the generalized robustness of steering gr of the input state
%   rho_AB when subjected to the optimal set of measurements opt_M found by 
%   the see-saw. The value gr is a lower bound for the generalized
%   robustness of steering of the state rho_AB when subjected to arbitrary
%   sets of N k-outcome measurements.
%
%   INPUT: rho_AB = bipartite quantum state
%               d = [dA dB], dimension vector, dA for Alice and dB for Bob
%               N = number of measurements
%               k = number of outcomes of each measurement      
%
%   OUTPUT: opt_M = optimal set of measurements opt_M and 
%              gr = lower bound for the generalized robustness 

dA = d(1);

% prelocates the variable M which stores the set of measurements 
% each matrix M(:,:,i,j) corresponds to effect j of measurement i
M_ax  = zeros(dA,dA,N,k);

eta_jm  = 1;
while eta_jm>0.99
      for i=1:N
          % assigns a random POVM to each measurement M(:,:,i,:) -- see RandomPOVM
          M_ax(:,:,i,:) = RandomPOVM(dA,k);
      end
      % checks if the set M is incompatible -- see jm
      % the set is incompatible if eta_jm < 1
      [~, eta_jm] = jm(M_ax);
end
% while loop ends when a set of incompatible measurements M is found

% gap is the difference between the values of gr in two subsequent iterations
gap = 1;

% when gap is smaller than precision the see-saw is halted.  
precision = 10^(-4);

count = 0;
old_J = 1;

%tic;
while gap>precision
    
    % first SDP of the see-saw optimizes over steering inequalities for
    % fixed state rho_AB and set of measurements M -- see GR_sdp1
    [F_ax, J, flag] = GR_sdp1(M_ax,rho_AB);
    
    if flag.problem~=0
        GR_sdp1_problem = flag.problem % will print flag on screen if SDP fails
        info            = flag.info % will print flag info on screen if SDP fails
    end
    
    % second SDP of the see-saw optimizes over sets of measurements for
    % fixed state rho_AB and coefficients F_ax of fized steering inequality
    % -- see GR_sdp2
    [M_ax, flag] = GR_sdp2(F_ax,rho_AB);  
    
    if flag.problem~=0
        GR_sdp2_problem = flag.problem % will print flag on screen if SDP fails
        info            = flag.info % will print flag info on screen if SDP fails
    end
    
    check_povm(M_ax); % checks if all measurements in the set M output of GR_sdp2 are valid POVMs
    
    gap     = abs(old_J - J);
    old_J   = J;
    count   = count + 1;
end

opt_M   = M_ax; % optimal set of measurements opt_M found by the see-saw
[gr, ~] = gr_ste(rho_AB,opt_M); % calculates the generalized robustness gr of the input state rho_AB and optimal set of measurements opt_M

%iterations = count; % number of iterations that were necessary to achieve gap<=precision
%time       = toc;   % time necessary to achieve gap<=precision

end

function M = RandomPOVM(dA,k)
%RandomPOVM Generates a random POVM tensor given dimension dA and number of
%   outcomes k
%
%   INPUT: dA = dimension of Hilbert space
%           k = number of outcomes in the POVM
%
%   OUTPUT: M = random k-outcome POVM in dimension dA

rng('shuffle');

T             = RandomUnitary(k);
T             = abs(T.^2);
prob          = sparse(k,1);
r             = k*rand;
prob(ceil(r)) = 1;
prob          = T*prob;

% M(:,:,i) is effect i of measurement M
M = zeros(dA,dA,k);

for i=1:k-1
    U = RandomUnitary(dA);
%     T = RandomUnitary(d);
%     T = abs(T.^2);
    p = sparse(dA,1);
    r = dA*rand;
    p(ceil(r)) = 1;
%     M(:,:,i) = U*(diag(p))*U'/no;
    M(:,:,i) = U*(diag(p))*U'*prob(i);
end

povm_ok = 0;
n       = 1;
while ~povm_ok
    if min(real(eig(eye(dA)-sum(M,3))))>-1e-6
       povm_ok = 1;
    else
        U = RandomUnitary(dA);
%         T = RandomUnitary(d);
%         T = abs(T.^2);

        p = sparse(dA,1);
        r = dA*rand;
        p(ceil(r)) = 1;
%         M(:,:,mod(n,no-1)+1) = U*(diag(T*p))*U'*prob(mod(n,no-1)+1);
        M(:,:,mod(n,k-1)+1) = U*(diag(p))*U'*prob(mod(n,k-1)+1);
        n = n + 1;
    end
end

M(:,:,k) = eye(dA) - sum(M,3);

end

function [joint_M, eta_jm] = jm(M)
%jm Is the SDP that calculates the white noise robustness of the
%   incompatibility of a set of measurements M.If eta_jm = 1 then M is
%   compatible. If eta_jm < 1 than M is incompatible.
%
%   INPUT:        M = set of measurements
%
%   OUTPUT: joint_M = joint POVM for the set of measurements M depolarized by
%            eta_jm = critical visibility

dA = size(M,1);
N  = size(M,3);
k  = size(M,4);
D  = zeros(N,k,k^N);

yalmip('clear');

sdpvar eta_jm
joint_M = sdpvar(dA,dA,k^N,'hermitian','complex'); % joint POVM

% generates all deterministic probability distributions D for N inputs and
% k possible outputs.
F = [];
for l=1:k^N
    F = F + [joint_M(:,:,l)>=0]; % positivity constraints for the effects of the joint POVM joint_M
    string = dec2base(l-1,k,N);
    for i=1:N
        c = str2double(string(i));
        D(i,c+1,l) = 1;
    end
end

for i=1:N
    for j=1:k
        cg_ax = zeros(dA,dA);
        for l=1:k^N
            cg_ax = cg_ax + D(i,j,l)*joint_M(:,:,l);
        end
        % constraint that the depolarized set of measurements can be
        % obtained via coarse-graining of the joint POVM
        F = F + [cg_ax==eta_jm*M(:,:,i,j)+((1-eta_jm)/dA)*trace(M(:,:,i,j))*eye(dA)];
    end
end

% maximizes the visibility eta_jm subjected to constraints F
SOLUTION = solvesdp(F,-eta_jm,sdpsettings('solver','mosek','verbose',0));

joint_M = double(joint_M);
eta_jm  = double(eta_jm);

end

function [F_ax, J, flag] = GR_sdp1(M_ax,rho_AB)
%GR_sdp1 Is the dual formulation of the SDP that calculates the generalized
%   robustness of steering of a quantum state rho_AB subjected to local set
%   of measurements M. It outputs the coefficients F_ax of steering
%   inequality whose value obtained by the input state and measurements is
%   precisely their generalized robustness of steering.
%
%   INPUT: rho_AB = quantum state 
%               M = set of measurements
%
%   OUTPUT:  F_ax = coefficients of optimal steering inequality 
%               J = generalized robustness
%            flag = exit flag of the SDP solution

dA = size(M_ax,1);
dB = size(rho_AB,1)/dA;
N  = size(M_ax,3);
k  = size(M_ax,4);
  
yalmip('clear');

% variables are the coefficients of a steering inequality
F_ax = sdpvar(dB,dB,N,k,'hermitian','complex');

F = [];
D = zeros(N,k,k^N); % deterministic probability distributions
for l=1:k^N
    U = zeros(dB,dB);
    string = dec2base(l-1,k,N);
    for i=1:N
        c = str2double(string(i));
        D(i,c+1,l) = 1;
        for j=1:k       
            U = U + D(i,j,l)*F_ax(:,:,i,j);
        end
    end
    F = F + [eye(dB)-U>=0]; % SDP constraint
end

tr_Fsig  = 0;

for i=1:N
    for j=1:k
        sig_ax  = PartialTrace(kron(M_ax(:,:,i,j),eye(dB))*rho_AB,1,[dA dB]);
        tr_Fsig = tr_Fsig + trace(F_ax(:,:,i,j)*sig_ax);
        F = F + [F_ax(:,:,i,j)>=0]; % SDP constraint 
    end
end

J = tr_Fsig - 1; % generalized robustness of steering

% maximizes generalized robustness of steering J subjected to constraints F
flag = solvesdp(F,-J,sdpsettings('solver','mosek','verbose',0));

F_ax = double(F_ax);
J    = double(J);

end

function [M_ax, flag] = GR_sdp2(F_ax,rho_AB)
%GR_sdp2 Finds the set M of N k-outcome of measurements that, when locally
%   applied on input state rho_AB, generates the assemblage that maximally
%   violates the steering inequality defined by the input coefficients F_ax
%
%   INPUT: rho_AB = quantum state 
%            F_ax = coefficients of steering inequality
%
%   OUTPUT:  M    = set of measurements that maximally violates the inequality
%            flag = exit flag of the SDP solution.

dB = size(F_ax,1);
dA = size(rho_AB,1)/dB;
N  = size(F_ax,3);
k  = size(F_ax,4);

yalmip('clear')

% variable is a set of measurements
M_ax = sdpvar(dA,dA,N,k,'hermitian','complex');

sum_M = sum(M_ax,4);

tr_Fsig  = 0;

F = [];
for i=1:N
    for j=1:k
        tr_Fsig = tr_Fsig + trace(F_ax(:,:,i,j)*PartialTrace(kron(M_ax(:,:,i,j),eye(dB))*rho_AB,1,[dA dB]));
        F = F + [M_ax(:,:,i,j)>=0]; % constraint that the effects of the measurements are positive semidefinite
    end
    F = F + [sum_M(:,:,i)==eye(dA)]; % constraint that the effects of the measurements sum to the identity
end

J = tr_Fsig; % expression of the steering inequality 

% maximizes the violation J of the steering inequality defined by F_ax
% subjected to the constraints that the variable M is a valid set of POVMs
flag = solvesdp(F,-J,sdpsettings('solver','mosek','verbose',0));

for i=1:N
    for j=1:k
        M_ax(:,:,i,j) = double(M_ax(:,:,i,j));
        M_ax(:,:,i,j) = (M_ax(:,:,i,j)+M_ax(:,:,i,j)')/2; 
    end
end

end

function [gr, solution] = gr_ste(rho_AB, M_ax)
%gr_ste Calculates the generalized robustness of steering gr of the quantum
%   state rho_AB subjected to local measurements M
%
%   INPUT:   rho_AB = quantum state 
%                 M = set of measurements
%  
%   OUTPUT:      gr = generalized robustness
%          solution = exit flag of the SDP solution

dA = size(M_ax,1);
dB = size(rho_AB,1)/dA;
N  = size(M_ax,3);
k  = size(M_ax,4);
D  = zeros(N,k,k^N);
   
yalmip('clear');

% variable is assemblage in the LHS model
sig_loc = sdpvar(dB,dB,k^N,'hermitian','complex'); 

% generates deterministic probability distribution
F = [];
for l=1:k^N
    F = F + [sig_loc(:,:,l)>=0]; % positivity constraint on the elements of the assemblage in the LHS model
    string = dec2base(l-1,k,N);
    for i=1:N
        c = str2double(string(i));
        D(i,c+1,l) = 1;
    end
end

for i=1:N
    for j=1:k
        sig_ax = PartialTrace(kron(M_ax(:,:,i,j),eye(dB))*rho_AB,1,[dA dB]);
        uns_ax = zeros(dB,dB);
        for l=1:k^N
            uns_ax = uns_ax + D(i,j,l)*sig_loc(:,:,l); % constructs LHS model
        end
        F = F + [uns_ax-sig_ax>=0]; % SDP constraint
    end
end

gr = trace(sum(sig_loc,3)) - 1;

% minimizes the amount of generalized noise that is necessary to apply to
% the assemblage generated by input state and measurements in order for it
% to accept a LHS model
solution = solvesdp(F, gr, sdpsettings('solver','mosek','verbose',0))

gr = double(gr);

end

function test = check_povm(M)
%check_povm tests if tensor M is a valid POVM. Tensor M, M(:,:,i) is the i 
%   effect of POVM M
%
%   INPUT:  M = POVM
%
%   OUTPUT: test==1 M is valid POVM
%           test==0 M is NOT valid POVM

dA = size(M,1);
N  = size(M,3);
k  = size(M,4);

test = 1;
for i=1:N
    for j=1:k
        eig_M = eig(M(:,:,i,j));
            for d=1:dA
                if eig_M(d)<-10^(-4) % checks positivity of the effects of the POVM
                   test = 0;
                   error('POVM NOT VALID: negative eigenvalue')
                end
            end
    end
    
    sum_M = sum(M(:,:,i,:),4);
    
    if norm(sum_M-eye(dA))>10^(-4) % checks if effects of the POVM sum to identity
        test = 0;
        error('POVM NOT VALID: effects do not sum to Id')
    end
end

end
