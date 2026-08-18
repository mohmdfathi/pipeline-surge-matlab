function result = moc_solver(g,pipe,Q0,UpstreamHead0,DownstreamHead0, ...
    fEquivalent,headAfterValve,waveSpeed,closureTime,gridSize,endTime,valveLaw)
%MOC_SOLVER The function calculates transient pressure head and discharge after
%   valve closure using the Method of Characteristics (MOC).
%
%   Inputs:
%       g                     Gravitational acceleration [m/s^2]
%       pipe                  Structure containing pipe.L, pipe.D and pipe.area
%       Q0                    Initial steady discharge [m^3/s]
%       UpstreamHead0         Initial head at the upstream pipe end [m]
%       DownstreamHead0       Initial head immediately before the valve [m]
%       fEquivalent           Equivalent Darcy friction factor [-]
%       headAfterValve        Constant head immediately after the valve [m]
%       waveSpeed             Pressure-wave speed [m/s]
%       closureTime           Valve closure time [s]
%       gridSize              Number of computational pipe reaches [-]
%       endTime               Simulation end time [s]
%       valveLaw              Function handle defining relative valve opening
%
%   Output:
%       result.head           Pressure head at each node and time step [m]
%       result.flow           Discharge at each node and time step [m^3/s]
%       result.t              Time vector [s]
%       result.x              Pipe position vector [m]

L = pipe.L;
D = pipe.D;
A = pipe.area;

N = gridSize;

dx = L / N;
dt = dx / waveSpeed;

time = 0:dt:endTime;
if time(end) < endTime
    time(end+1) = time(end) + dt;
end
nt = numel(time);

% MOC characteristic and friction coefficients
B = waveSpeed / (g*A);
R = fEquivalent*dx / (2*g*D*A^2);

H = zeros(N+1,nt);
Q = zeros(N+1,nt);

% Initial steady-state head and discharge distributions
H(:,1) = linspace(UpstreamHead0,DownstreamHead0,N+1).';
Q(:,1) = Q0;

% Fully open valve coefficient from the initial operating condition
headAcrossValve0 = DownstreamHead0 - headAfterValve;
Cv = Q0 / sqrt(headAcrossValve0);

for n = 1:nt-1

    % Interior-node characteristic equations
    Cp = H(1:N-1,n) + B*Q(1:N-1,n) ...
        - R*Q(1:N-1,n).*abs(Q(1:N-1,n));

    Cm = H(3:N+1,n) - B*Q(3:N+1,n) ...
        + R*Q(3:N+1,n).*abs(Q(3:N+1,n));

    H(2:N,n+1) = 0.5*(Cp + Cm);
    Q(2:N,n+1) = (Cp - Cm) / (2*B);

    % Upstream constant-head reservoir boundary
    CmUp = H(2,n) - B*Q(2,n) + R*Q(2,n)*abs(Q(2,n));
    H(1,n+1) = UpstreamHead0;
    Q(1,n+1) = (UpstreamHead0 - CmUp) / B;

    % Downstream valve boundary
    CpDown = H(N,n) + B*Q(N,n) - R*Q(N,n)*abs(Q(N,n));

    if closureTime == 0
        tau = 0;
    else
        tau = valveLaw(time(n+1)/closureTime);
    end

    % With y=sqrt(H-Hd), compatibility gives
    % y^2 + B*tau*Cv*y + Hd-Cp = 0.
    b = B*tau*Cv;
    y = 0.5*(-b + sqrt(b^2 + 4*max(CpDown-headAfterValve,0)));

    Q(N+1,n+1) = tau*Cv*y;
    H(N+1,n+1) = CpDown - B*Q(N+1,n+1);
end

result.head = H;
result.flow = Q;
result.t = time;
result.x = linspace(0,pipe.L,N+1).';

end