clear all;
% Initialiation
syms theta_x theta_y bg_x bg_y bg_z u_x u_y u_z n_x n_y n_z m_x m_y m_z;

% Form f
A_1=[1 sin(theta_x)*tan(theta_y) cos(theta_x)*tan(theta_y)
    0 cos(theta_x) -sin(theta_x)];

A_2=[u_x-bg_x-n_x
    u_y-bg_y-n_y
    u_z-bg_z-n_z];

f=A_1*A_2;

% State vector
X=[theta_x theta_y bg_x bg_y bg_z];

% Input vector
U=[u_x u_y u_z];

% Noise Vectors
V=[n_x n_y n_z];
W=[m_x m_y m_z];

% Form h
h = 9.8*[-sin(theta_y); cos(theta_y)*sin(theta_x); cos(theta_y)*cos(theta_x)];

% Initial State
theta_x=0; theta_y=0; bg_x=0.1; bg_y=0.1; bg_z=0.1; u_=0; u_y=0; u_z=0;
X_init=[0; 0; 0.1; 0.1; 0.1];

for i=1:150
    % Generate Noise
    n_x=random('norm',0,sqrt(0.00015)); n_y=random('norm',0,sqrt(0.00015)); n_z=random('norm',0,sqrt(0.00015));
    m_x=random('norm',0,sqrt(0.001)); m_y=random('norm',0,sqrt(0.001)); m_z=random('norm',0,sqrt(0.001));
    
    % linearize matrices
    A=eval(jacobian(f,X));
    A=[A;zeros(3,5)];
    
    B=eval(jacobian(f,U));
    B=[B;zeros(3,3)];
    
    C=eval(jacobian(h,X));
    D=eval(jacobian(h,U));
    
    G=eval(jacobian(f,V));
    G=[G;zeros(3)];
    
    H=zeros(3);
    
    % Differential Ricatti Equation
    [T, P] = ode45(@(t,P)mRiccati(t, P, A, (1/0.001)*(C')*C, 0.00015*G*G'), [0 10], eye(5));
    P=[P(end,1:5);P(end,6:10);P(end,11:15);P(end,16:20);P(end,21:25)];
    
    % Gain
    Gain=P*C'*(1/0.001);
    
    % Trajectory
    t1=i*pi/10;
    t2=i*pi/60;
    y(:,i) = 9.8*[-sin(t2)+m_x; cos(t2)*sin(t1)+m_y; cos(t1)*cos(t2)+m_z];
    
    % Calculated Input
    U_new=[pi/10+n_x+bg_x pi/60*cos(t1)+n_y+bg_y pi/60*sin(t2)+n_z+bg_z]';
    
    % Update states
    [t, X_new] = ode45(@(t,X_new)f1(t, A, B, C, Gain, U_new, X_new, y(:,i)), [0 10], [theta_x; theta_y; bg_x; bg_y; bg_z]);
    l=length(X_new);
    X1=X_new(l,:)';
    theta_x=X1(1); theta_y=X1(2); bg_x=X1(3); bg_y=X1(4); bg_z=X1(5);
    X2(i,:)=X1;
    
    % Output
    y_out(:,i)=C*X1;
    
end

%%
t=1:1:150;
figure(1)
subplot(5,1,1)
plot(t,X2(:,1),t,sawtooth(t*pi/10-pi));title('Theta_x(pitching with f=1/20)');legend('Estimated','Real');
subplot(5,1,2)
plot(t,X2(:,2),t,sawtooth(t*pi/60-pi));title('Theta_y(rolling with f=1/120)');legend('Estimated','Real');
subplot(5,1,3)
plot(t,X2(:,3),[1 150],[0.1 0.1]);title('Constant 0.1');legend('Estimated','Real');
subplot(5,1,4)
plot(t,X2(:,4),[1 150],[0.1 0.1]);title('Constant 0.1');legend('Estimated','Real');
subplot(5,1,5)
plot(t,X2(:,5),[1 150],[0.1 0.1]);title('Constant 0.1');legend('Estimated','Real');

%%
t=1:1:150;
figure(2)
subplot(3,1,1)
plot(t,y_out(1,:),t,y(1,:));title('First Output');legend('Estimated','Real');
subplot(3,1,2)
plot(t,y_out(2,:),t,y(2,:));title('Second Output');legend('Estimated','Real');
subplot(3,1,3)
plot(t,y_out(3,:),t,y(3,:));title('Third Output');legend('Estimated','Real');