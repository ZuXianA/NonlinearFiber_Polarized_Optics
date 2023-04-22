clear; clc;

viz = 2;    % 题号选择，记得修改保存文件名

s1 = [1 0 0];
s2 = [0 1 0];
s3 = [0 0 1];

if viz == 1
    qwps1 = SingleWaveplate(s3, 359, 'qwp');
    PolarizationTraceShow(qwps1', s3);
elseif viz == 2
    qwps1 = SingleWaveplate(s1, 359, 'hwp');
    PolarizationTraceShow(qwps1', s1);
elseif viz == 3
    qwp_hwp = TwoWaveplate(s3, 179);
    PolarizationTraceShow(qwp_hwp', s3);
elseif viz == 4
    qwp_hwp_qwp = TwoWaveplate(s1, 179);
    PolarizationTraceShow(qwp_hwp_qwp', s1)
else
    error('The number of viz is error.')
end


function rot = ROT(theta)  % 定义一个旋转器，庞加莱球上是绕着 S3 右旋 2*theta 角
rot = [cos(deg2rad(2*theta)), -sin(deg2rad(2*theta)), 0;...
    sin(deg2rad(2*theta)), cos(deg2rad(2*theta)), 0;...
    0, 0, 1];
end

function rotinv = ROTINV(mat)   % 定义一个旋转器的逆，用于计算米勒矩阵 R
rotinv = inv(mat);
end

function rwp = RWP(delta)   % 定义零方位角相位延迟器，庞加莱球上是绕着 S1 右旋 delta 角
rwp = [1, 0, 0;...
    0, cos(deg2rad(delta)), -sin(deg2rad(delta));...
    0, sin(deg2rad(delta)), cos(deg2rad(delta))];
end


function out = SingleWaveplate(state, max_theta, waveplate)
if strcmpi(waveplate, 'qwp')
    wp = RWP(90);
elseif strcmpi(waveplate, 'hwp')
    wp = RWP(180);
else
    error('Error, please re-input waveplate.')
end

out = [];
theta_range = 0:1:max_theta;

for i = theta_range
    rot = ROT(i);
    res = rot * wp * ROTINV(rot);  % R_wp,theta = R_theta * R_wp,0 * R_thetainv
    res(abs(res) < 1.3e-16) = 0;

    out = [out, res * state'];  % OUT = R_wp,theta * INPUT
end
end


function out = TwoWaveplate(state, max_theta)

out = [];
qwp = RWP(90);      % 1/4
hwp = RWP(180);     % 1/2
theta_range = 0:1:max_theta;

for i = theta_range
        rot1 = ROT(i);
        for j = theta_range
            rot2 = ROT(j);
            % 多个波片的级联，实际上应该是矩阵左乘
            % R_wp,theta = R_wp,theta1 * R_wp,theta0
            res = rot2 * hwp * ROTINV(rot2) * rot1 * qwp * ROTINV(rot1);
            res(abs(res) < 1.3e-16) = 0;

            out = [out, res * state'];  % OUT = R_wp,theta * INPUT
        end
end
end


function out = ThreeWaveplate(state, max_theta)

out = [];
qwp1 = RWP(90); % 1/4 波片
hwp = RWP(180); % 1/2 波片
qwp2 = RWP(90); % 1/4 波片
theta_range = 0:1:max_theta;

for i = theta_range
    rot1 = ROT(i);
    for j = theta_range
        rot2 = ROT(j);
        for k = theta_range
            rot3 = ROT(k);
            % 多个波片的级联，实际上应该是矩阵左乘
            % R_wp,theta = R_wp,theta2 * R_wp,theta1 * R_wp,theta0
            res = rot3 * qwp2 * ROTINV(rot3) * rot2 * hwp * ROTINV(rot2) * rot1 * qwp1 * ROTINV(rot1);
            res(abs(res) < 1.3e-16) = 0;

            out = [out, res * state'];  % OUT = R_wp,theta * INPUT
        end
    end
end
end


