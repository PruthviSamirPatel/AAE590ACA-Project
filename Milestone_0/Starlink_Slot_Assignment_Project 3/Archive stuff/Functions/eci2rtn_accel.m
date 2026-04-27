function a_rtn = eci2rtn_accel(r_eci, v_eci, a_eci)
% Rotate inertial acceleration vector into RTN frame

    Rhat = r_eci / norm(r_eci);

    hvec = cross(r_eci, v_eci);
    Nhat = hvec / norm(hvec);

    That = cross(Nhat, Rhat);

    Q = [Rhat, That, Nhat];   % columns are RTN basis in inertial frame

    a_rtn = Q' * a_eci;
end