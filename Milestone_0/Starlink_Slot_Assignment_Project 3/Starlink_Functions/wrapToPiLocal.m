function ang = wrapToPiLocal(ang)
%WRAPTOPILOCAL Wrap angle to [-pi, pi).
    ang = mod(ang + pi, 2*pi) - pi;
end
