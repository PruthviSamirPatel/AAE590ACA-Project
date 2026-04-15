function C_eci2rtn = eci2rtn_dcm(r, v)

rhat = r / norm(r);
h = cross(r, v);
hhat = h / norm(h);
that = cross(hhat, rhat);

C_eci2rtn = [rhat.'; that.'; hhat.'];
end