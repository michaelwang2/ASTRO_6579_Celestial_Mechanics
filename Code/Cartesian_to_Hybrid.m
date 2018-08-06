function z0 = Cartesian_to_Hybrid(rvec, vvec, mu)
% [a0, e0, Omega0, I0, w0, nu0] = pos_vel2orb_ele(rvec, vvec, mu);
[a0 ,e0 ,E0 ,I0 , w0 , Omega0 ,P , tau ,A , B] = vec2orbElem (rvec, vvec, mu);
nu0 = 2*atan2(sqrt((1+e0)/(1-e0))*tan(E0/2), 1);
r0 = norm(rvec);
rd0 = vvec'*(rvec./norm(rvec));
h0 = norm(cross(rvec, vvec));
th0 = w0 + nu0;
z0 = [r0; rd0; h0; th0; I0; Omega0];
end