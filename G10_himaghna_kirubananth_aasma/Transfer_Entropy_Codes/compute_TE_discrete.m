function TE = compute_TE_discrete(Yf, Yp_idx, Xp_idx)

joint_xyz = accumarray([Yf, Yp_idx, Xp_idx], 1);
joint_yz  = accumarray([Yp_idx, Xp_idx], 1);
joint_xy  = accumarray([Yf, Yp_idx], 1);
joint_y   = accumarray(Yp_idx, 1);

P_xyz = joint_xyz / sum(joint_xyz(:));
P_yz  = joint_yz / sum(joint_yz(:));
P_xy  = joint_xy / sum(joint_xy(:));
P_y   = joint_y / sum(joint_y(:));

TE = 0;

for i = 1:numel(P_xyz)
    if P_xyz(i) > 0
        [a,b,c] = ind2sub(size(P_xyz), i);
        p = P_xyz(i);
        p1 = p / P_yz(b,c);
        p2 = P_xy(a,b) / P_y(b);
        TE = TE + p * log2(p1/p2);
    end
end

end