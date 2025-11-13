function dst_mask = resample_mask_nn(src_mask, A_src, dst_size, A_dst)
    % Nearest-neighbor resample of a logical 3D mask from src grid to dst grid via RAS.
    % A_src, A_dst: voxel(0-based)->RAS (4x4). dst_size = [M N P].
    
    % Build all dst voxel centers (1-based), convert to 0-based for A_* usage
    [M,N,P] = deal(dst_size(1), dst_size(2), dst_size(3));
    [j_dst,i_dst,k_dst] = meshgrid(1:N, 1:M, 1:P);   % note: MATLAB order (row=i, col=j, slice=k)
    
    % dst (1-based) -> 0-based
    V0_dst = [j_dst(:)-1, i_dst(:)-1, k_dst(:)-1, ones(numel(j_dst),1)];
    
    % 0-based voxel -> RAS
    RAS = (V0_dst * A_dst.').';         % 4 x numel
    RAS = RAS(1:3,:).';                 % Nx3
    
    % RAS -> src 0-based voxel
    V0_src_h = (A_src \ [RAS, ones(size(RAS,1),1)].').';  % Nx4
    V0_src   = V0_src_h(:,1:3);
    
    % back to src 1-based voxel indices (nearest neighbor)
    I_src = round(V0_src(:,2) + 1);
    J_src = round(V0_src(:,1) + 1);
    K_src = round(V0_src(:,3) + 1);
    
    % clamp to bounds
    [Ms, Ns, Ps] = size(src_mask);
    I_src = min(max(I_src,1), Ms);
    J_src = min(max(J_src,1), Ns);
    K_src = min(max(K_src,1), Ps);
    
    % sample
    lin = sub2ind([Ms,Ns,Ps], I_src, J_src, K_src);
    dst_mask = reshape(src_mask(lin), [M,N,P]);
end