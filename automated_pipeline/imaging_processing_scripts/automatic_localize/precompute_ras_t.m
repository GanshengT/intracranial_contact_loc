function [Xras, t_field] = precompute_ras_t(Avox2ras0, start_traj, U, M,N,P)
    [I,J,K] = ndgrid(1:M,1:N,1:P);
    X = Avox2ras0(1:3,1:3)*[J(:)'-1; I(:)'-1; K(:)'-1] + Avox2ras0(1:3,4);
    X = reshape(X.',[M,N,P,3]);
    Xras = X;
    dX = X - reshape(start_traj,1,1,1,3);
    t_field = dX(:,:,:,1)*U(1) + dX(:,:,:,2)*U(2) + dX(:,:,:,3)*U(3);
end