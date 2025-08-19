function out = RAS_function(i,j,k,Avox2ras0)
    M = [i(:) j(:) k(:) ones(numel(i),1)] * Avox2ras0.';
    out = M(:,1:3);
end