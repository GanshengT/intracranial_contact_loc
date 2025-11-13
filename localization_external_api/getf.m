function v = getf(S, name, dflt)
if isstruct(S) && isfield(S, name) && ~isempty(S.(name))
    v = S.(name);
else
    v = dflt;
end
end