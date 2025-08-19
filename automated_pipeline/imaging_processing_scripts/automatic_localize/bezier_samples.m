function [t,B] = bezier_samples(E, C, D, M)
t = linspace(0,1,M).';
B = (1-t).^2 .* E + 2*(1-t).*t .* C + t.^2 .* D;
end