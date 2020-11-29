function [out] = tal2mni(in)
    trans = [1/0.88 0 0 0.8/0.88;1/0.97 0.97 0 3.32/0.97;0 -0.05/0.88 1/0.88 0.44/0.88;0 0 0 1];
    tmp = trans * [in 1]';
    out = tmp(1:3)';
end