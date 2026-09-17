function [in_mat] = in_range(matrix,range)
    in_mat = matrix>=range(1) & matrix<=range(2);
end