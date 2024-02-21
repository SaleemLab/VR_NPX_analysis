function g = grayscale_based_on_num_input(n)


g = 0.3:(0.7/n):1;
g = repmat(g,[3 1]);
g = transpose(g);