function ehat = myFilter(A, C, y)
    ehat = filter(A, C, y);
    ehat = ehat(length(A):end);
end