def nloop(arr, n):
    if n == 0:
        yield tuple()
    else:
        for i in arr:
            for l in nloop(arr,n-1):
                yield (i,) + l
