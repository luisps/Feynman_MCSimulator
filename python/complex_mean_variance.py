def c_mean (zList):
    M = len(zList)
    m = sum (z for z in zList)/ M
    return m

def c_var (zList):
    M = len (zList)
    m = c_mean (zList)
    mReal = m.real
    mIm = m.imag 
    varReal = sum((z.real - mReal) ** 2 for z in zList) / M 
    varIm = sum((z.imag - mIm) ** 2 for z in zList) / M
    var = varReal + varIm
    return var, varReal, varIm

def c_var_true_mean (zList, true_mean):
    M = len (zList)
    varReal = sum((z.real - true_mean.real) ** 2 for z in zList) / M 
    varIm = sum((z.imag - true_mean.imag) ** 2 for z in zList) / M
    var = varReal + varIm
    return var, varReal, varIm

    