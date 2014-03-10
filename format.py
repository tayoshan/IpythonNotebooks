def checkPD(data):
    dtype = type(data)
    if dtype == 'pandas.core.frame.DataFrame':
        return True
    else:
        return False