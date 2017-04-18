
def family_predict(alist):
    cnum = 0
    enum = 0
    hnum = 0
    for i in alist:
        if i == 'C':
            cnum += 1
        elif i == 'E':
            enum += 1
        else:
            hnum += 1

    return (hnum - enum) / (len(alist) - cnum)


def family_predict1(alist,window):

    result1 = ''
    for i in range(len(alist)):
        start = i- window if (i-window)>0 else 0
        end = i+ window if (i+window)<len(alist) else len(alist)
        t1 = alist[start:end].count('H')
        t2 = alist[start:end].count('E')
        t3 = alist[start:end].count('C')

        if t1 == max(t1,t2,t3):
            result1 += 'H'
        elif t2 == max(t1,t2,t3):
            result1 += 'E'
        else:
            result1 += 'C'

    print(result1)

    return family_predict(result1)


family_predict1('CEEHEEEEECCCHCCCCCCCCCCCCCCCCEEEHHCCE',1)