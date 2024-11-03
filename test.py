lis = [1,2,3,7,8,9]
firstcoord = lis.index(3)+1
lastcoord = lis.index(7)
sexlist = [4,5,6]
lis[firstcoord:firstcoord] = sexlist
print(lis)