# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Thu Oct 17 21:11:25 2024
"""

def longestDiverseString(a: int, b: int, c: int) -> str:
    r=""
    lastIDcount=0
    abc=[a,b,c]
    alpha=['a','b','c']
    while abc:
        maxID=max(range(len(abc)), key=lambda x: abc[x])
        print("abc=",abc)
        print("maxID=",maxID)
        
        if len(abc)==0:
            print('kk')
            return r
        
        elif abc[maxID]==0:
            print('2')
            del abc[maxID]
            del alpha[maxID]
            
        elif lastIDcount<2:
            print('3')
            r+=alpha[maxID]
            abc[maxID]-=1
            lastIDcount+=1
            
        else:
            temp=abc.copy()
            tempalpha=alpha.copy()
            del temp[maxID]
            del tempalpha[maxID]
            print(temp, tempalpha)
            
            if len(temp)!=0:
                print('41')
                secondMaxID=max(range(len(temp)), key=lambda x: temp[x])
                r+=alpha[secondMaxID]
                lastIDcount=1
            else:
                print('42')
                return r
                
                
    
print(longestDiverseString(a = 7, b = 1, c = 0))