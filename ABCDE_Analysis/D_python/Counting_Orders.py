# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 02:31:41 2023

@author: Shihab
"""
import matplotlib.pyplot as plt

one={1: 260.75, 2: 39.0, 3: 184.75, 4: 325.75, 5: 329.5, 6: 317.25, 7: 218.25, 8: 136.0, 9: 87.5, 10: 355.5, 11: 431.25, 12: 628.25, 13: 298.75, 14: 377.0, 15: 305.75, 16: 587.5, 17: 702.0, 18: 688.0, 19: 421.0, 20: 576.25, 21: 150.75, 22: 1000, 23: 315.5, 24: 112.5, 25: 190.75, 26: 477.75, 27: 386.25, 28: 803.25, 29: 119.0, 30: 163.25, 31: 608.25, 32: 282.5, 33: 408.0, 34: 338.5, 35: 568.5, 36: 443.5, 37: 355.0, 38: 185.0, 39: 323.75, 40: 112.5, 41: 39.5, 42: 571.0, 43: 227.75, 44: 426.75, 45: 33.25, 46: 334.0, 47: 24.75, 48: 242.75, 49: 376.75, 50: 624.75}

two={1: 1000, 2: 1000, 3: 713.75, 4: 65.25, 5: 1000, 6: 802.5, 7: 198.0, 8: 450.5, 9: 0, 10: 416.75, 11: 553.75, 12: 588.25, 13: 448.25, 14: 548.25, 15: 744.0, 16: 295.25, 17: 340.25, 18: 575.25, 19: 594.75, 20: 1000, 21: 1000, 22: 1000, 23: 26.75, 24: 332.0, 25: 640.5, 26: 718.75, 27: 1000, 28: 959.25, 29: 137.75, 30: 15.0, 31: 243.75, 32: 1000, 33: 210.0, 34: 871.5, 35: 585.5, 36: 1000, 37: 1000, 38: 130.0, 39: 614.5, 40: 1000, 41: 1000, 42: 1000, 43: 255.0, 44: 1000, 45: 290.75, 46: 180.25, 47: 530.25, 48: 132.0, 49: 8.25, 50: 520.0}

three= {1: 958.0, 2: 470.0, 3: 0, 4: 0, 5: 827.5, 6: 0, 7: 0, 8: 545.0, 9: 802.5, 10: 425.0, 11: 398.0, 12: 426.75, 13: 198.5, 14: 588.75, 15: 548.25, 16: 742.0, 17: 0, 18: 296.0, 19: 254.5, 20: 337.75, 21: 638.75, 22: 570.25, 23: 447.0, 24: 26.75, 25: 0, 26: 715.75, 27: 210.0, 28: 136.75, 29: 16.0, 30: 246.5, 31: 0, 32: 63.5, 33: 0, 34: 0, 35: 588.0, 36: 425.5, 37: 0, 38: 126.25, 39: 611.25, 40: 719.0, 41: 0, 42: 299.25, 43: 179.5, 44: 130.75, 45: 543.0, 46: 10.5, 47: 0, 48: 0, 49: 0, 50: 0}

four={1: 952.5, 2: 470.0, 3: 1000, 4: 1000, 5: 828.0, 6: 1000, 7: 1000, 8: 545.5, 9: 803.25, 10: 425.5, 11: 398.5, 12: 416.5, 13: 197.0, 14: 586.5, 15: 548.0, 16: 740.25, 17: 1000, 18: 294.5, 19: 254.5, 20: 336.5, 21: 639.25, 22: 570.75, 23: 446.0, 24: 24.0, 25: 1000, 26: 713, 27: 210.25, 28: 137.0, 29: 14.25, 30: 243.75, 31: 1000, 32: 63.75, 33: 1000, 34: 1000, 35: 588.75, 36: 374.75, 37: 1000, 38: 126.75, 39: 610.25, 40: 719.25, 41: 1000, 42: 294.5, 43: 176.75, 44: 131.5, 45: 530.5, 46: 9.0, 47: 1000, 48: 1000, 49: 1000, 50: 1000}

bond = [one,two,three,four]

for key in range(1,51):
    for i in range(4):
        if key in bond[i] and bond[i][key] in [0,1000]:
            bond[0].pop(key)
            bond[1].pop(key)
            bond[2].pop(key)
            bond[3].pop(key)
            
for x in bond:
    print(len(x))
print('-'*10)

beads  = list(zip(one.values(),two.values(),three.values(),four.values()))
series = list(zip([1]*50,[2]*50,[3]*50,[4]*50))
serial = []

for s,b in zip(series,beads):
    s = sorted(s,key=lambda x: b[s.index(x)])
    serial.append(s)


for s in serial:
    s.remove(1)
    s.remove(2)
print(len(serial),'-------',serial)
print('-'*10)
# print('{}-->{}-->{}'.format(*[4,3,2]),round(serial.count([4,3,2])*100/len(one),1),'%')
# print([4,2,3],round(serial.count([4,2,3])*100/len(one),1),'%')
# print(round(serial.count([3,2,4])*100/len(one),1),'%')
# print('{}-->{}-->{}'.format(*[2,4,3]),round(serial.count([2,4,3])*100/len(one),1),'%')
# print('{}-->{}-->{}'.format(*[3,4,2]),round(serial.count([3,4,2])*100/len(one),1),'%')
# print('{}-->{}-->{}'.format(*[2,3,4]),round(serial.count([2,3,4])*100/len(one),1),'%')

print(round(serial.count([4,3])*100/len(one),1),'%')
print(round(serial.count([3,4])*100/len(one),1),'%')



#plot bonds
label = ['Bond-1','Bond-2','Bond-3','Bond-4']
i=1
for bond in [two,three,four]:
    print(bond)
    plt.plot(bond.keys(), bond.values(),label=label[i])
    i+=1
plt.legend()
plt.xlabel('Bond')
plt.ylabel('Time (ps) when bond order changes')