import itertools
import numpy as np
import json

lists=[[0,1] for _ in range(18)]
print(lists)

print(len(list(itertools.product(*lists))))

synd_list = [[0,1] for _ in range(8)]
correction_table = {}

for ele in list(itertools.product(*synd_list)): #create empty lookup table for all possible syndromes
    key=str(np.array(ele)).strip('[,]')
    #print(key)
    correction_table[key]=['1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1', 18]
T = np.matrix([    #acts on error vector, 18 ele colum labelled x0->x8->z0->z8 to return
                  # 8 element syndrom row vector [z0,x1,z2,x3,x4,z5,x6,z7] as in ascii diagram in surface code file

    [1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,1,1,0,1,1,0,0,0,0],
    [0,1,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0],
    [0,0,0,1,1,0,1,1,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,1],
    [0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0],
            ])
def calc_weight(error_vec):
    '''
    calculate number of errors a string represents in the
      coin the depolarising noise model count y errors as same weight as a single x or z
      (don't count an x AND  z error)'''
    weight=sum(error_vec)
    for i in range(9):
        if vec[i] == 1 and vec[i + 9] == 1:
            weight -= 1

    return weight
print(len(list(itertools.product(*lists))))

# vec = [1] + 17*[0] # quick check T was giving correct weight corrections on q1 (see comment below)
# avec = np.array(vec)
# weight = calc_weight(avec)
# print(weight)
# print(np.matmul(T,avec) %2 )

for element in list(itertools.product(*lists))[:]:
    vec=np.array(element)
    #print(vec)
    weight=calc_weight(vec)
    #print('weight = {}'.format(weight))
    syndrome=np.matmul(T,vec) %2
    label=str(syndrome)
    label = label.strip('[,]')
    #print(label)
    #print(correction_table[label])
    if correction_table[label][1] > weight: #start with a max weight default
        correction_table[label] = ([int(v) for v in vec], int(weight))

  # combo of two bits below gave weight 2 corrections for 10000000, not sure why
    # for tab_ele in correction_table[label]:
    #     if tab_ele[1] > weight:
    #         del tab_ele
        # if len(tab_ele) == 0:
    # if len(correction_table[label]) == 0:
    #     correction_table[label].append(([int(v) for v in vec], int(weight)))
        # elif correction_table[label][0][1] == weight:
        #     correction_table[label].append((vec, weight))
#
#
for key in correction_table.keys():
    print((key,correction_table[key]))
#######

with open('correction_table_depolarising.json','w') as file:
    json.dump(correction_table,file)

with open('correction_table_depolarising.json','r') as infile:
    table=json.load(infile)
#
# for key in table.keys():
#     print((key,table[key]))
