#   What is included:
#       X-remove taper of a Hamiltonian.
#       Construction of ansatz for X-removed Hamiltonian.
import numpy as np
from qiskit import QuantumCircuit
from qiskit.opflow import I, X, Y, Z
from qiskit.circuit.library import XXPlusYYGate
from qiskit.opflow import PauliSumOp
from qiskit.circuit import ParameterVector
from itertools import combinations

## find the ZZ_symmetry
## input : hamiltonian
## ouput : ZZ generator
def find_ZZ_sym(hamiltonian):
    generator = []
    ham_str = []
    for ham in hamiltonian:
            str_ham = str(ham.primitive.paulis[0])
            str_ham = list(str_ham)
            ham_str.append(str_ham)
    for Z in ham_str:
        if Z == ['I']*hamiltonian.num_qubits:
            continue
        count = 0
        Bool = 1
        for string in ham_str:
            for i,word in enumerate(string):
                if (word != Z[i] and word != 'I' and Z[i] != 'I'):
                    count+=1
            if count%2 != 0:
                Bool = 0
        if Bool:
            g_word =''
            for word in Z:
                g_word = g_word + word
            generator.append(g_word)
    pairity = ['Z'*int(hamiltonian.num_qubits/2)+'I'*int((hamiltonian.num_qubits+1)/2),'I'*int((hamiltonian.num_qubits+1)/2)+'Z'*int(hamiltonian.num_qubits/2),'Z'*hamiltonian.num_qubits]
    generator_pair = []
    for Z in pairity: 
        count = 0
        Bool = 1
        for string in ham_str:
            for i,word in enumerate(string):
                if (word != Z[-i-1] and word != 'I' and Z[-i-1] != 'I'):
                    count+=1
            if count%2 != 0:
                Bool = 0
        if Bool:
            g_word =''
            for word in Z:
                g_word = g_word + word
            generator_pair.append(g_word)
    def dot_product(set1,set2):
        res = ''
        for index in range(len(set1)):
            if(set1[index] == 'Z' and set2[index] == 'I') or (set2[index] == 'Z' and set1[index] == 'I'):
                res = res + 'Z'
            else:
                res = res + 'I'
        return res

    def find_non_triv(full_set):
        for part in full_set:
            for set in full_set:
                res = dot_product(part,set)
                try:
                    full_set.remove(res)
                except:
                    continue
        return full_set
    
    generator = find_non_triv(generator)
    generator_pair = find_non_triv(generator_pair)
    pair_term = []
    for index, j in enumerate(generator_pair):
        if index == 1:
            pair_term.append(j)
            continue
        temp = j
        for i in generator:
            temp = dot_product(temp,i)
        pair_term.append(temp)
    generator = generator+pair_term
    generator = find_non_triv(generator)
    return generator

## Find the sector which gonna be taper off
## input : generator, pairity remove select
## ouput : taper off list
def Z_point_select(generator,slect = 0):
    result = []
    for Z in generator:
        temp = []
        index = 0
        index_word = 0
        bool = False
        for i in range(len(Z)):
            if Z[-i-1] == 'Z':
                index_word = index
                bool = True
                for Z_in in generator:
                    if Z != Z_in:
                        if Z_in[-index-1] == 'Z':
                            bool = False 
            if bool:
                temp.append(index_word)
            index += 1
        result.append(temp)
    res = []
    for i in result:
        res.append(i[slect])
    return res

## Find the sign of taper offed state eigen value
## input : generator, electron number
## ouput : signs of tapered state list
def joint_mat(generator,num_elec = (4,5)):
    result = []
    for i in generator:
        part1 = i[0:int(len(generator[0])/2)]
        part2 = i[int(len(generator[0])/2):]
        res = 1
        for i in range(len(part1)):
            if part1[-1-i] == 'Z' and i+1<=num_elec[1]:
                res*=-1
        for i in range(len(part2)):
            if part2[-1-i] == 'Z' and i+1<=num_elec[0]:
                res*=-1    
        result.append(res)
    return result

## Assign the sign in tapered hamiltonian
def assign_sign(joint_mat,Z_index,taper_ham):
    sign_list = []
    index = 0
    for i,sign in enumerate(joint_mat):
        if (sign == -1 and index == 0):
            for ham in taper_ham:
                if str(ham.primitive.paulis[0])[-Z_index[i]-1] == 'X':
                    sign_list.append(-1)
                else:
                    sign_list.append(1)
        if (sign == -1 and index != 0):
            for j,ham in enumerate(taper_ham):
                if str(ham.primitive.paulis[0])[-Z_index[i]-1] == 'X':
                    sign_list[j] *= -1
        if len(sign_list)!= 0:
            index +=1
    
    sum_list = 0
    for i in sign_list:
        sum_list += i
    if sum_list == len(sign_list):
        for i in range(len(taper_ham)):
            sign_list.append(1)
    result = (I^taper_ham.num_qubits)*0

    for i,ham in enumerate(taper_ham):
        result += sign_list[i]*ham

    return result.reduce()

## tapering method
def taper(hamiltonian,generator,num_elec,joint_mat1 = [],select = 0):
    ## taper 위치 계산
    Z_index = Z_point_select(generator,select)
    ## t + Z unitary 계산 및 UHU 계산
    Unitary = (I^hamiltonian[0].num_qubits)*1
    for index,gen in enumerate(generator):
        temp = list('I'*hamiltonian[0].num_qubits)
        temp[Z_index[index]] = 'X'
        temp.reverse()
        Op = ''
        for word1 in temp:
            Op = Op+word1
        Unitary = Unitary@PauliSumOp.from_list([(gen,np.sqrt(1/2)),(Op,np.sqrt(1/2))])
    taperop1 = (Unitary.adjoint()@hamiltonian@Unitary)
    taperop1 = taperop1.reduce()
    if len(joint_mat1) == 0:
        taperop1 = assign_sign(joint_mat=joint_mat(generator,num_elec),Z_index=Z_index,taper_ham= taperop1)
    else:
        taperop1 = assign_sign(joint_mat=joint_mat1,Z_index=Z_index,taper_ham= taperop1)
    taperop1 = taperop1.reduce()
    taper1 = (I^(hamiltonian[0].num_qubits-len(generator)))*0
    for op in taperop1:
        str_op = list(str(op.primitive.paulis[0]))
        for Z in Z_index:
            str_op[-Z-1] = ''
        temp = ''
        for word in str_op:
            temp = temp + word
        taper1+=PauliSumOp.from_list([(temp,op.coeffs[0])])
    return taper1.reduce()


## This sector is reconstruction of hamiltonian remove X and taper again
## remove the sigle X gate term in hamiltonian
def Single_X_remove(hamiltonian):
    X_taperop = PauliSumOp.from_list([("I"*hamiltonian[0].num_qubits,0)])
    for op in hamiltonian:
        temp = op.primitive.paulis.x[0]
        count = 0
        for i in temp:
            if i:
                count+=1
        if count%2 == 0:
            X_taperop += op
    X_taperop = X_taperop.reduce()
    return X_taperop

def Single_X_get(hamiltonian):
    X_taperop = PauliSumOp.from_list([("I"*hamiltonian[0].num_qubits,0)])
    for op in hamiltonian:
        temp = op.primitive.paulis.x[0]
        count = 0
        for i in temp:
            if i:
                count+=1
        if count%2 == 1:
            X_taperop += op
    X_taperop = X_taperop.reduce()
    return X_taperop
## find the single exitation terms in hamiltonian
def find_single_up(hamiltonian,pauli,hatree = []):
    single_list_y = []
    for op in hamiltonian:
        if str(op.primitive.paulis).count(pauli) == 2:
            temp = []
            tp = []
            for i,string in enumerate(op.primitive.paulis[0]):
                if str(string) == pauli:
                    temp.append(i)
            tp.append(temp[0]*10+temp[1])
            res = 1
            if len(hatree)!=0:
                for i in range(temp[0],temp[1]):
                    if hatree[i] == 1 and str(op.primitive.paulis[0])[i]=='Z':
                        res*=1
            tp.append(op.coeffs[0])
            single_list_y.append(tp)  
    single_list_y.sort()
    np_list_y = []
    for i in single_list_y:
        if len(np_list_y) == 0:
            np_list_y.append(i)
        else:
            if np_list_y[-1][0] == i[0]:
                np_list_y[-1][1] += i[1]
            else:
                np_list_y.append(i)
    return np_list_y
## remove the only XX term exitation
def remove_single_up(hamiltonian,pauli,remove_list):
    result = (I^hamiltonian.num_qubits)*0
    for op in hamiltonian:
        if str(op.primitive.paulis).count(pauli) == 2:
            temp = []
            for i,string in enumerate(op.primitive.paulis[0]):
                if str(string) == pauli:
                    temp.append(i)
            if remove_list.count(temp[0]*10+temp[1])>0.5:
                continue
        result += op
    return result
def find_only_XX(hamiltonian):
    X_list = find_single_up(hamiltonian,'X')
    Y_list = find_single_up(hamiltonian,'Y')
    remove_list = []
    for index in X_list:
        bool = False
        for indey in Y_list:
            if index[0] == indey[0]:
                bool = True
        if not(bool):
            remove_list.append(index[0])
    result = remove_single_up(hamiltonian,'X',remove_list)
    return result.reduce()

## taper with X remove
## reps is how many do X remove taper if reps is one, taper the qubits and remove proper X terms if reps > 1, taper more qubits
## output : tapered hamiltonian and Fock state list
def taper_X_remove(hamiltonian,reps,electron):
    generator = find_ZZ_sym(hamiltonian)
    taper_ham = taper(hamiltonian,generator,electron)
    tapped_term = Z_point_select(generator)
    all_term = [i for i in range(hamiltonian[0].num_qubits)]
    part_term = []
    
    for i in all_term:
        if i not in tapped_term:
            part_term.append(i)
    temp = taper_ham
    part_down = []
    part_up = []

    down_num = 0
    up_num = 0

    for i in part_term:
        if i < hamiltonian[0].num_qubits/2:
            part_down.append(i)
            if i < electron[0]:
                down_num+=1
        else:
            part_up.append(i)
            if i < electron[1]+hamiltonian[0].num_qubits/2:
                up_num += 1
    for i in range(reps-1):
        if part_down[0] <= part_up[0]-hamiltonian[0].num_qubits/2:
                if (up_num+down_num) %2 == 0:
                    sign = 1
                else:
                    sign = -1
                temp = Single_X_remove(temp)
                temp = find_only_XX(temp)
                temp_generator = ['Z'*temp.num_qubits]
                tapped_term = Z_point_select(temp_generator)
                elec = tuple([down_num,up_num])
                temp = taper(temp,temp_generator,electron,joint_mat1=[sign])
                if part_down[0] < electron[0]:
                    down_num -= 1
                del part_down[0]
        else:
                if (up_num+down_num) %2 == 0:
                    sign = 1
                else:
                    sign = -1
                temp = Single_X_remove(temp)
                temp = find_only_XX(temp)
                temp_generator = ['Z'*temp.num_qubits]
                tapped_term = Z_point_select(temp_generator)
                elec = tuple([down_num,up_num])
                temp = taper(temp,temp_generator,electron,select=len(part_down),joint_mat1=[sign])
                if part_up[0] < electron[1]+hamiltonian[0].num_qubits/2:
                    up_num -= 1
                del part_up[0]
    temp = Single_X_remove(temp)
    temp = find_only_XX(temp)
    temp = find_only_XX(temp)
    full_part = part_down+part_up
    Hatree = []
    for i in full_part:
        if (i < hamiltonian[0].num_qubits/2 and i < electron[0]) or (i >= hamiltonian[0].num_qubits/2 and i < electron[1]+hamiltonian[0].num_qubits/2):
            Hatree.append(1)
        else:
            Hatree.append(0)
    return temp,Hatree

## construct ansatz at X_tapered hamilitonian
def Givens_rotation(para):
    qc = QuantumCircuit(2)
    qc.h(0)
    qc.cx(0,1)
    qc.ry(para/2,0)
    qc.ry(para/2,1)
    qc.cx(0,1)
    qc.h(0)
    return qc.to_gate(label=f"Givens")

def EF_ansatz(hamiltonian,Hatree,Fraction=1/10,reverse = False):
    single_up_list = find_single_up(hamiltonian,'Y',Hatree)
    max = 0
    for part in single_up_list:
        if max<abs(part[1]):
            max = abs(part[1])
    select = []
    sign = []
    for part in single_up_list:
        if abs(part[1])>max*Fraction:
            select.append(part[0])
            sign.append(part[1])
    paras = ParameterVector('T',len(select)*2)
    AF_paras = ParameterVector('U',len(select)*2)
    qc = QuantumCircuit(hamiltonian.num_qubits)
    for i, temp in enumerate(Hatree):
        if temp == 1:
            qc.x(i)
    late_sector = []
    index = 0
    if not(reverse): 
        for i,sector in enumerate(select):
            if Hatree[int(sector/10)] == 1 and Hatree[int(sector%10)] == 1:
                late_sector.append(sector)
                continue
            if sign[i]<0:
                qc.append(Givens_rotation(paras[index]),[int(sector/10),sector%10])
            else:
                qc.append(Givens_rotation(paras[index]),[sector%10,int(sector/10)])
            index += 1 
        index = 0 
        for i,sector in enumerate(late_sector):
            qc.append(Givens_rotation(AF_paras[index]),[int(sector/10),sector%10])
            index += 1
    else:
        for i,sector in enumerate(select):
            if Hatree[int(sector/10)] == 1 and Hatree[int(sector%10)] == 1:
                late_sector.append(sector)
                continue
            if sign[i]>0:
                qc.append(Givens_rotation(paras[index]),[int(sector/10),sector%10])
            else:
                qc.append(Givens_rotation(paras[index]),[sector%10,int(sector/10)])
            index += 1 
        index = 0 
        for i,sector in enumerate(late_sector):
            qc.append(Givens_rotation(AF_paras[index]),[int(sector/10),sector%10])
            index += 1
    return qc

## Grouping Sector

import networkx as nx
from collections import defaultdict
## check if two Pauli terms are single pauli commute
def is_single_commute(Ham1, Ham2):
    str1 = str(Ham1.primitive.paulis[0])
    str2 = str(Ham2.primitive.paulis[0])
    for i in range(len(str1)):
        if str1[i] != 'I' and str2[i] != 'I':
            if str1[i] != str2[i]:
                return False
    return True


## Make the graph which is single commute
def Graph_made(Ham):
    num=len(Ham)
    Mat = []
    for i in range(num):
        line = []
        for j in range(num):
            if i==j:
                input = 0
            else:
                if is_single_commute(Ham[i],Ham[j]):
                    input = 1
                else:
                    input = 0
            line.append(input)
        Mat.append(line)
    return Mat

## grouping the hamiltonian
def group_ham(hamiltonian):
    result = []
    Mat = Graph_made(hamiltonian)
    G = nx.Graph()
    for i in range(len(Mat)):
        G.add_node(i)

    for i in range(len(Mat)):
        for j in range(i + 1, len(Mat)):
            if Mat[i][j] == 1:
                G.add_edge(i, j)
    C = nx.complement(G)
    groups = nx.coloring.greedy_color(C, strategy="largest_first")
    num_groups = len(set(groups.values()))
    for i in range(num_groups):
        temp = []
        for term, group_id in groups.items():
            if group_id == i:
                temp.append(hamiltonian[term])
        result.append(temp)
    return result

## extract the file's text as hamiltonian

def extract(text,qubits_num =12):
    global P_string
    temp = 0*(I^qubits_num)
    lines = text.split('\n')
    for line in lines:
        P_string = 0*I
        value = ""
        index = 0
        I_index = 0
        for word in line:
            if word != "*" and index==0:
                value+=word
            else:
                index +=1
            if index >=1:
                if I_index == 0:
                    if word == "I":
                        P_string += I
                        I_index+=1
                        
                    elif word == "X":
                        P_string += X
                        I_index+=1
                        
                    elif word == "Y":
                        P_string += Y
                        I_index+=1
                        
                    elif word == "Z":
                        P_string += Z
                        I_index+=1
                    else:
                        I_index+=0
                else:
                    if word == "I":
                        P_string = P_string^I
                        I_index+=1
                    elif word == "X":
                        P_string = P_string^X
                        I_index+=1
                    elif word == "Y":
                        P_string = P_string^Y
                        I_index+=1
                    elif word == "Z":
                        P_string = P_string^Z
        try:
            value = value.replace(" ", "")
            value = float(value)        
            temp = (temp) + ((P_string)*(value))
        except:
            continue
    return temp.reduce()
def Opeartation_string(PauliOPS):
    result = 'I'*PauliOPS[0].num_qubits
    result = list(result)
    for OP in PauliOPS:
        temp = str(OP.primitive.paulis[0])
        index = 0
        for word in temp:
            if word !='I' and word !='Z':
                result[index] = word
            index+=1
    return result
def operation_circ(ansatz,g_ham,para=None):
    result_list = []
    ansatz_temp = ansatz.copy()
    temp_index = 0
    for i in ansatz_temp.parameters:
        ansatz_temp = ansatz_temp.bind_parameters({i:para[temp_index]})
        temp_index+=1
    for part_ham in g_ham:
        temp = ansatz_temp.copy()
        temp.barrier(temp.qubits)
        op_str = Opeartation_string(part_ham)
        for i in range(len(op_str)):
            if op_str[-(i+1)]=='X':
                temp.h(i)
            elif op_str[-(i+1)]=='Y':
                temp.sdg(i)
                temp.h(i)
        if temp.num_clbits == 0:
            temp.measure_all()
        result_list.append(temp)
    return result_list

def expectation_value_circ(para,ansatz,g_ham,sampler,cb = None):
    circ_list = operation_circ(ansatz,g_ham,para)
    job = sampler.run(circ_list)
    result_list = job.result()
    index = 0
    exp_value = 0
    g_str = []
    for i in range(len(g_ham)):
        temp_list=[]
        for part_of_group in g_ham[i]:
            P_str = str(part_of_group.primitive.paulis[0])
            coeffs = part_of_group.coeffs
            temp_list.append([P_str,coeffs])
        g_str.append(temp_list)
    for index in range(len(g_ham)):
        for sampler_result in result_list.quasi_dists[index].keys():
            prob = result_list.quasi_dists[index][sampler_result]
            result_str = bin(sampler_result)
            result_str = result_str[2:]
            for part_of_group in g_str[index]:
                res = 1
                coeffs = part_of_group[1]
                P_str = part_of_group[0]
                for P in range(len(result_str)):
                    if P_str[-(P+1)]!='I' and result_str[-(P+1)]=='1':
                        res*=-1
                exp_value+=coeffs*res*prob
    try:
        cb(exp_value[0].real,para)
    except:
        None
    return exp_value[0].real


def hatree_value(hatree,hamiltonian):
    #Bloch sphere cos(t/2)+e^(ip)sin(t/2)
    hatree.reverse()
    g_str = []
    for P in hamiltonian:
        P_str = str(P.primitive.paulis[0])
        coeffs = P.coeffs
        g_str.append([P_str,coeffs])

    exp_value = 0
    for part_of_group in g_str:
        res = 1
        coeffs = part_of_group[1]
        P_str = part_of_group[0]
        temp = 0
        for P in P_str:
            if P == "I":
                res *= 1
            elif P == "X":
                res *= (np.cos((hatree[temp]*np.pi-np.pi/2)/2))**2-(np.sin((hatree[temp]*np.pi-np.pi/2)/2))**2
            elif P == "Y":
                res *= 0
            elif P == "Z":
                res *= (np.cos(hatree[temp]*np.pi/2))**2-(np.sin(hatree[temp]*np.pi/2))**2
            temp+=1
        exp_value+=coeffs*res
    hatree.reverse()
    return exp_value

def Z_exchange(hamiltonian):
    Z_inverse = (I^hamiltonian[0].num_qubits)*0
    P_str = []
    for op in hamiltonian:
        P_str.append(str(op.primitive.paulis[0]))

    for index,op in enumerate(P_str):
        res = 1
        for word in op:
            if word == 'Z':
                res*=-1
        Z_inverse += PauliSumOp.from_list([(op,hamiltonian[index].coeffs[0]*res)])
    return Z_inverse.reduce()




'''
readout error mitigation implementation on expectation value circ function.
'''

from qiskit.result.distributions.quasi import QuasiDistribution



def extract_ro_noise(noise_model):
    ro_noises = {}
    for noise in noise_model.to_dict()["errors"]:
        if noise["type"] == "roerror":
            ro_noises[noise["gate_qubits"][0][0]] = noise["probabilities"]
    return ro_noises

def generate_corr_mat(qubits:list, ro_noises:dict):
    err_mat = ro_noises[qubits[0]]
    for i in range(1,len(qubits)):
        err_mat = np.kron(err_mat, ro_noises[qubits[i]])
    corr_mat = np.linalg.inv(err_mat)
    return corr_mat

def ro_err_mit(quasi_dists, qubits:list, ro_noises:dict):
    dists_list = []
    for i in range(2**len(qubits)):
        if i in quasi_dists.keys():
            dists_list.append(quasi_dists[i])
        else:
            dists_list.append(0)
    corr_mat = generate_corr_mat(qubits, ro_noises)
    dists_list = np.dot(dists_list, corr_mat)
    corr_dists = {}
    for i in range(2**len(qubits)):
        if dists_list[i] != 0:
            corr_dists[i] = dists_list[i]
    return QuasiDistribution(corr_dists, quasi_dists.shots)

def expectation_value_circ_err_mit(para,ansatz,g_ham,sampler, qubits:list=None, noise_model=None, err_mit = True, nearest_prob_dists = False, cb = None):
    if qubits == None:
        qubits = list(range(len(ansatz.num_qubits)))
    ro_noises = extract_ro_noise(noise_model)
    circ_list = operation_circ(ansatz,g_ham,para)
    sampler.set_options(transpile_options=dict(initial_layout=qubits))
    job = sampler.run(circ_list)
    result_list = job.result()
    index = 0
    exp_value = 0
    g_str = []
    for i in range(len(g_ham)):
        temp_list=[]
        for part_of_group in g_ham[i]:
            P_str = str(part_of_group.primitive.paulis[0])
            coeffs = part_of_group.coeffs
            temp_list.append([P_str,coeffs])
        g_str.append(temp_list)
    for index in range(len(g_ham)):
        quasi_dists = result_list.quasi_dists[index]
        if err_mit:
            quasi_dists = ro_err_mit(quasi_dists, qubits, ro_noises)
        if nearest_prob_dists:
            quasi_dists = quasi_dists.nearest_probability_distribution()
        for sampler_result in quasi_dists.keys():
            prob = quasi_dists[sampler_result]
            result_str = bin(sampler_result)
            result_str = result_str[2:]
            for part_of_group in g_str[index]:
                res = 1
                coeffs = part_of_group[1]
                P_str = part_of_group[0]
                for P in range(len(result_str)):
                    if P_str[-(P+1)]!='I' and result_str[-(P+1)]=='1':
                        res*=-1
                exp_value+=coeffs*res*prob
    try:
        cb(exp_value[0].real,para)
    except:
        None
    return exp_value[0].real