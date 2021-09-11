import random
import numpy as np
import pypartitions
import json
from math import sqrt, pi
from projectq.ops import All, Measure, X, Y, Z, Rx, Ry, Rxx


#already compiled out all single qubit gates on ancillas, and the hadamards at beginning of x type stabs
#to do this required setting v in all x type cnots to same value (taken as +1), and no remaining gates
# in x type entangling ops have the v sign remaining, so v not a parameter for x type entangling op
# in z type stabs, v must alternate between steps in same stabilizer, so we have v as a param, but with a default
# still need to assign s to cancel some remaining Rx(-s*pi/2) gates, and see if reduced circuit produced
# (see New J. Phys. 20 (2018) 043038)
def process_syndrome_on_fly(prev_syndrome,syndrome):

    flips = (prev_syndrome - syndrome) % 2

    ft_syndrome = np.zeros(len(syndrome[0]),dtype=np.int8)

def insert_2q_error(qubit1,qubit2):
    rand = random.random()
    if rand<1/15:
        #I|qubit1
        X|qubit2
    if 1/15<rand<(2*1/15):
        # I|qubit1
        Y|qubit2
    if (2*1/15)<rand<(3*1/15):
        # I|qubit1
        Z|qubit2
    if (3*1/15)<rand<(4*1/15):
        X|qubit1
        #I|qubit2
    if (4*1/15)<rand<(5*1/15):
        X|qubit1
        X|qubit2
    if (5*1/15)<rand<(6*1/15):
        X|qubit1
        Y | qubit2
    if (6*1/15)<rand<(7*1/15):
        X|qubit1
        Z | qubit2
    if (7*1/15)<rand<(8*1/15):
        Y|qubit1
        #I|qubit2
    if (8*1/15)<rand<(9*1/15):
        Y|qubit1
        X|qubit2
    if (9*1/15)<rand<(10*1/15):
        Y|qubit1
        Y|qubit2
    if (10*1/15)<rand<(11*1/15):
        Y|qubit1
        Z|qubit2
    if (11*1/15)<rand<(12*1/15):
        Z|qubit1
        # I|qubit2
    if (12*1/15)<rand<(13*1/15):
        Z|qubit1
        X|qubit2
    if (13*1/15)<rand<(14*1/15):
        Z|qubit1
        Y|qubit2
    if (14*1/15)<rand<(15*1/15):
        Z|qubit1
        Z|qubit2

def insert_1q_error(qubit):
    rand = random.random()
    if rand < 1 / 3:
        X | qubit
    if 1 / 3 < rand < (2 * 1 / 3):  # depolarizing noise model
        Y | qubit
    if (2 * 1 / 3) < rand < (3 * 1 / 3):
        Z | qubit
def x_type_entangling(dataq,ancillaq,s=1,p2q=0,p1q=0):
    Rxx(s*pi/2) | (dataq,ancillaq)  # diff in project q def vs New J. Phys. 20 (2018) 043038 ms(chi) def by factor of 2
    if random.random()<p2q:
        insert_2q_error(dataq,ancillaq)
    Rx(-s*pi/2) | dataq
    if random.random()<p1q:
        insert_1q_error(dataq)

def z_type_entangling(dataq,ancillaq,v=1,s=1,p2q=0,p1q=0):
    Ry(v*pi/2) | dataq
    if random.random()<p1q:
        insert_1q_error(dataq)
    Rxx(s * pi / 2) | (dataq, ancillaq)# diff in project q def vs New J. Phys. 20 (2018) 043038 ms(chi) def by factor of 2

    if random.random()<p2q:
        insert_2q_error(dataq,ancillaq)
    Rx(-s * pi / 2) | dataq
    if random.random()<p1q:
        insert_1q_error(dataq)
    Ry(-v * pi / 2) | dataq
    if random.random()<p1q:
        insert_1q_error(dataq)

def z_type_entangling_no_Rx(dataq,ancillaq,v=1,s=1,p2q=0,p1q=0):
    Ry(v*pi/2) | dataq
    if random.random()<p1q:
        insert_1q_error(dataq)
    Rxx(s * pi / 2) | (dataq, ancillaq)# diff in project q def vs New J. Phys. 20 (2018) 043038 ms(chi) def by factor of 2

    if random.random()<p2q:
        insert_2q_error(dataq, ancillaq)
    Ry(-v * pi / 2) | dataq
    if random.random()<p1q:
        insert_1q_error(dataq)

def x_type_entangling_no_Rx(dataq,ancillaq, v=1, s=1, p2q=0):
    Rxx(-1*v*s*pi / 2) | (dataq, ancillaq)# angle dif by fac 2 in project q Rxx def
    # vs New J. Phys. 20 (2018) 043038 ms(chi) def
    if random.random() < p2q:
        insert_2q_error(dataq, ancillaq)
    # replace ent step with Rxx of correct sign, manually removing Rx that cancels
    # (as project q doesnt recognise this cancels as it cant comm Rx and Rxx)

def stabiliser_timestep_1(data,ancilla, old_order=False, p1q=0, p2q=0):
    if old_order:
        x_type_entangling(data[1], ancilla[1], s=1, p1q=p1q, p2q=p2q)
        x_type_entangling_no_Rx(data[3], ancilla[4], s=1, p2q=p2q)
        x_type_entangling_no_Rx(data[5], ancilla[6], s=1, p2q=p2q)

    else:
        x_type_entangling_no_Rx(data[4], ancilla[1], s=1, p2q=p2q)
        x_type_entangling(data[8], ancilla[6], s=1, p1q=p1q, p2q=p2q)
        x_type_entangling(data[6], ancilla[4], s=1, p1q=p1q, p2q=p2q)

def stabiliser_timestep_2(data,ancilla, old_order=False, p1q=0, p2q=0):
    if old_order:
        x_type_entangling_no_Rx(data[4], ancilla[1], s=-1, p2q=p2q)
        x_type_entangling(data[6], ancilla[4], s=-1, p1q=p1q, p2q=p2q)
        x_type_entangling(data[8], ancilla[6], s=-1, p1q=p1q, p2q=p2q)
    else:
        x_type_entangling(data[1], ancilla[1], s=-1, p1q=p1q, p2q=p2q)
        x_type_entangling_no_Rx(data[5], ancilla[6], s=-1, p2q=p2q)
        x_type_entangling_no_Rx(data[3], ancilla[4], s=-1, p2q=p2q)


def stabiliser_timestep_3(data,ancilla, old_order=False, p1q=0, p2q=0):
    if old_order:
        x_type_entangling(data[0], ancilla[1], s=1, p1q=p1q, p2q=p2q)
        x_type_entangling_no_Rx(data[4], ancilla[6], s=1, p2q=p2q)
        x_type_entangling(data[2], ancilla[3], s=1, p1q=p1q, p2q=p2q)
    else:
        x_type_entangling_no_Rx(data[3], ancilla[1], s=1, p2q=p2q)
        x_type_entangling(data[7], ancilla[6], s=1, p1q=p1q, p2q=p2q)
        x_type_entangling_no_Rx(data[5], ancilla[3], s=1, p2q=p2q)


def stabiliser_timestep_4(data, ancilla, old_order=False, p1q=0, p2q=0):
    if old_order:
        x_type_entangling_no_Rx(data[3],ancilla[1], s=-1, p2q=p2q)
        x_type_entangling(data[7], ancilla[6], s=-1, p1q=p1q, p2q=p2q)
        x_type_entangling_no_Rx(data[5], ancilla[3], s=-1, p2q=p2q)
    else:
        x_type_entangling(data[0], ancilla[1], s=-1, p1q=p1q, p2q=p2q)
        x_type_entangling_no_Rx(data[4], ancilla[6], s=-1, p2q=p2q)
        x_type_entangling(data[2], ancilla[3], s=-1, p1q=p1q, p2q=p2q)
def stabiliser_timestep_5(data,ancilla, old_order=False, p1q=0, p2q=0):
    if old_order:
        z_type_entangling(data[2], ancilla[2], s=1, p1q=p1q, p2q=p2q)
        #z_type_entangling(data[4], ancilla[5], s=1)
        z_type_entangling_no_Rx(data[4], ancilla[5], s=1, p1q=p1q, p2q=p2q)
        #replace z_type_ent with z_ent with Rx manually removed
        # (as project q doesnt recognise this cancels with z ent in step 4, as it cant comm Rx and Rxx)
        z_type_entangling(data[8], ancilla[7], s=1, p1q=p1q, p2q=p2q)
    else:
        z_type_entangling_no_Rx(data[1], ancilla[2], s=1, p1q=p1q, p2q=p2q)
        z_type_entangling(data[3], ancilla[5], s=1, p1q=p1q, p2q=p2q)
        z_type_entangling_no_Rx(data[7], ancilla[7], s=1, p1q=p1q, p2q=p2q)

def stabiliser_timestep_6(data,ancilla, old_order=False, p1q=0, p2q=0):
    if old_order:
        #z_type_entangling(data[1], ancilla[2], s=-1)
        z_type_entangling_no_Rx(data[1], ancilla[2], s=-1, p1q=p1q, p2q=p2q)
        #replace z_type_ent with z_ent with Rx manually removed
        # (as project q doesnt recognise this cancels with z ent in step 4, as it cant comm Rx and Rxx)
        z_type_entangling(data[3], ancilla[5], s=-1, p1q=p1q, p2q=p2q)
        #z_type_entangling(data[7], ancilla[7], s=-1)
        z_type_entangling_no_Rx(data[7], ancilla[7], s=-1, p1q=p1q, p2q=p2q)
        # replace z_type_ent with z_ent with Rx manually removed
        # (as project q doesnt recognise this cancels with z ent in step 4, as it cant comm Rx and Rxx)
    else:
        z_type_entangling(data[2], ancilla[2], s=-1, p1q=p1q, p2q=p2q)
        z_type_entangling_no_Rx(data[4], ancilla[5], s=-1, p1q=p1q, p2q=p2q)
        z_type_entangling(data[8], ancilla[7], s=-1, p1q=p1q, p2q=p2q)

def stabiliser_timestep_7(data,ancilla, old_order=False, p1q=0, p2q=0):
    if old_order:
        #z_type_entangling(data[1], ancilla[0], s=1)
        z_type_entangling_no_Rx(data[1], ancilla[0], s=1, p1q=p1q, p2q=p2q)
        # replace z_type_ent with z_ent with Rx manually removed
        # (as project q doesnt recognise this cancels with z ent in step 4, as it cant comm Rx and Rxx)
        z_type_entangling(data[5], ancilla[2], s=1, p1q=p1q, p2q=p2q)
        #z_type_entangling(data[7], ancilla[5], s=1)
        z_type_entangling_no_Rx(data[7], ancilla[5], s=1, p1q=p1q, p2q=p2q)
        # replace z_type_ent with z_ent with Rx manually removed
        # (as project q doesnt recognise this cancels with z ent in step 4, as it cant comm Rx and Rxx)
    else:
        z_type_entangling_no_Rx(data[4], ancilla[2], s=1, p1q=p1q, p2q=p2q)
        z_type_entangling(data[6], ancilla[5], s=1, p1q=p1q, p2q=p2q)
        z_type_entangling(data[0], ancilla[0], s=1, p1q=p1q, p2q=p2q)
def stabiliser_timestep_8(data,ancilla, old_order=False, p1q=0, p2q=0):
    if old_order:
        z_type_entangling(data[0], ancilla[0], s=-1, p1q=p1q, p2q=p2q)
        #z_type_entangling(data[4], ancilla[2], s=-1)
        z_type_entangling_no_Rx(data[4], ancilla[2], s=-1, p1q=p1q, p2q=p2q)
        #replace z_type_ent with z_ent with Rx manually removed
        # (as project q doesnt recognise this cancels with z ent in step 4, as it cant comm Rx and Rxx)
        z_type_entangling(data[6], ancilla[5], s=-1, p1q=p1q, p2q=p2q)
    else:
        z_type_entangling(data[5], ancilla[2], s=-1, p1q=p1q, p2q=p2q)
        z_type_entangling_no_Rx(data[7], ancilla[5], s=-1, p1q=p1q, p2q=p2q)
        z_type_entangling_no_Rx(data[1], ancilla[0], s=-1, p1q=p1q, p2q=p2q)

def stabilizer_cycle(data, ancilla, eng, reset=True, old_order=False, p1q=0, p2q=0):
    '''

    :param data:
    :param ancilla:
    :param eng:
    :param reset:
    :param old_order: stabilizer measure order prevents 'hook' errors (see fowler 1% thresh paper)
    this ordering got mixed up as we have reflected the surface code across the diagonal
    (log ops perpendicular to convention) so rewrote cycle with correct order (kept old ordering to compare performance)
    :return:
    '''
    if len(data)!=9:
        raise Exception('data qubit register does not correspond to the surface 17 QEC code')
    stabiliser_timestep_1(data, ancilla, old_order=old_order, p1q=p1q, p2q=p2q)
    stabiliser_timestep_2(data, ancilla, old_order=old_order, p1q=p1q, p2q=p2q)
    stabiliser_timestep_3(data, ancilla, old_order=old_order, p1q=p1q, p2q=p2q)
    stabiliser_timestep_4(data, ancilla, old_order=old_order, p1q=p1q, p2q=p2q)
    stabiliser_timestep_5(data, ancilla, old_order=old_order, p1q=p1q, p2q=p2q)
    stabiliser_timestep_6(data, ancilla, old_order=old_order, p1q=p1q, p2q=p2q)
    stabiliser_timestep_7(data, ancilla, old_order=old_order, p1q=p1q, p2q=p2q)
    stabiliser_timestep_8(data, ancilla, old_order=old_order, p1q=p1q, p2q=p2q)

    All(Measure) | ancilla
    eng.flush()
    syndrome_t = [int(q) for q in ancilla]
    if reset:
        for a in ancilla: #reset the ancillas to 0 at end of stab round (allow for repeat rounds)
            if int(a) == 1:
                X | a
    return syndrome_t

def get_results_log_e(rounds, runs, incorrect_count, time):
    res = {
        "runs": runs,
        "rounds": rounds,
        "incorrect_count": incorrect_count,
        "time_taken": time,
          }
    return res

def get_results_log_e_run_til_fail(rounds, time, p1q, p2q):
    res = {
        "rounds_til_fail": rounds,
        "p1q_error": p1q,
        "p2q_error": p2q,
        "time_taken": time,
          }
    return res