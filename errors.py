#Introduce insertions, substitutions and deletions and see what happens to the autocorrelation of barker with itself.

import numpy as np
import random 
import matplotlib.pyplot as plt

#function assumes input vector contains 1 and -1s and introduces errors - insertion or deletion or substitution - to the signal
def introduce_errors(signal, errorType):

    length = len(signal)

    if errorType == 'deletion':
        if length <= 1:
            return signal
        index = random.randint(0, length - 1)
        signal_with_error = np.delete(signal,index)

    elif errorType == 'substitution':
        if length == 0:
            return signal
        index = random.randint(0, length - 1)
        if signal[index]==1:
            new_element = -1
        elif signal[index]==-1:
            new_element = 1 
        signal_with_error = signal.copy()
        signal_with_error[index] = new_element

    elif errorType == 'insertion':
        index = random.randint(0, length)
        new_element = random.choice([-1, 1])
        signal_with_error = np.insert(signal,index,new_element)

    else:
        raise ValueError("Invalid error type specified. Choose from 'deletion', 'substitution', or 'insertion'.")

    return signal_with_error

signal = np.random.randint(0, 2, 13)
signal[np.where(signal == 0)] = -1

signal_w_errors = introduce_errors(signal, 'substitution')
print("original signal:", signal, " length of original signal:", len(signal))
print("signal w error:", signal_w_errors, " length of signal w errors", len(signal_w_errors))
