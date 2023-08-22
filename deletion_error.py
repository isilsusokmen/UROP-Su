#INPUT - error rate, confusion matrix 
#OUTPUT - table with index of barker and peak value without error, error rate, and index of barker again 
from tabulate import tabulate
import numpy as np 
import matplotlib.pyplot as plt
import random
import pandas as pd

barker7 = [1, 1j, -1, 1j, -1, 1j, 1]
barker11 = [1, 1j, -1, 1j, -1, 0-1j, -1, 1j, -1, 1j, 1]
barker13 = [1, 1j, -1, 0-1j, 1, 0-1j, 1, 0-1j, 1, 0-1j, -1, 1j, 1]
elements = [1,-1,1j,-1j]
vector_length = 60

def check_barker_occurrence(vector):
    barker_patterns = [barker7, barker11, barker13]
    for pattern in barker_patterns:
        count = 0
        for i in range(len(vector) - len(pattern) + 1):
            if vector[i:i + len(pattern)] == pattern:
                count += 1
                if count > 1:
                    return False
    return True

def inject_barker(vector):
    barker_choice = random.choice([barker7, barker11, barker13])
    random_index = np.random.randint(0, len(vector) + 1)
    vector = np.insert(vector, random_index, barker_choice)
    return vector, barker_choice

def generate_vector_with_single_barker():
    while True:
        randomVector = [random.choice(elements) for _ in range(vector_length)]
        if check_barker_occurrence(randomVector):  # it returns True if no Barker codes exist
            randomVector, barker_choice = inject_barker(randomVector)
            return randomVector, barker_choice

randomVector, barker_choice = generate_vector_with_single_barker()

#finds barker without error 
findBarkerInSignal = np.correlate(randomVector, barker_choice, mode='full') 
peak_index = np.argmax(findBarkerInSignal)
peak_value = findBarkerInSignal[peak_index]
plt.title('Convolution of Random Signal with Barker')
plt.xlabel("Time")
plt.ylabel("Amplitude")
plt.plot(findBarkerInSignal)
plt.plot(peak_index, peak_value, 'ro', label='Peak')
plt.legend()
plt.show()

def deletion_error(vector,error_rate):
    
    barker_seq_w_errors=vector.copy()
    index = 0
    
    while index < len(barker_seq_w_errors):
        # Check if error will occur for the element
        if random.random() < error_rate:
            # Remove the element at the current index
            barker_seq_w_errors.pop(index)
        else:
            index += 1  # Move to the next index

    return barker_seq_w_errors

error_rate = 0.8

#find barker with errors 
barker_seq_w_errors = deletion_error(randomVector,error_rate)
findBarkerInSignal_w_error = np.correlate(barker_seq_w_errors, barker_choice, mode='full') 
peak_index2 = np.argmax(findBarkerInSignal_w_error)
peak_value2 = findBarkerInSignal_w_error[peak_index2]
plt.title('Convolution of Random Signal with Barker')
plt.xlabel("Time")
plt.ylabel("Amplitude")
plt.plot(findBarkerInSignal_w_error)
plt.plot(peak_index2, peak_value2, 'ro', label='Peak')
plt.legend()
plt.show()

#comparing barker detection in sequences with/without error 
table = [["Barker index", peak_index, peak_index2],
         ["Barker peak value", f"{peak_value.real}+{peak_value.imag}j", f"{peak_value2.real}+{peak_value2.imag}j"],
         ["Error Rate", 0, error_rate]]
print(tabulate(table, headers=["", "Without Error", "With Deletion Error"]))