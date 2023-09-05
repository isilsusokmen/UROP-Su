
#FUNCTION
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
vector_length = 200

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
barker_choice = np.array(barker_choice)

#finds barker without error 
findBarkerInSignal = np.correlate(randomVector, barker_choice, mode='full') 
peak_index = np.argmax(findBarkerInSignal)
peak_value = findBarkerInSignal[peak_index]
plt.title('Convolution of Random Signal with Barker')
plt.xlabel("Time")
plt.ylabel("Amplitude")
plt.plot(np.real(findBarkerInSignal), label='Real Part')  # Plot the real part of the complex array
plt.plot(np.imag(findBarkerInSignal), label='Imaginary Part')  # Plot the imaginary part of the complex array
plt.plot(peak_index, np.real(peak_value), 'ro', label='Peak')
plt.legend()
plt.show()


#function assumes probability matrix is normalized 
def substitution_error(vector,error_rate, conf_matrix):
    
    barker_seq_w_errors=vector.copy()
    #mapping system 
    nucelotide_to_barker = {'A': 1, 'C': 1j, 'G':0-1j, 'T':-1}
    barker_to_nucleotide = {1: 'A', 1j: 'C', 0-1j: 'G', -1: 'T'}

    for i in range(len(vector)):
        
        #check if error will occur for the nucelotide 
        if np.random.rand() < error_rate:
            barker_value = vector[i]
        
            #use confusion matrix for nucelotide substitution probabibilty 
            probabilities = confusion_matrix[list(barker_to_nucleotide.keys()).index(barker_value)]
            new_value = np.random.choice(list(barker_to_nucleotide.keys()), p=probabilities)
            
            barker_seq_w_errors[i] = new_value
    
    return barker_seq_w_errors

#the confusion matrix found in the course (A,C,G,T)- order mantained 
confusion_matrix = [
    [0.80, 0.0, 0.0, 0.2],
 [0.05, 0.10, 0.60, 0.25],
 [0.10, 0.25, 0.60, 0.05],
 [0.2, 0.0, 0.0, 0.80]
]
error_rate = 0.2


#find barker with errors 
barker_seq_w_errors = substitution_error(randomVector,error_rate,confusion_matrix)
barker_seq_w_errors = np.array(barker_seq_w_errors)
findBarkerInSignal_w_error = np.correlate(barker_seq_w_errors, barker_choice, mode='full') 
peak_index2 = np.argmax(findBarkerInSignal_w_error)
peak_value2 = findBarkerInSignal_w_error[peak_index2]
plt.title('Convolution of Signal w Susbstitution error with Barker')
plt.xlabel("Time")
plt.ylabel("Amplitude")
plt.plot(np.real(findBarkerInSignal_w_error), label='Real Part')  # Plot the real part of the complex array
plt.plot(np.imag(findBarkerInSignal_w_error), label='Imaginary Part')  # Plot the imaginary part of the complex array
plt.plot(peak_index2, np.real(peak_value2), 'ro', label='Peak')
plt.legend()
plt.show()

table = [["Barker index", peak_index, peak_index2],
         ["Barker peak value", f"{peak_value.real}+{peak_value.imag}j", f"{peak_value2.real}+{peak_value2.imag}j"],
         ["Error Rate", 0, error_rate]]
print(tabulate(table, headers=["", "Without Error", "With Substitution Error"]))
