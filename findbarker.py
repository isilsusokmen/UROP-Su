import numpy as np 
import matplotlib.pyplot as plt

#barker sequences that will be implemented later 
barker7 = np.array([1,1,1,-1,-1,1,-1])
barker11 = np.array([+1, +1, +1, -1, -1, -1, +1, -1, -1, +1, -1])
barker13 = np.array([+1, +1, +1, +1, +1, -1, -1, +1, +1, -1, +1, -1, +1])

b7_str = ''.join(str(barker7))
b11_str = ''.join(str(barker11))
b13_str = ''.join(str(barker13))

#we are currently working with vectors, randomise a vector
found_barker = True 

#at this stage checking if barker sequence exists in random signal by converting to string - this can also be done w convolution 
while found_barker:
    randomVector = np.random.randint(0, 2, 13)
    randomVector[np.where(randomVector == 0)] = -1
    vct_str = ''.join(str(randomVector))
    
    if b7_str and b11_str and b13_str not in vct_str:
        found_barker = False
             
#Generate a random index and insert barker code at this location
random_index = np.random.randint(0, len(randomVector) + 1)
result_b7 = np.insert(randomVector, random_index, barker7)
result_b11 = np.insert(randomVector, random_index, barker11)
result_b13 = np.insert(randomVector, random_index, barker13)
    
#producing correlation graphs to find peak index and value to locate barker code 
findBarkerInSignal7 = np.correlate(result_b7, barker7, mode='full') 
peak_index7 = np.argmax(findBarkerInSignal7)
peak_value7 = findBarkerInSignal7[peak_index7]
plt.title('Convolution of Random Signal with Barker7')
plt.xlabel("Time")
plt.ylabel("Amplitude")
plt.plot(findBarkerInSignal7)
plt.plot(peak_index7, peak_value7, 'ro', label='Peak')
plt.legend()
plt.show()

print("Peak index:", peak_index7)
print("Peak value:", peak_value7)
plt.show()

findBarkerInSignal11 = np.correlate(result_b11, barker11, mode='full') 
peak_index11 = np.argmax(findBarkerInSignal11)
peak_value11= findBarkerInSignal11[peak_index11]
plt.title('Convolution of Random Signal with Barker11')
plt.xlabel("Time")
plt.ylabel("Amplitude")
plt.plot(findBarkerInSignal11)
plt.plot(peak_index11, peak_value11, 'ro', label='Peak')
plt.legend()
plt.show()

print("Peak index:", peak_index11)
print("Peak value:", peak_value11)
plt.show()

findBarkerInSignal13 = np.correlate(result_b13, barker13, mode='full') 
peak_index13 = np.argmax(findBarkerInSignal13)
peak_value13= findBarkerInSignal13[peak_index13]
plt.title('Convolution of Random Signal with Barker13')
plt.xlabel("Time")
plt.ylabel("Amplitude")
plt.plot(findBarkerInSignal13)
plt.plot(peak_index13, peak_value13, 'ro', label='Peak')
plt.legend()
plt.show()

print("Peak index:", peak_index13)
print("Peak value:", peak_value13)
plt.show()

