import time
import pandas as pd
import matplotlib.pyplot as plt
import math
import numpy as np

#load in data

blood_glucose = pd.read_csv(r"C:\Users\alexa\Downloads\New_BS.csv")

#return list of keys

def getKeys(dictionary):
    lst = []
    for i in dictionary:
        lst.append(i)
        
    return lst

#return a list of numbers seperated by anything non numerical from a string

def gettimeDataArray(str):
    sub = ""
    subs = []
    for c in str:
        if c in ['0','1','2','3','4','5','6','7','8','9']:
            sub+=c
        else:
            subs.append(int(sub))
            sub = ""
            
    if sub!="":
        subs.append(int(sub))
    return subs

print(blood_glucose.dtypes)



#seperate the data into two list and discard the date part of the time data

raw_time = blood_glucose["Device Timestamp"]
raw_levels = blood_glucose["Historic Glucose mmol/L"]

day_times = []
levels = []

#which data points to use

N = 63152
#N = 55750

start = 55750
#start = 0

length = N-start

for i in range(start,N):
    Time = gettimeDataArray(raw_time[i])
    day_times.append(Time[-1]/60+Time[-2])
    levels.append(raw_levels[i])
    
    
t = float("-inf")

#state_boundaries defines how we will break the continuos glucose levels into discrete states to use markov chains
# <4 corresponds to hypo(-glycemia), >10 corresponds to varying degrees of hyper

state_boundaries = [4,7,10,15,20]

num_bounds = len(state_boundaries)

num_states = num_bounds+1

chunk_size = 1

chunks = 24

chunk_data = []

chunk_average = [0]

transition_sums = []

transition_nums = []

#setup our 0 intialised transition counters which tracks how many times each state went to each other state at each time-chunk
#along with how many times that state went to any given state
#dividing the former by the latter will give us our probabilities of transitioning from one state to another

for i in range(chunks):
    chunk_data.append([])
    transition_num1 = []
    transition_sum1 = []
    for j in range(num_states):
        transition_num1.append(0)
        transition_sum2 = []
        for k in range(num_states):
            transition_sum2.append(0)
        transition_sum1.append(transition_sum2)
    transition_sums.append(transition_sum1)
    transition_nums.append(transition_num1)

#converts from our continuos glucose levels to the state space using the state_boundaries we defined earlier
#(n boundaries gives n+1 states)

def getState(num):
    state = num_bounds
    for i in range(num_bounds):
        if num<state_boundaries[i]:
            state = i
            break
           
    return state

#combining our glucose data into hour long averaged chunks and converting to states

consecutive_state_chunks = []
consecutive_state_chunk = []
        
for i in range(0,length):
    prev_time = t
    t = day_times[i]
    
    if prev_time>t:
        
        for j in range(chunks):
            l = len(chunk_data[j])
            if l==0:
                if len(consecutive_state_chunk)>1:
                    consecutive_state_chunks.append(consecutive_state_chunk)
                    consecutive_state_chunk = []
                else:
                    consecutive_state_chunk = []  
            else:
                state = getState(sum(chunk_data[j])/l)
                consecutive_state_chunk.append((j,state))
                
        chunk_data = []
        
        for i in range(chunks):
            chunk_data.append([])
        
    l = levels[i]
    if l*0 != 0:
        continue
    j = math.floor(t/chunk_size)
    chunk_data[j].append(l)

#populating our trasition counter
    
for chunk in consecutive_state_chunks:
    prev_state = chunk[0][1]    
    for i in range(1,len(chunk)):
        temp = chunk[i]
        Time = temp[0]  
        state = temp[1]
        
        t = (Time-1)%chunks
        
        transition_sums[t][prev_state][state] +=1
        transition_nums[t][prev_state] +=1
        prev_state = state


#setting up a matrix to store our transition probabilities
#labeled temp because it can be converted to a proper np matrix to update distributions of what state it is in
#however we dont make use of this feature

#since we are actually creating 24 different transition matrices (one for each chunk) and we're working with a small data set
#not every transtion probability will be known so we infer it by averaging the known rows from the other chunks
#failing to have any information on how that state evolves we assume it transitions to every state with equal probability

trans_matrixes_temp = []

average_replaced_lists = []
average_sums = []
average_nums = []

for i in range(num_states):
    average_sums.append([0]*num_states)  
    average_nums.append(0)
    average_replaced_lists.append([])

for i in range(chunks):
    trans_matrix = []
    for j in range(num_states):
        prob_sums = transition_sums[i][j]
        prob_num = transition_nums[i][j]

        if prob_num==0:
            probs=[1/num_states]*num_states
            average_replaced_lists[j].append(i)
        else:
            probs = []
            for k in range(num_states):
                p = prob_sums[k]/prob_num
                probs.append(p)
                average_sums[j][k]+=p
            average_nums[j]+=1
                 
        trans_matrix.append(probs)
        
    trans_matrixes_temp.append(trans_matrix)
    
for i in range(num_states):
    num = average_nums[i]
    
    if num==0:
        continue    
    
    temp_prob = []
    
    for j in average_sums[i]:
        temp_prob.append(j/num)    
    
    for j in range(chunks):
        if j in average_replaced_lists[i]:
            trans_matrixes_temp[j][i] = temp_prob

#working out the most likely states and times to lead to a hypo in the next chunk
        
top_hypo = {}

for i in range(chunks):
    for j in range(1,num_states):
        p = trans_matrixes_temp[i][j][0]
        if p in top_hypo:
            top_hypo[p].append([i,j])
        else:
            top_hypo[p] = [[i,j]]
          
keys = getKeys(top_hypo)
keys.sort(reverse = True)
for k in keys:
    print(k, end=" ")
    print(top_hypo[k])

#... a hyper in the next chunk
    
top_hyper = {}

for i in range(chunks):
    for j in range(0,num_states-2):
        p1 = trans_matrixes_temp[i][j][num_states-1]
        p2 = trans_matrixes_temp[i][j][num_states-2]
        p=p1+p2
        if p in top_hyper:
            top_hyper[p].append([i,j])
        else:
            top_hyper[p] = [[i,j]]
          
keys = getKeys(top_hyper)
keys.sort(reverse = True)
for k in keys:
    print(k, end=" ")
    print(top_hyper[k])

#plotting the pobability that state 1 (10-15 mmol/L) leads to subsequent hypo in various time chunks
        
x = []
y = []

for i in range(chunks):
    p = trans_matrixes_temp[i][1][0]
    x.append(i)
    y.append(p)

x_hyper = []
y_hyper = []

for i in range(chunks):
    p=0    
    p += trans_matrixes_temp[i][3][4]
    p += trans_matrixes_temp[i][3][5]
    x_hyper.append(i)
    y_hyper.append(p)

bar_width = .43

x_positions1 = []
x_positions2 = []
for i in range(0,len(x)):
    x_positions1.append(x[i] - bar_width / 2)  # Left bar positions
    x_positions2.append(x[i] + bar_width / 2)  # Right bar positions

# Plotting the bars
plt.figure(figsize=(10, 6))
plt.bar(x_positions1, y, width=bar_width, label='Probability of Hypo (<4 mmol/l)')
plt.bar(x_positions2, y_hyper, width=bar_width, label='Probability of Hyper (>15 mmol/l)')

# Adding labels and customizing the x-ticks
plt.xlabel('Hour of the day')
plt.ylabel('Probability')
plt.title('Markov Chain Predicted Events')
plt.xticks(x, labels=x)  # Ensure x-axis ticks are at the center of each pair
plt.legend()
plt.show()

time.sleep(100000)

trans_matrixes = []

for i in range(chunks):
    trans_matrixes.append(np.array(trans_matrixes_temp[i]))    
    

prob = np.array([1,0,0,0,0,0,0])

identity = np.array([[1,0,0,0,0,0,0],
                     [0,1,0,0,0,0,0],
                     [0,0,1,0,0,0,0],
                     [0,0,0,1,0,0,0],
                     [0,0,0,0,1,0,0],
                     [0,0,0,0,0,1,0],
                     [0,0,0,0,0,0,1]])

t = 0

print(prob)
    
"""for i in range(500):
    prob = prob@trans_matrixes[t]
    print(t, end=" ")
    print(prob)
    t=(t+1)%chunks"""
   
for t in range(chunks):
    print(trans_matrixes[t])
