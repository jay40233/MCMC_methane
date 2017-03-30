#The code is created under python 3.5.2
import sys
import getopt
import os
import math
import operator
import random
#import matplotlib.pyplot as plt
from timeit import default_timer
#define constant parameters
K_r = 367
K_theta = 35
r0 = 1.08
theta0 = 109.5
PI = math.pi

#read the methane PDB
def readPDB(fileName):
  coord = []
  f = open(fileName)
  for line in f:
    temp = line.split()
    if temp[2] != 'C':
      coord.append(float(temp[5]))
      coord.append(float(temp[6]))
      coord.append(float(temp[7]))
  f.close()
  return coord

#determine if the energy calculated is accepted or not
def accept_or_not(thresold):
  x = random.uniform(0,1)
  if x<=thresold:
    return 1
  return 0

#calculate the bond length
def bond_length(x,y,z):
  return math.sqrt(x**2+y**2+z**2)

#determine the direction for permutation
def p_or_n():
  return random.choice([1,-1])

#choose next parameter to permutate, on cartesian
def pick_next_permutation():
  return random.randrange(12)

#randomly decide the stepsize
def stepsize(step):
  return random.uniform(0,1)*step

#calculate the energy
def give_me_sum(para):
  sum = 0;
  for j in range(0,10):
    if j<4:
      sum += 1.0/2*K_r*(para[j]-r0)*(para[j]-r0)
    else:
      sum += 1.0/2*K_theta*(para[j]-theta0)*(para[j]-theta0)
  return sum

#convert from cartesian to internal coordination
def car2in(para, coord):
  d_map = {4:[0,1], 5:[0,2], 6:[0,3], 7:[1,2], 8:[1,3], 9:[2,3]}
  for i in range(0,10):
    if i<=3:
      para[i] = bond_length(coord[3*i],coord[3*i+1],coord[3*i+2])
    else:
      j, k = d_map[i][0], d_map[i][1]
      dot = coord[3*j]*coord[3*k]+coord[3*j+1]*coord[3*k+1]+coord[3*j+2]*coord[3*k+2]
      para[i] = math.acos(dot/(para[j]*para[k]))*180.0/PI
  return para

def main():
  (options, args) = getopt.getopt(sys.argv[1:], 'fbm')
  if not args:
    print("Usage: python MCMC.py [PDB file]")
    return 0
  #define parameters
  #-------------------------------------------------------
  E_tol = 0.1              	#The energy threshold, if the minimum energy drops below this point the searching will be terminated
  step = 0.05               #The max range a single step will reach
  ite_max = 100000          #Max allowable number of iterations
  #-------------------------------------------------------
  trail = {}
  while 1:
    #kT value that can be adjusted
    kT = float(input("Please enter kT value in kcal/mol (default = 0.6 kcal/mol), enter -1 to terminate:"))
    if kT==-1:
      break
    start_time = default_timer()
    para = [0]*10             #Initialize the internal parameters' list (bond length and bond angle)
    coord = readPDB(args[0])  #the order is A B C D AB AC AD BC BD CD
    accept = 0
    min_set_para = [0]*10
    min_set_coord = [0]*12
    ite = 0
    record = []
    para = car2in(para,coord)  
    E = give_me_sum(para)
    min = E
    print "Initial energy: ",E,"\n"
    while ite<=ite_max and min>E_tol:
      index = pick_next_permutation()
      pn = p_or_n()
      d = stepsize(step)
      coord[index] = coord[index]+ pn * d;
      para = car2in(para,coord);    
  
      if give_me_sum(para)<=E:
        E = give_me_sum(para)
        if E<min:
          min = E
          min_set_para=para[:]
          min_set_coord = coord[:]
      else:
        accept = accept + math.exp(-(give_me_sum(para)-E)/kT)
        if accept_or_not(math.exp(-(give_me_sum(para)-E)/kT)):
          E = give_me_sum(para)
        else:
          coord[index] = coord[index]- pn * d
      if ite%1000==0:
         record.append(min)
      ite +=1
    trail[kT] = record
    duration_time = default_timer() - start_time
    #--------------------------------------Show result
    print "\n===============RESULT================\n"
    print "parameters in internal coordination:\n"
    print 'Bond lengths:' ,para[:4]
    print 'Bond angles:' , para[4:]
    print "\n====================================\n"
    print "parameters in cartesian coordination:\n"
    print 'Atom 1:', coord[:3]
    print 'Atom 2:', coord[3:6]
    print 'Atom 3:', coord[6:9]
    print 'Atom 4:', coord[9:]
    print "\n====================================\n"
    print "Minimum Energy = ",min
    print "Execution time: ", duration_time," sec"
    print "Done at the", ite-1,"th step"
    print "\n====================================\n"
  """
  ax=plt.gca()
  plt.xlim(0,ite_max)
  for key in trail:
    x_axis = [1000*i for i in range(len(trail[key]))]
    plt.plot(x_axis,trail[key], label='kT = '+ str(key))
  plt.title('Minimum enery in MCMC')
  plt.grid()
  plt.xlabel('Iteration')
  plt.ylabel('Lowest energy')
  plt.legend()
  plt.show()
  """
  return 0

if __name__ == "__main__":
    main()