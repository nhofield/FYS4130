import numpy as np
import matplotlib.pyplot as plt

dim = 1
#np.random.seed(10)
lattice_length = 64
mc_cycles = 1000
J = 1.
#beta = 1.61
#beta = 0.61
beta = 2.0
lattice_size = lattice_length**dim
spin_possibilities = [1, -1]

lattice = np.random.choice(spin_possibilities, (1, lattice_length))[0]
lattice_indices = range(lattice_length)
lattice_initial = np.copy(lattice)
print "Start Lattice"
print lattice

#Periodic Boundary conditions
periodic_conds = lambda x: (x % lattice_length) - \
lattice_length*((x % lattice_length) // ((lattice_length + 1)//2))

T_eq = lambda beta: 216.*np.exp(-1.62*beta) #Equilibration time formula from fitting
alpha = -np.log( np.tanh(beta*J) )
C = np.vectorize(lambda r_value: np.exp( -alpha*r_value ))
#T_Equilib = int(np.round(T_eq(beta)))
T_Equilib = 0
#r_value = 6
print "Equilibration time = ", T_Equilib

def wolff_update(lattice, mc_cycle, r_val):

    #pick random site
    i_rand = np.random.choice(lattice_indices)
    SITE = i_rand
    SITES = [SITE]

    #Create stack and cluster
    lattice[SITE] = lattice[SITE]*(-1)
    stack_old = []

#    print "Initial lattice While loop"
#    print lattice
#    print lattice
#    print "Initial SITE ", SITE
    #Add neighbours to stack
    add_neighbours(SITE, stack_old)

#    print "Stack_old initial"
#    print stack_old
    while len(stack_old) != 0:    #Create cluster

        current_stack_list_size = len( stack_old[-1] )
#        print current_stack_list_size, "s len"
        counter = 0
        SITE_current = SITES[-1]
        rand_index = np.random.choice( range(current_stack_list_size) )
        current_neighbour = stack_old[-1][rand_index]
#        print current_neighbour, "Current neighbour"
#        print SITE_current, "Current SITE"
        if lattice[SITE_current] != lattice[ current_neighbour ]:

            if np.random.rand() <= 1.-np.e**(-2.*beta):
#                print "New SITE"
                new_SITE = current_neighbour
                stack_old[-1].pop( stack_old[-1].index(current_neighbour) )
                if stack_old[-1] == []:
                    stack_old.pop(-1)
                    SITES.pop(-1)
                else:
                    pass
                SITES.append(new_SITE)
                lattice[new_SITE] = lattice[new_SITE]*(-1)
                add_neighbours(new_SITE, stack_old)


            else: # | Antialigned, but not accepted

#                print "Antialigned, but not accepted"
                stack_old[-1].pop( stack_old[-1].index(current_neighbour) )
                if stack_old[-1] == []:
                    stack_old.pop(-1)
                    SITES.pop(-1)
                else:
                    pass


        else: #| Not antialigned |
#            print "Not antialigned"
            stack_old[-1].pop( stack_old[-1].index(current_neighbour) )
            if stack_old[-1] == []:
                stack_old.pop(-1)
                SITES.pop(-1)
            else:
                pass


#        print "Completed"
#        print "Current stack_old"
#        print stack_old
#        print "Current SITES"
#        print SITES

    energy, m = extract_data(lattice, J, beta)
#    plt.imshow(lattice)
#    print "Lattice after while loop complete"
#    print lattice
#    plt.show()
#    print "Completed ", counter , " iterations of while loop"
    return lattice, energy, m, lattice[0]*lattice[r_val]



def add_neighbours(SITE, stack):
    list_add = []
    list_add.append( periodic_conds( SITE + 1 ) )
    list_add.append( periodic_conds( SITE - 1 ) )
    stack.append(list_add)



def correlation_of_spins(spin_0, spin_r):
    corr = 0
    for i in range(lattice_length):
        corr += spin_0*spin_r



    return corr

def extract_data(lattice, J, beta):
    H = 0.  #Energy
    m = 0.  #Magnetization
    for i in range(lattice_length):
        sigma_i = lattice[i]
        H += -J*sigma_i*lattice[ periodic_conds(i + spin_possibilities[0] ) ]
        H += -J*sigma_i*lattice[ periodic_conds(i + spin_possibilities[1] ) ]
        m += sigma_i

    return H, m/float(lattice_size)

energy = np.zeros(mc_cycles)
m = np.zeros(mc_cycles)


#for k in range(mc_cycles):
#    print "k = ", k
#    lattice_equi, energy[k], m[k], correlation_0r[k] = wolff_update(lattice, mc_cycles, 1)
correlation_0r = np.zeros(lattice_length-1)
for r_vals in range(1, lattice_length):

    correlation_0r_temp = np.zeros(mc_cycles)
    for k in range(mc_cycles):
        lattice_equi, energy[k], m[k], correlation_0r_temp[k] = wolff_update(lattice, mc_cycles, r_vals)
    correlation_0r[r_vals-1] = np.mean(correlation_0r_temp)

print correlation_0r

print correlation_0r
print "Lattice before"
print lattice_initial
print "Lattice after"
print lattice_equi
slic = 1

print "-----------------"
print "Analytical Values"
print "-----------------"
m_analytical = lambda J, beta_val: (np.tanh(J*beta_val))**lattice_size
print "m_analytical = ", m_analytical(J, beta)
print "m_wolff = ", m[-1]


x_val = np.array(range(mc_cycles))[T_Equilib::slic]
#print x_val
#print x_val, "x val"
p = np.polyfit(x_val, energy[T_Equilib::slic], 1)
plt.figure(figsize=(10,6))
plt.subplot(1, 2, 1)
#plt.plot(x_val, p[1] + p[0]*x_val)
plt.plot(x_val, energy[T_Equilib::slic], label="Energy")
plt.grid(True)
plt.legend()
plt.subplot(1, 2, 2)
print extract_data(lattice, J, beta)[1], " m "
plt.plot(x_val, m[T_Equilib::slic], label="Magnetization")
plt.grid(True)
plt.legend()

plt.figure(figsize=(10,6))
x_corr_vals = np.array(range( len(correlation_0r) ) )
plt.scatter( x_corr_vals , correlation_0r, label = "Numerical" )
plt.plot( x_corr_vals, C(x_corr_vals) , label="C(r)")
plt.ylabel(r"$<\sigma_0 \sigma_r>$")
plt.xlabel("MC-Cycles")
plt.grid(True)
plt.legend()
plt.show()
