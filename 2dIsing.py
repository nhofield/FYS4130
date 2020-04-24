import numpy as np
import matplotlib.pyplot as plt

dim = 2
#np.random.seed(10)
lattice_length = 16
mc_cycles = 100
J = 1.

lattice_size = lattice_length**dim
spin_possibilities = [1, -1]

lattice = np.random.choice(spin_possibilities, (lattice_length, lattice_length))
lattice_indices = range(lattice_length)
lattice_initial = np.copy(lattice)
print "Start Lattice"
print lattice

#Periodic Boundary conditions
periodic_conds = lambda x: (x % lattice_length) - \
lattice_length*((x % lattice_length) // ((lattice_length + 1)//2))

T_eq = lambda beta: 216.*np.exp(-1.62*beta) #Equilibration time formula from fitting
#T_Equilib = int(np.round(T_eq(beta)))
T_Equilib = 0
print "Equilibration time = ", T_Equilib

def wolff_update(lattice, mc_cycle, beta):

    #pick random site
    i_rand = np.random.choice(lattice_indices)
    j_rand = np.random.choice(lattice_indices)
    SITE = (i_rand, j_rand)
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
#    if np.all(lattice == -1):
#        lattice = lattice*(-1)
#    else:
#        pass
    m, mm = extract_data(lattice, J, beta)
#    plt.imshow(lattice)
#    print "Lattice after while loop complete"
#    print lattice
#    plt.show()
#    print "Completed ", counter , " iterations of while loop"
    return m, mm, mmmm



def add_neighbours(SITE, stack):
    list_add = []
    list_add.append( (periodic_conds( SITE[0] ), periodic_conds( SITE[1] + 1 ) ))
    list_add.append( (periodic_conds(SITE[0] - 1), periodic_conds(SITE[1]) ) )
    list_add.append( (periodic_conds( SITE[0] ), periodic_conds( SITE[1] - 1 ) ))
    list_add.append( (periodic_conds(SITE[0] + 1), periodic_conds(SITE[1]) ) )
    stack.append(list_add)


def extract_data(lattice, J, beta):
    m = 0.  #Magnetization
    mm = 0.
    mmmm = 0.
    for i in range(lattice_length):
        for j in range(lattice_length):
            sigma_ij = lattice[i,j]

            m += sigma_ij

    mm = m**2
    mmmm = mm**2
    return m/float(lattice_size), mm/float(lattice_size**2), mmmm/float(lattice_size)

energy = np.zeros(mc_cycles)
beta_0 = 0.001
beta_f = 3.0
number_of_betas = 100
beta = np.linspace(beta_0 , beta_f, number_of_betas)


beta_iterations = 0

m = np.zeros(number_of_betas)
mm = np.zeros(number_of_betas)

for i in range(number_of_betas):
    m_temp = np.zeros(mc_cycles)
    mm_temp = np.zeros(mc_cycles)

    print "beta = ", beta[i]
    for k in range(mc_cycles):
#        print "k = ", k
        m_temp[k], mm_temp[k] = wolff_update(lattice, mc_cycles, beta[i])
    print mm_temp
    m[i] = np.mean( m_temp )
    mm[i] = np.mean( mm_temp )
    beta_iterations += 1
    print "i = ", beta_iterations


lattice_lengths = np.array([8, 16, 32])

for i in range( len(lattice_lengths) ):



#print "Lattice before"
#print lattice_initial
#print "Lattice after"
#print lattice_equi
x_val = range(mc_cycles)
#plt.figure(figsize=(8,6))
#plt.plot(x_val, energy, label="Energy")
#plt.grid(True)
#plt.legend()

plt.figure( figsize = ( 8, 6 ) )
plt.subplot(1, 2, 1)
plt.plot(beta, m, label=r"$<m>$")
plt.xlabel(r"$\beta$")
plt.ylabel(r"$<m>$")
plt.grid(True)
plt.legend()
plt.subplot(1, 2, 2)
plt.plot(beta, mm, label=r"$<m^2>$")
plt.xlabel(r"$\beta$")
plt.ylabel(r"$<m^2>$")
plt.legend()
plt.grid(True)
plt.show()
