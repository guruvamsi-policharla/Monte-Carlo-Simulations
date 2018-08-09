#Includes
include("ising_aux.jl")
using Distributions
using Plots

#Code
Tmin = 0.001
Tchange = 0.1
Tmax = 5.0
mcs = 100000
N = 4


lat = initialise(N)

norm=(1.0/float(mcs*N*N))

Temperature = Tmin:Tchange:Tmax
E_vec = zeros(length(Temperature),1)
M_vec = zeros(length(Temperature),1)
acc_vec = zeros(length(Temperature),1)

count = 1
for T in Temperature

    transient_results(lat,1000,T)
    M = total_mag(lat)
    E = total_energy(lat)

    etot=0;
    mabstot=0;
    acc_rat = 0;

    for i in 1:mcs
        for j in 1:N*N
            x = rand(1:N)
            y = rand(1:N)
            E_0 = energy_pos(x,y,lat)
            M_0 = lat[x,y]
            if(test_flip(x,y,lat,T))
                acc_rat = acc_rat + 1
                E = E + energy_pos(x,y,lat) - E_0
                M = M + lat[x,y] - M_0
            end
        end
        etot= etot + E
        mabstot = mabstot + vecdot(M,M)
    end

    acc_rat = acc_rat*norm
    M_vec[count] = mabstot*norm
    E_vec[count] = etot*norm
    acc_vec[count] = acc_rat
    println(T," ",M_vec[count])
    count = count + 1
end

#jackknife analysis
