include("ising_aux.jl")
Tmin = 0.001
Tchange = 0.1
Tmax = 5.0
Temperature = Tmin:Tchange:Tmax

N = 4

E_vec = zeros(length(Temperature),1)
Mabs_vec = zeros(length(Temperature),1)
count = 1

for T in Temperature
    e_tot = 0.0
    m_tot = 0.0
    part_tot = 0.0
    for i = 1:2^(N*N)
        lat = bin_init(i-1,N)
        e = total_energy(lat)
        mabs = abs.(total_mag(lat))
        m_tot = m_tot + mabs*exp(-e/T)
        e_tot = e_tot + e*exp(-e/T)
        part_tot = part_tot + exp(-e/T)
    end
    E_vec[count] = Float64(e_tot/part_tot/16)
    Mabs_vec[count] = Float64(m_tot/part_tot/16)
    count = count + 1
    println(T," ",e_tot/part_tot/16)
end
