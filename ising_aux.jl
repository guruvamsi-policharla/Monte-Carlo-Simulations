function initialise(N::Int)
    """ Initialising the lattice with random values """
    #lat = Array{Float64, 3}(N, N)
    lat = Array{Vector{Float64},2}(N, N);
    for i = 1:N
        for j = 1:N
            x = [rand(Uniform(0,1)), rand(Uniform(0,1)), rand(Uniform(0,1))]
            x = x/sqrt(vecdot(x,x))
            lat[i,j] = convert(Vector,x)
        end
    end
    return lat
end

function energy_pos(x, y, lat, a = [0,0,0])
    """ Calculating energy of the position """
    N = size(lat,1)
    if(y==N)
        up=1;
    else
        up=y+1;
    end

    if(y==1)
        down=N;
    else
        down=y-1;
    end

    if(x==1)
        left=N;
    else
        left=x-1;
    end

    if(x==N)
        right=1;
    else
        right=x+1;
    end

    if a == [0,0,0]
        energy = -1*vecdot(lat[x,y],(lat[left,y]+lat[right,y]+lat[x,up]+lat[x,down]));
        return energy
    else
        energy = -1*vecdot(a,(lat[left,y]+lat[right,y]+lat[x,up]+lat[x,down]));
        return energy
    end
end

function test_flip(x, y, lat, T)
""" Checks whether energy allows for a flip or not """
    theta = rand(Uniform(0,2*pi))
    phi = rand(Uniform(0,2*pi))
    a = convert(Vector,[sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)])
    de = -energy_pos(x,y,lat) + energy_pos(x,y,lat,a);

    if(de<0)
        lat[x,y] = a
        return true
    elseif(rand()<exp(-de/T))
        lat[x,y] = a
        return true
    else
        return false
    end
end

'''
function flip(x,y,lat)
    N = size(lat,1)
    theta = rand(Uniform(0,2*pi))
    phi = rand(Uniform(0,2*pi))

    #Flipping randomly
    lat[x,y] = [sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)]
end
'''

function transient_results(lat, transient::Int, T)
    """Takes lat as input and removes initial transients by running for transient number of steps"""
    N = size(lat,1)
    for i = 1:transient
        for j = 1:N*N
                x = rand(1:N)
                y = rand(1:N)
                test_flip(x,y,lat,T)
        end
    end
end

function total_mag(lat)
    return sum(lat)
end

function total_energy(lat)
    e = 0.0
    for i = 1:size(lat,1)
        for j = 1:size(lat,2)
            e = e + energy_pos(i,j,lat)
        end
    end
    return e/2
end

function bin_init(val,N)
    val_s = bin(val,N*N)
    lat = zeros(N*N,1)
    for i = 1:N*N
        if (val_s[i] == '0')
            lat[i] = -1
        elseif (val_s[i] == '1')
            lat[i] = 1
        end
    end
    lat = reshape(lat,(4,4))
    return lat
end
