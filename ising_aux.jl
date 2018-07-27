function initialise(N::Int)
    """ Initialising the lattice with random values """
lat = zeros(N,N)
    for i = 1:N
        for j = 1:N
            if (rand()>=0.5)
                lat[i,j] = 1
            else
                lat[i,j] = -1
            end
        end
    end
    return lat
end

function energy_pos(x, y, lat)
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

    energy = -1*lat[x,y]*(lat[left,y]+lat[right,y]+lat[x,up]+lat[x,down]);
    return energy
end

function test_flip(x, y, lat, T)
""" Checks whether energy allows for a flip or not """
    de = -energy_pos(x,y,lat);

    if(de<0)
        flip = true;
    elseif(rand()<exp(-2*de/T))
        flip = true;
    else
        flip = false;
    end
    return flip
end

function transient_results(lat, transient::Int, T)
    """Takes lat as input and removes initial transients by running for transient number of steps"""
    N = size(lat,1)
    for i = 1:transient
        for j = 1:N*N
            x = rand(1:N)
            y = rand(1:N)
            if(test_flip(x,y,lat,T) == true)
                lat[x,y] = -lat[x,y];
            end
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
