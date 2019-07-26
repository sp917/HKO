begin

    a = addfile("../data/X.nc", "r")

    U = a->X

    power = wave_number_spc(U,"G")

    newfile = addfile("../data/P.nc", "c")

    newfile->P = power

end
