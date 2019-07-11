load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin	

	year = sprinti("%0.4i",yy)
	day = sprinti("%0.2i",dd)
	month = sprinti("%0.2i",mm)
	hour = sprinti("%0.2i",hh)
	t = sprinti("%0.2i",tt)
	
	prefix="UV_interp_"
	date=str_concat((/ year,"-",day,"-",month,"_",hour,":00:00.t",t /))
	suffix=".nc"

	filename = str_concat((/prefix,date,suffix/))

	dir = "/home/sp917/data/"

	filename_full = str_concat((/dir,filename/))

	fin = addfile(filename_full, "r")

	U = wrf_user_getvar(fin,"U",-1)
	V = wrf_user_getvar(fin,"V",-1)	
	
	powerG = wave_number_spc(U,"G")
	powerF = wave_number_spc(U,"F")
	
	printVarSummary(powerG)
	printVarSummary(powerF)

	wks = gsn_open_wks("png", "testplot")

	res = True
        res@trXLog = True
        res@trYLog = True
        res@tiXAxisString = "wavenumber (km)"
        res@tiYAxisString = "Spectrum"

 	sizes = dimsizes(powerG)

	print(sizes)
	k = new(sizes(2),float)

	do i=0,sizes(2)-1
		k(i) = i+1
	end do		

	Kolmogorov = a*(k^(-5.0/3.0))

	data = new( (/	4, sizes(2)/) , float)
	
	res@xyMarkerColors = (/"blue", "red", "green", "orange"/)	

	data(0,:) = Kolmogorov
	data(3,:) = powerF(0,0,:)
	data(2,:) = powerG(0,0,:)
	data(1,:) = powerG(0,0,:) - powerF(0,0,:)

	plot = gsn_csm_xy(wks, k, data(1,:), res)

end

