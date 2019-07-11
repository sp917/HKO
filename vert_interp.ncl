load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin	

	year = sprinti("%0.4i",yy)
	day = sprinti("%0.2i",dd)
	month = sprinti("%0.2i",mm)
	hour = sprinti("%0.2i",hh)
	t = sprinti("%0.2i",tt)
	
	prefix="wrfout_xhka_"
	date=str_concat((/ year,"-",day,"-",month,"_",hour,":00:00.t",t /))
	suffix=".nc"

	filename = str_concat((/prefix,date,suffix/))

	dir = "/home/sp917/data/"

	filename_full = str_concat((/dir,filename/))

	fin = addfile(filename_full, "r")

	U = wrf_user_getvar(fin,"U",-1)
	V = wrf_user_getvar(fin,"V",-1)
	U10 = wrf_user_getvar(fin,"U10",-1)
	V10 = wrf_user_getvar(fin,"V10",-1)
	
	U_unstagger = wrf_user_unstagger(U,U@stagger)
	V_unstagger = wrf_user_unstagger(V,V@stagger)
;	U10_unstagger = wrf_user_unstagger(U10,U10@stagger)
;	V10_unstagger = wrf_user_unstagger(V10,V10@stagger)

	Udims = dimsizes(U_unstagger)
	nz = Udims(1)
 
	vert_coord = "ght_agl"
	interp_levels = fspan(0,16.6,nz)
	
	opts = True
	opts@extrapolate = True
	opts@field_type = "ght"

	U_intrp = wrf_user_vert_interp(fin, U_unstagger,vert_coord,interp_levels,opts)
	V_intrp = wrf_user_vert_interp(fin, V_unstagger,vert_coord,interp_levels,opts)

	xvals = wrf_user_getvar(fin,"XLONG",-1)
	yvals = wrf_user_getvar(fin,"XLAT",-1)

	newfile_name = str_concat((/ "UV_interp_",date,".nc"  /))
	newfile_name_full = str_concat((/dir,newfile_name/))

	system(str_concat((/ "/bin/rm ", newfile_name_full /)))

	newfile = addfile(newfile_name_full, "c")

	speed = sqrt(U_intrp^2 + V_intrp^2)
	ke = 0.5*(U_intrp^2 + V_intrp^2)
	speed10 = sqrt(U10^2 + V10^2)
	ke10 = 0.5*(U10^2 + V10^2)
	
	newfile->U = U_intrp
	newfile->V = V_intrp 
	newfile->XLONG = xvals
	newfile->XLAT = yvals
	newfile->U10 = U10 ;_unstagger
	newfile->V10 = V10 ;_unstagger
	newfile->speed = speed
	newfile->ke = ke
	newfile->speed10 = speed10
	newfile->ke10 = ke10

	U_spectrum = wave_number_spc(U_intrp,"F")
	V_spectrum = wave_number_spc(V_intrp,"F")
	speed_spectrum = wave_number_spc(speed,"F")
	ke_spectrum = wave_number_spc(ke,"F")
	U10_spectrum = wave_number_spc(U10, "F") ;_unstagger,"F")
	V10_spectrum = wave_number_spc(V10, "F") ;_unstagger,"F")
	speed10_spectrum = wave_number_spc(speed10,"F")
	ke10_spectrum = wave_number_spc(ke10,"F")
		
	newfile->U_spectrum = U_spectrum
	newfile->V_spectrum = V_spectrum
	newfile->U10_spectrum = U10_spectrum
	newfile->V10_spectrum = V10_spectrum
	newfile->speed_spectrum = speed_spectrum
	newfile->ke_spectrum = ke_spectrum
	newfile->speed10_spectrum = speed10_spectrum
	newfile->ke10_spectrum = ke10_spectrum
	
end

