load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"

begin

	dir = "/home/sp917/"
	datadir = str_concat((/dir,"data/"/))
	plotdir = str_concat((/dir,"plots/"/))

	year = sprinti("%0.4i",yy)
        day = sprinti("%0.2i",dd)
        month = sprinti("%0.2i",mm)
        hour = sprinti("%0.2i",hh)
        t = sprinti("%0.2i",tt)

	if (var .eq. 1) then
		variable = "U"
	else if (var .eq. 2) then
		variable = "V"
	else if (var .eq. 3) then
		variable = "ke"
	else if (var .eq. 4) then
		variable = "speed"
	else
		variable = "U"
	end if
	end if
	end if
	end if

	prefix=str_concat((/"Spectrum_",variable, "_" /))
        date=str_concat((/ year,"-",day,"-",month,"_",hour,":00:00.t",t /))
        suffix=".nc"

        filename = str_concat((/prefix,date,suffix/))

	filename_full = str_concat((/datadir,filename/))

        fin = addfile(filename_full, "r")

	S = fin->spectrum
	K = fin->k

	Kplot = K(1:)
	Splot = S(0,ll,1:)
	
	Kolmogorov = aa*(Kplot^(-5.0/3.0))

	wks = gsn_open_wks("png", plotdir + "spectrum_" + variable + "_layer_" + tostring(ll) + "_" + date)	

	res = True
	res@gsnDraw = False
	res@gsnFrame = False
	res@trXLog = True
	res@trYLog = True
	res@tiXAxisString = "wavenumber (km)"
	res@tiYAxisString = "Spectrum"

	
	
	plot1 = gsn_xy(wks,Kplot,Splot,res)
	plot2 = gsn_xy(wks,Kplot,Kolmogorov,res)

	overlay(plot1,plot2)

	draw(plot1)
	frame(wks)
	


end
