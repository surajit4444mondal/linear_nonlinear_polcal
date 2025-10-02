from casatasks import polcal,applycal,setjy,gaincal


S1='unpolarized_source_I10.ms'
S2='polarized_source_I50_Q1_U5_V0.ms'
S3='polarized_source_I50_Q5_U1_V0.ms'
S4='polarized_source_I50_Q4_U15_V0.ms'
S5='polarized_source_I50_Q15_neg_U5_V0.ms'

setjy(vis=S1,standard='manual',fluxdensity=[10,0,0,0],usescratch=True)
setjy(vis=S2,standard='manual',fluxdensity=[50,1,5,0],usescratch=True)

gaincal(vis=S1,caltable=S1[:-3]+".gcal",refant='0')
gain_cal=S1[:-3]+".gcal"

polcal(vis=S1,caltable=S1[:-3]+".leak",poltype='Df',refant='0',gaintable=gain_cal)

leak_cal=S1[:-3]+".leak"

polcal(vis=S2,caltable=S2[:-3]+".crossphase",poltype='Xf',refant='0',gaintable=[gain_cal,leak_cal])

crossphase_cal=S2[:-3]+".crossphase"

for msname in [S1,S2,S3,S4,S5]:
    applycal(vis=msname,gaintable=[gain_cal,leak_cal,crossphase_cal],calwt=[False,False,False])
    
