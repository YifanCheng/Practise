
PROGRAM two_layer_adv_diff
!
!------------------------------------------------------------------------- 
!           Practice 2-layer Sream Temperature  two_layer model
!                 by Ryan Niemeyer, January 2016
!                 modified by Yifan Cheng, Jan 24 2016		
!     based on Strzepek et al. (2015) equations 6 & 7, which is based
!      off of Chapra 1997
!
!	Update the equations calculated here and take advection into consideration
!
!-------------------------------------------------------------------------

implicit none

!-------------------------------------------------------------------------
!    define variables
!-------------------------------------------------------------------------

! 1-dim. real array for 10 years
real, dimension(365*10) :: flow_in,flow_out,flow_ex, flow_Tin
real, dimension(365*10) ::  temp_epil,temp_hypo, temp_out_tot
real, dimension(365*10) :: temp_change_ep, temp_change_hyp, energy
real, dimension(365*10) :: energy_tot, diffusion_tot,advection_tot , T_in_tot, T_out_tot

REAL, PARAMETER :: Pi = 3.1415927
real, parameter :: v_t = 2.1 !diffusion coefficient (m2/day) - based on Snodgrass, 1974
real  :: flow_constant
integer  :: i
real :: x1, x2, x3

real  :: depth_total, depth_e, depth_h, width, length, volume_e_x, outflow_x
real :: energy_x, volume_h_x, area, density, heat_c, temp_change
real  :: flow_in_hyp_x, flow_in_epi_x, flow_out_epi_x, flow_out_hyp_x
real  :: epix, hypox, dif_epi_x, dif_hyp_x, x, flow_epi_x, flow_hyp_x
real  :: adv_epi_x, adv_hyp_x

depth_total = 30
width = 200
length = 17000
area = width*length

!-------------------------------------------------------------------------
!     generate flow and energey temporal variables
!-------------------------------------------------------------------------

! --------------------------- flow -------------------------

! --------- constant flow paramaeter -------
! constant  to change day of
!  flow to go "up and down" with sin wave
! with 0 and 365 being "0" point
flow_constant = 365/(2*Pi)

! generates flow eacy day as a sin wave with the peak flow on April 1
! at 90000 cfs, and lowest point is 30000 cfs on October 1
! days are calendar year (day = 1 = January 1)


do  i=0,10*365

   ! ------ get flow in to the reservoir for the year ------
   flow_in(i) =  sin( i/flow_constant)
   flow_in(i) = (flow_in(i) + 2)*30000   ! gets flow to vary from 30000 to 90000 cfs

   ! ------ get flow in to the reservoir for the year ------
   flow_out(i) =  sin( i/flow_constant)
   flow_out(i) = (flow_out(i) + 2)*30000   ! gets flow to vary from 30000 to 90000

   ! ------ get flow exchange between the two layers for the year ------
   ! ------ here we assume the water exchange between the two layers equal to outflow -----
   flow_ex(i) = flow_out(i)

end  do

! print *, flow(366:456)

! --------------------------- flow - temperature   -------------------------

do  i=1,10*365
   flow_Tin(i)  =  cos((i/flow_constant)+ Pi)
   flow_Tin(i) = (flow_Tin(i) + 1.5)*10   ! gets temperature to vary form 5 to 25C
end  do

! print *, "Flow temperature from day 1 to 200"
! print *, flow_Tin(1:200)

! print *, "Flow temperature from day 3000 to 3200"
! print *, flow_Tin(3000:3200)

! --------------------------- energy  -------------------------

! --------- constant flow paramaeter -------
!  energy in W/m2
do  i=0,10*365
   energy(i) =  cos((i/flow_constant)+ Pi)
  ! gets net energy (positive downward) to vary from:
  ! " +0.7)*120":  -36 to 204 W/m2
  ! " +0.5)*120":  -60 to 180 W/m2
  ! converts energy from W/m2 to Joules/m2 * day
   energy(i) =( (energy(i) + 0.5)*120) * (24*60*60)
end  do

! print *, "energy from 2000 to 2100"
! print *, energy(2000:2100)

!-------------------------------------------------------------------------
!              Calculate two-layer stream temperature   
!-------------------------------------------------------------------------

density = 1000 !  density of water in kg / m3
heat_c = 4180  !  heat capacity of water in joules/ kg * C

! initial variables
depth_total = 30
depth_e = 5
depth_h = 15
volume_e_x = area*depth_e
volume_h_x = area*depth_h
temp_epil(1) = 5 ! starting epilimnion temperature at 5 C
temp_hypo(1) = 5 ! starting hypolimnion temperature at 5 C

! start at 2, because need to have an epil and hypo temperature to start out
do  i=2,10*365

  ! calculate incoming net energy (J/m2*day) to epilimnion
  energy_x = (energy(i)*area)/(density*heat_c*volume_e_x)

  ! divide incoming flow into reservoir to epilimnion and hypolimnion
  flow_in_hyp_x = flow_ex(i)
  flow_in_epi_x = flow_in(i)

  flow_out_hyp_x = flow_out(i)
  flow_out_epi_x = flow_ex(i)

  ! calculate temperature change due to diffusion
  ! NOTE: don't need to multiply by heat capacity or density of water because
  !       energy component is divided by those two
  dif_epi_x  = ( v_t * area * (temp_hypo(i-1) - temp_epil(i-1)))/volume_e_x * 24 * 60 * 60
  dif_hyp_x  = ( v_t * area * (temp_epil(i-1) - temp_hypo(i-1)))/volume_h_x * 24 * 60 * 60

  !calculate temperature change due to advection
  ! 
  adv_epi_x = ( flow_out_epi_x * (temp_hypo(i-1)-temp_epil(i-1)) )/volume_e_x*24*60*60
  adv_hyp_x = ( flow_in_hyp_x * (temp_epil(i-1)-temp_hypo(i-1)) )/volume_e_x*24*60*60

  ! calculate change in EPILIMNION  temperature (celsius)
  flow_epi_x = flow_in(i)*(flow_Tin(i)-temp_epil(i-1))/volume_e_x*24*60*60
  temp_change_ep(i) = flow_epi_x + energy_x + dif_epi_x + adv_epi_x

  ! save each temperature change component
  energy_tot(i) = energy_x
  diffusion_tot(i) = dif_epi_x
  T_in_tot(i) = flow_in_epi_x*flow_Tin(i)/volume_e_x*24*60*60
  T_out_tot(i) = flow_out_epi_x*temp_epil(i-1)/volume_e_x*24*60*60
  ! update epilimnion volume for next time step
  volume_e_x = volume_e_x + (flow_in_epi_x - flow_out_epi_x)*24*60*60
  temp_epil(i) = temp_epil(i-1) +  temp_change_ep(i)

  ! calculate change in HYPOLIMNION  temperature (celsius)
  temp_change_hyp(i) = adv_hyp_x  +  dif_hyp_x   ! diffusion 

  ! update hypolimnion volume for next time step
  volume_h_x = volume_h_x + (flow_in_hyp_x - flow_out_hyp_x)*24*60*60
  temp_hypo(i) = temp_hypo(i-1) +  temp_change_hyp(i)

  ! calculate combined (hypo. and epil.) temperature of outflow
  outflow_x = flow_out_hyp_x
  temp_out_tot(i) = temp_hypo(i)

 if (i==2 .or. i==180 .or. i==370 .or. i==550 .or. i==720 .or. i==3649) then
!  print *, "run: ", i
!  print *,"energy of incoming radiation " , energy(i)/(60*60*24)
!  print *, "energy - joules/day*m2", energy(i)
!  print *, "area:                 ", area
!  print *, "numeriator  of energy equation: ", x1
!  print *, "heat capacity  ", heat_c
!  print *, "density   ", density
!  print *, "volume of epilimnion: ", volume_e_x
!  print *, "denominator of energy equation: ", x2
!  print *, "volume of epilimnion: ", volume_e_x
!   print *, "change in volume - epilim.: ", flow_in_epi_x -  flow_out_epi_x
!   print *, "change in volume - hypolim.: ", flow_in_hyp_x - flow_out_hyp_x
!   print *, "depth of epilimnion: ", volume_e_x/area
!   print *, "depth of hypolimnion: ", volume_h_x/area
!   print *, "energy temp change is : ",energy_x
!   print *, "diffusion temp change  is: ", dif_epi_x
!   print *, "flow epilim.  temp change is: ", flow_epi_x
!   print *, "flow hypolim.  temp change is: ", flow_hyp_x
!   print *, "depth of epilimnion: ", volume_e_x/area 
!   print *, "temperature change in epilimnion:  ", temp_change_ep(i)
!   print *, "outflow temperature from epilimnion is: ", temp_epil(i)
!   print *, "depth of hypolimnion: ", volume_h_x/area
!  print *, "temperature change in hypolimnion:  ", temp_change_hyp(i)
!  print *, "outflow temperature from hypolimnion is: ", temp_hypo(i)
!   print *, "outflow (combined)  temperature is: ", temp_out_tot(i)
!  print *, "  "
 end if

end  do

!-------------------------------------------------------------------------
!         Calculations to print at the end 
!-------------------------------------------------------------------------

 print *,"sum of 1:720 temperature change in epilim.: ", sum(temp_change_ep(1:720))
!  print *,"sum of 1:720 temperature change in hypolim.: ",  sum(temp_change_hyp(1:720))

 print *, "sum of 1:720 temp change due to flow in: ", sum(T_in_tot(1:720))
 print *, "sum of 1:720 temp change due to flow out: ", sum(T_out_tot(1:720))
 print *, "sum of 1:720 temp change due to energy: ", sum(energy_tot(1:720))
 print *, "sum of 1:720 temp change due to diffusion: ", sum(diffusion_tot(1:720))


end program two_layer_adv_diff
