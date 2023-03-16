; These files are loaded by default in NCL V6.2.0 and newer
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;**************************************
begin
;****************************************

    f1 = addfile("../all_data/extend_h4_full.cam_nov_to_mar.nc","r")
    f2 = addfile("../all_data/ex_02.cam_may_to_sep.nc","r")
    f3 = addfile("../all_data/minu_02.cam_may_to_sep.nc","r")
     f4 = addfile("../all_data/extend_h4_full.cam.ana.ave.nc","r") 
    f5 = addfile("../monnnnsooon/extend_h4_full.cam_may_to_sep.nc","r")  
    f6  = addfile("../vege_data/ex_02.cam_nov_to_mar_new.nc","r")
    f7 = addfile("../all_data/minu_02.cam_nov_to_mar.nc","r")

   precp_minu_sum = f3->PRECL(0,:,:) + f3->PRECC(0,:,:)
   precp_minu_win = f7->PRECL(0,:,:) + f7->PRECC(0,:,:)
   prep_mis3_win = f1->PRECL(0,:,:) + f1->PRECC(0,:,:)
    prep_mis3_sum = f5->PRECL(0,:,:) + f5->PRECC(0,:,:)
 
   sum_diff = (precp_minu_sum - prep_mis3_sum)* 60*60*24*1000
   win_diff = (precp_minu_win - prep_mis3_win)* 60*60*24*1000

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   ts_mis3 = f5->PRECL(0,:,:) + f5->PRECC(0,:,:)
   ts_e02 = f2->PRECL(0,:,:) + f2->PRECC(0,:,:)
   ts_m02 = f6->PRECL(0,:,:) + f6->PRECC(0,:,:)

   ts_mis3 = ts_mis3* 60*60*24*1000
   ts_e02 = ts_e02* 60*60*24*1000
   ts_m02 = ts_m02* 60*60*24*1000
   ppp = f1->PRECL(0,:,:) + f1->PRECC(0,:,:)
   ppp = ppp* 60*60*24*1000

  ts_diff_e = ts_e02 - ts_mis3
  ts_diff_m = ts_m02 - ppp


   ccc = f1->TS(0,:,:)

  copy_VarCoords(ccc,ts_diff_e)    ;positive +0.2 summer
  copy_VarCoords(ccc,ts_diff_m)    ;positive +0.2 winter
   copy_VarCoords(ccc,sum_diff)    ;negative -0.2 summer
   copy_VarCoords(ccc,win_diff)    ;negative -0.2 winter

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


     lgm_c_sp_sum_pp = f5->PRECC(:,:,:) + f5->PRECL(:,:,:)
   lgm_c_sp_win_pp = f1->PRECC(:,:,:) + f1->PRECL(:,:,:)
   lgm_c_sp_ann = f4->PRECC(:,:,:) + f4->PRECL(:,:,:)

       avg_lgm_c_sp_sum_pp = dim_avg_n_Wrap(lgm_c_sp_sum_pp,0)
  avg_lgm_c_sp_win_pp = dim_avg_n_Wrap(lgm_c_sp_win_pp,0)
  avg_lgm_c_sp_ann = dim_avg_n_Wrap(lgm_c_sp_ann,0)

;;;;;;;;;;;;;;;;AVERGAE OVER SUMMER WINTER AND ANNUAL;;;;;;;;;;;;;


  lgm_c_sp_sum_total = avg_lgm_c_sp_sum_pp * 60*60*24*1000*150
  lgm_c_sp_win_total = avg_lgm_c_sp_win_pp * 60*60*24*1000*150
  lgm_c_sp_annul_total = avg_lgm_c_sp_ann * 60*60*24*1000*365

  lgm_c_sp_diff =  lgm_c_sp_sum_total - lgm_c_sp_win_total
  lgm_c_sp_sum_ratio = lgm_c_sp_sum_total/lgm_c_sp_annul_total
  lgm_c_sp_win_ratio  = lgm_c_sp_win_total/lgm_c_sp_annul_total


;;;;;;;;;;;;;;;;;;;CALCULATE SUM WIN DIFF AND RATIO;;;;;;;;;;;;;;;

 fun_array = new((/48,96/),"float")
 fun_array(0:23,0:95) = -1.0
 fun_array(24:47,0:95) = 1.0 
  
 abs_lgm_c_sp_diff = lgm_c_sp_diff*fun_array

  lgm_c_sp_ratio = new((/48,96/),"float")
  lgm_c_sp_ratio(0:23,0:95) = lgm_c_sp_win_ratio(0:23,0:95)
  lgm_c_sp_ratio(24:47,0:95) = lgm_c_sp_sum_ratio(24:47,0:95)

;;;;;;;;;;;WELL,SORT THE RATIO AND DIFFERENCE;;;;;;;;;;;;;;;;;;;;

  lgm_c_sp_gm  = new((/48,96/),"float")


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

 
  ff1 = addfile("../monthly_data/minu_02.cam.jun.nc","r")
  ff2 = addfile("../monthly_data/minu_02.cam.jul.nc","r")
  ff3 = addfile("../monthly_data/minu_02.cam.aug.nc","r") 
 
  ts = ff1->PRECC(0,:,:) + ff1->PRECL(0,:,:) + ff2->PRECC(0,:,:) + ff2->PRECL(0,:,:) + ff3->PRECC(0,:,:) + ff3->PRECL(0,:,:)
  printVarSummary(ts)

 precip =  ts /3    ;convert to mm/day   
 
      do i = 0,47
  do j = 0,95

   if(lgm_c_sp_ratio(i,j) .gt. 0.55 .and. abs_lgm_c_sp_diff(i,j) .gt. 300) then
; lgm_c_sp_gm(i,j) = 1 
   precip(i,j) = precip(i,j) *1000 * 60*60*24
    else
; lgm_c_sp_gm(i,j) = lgm_c_sp_gm@_FillValue
    precip(i,j) = lgm_c_sp_gm@_FillValue
     end if

    end do
   end do

   copy_VarCoords(ccc,lgm_c_sp_gm)
   copy_VarCoords(ccc,precip)

;  g_pre = precip *1000 * 60*60*24 *lgm_c_sp_gm

 ; print(max(g_pre))

; printVarSummary(g_pre)
;  copy_VarCoords(ccc,g_pre)

;   g_pre@_FillValue = 0



  gw = f1->gw                            ; get gaussian weights for ave


printVarSummary(gw)
;****************************************
; calculate averages
;****************************************
; arguments:
; data
; weights in y direction
; weights in x direction (we just use 1.0 since there are none)
; option (0 means use missing values, 1 means the ave will be missing if
;    ANY data in the array is missing.

  globav  = wgt_areaave(precip,gw,1.0,0)
  nh = wgt_areaave(precip({0:90},:),gw({0:90}),1.0,0)
  sh  = wgt_areaave(precip({-90:0},:),gw({-90:0}),1.0,0)
 
;   print("global mean = "+ globav)
;   print("HN mean = "+ nh)
;   print("SN mean = "+ sh)
 

   lat = ccc&lat
   lon =ccc&lon

   nlat = dimsizes(lat)
    gwty  = new((/nlat/),double)

gwty  = latGauWgt(nlat, "lat", "gaussian weights", "")


   globav1  = wgt_areaave(precip({-90:90},:),gwty({-90:90}),1.0,0)
  nh1 = wgt_areaave(precip({0:90},:),gwty({0:90}),1.0,0)
  sh1  = wgt_areaave(precip({-90:0},:),gwty({-90:0}),1.0,0)
 
   print("global mean = "+ globav1)
   print("HN mean = "+ nh1) 
   print("SN mean = "+ sh1) 


  wks = gsn_open_wks("pdf","ttttt")         ; send graphics to PNG file
  gsn_define_colormap(wks,"testcmap")  

  res = True
  res@mpFillOn     = False
  res@mpOutlineOn  = True                     ; turn the map outline on


;p1 = gsn_csm_contour_map(wks,g_pre,res)

end





