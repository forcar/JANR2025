      subroutine maid_lee(a1,a2,a3,a4,a5,
     1           theory_opt,channel_opt,resonance_opt,iopt,
     1           pol_elec,pol_targ,sigu,sigt,
     1           sigl,sigi,sige,siget,sigpt,g)

      implicit none
      
      real a1,a2,a3,a4,a5

      include 'mpintp.inc'
      include 'spp.inc'

      integer theory_opt, channel_opt, resonance_opt, n_call,iopt
      real pol_targ,pol_elec
      real sigu,sigt,sigl,sigi,g
      real sige,sigpt,siget
      character*132 fname

      data n_call /0/

      q2          = a1
      w           = a2
      epsilon     = a3
      csthcm      = a4
      phicm       = a5
      
      intopt      = iopt
      
      if (n_call.ne.0) goto 100

      method_helicity = 1

      if (theory_opt.eq.1) then
        if (channel_opt.eq.1) then
          call revinm('CLAS_PARMS',
     1    'spp_tbl/hlee-PPpi.tbl',data_file)
        elseif (channel_opt.eq.2) then
          call revinm('CLAS_PARMS',
     1    'spp_tbl/hlee-NPpi.tbl',data_file)
        elseif (channel_opt.eq.3) then 
          call revinm('CLAS_PARMS',
     1    'spp_tbl/hlee-PNpi.tbl',data_file)
        endif
      endif
      
      if (theory_opt.eq.2) then
        if (channel_opt.eq.1) then
          call revinm('CLAS_PARMS',
     1    'spp_tbl/maid98-PPpi_S.tbl',data_file)
        elseif (channel_opt.eq.2) then
          call revinm('CLAS_PARMS',
     1    'spp_tbl/maid98-NPpi.tbl',data_file)
        elseif (channel_opt.eq.3) then 
          call revinm('CLAS_PARMS',
     1    'spp_tbl/maid98-PNpi_S.tbl',data_file)
        endif
      endif
      

      if (theory_opt.eq.3) then
        if (channel_opt.eq.1) then
          call revinm('CLAS_PARMS',
     1    'spp_tbl/maid03-PPpi.tbl',data_file)
        elseif (channel_opt.eq.2) then
          call revinm('CLAS_PARMS',
     1    'spp_tbl/maid03-NPpi.tbl',data_file)
        elseif (channel_opt.eq.3) then 
          call revinm('CLAS_PARMS',
     1    'spp_tbl/maid03-PNpi.tbl',data_file)
        endif
      endif
      

      if (theory_opt.eq.4) then
        if (channel_opt.eq.1) then
          call revinm('CLAS_PARMS',
     1    'spp_tbl/maid07-PPpi.tbl',data_file)
        elseif (channel_opt.eq.2) then
          call revinm('CLAS_PARMS',
     1    'spp_tbl/maid07-NPpi.tbl',data_file)
        elseif (channel_opt.eq.3) then 
          call revinm('CLAS_PARMS',
     1    'spp_tbl/maid07-PNpi.tbl',data_file)
        endif
      endif
      

      if (theory_opt.eq.5) then
        if (channel_opt.eq.1) then
          call revinm('CLAS_PARMS',
     1    'spp_tbl/dmt-PPpi_all_on.tbl',data_file)
        elseif (channel_opt.eq.3) then 
          call revinm('CLAS_PARMS',
     1    'spp_tbl/dmt-PNpi_all_on.tbl',data_file)
        endif
      endif
      
      if (theory_opt.eq.6) then
        if (channel_opt.eq.4) then
          call revinm('CLAS_PARMS',
     1    'spp_tbl/etamaid-PPeta.tbl',data_file)
        endif
      endif
            
      if (theory_opt.ge.7) then
        print *, 'Enter spp_tbl/table: '
        read(5,'(a)') fname
        call revinm('CLAS_PARMS','spp_tbl/'//fname,data_file)
      endif      
            
      write(6,*) 'Enter max pi-N angular momentum (0-5): '

      read(5,*) mwave_L   

      write(6,*) 'Reading multipoles from ',data_file

      call read_sf_file(data_file,2)

      n_call = n_call + 1
      
100   call xsection
      
      sigu = sigma_t
      sigl = sigma_l
      sigt = sigma_tt
      sigi = sigma_lt  
     
      sigpt  = sigma_t_y+sigma_l_y
      sigpt  = sigpt+sigma_lt_x+sigma_lt_y+sigma_lt_z
      sigpt  = sigpt+sigma_tt_x+sigma_tt_y+sigma_tt_z      
      
      siget  = sigma_ltp_x+sigma_ltp_y+sigma_ltp_z
      siget  = siget+sigma_ttp_x+sigma_ttp_z
      
      sige   =  pol_elec*sigma_ltp  
      sigpt  =  pol_targ*sigpt
      siget  =  pol_elec*pol_targ*siget
      
      g      = -sigma_tt_z/sigma_t
      
      return
      end
