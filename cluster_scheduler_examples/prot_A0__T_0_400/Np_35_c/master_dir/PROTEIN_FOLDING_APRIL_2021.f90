! 			PROTEIN MEMBRANE
! 			  April 2020

!DR: (27.11.23)
! This version of VB program was modified by me to allow for command-line input (in particular random seed input).
! Moreover "num_max_aminoacids" was reduced form 21 to 20 to allow for the employment of aapot.dat file from 2016
! (without membrane bead and without removing the two errors in the MJ matrix of interactions).
! No other modifications wrt the original version inherited from J Bauermann.

!DR: (14.12.23)
! There is an additional modification by me in subroutine protein_rotation() to avoid flux coming from rounding
! of the center of mass position

!DR: (26.01.24)
! Modified the "flux" flag option to resemble 2016 code version
! Added the 'status' saving part
! Added the 'Bonds_DR' part

!DR: (28.02.24)
! Added saving of the cooperative bonds in the 'Bonds_DR' part [Nim, Nim_phob]
! Added saving of the number of w-p bonds of type PHI and PHO in the thermo output [bond_phob, bond_phil]

!DR: (16.04.24)
! Corrected LJ neighbor

!DR: (17.04.24)
! Corrected LJ neighbor

! Water neighbors: "0" down, "1" right, "2" up, "3" left
module	common_variables_folding

   REAL(KIND=8) :: DR_HB_f, DR_HB_beta = 30.0, DR_HB_l_threshold = 2.0d0**(1.0d0/3.0d0) !+10.0d0

   INTEGER, parameter :: rings=3, num_max_aminoacids=20, q=6    !number of shells and potts variables

   ! the number of amino acids are 21, the one with id 21 is a fake particle used to represent the membrane beads

   INTEGER	:: L, half_L, num_site, HB_bulk_local, HB_phil_local, HB_phob_local=0, HB_mix_local=0, number_species
   INTEGER :: Nim_local, Nim_phob_local, round_part, number_proteins, max_length, tot_prot_length
   INTEGER :: HB_bulk, HB_phob, HB_phil, HB_mix, bond_phob, bond_phil, Nim, Nim_phob, cp, cp_ran, pp_cp
   INTEGER :: flux
   INTEGER(KIND=8) :: acc_rotation, trial_rotation, acc_pivot, trial_pivot, acc_volume, trial_volume
   INTEGER(KIND=8) :: acc_crankshaft, trial_crankshaft, acc_jump, trial_jump, acc_shift, trial_shift
   INTEGER(KIND=8) :: acc_protein_rotation, trial_protein_rotation
   REAL(KIND=8) :: Temp, Press, lsample
   REAL(KIND=8) :: jj_bulk, jj_phob, jj_phil, jj_mix, jsigma, jsigma_phob
   REAL(KIND=8) :: jj_bulk_eff, jj_phob_eff, jj_phil_eff, jj_mix_eff, E_LJ, cp_max
   REAL(KIND=8) :: latt_size, latt_size2, latt, latt2, v12, v6, dL2
   REAL(KIND=8) :: energy, vol, vHB_bulk, vHB_phob, vHB_phil, vHB_mix, EPS, ratio

   INTEGER(KIND=8), dimension(0:3) :: cx, cy, neig_arm
   REAL(KIND=8), dimension(-num_max_aminoacids:6,-num_max_aminoacids:6) :: matrix

   INTEGER, dimension (:), allocatable :: protein_length, protein_move
   INTEGER, dimension (:), allocatable :: native_contacts_per_protein, non_native_contacts_per_protein, inter_contacts_per_protein
   INTEGER, dimension (:), allocatable :: total_contacts_per_protein, cp_max_per_protein
   INTEGER, dimension (:), allocatable :: native_contacts_per_species, non_native_contacts_per_species, inter_contacts_per_species
   INTEGER, dimension (:), allocatable :: total_contacts_per_species, cp_max_per_species


   INTEGER, dimension (:,:), allocatable :: protein_state, px, py
   INTEGER, dimension (:), allocatable :: pbc, rx, ry, protein_species
   INTEGER, dimension (:,:), allocatable :: displacement
   INTEGER, dimension (:,:), allocatable :: surface, surface_old,  protein_old, protein_aux
   INTEGER, dimension (:,:,:), allocatable :: status, status_old, protein, neigh, neigh_run


end module	common_variables_folding


module  random_variables_folding

   INTEGER(KIND=8) :: ic
   INTEGER, parameter :: ip=1279
   INTEGER(KIND=8), dimension(1:ip) :: ix

end module random_variables_folding

!*********
!* START *
!*********

program protein_design_folding

   use common_variables_folding

   implicit none
   INTEGER :: init_time, sample, conf_sample !initial MC step to equilibrate, final MC step in which I sample, timesteps between two samples
   INTEGER :: i, j, x, y, n, particle, time, st, nn
   INTEGER :: run_time, choose_prot_move, random_steps
   INTEGER(KIND=8) :: i_dran, iseed
   REAL(KIND=8) :: compress_bulk, compress_phob, compress_phil, prot_ene
   CHARACTER*100 :: word
!DR: added this variable
   CHARACTER(len=32) :: arg
   
   


   open(1133,file="hor_HB.dat",status="replace")
   open(1134,file="ver_HB.dat",status="replace")
   open(1135,file="coop_bonds.dat",status="replace")
   open(1136,file="coop_bonds_PHO.dat",status="replace")
   close(1133)
   close(1134)
   close(1135)
   close(1136)

   open(10,file="input_data_folding",status="unknown")
   read(10,*) word, L
   read(10,*) word, Temp
   read(10,*) word, Press
   read(10,*) word, number_proteins
   read(10,*) word, init_time
   read(10,*) word, run_time
   read(10,*) word, sample
   read(10,*) word, conf_sample
   read(10,*) word, iseed
   read(10,*) word, number_species
   read(10,*) word, flux				! this variable establishes if the protein is shifted isotropically (0) or if it always move of a single step in a fixed direction (1) (Y positive)
   close(10)



!DR: read commandline first argument
   CALL getarg(1, arg)
!DR: convert argument into integer and assign it to
   read(arg,*) iseed
!DR: print to check consistency of input
   open(10,file="output_dr.txt",status="unknown")
   write(10,'(A25,I20)') 'L: ', L
   write(10,'(A25,F20.10)') 'Temp: ', Temp
   write(10,'(A25,F20.10)') 'Press: ', Press
   write(10,'(A25,I20)') 'number_proteins: ', number_proteins
   write(10,'(A25,I20)') 'init_time: ', init_time
   write(10,'(A25,I20)') 'run_time: ', run_time
   write(10,'(A25,I20)') 'sample: ', sample
   write(10,'(A25,I20)') '(random seed) iseed: ', iseed
   write(10,'(A25,I20)') 'number_species: ', number_species
   write(10,'(A25,I20)') 'flux: ', flux
   close(10)

   if ( ( flux .ne. 0 ) .AND. ( flux .ne. 1 ) )	then
      write(*,*) "Error! The flux flag must be 0 or 1!"
      call exit()
   endif

   open(10,file="water_parameters",status="unknown")
   read(10,*) word,  jj_bulk
   read(10,*) word,  jj_phob
   read(10,*) word,  jj_phil
   read(10,*) word,  jsigma
   read(10,*) word,  jsigma_phob
   read(10,*) word,  vHB_bulk
   read(10,*) word,  vHB_phob
   read(10,*) word,  vHB_phil
   read(10,*) word,  compress_bulk
   read(10,*) word,  compress_phob
   read(10,*) word,  compress_phil
   close(10)
   !DR: print to check consistency of input
   open(10,file="output_dr_2.txt",status="unknown")
   write(10,'(A25,F20.10)') 'jj_bulk: ', jj_bulk
   write(10,'(A25,F20.10)') 'jj_phob: ', jj_phob
   write(10,'(A25,F20.10)') 'jj_phil: ', jj_phil
   write(10,'(A25,F20.10)') 'jsigma: ', jsigma
   write(10,'(A25,F20.10)') 'jsigma_phob: ', jsigma_phob
   write(10,'(A25,F20.10)') 'vHB_bulk: ', vHB_bulk
   write(10,'(A25,F20.10)') 'vHB_phob: ', vHB_phob
   write(10,'(A25,F20.10)') 'vHB_phil: ', vHB_phil
   write(10,'(A25,F20.10)') 'compress_bulk: ', compress_bulk
   write(10,'(A25,F20.10)') 'compress_phob: ', compress_phob
   write(10,'(A25,F20.10)') 'compress_phil: ', compress_phil
   close(10)

   vHB_bulk=vHB_bulk-compress_bulk*Press
   vHB_phob=vHB_phob-compress_phob*Press
   vHB_phil=vHB_phil-compress_phil*Press
   vHB_mix=0.5d0*vHB_phob+0.5d0*vHB_phil
   jj_mix=0.5d0*jj_phob+0.5*jj_phil
   jj_bulk_eff=jj_bulk-Press*vHB_bulk
   jj_phob_eff=jj_phob-Press*vHB_phob
   jj_phil_eff=jj_phil-Press*vHB_phil
   jj_mix_eff=jj_mix-Press*vHB_mix

   num_site=L*L
   half_L=L/2
   dL2=dble(L)*dble(L)
   latt_size=1.0d0 ! lattice spacing
   latt_size2=latt_size**2
   latt=dble(L)
   latt2=latt*latt
   lsample=0.01d0;
   EPS=0.5
   
   

888 format	(I12,4I6)

! time, energy, prot_ ene, vol, HB_bulk, HB_phob, HB_phil, HB_mix, Nim, Nim_phob, bond_phob, bond_phil, cp_ran/cp_max, cp/cp_max

   call dran_ini(iseed) ! Inizialisation of the random number generator
   call initial_setting() ! allocate arrays and Generation of the initial conditions
   call protein_surface()
   call potential_energy()

   open(20,file="Data_folding.dat",status="unknown")
   open(69,file="Protein_Configurations.dat",status="unknown")

   !DR:
   open(130,file="rx.dat",status="unknown")
   open(131,file="ry.dat",status="unknown")
   do i=0,num_site-1
      write(130,'(I4)',advance='no') rx(i)
      write(131,'(I4)',advance='no') ry(i)
   enddo
   close(130)
   close(131)
   open(130,file="status0.dat",status="unknown")
   open(131,file="status1.dat",status="unknown")
   open(132,file="status2.dat",status="unknown")
   open(133,file="status3.dat",status="unknown")
   open(13342,file="latt_size_parameters.dat",status="unknown")
   write(13342,*) DR_HB_beta, DR_HB_l_threshold
   close(13342)
   open(13342,file="latt_size.dat",status="unknown")

!***************************************
!*  EQUILIBRATE AND SAMPLE THE SYSTEM  *
!***************************************

   call protein_surface()
   call potential_energy()
   latt_size=1.0d0 ! lattice spacing
   latt_size2=latt_size**2
   DR_HB_f=1.0d0/( dexp((latt_size-DR_HB_l_threshold)*DR_HB_beta)  +1)
   
   

   do time=1,init_time+run_time
   	

	! if (mod(time,sample) .eq. 0)	then  
	! 	write(*,*) time,latt_size,DR_HB_f
	! endif

      ! FIRST TRY A GLOBAL MOVE OF THE PROTEIN
      choose_prot_move = i_dran(4)
      if (choose_prot_move .eq. 1) then
         trial_pivot = trial_pivot + number_proteins
         call pivot()
      else
         if (choose_prot_move .eq. 2)	then
            trial_shift = trial_shift + number_proteins
            call protein_shift()
         else
            if (choose_prot_move .eq. 3)	then
               trial_crankshaft = trial_crankshaft + number_proteins
               call crankshaft()
            else
               trial_protein_rotation = trial_protein_rotation + number_proteins
               call protein_rotation()
            endif
         endif
      endif

      ! THAN DO A RANDOM SERIES OF LOCAL MOVE OF WATER AND PROTEINS TO EQUILIBRATE THE SYSTEM
      random_steps = i_dran(8*num_site)	! 4*num_sites are the degrees of freedom of the system
      do i=1,random_steps
         particle=i_dran(num_site)-1
         x=rx(particle)
         y=ry(particle)

         st = status(x,y,0)
         if (st .gt. 0) then		! it is a water molecule. The variable st is used as flag and must be redifined for a water molecule
            n=i_dran(4)-1
            st=status(x,y,n) !!!!!!!!!!!!!!!!!!!!!!!!!
            trial_rotation=trial_rotation+1
            call rotation(x,y,n,st)
         else
            trial_jump=trial_jump+1
            call local_protein_jump(x,y,st)
         endif
      enddo ! i

      ! FINALLY UPDATE THE VOLUME
      call potential_energy() !DR added this
      trial_volume = trial_volume + 1
      call volume()
      latt_size=latt_size2**0.5d0
   	DR_HB_f=1.0d0/( dexp((latt_size-DR_HB_l_threshold)*DR_HB_beta)  +1)

      if (time .gt. init_time)	then
         ! SAVE THE DATA
         if (mod(time,sample) .eq. 0)	then
            energy=0.d0
            call protein_energy()
            prot_ene=energy
            call Bonds()
            call potential_energy()
            !energy=EPS*(v12-v6)+energy-jj_bulk*DR_HB_f*HB_bulk-jj_phob*DR_HB_f*HB_phob &
            energy=v12-v6+energy-jj_bulk*DR_HB_f*HB_bulk-jj_phob*DR_HB_f*HB_phob &
            -jj_phil*DR_HB_f*HB_phil-jj_mix*DR_HB_f*HB_mix &
            -jsigma*Nim -jsigma_phob*Nim_phob
            vol=latt2+HB_bulk*vHB_bulk+HB_phob*vHB_phob+HB_phil*vHB_phil+HB_mix*vHB_mix
            778 format 	(I10,3E20.10,8I7)		
            write(20,778) time, energy, prot_ene, vol, HB_bulk, HB_phob, HB_phil, HB_mix, Nim, Nim_phob, bond_phob, bond_phil!, cp_ran/cp_max, cp/cp_max
            write(13342,*) time,latt_size,DR_HB_f
            call flush(20)
            call flush(13342)
         endif
         

         ! SAVE THE CONFIFURATION
         if (mod(time,conf_sample) .eq. 0)	then
            do n=1,number_proteins
               do j=1,protein_length(n)
                  x=protein(n,j,1)
                  y=protein(n,j,2)
                  write(69,888) time, x,y, status(x,y,0), j
               enddo
            enddo

            !DR: saving the status
            do i=0,num_site-1
               write(130,'(I4)',advance='no') status(rx(i),ry(i),0)
               write(131,'(I4)',advance='no') status(rx(i),ry(i),1)
               write(132,'(I4)',advance='no') status(rx(i),ry(i),2)
               write(133,'(I4)',advance='no') status(rx(i),ry(i),3)
            enddo
            write(130,*) ''
            write(131,*) ''
            write(132,*) ''
            write(133,*) ''

            call flush(69)
            call flush(130)
            call flush(131)
            call flush(132)
            call flush(133)


            call protein_surface()
            call Bonds_DR()

         endif
      endif

   enddo


   close(20)
   close(69)
   !DR:
   close(130)
   close(131)
   close(132)
   close(133)
   close(13342)


   open (25, file="Final_Configuration.dat",status="unknown")
   do i=0, num_site-1
      x = rx(i)
      y = ry(i)
      write(25,*) x, y, status(x,y,0), status(x,y,1), status(x,y,2), status(x,y,3)
   enddo
   close(25)

   open (25, file="Final_Protein_Conformation.dat",status="unknown");
   do n=1,number_proteins
      do i=1,protein_length(n)
         write(25,*) protein(n,i,1), protein(n,i,2)
      enddo
   enddo
   close(25)

   open(25,file="Acceptance.dat",status="unknown")
   ratio=acc_rotation
   ratio=ratio/trial_rotation
   write(25,*) "Acceptance of Rotation move = ", ratio
   ratio=acc_pivot
   ratio=ratio/trial_pivot
   write(25,*) "Acceptance of Pivot move = ", ratio
   ratio=acc_jump
   ratio=ratio/trial_jump
   write(25,*) "Acceptance of Jump move = ", ratio
   ratio=acc_crankshaft
   ratio=ratio/trial_crankshaft
   write(25,*) "Acceptance of Crankshaft move = ", ratio
   ratio=acc_shift
   ratio=ratio/trial_shift
   write(25,*) "Acceptance of Shift move = ", ratio
   ratio=acc_protein_rotation
   ratio=ratio/trial_protein_rotation
   write(25,*) "Acceptance of Protein Rotation move = ", ratio
   ratio=acc_volume
   ratio=ratio/trial_volume
   write(25,*) "Acceptance of Volume move = ", ratio

   call flush(25)
   close(25)

end program protein_design_folding


! ******** SUBROUTINES ***********

! *******************
! * INITIAL SETTING *
! *******************

SUBROUTINE initial_setting()

   use common_variables_folding

   implicit none
   INTEGER :: n, m, i, j, p, k, sigma, flag, x, y, dx, dy, r2, initial_shift, mmm
   INTEGER(KIND=8) :: i_dran
   REAL(KIND=8) :: fl, dist
   LOGICAL :: file1_exists, file2_exists, file3_exists


   allocate (protein_length(1:number_proteins))
   allocate (protein_species(1:number_proteins))
   allocate (protein_move(1:number_proteins))

   allocate (native_contacts_per_protein(1:number_proteins))
   allocate (non_native_contacts_per_protein(1:number_proteins))
   allocate (inter_contacts_per_protein(1:number_proteins))
   allocate (total_contacts_per_protein(1:number_proteins))

   allocate (native_contacts_per_species(1:number_species))
   allocate (non_native_contacts_per_species(1:number_species))
   allocate (inter_contacts_per_species(1:number_species))
   allocate (total_contacts_per_species(1:number_species))

   allocate (cp_max_per_protein(1:number_proteins))

   max_length = 0
   tot_prot_length=0
   open (66, file="protein_species.dat",status="unknown");		!	THIS FILE CONTAINS THE INDEX IDENTIFYING FOR THE PROTEIN SPECIES. THIS NUMBER IS AN INTEGER AND SHOULD START FROM 1
   do i=1,number_proteins
      read(66,*) n
      protein_species(i) = n
   enddo
   close(66)

   open (66, file="protein_moves.dat",status="unknown");		!	THIS FILE CONTAINS A FLAG THAT ESTABLISHES IF THE PROTEIN WILL BE MOVED OR NOT. THIS FLAG CAN BE 0 (PROTEIN FIXED) OR 1 (MOBILE PROTEIN)
   do i=1,number_proteins
      read(66,*) n
      protein_move(i) = n
      if ( ( n .ne. 0 ) .AND. ( n .ne. 1 ) )	then
         write(*,*) "Error! The move flag must be 0 or 1!"
         call exit()
      endif
   enddo
   close(66)

   open (66, file="protein_lengths.dat",status="unknown");	!	THIS FILE CONTAINS THE NUMBER OF AMINO ACIDS OF EACH PROTEINS
   do i=1,number_proteins
      read(66,*) n

      if ((n .ge. L) .AND. (protein_move(i) .eq. 1 ))	then
         write(*,*) "Error! The length of a mobile protein cannot exceed the box size!"
         write(*,*) "Check the protein ",i
         call exit()
      endif
      if ((n .gt. L) .AND. (protein_move(i) .eq. 0 ))	then
         write(*,*) "Error! The length of a fixed protein cannot exceed the box size!"
         write(*,*) "Check the protein ",i
         call exit()
      endif

      protein_length(i) = n
      tot_prot_length=tot_prot_length+n
      if ( n .gt. max_length) then
         max_length = n
      endif
   enddo
   close(66)

   mmm = max_length*number_proteins
   allocate (pbc(-L+1:(2*L)))
   allocate (rx(0:(L*L-1)))
   allocate (ry(0:(L*L-1)))
   allocate (neigh(1:number_proteins , 1:max_length , 1:max_length))
   allocate (neigh_run(1:number_proteins , 1:max_length , 1:max_length))
   allocate (protein(1:number_proteins , 1:max_length , 1:2))
   allocate (protein_old(1:max_length , 1:2))
   allocate (protein_aux(1:max_length , 1:2))
   allocate (surface(1:L,1:L))
   allocate (surface_old(1:L,1:L))
   allocate (status(1:L,1:L,0:3))
   allocate (status_old(1:L,1:L,0:3))
   allocate (protein_state(1:number_proteins , 1:max_length))
   allocate (px(1:number_proteins , 1:max_length))
   allocate (py(1:number_proteins , 1:max_length))


   acc_rotation=0; trial_rotation=0; acc_pivot=0; trial_pivot=0; acc_volume=0; trial_volume=0;
   acc_crankshaft=0; trial_crankshaft=0; acc_jump=0; trial_jump=0;

   do n=1,L !periodic boundary conditions
      pbc(n)=n
      pbc(L+n)=n
      pbc(1-n)=L+1-n
   enddo

   do m=1,L ! coordinates of the lattice
      do p=1,L
         rx((m-1)*L+p-1)=m
         ry((m-1)*L+p-1)=p
      enddo
   enddo

   do i=0,num_site-1
      m=rx(i)
      p=ry(i)
      do k=0,3
         status(m,p,k)=i_dran(q) ! random initial Potts variables
      enddo
   enddo

   cx = (/0,1,0,-1/)
   cy = (/-1,0,1,0/)
   neig_arm = (/2,3,0,1/)

   round_part=0
   x=-rings
   do m=1,2*rings+1
      y=-rings
      do n=1,2*rings+1
         dist= sqrt( dble(x)**2+dble(y)**2)
         if ((dist .le. rings) .AND. ((x**2+y**2) .ne. 0))	round_part=round_part+1
         y=y+1
      enddo
      x=x+1
   enddo

   allocate (displacement(1:round_part,1:2))

   x=-rings
   i=1
   do m=1,2*rings+1 ! A square with side equal to "2*ring+1" contains a circle of radius "ring"
      y=-rings
      do n=1,2*rings+1
         dist= sqrt( dble(x)**2+dble(y)**2)
         if ((dist .le. rings) .AND. ((x**2+y**2) .ne. 0))		then
            displacement(i,1)=x
            displacement(i,2)=y
            i=i+1
         endif
         y=y+1
      enddo
      x=x+1
   enddo

   initial_shift = L/number_proteins

   inquire( file="Final_Configuration.dat", exist=file1_exists )
   inquire( file="Final_Protein_Conformation.dat", exist=file2_exists )

   if ( (file1_exists) .AND. (file2_exists) )	then
      write(*,*) "Found the files Final_Configuration.dat and Final_Protein_Conformation.dat of the previous simulation."
      write(*,*) "The new run will continue from such files."
      open (87, file="Final_Configuration.dat",status="unknown")
      do i=0, num_site-1
         read(87,*) x, y, status(x,y,0), status(x,y,1), status(x,y,2), status(x,y,3)
      enddo
      close(87)

      open (87, file="Final_Protein_Conformation.dat",status="unknown");
      do n=1,number_proteins
         do i=1,protein_length(n)
            read(87,*) protein(n,i,1), protein(n,i,2)
         enddo
      enddo
      close(87)
   else
      inquire( file="Starting_Protein_Conformation.dat", exist=file3_exists )
      if (file3_exists)	then
         write(*,*) "Found the file Starting_Protein_Conformation.dat."
         write(*,*) "The initial configuration will be random but the proteins will be placed according to such a file."
         open (87, file="Starting_Protein_Conformation.dat",status="unknown")
         do n=1,number_proteins
            do i=1,protein_length(n)
               read(87,*) protein(n,i,1), protein(n,i,2)
               x=protein(n,i,1)
               y=protein(n,i,2)
               do j=0,3
                  status(x,y,j)=-1
               enddo
            enddo
         enddo
      else
         write(*,*) "No initial conformation found. The initial configuration will be random."
         write(*,*) "The protein will be placed stretched and equally shifted one from the other."
         do n=1,number_proteins
            do i=1,protein_length(n)
               if (MOD(n,2) .eq. 0)	then
                  protein(n,i,1)=i
               else
                  protein(n,i,1)=L+1-i
               endif
               protein(n,i,2)=(n-1)*initial_shift + 1
            enddo
         enddo
      endif
   endif

   open(13,file="input_sequences.dat",status="unknown")
   do n=1,number_proteins
      do i=1,protein_length(n)
         read(13,*) flag, sigma
         if ((flag .ne. 1) .AND. (i .eq. 1))	then
            write(*,*) "Error in the input sequences!!!!"
            call exit()
         endif
         do k=0,3
            status(protein(n,i,1),protein(n,i,2),k)=sigma
         enddo
         protein_state(n,i)=sigma
      enddo
   enddo
   close(13)

   open (87, file="target_structures.dat",status="unknown")
   do n=1,number_proteins
      do i=1,protein_length(n)
         read(87,*) x, y
         px(n,i)=x
         py(n,i)=y
      enddo
   enddo
   close(87)

   neigh=0
   cp_max=0
   cp_max_per_protein = 0
   do n=1,number_proteins
      do i=1,protein_length(n)-1
         do j=i+1,protein_length(n)
            dx=abs(px(n,i)-px(n,j))
            dy=abs(py(n,i)-py(n,j))
            if (dx .ge. half_L)	dx=L-dx
            if (dy .ge. half_L)	dy=L-dy
            r2=dx*dx+dy*dy
            if ((r2 .eq. 1) .AND. (j .ne. i+1))	then
               neigh(n,i,j)=1
               neigh(n,j,i)=1
               cp_max_per_protein(n) = cp_max_per_protein(n) + 1
               cp_max=cp_max+1
            endif
         enddo
      enddo
   enddo

   open (87, file="aapot_water.dat",status="unknown");
   matrix=0
   do i=1,num_max_aminoacids	! 	ADD THE RESIDUE-WATER INTERACTION VALUES
      read(87,*) fl
      do j=1,6
         matrix(-i,j)=fl
         matrix(j,-i)=fl
      enddo
   enddo

   do j=1,num_max_aminoacids
      do i=j,num_max_aminoacids
         read(87,*) fl
         matrix(-i,-j)=fl
         matrix(-j,-i)=fl
      enddo
   enddo
   close(87)

   open (87, file="Initial_Configuration.dat",status="unknown")
   do i=0, num_site-1
      x = rx(i)
      y = ry(i)
      write(87,*) x, y, status(x,y,0), status(x,y,1), status(x,y,2), status(x,y,3)
   enddo
   close(87)

   open (87, file="Initial_Protein_Conformation.dat",status="unknown");
   do n=1,number_proteins
      do i=1,protein_length(n)
         write(87,*) protein(n,i,1), protein(n,i,2)
      enddo
   enddo
   close(87)

end SUBROUTINE




! ***********************
! * LJ  N.N. ENERGY     *
! ***********************
!  LJ energy of the nearest neighbours

SUBROUTINE potential_energy ()

   use common_variables_folding

   implicit none
   INTEGER :: i, n, x, y, dx, dy, vx, vy
   REAL(KIND=8) :: aux, aux3, r2

   v12=0.d0 ! initial lj energy, V(r)=EPS*(1/r**12-1/r**6), of all system
   v6=0.d0

   do n=0,num_site-1
      x=rx(n)
      y=ry(n)
      if (status(x,y,0) .gt. 0) then
         do i=1,round_part
            vx=pbc(x+displacement(i,1))
            vy=pbc(y+displacement(i,2))
            if (status(vx,vy,0) .gt. 0) then
               dx=x-vx ! horizontal distance
               dy=y-vy ! vertical distance
               if (iabs(dx) .gt. half_L) dx=dx-sign(L,dx)
               if (iabs(dy) .gt. half_L) dy=dy-sign(L,dy)
               r2=dble(dx**2+dy**2)*latt_size2
               aux=1.d0/r2
               aux3=aux*aux*aux
               v6=v6+aux3
               v12=v12+aux3*aux3
            endif
         enddo
      endif
   enddo

   v6=v6/2.d0
   v12=v12/2.d0

end SUBROUTINE


! **********************
! *    PROTEIN JUMP    *
! **********************

SUBROUTINE local_protein_jump(x,y,st)

   use common_variables_folding

   implicit none
   INTEGER :: x, y, st ! dummy variables
   INTEGER :: newx, newy, amin=-1
   INTEGER :: prax, prcx, pray, prcy, next_x, next_y, next_st, n, prot=-1
   INTEGER :: distance1, distance2, delta1, delta2, choice
   INTEGER :: deltax1, deltax2, deltay1, deltay2
   INTEGER(KIND=8) :: i_dran
   REAL(KIND=8) :: delta_h, metropolis, z, dran_u
   REAL(KIND=8) :: energy_old, energy_new, v6old, v6new, v12old, v12new, deltav12, deltav6, delta_lj
   INTEGER :: flag_move, flag_change
   INTEGER :: u, i

   flag_move=0

   do n=1,number_proteins
      do i=1,protein_length(n)
         if ((protein(n,i,1) .eq. x) .AND. (protein(n,i,2) .eq. y))	then
            amin=i
            prot=n
         endif
      enddo
   enddo

   if ((amin .eq. -1) .OR. (prot .eq. -1))	then
      write(*,*)	"Errore, Local_protein_jump, amin=-1 or prot=-1"
      write(*,*) "x=",x,"y=",y,"st=",st
      call exit()
   endif

   if ( protein_move(prot) .eq. 1 )	then	!	Check if the residue belongs to a fixed or mobile protein
      if (amin .eq. 1) then	! Its the first amino acid

         choice=i_dran(4)
         select case (choice)
          case (1)
            newx=pbc(protein(prot,2,1)+1)
            newy=protein(prot,2,2)
            if (status(newx,newy,0) .gt. 0) then
               if (newx .ne. x) then ! check if it's not the initial position
                  flag_move=1
               endif
            endif
          case (2)
            newx=pbc(protein(prot,2,1)-1)
            newy=protein(prot,2,2)
            if (status(newx,newy,0) .gt. 0) then
               if (newx .ne. x) then ! check if it's not the initial position
                  flag_move=1
               endif
            endif
          case (3)
            newx=protein(prot,2,1)
            newy=pbc(protein(prot,2,2)+1)
            if (status(newx,newy,0) .gt. 0) then
               if (newy .ne. y) then ! check if it's not the initial position
                  flag_move=1
               endif
            endif
          case (4)
            newx=protein(prot,2,1)
            newy=pbc(protein(prot,2,2)-1)
            if (status(newx,newy,0) .gt. 0) then
               if (newy .ne. y) then ! check if it's not the initial position
                  flag_move=1
               endif
            endif
          case default
            write(*,*) "ERRORE NELLA ROUTINE protein_jump"
         end select
      else
         if (amin .eq. protein_length(prot)) then ! Its the last amino acid
            choice=i_dran(4)
            select case (choice)
             case (1)
               newx=pbc(protein(prot,protein_length(prot)-1,1)+1)
               newy=protein(prot,protein_length(prot)-1,2)
               if (status(newx,newy,0) .gt. 0) then
                  if (newx .ne. x) then ! check if it's not the initial position
                     flag_move=1
                  endif
               endif
             case (2)
               newx=pbc(protein(prot,protein_length(prot)-1,1)-1)
               newy=protein(prot,protein_length(prot)-1,2)
               if (status(newx,newy,0) .gt. 0) then
                  if (newx .ne. x) then ! check if it's not the initial position
                     flag_move=1
                  endif
               endif
             case (3)
               newx=protein(prot,protein_length(prot)-1,1)
               newy=pbc(protein(prot,protein_length(prot)-1,2)+1)
               if (status(newx,newy,0) .gt. 0) then
                  if (newy .ne. y) then ! check if it's not the initial position
                     flag_move=1
                  endif
               endif
             case (4)
               newx=protein(prot,protein_length(prot)-1,1)
               newy=pbc(protein(prot,protein_length(prot)-1,2)-1)
               if (status(newx,newy,0) .gt. 0) then
                  if (newy .ne. y) then ! check if it's not the initial position
                     flag_move=1
                  endif
               endif
             case default
               write(*,*) "ERRORE NELLA ROUTINE PROTEIN_jump"
            end select
         else ! Its an internal amino acid
            prax=protein(prot,amin-1,1)! x coordinate of the i-1 amino acid
            prcx=protein(prot,amin+1,1)! x coordinate of the i+1 amino acid
            pray=protein(prot,amin-1,2)! y coordinate of the i-1 amino acid
            prcy=protein(prot,amin+1,2)! y coordinate of the i+1 amino acid

            if (pray .ne. prcy) then 	! Check if they are in an orizzontal line
               if (prax .ne. prcx) then 	! Check if they are on a vertical line
                  choice=i_dran(4)
                  select case (choice)
                   case (1) ! try to move the amino acid on the right
                     newx=pbc(protein(prot,amin,1)+1)
                     newy=pbc(protein(prot,amin,2)+1)
                     if (status(newx,newy,0) .gt. 0) then
                        delta1=newx-protein(prot,amin-1,1)
                        deltax1=min0((L-iabs(delta1)),iabs(delta1))
                        delta1=newy-protein(prot,amin-1,2)
                        deltay1=min0((L-iabs(delta1)),iabs(delta1))
                        distance1=deltax1+deltay1
                        delta2=newx-protein(prot,amin+1,1)
                        deltax2=min0((L-iabs(delta2)),iabs(delta2))
                        delta2=newy-protein(prot,amin+1,2)
                        deltay2=min0((L-iabs(delta2)),iabs(delta2))
                        distance2=deltax2+deltay2
                        if (distance1 .eq. 1) then
                           if (distance2 .eq. 1) then ! possible jump
                              flag_move=1
                           endif
                        endif
                     endif
                   case (2)	! try to move the amino acid on the left
                     newx=pbc(protein(prot,amin,1)-1)
                     newy=pbc(protein(prot,amin,2)+1)
                     if (status(newx,newy,0) .gt. 0) then
                        delta1=newx-protein(prot,amin-1,1)
                        deltax1=min0((L-iabs(delta1)),iabs(delta1))
                        delta1=newy-protein(prot,amin-1,2)
                        deltay1=min0((L-iabs(delta1)),iabs(delta1))
                        distance1=deltax1+deltay1
                        delta2=newx-protein(prot,amin+1,1)
                        deltax2=min0((L-iabs(delta2)),iabs(delta2))
                        delta2=newy-protein(prot,amin+1,2)
                        deltay2=min0((L-iabs(delta2)),iabs(delta2))
                        distance2=deltax2+deltay2
                        if (distance1 .eq. 1) then
                           if (distance2 .eq. 1) then ! possible jump
                              flag_move=1
                           endif
                        endif
                     endif
                   case (3) 	! try to move the amino acid up
                     newx=pbc(protein(prot,amin,1)+1)
                     newy=pbc(protein(prot,amin,2)-1)
                     if (status(newx,newy,0) .gt. 0) then
                        delta1=newx-protein(prot,amin-1,1)
                        deltax1=min0((L-iabs(delta1)),iabs(delta1))
                        delta1=newy-protein(prot,amin-1,2)
                        deltay1=min0((L-iabs(delta1)),iabs(delta1))
                        distance1=deltax1+deltay1
                        delta2=newx-protein(prot,amin+1,1)
                        deltax2=min0((L-iabs(delta2)),iabs(delta2))
                        delta2=newy-protein(prot,amin+1,2)
                        deltay2=min0((L-iabs(delta2)),iabs(delta2))
                        distance2=deltax2+deltay2
                        if (distance1 .eq. 1) then
                           if (distance2 .eq. 1) then ! possible jump
                              flag_move=1
                           endif
                        endif
                     endif
                   case (4)	! try to move the amino acid down
                     newx=pbc(protein(prot,amin,1)-1)
                     newy=pbc(protein(prot,amin,2)-1)
                     if (status(newx,newy,0) .gt. 0) then
                        delta1=newx-protein(prot,amin-1,1)
                        deltax1=min0((L-iabs(delta1)),iabs(delta1))
                        delta1=newy-protein(prot,amin-1,2)
                        deltay1=min0((L-iabs(delta1)),iabs(delta1))
                        distance1=deltax1+deltay1
                        delta2=newx-protein(prot,amin+1,1)
                        deltax2=min0((L-iabs(delta2)),iabs(delta2))
                        delta2=newy-protein(prot,amin+1,2)
                        deltay2=min0((L-iabs(delta2)),iabs(delta2))
                        distance2=deltax2+deltay2
                        if (distance1 .eq. 1) then
                           if (distance2 .eq. 1) then ! possible jump
                              flag_move=1
                           endif
                        endif
                     endif
                   case default
                     write(*,*) "ERRORE NELLA ROUTINE PROTEIN_jump"
                  end select
               endif !end of horizontal check
            endif !end of vertical check
         endif !end of terminal acid check
      endif !end of first acid check

      if (flag_move .eq. 1) then
         energy_old=0.d0	! the energy of the adjacent residues will be counted for energy_old and energy_new
         next_x=pbc(x+1);		next_y=y;			next_st=status(next_x,next_y,0);
         energy_old=energy_old+matrix(st,next_st)
         next_x=pbc(x-1);		next_y=y;			next_st=status(next_x,next_y,0);
         energy_old=energy_old+matrix(st,next_st)
         next_x=x;				next_y=pbc(y+1);	next_st=status(next_x,next_y,0);
         energy_old=energy_old+matrix(st,next_st)
         next_x=x;				next_y=pbc(y-1);	next_st=status(next_x,next_y,0);
         energy_old=energy_old+matrix(st,next_st)

         call molecule_potential_energy(newx,newy,v12old,v6old) ! Evaluate the OLD LJ energy of the water molecule in position (newx, newy)

         call local_Bonds(x,y)	! calculate the water bonds of the sorrounding molecules
         energy_old=energy_old-jj_bulk_eff*DR_HB_f*HB_bulk_local-jj_phob_eff*DR_HB_f*HB_phob_local & 
         -jj_phil_eff*DR_HB_f*HB_phil_local &
         -jj_mix_eff*DR_HB_f*HB_mix_local -jsigma*Nim_local -jsigma_phob*Nim_phob_local


         !write(*,*) "BEFORE LOCAL JUMP"
         !do n=1,number_proteins
         !	do i=1,protein_length(n)
         !		write(*,*) protein(n,i,1), protein(n,i,2), n
         !	enddo
         !enddo


         do u=0,3
            flag_change=status(x,y,u)	! change of status ...
            status(x,y,u)=status(newx,newy,u)
            status(newx,newy,u)=flag_change
         enddo
         protein(prot,amin,1)=newx						! ... protein ...
         protein(prot,amin,2)=newy

         call protein_surface()

         energy_new=0.d0	! the energy of the adjacent residues will be counted for energy_new and energy_new

         next_x=pbc(newx+1);		next_y=newy;			next_st=status(next_x,next_y,0);
         energy_new=energy_new+matrix(st,next_st)
         next_x=pbc(newx-1);		next_y=newy;			next_st=status(next_x,next_y,0);
         energy_new=energy_new+matrix(st,next_st)
         next_x=newx;			next_y=pbc(newy+1);		next_st=status(next_x,next_y,0);
         energy_new=energy_new+matrix(st,next_st)
         next_x=newx;			next_y=pbc(newy-1);		next_st=status(next_x,next_y,0);
         energy_new=energy_new+matrix(st,next_st)

         call local_Bonds(x,y)	! calculate the water bonds of the sorrounding molecules
         energy_new=energy_new-jj_bulk_eff*DR_HB_f*HB_bulk_local-jj_phob_eff*DR_HB_f*HB_phob_local & 
         -jj_phil_eff*DR_HB_f*DR_HB_f*HB_phil_local &
            -jj_mix_eff*DR_HB_f*HB_mix_local -jsigma*Nim_local -jsigma_phob*Nim_phob_local

         call molecule_potential_energy(x,y,v12new,v6new) ! Evaluate the NEW LJ energy of the water molecule in position (x,y)
         deltav12=dble(v12new-v12old)
         deltav6=dble(v6new-v6old)
         delta_lj=EPS*(deltav12-deltav6)

         delta_h=delta_lj+energy_new-energy_old
         acc_jump=acc_jump+1
         if (delta_h .gt. 0.d0) then ! metropolis
            z=dran_u()
            metropolis=dexp(-delta_h/Temp)
            if (z .gt. metropolis) then ! secondary
               protein(prot,amin,1)=x	! move not accepted, re-change the status
               protein(prot,amin,2)=y
               do u=0,3
                  flag_change=status(x,y,u)
                  status(x,y,u)=status(newx,newy,u)
                  status(newx,newy,u)=flag_change
               enddo
               call protein_surface()
               acc_jump=acc_jump-1
            endif ! end secondary metropolis
         endif ! end main metropolis
      endif !

      !write(*,*) "AFTER LOCAL JUMP"
      !	do n=1,number_proteins
      !		do i=1,protein_length(n)
      !			write(*,*) protein(n,i,1), protein(n,i,2), n
      !		enddo
      !	enddo

      do n=1,number_proteins
         do i=1,protein_length(n)-1
            x=protein(n,i+1,1)-protein(n,i,1)
            y=protein(n,i+1,2)-protein(n,i,2)
            if (iabs(x) .gt. 1) x=x-sign(L,x)
            if (iabs(y) .gt. 1) y=y-sign(L,y)
            if ((iabs(x) .gt. 1) .OR. (iabs(y) .gt. 1))	then
               write (*,*) "After local_jump, Broken protein , n = ",n
               call exit()
            endif
         enddo
      enddo
   endif

end SUBROUTINE


! ***************************
! * LJ  N.N. ENERGY CHECK   *
! ***************************

SUBROUTINE molecule_potential_energy(x_aux,y_aux,v12_aux,v6_aux)

   use common_variables_folding

   implicit none
   INTEGER:: x_aux, y_aux		! dummy variables
   REAL(KIND=8) :: v12_aux, v6_aux	! dummy variables
   REAL(KIND=8) :: r2, aux_ext, aux3_ext
   INTEGER:: i, dx_aux, dy_aux, vx_aux, vy_aux

   v12_aux=0.d0
   v6_aux=0.d0
   ! check the LJ energy of a single molecule

   do i=1,round_part
      vx_aux=pbc(x_aux+displacement(i,1)) ! coordinates of the site with x and y distances less/equal than "ring"
      vy_aux=pbc(y_aux+displacement(i,2))
      if (status(vx_aux,vy_aux,0) .gt. 0) then   !checks if the other molecule is a water molecule
         dx_aux=x_aux-vx_aux ! horizontal distance
         dy_aux=y_aux-vy_aux ! vertical distance
         if (iabs(dx_aux) .gt. half_L) dx_aux=dx_aux-sign(L,dx_aux)
         if (iabs(dy_aux) .gt. half_L) dy_aux=dy_aux-sign(L,dy_aux)
         r2=dble(dx_aux**2+dy_aux**2)*latt_size2
         aux_ext=1.d0/r2
         aux3_ext=aux_ext*aux_ext*aux_ext
         v6_aux=v6_aux+aux3_ext
         v12_aux=v12_aux+aux3_ext*aux3_ext
      endif
   enddo

end SUBROUTINE


! *************************
! *    PROTEIN SURFACE    *
! *************************

SUBROUTINE protein_surface()	! the surface elements are equal to 1, 2, or 3

   use common_variables_folding

   implicit none
   INTEGER :: m, x, y, st, n

   surface=0
   do n=1,number_proteins
      do m=1,protein_length(n) ! find the surface of the protein
         x=protein(n,m,1)
         y=protein(n,m,2)
         st=status(x,y,0)
         if (matrix(st,1) .ge. 0) then	! hydrophobic element - I check the sign of the residue-water interaction
            if (status(pbc(x+1),y,0) .gt. 0) then	! check if the new site is a water molecule of a protein monomer
               if (surface(pbc(x+1),y) .eq. 0)	then
                  surface(pbc(x+1),y)=1;		! hydrophobic interface
               elseif (surface(pbc(x+1),y) .eq. 2)	then
                  surface(pbc(x+1),y)=3;	! mix interface
               endif
            endif
            if (status(pbc(x-1),y,0) .gt. 0) then
               if (surface(pbc(x-1),y) .eq. 0)	then
                  surface(pbc(x-1),y)=1;
               elseif (surface(pbc(x-1),y) .eq. 2)	then
                  surface(pbc(x-1),y)=3;
               endif
            endif
            if (status(x,pbc(y+1),0) .gt. 0) then
               if (surface(x,pbc(y+1)) .eq. 0)	then
                  surface(x,pbc(y+1))=1;
               elseif (surface(x,pbc(y+1)) .eq. 2)	then
                  surface(x,pbc(y+1))=3;
               endif
            endif
            if (status(x,pbc(y-1),0) .gt. 0) then
               if (surface(x,pbc(y-1)) .eq. 0)	then
                  surface(x,pbc(y-1))=1;
               elseif (surface(x,pbc(y-1)) .eq. 2)	then
                  surface(x,pbc(y-1))=3;
               endif
            endif
         elseif (matrix(st,1) .lt. 0) then					!hydropilic element
            if (status(pbc(x+1),y,0) .gt. 0) then
               if (surface(pbc(x+1),y) .eq. 0)	then
                  surface(pbc(x+1),y)=2;			!hydrophilic interface
               elseif (surface(pbc(x+1),y) .eq. 1)	then
                  surface(pbc(x+1),y)=3;
               endif
            endif
            if (status(pbc(x-1),y,0) .gt. 0) then
               if (surface(pbc(x-1),y) .eq. 0)	then
                  surface(pbc(x-1),y)=2;
               elseif (surface(pbc(x-1),y) .eq. 1)	then
                  surface(pbc(x-1),y)=3;
               endif
            endif
            if (status(x,pbc(y+1),0) .gt. 0) then
               if (surface(x,pbc(y+1)) .eq. 0)	then
                  surface(x,pbc(y+1))=2;
               elseif (surface(x,pbc(y+1)) .eq. 1)	then
                  surface(x,pbc(y+1))=3;
               endif
            endif
            if (status(x,pbc(y-1),0) .gt. 0) then
               if (surface(x,pbc(y-1)) .eq. 0)	then
                  surface(x,pbc(y-1))=2;
               elseif (surface(x,pbc(y-1)) .eq. 1)	then
                  surface(x,pbc(y-1))=3;
               endif
            endif
         endif
      enddo
   enddo !surface found

end SUBROUTINE


! ******************
! *  LOCAL BONDS   *
! ******************
! Calculate the local number of HB and the surface bounds!

SUBROUTINE local_Bonds(x0,y0)  ! IMPORTANT - YOU EVER MUST CALL "protein_surface" SUBROUTINE BEFORE TO CALL THIS SUBROUTINE

   use common_variables_folding

   implicit none
   INTEGER :: x0, y0 ! dummy variables
   INTEGER :: st, m, n, i, x, y, k, neig_st, su, neig_su, neig_x, neig_y

   HB_bulk_local=0;	HB_phil_local=0;	HB_phob_local=0;	HB_mix_local=0;
   Nim_local=0;		Nim_phob_local=0;
   do m=-2,2
      do n=-2,2
         x=pbc(x0+m)
         y=pbc(y0+n)
         su=surface(x,y)
         do i=0,1
            st=status(x,y,i)
            neig_x=pbc(x+cx(i))
            neig_y=pbc(y+cy(i))
            neig_st=status(neig_x,neig_y,neig_arm(i))
            neig_su=surface(neig_x,neig_y)
            if ((st .eq. neig_st) .AND. (st .gt. 0))	then	! both are water molecules
               if ((su .eq. 0) .OR. (neig_su .eq. 0))	then
                  HB_bulk_local=HB_bulk_local+1
               else
                  if ((su .eq. 1) .AND. (neig_su .eq. 1))	then
                     HB_phob_local=HB_phob_local+1
                  else
                     if ((su .eq. 2) .AND. (neig_su .eq. 2))	then
                        HB_phil_local=HB_phil_local+1
                     else
                        HB_mix=HB_mix+1
                     endif
                  endif
               endif
            endif
         enddo

         st=status(x,y,0)
         if (st .gt. 0)	then	! it is a water molecule
            if ((su .ne. 1) .AND. (su .ne. 3))	then
               do i=0,2
                  do k=i+1,3
                     if (status(x,y,i) .eq. status(x,y,k))	Nim_local=Nim_local+1
                  enddo
               enddo
            else
               do i=0,2
                  do k=i+1,3
                     if (status(x,y,i) .eq. status(x,y,k))	Nim_phob_local=Nim_phob_local+1
                  enddo
               enddo
            endif
         endif
      enddo
   enddo

end SUBROUTINE


! *********
! * PIVOT *
! *********

SUBROUTINE pivot()

   use common_variables_folding

   implicit none
   INTEGER :: dir, node, x, y, x1, y1, i, k, dx, dy, rot, n, nn, useless_flag, j
   INTEGER :: xx, yy, pp, st, particle, rand_step
   INTEGER(KIND=8) :: i_dran
   REAL (KIND=8) :: energy_new, energy_old, delta_h, z, metropolis, dran_u


   do nn=1,number_proteins

      n=i_dran(number_proteins)

      if ( protein_move(n) .eq. 1 )	then	! Check if the protein is fixed or mobile

         ! directions: "0" down, "1" right, "2" up, "3" left
         node=i_dran(protein_length(n)-3)+1	! choose a random monomer excluding the two on extremes

         do j=1,protein_length(n)
            protein_aux(j,1)=protein(n,j,1)
            protein_aux(j,2)=protein(n,j,2)
         enddo
         rot=i_dran(3)	! choose if you want to turn the protein of 90° (1), 180° (2) or 270° (3)

         do i=node,protein_length(n)-1
            dx=protein(n,i+1,1)-protein(n,i,1)
            dy=protein(n,i+1,2)-protein(n,i,2)
            if (iabs(dx) .gt. 1) dx=dx-sign(L,dx)
            if (iabs(dy) .gt. 1) dy=dy-sign(L,dy)
            dir=(dx+2*dy)*rot*rot
            select case (dir)
             case (1)	! dx=1 dy=0 rot=1 ---> dx_new=0 dy_new=1
               protein_aux(i+1,1)=protein_aux(i,1)
               protein_aux(i+1,2)=pbc(protein_aux(i,2)+1)
             case (4)	! dx=1 dy=0 rot=2 ---> dx_new=-1 dy_new=0
               protein_aux(i+1,1)=pbc(protein_aux(i,1)-1)
               protein_aux(i+1,2)=protein_aux(i,2)
             case (9)	! dx=1 dy=0 rot=3 ---> dx_new=0 dy_new=-1
               protein_aux(i+1,1)=protein_aux(i,1)
               protein_aux(i+1,2)=pbc(protein_aux(i,2)-1)

             case (2)	! dx=0 dy=1 rot=1 ---> dx_new=-1 dy_new=0
               protein_aux(i+1,1)=pbc(protein_aux(i,1)-1)
               protein_aux(i+1,2)=protein_aux(i,2)
             case (8)	! dx=0 dy=1 rot=2 ---> dx_new=0 dy_new=-1
               protein_aux(i+1,1)=protein_aux(i,1)
               protein_aux(i+1,2)=pbc(protein_aux(i,2)-1)
             case (18)	! dx=0 dy=1 rot=3 ---> dx_new=1 dy_new=0
               protein_aux(i+1,1)=pbc(protein_aux(i,1)+1)
               protein_aux(i+1,2)=protein_aux(i,2)

             case (-1)	! dx=-1 dy=0 rot=1 ---> dx_new=0 dy_new=-1
               protein_aux(i+1,1)=protein_aux(i,1)
               protein_aux(i+1,2)=pbc(protein_aux(i,2)-1)
             case (-4)	! dx=-1 dy=0 rot=2 ---> dx_new=1 dy_new=0
               protein_aux(i+1,1)=pbc(protein_aux(i,1)+1)
               protein_aux(i+1,2)=protein_aux(i,2)
             case (-9)	! dx=-1 dy=0 rot=3 ---> dx_new=0 dy_new=1
               protein_aux(i+1,1)=protein_aux(i,1)
               protein_aux(i+1,2)=pbc(protein_aux(i,2)+1)

             case (-2)	! dx=0 dy=-1 rot=1 ---> dx_new=1 dy_new=0
               protein_aux(i+1,1)=pbc(protein_aux(i,1)+1)
               protein_aux(i+1,2)=protein_aux(i,2)

             case (-8)	! dx=0 dy=-1 rot=2 ---> dx_new=0 dy_new=1
               protein_aux(i+1,1)=protein_aux(i,1)
               protein_aux(i+1,2)=pbc(protein_aux(i,2)+1)

             case (-18)	! dx=0 dy=-1 rot=3 ---> dx_new=-1 dy_new=0
               protein_aux(i+1,1)=pbc(protein_aux(i,1)-1)
               protein_aux(i+1,2)=protein_aux(i,2)
             case default
               write(*,*) "ERRORE NELLA ROTAZIONE CON PIVOT"
               call exit()
            end select
         enddo


         do i=1,protein_length(n)-1	! check of the overlap with itself to establish if the rotation is possible
            do k=i+1,protein_length(n)
               dx=protein_aux(i,1)-protein_aux(k,1)
               dy=protein_aux(i,2)-protein_aux(k,2)
               if ((dx .eq. 0) .AND. (dy .eq. 0))	goto 87
            enddo
         enddo

         do i=node,protein_length(n)	! check of the overlap with other proteins to establish if the rotation is possible
            x = protein_aux(i,1)
            y = protein_aux(i,2)
            if (status(x,y,0) .lt. 0)	then	! possible overlap
               do j=1,number_proteins		! look all the proteins
                  if (j .ne. n)	then	! exclude the moving protein that has been already checked
                     do k=1,protein_length(j)
                        x1 = protein(j,k,1)
                        y1 = protein(j,k,2)
                        if ((x1 .eq. x) .AND. (y1 .eq. y))	goto 87		! overlap with the other proteins
                     enddo
                  endif
               enddo
            endif
         enddo

         status_old=status	! memorize the old configurations

         do j=1,protein_length(n)
            protein_old(j,1)=protein(n,j,1)
            protein_old(j,2)=protein(n,j,2)
         enddo

         energy=0.d0
         call protein_energy()
         call Bonds()
         call potential_energy()
         energy_old=energy+EPS*(v12-v6)-jj_bulk_eff*DR_HB_f*HB_bulk-jj_phob_eff*DR_HB_f*HB_phob-jj_phil_eff*DR_HB_f*HB_phil &
         -jj_mix_eff*DR_HB_f*HB_mix &
         -jsigma*Nim -jsigma_phob*Nim_phob

         ! FLIP THE CHAIN AND CHANGE THE STATUS
         ! note that I cannot simply exchange the states because the new protein structure could overlap with the old one

         do i=node+1,protein_length(n)
            x=protein(n,i,1)
            y=protein(n,i,2)
            !x1=protein_aux(i,1)
            !y1=protein_aux(i,2)
            !if (status(x1,y1,0) .gt. 0) then	! assigne new values to the old protein sites rotating the water molecules
            !	do k=0,3
            !		status(x,y,k)=status(x1,y1,MOD(k+rot,4))
            !	enddo
            !else
            do k=0,3
               status(x,y,k)=i_dran(q)	! assigne new random values to the old protein sites
            enddo
            !endif
         enddo

         do j=1,protein_length(n)	! define the new protein coordinates
            protein(n,j,1)=protein_aux(j,1)
            protein(n,j,2)=protein_aux(j,2)
         enddo

         do i=node+1,protein_length(n)	! define the new protein state
            x=protein(n,i,1)
            y=protein(n,i,2)
            do k=0,3
               status(x,y,k)=protein_state(n,i)
            enddo
         enddo

         call protein_surface() ! find the new surface

         energy=0.d0
         call protein_energy()
         call Bonds()
         call potential_energy()
         energy_new=energy+EPS*(v12-v6)-jj_bulk_eff*DR_HB_f*HB_bulk-jj_phob_eff*DR_HB_f*HB_phob & 
         -jj_phil_eff*DR_HB_f*HB_phil-jj_mix_eff*DR_HB_f*HB_mix &
         -jsigma*Nim -jsigma_phob*Nim_phob

         delta_h=energy_new-energy_old
         acc_pivot=acc_pivot+1


         if (delta_h .gt. 0)	then
            z=dran_u()
            metropolis=dexp(-delta_h/Temp)
            if (z .gt. metropolis) then
               status=status_old
               do j=1,protein_length(n)	! define the new protein coordinates
                  protein(n,j,1)=protein_old(j,1)
                  protein(n,j,2)=protein_old(j,2)
               enddo
               call protein_surface()
               acc_pivot=acc_pivot-1
            endif
         endif

87       useless_flag=0
      endif
   enddo

end SUBROUTINE


! *****************
! * PROTEIN SHIFT *
! *****************

SUBROUTINE protein_shift()

   use common_variables_folding

   implicit none
   INTEGER :: x, y, x1, y1, i, k, nn, n, useless_flag, j, dx, dy, flag_overlap, jj
   INTEGER :: xx, yy, st, particle, pp, rand_step
   INTEGER(KIND=8) :: i_dran
   REAL (KIND=8) :: energy_new, energy_old, delta_h, z, metropolis, dran_u

   do nn=1,number_proteins

      n=i_dran(number_proteins)

      if ( protein_move(n) .eq. 1 )	then	! Check if the protein is fixed or mobile
         do j=1,protein_length(n)
            protein_aux(j,1)=protein(n,j,1)
            protein_aux(j,2)=protein(n,j,2)
         enddo

         if ( flux .eq. 0 )	then
            dx = i_dran(L + 1) - (L/2 + 1)
            dy = i_dran(L + 1) - (L/2 + 1)
         else
            !DR: Original 2021 version:
            !dx = i_dran(3) - 2! In this case we assume the presence of a flux that forces the protein to move in a fixed direction
            !dy = 1
            !DR: modified to:
            dx = i_dran(2) - 1
            dy = i_dran(2) - 1
            !for similarity with 2016
         endif

         if ((dx .ne. 0) .OR. (dy .ne. 0))	then
            ! shift the auxiliar protein
            do j=1,protein_length(n)
               protein_aux(j,1) = pbc(protein_aux(j,1) + dx)
               protein_aux(j,2) = pbc(protein_aux(j,2) + dy)
            enddo

            !check if the protein is overlapping with itself or with a different protein
            do j=1,protein_length(n)
               x = protein_aux(j,1)
               y = protein_aux(j,2)
               if (status(x,y,0) .lt. 0)	then	!overlap with a protein ... establish which
                  flag_overlap = 0
                  do jj=1,protein_length(n)
                     x1 = protein(n,jj,1)
                     y1 = protein(n,jj,2)
                     if ((x .eq. x1) .AND. (y .eq. y1))	then	! the new protein conformation overlaps with the old one. Shift is possible
                        flag_overlap = 1
                     endif
                  enddo
                  if (flag_overlap .eq. 0)	then
                     goto 807
                  endif
               endif
            enddo

            status_old=status	! memorize the old configurations
            do j=1,protein_length(n)
               protein_old(j,1)=protein(n,j,1)
               protein_old(j,2)=protein(n,j,2)
            enddo

            energy=0.d0
            call protein_energy()
            call Bonds()
            call potential_energy()
            energy_old=energy+EPS*(v12-v6)-jj_bulk_eff*DR_HB_f*HB_bulk-jj_phob_eff*DR_HB_f*HB_phob-jj_phil_eff*DR_HB_f*HB_phil &
            -jj_mix_eff*DR_HB_f*HB_mix &
            -jsigma*Nim -jsigma_phob*Nim_phob

            ! SHIFT THE CHAIN AND CHANGE THE STATUS
            ! note that I cannot simply exchange the states because the new protein structure could overlap with the old one

            do i=1,protein_length(n)
               x=protein(n,i,1)
               y=protein(n,i,2)
               do k=0,3
                  status(x,y,k)=i_dran(q)	! assigne new random values to the old protein sites
               enddo
            enddo

            do j=1,protein_length(n)	! define the new protein coordinates
               protein(n,j,1)=protein_aux(j,1)
               protein(n,j,2)=protein_aux(j,2)
            enddo

            do i=1,protein_length(n)	! define the new protein state
               x=protein(n,i,1)
               y=protein(n,i,2)
               do k=0,3
                  status(x,y,k)=protein_state(n,i)
               enddo
            enddo

            call protein_surface() ! find the new surface

            energy=0.d0
            call protein_energy()
            call Bonds()
            call potential_energy()
            energy_new=energy+EPS*(v12-v6)-jj_bulk_eff*DR_HB_f*HB_bulk-jj_phob_eff*DR_HB_f*HB_phob & 
            -jj_phil_eff*DR_HB_f*HB_phil-jj_mix_eff*DR_HB_f*HB_mix &
            -jsigma*Nim -jsigma_phob*Nim_phob

            delta_h=energy_new-energy_old
            acc_shift=acc_shift+1
            if (delta_h .gt. 0)	then
               z=dran_u()
               metropolis=dexp(-delta_h/Temp)
               if (z .gt. metropolis) then
                  status=status_old
                  do j=1,protein_length(n)	! define the new protein coordinates
                     protein(n,j,1)=protein_old(j,1)
                     protein(n,j,2)=protein_old(j,2)
                  enddo
                  call protein_surface()
                  acc_shift=acc_shift-1
               endif
            endif

         endif

807      useless_flag=0
      endif
   enddo

end SUBROUTINE


! ********************
! * PROTEIN ROTATION *
! ********************

SUBROUTINE protein_rotation()

   use common_variables_folding

   implicit none
   INTEGER :: x, y, x1, y1, i, k, nn, n, useless_flag, j, theta, flag_overlap, jj
   INTEGER :: xx, yy, st, pp, rand_step, particle, cmx, cmy, dx, dy
   INTEGER(KIND=8) :: i_dran
   REAL (KIND=8) :: energy_new, energy_old, delta_h, z, metropolis, dran_u

   do nn=1,number_proteins

      n=i_dran(number_proteins)

      if ( protein_move(n) .eq. 1 )	then	! Check if the protein is fixed or mobile

         theta = i_dran(4) ! theta =1,2,3,4 -> rotation of 0,90,180,270 degrees

         if (theta .gt. 1)	then
            ! shift the protein out of the pbc to compute the cm
            protein_aux(1,1) = protein(n,1,1)
            protein_aux(1,2) = protein(n,1,2)
            do j=2,protein_length(n)
               dx = protein(n,j,1) - protein(n,j-1,1)
               dy = protein(n,j,2) - protein(n,j-1,2)

               if (dx .eq. 0)	then
                  protein_aux(j,1) = protein_aux(j-1,1)
               endif
               if (dx .eq. -1)	then
                  protein_aux(j,1) = protein_aux(j-1,1) - 1
               endif
               if (dx .eq. 1)	then
                  protein_aux(j,1) = protein_aux(j-1,1) + 1
               endif
               if (dx .gt. 1)	then
                  protein_aux(j,1) = protein_aux(j-1,1) - 1
               endif
               if (dx .lt. -1)	then
                  protein_aux(j,1) = protein_aux(j-1,1) + 1
               endif

               if (dy .eq. 0)	then
                  protein_aux(j,2) = protein_aux(j-1,2)
               endif
               if (dy .eq. -1)	then
                  protein_aux(j,2) = protein_aux(j-1,2) - 1
               endif
               if (dy .eq. 1)	then
                  protein_aux(j,2) = protein_aux(j-1,2) + 1
               endif
               if (dy .gt. 1)	then
                  protein_aux(j,2) = protein_aux(j-1,2) - 1
               endif
               if (dy .lt. -1)	then
                  protein_aux(j,2) = protein_aux(j-1,2) + 1
               endif

            enddo

            !DR: this was original code:
            ! cmx = 0
            ! cmy = 0
            ! do j=1,protein_length(n)
            ! 	cmx = cmx + protein_aux(j,1)
            ! 	cmy = cmy + protein_aux(j,2)
            ! enddo
            ! cmx = pbc(cmx/protein_length(n))
            ! cmy = pbc(cmy/protein_length(n))
            !DR: this is my modification, which amounts to choosing a random monomer as the pivot point for the rotation
            cmx = pbc(protein_aux(i_dran(protein_length(n)),1))
            cmy = pbc(protein_aux(i_dran(protein_length(n)),2))

            ! rotate the auxiliar protein
            if ( theta .eq. 2 )	then	! 90 degrees
               do j=1,protein_length(n)
                  x = pbc(protein(n,j,1) - cmx)
                  y = pbc(protein(n,j,2) - cmy)
                  protein_aux(j,1) = pbc(-y + cmx)	! x' = -y
                  protein_aux(j,2) = pbc( x + cmy)	! y' = x
               enddo
            endif

            if ( theta .eq. 3 )	then	! 180 degrees
               do j=1,protein_length(n)
                  x = pbc(protein(n,j,1) - cmx)
                  y = pbc(protein(n,j,2) - cmy)
                  protein_aux(j,1) = pbc(-x + cmx)	! x' = -x
                  protein_aux(j,2) = pbc(-y + cmy)	! y' = -y
               enddo
            endif

            if ( theta .eq. 4 )	then	! 270 degrees
               do j=1,protein_length(n)
                  x = pbc(protein(n,j,1) - cmx)
                  y = pbc(protein(n,j,2) - cmy)
                  protein_aux(j,1) = pbc( y + cmx)	! x' = y
                  protein_aux(j,2) = pbc(-x + cmy)	! y' = -x
               enddo
            endif

            do j=1,protein_length(n)-1
               xx=protein_aux(j+1,1)-protein_aux(j,1)
               yy=protein_aux(j+1,2)-protein_aux(j,2)
               if (iabs(xx) .gt. 1) xx=xx-sign(L,xx)
               if (iabs(yy) .gt. 1) yy=yy-sign(L,yy)
               if ((iabs(xx) .gt. 1) .OR. (iabs(yy) .gt. 1))	then
                  write (*,*) "During rotation, Broken protein "
                  call exit()
               endif
            enddo

            !check if the protein is overlapping with itself or with a different protein
            do j=1,protein_length(n)
               x = protein_aux(j,1)
               y = protein_aux(j,2)
               if (status(x,y,0) .lt. 0)	then	!overlap with a protein ... establish which
                  flag_overlap = 0
                  do jj=1,protein_length(n)
                     x1 = protein(n,jj,1)
                     y1 = protein(n,jj,2)
                     if ((x .eq. x1) .AND. (y .eq. y1))	then	! the new protein conformation overlaps with the old one. rotation is possible
                        flag_overlap = 1
                     endif
                  enddo
                  if (flag_overlap .eq. 0)	then
                     !	write(*,*) "REJECTED"
                     goto 807
                  endif
               endif
            enddo

            status_old=status	! memorize the old configurations
            do j=1,protein_length(n)
               protein_old(j,1)=protein(n,j,1)
               protein_old(j,2)=protein(n,j,2)
            enddo

            energy=0.d0
            call protein_energy()
            call Bonds()
            call potential_energy()
            energy_old=energy+EPS*(v12-v6)-jj_bulk_eff*DR_HB_f*HB_bulk-jj_phob_eff*DR_HB_f*HB_phob-jj_phil_eff*DR_HB_f*HB_phil & 
            -jj_mix_eff*DR_HB_f*HB_mix &
            -jsigma*Nim -jsigma_phob*Nim_phob

            ! ROTATE THE CHAIN AND CHANGE THE STATUS
            ! note that I cannot simply exchange the states because the new protein structure could overlap with the old one

            do i=1,protein_length(n)
               x=protein(n,i,1)
               y=protein(n,i,2)
               do k=0,3
                  status(x,y,k)=i_dran(q)	! assigne new random values to all the old protein sites
                  ! (some of them will be in case replaced by the protein_state)
               enddo
            enddo

            do j=1,protein_length(n)	! define the new protein coordinates
               protein(n,j,1)=protein_aux(j,1)
               protein(n,j,2)=protein_aux(j,2)
            enddo

            do i=1,protein_length(n)	! define the new protein state
               x=protein(n,i,1)
               y=protein(n,i,2)
               do k=0,3
                  status(x,y,k)=protein_state(n,i)
               enddo
            enddo

            call protein_surface() ! find the new surface

            energy=0.d0
            call protein_energy()
            call Bonds()
            call potential_energy()
            energy_new=energy+EPS*(v12-v6)-jj_bulk_eff*DR_HB_f*HB_bulk-jj_phob_eff*DR_HB_f*HB_phob & 
            -jj_phil_eff*DR_HB_f*HB_phil-jj_mix_eff*DR_HB_f*HB_mix &
            -jsigma*Nim -jsigma_phob*Nim_phob

            delta_h=energy_new-energy_old
            acc_protein_rotation=acc_protein_rotation+1
            if (delta_h .gt. 0)	then
               z=dran_u()
               metropolis=dexp(-delta_h/Temp)
               if (z .gt. metropolis) then
                  status=status_old
                  do j=1,protein_length(n)	! define the new protein coordinates
                     protein(n,j,1)=protein_old(j,1)
                     protein(n,j,2)=protein_old(j,2)
                  enddo
                  call protein_surface()
                  acc_protein_rotation=acc_protein_rotation-1
                  !	write(*,*) "REJECTED"
               endif
            endif

         endif

807      useless_flag=0
      endif
   enddo

end SUBROUTINE


! **************
! * CRANKSHAFT *
! **************

SUBROUTINE crankshaft()

   use common_variables_folding

   implicit none
   INTEGER :: node1, node2, x1, x2, y1, y2, dx, dy, flag_move, delta, i, j, k, flag, n, nn, useless_flag
   INTEGER :: xx, yy, pp, st, particle, rand_step
   INTEGER(KIND=8) :: i_dran
   REAL (KIND=8) :: energy_new, energy_old, delta_h, z, metropolis, dran_u

   do nn=1,number_proteins

      n=i_dran(number_proteins)

      if ( protein_move(n) .eq. 1 )	then	! Check if the protein is fixed or mobile

         node1=i_dran(protein_length(n)-3)
         node2=node1+3

         dx=protein(n,node1,1)-protein(n,node2,1)
         dy=protein(n,node1,2)-protein(n,node2,2)

         flag_move=0
         if (dx .eq. 0)	then
            do j=1,protein_length(n)
               protein_aux(j,1)=protein(n,j,1)
               protein_aux(j,2)=protein(n,j,2)
            enddo
            x1=protein(n,node1,1)
            protein_aux(node1+1,1)=pbc(2*x1 - protein_aux(node1+1,1))
            protein_aux(node1+2,1)=pbc(2*x1 - protein_aux(node1+2,1))
            delta=iabs(protein_aux(node1+1,1)-protein(n,node1+1,1))+iabs(protein_aux(node1+2,1)-protein(n,node1+2,1))
            if (delta .ne. 0) flag_move=1
         elseif (dy .eq. 0)	then
            do j=1,protein_length(n)
               protein_aux(j,1)=protein(n,j,1)
               protein_aux(j,2)=protein(n,j,2)
            enddo
            y1=protein(n,node1,2)
            protein_aux(node1+1,2)=pbc(2*y1 - protein_aux(node1+1,2))
            protein_aux(node1+2,2)=pbc(2*y1 - protein_aux(node1+2,2))
            delta=iabs(protein_aux(node1+1,2)-protein(n,node1+1,2))+iabs(protein_aux(node1+2,2)-protein(n,node1+2,2))
            if (delta .ne. 0) flag_move=1
         endif

         if (flag_move .eq. 1)	then

            delta=0
            do i=1,2	! check if the new protein overlaps with itself
               do j=1,protein_length(n)
                  if (j .ne. node1+i)	then
                     delta=iabs(protein_aux(node1+i,1)-protein_aux(j,1)) + iabs(protein_aux(node1+i,2)-protein_aux(j,2))
                     if (delta .eq. 0)	then
                        goto 666
                     endif
                  endif
               enddo
            enddo

            do i=node1,node2	! check of the overlap with other proteins to establish if the rotation is possible
               x1 = protein_aux(i,1)
               y1 = protein_aux(i,2)
               if (status(x1,y1,0) .lt. 0)	then	! possible overlap
                  do j=1,number_proteins		! look all the other proteins proteins
                     if (j .ne. n)	then	! exclude the moving protein that has been already checked
                        do k=1,protein_length(j)
                           x2 = protein(j,k,1)
                           y2 = protein(j,k,2)
                           if ((x1 .eq. x2) .AND. (y1 .eq. y2))	goto 666		! overlap with the other proteins
                        enddo
                     endif
                  enddo
               endif
            enddo

            energy=0.d0
            call protein_energy()
            call Bonds()
            call potential_energy()
            energy_old=energy+EPS*(v12-v6)-jj_bulk_eff*DR_HB_f*HB_bulk-jj_phob_eff*DR_HB_f*HB_phob-jj_phil_eff*DR_HB_f*HB_phil &
            -jj_mix_eff*DR_HB_f*HB_mix-jsigma*Nim -jsigma_phob*Nim_phob

            do j=1,protein_length(n)
               protein_old(j,1)=protein(n,j,1)
               protein_old(j,2)=protein(n,j,2)
            enddo

            status_old=status

            do i=1,2
               x1=protein(n,node1+i,1)
               y1=protein(n,node1+i,2)
               x2=protein_aux(node1+i,1)
               y2=protein_aux(node1+i,2)
               do j=0,3
                  flag=status(x1,y1,j)
                  status(x1,y1,j)=status(x2,y2,j)
                  status(x2,y2,j)=flag
               enddo
            enddo

            do j=1,protein_length(n)
               protein(n,j,1)=protein_aux(j,1)
               protein(n,j,2)=protein_aux(j,2)
            enddo

            call protein_surface()

            energy=0.d0
            call protein_energy()
            call Bonds()
            call potential_energy()
            energy_new=energy+EPS*(v12-v6)-jj_bulk_eff*DR_HB_f*HB_bulk-jj_phob_eff*DR_HB_f*HB_phob-jj_phil_eff*DR_HB_f*HB_phil &
            -jj_mix_eff*DR_HB_f*HB_mix-jsigma*Nim -jsigma_phob*Nim_phob

            delta_h=energy_new-energy_old
            acc_crankshaft=acc_crankshaft+1
            if (delta_h .gt. 0)	then
               z=dran_u()
               metropolis=dexp(-delta_h/Temp)
               if (z .gt. metropolis) then
                  status=status_old
                  do j=1,protein_length(n)
                     protein(n,j,1)=protein_old(j,1)
                     protein(n,j,2)=protein_old(j,2)
                  enddo
                  call protein_surface()
                  acc_crankshaft=acc_crankshaft-1
               endif
            endif
         endif
666      useless_flag = 0
      endif
   enddo

end SUBROUTINE


! ******************
! *  TOTAL BONDS   *
! ******************
! Calculate the total number of HB and the surface bounds!

SUBROUTINE Bonds()  ! IMPORTANT - YOU EVER MUST CALL "protein_surface" BEFORE TO CALL THIS SUBROUTINE

   use common_variables_folding

   implicit none
   INTEGER :: st, m, i, x, y, k, neig_st, su, neig_su, neig_x, neig_y

   HB_bulk=0;	HB_phil=0;	HB_phob=0;	HB_mix=0;	bond_phil=0;	bond_phob=0;	Nim=0;	Nim_phob=0;
   do m=0,num_site-1
      x=rx(m)
      y=ry(m)
      su=surface(x,y)
      do i=0,1
         st=status(x,y,i)
         neig_x=pbc(x+cx(i))
         neig_y=pbc(y+cy(i))
         neig_st=status(neig_x,neig_y,neig_arm(i))
         neig_su=surface(neig_x,neig_y)
         if ((st .eq. neig_st) .AND. (st .gt. 0))	then	! both are water molecules
            if ((su .eq. 0) .OR. (neig_su .eq. 0))	then
               HB_bulk=HB_bulk+1
            else
               if ((su .eq. 1) .AND. (neig_su .eq. 1))	then
                  HB_phob=HB_phob+1
               else
                  if ((su .eq. 2) .AND. (neig_su .eq. 2))	then
                     HB_phil=HB_phil+1
                  else
                     HB_mix=HB_mix+1
                  endif
               endif
            endif
         else
            if (st*neig_st .le. 0)	then	! one site is residue and the other is a water molecule
               if (matrix(st,neig_st) .ge. 0)	then	! hydrophobic residue
                  bond_phob=bond_phob+1
               else	! one site is a hydrophilic residue
                  bond_phil=bond_phil+1
               endif
            endif
         endif
      enddo
      st=status(x,y,0)
      if (st .gt. 0)	then	! it is a water molecule
         if ((su .ne. 1) .AND. (su .ne. 3))	then
            do i=0,2
               do k=i+1,3
                  if (status(x,y,i) .eq. status(x,y,k))	Nim=Nim+1
               enddo
            enddo
         else
            do i=0,2
               do k=i+1,3
                  if (status(x,y,i) .eq. status(x,y,k))	Nim_phob=Nim_phob+1
               enddo
            enddo
         endif
      endif
   enddo ! no. of particles

end SUBROUTINE


SUBROUTINE Bonds_DR()

   use common_variables_folding

   implicit none
   INTEGER(KIND=8) :: st, m, i, x, y, k, neig_st, su, neig_su, neig_x, neig_y

   HB_bulk=0;	HB_phil=0;	HB_phob=0;	HB_mix=0;	bond_phil=0;	bond_phob=0;	Nim=0;	Nim_phob=0

   open(1133,file="hor_HB.dat",status="unknown",access="append")
   open(1134,file="ver_HB.dat",status="unknown",access="append")
   open(1135,file="coop_bonds.dat",status="unknown",access="append")
   open(1136,file="coop_bonds_PHO.dat",status="unknown",access="append")

   do m=0,num_site-1
      x=rx(m)
      y=ry(m)
      su=surface(x,y)

      i=0
      st=status(x,y,i)
      neig_x=pbc(x+cx(i))
      neig_y=pbc(y+cy(i))
      neig_st=status(neig_x,neig_y,neig_arm(i))
      neig_su=surface(neig_x,neig_y)
      if ((st .gt. 0) .AND. (neig_st .gt. 0)) then		! both are water molecules
         if (st .eq. neig_st)	then !DR: there is HB, check wich type according to surface
            if ((su .eq. 0) .OR. (neig_su .eq. 0))	then !DR: at least one of the water molecules is bulk
               write(1134,'(I3)', advance="no") 1 ! bulk HB
            else !DR:
               if ((su .eq. 1) .AND. (neig_su .eq. 1))	then
                  write(1134,'(I3)', advance="no") 2 ! phob HB
               else
                  if ((su .eq. 2) .AND. (neig_su .eq. 2))	then
                     write(1134,'(I3)', advance="no") 3 ! phil HB
                  else
                     write(1134,'(I3)', advance="no") 4 ! mix HB
                  endif
               endif
            endif
         else !there is no HB here
            write(1134,'(I3)', advance="no") 0 ! no HB
         endif
      else
         write(1134,'(I3)', advance="no") 0 ! no HB
         if (st*neig_st .le. 0)	then	! one site is residue and the other is a water molecule
            if (matrix(st,neig_st) .ge. 0)	then	! one site is a hydrophobic residue
               bond_phob=bond_phob+1
            else	! one site is a hydrophilic residue
               bond_phil=bond_phil+1
            endif
         endif
      endif


      i=1
      st=status(x,y,i)
      neig_x=pbc(x+cx(i))
      neig_y=pbc(y+cy(i))
      neig_st=status(neig_x,neig_y,neig_arm(i))
      neig_su=surface(neig_x,neig_y)
      if ((st .gt. 0) .AND. (neig_st .gt. 0)) then		! both are water molecules
         if (st .eq. neig_st)	then ! there is HB, check wich type according to surface
            if ((su .eq. 0) .OR. (neig_su .eq. 0))	then
               write(1133,'(I3)', advance="no") 1 ! bulk HB
            else
               if ((su .eq. 1) .AND. (neig_su .eq. 1))	then
                  write(1133,'(I3)', advance="no") 2 ! phob HB
               else
                  if ((su .eq. 2) .AND. (neig_su .eq. 2))	then
                     write(1133,'(I3)', advance="no") 3 ! phil HB
                  else
                     write(1133,'(I3)', advance="no") 4 ! mix HB
                  endif
               endif
            endif
         else !there is no HB here
            write(1133,'(I3)', advance="no") 0 ! no HB
         endif
      else
         write(1133,'(I3)', advance="no") 0 ! no HB
         if (st*neig_st .le. 0)	then	! one site is residue and the other is a water molecule
            if (matrix(st,neig_st) .ge. 0)	then	! one site is a hydrophobic residue
               bond_phob=bond_phob+1
            else	! one site is a hydrophilic residue
               bond_phil=bond_phil+1
            endif
         endif
      endif

      Nim=0; Nim_phob=0
      st=status(x,y,0)
      if (st .gt. 0)	then	! it is a water molecule
         if ((su .ne. 1) .AND. (su .ne. 3))	then
            do i=0,2
               do k=i+1,3
                  if (status(x,y,i) .eq. status(x,y,k)) Nim=Nim+1
               enddo
            enddo
         else
            do i=0,2
               do k=i+1,3
                  if (status(x,y,i) .eq. status(x,y,k)) Nim_phob=Nim_phob+1
               enddo
            enddo
         endif
      endif
      write(1135,'(I3)', advance="no") Nim
      write(1136,'(I3)', advance="no") Nim_phob
      
   enddo ! no. of particles

   write(1133,'(A1)') ""
   write(1134,'(A1)') ""
   write(1135,'(A1)') ""
   write(1136,'(A1)') ""
   call flush(1133)
   call flush(1134)
   call flush(1135)
   call flush(1136)
   close(1133)
   close(1134)
   close(1135)
   close(1136)
end SUBROUTINE


! ***************
! *  ROTATION   *
! ***************

SUBROUTINE rotation(x,y,arm,st)

   use common_variables_folding

   implicit none
   INTEGER :: x, y, arm, st ! dummy variables
   INTEGER :: i, su, new_st, neig_x, neig_y, neig_st, neig_su, kkk, delta_jsigma
   INTEGER(KIND=8) :: i_dran
   REAL(KIND=8):: deltah, z, metropolis, dran_u

72 new_st=i_dran(q)

   if (new_st .eq. st) then
      goto 72
   endif

   su = surface(x,y)
   neig_x=pbc(x+cx(arm)) ! neigbohring particle in the arm's direction
   neig_y=pbc(y+cy(arm))
   neig_st=status(neig_x,neig_y,neig_arm(arm))
   neig_su = surface(neig_x,neig_y)

   deltah=0.d0
   if (neig_st .gt. 0)	then	! the neighbour is a water molecule
      if (st .eq. neig_st)	then
         if ((su .eq. 0) .OR. (neig_su .eq. 0))	then
            deltah=jj_bulk_eff
         else
            if ((su .eq. 1) .AND. (neig_su .eq. 1))	then
               deltah=jj_phob_eff
            else
               if ((su .eq. 2) .AND. (neig_su .eq. 2))	then
                  deltah=jj_phil_eff
               else
                  deltah=jj_mix_eff
               endif
            endif
         endif
      else
         if (new_st .eq. neig_st)	then
            if ((su .eq. 0) .OR. (neig_su .eq. 0))	then
               deltah=-jj_bulk_eff
            else
               if ((su .eq. 1) .AND. (neig_su .eq. 1))	then
                  deltah=-jj_phob_eff
               else
                  if ((su .eq. 2) .AND. (neig_su .eq. 2))	then
                     deltah=-jj_phil_eff
                  else
                     deltah=-jj_mix_eff
                  endif
               endif
            endif
         endif
      endif
   endif

   delta_jsigma=0 ! jsigma interaction
   do i=1,3
      kkk=status(x,y,mod(arm+i,4)) ! the other three arms of the molecule
      if (st .eq. kkk) delta_jsigma=delta_jsigma-1
      if (new_st .eq. kkk) delta_jsigma=delta_jsigma+1
   enddo

   if ((su .ne. 1) .AND. (su .ne. 3))	then	! the molecule is at the hydrophilic interface
      deltah=deltah-dble(delta_jsigma)*jsigma
   else
      deltah=deltah-dble(delta_jsigma)*jsigma_phob
   endif

   if (deltah .le. 0.d0) then ! metropolis
      status(x,y,arm)=new_st
      acc_rotation=acc_rotation+1
   else
      z=dran_u()
      metropolis=dexp(-deltah/Temp)
      if (z .lt. metropolis) then
         status(x,y,arm)=new_st
         acc_rotation=acc_rotation+1
      endif
   endif ! end metropolis

end SUBROUTINE



! ********************************
! *   PROTEIN ENERGY ENTHALPY    *
! ********************************

SUBROUTINE 	protein_energy()

   use common_variables_folding

   implicit none
   INTEGER(KIND=8) :: i, st, neig_st, x, y, n
   REAL (KIND=8) :: ene

   do n=1,number_proteins
      do i=1,protein_length(n)
         st=status(protein(n,i,1),protein(n,i,2),0)

         x=pbc(protein(n,i,1)+1)
         y=protein(n,i,2)
         neig_st=status(x,y,3)	! status of the right site
         if (neig_st .gt. 0)	then	! it is a water molecule
            ene=matrix(st,neig_st)	! water-protein energy
            energy=energy+ene
         else
            ene=matrix(st,neig_st)
            energy=energy+ene/2.d0
         endif

         x=protein(n,i,1)
         y=pbc(protein(n,i,2)+1)
         neig_st=status(x,y,0)	! status of the upper site
         if (neig_st .gt. 0)	then	! it is a water molecule
            ene=matrix(st,neig_st)
            energy=energy+ene
         else
            ene=matrix(st,neig_st)
            energy=energy+ene/2.d0
         endif

         x=pbc(protein(n,i,1)-1)
         y=protein(n,i,2)
         neig_st=status(x,y,1)	! status of the left site
         if (neig_st .gt. 0)	then	! it is a water molecule
            ene=matrix(st,neig_st)
            energy=energy+ene
         else
            ene=matrix(st,neig_st)
            energy=energy+ene/2.d0
         endif

         x=protein(n,i,1)
         y=pbc(protein(n,i,2)-1)
         neig_st=status(x,y,2)	! status of the down site
         if (neig_st .gt. 0)	then	! it is a water molecule
            ene=matrix(st,neig_st)
            energy=energy+ene
         else
            ene=matrix(st,neig_st)
            energy=energy+ene/2.d0
         endif
      enddo

      ! subtract the energy of the amonoacids along the chain
      do i=1,protein_length(n)-1
         st=status(protein(n,i,1),protein(n,i,2),0)
         neig_st=status(protein(n,i+1,1),protein(n,i+1,2),0)
         ene=matrix(st,neig_st)
         energy=energy-ene
      enddo
   enddo

end SUBROUTINE



! ***************
! *   VOLUME    *
! ***************

SUBROUTINE volume()

   use common_variables_folding

   implicit none
   REAL(KIND=8) :: v6new, v12new, aux, aux3, deltav, extra, deltal, deltah, new_latt2, z, metropolis, dran_u
   REAL(KIND=8) :: DR_latt,DR_HB_f_new,DR_HB_extra_contribution

5  deltal=lsample*(2.d0*dran_u()-1.d0) ! -lsample < deltal < lsample

   new_latt2=latt2*(1.d0+deltal)**2 ! it's the volume without hb's that is updated

   DR_latt=(1.d0+deltal)*latt_size
   DR_HB_f_new=1.0d0/( dexp((DR_latt-DR_HB_l_threshold)*DR_HB_beta)  +1)

   if(new_latt2/dL2 .le. 1.d0) then	!  volume per particle greater than 1
      goto 5
   endif

   aux=latt2/new_latt2
   aux3=aux*aux*aux
   v6new=v6*aux3
   v12new=v12*aux3*aux3
   deltav=EPS*(v12new-v6new-v12+v6)
   extra=num_site*Temp*dlog(aux)

   DR_HB_extra_contribution=-jj_bulk*(DR_HB_f_new-DR_HB_f)*HB_bulk-jj_phob*(DR_HB_f_new-DR_HB_f)*HB_phob &
   -jj_phil*(DR_HB_f_new-DR_HB_f)*HB_phil-jj_mix*(DR_HB_f_new-DR_HB_f)*HB_mix

   deltah=deltav+extra+Press*(new_latt2-latt2)  &
   +DR_HB_extra_contribution

   ! write(*,*) DR_latt,latt, DR_HB_f_new, DR_HB_f, DR_HB_extra_contribution
   ! deltah=deltav+extra+Press*(new_latt2-latt2)

   if (deltah .le. 0.d0) then ! metropolis
      latt2=new_latt2
      latt=(1.d0+deltal)*latt
      latt_size2=latt2/dL2
      v12=v12new
      v6=v6new
      acc_volume=acc_volume+1
   else
      z=dran_u()
      metropolis=dexp(-deltah/Temp)
      if (z .lt. metropolis) then
         latt2=new_latt2
         latt=(1.d0+deltal)*latt
         latt_size2=latt2/dL2
         v12=v12new
         v6=v6new
         acc_volume=acc_volume+1
      endif
   endif

end SUBROUTINE


!*****************************
!  RANDOM NUMBER GENERATOR   *
!*****************************

SUBROUTINE dran_ini(iseed0)

   use random_variables_folding

   implicit none
   INTEGER(KIND=8) :: iseed0  ! dummy variable
   INTEGER(KIND=8):: i, j
   INTEGER(KIND=8), parameter :: nbit=31 !np=14
   REAL(KIND=8) :: dseed, rand_xx

   dseed=iseed0
   do i=1,ip
      ix(i)=0
      do j=0,nbit-1
         if (rand_xx(dseed) .lt. 0.5d0) ix(i)=ibset(ix(i),j)
      enddo
   enddo
   ic=0

end SUBROUTINE dran_ini

INTEGER(KIND=8) FUNCTION i_dran(n)

   use random_variables_folding

   implicit none
   INTEGER :: n !dummy variable
   INTEGER(KIND=8) :: i_ran
   INTEGER, parameter ::  iq=418, is=ip-iq

   ic = ic + 1
   if (ic .gt. ip)	ic=1
   if (ic .gt. iq)	then
      ix(ic)=ieor(ix(ic),ix(ic-iq))
   else
      ix(ic)=ieor(ix(ic),ix(ic+is))
   endif
   i_ran=ix(ic)
   if (n .gt. 0) i_dran=mod(i_ran,n) + 1

end FUNCTION i_dran

REAL(KIND=8) FUNCTION dran_u()

   use random_variables_folding

   implicit none
   INTEGER(KIND=8), parameter ::  iq=418, is=ip-iq
   REAL(KIND=8), parameter :: rmax=2147483647.0

   ic = ic + 1
   if (ic .gt. ip) ic = 1
   if (ic .gt. iq)	then
      ix(ic)=ieor(ix(ic),ix(ic-iq))
   else
      ix(ic)=ieor(ix(ic),ix(ic+is))
   endif
   dran_u=dfloat(ix(ic))/rmax

end FUNCTION dran_u

REAL(KIND=8) FUNCTION rand_xx(dseed)

   implicit none
   REAL(KIND=8) :: dseed  ! dummy variable
   REAL(KIND=8), parameter ::xm=2.d0**32, rm=1.d0/xm, a=69069.d0, c=1.d0

   dseed=mod(dseed*a+c,xm)
   rand_xx=dseed*rm

end FUNCTION rand_xx

SUBROUTINE StripSpaces(string)
   character(len=*) :: string
   integer :: stringLen
   integer :: last, actual

   stringLen = len (string)
   last = 1
   actual = 1

   do while (actual < stringLen)
      if (string(last:last) == ' ') then
         actual = actual + 1
         string(last:last) = string(actual:actual)
         string(actual:actual) = ' '
      else
         last = last + 1
         if (actual < last) &
            actual = last
      endif
   end do

end SUBROUTINE
