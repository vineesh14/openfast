@test
subroutine test_BD_FEinternalForceQPweights()
   ! test branches
   ! - p = 1, invalid value
   ! - p = 2, boundaries only
   ! - p = 5, odd number
   ! - p = 6, even number
   ! - p = 97, large, prime number

   use pFUnit_mod
   use BeamDyn_Subs
   use NWTC_Num
   use test_tools
   
   implicit none
   
   integer(IntKi)             :: idx_FE, idx_qp, nelem, i
   type(BD_ParameterType)     :: p
   real(BDKi), allocatable    :: baseline_QPtWghtIntForceFE(:,:,:,:), baseline_QPrangeOverlapFE(:,:,:)
   logical,    allocatable    :: baseline_FEoutboardOfQPt(:,:)
   integer(IntKi)             :: ErrStat
   character                  :: ErrMsg
   character(1024)            :: testname
   real(BDKi)                 :: tolerance
   
   ! initialize NWTC_Num constants
   call SetConstants()
   
   tolerance = 1e-10
  
   
   ! --------------------------------------------------------------------------
   ! Gaussian quadrature has QP's inboard of the ends.  This changes the weightings distributions.
   testname = "Gauss quadrature (6 FE, 6 QP, 1 element): "
   p%nqp             = 6
   p%nodes_per_elem  = 6
   p%elem_total      = 1
   call SetupArrays
      ! setup the basic input arrays required
   p%GLL_Nodes       = (/   -1.0000_BDKi,  -0.7651_BDKi,  -0.2852_BDKi,          0.2852_BDKi,   0.7651_BDKi,   1.0000_BDKi /)
   p%QPtN            = (/       -0.9325_BDKi,  -0.6612_BDKi,  -0.2386_BDKi,   0.2386_BDKi,   0.6612_BDKi,   0.9325_BDKi /)
      ! expected results
   baseline_FEoutboardOfQPt = reshape( (/ .false.,.true. ,.false.,.false.,.false.,.false.,&
                                          .false.,.false.,.true. ,.false.,.false.,.false.,&
                                          .false.,.false.,.false.,.false.,.false.,.false.,&
                                          .false.,.false.,.false.,.true. ,.false.,.false.,&
                                          .false.,.false.,.false.,.false.,.true. ,.false.,&
                                          .false.,.false.,.false.,.false.,.false.,.true. /), &
                                        (/ p%nodes_per_elem, p%nqp /) )
   baseline_QPrangeOverlapFE = reshape( (/ 1, 1,&
                                           1, 2,&
                                           2, 4,&
                                           4, 5,&
                                           5, 5,&
                                           0, 0 /), &
                                        (/ 2, p%nodes_per_elem, p%elem_total /) )
   baseline_QPtWghtIntForceFE = reshape((/ 0.0675000000000000_BDKi   ,0.0516453372650202_BDKi   ,0.0000000000000000_BDKi   ,0.0000000000000000_BDKi   ,0.0000000000000000_BDKi   ,0.0000000000000000_BDKi   ,&
                                           0.1157546627349797_BDKi   ,0.0000000000000000_BDKi   ,0.0000000000000000_BDKi   ,0.0000000000000000_BDKi   ,0.0000000000000000_BDKi   ,0.0000000000000000_BDKi   ,&
                                           0.0000000000000000_BDKi   ,0.0840046627349797_BDKi   ,0.1672692853762423_BDKi   ,0.0000000000000000_BDKi   ,0.0000000000000000_BDKi   ,0.0000000000000000_BDKi   ,&
                                           0.0198953372650202_BDKi   ,0.2087307146237576_BDKi   ,0.0000000000000000_BDKi   ,0.0000000000000000_BDKi   ,0.0000000000000000_BDKi   ,0.0000000000000000_BDKi   ,&
                                           0.0000000000000000_BDKi   ,0.0000000000000000_BDKi   ,0.0440307146237576_BDKi   ,0.2386000000000000_BDKi   ,0.0025692853762423_BDKi   ,0.0000000000000000_BDKi   ,&
                                           0.0000000000000000_BDKi   ,0.0025692853762423_BDKi   ,0.2386000000000000_BDKi   ,0.0440307146237576_BDKi   ,0.0000000000000000_BDKi   ,0.0000000000000000_BDKi   ,&
                                           0.0000000000000000_BDKi   ,0.0000000000000000_BDKi   ,0.0000000000000000_BDKi   ,0.0000000000000000_BDKi   ,0.2087307146237576_BDKi   ,0.0198953372650202_BDKi   ,&
                                           0.0000000000000000_BDKi   ,0.0000000000000000_BDKi   ,0.0000000000000000_BDKi   ,0.1672692853762423_BDKi   ,0.0840046627349797_BDKi   ,0.0000000000000000_BDKi   ,&
                                           0.0000000000000000_BDKi   ,0.0000000000000000_BDKi   ,0.0000000000000000_BDKi   ,0.0000000000000000_BDKi   ,0.0000000000000000_BDKi   ,0.1157546627349797_BDKi   ,&
                                           0.0000000000000000_BDKi   ,0.0000000000000000_BDKi   ,0.0000000000000000_BDKi   ,0.0000000000000000_BDKi   ,0.0516453372650202_BDKi   ,0.0000000000000000_BDKi   ,&
                                           0.0000000000000000_BDKi   ,0.0000000000000000_BDKi   ,0.0000000000000000_BDKi   ,0.0000000000000000_BDKi   ,0.0000000000000000_BDKi   ,0.0000000000000000_BDKi   ,&
                                           0.0000000000000000_BDKi   ,0.0000000000000000_BDKi   ,0.0000000000000000_BDKi   ,0.0000000000000000_BDKi   ,0.0000000000000000_BDKi   ,0.0675000000000000_BDKi   /), &
                                        (/  p%nqp, 2, p%nodes_per_elem, p%elem_total /) )

   call BD_FEinternalForceQPweights(p, ErrStat, ErrMsg)

      ! check mapping flags
   do idx_qp = 1, p%nqp
      do idx_FE = 1, p%nodes_per_elem
         @assertEqual(baseline_FEoutboardOfQPt(idx_FE,idx_qp), p%FEoutboardOfQPt(idx_FE,idx_qp), testname)
      end do
   end do

      ! check QP mapping ranges
   do nelem = 1, p%elem_total
      do idx_FE = 1, p%nodes_per_elem
         do i=1,2
            @assertEqual(baseline_QPrangeOverlapFE(i,idx_FE,nelem), p%QPrangeOverlapFE(i,idx_FE,nelem), testname)
         end do
      end do
   end do

      ! check QPt weightings
   do nelem = 1, p%elem_total
      do idx_FE = 1, p%nodes_per_elem
         do i=1,2
            do idx_qp=1,p%nqp
               @assertEqual(baseline_QPtWghtIntForceFE(idx_qp,i,idx_FE,nelem), p%QPtWghtIntForceFE(idx_qp,i,idx_FE,nelem), tolerance, testname)
            enddo
         enddo
      enddo
   enddo

   call CleanUp


   ! --------------------------------------------------------------------------
   ! Trapezoidal quadrature.  This is the 5MW blade model.
   testname = "Trapezoidal quadrature 5MW blade (6 FE, 49 QP, 1 element): "
   p%nqp             = 49
   p%nodes_per_elem  = 6
   p%elem_total      = 1
   call SetupArrays
      ! setup the basic input arrays required
   p%GLL_Nodes       = (/   -1.0000_BDKi,  -0.7651_BDKi,  -0.2852_BDKi,          0.2852_BDKi,   0.7651_BDKi,   1.0000_BDKi /)
   p%QPtN            = (/ -1.0000_BDKi,-0.9935_BDKi,-0.9610_BDKi,-0.9285_BDKi,-0.8959_BDKi,-0.8634_BDKi,-0.8309_BDKi, &
                          -0.7984_BDKi,-0.7659_BDKi,-0.7333_BDKi,-0.7008_BDKi,-0.6683_BDKi,-0.6358_BDKi,-0.6033_BDKi, &
                          -0.5707_BDKi,-0.5382_BDKi,-0.5057_BDKi,-0.4732_BDKi,-0.4081_BDKi,-0.3431_BDKi,-0.2780_BDKi, &
                          -0.2130_BDKi,-0.1480_BDKi,-0.0829_BDKi,-0.0179_BDKi, 0.0472_BDKi, 0.1122_BDKi, 0.1772_BDKi, &
                           0.2423_BDKi, 0.3073_BDKi, 0.3724_BDKi, 0.4374_BDKi, 0.5024_BDKi, 0.5675_BDKi, 0.6325_BDKi, &
                           0.6976_BDKi, 0.7626_BDKi, 0.7951_BDKi, 0.8276_BDKi, 0.8602_BDKi, 0.8764_BDKi, 0.8927_BDKi, &
                           0.9089_BDKi, 0.9252_BDKi, 0.9415_BDKi, 0.9577_BDKi, 0.9740_BDKi, 0.9902_BDKi, 1.0000_BDKi /)
      ! expected results
   baseline_FEoutboardOfQPt = reshape( (/  .true., .false., .false., .false., .false., .false., .false., &
                                          .false., .false., .false., .false., .false., .false., .false., &
                                          .false., .false., .false., .false., .false., .false., .false., &
                                          .false., .false., .false., .false., .false., .false., .false., &
                                          .false., .false., .false., .false., .false., .false., .false., &
                                          .false., .false., .false., .false., .false., .false., .false., &
                                          .false., .false., .false., .false., .false., .false., .false., &
                                           .true., .false., .false., .false., .false., .false., .false., &
                                          .false., .false., .false., .false., .false., .false., .false., &
                                          .false., .false., .false., .false., .false., .false., .false., &
                                          .false., .false., .false., .false., .false., .false., .false., &
                                          .false., .false., .false., .false., .false., .false., .false., &
                                          .false., .false., .false., .false., .false., .false., .false., &
                                          .false., .false., .false., .false., .false., .false., .false., &
                                          .false., .false., .false., .false., .false., .false., .false., &
                                          .false., .false., .false., .false., .false., .false., .false., &
                                          .false., .false., .false., .false.,  .true., .false., .false., &
                                          .false., .false., .false., .false., .false., .false., .false., &
                                          .false., .false., .false., .false., .false., .false., .false., &
                                          .false., .false., .false., .false., .false., .false., .false., &
                                          .false., .false., .false., .false., .false., .false., .false., &
                                          .false., .false., .false., .false., .false., .false., .false., &
                                          .false., .false., .false., .false., .false., .false., .false., &
                                          .false., .false., .false., .false., .false., .false., .false., &
                                          .false., .false., .false.,  .true., .false., .false., .false., &
                                          .false., .false., .false., .false., .false., .false., .false., &
                                          .false., .false., .false., .false., .false., .false., .false., &
                                          .false., .false., .false., .false., .false., .false., .false., &
                                          .false., .false., .false., .false., .false., .false., .false., &
                                          .false., .false., .false., .false., .false., .false., .false., &
                                          .false., .false., .false., .false., .false., .false., .false., &
                                          .false., .false., .false.,  .true., .false., .false., .false., &
                                          .false., .false., .false., .false., .false., .false., .false., &
                                          .false., .false., .false., .false., .false., .false., .false., &
                                          .false., .false., .false., .false., .false., .false., .false., &
                                          .false., .false., .false., .false., .false., .false., .false., &
                                          .false., .false., .false., .false., .false., .false., .false., &
                                          .false., .false., .false., .false., .false., .false., .false., &
                                          .false., .false., .false., .false., .false., .false., .false., &
                                          .false., .false., .false., .false., .false., .false., .false., &
                                          .false., .false., .false., .false., .false., .false., .false., &
                                          .false., .false., .false., .false., .false., .false.,  .true. /), &
                                        (/ p%nodes_per_elem, p%nqp /) )
   baseline_QPrangeOverlapFE = reshape( (/ 1,   9,&
                                           9,  20,&
                                          20,  29,&
                                          29,  37,&
                                          37,  48,&
                                           0,   0/), &
                                        (/ 2, p%nodes_per_elem, p%elem_total /) )
   baseline_QPtWghtIntForceFE =reshape((/ 0.000000000000_BDKi, 0.003249999999_BDKi, 0.016250000000_BDKi, 0.016250000000_BDKi, 0.016300000000_BDKi, 0.016250000000_BDKi, 0.016250000000_BDKi, &
                                          0.016250000000_BDKi, 0.016250000000_BDKi, 0.000009815950_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.003249999999_BDKi, 0.016250000000_BDKi, 0.016250000000_BDKi, 0.016300000000_BDKi, 0.016250000000_BDKi, 0.016250000000_BDKi, 0.016250000000_BDKi, &
                                          0.016250000000_BDKi, 0.000790184049_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.016290184049_BDKi, 0.016250000000_BDKi, 0.016250000000_BDKi, 0.016250000000_BDKi, 0.016250000000_BDKi, &
                                          0.016300000000_BDKi, 0.016250000000_BDKi, 0.016250000000_BDKi, 0.016250000000_BDKi, 0.032550000000_BDKi, 0.032500000000_BDKi, 0.025748156682_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.015509815950_BDKi, 0.016250000000_BDKi, 0.016250000000_BDKi, 0.016250000000_BDKi, 0.016250000000_BDKi, 0.016300000000_BDKi, &
                                          0.016250000000_BDKi, 0.016250000000_BDKi, 0.016250000000_BDKi, 0.032550000000_BDKi, 0.032500000000_BDKi, 0.032151843317_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.006801843317_BDKi, &
                                          0.032500000000_BDKi, 0.032500000000_BDKi, 0.032550000000_BDKi, 0.032500000000_BDKi, 0.032550000000_BDKi, 0.032500000000_BDKi, 0.032500000000_BDKi, &
                                          0.032550000000_BDKi, 0.014157000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000398156682_BDKi, 0.032500000000_BDKi, &
                                          0.032500000000_BDKi, 0.032550000000_BDKi, 0.032500000000_BDKi, 0.032550000000_BDKi, 0.032500000000_BDKi, 0.032500000000_BDKi, 0.032550000000_BDKi, &
                                          0.028743000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.018343000000_BDKi, 0.032550000000_BDKi, 0.032500000000_BDKi, 0.032500000000_BDKi, 0.032550000000_BDKi, 0.032500000000_BDKi, &
                                          0.032550000000_BDKi, 0.032500000000_BDKi, 0.000096153846_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.003757000000_BDKi, 0.032550000000_BDKi, 0.032500000000_BDKi, 0.032500000000_BDKi, 0.032550000000_BDKi, 0.032500000000_BDKi, 0.032550000000_BDKi, &
                                          0.032500000000_BDKi, 0.002403846153_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.016153846153_BDKi, 0.016250000000_BDKi, 0.016300000000_BDKi, 0.008100000000_BDKi, 0.008150000000_BDKi, &
                                          0.008100000000_BDKi, 0.008149999999_BDKi, 0.008149999999_BDKi, 0.008100000000_BDKi, 0.008149999999_BDKi, 0.008100000000_BDKi, 0.004900000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.013846153846_BDKi, 0.016250000000_BDKi, 0.016300000000_BDKi, 0.008100000000_BDKi, 0.008150000000_BDKi, 0.008100000000_BDKi, &
                                          0.008149999999_BDKi, 0.008149999999_BDKi, 0.008100000000_BDKi, 0.008149999999_BDKi, 0.008100000000_BDKi, 0.004900000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, &
                                          0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi, 0.000000000000_BDKi  &
                                       /), (/ p%nqp, 2, p%nodes_per_elem, p%elem_total /))

   call BD_FEinternalForceQPweights(p, ErrStat, ErrMsg)

      ! check mapping flags
   do idx_qp = 1, p%nqp
      do idx_FE = 1, p%nodes_per_elem
         @assertEqual(baseline_FEoutboardOfQPt(idx_FE,idx_qp), p%FEoutboardOfQPt(idx_FE,idx_qp), testname)
      end do
   end do

      ! check QP mapping ranges
   do nelem = 1, p%elem_total
      do idx_FE = 1, p%nodes_per_elem
         do i=1,2
            @assertEqual(baseline_QPrangeOverlapFE(i,idx_FE,nelem), p%QPrangeOverlapFE(i,idx_FE,nelem), testname)
         end do
      end do
   end do

      ! check QPt weightings
   do nelem = 1, p%elem_total
      do idx_FE = 1, p%nodes_per_elem
         do i=1,2
            do idx_qp=1,p%nqp
               @assertEqual(baseline_QPtWghtIntForceFE(idx_qp,i,idx_FE,nelem), p%QPtWghtIntForceFE(idx_qp,i,idx_FE,nelem), tolerance, testname)
            enddo
         enddo
      enddo
   enddo

   call CleanUp

   ! --------------------------------------------------------------------------
contains
subroutine SetupArrays
      ! arrays given
   call AllocAry( p%GLL_Nodes, p%nodes_per_elem, 'GLL_Nodes in natural [-1,1] frame', ErrStat, ErrMsg)
   call AllocAry( p%QPtN, p%nqp,'p%QPtN', ErrStat, ErrMsg)
      ! arrays returned
   call AllocAry( p%FEoutboardOfQPt,             p%nodes_per_elem, p%nqp,        'p%FEoutboardOfQPt',                   ErrStat, ErrMsg )
   call AllocAry( p%QPrangeOverlapFE,         2, p%nodes_per_elem, p%elem_total, 'p%QPrangeOverlapFE  -- optimization', ErrStat, ErrMsg )
   call AllocAry( p%QPtWghtIntForceFE, p%nqp, 2, p%nodes_per_elem, p%elem_total, 'p%QPtWghtIntForceFE -- optimization', ErrStat, ErrMsg )
         ! baseline for arrays returned
      call AllocAry( baseline_FEoutboardOfQPt,             p%nodes_per_elem, p%nqp,        'baseline_FEoutboardOfQPt',                   ErrStat, ErrMsg )
      call AllocAry( baseline_QPrangeOverlapFE,         2, p%nodes_per_elem, p%elem_total, 'baseline_QPrangeOverlapFE  -- optimization', ErrStat, ErrMsg )
      call AllocAry( baseline_QPtWghtIntForceFE, p%nqp, 2, p%nodes_per_elem, p%elem_total, 'baseline_QPtWghtIntForceFE -- optimization', ErrStat, ErrMsg )
  
   end subroutine SetupArrays

   subroutine CleanUp
      deallocate(p%GLL_Nodes)
      deallocate(p%QPtN)
      deallocate(p%FEoutboardOfQPt)
      deallocate(p%QPrangeOverlapFE)
      deallocate(p%QPtWghtIntForceFE)
   end subroutine CleanUp

end subroutine
