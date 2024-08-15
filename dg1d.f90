program femDG
   implicit none 
     integer, parameter :: dp=selected_real_kind(32)
     integer :: nel, nen,ndof, np, nint, i, j, e, nm, npm
     integer :: n,l, ierror, ig, jg, no, kk
     
     real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)
     integer , dimension(:,:) , allocatable :: loc_glob
     
     real(dp), dimension(:,:),allocatable :: Kglob, Mglob
     real(dp), dimension(:,:),allocatable :: uloc
     
     real(dp), dimension(:),allocatable :: Fglob, coord, w, pos, det, deter1, sol_edp
     real(dp), dimension(:),allocatable :: Uexata, rf, xl
     real(dp), dimension(:),allocatable :: erro, numel
     
     real(dp), dimension(:,:,:),allocatable :: shl, shg, shlmt, shgmt
     
     real(dp) :: dx, xini, xfim, he, emax, erro_max
     real(dp) :: tt1,tt2
     real(dp) :: beta, beta0, eps, alfa, xx

     call cpu_time(tt1)
     kk=5
     allocate(erro(kk),numel(kk))
    
    do nm=1,kk
      nel       = 4*2**(nm-1)
      numel(nm) = nel 
      ndof = 1
      nen  = 2
      np   = nel*nen*ndof
      npm  = nel*(nen-1) + 1
      nint = nen + 1      

   allocate(Kglob(np,np),Fglob(np),loc_glob(nel,nen))
   allocate(coord(npm),w(nint),pos(nint), uloc(nel,nen), Uexata(np) )
   allocate(shl(2,nen,nint),shg(2,nen,nint), det(nint), xl(nen), deter1(nen), sol_edp(np))
   allocate(shlmt(2,nen,nen), shgmt(2,nen,nen), rf(nen))

      xini = 0.0_dp
      xfim = 1.0_dp    
      dx   = (xfim-xini)/(npm-1)

      beta0= 10.0_dp*(nen-1)**2 ! stabilization parameter
      beta = beta0/dx

      alfa = -1._dp  ! adjoint consistent parameter
      eps  = 1._dp  ! diffusive parameter

!==================================!
!== Vetor de coordenadas globais ==!
!==================================!
      do i=1,npm
            coord(i) = xini + (i-1)*dx
      enddo


!========================================! 
!== Constroi a matriz de conectividade ==!
!========================================!
       do e=1,nel
            do i=1,nen
                   loc_glob(e,i) = (e-1)*(nen-1) + i
            enddo
      enddo

!=============================================!
! ===   Monta a Matriz e Vetor Global  ====   !
!=============================================!
      Kglob = 0._dp
      Fglob = 0._dp


!======================!
! Funções forma locais !
!======================!
      call shl1d(shl,w,nint,nen,pos)
      call shl_mult(shlmt,nen,rf)


!================================================!
! Monta a Matriz de Rigidez                      !
!================================================!
      do n=1,nel
             do i=1,nen
                   no = loc_glob(n,i)
                   xl(i) = coord(no)
             enddo
!------------------------------------------------!
! Calcula as derivadas globais das funções forma !
!------------------------------------------------!
            call shglobal(xl,det,shl,shg,nint,nel,nen)
            call shglobal_mult(xl,deter1,shlmt,shgmt,nel,nen,rf)

        do l=1,nint
            xx = 0.0_dp
            do i = 1,nen
                xx = xx + xl(i)*shg(1,i,l)
            end do
            do i=1,nen
                ig = (n-1)*nen+i
                Fglob(ig) = Fglob(ig) + eps*pi*pi*sin(pi*xx)*shg(1,i,l)*w(l)*det(l)
            end do
        end do

        do i=1,nen
                ig = (n-1)*nen+i
                Uexata(ig) = sin(pi*xl(i))
            do j=1,nen
                jg = (n-1)*nen+j
                do l=1,nint
                    Kglob(ig,jg) = Kglob(ig,jg) + eps*shg(2,i,l)*shg(2,j,l)*w(l)*det(l)
                enddo
                if(n==1) then
                    !x0
                    Kglob(ig,jg) = Kglob(ig,jg) + eps*shgmt(2,j,1)*shgmt(1,i,1) &
                                                - alfa*eps*shgmt(1,j,1)*shgmt(2,i,1) &
                                                + beta*shgmt(1,j,1)*shgmt(1,i,1)
                    !x+
                    Kglob(ig,jg) = Kglob(ig,jg) - 0.5_dp*eps*shgmt(2,j,nen)*shgmt(1,i,nen) &
                                                + 0.5_dp*alfa*eps*shgmt(1,j,nen)*shgmt(2,i,nen) &
                                                + beta*shgmt(1,j,nen)*shgmt(1,i,nen)
                    !x+/x-
                    Kglob(ig,jg+nen) = Kglob(ig,jg+nen) - 0.5_dp*eps*shgmt(2,j,1)*shgmt(1,i,nen) &
                                                        - 0.5_dp*alfa*eps*shgmt(1,j,1)*shgmt(2,i,nen) &
                                                        - beta*shgmt(1,j,1)*shgmt(1,i,nen)

                else if(n==nel) then
                    !xN
                    Kglob(ig,jg) = Kglob(ig,jg) - eps*shgmt(2,j,nen)*shgmt(1,i,nen) &
                                                + alfa*eps*shgmt(1,j,nen)*shgmt(2,i,nen) &
                                                + beta*shgmt(1,j,nen)*shgmt(1,i,nen)
                    !x-
                    Kglob(ig,jg) = Kglob(ig,jg) + 0.5_dp*eps*shgmt(2,j,1)*shgmt(1,i,1) &
                                                - 0.5_dp*alfa*eps*shgmt(1,j,1)*shgmt(2,i,1) &
                                                + beta*shgmt(1,j,1)*shgmt(1,i,1)
                    !x-/x+
                    Kglob(ig,jg-nen) = Kglob(ig,jg-nen) + 0.5_dp*eps*shgmt(2,j,nen)*shgmt(1,i,1) &
                                                        + 0.5_dp*alfa*eps*shgmt(1,j,nen)*shgmt(2,i,1) &
                                                        - beta*shgmt(1,j,nen)*shgmt(1,i,1)

                else
                    !x+
                    Kglob(ig,jg) = Kglob(ig,jg) - 0.5_dp*eps*shgmt(2,j,nen)*shgmt(1,i,nen) &
                                                + 0.5_dp*alfa*eps*shgmt(1,j,nen)*shgmt(2,i,nen) &
                                                + beta*shgmt(1,j,nen)*shgmt(1,i,nen)

                    !x-
                    Kglob(ig,jg) = Kglob(ig,jg) + 0.5_dp*eps*shgmt(2,j,1)*shgmt(1,i,1) &
                                                - 0.5_dp*alfa*eps*shgmt(1,j,1)*shgmt(2,i,1) &
                                                + beta*shgmt(1,j,1)*shgmt(1,i,1)
                    !x+/x-
                    Kglob(ig,jg+nen) = Kglob(ig,jg+nen) - 0.5_dp*eps*shgmt(2,j,1)*shgmt(1,i,nen) &
                                                        - 0.5_dp*alfa*eps*shgmt(1,j,1)*shgmt(2,i,nen) &
                                                        - beta*shgmt(1,j,1)*shgmt(1,i,nen)

                    !x-/x+
                    Kglob(ig,jg-nen) = Kglob(ig,jg-nen) + 0.5_dp*eps*shgmt(2,j,nen)*shgmt(1,i,1) &
                                                        + 0.5_dp*alfa*eps*shgmt(1,j,nen)*shgmt(2,i,1) &
                                                        - beta*shgmt(1,j,nen)*shgmt(1,i,1)

                end if
            enddo

        enddo

    enddo ! fim do loop dos elementos
!write(*,"(8f16.8)") ((Kglob(i,j), j=1,np), i=1,np)

!======================================!
! Impoe as cond de contorno Dirichlet  !
!======================================!
      Kglob(:,1) = 0._dp;      Kglob(1,:) = 0._dp
      Kglob(:,np) = 0._dp;  Kglob(np,:) = 0._dp
      Kglob(1,1) = 1._dp;      Kglob(np,np) = 1._dp

  !  write(*,"(8f16.8)") ((Kglob(i,j), j=1,np), i=1,np)
  !  write(*,"(1f16.8)") (Fglob(j), j=1,np)

        !======================================!
        ! Impoe as cond de contorno Dirichlet  !
        !======================================!

        Fglob(1)  = 0.0_dp
        Fglob(np) = 0.0_dp


        !=============================!
        ! Resolve o sistema linear    !
        !=============================!
				call solve2(np,Kglob,Fglob,sol_edp,ierror)
   !   do n=1,nel
   !          do i=1,nen
   !             ig = (n-1)*nen+i
   !             uloc(n,i) = sol_edp(ig)
   !         enddo
   !   enddo

    !   do i=1,np
     !       write(*,*)Uexata(i), sol_edp(i)
      !  enddo
        !*****************************************!
        !    Deixa a solucao na forma continua    !
        !*****************************************!

  !        if (nen .eq. 2)then
  !              do n=1,nel-1
  !                      sol_edp(loc_glob(n,nen))=(uloc(n+1,1)+uloc(n,nen))/2.0_dp
  !              enddo
  !        else
  !            do n=1,nel-1
  !                     sol_edp(loc_glob(n,nen))=(uloc(n+1,1)+uloc(n,nen))/2.0_dp
  !            enddo
  !            do n=1,nel
  !                do i=2,nen-1
  !                     sol_edp(loc_glob(n,i)) = uloc(n,i)
  !                enddo
  !            enddo
  !        endif
  !          sol_edp(loc_glob(1,1)) = uloc(1,1)
  !          sol_edp(loc_glob(nel,nen)) = uloc(nel,nen)

            !=========================!
            !Calculo do Erro do maximo!
            !=========================!

            erro_max = 0._dp
            do i = 1,np
                  emax = Uexata(i) - sol_edp(i)
                  if(abs(emax) .gt. erro_max)then
                        erro_max = abs(emax)
                  endif
            enddo
            erro(nm)=erro_max

            !===================!
            !Calculo do Erro L2 !
            !===================!

            deallocate(Kglob,Fglob, loc_glob)
            deallocate(coord, w, pos, Uexata)
            deallocate(shl,shg, det, xl, deter1)
            deallocate(shlmt, shgmt, rf)
            deallocate(uloc,sol_edp)
  enddo


        do i=1,kk-1
            write(*,*) -(log(erro(i+1)) - log(erro(i)))/(log(real(numel(i+1)))-log(real(numel(i))))
        end do
  call cpu_time(tt2)
       write(*,*) "Tempo de Execucao =", abs(tt2-tt1)

            
 

end program



