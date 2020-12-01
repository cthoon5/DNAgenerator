program GenDNAcoordinate
!
!
  implicit none
  integer,parameter   :: nin=19, nin2=20, nout=21, maxmer=100
  integer             :: i, j, k
  character*80        :: line, arg
  character*3         :: SuC(maxmer),LinkC(maxmer),BaseC(maxmer)
  character*4         :: CheckC

  real*8  Thy(14,3), Ade(14,3), Gua(15,3), Cyt(12,3), mCy(15,3)
  real*8  SugarI(14,3), Sugar(14,3), Terminal(3), Linkage(4,3),LNA(16,3)
  integer nAtom, Is, iPO3, iP,nmer
  integer iO2(maxmer), iC2(maxmer), iC3(maxmer), iC5(maxmer), iO1(maxmer)
!
  common /coord_data/ Thy, Ade, Gua, Cyt, mCy,SugarI,Sugar,Terminal,Linkage,LNA
!

  call getdata
  
  call getarg(1,line) 
!
  open(nin,file=trim(line))
  open(nout,file="zmatrix.dat")
! 
  read(nin,'(a80)',end=1000) line 
  read(line,'(I6)',end=1000) nmer
  Do k=1, 3
    read(nin,'(a80)') line 
    read(line,'(1x,A4)') CheckC
    If (CheckC .eq. "Suga" ) then
     Do i=1, nmer
      read(nin,'(a80)') line 
      read(line,'(A3)') SuC(i)
     enddo
    else if (CheckC .eq. "Base" ) then
     Do i=1, nmer
      read(nin,'(a80)') line 
      read(line,'(A3)') BaseC(i)
     enddo
    else if (CheckC .eq. "Link" ) then
     Do i=1, nmer
      read(nin,'(a80)') line 
      read(line,'(A3)') LinkC(i)
     enddo
    endif
  Enddo 

   write(nout,'(A27)')"#ZMATRIX                   "
   write(nout,'(A27)')"#                          "
   write(nout,'(A27)')"H                          "
   write(nout,'(A27)')"O       1        0.95935   "
  nAtom=2 
  Do i=1, nmer
   if (SuC(i) .eq. "DNA") then
      if (i .eq. 1 ) then
        iPO3 =2; iP=1; Is=3;
        write(nout,10) "C", iPO3, Sugar(1,1),   iP, Sugar(1,2)
      else
        Is=nAtom+1
        write(nout,11) "C", iPO3, Sugar(1,1),   iP, Sugar(1,2),  iO2(i-1), Sugar(1,3)
      endif
      write(nout,11) "C", Is+0, Sugar(2,1), iPO3, Sugar(2,2),   iP, Sugar(2,3)
      write(nout,11) "O", Is+1, Sugar(3,1), Is+0, Sugar(3,2), iPO3, Sugar(3,3)
      write(nout,11) "C", Is+1, Sugar(4,1), Is+0, Sugar(4,2), iPO3, Sugar(4,3)
      write(nout,11) "O", Is+3, Sugar(5,1), Is+1, Sugar(5,2), Is+0, Sugar(5,3)
      write(nout,11) "C", Is+3, Sugar(6,1), Is+1, Sugar(6,2), Is+0, Sugar(6,3)
      write(nout,11) "C", Is+2, Sugar(7,1), Is+1, Sugar(7,2), Is+0, Sugar(7,3)
      write(nout,11) "H", Is+0, Sugar(8,1), iPO3, Sugar(8,2),   iP , Sugar(8,3)
      write(nout,11) "H", Is+0, Sugar(9,1), iPO3, Sugar(9,2),   iP, Sugar(9,3)
      write(nout,11) "H", Is+1, Sugar(10,1),Is+0, Sugar(10,2),iPO3, Sugar(10,3)
      write(nout,11) "H", Is+6, Sugar(11,1),Is+2, Sugar(11,2),Is+1, Sugar(11,3)
      write(nout,11) "H", Is+5, Sugar(12,1),Is+3, Sugar(12,2),Is+1, Sugar(12,3)
      write(nout,11) "H", Is+5, Sugar(13,1),Is+3, Sugar(13,2),Is+1, Sugar(13,3)
      write(nout,11) "H", Is+3, Sugar(14,1),Is+1, Sugar(14,2),Is+0, Sugar(14,3)
      nAtom=nAtom+14
      iO2(i)=Is+4; iC2(i)=Is+1; iC3(i)=Is+3; iC5(i)=Is+6; iO1(i)=Is+2; 
   else if (SuC(i) .eq. "LNA") then
      if (i .eq. 1 ) then
        iPO3 =2; iP=1; Is=3;
        write(nout,10) "C", iPO3, LNA(1,1),   iP, LNA(1,2)
      else
        Is=nAtom+1
        write(nout,11) "C", iPO3, LNA(1,1),   iP, LNA(1,2),  iO2(i-1), LNA(1,3)
      endif
      write(nout,11) "C", Is+0, LNA(2,1), iPO3, LNA(2,2),   iP, LNA(2,3)
      write(nout,11) "O", Is+1, LNA(3,1), Is+0, LNA(3,2), iPO3, LNA(3,3)
      write(nout,11) "C", Is+1, LNA(4,1), Is+0, LNA(4,2), Is+2, LNA(4,3)
      write(nout,11) "O", Is+3, LNA(5,1), Is+1, LNA(5,2), Is+0, LNA(5,3)
      write(nout,11) "C", Is+3, LNA(6,1), Is+1, LNA(6,2), Is+0, LNA(6,3)
      write(nout,11) "C", Is+2, LNA(7,1), Is+1, LNA(7,2), Is+0, LNA(7,3)
      write(nout,11) "C", Is+1, LNA(8,1), IS+0, LNA(8,2), Is+2, LNA(8,3)
      write(nout,11) "O", Is+5, LNA(9,1), Is+3, LNA(9,2), Is+1, LNA(9,3)
      write(nout,11) "H", Is+0, LNA(10,1),Is+1, LNA(10,2),Is+2, LNA(10,3)
      write(nout,11) "H", Is+0, LNA(11,1),Is+1, LNA(11,2),Is+2, LNA(11,3)
      write(nout,11) "H", Is+6, LNA(12,1),Is+2, LNA(12,2),Is+1, LNA(12,3)
      write(nout,11) "H", Is+5, LNA(13,1),Is+3, LNA(13,2),Is+1, LNA(13,3)
      write(nout,11) "H", Is+3, LNA(14,1),Is+1, LNA(14,2),Is+0, LNA(14,3)
      write(nout,11) "H", Is+7, LNA(15,1),Is+1, LNA(15,2),Is+0, LNA(15,3)
      write(nout,11) "H", Is+7, LNA(16,1),Is+1, LNA(16,2),Is+0, LNA(16,3)
      nAtom=nAtom+16
      iO2(i)=Is+4; iC2(i)=Is+1; iC3(i)=Is+3; iC5(i)=Is+6; iO1(i)=Is+2; 
   endif
   if (i .eq. nmer) then
     write(nout,11) "H", Is+4,Terminal(1), Is+3, Terminal(2), Is+1,Terminal(3)
     nAtom=nAtom+1
   else
     Is=nAtom+1
     write(nout,11) "P",iO2(i),Linkage(1,1),iC3(i), Linkage(1,2), iC2(i), Linkage(1,3)
     write(nout,11) "S", Is+0, Linkage(2,1), iO2(i), Linkage(2,2), iC3(i), Linkage(2,3)
     write(nout,11) "O", Is+0, Linkage(3,1), iO2(i), Linkage(3,2), iC3(i), Linkage(3,3)
     write(nout,11) "O", Is+0, Linkage(4,1), iO2(i), Linkage(4,2), iC3(i), Linkage(4,3)
     nAtom=nAtom+4
     iPO3=Is+3; iP=Is;
   endif 
  enddo

  Do i=1,nmer 
   if (BaseC(i) .eq. "Thy") then
     Is=nAtom+1
     write(nout,11) "N",iC5(i), Thy(1,1),iO1(i), Thy(1,2),iC2(i), Thy(1,3)
     write(nout,11) "C",  Is+0, Thy(2,1),iC5(i), Thy(2,2),IO1(i), Thy(2,3)
     write(nout,11) "O",  Is+1, Thy(3,1),  Is+0, Thy(3,2),iC5(i), Thy(3,3)
     write(nout,11) "N",  Is+1, Thy(4,1),  Is+0, Thy(4,2),iC5(i), Thy(4,3)
     write(nout,11) "C",  Is+3, Thy(5,1),  Is+1, Thy(5,2), Is+0, Thy(5,3)
     write(nout,11) "O",  Is+4, Thy(6,1),  Is+3, Thy(6,2), Is+1, Thy(6,3)
     write(nout,11) "C",  Is+4, Thy(7,1),  Is+3, Thy(7,2), Is+1, Thy(7,3)
     write(nout,11) "C",  Is+6, Thy(8,1),  Is+4, Thy(8,2), Is+3, Thy(8,3)
     write(nout,11) "C",  Is+6, Thy(9,1),  Is+4, Thy(9,2), Is+3, Thy(9,3)
     write(nout,11) "H",  Is+8, Thy(10,1), Is+6, Thy(10,2),Is+4, Thy(10,3)
     write(nout,11) "H",  Is+7, Thy(11,1), Is+6, Thy(11,2),Is+4, Thy(11,3)
     write(nout,11) "H",  Is+7, Thy(12,1), Is+6, Thy(12,2),Is+4, Thy(12,3)
     write(nout,11) "H",  Is+7, Thy(13,1), Is+6, Thy(13,2),Is+4, Thy(13,3)
     write(nout,11) "H",  Is+3, Thy(14,1) ,Is+1, Thy(14,2),Is+0, Thy(14,3)
     nAtom=nAtom+14
   else if (BaseC(i) .eq. "Ade") then
     Is=nAtom+1
     write(nout,11) "N",iC5(i), Ade(1,1),iO1(i), Ade(1,2),iC2(i), Ade(1,3)
     write(nout,11) "C",  Is+0, Ade(2,1),iC5(i), Ade(2,2),iO1(i), Ade(2,3)
     write(nout,11) "N",  Is+1, Ade(3,1),  Is+0, Ade(3,2),iC5(i), Ade(3,3)
     write(nout,11) "C",  Is+2, Ade(4,1),  Is+1, Ade(4,2), Is+0, Ade(4,3)
     write(nout,11) "C",  Is+3, Ade(5,1),  Is+2, Ade(5,2), Is+1, Ade(5,3)
     write(nout,11) "N",  Is+4, Ade(6,1),  Is+3, Ade(6,2), Is+2, Ade(6,3)
     write(nout,11) "N",  Is+4, Ade(7,1),  Is+3, Ade(7,2), Is+2, Ade(7,3)
     write(nout,11) "C",  Is+6, Ade(8,1),  Is+4, Ade(8,2), Is+3, Ade(8,3)
     write(nout,11) "N",  Is+7, Ade(9,1),  Is+6, Ade(9,2), Is+4, Ade(9,3)
     write(nout,11) "C",  Is+8, Ade(10,1), Is+7, Ade(10,2),Is+6, Ade(10,3)
     write(nout,11) "H",  Is+1, Ade(11,1), Is+0, Ade(11,2),Is+9, Ade(11,3)
     write(nout,11) "H",  Is+5, Ade(12,1), Is+4, Ade(12,2),Is+3, Ade(12,3)
     write(nout,11) "H",  Is+5, Ade(13,1), Is+4, Ade(13,2),Is+3, Ade(13,3)
     write(nout,11) "H",  Is+7, Ade(14,1) ,Is+6, Ade(14,2),Is+4, Ade(14,3)
     nAtom=nAtom+14
   else if (BaseC(i) .eq. "Gua") then
     Is=nAtom+1
     write(nout,11) "N",iC5(i), Gua(1,1),iO1(i), Gua(1,2),iC2(i), Gua(1,3)
     write(nout,11) "C",  Is+0, Gua(2,1),iC5(i), Gua(2,2),iO1(i), Gua(2,3)
     write(nout,11) "N",  Is+1, Gua(3,1),  Is+0, Gua(3,2),iC5(i), Gua(3,3)
     write(nout,11) "C",  Is+2, Gua(4,1),  Is+1, Gua(4,2), Is+0, Gua(4,3)
     write(nout,11) "C",  Is+3, Gua(5,1),  Is+2, Gua(5,2), Is+1, Gua(5,3)
     write(nout,11) "O",  Is+4, Gua(6,1),  Is+3, Gua(6,2), Is+2, Gua(6,3)
     write(nout,11) "N",  Is+4, Gua(7,1),  Is+3, Gua(7,2), Is+2, Gua(7,3)
     write(nout,11) "C",  Is+6, Gua(8,1),  Is+4, Gua(8,2), Is+3, Gua(8,3)
     write(nout,11) "N",  Is+7, Gua(9,1),  Is+6, Gua(9,2), Is+4, Gua(9,3)
     write(nout,11) "N",  Is+7, Gua(10,1), Is+6, Gua(10,2),Is+4, Gua(10,3)
     write(nout,11) "C",  Is+9, Gua(11,1), Is+7, Gua(11,2),Is+6, Gua(11,3)
     write(nout,11) "H",  Is+1, Gua(12,1), Is+0, Gua(12,2),Is+10, Gua(12,3)
     write(nout,11) "H",  Is+6, Gua(13,1), Is+4, Gua(13,2),Is+3, Gua(13,3)
     write(nout,11) "H",  Is+8, Gua(14,1) ,Is+7, Gua(14,2),Is+6, Gua(14,3)
     write(nout,11) "H",  Is+8, Gua(15,1) ,Is+7, Gua(15,2),Is+6, Gua(15,3)
     nAtom=nAtom+15
   else if (BaseC(i) .eq. "Cyt") then
     Is=nAtom+1
     write(nout,11) "N",iC5(i), Cyt(1,1),iO1(i), Cyt(1,2),iC2(i), Cyt(1,3)
     write(nout,11) "C",  Is+0, Cyt(2,1),iC5(i), Cyt(2,2),iO1(i), Cyt(2,3)
     write(nout,11) "O",  Is+1, Cyt(3,1),  Is+0, Cyt(3,2),iC5(i), Cyt(3,3)
     write(nout,11) "N",  Is+1, Cyt(4,1),  Is+0, Cyt(4,2),iC5(i), Cyt(4,3)
     write(nout,11) "C",  Is+3, Cyt(5,1),  Is+1, Cyt(5,2), Is+0, Cyt(5,3)
     write(nout,11) "N",  Is+4, Cyt(6,1),  Is+3, Cyt(6,2), Is+1, Cyt(6,3)
     write(nout,11) "C",  Is+4, Cyt(7,1),  Is+3, Cyt(7,2), Is+1, Cyt(7,3)
     write(nout,11) "C",  Is+6, Cyt(8,1),  Is+4, Cyt(8,2), Is+3, Cyt(8,3)
     write(nout,11) "H",  Is+7, Cyt(9,1),  Is+6, Cyt(9,2), Is+4, Cyt(9,3)
     write(nout,11) "H",  Is+6, Cyt(10,1), Is+4, Cyt(10,2),Is+3, Cyt(10,3)
     write(nout,11) "H",  Is+5, Cyt(11,1) ,Is+4, Cyt(11,2),Is+3, Cyt(11,3)
     write(nout,11) "H",  Is+5, Cyt(12,1) ,Is+4, Cyt(12,2),Is+3, Cyt(12,3)
     nAtom=nAtom+12
   else if (BaseC(i) .eq. "mCy") then
     Is=nAtom+1
     write(nout,11) "N",iC5(i), mCy(1,1),iO1(i), mCy(1,2),iC2(i), mCy(1,3)
     write(nout,11) "C",  Is+0, mCy(2,1),iC5(i), mCy(2,2),iO1(i), mCy(2,3)
     write(nout,11) "O",  Is+1, mCy(3,1),  Is+0, mCy(3,2),iC5(i), mCy(3,3)
     write(nout,11) "N",  Is+1, mCy(4,1),  Is+0, mCy(4,2),iC5(i), mCy(4,3)
     write(nout,11) "C",  Is+3, mCy(5,1),  Is+1, mCy(5,2), Is+0, mCy(5,3)
     write(nout,11) "N",  Is+4, mCy(6,1),  Is+3, mCy(6,2), Is+1, mCy(6,3)
     write(nout,11) "C",  Is+4, mCy(7,1),  Is+3, mCy(7,2), Is+1, mCy(7,3)
     write(nout,11) "C",  Is+6, mCy(8,1),  Is+4, mCy(8,2), Is+3, mCy(8,3)
     write(nout,11) "C",  Is+6, mCy(9,1),  Is+4, mCy(9,2), Is+3, mCy(9,3)
     write(nout,11) "H",  Is+8, mCy(10,1), Is+6, mCy(10,2),Is+4, mCy(10,3)
     write(nout,11) "H",  Is+7, mCy(11,1) ,Is+6, mCy(11,2),Is+4, mCy(11,3)
     write(nout,11) "H",  Is+7, mCy(12,1) ,Is+6, mCy(12,2),Is+4, mCy(12,3)
     write(nout,11) "H",  Is+7, mCy(13,1) ,Is+6, mCy(13,2),Is+4, mCy(13,3)
     write(nout,11) "H",  Is+5, mCy(14,1) ,Is+4, mCy(14,2),Is+3, mCy(14,3)
     write(nout,11) "H",  Is+5, mCy(15,1) ,Is+4, mCy(15,2),Is+3, mCy(15,3)
     nAtom=nAtom+15
   endif
  enddo

1000 continue
10 format(A1,2x,I6,5x,F10.5,I6,5x,F10.5)
11 format(A1,2x,I6,5x,F10.5,I6,5x,F10.5,I6,5x,F10.5)
2 format(3F20.10)
end program GenDNAcoordinate
!
SUBROUTINE getdata

  real*8  Thy(14,3), Ade(14,3), Gua(15,3), Cyt(12,3), mCy(15,3)
  real*8  SugarI(14,3), Sugar(14,3), Terminal(3), Linkage(4,3), LNA(16,3)
!
  common /coord_data/Thy,Ade,Gua,Cyt,mCy,SugarI,Sugar,Terminal,Linkage,LNA
!

Thy(1,1) = 1.48688 ; Thy(1,2) = 107.84238 ; Thy(1,3) = 232.31857 ;
Thy(2,1) = 1.37795 ; Thy(2,2) = 117.3744 ; Thy(2,3) = 262.14231 ;
Thy(3,1) = 1.21728 ; Thy(3,2) = 122.71194 ; Thy(3,3) = 0.03124 ;
Thy(4,1) = 1.38184 ; Thy(4,2) = 115.25611 ; Thy(4,3) = 179.9509 ;
Thy(5,1) = 1.3733 ; Thy(5,2) = 126.63024 ; Thy(5,3) = 0.20298 ; 
Thy(6,1) = 1.23702 ; Thy(6,2) = 120.80598 ; Thy(6,3) = 179.92276 ;
Thy(7,1) = 1.44265 ; Thy(7,2) = 114.29504 ; Thy(7,3) = 0.01731 ;
Thy(8,1) = 1.49828 ; Thy(8,2) = 116.1303 ; Thy(8,3) = 180.33442 ;
Thy(9,1) = 1.34855 ; Thy(9,2) = 120.76651 ; Thy(9,3) = 359.64172 ;
Thy(10,1) = 1.08982 ; Thy(10,2) = 119.49581 ; Thy(10,3) = 180.4412 ;
Thy(11,1) = 1.10075 ; Thy(11,2) = 109.57018 ; Thy(11,3) = 258.95033 ;
Thy(12,1) = 1.09885 ; Thy(12,2) = 109.49982 ; Thy(12,3) = 138.751 ;
Thy(13,1) = 1.10086 ; Thy(13,2) = 109.4565 ; Thy(13,3) = 18.74006 ;
Thy(14,1) = 1.00006 ; Thy(14,2) = 116.5931 ; Thy(14,3) = 180.15357 ;

Ade(1,1) = 1.4861 ; Ade(1,2) = 107.8767 ; Ade(1,3) = 232.37137 ;
Ade(2,1) = 1.37248 ; Ade(2,2) = 128.19626 ; Ade(2,3) = 82.37625 ;
Ade(3,1) = 1.29827 ; Ade(3,2) = 113.54738 ; Ade(3,3) = 180.25672 ;
Ade(4,1) = 1.38258 ; Ade(4,2) = 103.87701 ; Ade(4,3) = 359.3457 ;
Ade(5,1) = 1.41201 ; Ade(5,2) = 132.35294 ; Ade(5,3) = 179.82368 ;
Ade(6,1) = 1.32936 ; Ade(6,2) = 123.73525 ; Ade(6,3) = 0.60177 ;
Ade(7,1) = 1.34178 ; Ade(7,2) = 117.40815 ; Ade(7,3) = 180.56925 ;
Ade(8,1) = 1.32924 ; Ade(8,2) = 118.90405 ; Ade(8,3) = 0.03672 ;
Ade(9,1) = 1.31548 ; Ade(9,2) = 129.50401 ; Ade(9,3) = 0.32922 ;
Ade(10,1) = 1.3457 ; Ade(10,2) = 110.21199 ; Ade(10,3) = 359.76836 ;
Ade(11,1) = 1.08942 ; Ade(11,2) = 123.17791 ; Ade(11,3) = 180.29093 ;
Ade(12,1) = 1.02937 ; Ade(12,2) = 120.04243 ; Ade(12,3) = 0.01924 ;
Ade(13,1) = 1.03059 ; Ade(13,2) = 119.97759 ; Ade(13,3) = 179.94974 ;
Ade(14,1) = 1.08393 ; Ade(14,2) = 115.22735 ; Ade(14,3) = 180.36626 ;

Gua(1,1) = 1.4861 ; Gua(1,2) = 107.84262 ; Gua(1,3) = 232.34404 ;
Gua(2,1) = 1.38038 ; Gua(2,2) = 129.20352 ; Gua(2,3) = 82.30471 ;
Gua(3,1) = 1.30721 ; Gua(3,2) = 114.10223 ; Gua(3,3) = 179.86859 ;
Gua(4,1) = 1.39533 ; Gua(4,2) = 104.04606 ; Gua(4,3) = 0.15864 ;
Gua(5,1) = 1.41385 ; Gua(5,2) = 130.2538 ; Gua(5,3) = 179.85847 ;
Gua(6,1) = 1.22987 ; Gua(6,2) = 129.11388 ; Gua(6,3) = 359.98939 ;
Gua(7,1) = 1.40682 ; Gua(7,2) = 111.43082 ; Gua(7,3) = 179.88905 ;
Gua(8,1) = 1.37956 ; Gua(8,2) = 125.04687 ; Gua(8,3) = 359.68508 ;
Gua(9,1) = 1.33499 ; Gua(9,2) = 116.20837 ; Gua(9,3) = 180.25052 ;
Gua(10,1) = 1.33061 ; Gua(10,2) = 123.60091 ; Gua(10,3) = 0.75571 ;
Gua(11,1) = 1.35997 ; Gua(11,2) = 112.31204 ; Gua(11,3) = 359.44467 ;
Gua(12,1) = 1.08979 ; Gua(12,2) = 122.92407 ; Gua(12,3) = 179.89225 ;
Gua(13,1) = 1.00012 ; Gua(13,2) = 117.47555 ; Gua(13,3) = 179.72095 ;
Gua(14,1) = 1.03031 ; Gua(14,2) = 120.11494 ; Gua(14,3) = 180.49281 ;
Gua(15,1) = 1.03151 ; Gua(15,2) = 120.05082 ; Gua(15,3) = 0.5101 ;

Cyt(1,1) = 1.48688 ; Cyt(1,2) = 107.84238 ; Cyt(1,3) = 232.31857 ;
Cyt(2,1) = 1.39104 ; Cyt(2,2) = 118.05937 ; Cyt(2,3) = 262.19427 ;
Cyt(3,1) = 1.23975 ; Cyt(3,2) = 118.87896 ; Cyt(3,3) = 359.76613 ;
Cyt(4,1) = 1.51013 ; Cyt(4,2) = 113.7231 ; Cyt(4,3) = 180.46983 ;
Cyt(5,1) = 1.14692 ; Cyt(5,2) = 123.01232 ; Cyt(5,3) = 358.71583 ;
Cyt(6,1) = 1.32591 ; Cyt(6,2) = 115.43867 ; Cyt(6,3) = 180.74506 ;
Cyt(7,1) = 1.4301 ; Cyt(7,2) = 124.43426 ; Cyt(7,3) = 1.65854 ; 
Cyt(8,1) = 1.35977 ; Cyt(8,2) = 116.69092 ; Cyt(8,3) = 358.71904 ;
Cyt(9,1) = 1.0884 ; Cyt(9,2) = 119.48642 ; Cyt(9,3) = 180.52263 ;
Cyt(10,1) = 1.08906 ; Cyt(10,2) = 121.65677 ; Cyt(10,3) = 178.67667 ;
Cyt(11,1) = 1.03054 ; Cyt(11,2) = 119.94855 ; Cyt(11,3) = 0.23466 ;
Cyt(12,1) = 1.02935 ; Cyt(12,2) = 119.97346 ; Cyt(12,3) = 180.23909 ;

mCy(1,1) = 1.48688 ; mCy(1,2) = 107.84238 ; mCy(1,3) = 232.31857 ;
mCy(2,1) = 1.39104 ; mCy(2,2) = 118.05937 ; mCy(2,3) = 262.19427 ;
mCy(3,1) = 1.23975 ; mCy(3,2) = 118.87896 ; mCy(3,3) = 359.76613 ;
mCy(4,1) = 1.51013 ; mCy(4,2) = 113.7231 ; mCy(4,3) = 180.46983 ;
mCy(5,1) = 1.14692 ; mCy(5,2) = 123.01232 ; mCy(5,3) = 358.71583 ;
mCy(6,1) = 1.32591 ; mCy(6,2) = 115.43867 ; mCy(6,3) = 180.74506 ;
mCy(7,1) = 1.4301 ; mCy(7,2) = 124.43426 ; mCy(7,3) = 1.65854 ; 
mCy(8,1) = 1.49828 ; mCy(8,2) = 116.1303 ; mCy(8,3) = 180.33442 ;
mCy(9,1) = 1.35977 ; mCy(9,2) = 116.69092 ; mCy(9,3) = 358.71904 ;
mCy(10,1) = 1.0884 ; mCy(10,2) = 119.48642 ; mCy(10,3) = 180.52263 ;
mCy(11,1) = 1.10075 ; mCy(11,2) = 109.57018 ; mCy(11,3) = 258.95033 ;
mCy(12,1) = 1.09885 ; mCy(12,2) = 109.49982 ; mCy(12,3) = 138.751 ;
mCy(13,1) = 1.10086 ; mCy(13,2) = 109.4565 ; mCy(13,3) = 18.74006 ;
mCy(14,1) = 1.03054 ; mCy(14,2) = 119.94855 ; mCy(14,3) = 0.23466 ;
mCy(15,1) = 1.02935 ; mCy(15,2) = 119.97346 ; mCy(15,3) = 180.23909 ;

Sugar(1,1) = 1.44902 ; Sugar(1,2) = 118.7766 ; Sugar(1,3) = 327.86658 ;
Sugar(2,1) = 1.51335 ; Sugar(2,2) = 109.7751 ; Sugar(2,3) = 213.86141 ;
Sugar(3,1) = 1.46168 ; Sugar(3,2) = 108.67189 ; Sugar(3,3) = 279.20556 ;
Sugar(4,1) = 1.52752 ; Sugar(4,2) = 116.30496 ; Sugar(4,3) = 36.4485 ;
Sugar(5,1) = 1.42212 ; Sugar(5,2) = 112.2627 ; Sugar(5,3) = 156.35578 ;
Sugar(6,1) = 1.522 ; Sugar(6,2) = 102.8081 ; Sugar(6,3) = 273.68762 ;
Sugar(7,1) = 1.41517 ; Sugar(7,2) = 109.79169 ; Sugar(7,3) = 106.20667 ;
Sugar(8,1) = 1.09938 ; Sugar(8,2) = 109.31868 ; Sugar(8,3) = 333.94883 ;
Sugar(9,1) = 1.09925 ; Sugar(9,2) = 109.45423 ; Sugar(9,3) = 93.84885 ;
Sugar(10,1) = 1.09638 ; Sugar(10,2) = 107.88402 ; Sugar(10,3) = 159.06718 ;
Sugar(11,1) = 1.09982 ; Sugar(11,2) = 109.77529 ; Sugar(11,3) = 114.07815 ;
Sugar(12,1) = 1.10068 ; Sugar(12,2) = 111.43212 ; Sugar(12,3) = 83.93756 ;
Sugar(13,1) = 1.10089 ; Sugar(13,2) = 111.1722 ; Sugar(13,3) = 206.28999 ;
Sugar(14,1) = 1.10125 ; Sugar(14,2) = 111.01739 ; Sugar(14,3) = 32.97188 ;

SugarI(1,1) = 1.44896 ; SugarI(1,2) = 122.26604 ; SugarI(1,3) = 0.0000 ;
SugarI(2,1) = 1.51306 ; SugarI(2,2) = 109.7952 ; SugarI(2,3) = 207.96653 ; 
SugarI(3,1) = 1.46261 ; SugarI(3,2) = 108.58236 ; SugarI(3,3) = 279.09613 ; 
SugarI(4,1) = 1.52712 ; SugarI(4,2) = 116.35105 ; SugarI(4,3) = 36.26968 ;
SugarI(5,1) = 1.42244 ; SugarI(5,2) = 112.26805 ; SugarI(5,3) = 156.39098 ; 
SugarI(6,1) = 1.52123 ; SugarI(6,2) = 102.8632 ; SugarI(6,3) = 273.78762 ; 
SugarI(7,1) = 1.41491 ; SugarI(7,2) = 109.75327 ; SugarI(7,3) = 106.23302 ; 
SugarI(8,1) = 1.10116 ; SugarI(8,2) = 109.28538 ; SugarI(8,3) = 327.99751 ; 
SugarI(9,1) = 1.09927 ; SugarI(9,2) = 109.47937 ; SugarI(9,3) = 87.8973 ; 
SugarI(10,1) = 1.09583 ; SugarI(10,2) = 107.92177 ; SugarI(10,3) = 159.02499 ; 
SugarI(11,1) = 1.09977 ; SugarI(11,2) = 109.80514 ; SugarI(11,3) = 114.17688 ; 
SugarI(12,1) = 1.10066 ; SugarI(12,2) = 111.40749 ; SugarI(12,3) = 83.87241 ;
SugarI(13,1) = 1.10078 ; SugarI(13,2) = 111.23286 ; SugarI(13,3) = 206.31246 ; 
SugarI(14,1) = 1.1012 ; SugarI(14,2) = 110.99405 ; SugarI(14,3) = 33.07257 ;


LNA(1,1) = 1.44902 ; LNA(1,2) = 118.7766 ; LNA(1,3) = 327.86658;
LNA(2,1) = 1.51565 ; LNA(2,2) = 109.7751 ; LNA(2,3) = 175.0000;
LNA(3,1) = 1.42159 ; LNA(3,2) = 111.12936 ; LNA(3,3) = 320.000;
LNA(4,1) = 1.54051 ; LNA(4,2) = 117.73836 ; LNA(4,3) = 118.01846;
LNA(5,1) = 1.39391 ; LNA(5,2) = 112.03748 ; LNA(5,3) = 62.58704;
LNA(6,1) = 1.53941 ; LNA(6,2) = 90.03892 ; LNA(6,3) = 180.66502;
LNA(7,1) = 1.42739 ; LNA(7,2) = 106.0307 ; LNA(7,3) = 164.74641;
LNA(8,1) = 1.54026 ; LNA(8,2) = 116.66205 ; LNA(8,3) = 236.54801;
LNA(9,1) = 1.41639 ; LNA(9,2) = 103.07467 ; LNA(9,3) = 302.62074;
LNA(10,1) = 1.0985 ; LNA(10,2) = 108.22029 ; LNA(10,3) = 170.49375;
LNA(11,1) = 1.09762 ; LNA(11,2) = 108.17067 ; LNA(11,3) = 51.90955;
LNA(12,1) = 1.09405 ; LNA(12,2) = 109.81257 ; LNA(12,3) = 118.20627;
LNA(13,1) = 1.08713 ; LNA(13,2) = 119.22717 ; LNA(13,3) = 181.55563;
LNA(14,1) = 1.09552 ; LNA(14,2) = 111.23163 ; LNA(14,3) = 295.66211;
LNA(15,1) = 1.09353 ; LNA(15,2) = 111.17992 ; LNA(15,3) = 78.50705;
LNA(16,1) = 1.0932 ; LNA(16,2) = 112.24301 ; LNA(16,3) = 314.42593;


Linkage(1,1) = 1.56703 ; Linkage(1,2) = 103.04905 ; Linkage(1,3) = 164.19925 ;
Linkage(2,1) = 1.47853 ; Linkage(2,2) = 100.5589 ; Linkage(2,3) = 146.25452 ;
Linkage(3,1) = 1.48053 ; Linkage(3,2) = 127.11291 ; Linkage(3,3) = 11.83068 ;
Linkage(4,1) = 1.5959 ; Linkage(4,2) = 90.62256 ; Linkage(4,3) = 256.60228 ;

Terminal(1) = 0.95935;  Terminal(2) = 109.38220;  Terminal(3) = 51.44904;


 return
 end


