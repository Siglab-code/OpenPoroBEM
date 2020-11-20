!****************************************************************************************
!
!							 CONTROL Subroutine
!
!	AIM			:	controls the program
!	callED IN	:	HYBRID
!	callS		:	DIM1,INPUT,NumToStr,DIM2,DETNA,INIT1,INIT3,
!					BEM,CHARGE1,CHARGE2,CHARGE3,CHARGE4,
!					SEISMIC,SOLVE,UPDATE1,FSTRESS,UPDATE3,OUTPUT,DEL
!
!	VARIABLES	:
!
!	in block /a/, these variables are presented:
!
!		nnp		: total number of nodes in problem (FE &BE zones)
!		ifem	: FEM code [0,1];
!							[0]: program studies just the topographical effect by BEM
!							     (FE is not done),
!							[1]: program studies geotechnical effect in a 2D topographical
!							     structure by hybrid FEM/BEM (FE is done)
!		ibem	: BEM code [0,1,2,3];
!							[0]: program doesn't consider BEM analysis. except some
!								 simple examples IBEM must never be [0]
!							[1]: BEM_dry
!							[2]: BEM_saturate
!							[3]: BEM_unsaturate
!		idyn	: dynamic code or type of loading [0,1,2];
!							[0]: static
!							[1]: quasi-static
!							[2]: dynamic problem
!		ieaq	: earthquake code [0,1,2,3];
!							[0]: static case
!							[1]: accelaration
!							[2]: imposed deplacement
!							[3]: incident wave, in most of the time IEAQ=[3]
!		nr		: type newton raphson code [0,1] in FE calculation;
!							[0]: modified
!							[1]: normal
!
!	in block /a1/, these variables are presented:
!
!		ne8d	: number of drained (dry) elements in FE zone
!		ne8c	: number of consolidated (saturated) elements in FE zone
!		ne8u	: number of unsaturated elements in FE zone
!
!	in block /a2/, these variables are presented:
!
!		nbed1	: number of D1 BEM areas (finite areas), it can be at maximum 5 zones
!		nbed2	: number of D2 BEM areas (infinite areas), it can be at maximum 5 zones
!		nbed3	: number of D3 BEM areas (semi-infinite area), it can be at maximum 1 zone
!		nnpbed1	: number of nodes of each D1 BEM area; NNPBED1(1 to NBED1)
!		nnpbed2	: number of nodes of each D2 BEM area; NNPBED2(1 to NBED2)
!		nnpbed3	: number of nodes of one D3 BEM area
!		nbeint	: number of internal nodes in BE D3 zone
!		nnpen	: number of nodes of fictitious enclosing elements in D3 BEM area, it can be
!				  at maximum 500 nodes
!		npbmax	: maximum nodes in all of the BE zone + 1; for example if we consider D3 BE
!				  area NPBMAX = nnpbed3 + 1
!		n12gh	: in geomechanical problems there is more than one unknown per node, so the
!				  matrices in BEM are assambled according to the doF, rather than node
!				  number. for dry soil each node has [2] doF (ux,uy), for for saturated soil
!				  each node has [3] (ux,uy,pw) doF & for unsaturated soil each node has
!				  [4] doF(ux,uy,pw,pa). then total number of doF in BE=(npbmax-1)*kbem
!		n3gh	: in semi-infinite problem in BEM we have to consider a Fictitious enclosing
!				  surface to perform point collocation method. N3GH shows the total number
!				  of doF in (D3+enclosing surface) in BE technique=(npbmax-1+nnpen)*kbem
!		kfsol	: kind of fundamental solution in saturated zone [1,2];
!							[1]: numerical laplace inverse
!							[2]: analytical laplace inverse (incompressible composent)
!
!	in block /a3/, these variables are presented:
!
!		ne8d1	: this parameter is defined for storage: ne8d1=ne8d if ne8d=0 then ne8d1=1
!		ne8c1	: this parameter is defined for storage: ne8c1=ne8c if ne8c=0 then ne8c1=1
!		ne8u1	: this parameter is defined for storage: ne8u1=ne8u if ne8u=0 then ne8u1=1
!		nload1	: this parameter is defined for storage: nload1=nload if nload=0 then
!				  nload1=1
!		nbcx1	: this parameter is defined for storage: nbcx1=nbcx if nbcx=0 then nbcx1=1
!		nbcy1	: this parameter is defined for storage: nbcy1=nbcy if nbcy=0 then nbcy1=1
!		nbcw1	: this parameter is defined for storage: nbcw1=nbcw if nbcw=0 then nbcw1=1
!		nbca1	: this parameter is defined for storage: nbca1=nbca if nbca=0 then nbca1=1
!		m8d1	: this parameter is defined for storage: m8d1=m8d if m8d=0 then m8d1=1
!		m8c1	: this parameter is defined for storage: m8c1=m8c if m8c=0 then m8c1=1
!		m8u1	: this parameter is defined for storage: m8u1=m8u if m8u=0 then m8u1=1
!		nbed11	: this parameter is defined for storage: nbed11=nbed1 if nbed1=0 then
!				  nbed11=1
!		nbed21	: this parameter is defined for storage: nbed21=nbed2 if nbed2=0 then
!				  nbed21=1
!		nbed31	: this parameter is defined for storage: nbed31=nbed3 if nbed3=0 then
!				  nbed31=3
!		ne8dd	: this parameter is defined for storage in plasticity case: ne8dd=ne8d
!				  if idplst=0 then ne8dd=1
!		ne8cc	: this parameter is defined for storage in plasticity case: ne8cc=ne8c
!				  if icplst=0 then ne8cc=1
!		ne8uu	: this parameter is defined for storage in plasticity case: ne8uu=ne8u
!				  if iuplst=0 then ne8uu=1
!
!	in block /a4/, these variables are presented:
!
!		xencl: X-coordinates of enclosing elements' nodes (nodal data)
!		yencl: Y-coordinates of enclosing elements' nodes (nodal data)
!
!	in block /b/, these variables are presented:
!
!		m8d		: number of drained (dry) materials in FE & BE zones; in BE area the program
!				  is written for just one material (homogeneous space) then if we use a dry
!				  material in BE zone (by considering the dry materials in FE zone) the
!				  number of material of BE zone is always equal to M8D because the materials
!			      from "1" to "M8D-1" are used in FE zone
!		m8c		: number of consolidated (saturated) materials in FE & BE zones; in BE area
!				  the program is written for just one material (homogeneous space) then if
!				  we use a saturated material in BE zone (by considering the saturated
!				  materials in FE zone) the number of material of BE zone is always equal to
!				  M8C because the materials from "1" to "M8C-1" are used in FE zone
!		m8u		: number of unsaturated materials in FE & BE zones; in BE area the program
!				  is written for just one material (homogeneous space) then if we use an
!				  unsaturated material in BE zone (by considering the unsaturated materials
!				  in FE zone) the number of material of BE zone is always equal to M8U
!				  because the materials from "1" to "M8U-1" are used in FE zone
!		m8dep	: behaviour type of drained (dry) materials [1,2,3,4];
!								[1]: linear,
!								[2]: hyperbolic for statics,
!								[3]: hyperbolic for dynamics,
!								[4]: elastoplastic
!				  if we use a dry material in BE zone then the behaviour type of BE zone's
!				  material is: M8DEP(m8d)
!		m8cep	: behaviour type of saturated materials [1,2,3,4];
!								[1]: linear,
!								[2]: hyperbolic for statics,
!								[3]: hyperbolic for dynamics,
!								[4]: elastoplastic
!				  if we use a saturated material in BE zone then the behaviour type of BE
!				  zone's material is: M8CEP(m8c)
!???	m8uep	: behaviour type of unsaturated materials [1,2,3,4];
!								[1]: linear,
!								[2]: hyperbolic for statics,
!								[3]: hyperbolic for dynamics,
!								[4]: elastoplastic
!				  if we use an unsaturated material in BE zone then the behaviour type of BE
!				  zone's material is: M8UEP(m8u)
!
!	in block /c/, these variables are presented:
!
!		nload	: number of load steps; load is applied in "nload" steps, in each
!				  "iload=1:nload" the characteristics of loading such as 'intensity-nstep-...
!				  are differents
!		ntime	: number of time steps or total number of "nstep"s in all "nload"s. loading
!				  is divided into SUM[nstep]s, when the characteristics of "nstep"s are
!				  similar we put these "nsteps" in one group which is named iload
!		itmax	:
!		itime	:
!		iter	:
!		iload	:
!		time	:
!		dtime	:
!		dtimei	: minimal timestep size
!		dtimem	: maximal timestep size
!		xinc	:
!
!	in block /d/, these variables are presented:
!
!		nbcx	: number of (null) boundary conditions introduced for x-displacements
!		nbcy	: number of (null) boundary conditions introduced for y-displacements
!		nbcw	: number of (null) boundary conditions introduced for nodal water pressure
!		nbca	: number of (null) boundary conditions introduced for nodal air pressure
!		nbc		: total number of null boundary conditions; nbc = nbcx+nbcy+nbcw+nbca
!
!	in block /e/, these variables are presented:
!
!		gammaw	: water volumetric weight
!		gammaa	: air volumetric weight
!		grav	: gravity acceleration
!		atmp	: ambient atmospheric pressure (taken constant)
!
!	in block /g/, these variables are presented:
!
!		init	: code describing the type of initial conditions [0,1,2,3]. when we want to
!				  do a topographical effect (ifem=0) we continue by INIT=[0]
!							[0] if the results are taken from a previous analysis
!							[1] if initial normal stresses are considered equal to the
!								weight of the upper soil layers
!							[2] if initial stresses are computed for the steady state of
!							    the system
!							[3] if initial stresses are calculated as the weight of the
!								upper soil layers, considering that the system has a mere
!							    parallelogram shape, and that the x-direction of the
!							    space-frame corresponds to the physical horizontal
!		iprint	: initial data output print code (control code for the input and output files)
!				  [0,1];
!							[0] to avoid the outprint of the mesh file
!							[1] to print out the mesh file
!
!	in block /h/, these variables are presented:
!
!		mdof	: maximum degree of freedom by considering the boundary conditions
!		mdofn	: maximum degree of freedom without considering the boundary conditions
!		ls		:
!		ls1		:
!		isolv	: solution type indicator [1,3];
!  								[1]: for a symmetric computation,
!								[3]: for a non-symmetric computation problem; in most of the
!									 case is solved by ISOLV=[3]
!
!	in block /h1/, these variables are presented:
!
!???	ptoll1	: time-integration constant in FE formulation, for isolve=3,TETA>0.5
!???	ptoll3	: permissible tollerance in FE formulation
!
!	in block /i/, these variables are presented:
!
!		iterp	: number of iloads for which the results are printed out in the output file
!		itprt	: iload's numbers which are selected to print out the results;
!							[0] if no printout is expected
!		ne8do: number of drained elements that we need their outputs, it must be at maximum
!			   10 elms.
!		ne8co: number of saturated elements that we need their outputs
!		ne8uo: number of unsaturated elements that we need their outputs
!		ie8dout	: numbers of the nodes forming the current drained element that we need
!				  their outputs, IE8doUT (1:NE8do)
!		ie8cout	: numbers of the nodes forming the current saturated element that we need
!				  their outputs, IE8COUT (1:NE8CO)
!		ie8uout	: numbers of the nodes forming the current unsaturated element that we need
!				  their outputs, IE8UOUT (1:NE8UO)
!		nnpo: number of nodes that we need their outputs, it must be at maximum 80 nodes.
!		inpout	: nodes' numbers which we need their outputs, INPOUT (1:nnpo)
!
!	in block /chrg1/, these variables are presented:
!
!		nstep	:
!		istep	:
!
!	in block /chrg2/, these variables are presented:
!
!***	npw		: number of nodes on which water pressure is imposed for the current loadstep (NPW <= 1000)
!***	npa		: number of nodes on which air pressure is imposed for the current loadstep (NPA <= 1000)
!***	nnpw	: node numbers of the nodes on which water pressure is imposed for the current loadstep
!***	nnpa	: node numbers of the nodes on which air pressure is imposed for the current loadstep
!***	vnpw	: water pressure values imposed on the preceding nodes for the current loadstep
!***	vnpa	: air pressure values imposed on the preceding nodes for the current loadstep
!
!	in block /chrg3/, these variables are presented:
!
!		nstime:
!		neaq:
!		nneaq(500):
!		uxeaq:
!		uyeaq:
!		axeaq:
!		ayeaq:
!		uxeq(3,500):
!		uyeq(3,500):
!
!	in block /ini1/, these variables are presented:
!
!		vnpwi	: initial water pressure values imposed on the preceding nodes
!		vnpai	: initial air pressure values imposed on the preceding nodes
!		nnpwi	: nodes' numbers on which an initial water-p is imposed
!		nnpai	: nodes' numbers on which an initial air-p is imposed
!
!	in block /ini2/, these variables are presented:
!
!		npwi	: number of nodes on which an initial water pressure is imposed which can
!				  be at max. 400 nodes
!		npai	: number of nodes on which an initial gas pressure is imposed which can be
!				  at max. 400 nodes
!		ndispi	: number of nodes on which an initial displacement is imposed which can be
!				  at max. 400 nodes
!		nveli	: number of nodes on with an initial velocity is imposed which can be at
!				  max. 400 nodes
!		nacci	: number of nodes on with an initial acceleration is imposed which can be
!				  at max. 400 nodes
!
!	in block /ini3/, these variables are presented:
!
!		nndispi	: nodes' numbers on which an initial disp. is imposed
!		nnveli	: nodes' numbers on which an initial vel. is imposed
!		nnacci	: nodes' numbers on which an initial acc. is imposed
!
!	in block /ini4/, these variables are presented:
!
!		vxdispi	: initial displacement values imposed on the preceding nodes (X-dir)
!		vydispi	: initial displacement values imposed on the preceding nodes (Y-dir)
!
!	in block /ini5/, these variables are presented:
!
!		vxveli	: initial velocity values imposed on the preceding nodes (X-dir)
!		vyveli	: initial velocity values imposed on the preceding nodes (Y-dir)
!
!	in block /ini6/, these variables are presented:
!
!		vxacci	: initial acc. values imposed on the preceding nodes (X-dir)
!		vyacci	: initial acc. values imposed on the preceding nodes (Y-dir)
!
!
!		L73	 -> S1 ; L74 -> S2 ;
!		Lend -> storage required for arrays
!		Lsx -> storage for stiff
!		Lmax -> storage allowed
!
!
!****************************************************************************************
!	
	subroutine control(a,lmax)
!
	implicit double precision (a-h,o-z)
	character (72) :: infile, outfile, out*4
!
	common /a/ nnp,ifem,ibem,idyn,ieaq,nr
	common /a1/ ne8d,ne8c,ne8u
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /a3/ ne8d1,ne8c1,ne8u1,ne3u1,nload1,nbcx1,nbcy1,nbcw1,nbca1,m8d1,m8c1,m8u1,&
				nbed11,nbed21,nbed31,nbeint1,ne8dd,ne8cc,ne8uu
	common /a4/ xencl(500),yencl(500)
	common /b/ m8d,m8c,m8u,m8dep(10),m8cep(10),m8uep(10)
	common /c/ nload,ntime,itmax,itime,iter,iload,time,dtime,dtimei,dtimem,xinc
	common /d/ nbcx,nbcy,nbcw,nbca,nbc
	common /e/ gammaw,gammas,gammaa,grav,atmp
	common /g/ init,iprint
	common /h/ mdof,mdofn,ls,ls1,isolv
	common /h1/ ptoll1,ptoll3
	common /i/iterp,itprt(100),ne8do,ne8co,ne8uo,ie8dout(10),ie8cout(10),ie8uout(10),nnpo,&
			  inpout(100)
	common /chrg1/ nstep,istep
	common /chrg2/ npw,npa,nnpw(200),nnpa(200),vnpw(200),vnpa(200)
	common /chrg3/ nstime,neaq,nneaq(500),uxeaq,uyeaq,axeaq,ayeaq,uxeq(3,500),uyeq(3,500)
	common /ini1/ vnpwi(400),vnpai(400),nnpwi(1000),nnpai(400)
	common /ini2/ npwi,npai,ndispi,nveli,nacci
	common /ini3/ nndispi(400),nnveli(400),nnacci(400)
	common /ini4/ vxdispi(400),vydispi(400)
	common /ini5/ vxveli(400),vyveli(400)
	common /ini6/ vxacci(400),vyacci(400)
!
	dimension :: a(lmax)
!
!
	write (*,*)
	write (*,*) ' Dynamic Analysis of Dry/Saturated/Unsaturated Porous Media &
										With FEM / BEM		'
	write (*,*)
	write (*,*) '				OpenPoroBEM (Version 2020)				'
	write (*,*)
!
	write (*,*) 'name of input file: '
	read (*,*) infile
!
	write (*,*) 'name of output file: '
	read (*,*) outfile
!
	open (UNIT=4, FILE=infile, STATUS='OLD')
	open (UNIT=7, FILE=outfile, STATUS='REPLACE')
!
!
!	----- call subroutine to read basic data and allocate storage (subroutine dim1) -----
!
	write (*,*)
	write (*,*) 'dim1 started'
	call dim1(l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14,l15,l16,l17,l18,l19,l20,l21,&
			  l22,l23,l24,l25,l26,l27,l28,l29,l30,l31,l32,l33,l34,l35,l36,l37,l38,l39,l40)
	write (*,*) '			dim1 concluded'
!
!
!	--------- call subroutine to read mesh and material data (subroutine input) ---------
!
	write (*,*) 'input started'
	call input(a(l1),a(l2),a(l3),a(l8),a(l9),a(l10),a(l11),a(l12),a(l13),a(l14),a(l15),&
			   a(l16),a(l17),a(l4),a(l5),a(l6),infile)
	write (*,*) '			input concluded'
!
!
!	---------------- opening output files for selected nodes & elements -----------------
!
	j=9
!   nodes
	if (nnpo.eq.0) goto 8
	do i=1,nnpo
		j=j+i
		call NumToStr(out,inpout(i))
		open (j,file='N'//out//'.txt')
	enddo
	j=j+1
	open (j,file='NodesOut.txt')
!
8	if (nbeint.eq.0) goto 10
	do i=1,nbeint
		j=j+i
		call NumToStr(out,i)
		open (j,file='Nint'//out//'.txt')
	enddo
	j=j+1
	open (j,file='NodesINTOut.txt')
!
10	continue
!
!	dry elements
	if (ne8do.eq.0) goto 20
	do i=1,ne8do
		j=j+i
		call NumToStr(out,ie8dout(i))
		open (j,file='DE'//out//'.txt')
		write (j,4)
	enddo
20	continue
!
!	consolidation elements
	if (ne8co.eq.0) goto 30
	do i=1,ne8co
		j=j+i
		call NumToStr(out,ie8cout(i))
		open (j,file='CE'//out//'.txt')
		write (j,5)
	enddo
30	continue
!
!
!	---- call subroutine for further allocation of dynamic storage (subroutine dim2) ----
!
	write (*,*) 'dim2 started'
	call dim2(a(l1),a(l2),a(l3),a(l7),l40,l41,l42,l43,l44,l45,l46,l47,l48,l49,l50,l51,l52,&
			  l53,l54,l55,l56,l57,l58,l59,l60,l61,l62,l63,l64,l65,l66,l67,l68,l69,l70,l71,&
			  l72,l73,a(l4),a(l5),a(l6))
	write (*,*) '			dim2 concluded'
!
!
!	---------- call subroutine to determine na vector (subroutine detna) ----------
!
	write (*,*) 'detna started'
	call detna(a(l1),a(l2),a(l3),a(l7),a(l9),a(l10),a(l11),a(l12),a(l44),a(l50),a(l4),a(l5),&
			   a(l6))
	write (*,*) '			detna concluded'
!
!	---------------- check storage ----------------
!
	l74=l73+ls	! S1
	l75=l74+ls1	! S2
	if (ieaq.ne.2) then
		l76=l75+1
		l77=l76+1
	else
        l76=l75+ls
		l77=l76+ls1
	endif
	l78=l77+mdof
	lend=l78+mdof
	lsx=ls+ls1
	write (7,3) lend,lsx
	if (lend.le.lmax) goto 50
	write (7,1) lend,lmax
	stop
50	continue
!
!
!	----------------------------------- initializing -----------------------------------
!
	if (init.gt.3) goto 1100
	write (*,*)
	write (*,*)'********************** initializing **********************'
	write (*,*)
!
!
!	equaling all displacement, velocity, acceleration and stress arrays to zero
!
	write (*,*)'init1 started'
	call init1(a(l1),a(l2),a(l3),a(l18),a(l19),a(l20),a(l24),a(l25),a(l26),a(l27),&
			   a(l28),a(l29),a(l21),a(l22),a(l23),a(l39),a(l40),a(l41),a(l46),a(l51),&
			   a(l4),a(l5),a(l6),a(l44),a(l45),a(l47),a(l48),a(l49),a(l42),a(l43),&
			   a(l77),a(l78))
	write (*,*) '			init1 concluded'
!
!
!   inputs initial values of displacements, velocities, accelerations and pore pressurs
!
	write (*,*)'init3 started'
	call init3(a(l7),a(l47),a(l48),a(l49),a(l39),a(l17),a(l4),a(l5),a(l6),a(l13),a(l45),&
			   a(l46),a(l44))
	write (*,*) '			init3 concluded'
	write (*,*)
	write (*,*)'**********************************************************'
!
!
	if (init.eq.0) goto 300	! it means that if the results are taken from a previous analysis
	if (ifem.eq.0) goto 300
!
!
!*********************************************************************************************
!	------------------- evaluating overburden stresses of the FEM area ------------------
!	dtime=1.0	! to calculate the initial stresses the min timestep is equalized to 1
!	write (*,*)'****************** evaluation of fem overburden stresses ******************'
!	goto (100,100,200) init
!100	continue
!	if (ibem.eq.0) goto 150
!	write (*,*)'bem started'
!	iter=1
!	itime=1
!	idyn1=idyn
!	idyn=0
!	call bem(a(l7),a(l48),a(l13),a(l14),a(l15),a(l16),a(l71),a(l72),a(l4),a(l5),a(l6),a(l41),a(l45),&
!			 a(l50),a(l51),a(l52),a(l53),a(l54),a(l55),a(l56),a(l57),a(l58),a(l59),a(l60),a(l61),&
!			 a(l62),a(l63),a(l64),a(l65),a(l66),a(l67),a(l68),a(l69),a(l70))
!	idyn=idyn1
!150	if (ifem.eq.0) goto 300
!	write (*,*)'fem started'
!	call fem(a(l1),a(l2),a(l6),a(l11),a(l12),a(l13),a(l14),a(l15),a(l18),a(l19),a(l20),a(l21),a(l16),&
!			 a(l17),a(l28),a(l31),a(l32),a(l34),a(l38),a(l61),a(l62),a(l35),a(l36),a(l37),a(l33),a(l22),&
!			 a(l23),a(l24),a(l25),a(l26),a(l27),a(l65),a(l66))
!	write (*,*)'solve started'
!	call solve(a(l61),a(l62),a(l28),a(l38))
!   write (*,*)'update1 started'
!	call update1(a(l28),a(l32),a(l34),a(l33))
!   write (*,*)'stress started'
!	call fstress(a(l1),a(l2),a(l6),a(l11),a(l12),a(l13),a(l14),a(l15),a(l18),a(l19),a(l20),a(l21),a(l16),&
!				 a(l17),a(l32),a(l35),a(l34),a(l22),a(l23),a(l24),a(l25),a(l26),a(l27))
!200	write (*,*)'init2 started'
!	call init2(a(l1),a(l2),a(l11),a(l12),a(l13),a(l14),a(l15),a(l18),a(l19),a(l20),a(l21),a(l16),a(l17))
!	           ,a(l28),a(l29),a(l34),a(l6))
!*********************************************************************************************
!
!
300	continue
!
!
!	-------------------------------- start time increment -------------------------------
!	-------------------------------------------------------------------------------------
!	-------------------------------------------------------------------------------------
!
	time=0.0d0
	dtime=dtimei
	iload=0
!
	write (*,*)
	write (*,*)'************************ loading *************************'
	write (*,*)
	do 1000 itime=1,ntime
!!!!		write (*,6) itime,ntime
!
!
!	check for loading; zero all internal nodal loads
!
		k=0
		if (nload.eq.0) goto 410
		write (*,*) 'charge1 started'
		call charge1 (a(l8),a(l43),nload,itime,k,i)
		write (*,*) '			charge1 concluded'
		if (k.eq.0) goto 410
		iload=iload+1
		write (7,2) i
410		continue
!
!
!	define new increments of external nodal loads, flows (and seismic effects)
!
		if (k.eq.0) goto 420
		write (*,*)'charge2 started'
		call charge2 (a(l7),a(l8),a(l13),a(l41),a(l47))
		write (*,*) '			charge2 concluded'
420		continue
!
!
!	update time and istep
!
		time=time+dtime
!
		if (k.eq.1) istep=nstep
		istep=istep-1
!
!
!	apply the sum of all external nodal loads
!
		if (istep.lt.0) goto 430
		write (*,*)'charge3 started'
		call charge3(a(l41),a(l42))
		write (*,*) '			charge3 concluded'
430		continue
!
!
!	-------------------------------- start iteration -------------------------------
!
		do 500 iter=1,itmax
!
			write (*,*)
			write (*,*)'**********************************************************'
			write (*,*)
			write (*,*)'itime= ',itime,'  iload= ',iload,'  iter= ',iter
			write (7,*)'itime= ',itime,'  iload= ',iload,'  iter= ',iter
!
!
!	assemble bem stiffness matrices & load vectors
!
			if (ibem.eq.0) goto 440
			if (iter.gt.1.and.nr.eq.1) goto 440
			write (*,*)'bem started'
			call bem(a(l7),a(l50),a(l13),a(l14),a(l15),a(l16),a(l17),a(l73),a(l74),a(l4),a(l5),&
					a(l6),a(l39),a(l40),a(l43),a(l47),a(l52),a(l53),a(l54),a(l55),a(l56),a(l57),&
					a(l58),a(l59),a(l60),a(l61),a(l62),a(l63),a(l64),a(l65),a(l66),a(l67),a(l68),&
					a(l69),a(l70),a(l71),a(l72))
			write (*,*) '		bem concluded'
440			continue
!
!
!	assemble fem stiffness matrices & load vectors
!
!      write (*,*)'fem started'
!      call fem(a(l1),a(l2),a(l6),a(l11),a(l12),a(l13),a(l14),
!     $a(l15),a(l18),a(l19),a(l20),a(l21),a(l16),a(l17),a(l28),
!     $a(l31),a(l32),a(l34),a(l38),a(l61),a(l62),a(l35),a(l36),a(l37),
!     $a(l33),a(l22),a(l23),a(l24),a(l25),a(l26),a(l27),a(l65),a(l66))
!
!
460			continue
!
!
!	apply the sum of all (ext. & int.) nodal loads
!
			write (*,*)'charge4 started'
			call charge4(a(l40),a(l42),a(l43),a(l77),a(l78))
			write (*,*) '		 charge4 concluded'
!
!
!	take account sismic event
!
			if (ieaq.eq.2) then
				write (*,*) 'seismic started'
				call seismic(a(l7),a(l50),a(l73),a(l74),a(l75),a(l76),a(l40),iter)
				write (*,*) '		 seismic concluded'
			endif
!
!
!	solve for displacements and pore pressures
!
			write (*,*)'solve started'
			call solve(a(l73),a(l74),a(l40),a(l50))
			write (*,*) '		 solve concluded'
!
!
!	update displacement arrays after each iteration
!
			write (*,*)'update1 started'
			call update1(a(l40),a(l44),a(l46),a(l45))
			write (*,*) '		 update1 concluded'
!
!
!	calcul stress arrays after each iteration
!
			write (*,*)'stress started'
			call fstress(a(l1),a(l2),a(l3),a(l4),a(l7),a(l13),a(l15),a(l16),a(l17),a(l39),&
						 a(l70),a(l55),a(l61),a(l18),a(l19),&
						 a(l20),a(l24),a(l25),a(l26),a(l27),a(l28),a(l29),a(l21),a(l22),&
						 a(l23),a(l44),a(l47),a(l46),a(l30),a(l31),a(l32),a(l33),a(l34),&
						 a(l35),a(l36),a(l37),a(l38))
			write (*,*) '		 fstress concluded'
!
!	call subroutine to control the tollerance
!
			if (nr.eq.0) goto 550
			write (*,*)'toll started'
			call toll(a(l44),a(l46),a(l1),a(l2),a(l3),a(l7),a(l77),a(l78),tol1,tol2,tol3,&
					  tol4,tol5)
			write (7,*) tol1,tol4
!			if (tol1.le.ptoll1 .and. tol4.le.ptoll3) goto 550
			if (tol1.le.ptoll1) goto 550
!
!
500		continue	! the end of the iteration loop
!
!	-------------------------------------------------------------
!
550		continue
!
!
!	update displacement, velocity, acceleration and stress arrays after each increment
!
		write (*,*)'update3 started'
		call update3(a(l47),a(l48),a(l49),a(l18),a(l19),a(l20),a(l27),a(l28),a(l29),a(l46))
		write (*,*)'		update3 concluded'
!
!
!	solution at any interior point in BE zone
!
		if (nbeint.eq.0) goto 551
		call intdisp(a(l7),a(l13),a(l14),a(l15),a(l16),a(l17),a(l4),a(l5),a(l6),a(l55),a(l56),&
					 a(l57),a(l61),a(l62),a(l63),a(l70),a(l71),a(l72),a(l47),a(l39),&
					 a(l64),a(l65),a(l66),a(l67),a(l68),a(l69))
!
551		continue
!
!	output the results
!
		write (*,*)'output started'
		call output(a(l1),a(l2),a(l3),a(l7),a(l13),a(l18),a(l19),a(l20),a(l46),a(l51),a(l39))
		write (*,*)'output concluded'

!
!	update time step
!
		dtime=dtime*xinc
		if (dtime.gt.dtimem) dtime=dtimem
		if (dtime.lt.dtimei) dtime=dtimei
!
!
!   -------------------------------------------------------------
1000 continue
1100 continue
!
!
!	----------------- delete the temporary files ----------------
!
	call del(ntime,nbed1,nbed2,nbed3,nbeint)
!
!
!	------------------------- formats -------------------------
!
1	FORMAT (//10x,'storage exceeded'//10x,'storage required = ',i7//&
			  10x,'storage allowed  = ',i7)
2	FORMAT (//10x,'loading step number = ',i7//)
3	FORMAT (//10x,'storage for arrays = ',i7/&
              10x,'storage for stiff  = ',i7)
4	FORMAT (//10x,'Stress and strain at the center of the element',//,&
				  'time   sigma_x   sigma_y   sigma_xy   sigma_m',&
				  '   tau   epxilon_x epxilon_y   gamma_xy   gamma',//)
5	FORMAT (//10x,'Stress and strain at the center of the element',//,&
				  'time   sigma_x   sigma_y   sigma_xy  sigma_m   tau',&
				  '   pression  epxilon_x   epxilon_y   gamma-xy   gamma',//)
!!!!6	FORMAT (/20x, 'Time step is:'i6', From:'i8)
!

	return
	end
!****************************************************************************************
!
!							           DMATL subroutine
!
!	MODEL ELASTICITY LINEAR.
!
!	this subroutine is called in: STR8ND
!
!	variables used are:
!
!	INPUT	:
!
!	OUTPUT	:
!
!
!****************************************************************************************
!
	subroutine dmatl(kode,mt,d,sm8d,sm8c)
!
 	implicit double precision (a-h,o-z)
!
	common /a3/ ne8d1,ne8c1,ne8u1,ne3u1,nload1,nbcx1,nbcy1,nbcw1,nbca1,m8d1,m8c1,m8u1,nbed11,nbed21,nbed31,nbeint1,ne8dd,ne8cc,ne8uu
!
	dimension d(3,3),sm8d(60,m8d1),sm8c(60,m8c1)
!
!
	if (kode.eq.1) then
!	module Young
		et=sm8d(3,mt)
!	module de compressibilit�
		bt=sm8d(7,mt)
	endif
!
	if (kode.eq.2) then
		et=sm8c(3,mt)
		bt=sm8c(7,mt)
	endif
!
	d1=3.0*bt/(9.0*bt-et)
	d2=3.0*bt+et
      d3=3.0*bt-et
      d(1,1)=d1*d2
      d(2,2)=d(1,1)
      d(1,2)=d1*d3
      d(2,1)=d(1,2)
      d(1,3)=0.0
      d(2,3)=0.0
      d(3,1)=0.0
      d(3,2)=0.0
      d(3,3)=d1*et
!
	return
	end
!c
!!!!!!!c*************************************!!!!!!c*********************************
!c
!!!!cMODEL ELASTICITY HYPERBOLIC (DUNCAN & CHEN)
!c
!!!!!!!c*************************************!!!!!!c*********************************
!c
      subroutine dmath1(kode,mt,d,sm8d,sm8c,sstr)
!c
      implicit double precision (a-h,o-z)
      common /a3/ne8d1,ne8c1,nload1,nbcx1,nbcy1,nbcw1,m8d1,m8c1,nbed11,nbed21,nbed31,ne8dd,ne8cc
      common /e/gammaw,grav,atmp
      common /g/init,iprint
      common /c/nload,ntime,itmax,itime,iter,iload,time,dtime,dtimei,dtimem,beta,xinc
!c
      dimension d(3,3),sstr(6),sm8d(60,m8d1),sm8c(60,m8c1)
      data einit/48.d6/,  binit/40.d6/
!c
      if (init.eq.0) goto 50
        et=einit ; bt=binit
        goto 1000
50    continue
      sx=-sstr(1)    ; sy=-sstr(2)
      sxy=-sstr(3)
      if (kode.eq.1) then
        phi=sm8d(1,mt) ; c=sm8d(2,mt)
        xkl=sm8d(3,mt) ; xku=sm8d(4,mt) ; xkb=sm8d(7,mt)
        xn=sm8d(5,mt)  ; xm=sm8d(8,mt)
        rf=sm8d(6,mt)  ; ef=sm8d(9,mt)
        ten=sm8d(12,mt)
      endif
      if (kode.eq.2) then
        phi=sm8c(1,mt) ; c=sm8c(2,mt)
        xkl=sm8c(3,mt) ; xku=sm8c(4,mt) ; xkb=sm8c(7,mt)
        xn=sm8c(5,mt)  ; xm=sm8c(8,mt)
        rf=sm8c(6,mt)  ; ef=sm8c(9,mt)
        ten=sm8c(12,mt)
      endif
!c
!!!c check for tensile failure
!c
100   call princp(sx,sy,sxy,s1,s3)
      smean=s3
      if (s3.le.0.d0) stop
!!!cap1=0.1*atmp
!!!ctens=s3+ten
!!!cif (s3.gt.0.d0) goto 300
!!!cif (tens.gt.0.d0) goto 200
!!!csmean=ap1
!!!cbt=xkb*atmp*(smean/atmp)**xm
!!!cet=0.01*bt
!!!csr=1.0000
!!!cgoto 1000
!c200   smean=tens
!c
300   sr=(1.d0-dsin(phi))*(s1-smean)/2.d0/(c*dcos(phi)+smean*dsin(phi))
!!!cif (sr.ge.0.99) stop
!!!caa=0.75d0
!!!cif (sr.le.sstr(4)*aa) then
!!!!cet=xku*atmp*(smean/atmp)**xn
!!!celseif (sr.ge.sstr(4)) then
        et=xkl*atmp*((smean/atmp)**xn)*(1.d0-rf*sr)**2
!!!celse
!!!!ce1=xku*atmp*(smean/atmp)**xn
!!!!ce2=xkl*atmp*((smean/atmp)**xn)*(1.d0-rf*sr)**2
!!!!cet=e2+(e1-e2)*(1.d0-sr/sstr(4))/(1.d0-aa)
!!!cendif
!!!cbt=xkb*atmp*(smean/atmp)**xm
      xnuy=0.3 ; bt=et/3.d0/(1-2.d0*xnuy)
      btmax=17.d0*et ; btmin=et/3.d0
      if (bt.lt.btmin) bt=btmin
      if (bt.gt.btmax) bt=btmax
1000  continue
!c
!!c form d matrix
!c
      d1=3.0*bt/(9.0*bt-et)
      d2=3.0*bt+et
      d3=3.0*bt-et
      d(1,1)=d1*d2
      d(2,2)=d(1,1)
      d(1,2)=d1*d3
      d(2,1)=d(1,2)
      d(1,3)=0.0
      d(2,3)=0.0
      d(3,1)=0.0
      d(3,2)=0.0
      d(3,3)=d1*et
!c
      return
      end
!c
!!!!!!!c*************************************!!!!!!c*********************************
!c
!!!cMODEL ELASTICITY HYPERBOLIC (HARDIN & DRNEVICH)
!c
!!!!!!!c*************************************!!!!!!c*********************************
!c
      subroutine dmath2(kode,max,mt,d,sm8d,sm8c,sig,epi)
!!c kode=1:sec kode=2:saturated
!!c max=1:Gmax max=2:G tangeant max=3:G secant
!c
      implicit double precision (a-h,o-z)
      common /a3/ne8d1,ne8c1,nload1,nbcx1,nbcy1,nbcw1,m8d1,m8c1,nbed11,nbed21,nbed31,ne8dd,ne8cc
      common /e/gammaw,grav,atmp
      common /g/init,iprint
!c
      dimension d(3,3),sig(6),sm8d(60,m8d1),sm8c(60,m8c1),epi(3)
      data ginit/18.d6/,  binit/48.d6/
!!!cdata gomax/1.d6/, tomax/1.d6/
!c
      if (init.eq.0) goto 50
        gt=ginit ; bt=binit
        goto 100
50    sx=-sig(1) ; sy=-sig(2)
      sm=(sx+sy)/2.d0
      if (sm.lt.45.d3) sm=45.d3
!!!cgamma=epi(3)
      gamma=dsqrt((epi(1)-epi(2))**2+epi(3)**2)
      if (kode.eq.1) then
        ocr=sm8d(3,mt) ; zk=sm8d(4,mt)
        zn=sm8d(5,mt)  ; e=sm8d(6,mt)
        zkb=sm8d(7,mt) ; zm=sm8d(8,mt)
        phi=sm8d(1,mt) ; c=sm8d(2,mt)
      endif
      if (kode.eq.2) then
        ocr=sm8c(3,mt) ; zk=sm8c(4,mt)
        zn=sm8c(5,mt)  ; e=sm8c(6,mt)
        zkb=sm8c(7,mt) ; zm=sm8c(8,mt)
        phi=sm8c(1,mt) ; c=sm8c(2,mt)
      endif
      gmax=625.d0*(ocr**zk)*atmp*((sm/atmp)**zn)/(0.3d0+0.7d0*e**2)
!!!ctmax=dsqrt((sm*dsin(phi)+c*dcos(phi))**2-((sy-sx)/2.d0)**2)
      tmax=sm*dsin(phi)+c*dcos(phi)
!!!cif (gmax.lt.gomax) gmax=gomax
!!!cif (tmax.lt.tomax) tmax=tomax
      if (max.eq.1) gt=gmax
      if (max.eq.2) gt=gmax/((1.d0+abs(gamma)*gmax/tmax)**2)
      if (max.eq.3) gt=gmax/(1.d0+abs(gamma)*gmax/tmax)
      xnuy=1.d0/3.d0 ; bt=gmax*2.d0*(1.d0+xnuy)/3.d0/(1-2.d0*xnuy)
!!!cxnuy=1.d0/3.d0 ; bt=gt*2.d0*(1.d0+xnuy)/3.d0/(1-2.d0*xnuy)
!!!cbt=zkb*atmp*(sm/atmp)**zm
      btmax=gmax*49.667d0 ; btmin=gmax*0.667d0
!!!cbtmax=gt*49.667d0 ; btmin=gt*0.667d0
      if (bt.lt.btmin) bt=btmin
      if (bt.gt.btmax) bt=btmax
!c
!!c form d matrix
100   continue
      d(1,1)=bt+gt*4.d0/3.d0
      d(2,2)=d(1,1)
      d(1,2)=bt-gt*2.d0/3.d0
      d(2,1)=d(1,2)
      d(1,3)=0.0
      d(2,3)=0.0
      d(3,1)=0.0
      d(3,2)=0.0
      d(3,3)=gt
!c
      return
      end
!c

!!!!!!c*************************************!!!!!!c*********************************
      subroutine princp(sx,sy,sxy,s1,s3)
!!!!!!c*************************************!!!!!!c*********************************

!!c calculates principal stresses

      implicit double precision (a-h,o-z)

      cntr=(sx+sy)/2.d0
      hleg=(sx-sy)/2.d0
      rd=dsqrt(hleg*hleg+sxy*sxy)
      s1=cntr+rd
      s3=cntr-rd

      return
      end
!****************************************************************************************
!
!							           DEL subroutine
!
!	This subroutine deletes the temporary 'gh' files.
!
!	this subroutine is called in: CONTROL
!
!
!	variables used are:
!
!
!	INPUT	:
!
!
!	OUTPUT	:
!
!
!****************************************************************************************
!
	subroutine del(ntime,nbed1,nbed2,nbed3,nbeint)
!
	integer ntime,i
	character b*5
	logical fil
!
	if (nbed1.eq.0) goto 10
	do i=1,ntime
		call NumToStr(b,i)
		INQUIRE (file='ghd1'//b,exist=fil)
		if (fil) then
			open (2,file='ghd1'//b)
			close (2,status='delete')
		endif
	enddo
!
10	if (nbed2.eq.0) goto 20
	do i=1,ntime
		call NumToStr(b,i)
		INQUIRE (file='ghd2'//b,exist=fil)
		if (fil) then
			open (2,file='ghd2'//b)
			close (2,status='delete')
		endif
	enddo
!
20	if (nbed3.eq.0) goto 25
	do i=1,ntime
		call NumToStr(b,i)
		INQUIRE (file='ghd3'//b,exist=fil)
		if (fil) then
			open (2,file='ghd3'//b)
			close (2,status='delete')
		endif
	enddo
!
25	if (nbeint.eq.0) goto 30
	do i=1,ntime
		call NumToStr(b,i)
		INQUIRE (file='ghd1int'//b,exist=fil)
		if (fil) then
			open (2,file='ghd1int'//b)
			close (2,status='delete')
		endif
	enddo
!
30	return
	end

!****************************************************************************************
!
!							           DETNA subroutine
!
!	This subroutine determines the NA vector (this locates the diagonal coefficient for
!	skyline storage) and the length of the stiffness vector.
!	This subroutine is called in: CONTROL
!	This subroutine calls: -
!
!	variables used are:
!
!		nnp	: total number of nodes in problem (FE &BE zones)
!		ne8d	: number of drained (dry) elements in FE zone
!		ne8c	: number of consolidated (saturated) elements in FE zone
!		ne8u	: number of unsaturated elements in FE zone
!		nbed1	: number of D1 BEM areas (finite areas), it can be at maximum 5 zones
!		nbed2	: number of D2 BEM areas (infinite areas), it can be at maximum 5 zones
!		nbed3	: number of D3 BEM areas (semi-infinite area), it can be at maximum 1 zone
!		nbcx	: number of (null) boundary conditions introduced for x-displacements
!		nbcy	: number of (null) boundary conditions introduced for y-displacements
!		nbcw	: number of (null) boundary conditions introduced for nodal water pressure
!		nbca	: number of (null) boundary conditions introduced for nodal air pressure
!		ie3d1	: ie3d1(1:nnpbed1,1:nbed): nodes' number which form the current D1 BEM area
!										   Numbering direction for external surface close
!										   (finite) domain is ANTICLOCKWISE
!				  ie3d1(nnpbed1(i)+1,1:nbed): it shows the type of just one material which
!											  is used:
!							[1]: dry boundary elm
!							[2]: saturated boundary elm
!							[3]: unsaturated boundary elm
!		ie3d2	: ie3d2(1:nnpbed2,1:nbed): nodes' number which form the current D2 BEM area
!										   Numbering direction for internal surface open
!										   (infinite) domain is CLOCKWISE
!				  ie3d2(nnpbed2(i)+1,1:nbed): it shows the type of just one material which
!											  is used:
!							[1]: dry boundary elm
!							[2]: saturated boundary elm
!							[3]: unsaturated boundary elm
!		ie3d3	: ie3d3(1:nnpbed3): nodes' number which form the current D3 BEM area
!										   Numbering direction for semi-infinite domain is
!										   ANTICLOCKWISE
!				  ie3d3(nnpbed3+1): it shows the type of just one material which
!											  is used:
!							[1]: dry boundary elm
!							[2]: saturated boundary elm
!							[3]: unsaturated boundary elm
!		ie8d	:
!		ie8c	:
!		ie8u	:
!		id		: degrees of freedom of each node from "1" to "mdof"; for the nodes with
!				  the boundary conditions, each node takes from "mdof+1" to "mdof+nbc"
!				  [X,Y]; which are free for the dry elm &
!				  [X,Y,Pw]; which are free for the saturated elm
!				  [X,Y,Pw,Pa]; which are free for the unsaturated elm
!		mdof	: maximum degree of freedom by considering the boundary conditions
!		kodex	: numbers of nodes in which zero displacement in X-dir is imposed
!		kodey	: numbers of nodes in which zero displacement in Y-dir is imposed
!		kodew	: numbers of nodes in which zero pore water pressure is imposed
!		kodea	: numbers of nodes in which zero pore gas pressure is imposed
!		isolv	: solution type indicator [1,3];
!  								[1]: for a symmetric computation,
!								[3]: for a non-symmetric computation problem; in most of the
!									 case is solved by ISOLV=[3]
!		di		:
!		na		: this locates the diagonal coefficient for skyline storage. we economise
!				  the space by neglecting the stockage of certain coefs which are zero
!				  because of the absence of the connectivity of doF. we delete just the zeros
!				  lied to the non connectivity. therefore na(i) shows the i-th coef to stock.
!		ls		: lenght of the global rigidity matrix in symmetric case. if the matrix is
!				  not symmetric, 2 triangles must be filled.
!		ls1		: if the problem is non symmetric then this partie must be taken into account
!				  then for this partie there is the same number of coefs to stocke like ls.
!
!	INPUT	: ie8d,ie8c,ie8u,ie3d1,ie3d2,ie3d3,kodex,kodey,kodew,kodea,id,mdof
!
!	OUTPUT	: id,di,na,ls,ls1
!
!****************************************************************************************
!
	subroutine detna(ie8d,ie8c,ie8u,id,kodex,kodey,kodew,kodea,di,na,ie3d1,ie3d2,ie3d3)
!
	implicit double precision (a-h,o-z)
!
	common /a/ nnp,ifem,ibem,idyn,ieaq,nr
	common /a1/ ne8d,ne8c,ne8u
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /a3/ ne8d1,ne8c1,ne8u1,ne3u1,nload1,nbcx1,nbcy1,nbcw1,nbca1,m8d1,m8c1,m8u1,&
				nbed11,nbed21,nbed31,nbeint1,ne8dd,ne8cc,ne8uu
	common /d/ nbcx,nbcy,nbcw,nbca,nbc
	common /h/ mdof,mdofn,ls,ls1,isolv
!
	dimension ie8d(9,ne8d1),ie8c(9,ne8c1),ie8u(9,ne8u1),id(4,nnp),kodex(nbcx1),kodey(nbcy1),&
	          kodew(nbcw1),kodea(nbca1),na(mdof),di(mdofn),ie3d1(npbmax,nbed11),&
			  ie3d2(npbmax,nbed21),ie3d3(npbmax)
!
!     allocate locations for fixed X-displacements Ux=0
!
	n=mdof
	if (nbcx.eq.0) goto 175
	do 170 i=1,nbcx
		ii=kodex(i)
		n=n+1
		di(n)=0.0
170		id(1,ii)=n
175	continue
!
!     allocate locations for fixed Y-displacements Uy=0
!
	if (nbcy.eq.0) goto 185
	do 180 i=1,nbcy
		ii=kodey(i)
		n=n+1
		di(n)=0.0
180		id(2,ii)=n
185	continue
!
!     allocate locations for fixed W pressions Pw=0
!
	if (nbcw.eq.0) goto 195
	do 190 i=1,nbcw
		ii=kodew(i)
		n=n+1
		di(n)=0.0
190		id(3,ii)=n
195	continue
!
!     allocate locations for fixed A pressions Pa=0
!
	if (nbca.eq.0) goto 197
	do 196 i=1,nbca
		ii=kodea(i)
		n=n+1
		di(n)=0.0
196		id(4,ii)=n
197	continue
!
!     allocate degrees of freedom
!
	n=0
	do 200 i=1,nnp
		do 200 j=1,4
		if (id(j,i).ne.0) goto 200
		n=n+1
		id(j,i)=n
200	continue
!
	do 210 i=1,mdof
210		na(i)=i
!
	if (ne8d.eq.0) goto 235
		do 230 m=1,ne8d
		imin=mdof
			do 220 i=1,8
			nnod=ie8d(i,m)
			do 220 j=1,4
220			if (imin.gt.id(j,nnod)) imin=id(j,nnod)
		do 230 i=1,8
		nnod=ie8d(i,m)
		do 230 j=1,4
		ndf=id(j,nnod)
		if (ndf.gt.mdof) goto 230
		if (na(ndf).gt.imin) na(ndf)=imin
230		continue
235	continue
!
	if (ne8c.eq.0) goto 260
		do 250 m=1,ne8c
		imin=mdof
			do 240 i=1,8
			nnod=ie8c(i,m)
			do 240 j=1,4
240			if (imin.gt.id(j,nnod)) imin=id(j,nnod)
		do 250 i=1,8
		nnod=ie8c(i,m)
		do 250 j=1,4
		ndf=id(j,nnod)
		if (ndf.gt.mdof) goto 250
		if (na(ndf).gt.imin) na(ndf)=imin
250		continue
260	continue
!
	if (ne8u.eq.0) goto 261
		do 251 m=1,ne8u
		imin=mdof
			do 241 i=1,8
			nnod=ie8u(i,m)
			do 241 j=1,4
241			if (imin.gt.id(j,nnod)) imin=id(j,nnod)
		do 251 i=1,8
		nnod=ie8u(i,m)
		do 251 j=1,4
		ndf=id(j,nnod)
		if (ndf.gt.mdof) goto 251
		if (na(ndf).gt.imin) na(ndf)=imin
251		continue
261	continue
!
	if (nbed1.eq.0) goto 290
	do 285 m=1,nbed1
		nad1=mdof
		do 265 i=1,nnpbed1(m)
			nnod=ie3d1(i,m)
			do 265 j=1,4
				if (nad1.gt.id(j,nnod)) nad1=id(j,nnod)
265		continue
		do 270 i=1,nnpbed1(m)
			nnod=ie3d1(i,m)
			do 270 j=1,4
				ndf=id(j,nnod)
				if (ndf.gt.mdof) goto 270
				if (na(ndf).gt.nad1) na(ndf)=nad1
270		continue
285	continue
290	continue
!
	if (nbed2.eq.0) goto 340
	do 335 m=1,nbed2
		nad2=mdof
		do 315 i=1,nnpbed2(m)
			nnod=ie3d2(i,m)
			do 315 j=1,4
				if (nad2.gt.id(j,nnod)) nad2=id(j,nnod)
315		continue
	do 320 i=1,nnpbed2(m)
		nnod=ie3d2(i,m)
		do 320 j=1,4
			ndf=id(j,nnod)
			if (ndf.gt.mdof) goto 320
			if (na(ndf).gt.nad2) na(ndf)=nad2
320	continue
335	continue
340	continue
!
	if (nbed3.eq.0) goto 355
		nad3=mdof
		do 345 i=1,nnpbed3
		nnod=ie3d3(i)
		do 345 j=1,4
		if (nad3.gt.id(j,nnod)) nad3=id(j,nnod)
345		continue
		do 350 i=1,nnpbed3
		nnod=ie3d3(i)
		do 350 j=1,4
		ndf=id(j,nnod)
		if (ndf.gt.mdof) goto 350
		if (na(ndf).gt.nad3) na(ndf)=nad3
350		continue
355	continue
!
!     calculate na vector
!
	i=1
	do 360 j=2,mdof
		i=i+j-na(j)+1
		na(j)=i
360	continue
	ls=na(mdof)
	ls1=1
	if (isolv.eq.3) ls1=ls
!
	return
	end

!****************************************************************************************
!
!							           DIM1 subroutine
!
!	AIM			: reads basic data and allocates dynamic storage from l1 until l39.
!	callED IN	: CONTROL
!	callS		: -
!
!	VARIABLES	:
!
!		title	: title of the current problem
!
!	in block /stability/, these variables are presented:
!
!		anew1	: Gamma coefficient in Newmark direct integration method. this coef must be
!				  [>=0.5] because that the convergence be inconditionally stable.
!				  In this program, it is chosen 0.5. this choice corresponds to the rule of
!				  the trapezoid
!		anew2	: Beta coefficient in Newmark direct integration method. this coef must be
!				  [>=0.25] because that the convergence be inconditionally stable.
!				  In this program, it is chosen 0.25. this choice corresponds to the rule of
!				  the trapezoid
!		wil		: TETA coefficient in TETA-Wilson method for the stability of numerical
!				  solution [>1]. the coefficient plays a numerical damping role of artificial
!				  wave propagation due to the perturbations. When TETA=[1] it is meaning that
!				  TETA-Wilson method is not considered
!		intb	: interpolation function of time [1,2,3];
!								[1]: constant: two fields of displacement and stress (u,t)
!											   are supposed to be constant in each time
!											   interval
!								[2]: linear	 : two fields of displacement and stress (u,t)
!											   vary linearly in each time interval
!								[3]: mixte	 : displacement field (u) vary linearly while
!											   stress field (t) is constant in each time
!											   interval
!
!	in block /Kcor/, these variables are presented:
!
!		kcorf	: K-correction parameters: K-correction code [0,1]; [0]: wrobel method ;
!				  [1]: time truncated method
!		mtime	: K-correction parameters: cutting point of near history. Time integration
!				  is limited to the number of time steps (MTIME) in time truncation method
!		rmtol	: K-correction parameters: K-correction tolerance
!		c1		:
!		c2		:
!
!	in block /ocqm/, these variables are presented:
!
!		EPScqm	:
!		Icqm	: code of CQM [0,1]
!		Lcqm	:
!		NPcqm	:
!
!		idplst	: plasticity code for dry soil
!		icplst	: plasticity code for saturated soil
!		iuplst	: plasticity code for unsaturated soil
!
!
!	INPUT	:	nnp,ifem,ibem
!				ne8d,ne8c,ne8u
!				nbed1,nbed2,nbed3
!				nnpbed1,nnpbed2,nnpbed3,nnpen
!				m8d,m8c,m8u
!				m8dep,m8cep,m8uep
!				idyn,ieaq,nload,ntime,isolv,nr,dtimei,dtimem,ptoll1,ptoll3
!				kcorf,mtime,rmtol
!				nbcx,nbcy,nbcw,nbca
!				gammaw,gammaa,grav,atmp
!				anew1,anew2,wil,intb
!
!	OUTPUT	:	ne8d1,ne8c1,ne8u1,nbed11,nbed21,nbed31,nload1,nbcx1,nbcy1,nbcw1,nbca1
!				m8d1,m8c1,m8u1,ne8dd,ne8cc,ne8uu,idplst,icplst,iuplst,npbmax,n12gh,n3gh,
!				l1-l38
!
!****************************************************************************************
!
	subroutine dim1(l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14,l15,l16,l17,l18,l19,l20,&
					l21,l22,l23,l24,l25,l26,l27,l28,l29,l30,l31,l32,l33,l34,l35,l36,l37,l38,&
					l39,l40)
!
	implicit double precision (a-h,o-z)
	character (72) :: title, card*5
	logical strcomp
!
	common /a/ nnp,ifem,ibem,idyn,ieaq,nr
	common /a1/ ne8d,ne8c,ne8u
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /a3/ ne8d1,ne8c1,ne8u1,ne3u1,nload1,nbcx1,nbcy1,nbcw1,nbca1,m8d1,m8c1,m8u1,&
				nbed11,nbed21,nbed31,nbeint1,ne8dd,ne8cc,ne8uu
	common /b/ m8d,m8c,m8u,m8dep(10),m8cep(10),m8uep(10)
	common /c/ nload,ntime,itmax,itime,iter,iload,time,dtime,dtimei,dtimem,xinc
	common /d/ nbcx,nbcy,nbcw,nbca,nbc
	common /e/ gammaw,gammas,gammaa,grav,atmp
	common /h/ mdof,mdofn,ls,ls1,isolv
	common /h1/ ptoll1,ptoll3
	common /stability/ anew1,anew2,wil,intb
	common /Kcor/ kcorf,mtime,rmtol,c1,c2
	common /ocqm/ EPScqm,taw,Icqm,Lcqm,NPcqm
	common /comp/ icomp
!
!
	read (4,*) title
	write (7,30) title
!
!-------------------------------
!	    read BASIC data
!-------------------------------
!
	read (4,*) card
	if (.not.strcomp(card,'basic')) then
		write (*,*) '**ERROR**: THIS PROGRAM FIRSTLY readS THE BASIC data'
		stop
	endif
!
	read (4,*) nnp,ifem,ibem
	FEMcontrol: SELECT CASE (ifem)
					CASE (0)
						goto 1
					CASE (1)
						read (4,*) ne8d,ne8c,ne8u
				 end SELECT FEMcontrol
1	continue
	BEMcontrol: SELECT CASE (ibem)
				CASE (0)
					write (*,*) '**WARNING**: UNDER CONSTRUCTION! ACTUALLY, THE BEM IS &
											  NECESSARY TO RUN THE COMPLICATED PROBLEM'
					goto 2
				CASE DEFAULT
					read (4,*) nbed1,nbed2,nbed3,nbeint
					read (4,*) Icqm
					if (Icqm.eq.1) read (4,*) EPScqm,NPcqm
					if (nbed1.gt.0)  read (4,*) (nnpbed1(i),i=1,nbed1)
					if (nbed2.gt.0)  read (4,*) (nnpbed2(i),i=1,nbed2)
					if (nbed3.gt.0)  read (4,*) nnpbed3,nnpen
					if (ibem.eq.2)   read (4,*) kfsol,icomp
 		       end SELECT BEMcontrol
2	continue
!
	read (4,*) m8d,m8c,m8u
!
	read (4,*) idyn,ieaq,nload,ntime,isolv,nr,dtimei,dtimem,ptoll1,ptoll3
!
	read (4,*) kcorf,mtime,rmtol !	this line defines the parameters for K-Correction process
	read (4,*) nbcx,nbcy,nbcw,nbca
	nbc=nbcx+nbcy+nbcw+nbca
	read (4,*) gammaw,gammas,gammaa,grav,atmp
	read (4,*) anew1,anew2,wil,intb
!
!
!------------------------------------------
!	    define parameters for storage
!------------------------------------------
!
	ne8d1=ne8d
    ne8c1=ne8c
	ne8u1=ne8u
	if (ne8d.eq.0) ne8d1=1
	if (ne8c.eq.0) ne8c1=1
	if (ne8u.eq.0) ne8u1=1
    nbed11=nbed1
    nbed21=nbed2
    nbed31=nbed3
	ne3u1 =ne3u
	if (nbed1.eq.0) nbed11=1
	if (nbed2.eq.0) nbed21=1
	if (nbed3.eq.0) nbed31=1
	if (ne3u .eq.0) ne3u1=1
    nload1=nload
    nbcx1=nbcx
    nbcy1=nbcy
    nbcw1=nbcw
	nbca1=nbca
	if (nload.eq.0) nload1=1
	if (nbcx.eq.0) nbcx1=1
	if (nbcy.eq.0) nbcy1=1
	if (nbcw.eq.0) nbcw1=1
	if (nbca.eq.0) nbca1=1
	m8d1=m8d
	m8c1=m8c
	m8u1=m8u
	if (m8d.eq.0) m8d1=1
	if (m8c.eq.0) m8c1=1
	if (m8u.eq.0) m8u1=1
!
	ne8dd=ne8d
	ne8cc=ne8c
	ne8uu=ne8u
	idplst=0
	icplst=0
	iuplst=0
!
	if (m8d.eq.0) goto 5
	do i=1,m8d
		if (m8dep(i).eq.3) idplst=1
	enddo
5	continue
	if	(m8c.eq.0) goto 10
	do i=1,m8c
		if (m8cep(i).eq.3) icplst=1
	enddo
10	continue
	if	(m8u.eq.0) goto 15
	do i=1,m8u
		if (m8uep(i).eq.3) iuplst=1
	enddo
15	continue
	if (idplst.eq.0) ne8dd=1
	if (icplst.eq.0) ne8cc=1
	if (iuplst.eq.0) ne8uu=1
!
	if (ibem.eq.1) kbem=2 ! dry mat. [Ux,Uy]
	if (ibem.eq.2) kbem=3 ! saturated mat. [Ux,Uy,Pw]
	if (ibem.eq.3) kbem=4 ! unsaturated mat. [Ux,Uy,Pw,Pa]
!
	npbmax=1
	if (nbed1.eq.0) goto 20
	do i=1,nbed1
		if (npbmax.lt.nnpbed1(i)) npbmax=nnpbed1(i)
	enddo
20	if (nbed2.eq.0) goto 25
	do i=1,nbed2
		if (npbmax.lt.nnpbed2(i)) npbmax=nnpbed2(i)
	enddo
25	if (npbmax.lt.nnpbed3) npbmax=nnpbed3
	npbmax=npbmax+1
	n12gh=1
	n3gh=1
	if (ibem.eq.1) then		! dry elm. [X,Y]
		n12gh=2*(npbmax-1)
		n3gh=2*(npbmax-1+nnpen)
	endif
	if (ibem.eq.2) then		! saturated elm. [X,Y,Pw]
		n12gh=3*(npbmax-1)
		n3gh=3*(npbmax-1+nnpen)
	endif
	if (ibem.eq.3) then		! unsaturated elm. [X,Y,Pw,Pa]
		n12gh=4*(npbmax-1)
		n3gh=4*(npbmax-1+nnpen)
	endif
!
    nbeint1=nbeint
	if (nbeint.eq.0) nbeint1=1
!
26	continue
!
!
!------------------------------------------
!	    allocate dynamic storage
!------------------------------------------
!
	l1=1						! on stock 1
	l2=l1+(ne8d1*9+1)/2			! (*9), car il faut stocker les num�ros des 8 noeuds qui le constituent et
!									le num�ro du type de mat�riau choisi
!									(le num�ro de l'�l�ment lui-m�me est relatif;
!									il est affect� � l'int�rieur d'une subroutine)
!								  (+1), on arrondit le r�sultat de la division au nombre sup�rieur
!									(le reste de la division n'est pas pris en compte)
!								  (/2), on stock un integer(4 bytes) au lieu d'un real*8
!								  it gives us ie8d
	l3=l2+(ne8c1*9+1)/2			! it gives us ie8c
	l4=l3+(ne8u1*9+1)/2			! it gives us ie8u
	l5=l4+(nbed11*npbmax+1)/2	! it gives us ie3d1;
	l6=l5+(nbed21*npbmax+1)/2	! it gives us ie3d2;
    l7=l6+(nbed31*npbmax+1)/2	! it gives us ie3d3;
	l8=l7+(nnp*4+1)/2			! it gives us id; (*4), (x,y,Pw,Pa)
    l9=l8+(nload1*3+1)/2		! it gives us iconst;
	l10=l9+(nbcx1+1)/2			! it gives us kodex;
    l11=l10+(nbcy1+1)/2			! it gives us kodey;
    l12=l11+(nbcw1+1)/2			! it gives us kodew;
	l13=l12+(nbca1+1)/2			! it gives us kodea;
	l14=l13+nnp*2				! it gives us (x,y);
	l15=l14+nbeint1*2			! it gives us xint;
	l16=l15+m8d1*60				! it gives us sm8d; (*60)
	l17=l16+m8c1*60				! it gives us sm8c; (*60)
	l18=l17+m8u1*60				! it gives us sm8u; (*60)
	l19=l18+ne8d1*5*9			! it gives us sig8d; (*5[sx,sy,sxy,sm,-]*9)
	l20=l19+ne8c1*6*9			! it gives us sig8c; (*6[sx,sy,sxy,sm,-,pw]*9) 9:Gauss point
	l21=l20+ne8u1*7*9			! it gives us sig8u; (*7[sx,sy,sxy,sm,-,suc,pa]*9) 9:Gauss point
	l22=l21+ne8d1*5*9			! it gives us sig8do; (*5*9)
	l23=l22+ne8c1*6*9			! it gives us sig8co; (*6*9)
	l24=l23+ne8u1*7*9			! it gives us sig8uo; (*7*9)
	l25=l24+ne8d1*5*9			! it gives us sig8di; (*5*9)
	l26=l25+ne8c1*6*9			! it gives us sig8ci; (*6*9)
	l27=l26+ne8u1*7*9			! it gives us sig8ui; (*7*9)
	l28=l27+ne8d1*5*9			! it gives us sigmad; (*5*9)
	l29=l28+ne8c1*6*9			! it gives us sigmac; (*6*9)
	l30=l29+ne8u1*7*9			! it gives us sigmau; (*7*9)
	l31=l30+ne8dd*3*3			! it gives us dpad; ???
	l32=l31+ne8cc*3*3			! it gives us dpac; ???
	l33=l32+ne8uu*3*3			! it gives us dpau; ???
	l34=l33+ne8dd*4*10			! it gives us alfad; ???
	l35=l34+ne8cc*4*10			! it gives us alfac; ???
	l36=l35+ne8uu*4*10			! it gives us alfau; ???
	l37=l36+ne8dd*2				! it gives us iysd; ???
	l38=l37+ne8cc*2				! it gives us iysc; ???
	l39=l38+ne8uu*2				! it gives us iysu; ???
	l40=l39+5*500*nbed11			! it gives us sig3u; (*5[sx,sy,sxy,suc,pa])
!
!
30	FORMAT (1x,a72)
	return
	end
!****************************************************************************************
!
!							           DIM2 SUBROUTINE
!
!	This subroutine determines the degrees of freedom (id array) & allocates further
!	dynamic dynamic storage (from l39 until l71)
!	This subroutine is called in: CONTROL
!	This subroutine calls: -
!
!	variables used are:
!
!		nnp		: total number of nodes in problem (FE &BE zones)
!		ne8d	: number of drained (dry) elements in FE zone
!		ne8c	: number of consolidated (saturated) elements in FE zone
!		ne8u	: number of unsaturated elements in FE zone
!		nbed1	: number of D1 BEM areas (finite areas), it can be at maximum 5 zones
!		nbed2	: number of D2 BEM areas (infinite areas), it can be at maximum 5 zones
!		nbed3	: number of D3 BEM areas (semi-infinite area), it can be at maximum 1 zone
!		nnpbed1	: number of nodes of each D1 BEM area; NNPBED1(1 to NBED1)
!		nnpbed2	: number of nodes of each D2 BEM area; NNPBED2(1 to NBED2)
!		nnpbed3	: number of nodes of one D3 BEM area
!		ie3d1	: ie3d1(1:nnpbed1,1:nbed): nodes' number which form the current D1 BEM area
!										   Numbering direction for external surface close
!										   (finite) domain is ANTICLOCKWISE
!				  ie3d1(nnpbed1(i)+1,1:nbed): it shows the type of just one material which
!											  is used:
!							[1]: dry boundary elm
!							[2]: saturated boundary elm
!							[3]: unsaturated boundary elm
!		ie3d2	: ie3d2(1:nnpbed2,1:nbed): nodes' number which form the current D2 BEM area
!										   Numbering direction for internal surface open
!										   (infinite) domain is CLOCKWISE
!				  ie3d2(nnpbed2(i)+1,1:nbed): it shows the type of just one material which
!											  is used:
!							[1]: dry boundary elm
!							[2]: saturated boundary elm
!							[3]: unsaturated boundary elm
!		ie3d3	: ie3d3(1:nnpbed3): nodes' number which form the current D3 BEM area
!										   Numbering direction for semi-infinite domain is
!										   ANTICLOCKWISE
!				  ie3d3(nnpbed3+1): it shows the type of just one material which
!											  is used:
!							[1]: dry boundary elm
!							[2]: saturated boundary elm
!							[3]: unsaturated boundary elm
!		nbc		: total number of null boundary conditions; nbc = nbcx+nbcy+nbcw+nbca
!		npbmax	: maximum nodes in all of the BE zone + 1; for example if we consider D3 BE
!				  area NPBMAX = nnpbed3 + 1
!		n12gh	: in geomechanical problems there is more than one unknown per node, so the
!				  matrices in BEM are assambled according to the DOF, rather than node
!				  number. for dry soil each node has [2] DOF (ux,uy), for for saturated soil
!				  each node has [3] (ux,uy,pw) DOF & for unsaturated soil each node has
!				  [4] DOF(ux,uy,pw,pa). then total number of DOF in BE=(npbmax-1)*kbem
!		n3gh	: in semi-infinite problem in BEM we have to consider a Fictitious enclosing
!				  surface to perform point collocation method. N3GH shows the total number
!				  of DOF in (D3+enclosing surface) in BE technique=(npbmax-1+nnpen)*kbem
!		id		: degrees of freedom of each node from "1" to "mdof"; for the nodes with
!				  the boundary conditions, each node takes from "mdof+1" to "mdof+nbc"
!				  [X,Y]; which are free for the dry elm &
!				  [X,Y,Pw]; which are free for the saturated elm
!				  [X,Y,Pw,Pa]; which are free for the unsaturated elm
!		kbem	: number of degree of freedom per node in BE zone [2,3,4];
!							if IBEM=1 -> KBEM=2 [X,Y]
!							if IBEM=2 -> KBEM=3 [X,Y,Pw]
!							if IBEM=4 -> KBEM=4 [X,Y,Pw,Pa]
!		mdofn	: maximum degree of freedom without considering the boundary conditions
!		mdof	: maximum degree of freedom by considering the boundary conditions
!		ie8d	:
!		ie8c	:
!		ie8u	:
!		l39-l71	: dimensions of dynamic storage

!
!	INPUT	:	nnp,ne8d,ne8c,ne8u,nbed1,nbed2,nbed3,nnpbed1,nnpbed2,nnpbed3,n12gh,n3gh
!
!	OUTPUT	:	id,mdof,mdofn,l39-l71
!
!****************************************************************************************
!
	SUBROUTINE dim2(ie8d,ie8c,ie8u,id,l40,l41,l42,l43,l44,l45,l46,l47,l48,l49,l50,l51,l52,&
			  l53,l54,l55,l56,l57,l58,l59,l60,l61,l62,l63,l64,l65,l66,l67,l68,l69,l70,l71,&
			  l72,l73,ie3d1,ie3d2,ie3d3)
!
	IMPLICIT DOUBLE PRECISION (a-h,o-z)
!
	COMMON /a/ nnp,ifem,ibem,idyn,ieaq,nr
	COMMON /a1/ ne8d,ne8c,ne8u
	COMMON /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	COMMON /a3/ ne8d1,ne8c1,ne8u1,ne3u1,nload1,nbcx1,nbcy1,nbcw1,nbca1,m8d1,m8c1,m8u1,&
				nbed11,nbed21,nbed31,nbeint1,ne8dd,ne8cc,ne8uu
	COMMON /b/ m8d,m8c,m8u,m8dep(10),m8cep(10),m8uep(10)
	COMMON /c/ nload,ntime,itmax,itime,iter,iload,time,dtime,dtimei,dtimem,xinc
	COMMON /d/ nbcx,nbcy,nbcw,nbca,nbc
	COMMON /h/ mdof,mdofn,ls,ls1,isolv
!
	DIMENSION ie8d(9,ne8d1), ie8c(9,ne8c1), ie8u(9,ne8u1),id(4,nnp),ie3d1(npbmax,nbed11),&
	          ie3d2(npbmax,nbed21),ie3d3(npbmax)
!
!
	kk=4*nnp+1		! because for each node,at max, there are four dof [X,Y,Pw,Pa]
!					  for dry elm:			[X,Y]
!					  for saturated elm:	[X,Y,Pw]
!					  for unsaturated elm:	[X,Y,Pw,Pa]
!
	DO 100 i=1,nnp
		DO 100 j=1,4
100			id(j,i)=kk
!					  in the first stage, the ID matrix was constructed
!
	IF (ne8d.EQ.0) GO TO 350
	DO 300 i=1,ne8d
		DO 200 j=1,8
			nnod=ie8d(j,i)
			DO 200 k=1,2
!
200			   id(k,nnod)=0
300	CONTINUE
!     we work on the 8 first lines of the IE8D (i.e. on the node numbers of the elms)
!	  we work on the 2 first lines of the mat. of nodes [X,Y]. because we are in drained cond.
!	  in which in this cond. there is no influence on Pw & Pa which are fixed.
!	  IDEA: in each node in a drianed elm [X,Y] are free (libre) therefore, if there is no b.c.
!		    there are the degrees of freedom. to indicate it, we annule the corresponding ID
!		    coeff.
!
350	IF (ne8c.EQ.0) GO TO 600
	DO 500 i=1,ne8c
		DO 400 j=1,8
			nnod=ie8c(j,i)
			DO 400 k=1,3
400				id(k,nnod)=0
500	CONTINUE
!	  IDEA: in each node in a saturated elm [X,Y,Pw] are free (libre) therefore, if there is
!		    no b.c. there are the degrees of freedom. to indicate it, we annule the
!			corresponding ID coeff.
!
600	IF (ne8u.EQ.0) GO TO 700
	DO 800 i=1,ne8u
		DO 750 j=1,8
			nnod=ie8u(j,i)
			DO 750 k=1,4
750				id(k,nnod)=0
800	CONTINUE
!	  IDEA: in each node in an unsaturated elm [X,Y,Pw,Pa] are free (libre) therefore,
!			if there is no b.c., there are the degrees of freedom. to indicate it,
!			we annule the corresponding ID coeff.
!
700	IF (nbed1.EQ.0) GO TO 900
	DO i=1,nbed1
		DO j=1,nnpbed1(i)
			nnod=ie3d1(j,i)
			DO k=1,kbem
				id(k,nnod)=0
			END DO
		END DO
	END DO
900	CONTINUE
!
	IF (nbed2.EQ.0) GO TO 1200
	DO i=1,nbed2
		DO j=1,nnpbed2(i)
			nnod=ie3d2(j,i)
			DO k=1,kbem
				id(k,nnod)=0
			END DO
		END DO
	END DO
1200	CONTINUE
!
	IF (nbed3.EQ.0) GO TO 1500
	DO j=1,nnpbed3
		nnod=ie3d3(j)
		DO k=1,kbem
			id(k,nnod)=0
		END DO
	END DO
1500	CONTINUE
!
!	  calculate maximum degree of freedom
	mdof=0
	DO i=1,nnp
		DO j=1,4
			IF (id(j,i).EQ.0) mdof=mdof+1
		END DO
	END DO
	mdofn=mdof
	mdof=mdofn-nbc	! we must remove the boundary conditions. therefore we obtain the
!			          number of DoF.
!
	l41=l40+mdof
	l42=l41+mdof
    l43=l42+mdof
	l44=l43+mdof
	l45=l44+mdofn
	l46=l45+mdofn
	l47=l46+mdofn
	l48=l47+mdofn
	l49=l48+mdofn
	l50=l49+mdofn
	l51=l50+(mdof+1)/2
	l52=l51+(nnp+1)/2
	l53=l52+(n12gh*n12gh)*nbed11
	l54=l53+(n12gh*n12gh)*nbed21
	l55=l54+(n12gh*n12gh)
	l56=l55+(n12gh*n12gh)*nbed11
	l57=l56+(n12gh*n12gh)*nbed21
	l58=l57+(n12gh*n12gh)
	l59=l58+(n12gh*n12gh)*nbed11
	l60=l59+(n12gh*n12gh)*nbed21
	l61=l60+(n12gh*n12gh)
	l62=l61+(n12gh*n12gh)*nbed11
	l63=l62+(n12gh*n12gh)*nbed21
	l64=l63+(n12gh*n12gh)
	l65=l64+n12gh*ntime*nbed11
	l66=l65+n12gh*ntime*nbed21
	l67=l66+n12gh*ntime
	l68=l67+n12gh*ntime*nbed11
	l69=l68+n12gh*ntime*nbed21
	l70=l69+n12gh*ntime
	l71=l70+n12gh*ntime*nbed11
	l72=l71+n12gh*ntime*nbed21
	l73=l72+n12gh*ntime
!
	RETURN
	END

!****************************************************************************************
!
!							            DINVERSE subroutine
!
!	This subroutine calculates the inverse of a matrix Gd3**(-1).
!	This subroutine is called in: GHMATD3
!	This subroutine calls: DGECO,DGEFA,DGEDI
!		dgeco	: factors a double precision matrix by gaussian elimination and estimates
!				  the condition of the matrix.
!		dgefa	: factors a double precision matrix by gaussian elimination. dgefa is
!				  usually called by dgeco, but it can be called directly with a saving in
!				  time if  rcond  is not needed.
!				  (time for dgeco) = (1 + 9/n)*(time for dgefa) .
!		dgedi	: computes the inverse of a matrix using the factors computed by dgeco or dgefa
!
!	variables used are:
!
!		a		: matrix to inverse = gdinv
!		lda		: the leading dimension of the array a = n12gh
!		n		: the order of the matrix  a = n1
!		info	: [0,-1,k]
!						  [0] : normal value
!						  [-1]: bad condition of matrix : u(k,k) < e-8
!						  [k] : singularity (divise by zero):  u(k,k) =  0.0
!
!
!	INPUT	:
!
!	OUTPUT	: a (inverse of original matrix)
!
!
!****************************************************************************************
!
	subroutine dinverse(a,lda,n,info)
!
	integer lda,n,ipvt(n),info
	double precision a(lda,1),work(n)
!
!	call dgeco(a,lda,n,ipvt,info,z)
	call dgefa(a,lda,n,ipvt,info)
	if (info.gt.0) return
	call dgedi (a,lda,n,ipvt,work)
!
	return
	end
!
!
!****************************************************************************************
!							            DGEFA subroutine
!	This subroutine factors a double precision matrix by gaussian elimination. dgefa is
!	usually called by dgeco, but it can be called directly with a saving in time if  rcond
!	is not needed.
!	This subroutine is called in: DINVERSE
!	This subroutine calls: DSCAL,DAXPY
!	variables used are:
!		a		: an upper triangular matrix and the multipliers which were used to
!				  obtain it. the factorization can be written  a = l*u  where l  is a
!				  product of permutation and unit lower triangular matrices and  u  is
!				  upper triangular. = gdinv
!		lda		: the leading dimension of the array  a = n12gh
!		ipvt	: an integer vector of pivot indices. number of line which has the
!				  max. absolute value in column k.
!		info	: [0,-1,k]
!						  [0] : normal value
!						  [-1]: bad condition of matrix : u(k,k) < e-8
!						  [k] : singularity (divise by zero):  u(k,k) =  0.0
!				  this is not an error condition for this subroutine, but it does
!				  indicate that dgesl or dgedi will divide by zero if called.  use  rcond
!				  in dgeco for a reliable indication of singularity.
!****************************************************************************************
!
	subroutine dgefa(a,lda,n,ipvt,info)
!
	integer lda,n,ipvt(1),info
	double precision a(lda,1)
!
!	internal variables
	double precision t
	integer idamax,j,k,kp1,l,nm1
!
!	gaussian elimination with partial pivoting
!
	info = 0
	nm1 = n - 1
	if (nm1 .lt. 1) goto 70
	do 60 k = 1, nm1	! iterate over columns
		 kp1 = k + 1
!
!	find l = pivot index
		l = idamax(n-k+1,a(k,k),1) + k - 1	! n-k+1: number of composents in each column
!				         					  under the diameter + the diagonal comps to
!											  study the max. value in each column
		ipvt(k) = l
!
!	zero pivot implies this column already triangularized
!
		if (dabs(a(l,k)) .lt. 1.d-8) info = -1	! if the max. absolute value in column k
!												  is equal to 1.d-8 then info = -1
		if (a(l,k) .eq. 0.0d0) goto 40			! no unique solution exists, stop
!
!	interchange if necessary
		if (l .eq. k) goto 10
		t = a(l,k)
		a(l,k) = a(k,k)
		a(k,k) = t
10		continue
!
!	compute multipliers
		t = -1.0d0/a(k,k)
		call dscal(n-k,t,a(k+1,k),1)
!
!	row elimination with column indexing
		do 30 j = kp1, n
			t = a(l,j)
			if (l .eq. k) goto 20
			a(l,j) = a(k,j)
			a(k,j) = t
20			continue
			call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
30		continue
		goto 50
40		continue
		info = k
50		continue
60	continue
70	continue
!
	ipvt(n) = n
	if (dabs(a(n,n)) .lt. 1.d-8) info = -1
	if (a(n,n) .eq. 0.0d0) info = n
	return
	end
!
!	--------------------------------------------------------------------------
!
!	finds the index of element having max. absolute value.
!
!	n		: number of elements studied in each column to find the max. value
!			  which is the number of elements under the diameter + 1 (diameter)
!	dx		: diagonal element a(k,k)
!	incx	: increment ???????????
!	dmax	: max. absolute value in each column
!	idamax	: number of line which has the max. absolute value in column
!
	integer function idamax(n,dx,incx)
	double precision dx(1),dmax
	integer i,incx,ix,n
!
	idamax = 0
	if ( n .lt. 1 ) return
	idamax = 1
	if (n.eq.1) return
	if (incx.eq.1) goto 20
!
!	code for increment not equal to 1
	ix = 1
	dmax = dabs(dx(1))
	ix = ix + incx
	do 10 i = 2,n
		if (dabs(dx(ix)).le.dmax) goto 5
		idamax = i
		dmax = dabs(dx(ix))
5		ix = ix + incx
10	continue
	return
!
!	code for increment equal to 1
20	dmax = dabs(dx(1))
	do 30 i = 2,n
		if (dabs(dx(i)).le.dmax) goto 30
		idamax = i
		dmax = dabs(dx(i))
30	continue
	return
	end
!
!
!****************************************************************************************
!							            DGEDI subroutine
!	This subroutine computes inverse(u)
!	This subroutine is called in: DINVERSE
!	This subroutine calls: DSCAL,DAXPY,DSWAP
!****************************************************************************************
!
	subroutine dgedi(a,lda,n,ipvt,work)
!
	integer lda,n,ipvt(1)
	double precision a(lda,1),work(1)
!
!	internal variables
	double precision t
	integer i,j,k,kb,kp1,l,nm1
!
!	compute inverse(u)
!
	do 100 k = 1, n
		a(k,k) = 1.0d0/a(k,k)
		t = -a(k,k)
		call dscal(k-1,t,a(1,k),1)
		kp1 = k + 1
		if (n .lt. kp1) goto 90
		do 80 j = kp1, n
			t = a(k,j)
            a(k,j) = 0.0d0
			call daxpy(k,t,a(1,k),1,a(1,j),1)
80		continue
90		continue
100	continue
!
!	form inverse(u)*inverse(l)
!
	nm1 = n - 1
	if (nm1 .lt. 1) goto 140
	do 130 kb = 1, nm1
		k = n - kb
		kp1 = k + 1
		do 110 i = kp1, n
			work(i) = a(i,k)
            a(i,k) = 0.0d0
110		continue
		do 120 j = kp1, n
			t = work(j)
			call daxpy(n,t,a(1,j),1,a(1,k),1)
120		continue
		l = ipvt(k)
		if (l .ne. k) call dswap(n,a(1,k),1,a(1,l),1)
130	continue
140	continue
150	continue
	return
	end
!
!
!****************************************************************************************
!							            DSCAL subroutine
!	This subroutine scales a vector by a constant. it uses unrolled loops for increment
!	equal to one. in this subroutine
!	This subroutine is called in: DGEFA
!	variables used are:
!		n		: number of elements studied under the diameter in each column
!		da		: -1/a(k,k)
!		dx		: a(k+1,k)
!		incx	: increment xhich is equal to 1 ?????????
!****************************************************************************************
!
	subroutine dscal(n,da,dx,incx)
!
	double precision da,dx(1)
	integer i,incx,m,mp1,n,nincx
!
	if (n.le.0) return
	if (incx.eq.1) goto 20
!
!	code for increment not equal to 1
	nincx = n*incx
	do 10 i = 1,nincx,incx
		dx(i) = da*dx(i)
10	continue
	return
!
!	code for increment equal to 1
!	clean-up loop
!
20	m = MOD (n,5)
	if ( m .eq. 0 ) goto 40
	do 30 i = 1,m
		dx(i) = da*dx(i)
30	continue
	if ( n .lt. 5 ) return
40	mp1 = m + 1
	do 50 i = mp1,n,5
		dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
50	continue
	return
	end
!
!
!****************************************************************************************
!							            DAXPY subroutine
!	This subroutine constant times a vector plus a vector.
!	This subroutine is called in: DGEFA
!****************************************************************************************
!
	subroutine daxpy(n,da,dx,incx,dy,incy)
!
	double precision dx(1),dy(1),da
	integer i,incx,incy,m,mp1,n
!
	if (n.le.0) return
	if (da .eq. 0.0d0) return
	if (incx.eq.1.and.incy.eq.1) goto 20
!
!	code for unequal increments or equal increments not equal to 1
	ix = 1
	iy = 1
	if (incx.lt.0) ix = (-n+1)*incx + 1
	if (incy.lt.0) iy = (-n+1)*incy + 1
	do 10 i = 1,n
		dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
10	continue
	return
!
!	code for both increments equal to 1
!	clean-up loop
20	m = MOD (n,4)
	if ( m .eq. 0 ) goto 40
	do 30 i = 1,m
		dy(i) = dy(i) + da*dx(i)
30	continue
	if ( n .lt. 4 ) return
40	mp1 = m + 1
	do 50 i = mp1,n,4
		dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
50	continue
	return
	end
!
!
!****************************************************************************************
!							            DSWAP subroutine
!	This subroutine interchanges two vectors.
!	This subroutine is called in: DEGEDI
!****************************************************************************************
!
	subroutine dswap(n,dx,incx,dy,incy)
!
	double precision dx(1),dy(1),dtemp
	integer i,incx,incy,ix,iy,m,mp1,n
!
	if (n.le.0) return
	if (incx.eq.1.and.incy.eq.1) goto 20
!
!	code for unequal increments or equal increments not equal to 1
	ix = 1
	iy = 1
	if (incx.lt.0) ix = (-n+1)*incx + 1
	if (incy.lt.0) iy = (-n+1)*incy + 1
	do 10 i = 1,n
		dtemp = dx(ix)
        dx(ix) = dy(iy)
        dy(iy) = dtemp
        ix = ix + incx
        iy = iy + incy
10	continue
	return
!
!	code for both increments equal to 1
!	clean-up loop
20	m = MOD (n,3)
	if ( m .eq. 0 ) goto 40
	do 30 i = 1,m
		dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
30	continue
	if ( n .lt. 3 ) return
40	mp1 = m + 1
	do 50 i = mp1,n,3
		dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
        dtemp = dx(i + 1)
        dx(i + 1) = dy(i + 1)
        dy(i + 1) = dtemp
        dtemp = dx(i + 2)
        dx(i + 2) = dy(i + 2)
        dy(i + 2) = dtemp
50	continue
	return
	end
!
!	--------------------------------------------------------------------
!
!	takes the sum of the absolute values.
!
	double precision function dasum(n,dx,incx)
	double precision dx(1),dtemp
	integer i,incx,m,mp1,n,nincx
!
	dasum = 0.0d0
	dtemp = 0.0d0
	if (n.le.0) return
	if (incx.eq.1) goto 20
!
!	code for increment not equal to 1
	nincx = n*incx
	do 10 i = 1,nincx,incx
		dtemp = dtemp + dabs (dx(i))
10	continue
	dasum = dtemp
	return
!
!	code for increment equal to 1
!	clean-up loop
!
20	m = MOD (n,6)
	if ( m .eq. 0 ) goto 40
	do 30 i = 1,m
		dtemp = dtemp + dabs (dx(i))
30	continue
	if ( n .lt. 6 ) goto 60
40	mp1 = m + 1
	do 50 i = mp1,n,6
		dtemp = dtemp + dabs(dx(i)) + dabs(dx(i + 1)) + dabs(dx(i + 2)) + dabs(dx(i + 3)) +&
				dabs(dx(i + 4)) + dabs(dx(i + 5))
50	continue
60	dasum = dtemp
	return
	end
!
!	---------------------------------------------------------------------
!
!	forms the dot product of two vectors.
!
	double precision function ddot(n,dx,incx,dy,incy)
	double precision dx(1),dy(1),dtemp
	integer i,incx,incy,ix,iy,m,mp1,n
!
	ddot = 0.0d0
	dtemp = 0.0d0
	if (n.le.0) return
	if (incx.eq.1.and.incy.eq.1) goto 20
!
!	code for unequal increments or equal increments not equal to 1
	ix = 1
	iy = 1
	if (incx.lt.0) ix = (-n+1)*incx + 1
	if (incy.lt.0) iy = (-n+1)*incy + 1
	do 10 i = 1,n
		dtemp = dtemp + dx(ix)*dy(iy)
		ix = ix + incx
		iy = iy + incy
10	continue
	ddot = dtemp
	return
!
!	code for both increments equal to 1
!	clean-up loop
20	m = MOD(n,5)
	if ( m .eq. 0 ) goto 40
	do 30 i = 1,m
		dtemp = dtemp + dx(i)*dy(i)
30	continue
	if ( n .lt. 5 ) goto 60
40	mp1 = m + 1
	do 50 i = mp1,n,5
		dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) + dx(i + 2)*dy(i + 2) +&
		dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
50	continue
60	ddot = dtemp
	return
	end

!!c ****************************************************************

!!c CALCUL NUMERIQUE DE LA TRANSFORMATION LAPLACE INVERSE
!!c using the durbin formula in combination with the epsilon algorithm

!!c ****************************************************************

      subroutine dlainv(fun,t,c,epsre,epsab,max,rslt)

!!c inversion of laplace transform using the durbin formula
!!c in combination with the epsilon algorithm

!!c input parameters
!!!!cfun    - complex double precision
!!!!!c       the actual name for fun needs to be declared
!!!!!c       external in the driver program.

!!!!ct      - double precision
!!!!!c       value of the independent variable for which the
!!!!!c       inverse laplace transform has to be computed.
!!!!!c       t should be greater than zero.

!!!!c!!!c- double precision
!!!!!c       abscissa of convergence of the laplace transform

!!!!cepsre  - double precision
!!!!!c       relative accuracy requested

!!!!cepsab  - double precision
!!!!!c       absolute accuracy requested.
!!!!!c       the routine tries to satisfy the least stringent
!!!!!c       of both accuracy requirement.
!!!!cmax    - bound on the number of terms used in the durbin formula

!!c output parameters
!!!!crslt - double precision
!!!!!c       inverse laplace transform

!!!!cesterr - double precision
!!!!!c       estimate of the absolute error abs(f(t)-rslt)

!!!!cnum    - integer
!!!!!c       number of evaluations of fun

!!!!cier    - integer
!!!!!c       parameter giving information on the termination
!!!!!c       of the algorithm
!!!!!c       ier = 0 normal and reliable termination of the
!!!!!c               routine
!!!!!c       ier = 1 the computations are terminated because
!!!!!c               the bound on the number of evaluations
!!!!!c               of fun has been achieved.  this bound
!!!!!c               is equal to 8*max+5 where  max  is a
!!!!!c               number initialized in a data
!!!!!c               statement.  one can allow more function
!!!!!c               evaluations by increasing the value of
!!!!!c               max in the data-statement.
!!!!!c       ier = 2 the value of t is less than or equal
!!!!!c               to zero.

      double precision aim,ak,are,arg,bb,c,dabs,datan,dexp,dmax1,dsin,epsab,epsre,esterr,fim,fre,pid16,r,rslt,res3la,rex,si,t
      integer i,ier,k,kc,kk,ks,m,nex,nres,num
      double complex s,ff,fun
      dimension si(32),res3la(3),rex(52)

!!c the array si contains values of the sine and cosine functions
!!c required in the durbin formula.
      data si(8),si(16)/ 1.0d+00,0.0d+00/
!!c pid16 is equal to pi/16
      pid16 = datan(1.0d+00)/4.0d+00
      ak = 1.0d+00
      do 10 k=1,7
        si(k) = dsin(ak*pid16)
        ak = ak+1.0d+00
        kk = 16-k
        si(kk) = si(k)
10    continue
      do 20 k=17,32
        si(k) = -si(k-16)
20    continue

!!c test on validity of the input parameter t

      ier = 2
      rslt = 0.0d+00
      esterr = 1.0d+00
      num = 0
      if (t.le.0.0d+00) go to 999
      ier = 0
      nres = 0

!!c initialization of the summation of the durbin formula.

      arg = pid16/t
      are = c+2.0d+00/t
      aim = 0.0d+00
      bb = dexp(are*t)/(1.6d+01*t)
      s=dcmplx(are,aim)
      ff=fun(s)
      fre=dble(ff)
      fim=dimag(ff)
      num = 5
      r = 5.0d-01*fre
      nex = 0
      kc = 8
      ks = 0

!!c main loop for the summation

      do 40 i=1,max
        m = 8
        if (i.eq.1) m = 12
        do 30 k=1,m
          aim = aim+arg
          kc = kc+1
          ks = ks+1
          if (kc.gt.32) kc = 1
          if (ks.gt.32) ks = 1
          s=dcmplx(are,aim)
          ff=fun(s)
          fre=dble(ff)
          fim=dimag(ff)
          r = r+fre*si(kc)-fim*si(ks)
30      continue
        num = num+8
        nex = nex+1
        rex(nex) = r

!!c extrapolation using the epsilon algorithm

        if(nex.ge.3) call dqext(nex,rex,rslt,esterr,res3la,nres)
        if(nres.lt.4) go to 40

!!c computation of intermediate result and estimate of the absolute error

        rslt = rslt * bb
        esterr = esterr * bb
        if (esterr.lt.dmax1(epsab,epsre*dabs(rslt)).and.dabs(r*bb-rslt).lt.5.0d-01*dabs(rslt)) go to 999
40    continue

!!c set error flag in the case that the number of terms in the
!!c summation is equal to max

      ier = 1
!!!cwrite(*,*) 'DLAINV: bound on the number of evaluations has been
!!c $achieved'
999   return
      end

      subroutine dqext(n,epstab,rslt,abserr,res3la,nres)

!!c epsilon algorithm

!!c the routine determines the limit of a given sequence of
!!c approximations, by means of the epsilon algorithm
!!c of p. wynn.
!!c an estimate of the absolute error is also given.
!!c the condensed epsilon table is computed. only those
!!c elements needed for the computation of the next diagonal
!!c are preserved.


	  double precision abserr,dabs,delta1,delta2,delta3,dmax1,epmach,epsinf,epstab,error,err1,err2,err3,e0,e1,e1abs,e2,e3, &
	  				    oflow,res,rslt,res3la,ss,tol1,tol2,tol3,d1mach
      integer i,ib,ib2,ie,indx,k1,k2,k3,limexp,n,newelm,nres,num
      dimension epstab(52),res3la(3)
!c
!!c list of major variables
!!c -----------------------
!!c e0     - the 4 elements on which the
!!c e1       computation of a new element in
!!c e2       the epsilon table is based
!!c e3                 e0
!!!!!c        e3    e1    new
!!!!!c              e2
!!c newelm - number of elements to be computed in the new
!!!!!c    diagonal
!!c error  - error = abs(e1-e0)+abs(e2-e1)+abs(new-e2)
!!c rslt - the element in the new diagonal with least value
!!!!!c    of error

      oflow = d1mach(2)
      epmach = d1mach(4)
      nres = nres+1
      abserr = oflow
      rslt = epstab(n)
      if(n.lt.3) go to 100
      limexp = 50
      epstab(n+2) = epstab(n)
      newelm = (n-1)/2
      epstab(n) = oflow
      num = n
      k1 = n
      do 40 i = 1,newelm
        k2 = k1-1
        k3 = k1-2
        res = epstab(k1+2)
        e0 = epstab(k3)
        e1 = epstab(k2)
        e2 = res
        e1abs = dabs(e1)
        delta2 = e2-e1
        err2 = dabs(delta2)
        tol2 = dmax1(dabs(e2),e1abs)*epmach
        delta3 = e1-e0
        err3 = dabs(delta3)
        tol3 = dmax1(e1abs,dabs(e0))*epmach
        if(err2.gt.tol2.or.err3.gt.tol3) go to 10

!!c if e0, e1 and e2 are equal to within machine
!!c accuracy, convergence is assumed.
!!c rslt = e2
!!c abserr = abs(e1-e0)+abs(e2-e1)

        rslt = res
        abserr = err2+err3
!!c jump out of do-loop
        go to 100
   10   e3 = epstab(k1)
        epstab(k1) = e1
        delta1 = e1-e3
        err1 = dabs(delta1)
        tol1 = dmax1(e1abs,dabs(e3))*epmach

!!c if two elements are very close to each other, omit
!!c a part of the table by adjusting the value of n

        if(err1.le.tol1.or.err2.le.tol2.or.err3.le.tol3) go to 20
        ss = 1.0d+00/delta1+1.0d+00/delta2-1.0d+00/delta3
        epsinf = dabs(ss*e1)

!!c test to detect irregular behaviour in the table, and
!!c eventually omit a part of the table adjusting the value of n.

        if(epsinf.gt.1.0d-04) go to 30
   20   n = i+i-1
!!c jump out of do-loop
        go to 50

!!c compute a new element and eventually adjust the value of rslt.

   30   res = e1+1.0d+00/ss
        epstab(k1) = res
        k1 = k1-2
        error = err2+dabs(res-e2)+err3
        if(error.gt.abserr) go to 40
        abserr = error
        rslt = res
   40 continue

!!c shift the table.

   50 if(n.eq.limexp) n = 2*(limexp/2)-1
      ib = 1
      if((num/2)*2.eq.num) ib = 2
      ie = newelm+1
      do 60 i=1,ie
        ib2 = ib+2
        epstab(ib) = epstab(ib2)
        ib = ib2
   60 continue
      if(num.eq.n) go to 80
      indx = num-n+1
      do 70 i = 1,n
        epstab(i)= epstab(indx)
        indx = indx+1
   70 continue
   80 if(nres.ge.4) go to 90
      res3la(nres) = rslt
      abserr = oflow
      go to 100

!!c compute error estimate

   90 abserr = dabs(rslt-res3la(3))+dabs(rslt-res3la(2))+dabs(rslt-res3la(1))
      res3la(1) = res3la(2)
      res3la(2) = res3la(3)
      res3la(3) = rslt
  100 abserr = dmax1(abserr,5.0d+00*epmach*dabs(rslt))
      return
      end

!!!!!!c*********************************************************************
!!c double precision MACHINE CONSTANTS
!!!!!!c*********************************************************************

      double precision function d1mach(i)
      integer i

!!c d1mach( 1) = b**(emin-1), the smallest positive magnitude
!!c d1mach( 2) = b**emax*(1 - b**(-t)), the largest magnitude
!!c d1mach( 3) = b**(-t), the smallest relative spacing
!!c d1mach( 4) = b**(1-t), the largest relative spacing.
!!c d1mach( 5) = log10(b)

!!!cif (i.eq.1) d1mach = 2.225073858507201e-308
!!!cif (i.eq.2) d1mach = 1.797693134862316e+308
!!!cif (i.eq.3) d1mach = 2.22044049250313e-16
!!!cif (i.eq.4) d1mach = 2*2.22044049250313e-16
      if (i.eq.1) d1mach = tiny(d1mach)
      if (i.eq.2) d1mach = huge(d1mach)
      if (i.eq.3) d1mach = epsilon(d1mach)
      if (i.eq.4) d1mach = 2*epsilon(d1mach)

      return
      end
!****************************************************************************************
!
!							           DMATU subroutine
!
!	This subroutine calculates the D matrix for unsat. soil and aslo suction matrix
!	This subroutine is called in: GHMATD1
!	This subroutine calls: -
!
!	variables used are:
!				eload	: (El) elastic modulus in absence of suction Et=El+Es
!				evkb	: (K0) bulk modulud in absence of suction (initial value of B)
!				ten		:
!				emin	: minimum value of Young modulus
!				sati	: initial saturation degree
!				bmkw	: Cw=d(row)/row/d(pw), compressibility of water
!				perwa	: aw(m.s^-1) used in the formula of water permeability
!				perwal	: alpha_w used in the formula of water permeability
!				sru		: used in the formula of water permeability
!				perwm	: kw_max, maximum authorized water permeability
!				bmka	: Ca=d(roa)/roa/d(pa), compressibility of air
!				perab	: ba used in the formula of air permeability
!				peraal	: alpfa_a used in the formula of air permeability
!				xmua	: viscosity of air
!				peram	: ka_max, maximum authorized air permeability
!				evae	: ae used in the formula of the void ratio state surface
!				evbe	: be used in the formula of the void ratio state surface
!				evce	: evkb(adim.)*P_atm
!				eve0	: e0 initial void ratio, used in the formula of the void ratio
!						  state surface
!				evsigb	: maximum traction resistance
!				srbetaw	: beta_w sed in the formula of the saturation degree state surface
!				sucm	: ms used in the formula of Es=ms(Pa-Pw)
!				henry	: Henry coefficient
!
!	INPUT	:
!
!	OUTPUT	:
!
!****************************************************************************************
!
	subroutine dmatuelas(mtype,sige,sm8u)
!
	implicit double precision (a-h,o-z)
	common /a/ nnp,ifem,ibem,idyn,ieaq,nr
	common /a3/ ne8d1,ne8c1,ne8u1,nload1,nbcx1,nbcy1,nbcw1,nbca1,m8d1,m8c1,m8u1,nbed11,&
				nbed21,nbed31,ne8dd,ne8cc,ne8uu
	common /c/ nload,ntime,itmax,itime,iter,iload,time,dtime,dtimei,dtimem,xinc
	common /e/ gammaw,gammas,gammaa,grav,atmp
	common /umate/ d(3,3),dfs(3),xn,xkw,cww,cwg,xka,cgg,cgw,ev,sat,g1,rs,rw,raa,rmix,xm1,&
				   bt,et,xlambda,xmu,henry
!	common /Dumate/ Dxlambda,Dxmu,Drw,Draa,Dgammaw,Dgammaa,Dxkw,Dxka
!
	dimension ds(3),sm8u(60,m8u1),sige(5)
!
!
	sx =sige(1)
	sy =sige(2)
	sxy=sige(3)
	suc=sige(4)
	pg =sige(5)
!
	if (itime.eq.1) sy=pg
!
!	elastic data 1
!
	eload=sm8u(3,mtype)
	evkb =sm8u(7,mtype)
	evxm =sm8u(8,mtype)
!
!	elastic data 2
!
	ten =sm8u(10,mtype)
	emin=sm8u(11,mtype)
	sati=sm8u(14,mtype)
!
!	water flow data
!
	bmkw=sm8u(6,mtype)
	perwa=sm8u(15,mtype)
	perwal=sm8u(16,mtype)
	sru=sm8u(17,mtype)
	perwm=sm8u(22,mtype)
!
!	air flow data
!
	bmka=sm8u(5,mtype)
	perab=sm8u(19,mtype)
	peraal=sm8u(20,mtype)
	xmua=sm8u(21,mtype)
	peram=sm8u(23,mtype)
	henry=sm8u(34,mtype)
!
!	state surface data
!
	evae=sm8u(24,mtype)
	evbe=sm8u(25,mtype)
	evce=evkb*atmp
	eve0=sm8u(27,mtype)
	evsigb=sm8u(28,mtype)
!
	srbetaw=sm8u(31,mtype)
!
!	suction data
!
	sucm=sm8u(33,mtype)
!
10	continue
!
!	call subroutines to calculate void ratio,saturation,gradients and permeabilities
!
	if (suc.eq.0.d0) suc=dabs((1.d0-sati)/srbetaw)
!
	sx=sx-pg
	sy=sy-pg
	call STATEVT (sy,suc,evae,evbe,evce,eve0,evxm,evsigb)
	call STATEST (sy,suc,srbetaw,sati)
	call perwt	 (perwa,perwal,sru,perwm)
	call PERAT	 (perab,peraal,xmua,sru,peram)
!
!	check for tensile failure
!	check for shear failure
!
	et=eload*atmp+suc*sucm
	btmax=17.d0*et
    btmin=0.5d0*et
	if (bt.lt.btmin) bt=btmin
	if (bt.gt.btmax) bt=btmax
!
!	elasticity matrix
!
	d1=3.d0*bt/(9.d0*bt-et)
	d2=3.d0*bt+et
	d3=3.d0*bt-et
	d(1,1)=d1*d2
	d(2,2)=d(1,1)
	d(1,2)=d1*d3
	d(2,1)=d(1,2)
	d(1,3)=0.d0
	d(2,3)=0.d0
	d(3,1)=0.d0
	d(3,2)=0.d0
	d(3,3)=d1*et
!
!	form Ds^-1 & Fs
!
	do i=1,2
		ds(i)=xm1
	enddo
	ds(3)=0.d0
	do i=1,3
		temp=0.d0
		do j=1,3
			temp=temp+d(i,j)*ds(j)
		enddo
		dfs(i)=temp
	enddo
!
!	solid skeleton parameters
!
	xlambda=d1*d3
	xmu=d1*et
	if (xlambda.lt.0.D0.OR.xmu.lt.0.D0) then
		write (*,*) 'ERROR! LAME COEF IS NEGATIVE'
		stop
	endif
	xn=ev/(1.d0+ev)
	rs=gammas/grav
	rw=gammaw/grav
	raa=gammaa/grav
	uwt=xn*sat
	rmix=(1.d0-xn)*rs+uwt*rw+(xn-uwt)*raa
!
!	water phase parameters
!
	cww=xn*g1-xn*sat*bmkw
	cwg=-xn*g1
!
!	air phase parameters
!
	cgg=xn*g1-xn*(1.d0-sat)*bmka
	cgw=cwg
!
!	dimensionless variables
!
	gammam=rmix*grav
	Drw= rw/rmix
	Draa=raa/rmix
	Dgammaw=gammaw/gammam
	Dgammaa=gammaa/gammam
	Dxkw=xkw/DMIN1(xkw,xka)
	Dxka=xka/DMIN1(xkw,xka)
	Dxlambda=xlambda/(xlambda+2.D0*xmu)
	Dxmu=xmu/(xlambda+2.D0*xmu)
!
	return
	end
!
!
!*******************************************************************************************
!
!*******************************************************************************************

	subroutine STATEVT(sy,suc,a,b,c,e0,xm,sigb)
!
!	this subroutine calculates the void ratio ev, bulk modulus B and
!	suction modulus Dsuc^-1 from the new void ratio state surface
!
	implicit double precision (a-h,o-z)
	common /umate/ d(3,3),dfs(3),xn,xkw,cww,cwg,xka,cgg,cgw,ev,sat,g1,rs,rw,raa,rmix,xm1,&
				   bt,et,xlambda,xmu,henry
	common /e/ gammaw,gammas,gammaa,grav,atmp
!
!	if (suc.le.0.d0) suc=0.d0
	btinit=c
	sucn=suc/atmp
	syn=sy/atmp
	sigbn=sigb/atmp
	asig=a*sigbn
	bsuc=b*sucn
	if (bsuc.ge.asig) goto 100
	dnexp=c*(1.d0-xm)/atmp
	dnom1=a-b*sucn/sigbn
	dnom2=dnom1*syn+b*sucn
!
	if (xm.eq.1.) then
		ev=e0
		bt=(c/dnom1)*(dnom2**xm)
		xm1=(b*(1.-syn/sigbn))/(c*DNom2**xm)
		goto 500
	endif
!
	if (sucn.ne.0.d0) goto 200
	if (syn.eq.0.d0) then
		ev=e0
	else
		ev=(1.d0+e0)/(exp(((a*syn)**(1.-xm))/dnexp))-1.
	endif
	bt=c
	xm1=0.d0
	goto 500
!
200	continue
	ev=(1.d0+e0)/(exp ((dnom2)**(1.-xm)/dnexp))-1.d0
	bt=(c/dnom1)*(dnom2**xm)
	xm1=-b*(1.d0-syn/sigbn)/(c*DNom2**xm)
	goto 500
100	write (*,1) a,b,sigb,suc
	stop
500	continue
!
!	formats
!
1	FORMAT (5x,'state surface of void ratio is not correct'//,2x,4e10.3)
	return
	end subroutine
!
!
!*******************************************************************************************
!
	subroutine STATEST(sy,suc,b,sati)
!
!	this subroutine calculates the saturation and gradients
!
	implicit double precision (a-h,o-z)
	common /umate/ d(3,3),dfs(3),xn,xkw,cww,cwg,xka,cgg,cgw,ev,sat,g1,rs,rw,raa,rmix,xm1,&
				   bt,et,xlambda,xmu,henry
	common /e/ gammaw,gammas,gammaa,grav,atmp
!
!	if (suc.lt.0.d0) suc=0.d0
	cc=exp(b*suc)
!	cc=1.d0+b*suc
	sat=cc
	g1=b*sat
!
	if (cc.gt.1D0) sat=sati
!
	return
	end
!
!
!*******************************************************************************************
!
!*******************************************************************************************
!
	subroutine perwt(a,al,sru,pm)
!
!	this subroutine calculates the water permeability
!
	implicit double precision (a-h,o-z)
	common /umate/ d(3,3),dfs(3),xn,xkw,cww,cwg,xka,cgg,cgw,ev,sat,g1,rs,rw,raa,rmix,xm1,&
				   bt,et,xlambda,xmu,henry
	common /e/ gammaw,gammas,gammaa,grav,atmp
!
	if (sat.lt.0.d0) sat=0.d0
	if (sat.gt.1.d0) sat=1.d0
	xn=ev/(1.d0+ev)
	aa=sat-sru
	if (aa.gt.0.d0) goto 100
	xkw=0.d0
	if (xkw.eq.0.d0) xkw=10.d0**(-10.)
	goto 200
100	aa=aa/(1.d0-sru)
	aa=a*aa**3.5d0
	bb=10.d0**(al*ev)
	xkw=aa*bb
!	if (xkw.gt.pm) xkw=pm
	xkw=xkw/gammaw
!
200	continue
	return
	end
!
!
!*******************************************************************************************
!
!*******************************************************************************************
!
	subroutine PERAT(b,bet,xmua,sru,pm)
!
!	this subroutine calculates the air permeability
!
	implicit double precision (a-h,o-z)
	common /umate/ d(3,3),dfs(3),xn,xkw,cww,cwg,xka,cgg,cgw,ev,sat,g1,rs,rw,raa,rmix,xm1,&
				   bt,et,xlambda,xmu,henry
	common /e/ gammaw,gammas,gammaa,grav,atmp
!
	if (sat.lt.0.d0) sat=0.d0
	if (sat.gt.1.d0) sat=1.d0
	aa=1.d0-sat
	if (aa.gt.0.0001) goto 100
	xka=0.d0
	if (xka.eq.0.d0) xka=10.d0**(-15.)
	goto 200
100	aa=ev*aa
	aa=aa**bet
	xka=b*gammaa*aa/xmua
!	if (xka.gt.pm) xka=pm
	xka=xka/gammaa
!
200	continue
	return
	end
!
!*******************************************************************************************


!****************************************************************************************
!
!							           DMATU subroutine
!
!	This subroutine calculates the D matrix for unsat. soil and aslo suction matrix
!	This subroutine is called in: GHMATD1
!	This subroutine calls: -
!
!	variables used are:
!				eload	: (El) elastic modulus in absence of suction Et=El+Es
!				evkb	: (K0) bulk modulud in absence of suction (initial value of B)
!				ten		:
!				emin	: minimum value of Young modulus
!				sati	: initial saturation degree
!				bmkw	: Cw=d(row)/row/d(pw), compressibility of water
!				perwa	: aw(m.s^-1) used in the formula of water permeability
!				perwal	: alpha_w used in the formula of water permeability
!				sru		: used in the formula of water permeability
!				perwm	: kw_max, maximum authorized water permeability
!				bmka	: Ca=d(roa)/roa/d(pa), compressibility of air
!				perab	: ba used in the formula of air permeability
!				peraal	: alpfa_a used in the formula of air permeability
!				xmua	: viscosity of air
!				peram	: ka_max, maximum authorized air permeability
!				evae	: ae used in the formula of the void ratio state surface
!				evbe	: be used in the formula of the void ratio state surface
!				evce	: evkb(adim.)*P_atm
!				eve0	: e0 initial void ratio, used in the formula of the void ratio
!						  state surface
!				evsigb	: maximum traction resistance
!				srbetaw	: beta_w sed in the formula of the saturation degree state surface
!				sucm	: ms used in the formula of Es=ms(Pa-Pw)
!				henry	: Henry coefficient
!
!	INPUT	:
!
!	OUTPUT	:
!
!****************************************************************************************
!
!
!							            EXTINeq subroutine
!
!	This subroutine calculates the ehd1 & egd1 matrices in ordinary elements
!	This subroutine is called in: GHMATD3
!	This subroutine calls:		  EXTINeq1
!
!	variables used are:
!
!		nn			: it is equal to 3*kbem; where 3 is related to the number of nodes in
!					  each element then nn = [6, 9, or 12]
!		ehstd1		: H matrix in an ordinary element in static case
!		egstd1		: G matrix in an ordinary element in static case
!		ehd1		: H matrix in an ordinary element in dynamic case
!		egd1		: G matrix in an ordinary element in dynamic case
!		sehd1		: H matrix in each ordinary sub-elm in dynamic case
!		segd1		: G matrix in each ordinary sub-elm in dynamic case
!		sehstd1		: H matrix in each ordinary sub-elm in static case
!		segstd1		: G matrix in each ordinary sub-elm in static case
!		eItn		: ??????????
!		seItn		: ??????????
!
!		e11			: local coordinate of the first node in element [-1]
!		e22			: local coordinate of the third node in element [1]
!		ajab		: coefficient in tangential vector to 1D element (Vxsi) in x-direction
!					  Vxsi(x)=(ajab) xsi + bjab
!		bjab		: coefficient in tangential vector to 1D element (Vxsi) in x-direction
!					  Vxsi(x)=(ajab) xsi + bjab
!		cjab		: coefficient in tangential vector to 1D element (Vxsi) in y-direction
!					  Vxsi(y)=(cjab) xsi + djab
!		djab		: coefficient in tangential vector to 1D element (Vxsi) in y-direction
!					  Vxsi(y)=(cjab) xsi + djab
!		nsub		: number of subelements
!		deltae		: length of sub regions in each elemnt
!		itm			:
!		itime		: time increment; [1]: first time increment,
!									  [>1]: higher time increment
!		idyn		: dynamic code or type of loading [0,1,2]; [0]: static,
!															   [1]: quasi-static and,
!															   [2]: dynamic problem,
!
!
!	INPUT	: itm,itime,dtime,idyn,nodo,nsub,nptg,kbem,eItn
!
!	OUTPUT	: ehd1,egd1,ehstd1,egstd1
!
!****************************************************************************************
!
	subroutine extineq(itm,itime,dtime,idyn,nodo,nsub,nptg,ehd1,egd1,ehstd1,egstd1,&
					   eItn)
!
	implicit double precision (a-h,o-z)
!
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /mate/zlamda,zmuy,znuu,znu,zalpha,zrho,zrhof,zm,zk,zr,zn
	common /umate/ d(3,3),dfs(3),xn,xkw,cww,cwg,xka,cgg,cgw,ev,sat,g1,rs,rw,raa,rmix,xm1,&
				   bt,et,xlambda,xmu,henry
	common /jab/ajab,bjab,cjab,djab,xelm(3),yelm(3),xp,yp
!
	dimension ehd1(kbem,3*kbem),egd1(kbem,3*kbem),ehstd1(kbem,3*kbem),egstd1(kbem,3*kbem),&
			  sehstd1(kbem,3*kbem),segstd1(kbem,3*kbem),sehd1(kbem,3*kbem),segd1(kbem,3*kbem)
	dimension sItn(kbem,3*kbem),eItn(kbem,3*kbem)
!
!
!	------------------- initializing -------------------
!
	nn=3*kbem	! 3 nodes in the boundary elm.*[Ux,Uy,Pw,Pa] (we define it![2,3 or 4])
	do i=1,kbem	! load is applied to each D3 node. it can be inside of elemen studied or
				! outside (nodo). then in each node load can be in kbem direction
				! [Ux,Uy,Pw,Pa]
		do j=1,nn
			ehstd1(i,j)=0.d0
			egstd1(i,j)=0.d0
			ehd1(i,j)=0.d0
			egd1(i,j)=0.d0
			eItn(i,j)=0.d0
		enddo
	enddo
!
!	-------------- geometrical parameters --------------
!
	e11=-1.d0
	e22=1.d0
	ajab=xelm(3)-2.d0*xelm(2)+xelm(1)
	bjab=(xelm(3)-xelm(1))/2.d0
	cjab=yelm(3)-2.d0*yelm(2)+yelm(1)
	djab=(yelm(3)-yelm(1))/2.d0
	deltae=(e22-e11)/nsub
	ee1=e11
!
!	-- (ehstd1 & egstd1) and (ehd1 & egd1) in each element --
!
	do k=1,nsub
		ee2=ee1+deltae
!
!	-- (sehstd1 & segstd1) and (sehd1 & segd1) in each sub-region --
!
		call extineq1(itm,itime,dtime,nodo,nptg,ee1,ee2,sehd1,segd1,sehstd1,segstd1,&
					  sItn)
!
!	-- assembling (ehstd1 & egstd1) and (ehd1 & egd1) in each element --
!
		do 50 i=1,kbem
			do 40 j=1,nn
				if (itime.gt.1) goto 35
				ehstd1(i,j)=ehstd1(i,j)+sehstd1(i,j)
				egstd1(i,j)=egstd1(i,j)+segstd1(i,j)
				if (idyn.eq.0) goto 40
35				ehd1(i,j)=ehd1(i,j)+sehd1(i,j)
				egd1(i,j)=egd1(i,j)+segd1(i,j)
!				if(itime.eq.2) then
!				eItn(i,j)=eItn(i,j)+sItn(i,j)
!				endif
40			continue
50		continue
		ee1=ee2
	enddo
!
200	return
	end
!****************************************************************************************
!
!							            EXTINeq1 subroutine
!
!	This subroutine computes the sehstd1, segstd1, sehd1 & segd1 matrices that relate a
!	node (xp,yp) with a boundary sub-element using gauss quadrature. the sub-region g & h
!	matrices are made in this subroutine using Gauss quadrature WITHOUT CONSIDERING THE
!	SINGULARITY EFFECT (ordinary element) because the diagonal composants of gd3 & hd3 are
!	obtained in GHMATD3 by rigid body motion method.
!
!	This subroutine is called in: EXTINeq
!	This subroutine calls:		  SAICOR,ELSTFSOL,SATPOELSTFSOL,UNSATPOELSTFSOL,ELDYNFSOL,
!								  SATCOEFFICIENTSLAP,SATDYNFSOLGLAP,SATDYNFSOLTLAP,
!								  SATCOEFFICIENTSANAL,SATDYNFSOLGALT,SATDYNFSOLTALT,
!								  UNSATCOEFFICIENTSLAP,UNSATDYNFSOLGLAP,UNSATDYNFSOLTLAP
!
!	variables used are:
!
!		nptg	: number of Gauss points which is determined by the Stroud and Secrest method
!		gi		: the Gauss points which can be between 3 & 10
!		wi		: the Gauss weighting functions
!		nn		: it is equal to 3*kbem; where 3 is related to the number of nodes in each
!				  element then nn = [6, 9, or 12]
!		sehd1	: H matrix in each ordinary sub-elm in dynamic case
!		segd1	: G matrix in each ordinary sub-elm in dynamic case
!		sehstd1	: H matrix in each ordinary sub-elm in static case
!		segstd1	: G matrix in each ordinary sub-elm in static case
!		sItn	: ??????????
!
!		xja2	: second jacobian; it shows 1/nsub=Jacob_bar in the local coordiante system
!				  of the concerning sub-region. the length of the concerning subregions is
!				  equal to "deltae" obtained in EXTINeq where is (e22-e11)/nsub=2/nsub.
!				  xja2=(ee2-ee1)/2=deltae/2=1/nsub
!
!		e		: evaluating coef; tha Gauss points have to be defined in a local coordinate
!				  system of each sub-region. then one must change the global coordinate of
!				  the Gauss points which are between -1 & 1 to local coordinate system in each
!				  sub-region which is between ee1 & ee2. like a linear elm:
!				  X(xsi) = 0.5*(1-xsi)*xe(1) + 0.5 * (1+xsi)*xe(2)
!
!		f(1)	: shape function at first point
!		f(2)	: shape function at scond point
!		f(3)	: shape function at third point
!		xco		: x-component of integration point along the element
!		yco		: y-component of integration point along the element
!		dist(1)	: distance from collocation point to the gauss integration
!				  points on the x-plane; = xco-xp
!		dist(2)	: distance from collocation point to the gauss integration
!				  points on the y-plane; =yco-yp
!
!		xja1	: first jacobian
!				  sqrt[ (d y/d xsi)^2 + (d x/d xsi)^2 ] where
!					 Vxsi(x)=(d x/d xsi)=(xsi-0.5)xelm(1)-2 xsi xelm(2)+ (xsi+0.5)xelm(3)
!							=[xelm(3)-2.d0*xelm(2)+xelm(1)] xsi + (xelm(3)-xelm(1))/2.d0
!							=(ajab) xsi + bjab
!					 Vxsi(y)=(d y/d xsi)=(xsi-0.5)yelm(1)-2 xsi yelm(2)+ (xsi+0.5)yelm(3)
!							=[yelm(3)-2.d0*yelm(2)+yelm(1)] xsi + (yelm(3)-yelm(1))/2.d0
!							=(cjab) xsi + djab
!				  where
!					 (xsi-0.5)= d N1/d xsi
!					 -2 xsi   = d N2/d xsi
!					 (xsi-0.5)= d N3/d xsi
!
!		eta(1)	: component of the unit normal to the element in direction 1
!				  V3=Vxsi*Vz=[d y/d xsi, -d x/d xsi, 0] then
!				  eta(1)=Vxsi(y)=(e*cjab+djab)/xja1
!		eta(2)	: component of the unit normal to the element in direction 2
!				  V3=Vxsi*Vz=[d y/d xsi, -d x/d xsi, 0] then
!				  eta(2)= -Vxsi(x)=(e*ajab+bjab)/xja1
!
!		ra		: distance from collocation point to the gauss integration points on the
!				  boundary elements
!		rd(1)	: radius derivative dr/dx1; r,x
!		rd(2)	: radius derivative dr/dx2; r,y
!		rdn		: normal derivative of the distance r with respect to the normal direction
!				  n; dr/dn=r,i*ni=r,x*nx+r,y*ny
!		coef	: xja1*xja2*wi(i)
!
!		itime	: time increment; [1]: first time increment,
!								  [>1]: higher time increment
!		ibem	: BEM code, [0,1,2,3]; [0]: BE is not done,
!									   [1]: BEM_dry,
!									   [2]: BEM_saturate,
!									   [3]: BEM_unsaturate
!		idf		: identifier the colonel place in segstd1 & sehstd1


!		nodo	: code of singularity of element [0,1,2,3]
!									[0]: no singularity, loading point is outside of the
!										 boundary element studied
!									[1]: singularity in the first node of element
!									[2]: singularity in the second node of element
!									[3]: singularity in the third node of element
!		ee1		: first node in each sub-region
!		ee2		: second node in each sub-region
!
!
!
!	INPUT	: itm,itime,dtime,nodo,nptg,ee1,ee2,kbem,sItn
!
!	OUTPUT	: sehd1,segd1,sehstd1,segstd1
!
!
!****************************************************************************************
!
	subroutine extineq1(itm,itime,dtime,nodo,nptg,ee1,ee2,sehd1,segd1,sehstd1,segstd1,&
						sItn)
!
	implicit double precision (a-h,o-z)
!
	common /a/ nnp,ifem,ibem,idyn,ieaq,nr
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /mate/zlamda,zmuy,znuu,znu,zalpha,zrho,zrhof,zm,zk,zr,zn
	common /Dmate/cp,zt,zx,dzrho,dzrhof,dzlamda,dzmuy,dzk
	common /umate/ d(3,3),dfs(3),xn,xkw,cww,cwg,xka,cgg,cgw,ev,sat,g1,rs,rw,raa,rmix,xm1,&
				   bt,et,xlambda,xmu,henry
	common /jab/ajab,bjab,cjab,djab,xelm(3),yelm(3),xp,yp
	common /point/dist(2),ra,rd(2),rdn,eta(2)
	common /Dpoint/ ddist(2),dra
	common /Kcor/kcorf,mtime,rmtol,c1,c2
	common /ocqm/ EPScqm,taw,Icqm,Lcqm,NPcqm
!
	dimension sehd1(kbem,3*kbem),segd1(kbem,3*kbem),sehstd1(kbem,3*kbem),&
			  segstd1(kbem,3*kbem),idf(kbem),stg(kbem,kbem),stf(kbem,kbem),&
			  dg(kbem,kbem),df(kbem,kbem),f(3),gi(10),wi(10)
!
!	determine of factors integral
!
	dimension g10(10),w10(10),g9(9),w9(9),g8(8),w8(8),g7(7),w7(7),g6(6),w6(6),g5(5),w5(5),&
			  g4(4),w4(4),g3(3),w3(3)
	dimension sai(4,4),sItn(kbem,3*kbem) !?????
!
!
	save g10,w10,g9,w9,g8,w8,g7,w7,g6,w6,g5,w5,g4,w4,g3,w3
!
	data g10 / 0.973906528517172,-0.973906528517172,&
			   0.865063366688985,-0.865063366688985,&
			   0.679409568299024,-0.679409568299024,&
			   0.433395394129247,-0.433395394129247,&
			   0.148874338981631,-0.148874338981631 /
	data w10 / 0.066671344308688,0.066673344308688,&
			   0.149451349150581,0.149451349150581,&
			   0.219086362515982,0.219086362515982,&
			   0.269266719309996,0.269266719309996,&
			   0.295524224714753,0.295524224714753 /
!
	data g9 / 0.968160239507626,-0.968160239507626,&
			  0.836031107326636,-0.836031107326636,&
			  0.613371432700590,-0.613371432700590,&
			  0.324253423403809,-0.324253423403809,&
			  0.d0 /
	data w9 / 0.081274388361574,0.081274388361574,&
			  0.180648160694857,0.180648160694857,&
			  0.260610696402935,0.260610696402935,&
			  0.312347077040003,0.312347077040003,&
			  0.330239355001260 /
!
	data g8 / 0.960289856497536,-0.960289856497536,&
			  0.796666477413627,-0.796666477413627,&
			  0.525532409916329,-0.525532409916329,&
			  0.183434642495650,-0.183434642495650 /
	data w8 / 0.101228536290376,0.101228536290376,&
			  0.222381034453374,0.222381034453374,&
			  0.313706645877887,0.313706645877887,&
			  0.362683783378362,0.362683783378362 /
!
	data g7 / 0.949107912342759,-0.949107912342759,&
			  0.741531185599394,-0.741531185599394,&
			  0.405845151377397,-0.405845151377397,&
			  0.d0 /
	data w7 / 0.129484966168870,0.129484966168870,&
			  0.279705391489277,0.279705391489277,&
			  0.381830050505119,0.381830050505119,&
			  0.417959183673469 /
!
	data g6 / 0.932469514203152,-0.932469514203152,&
			  0.661209386466265,-0.661209386466265,&
			  0.238619186083197,-0.238619186083197 /
	data w6 / 0.171324492379170,0.171324492379170,&
			  0.360761573048139,0.360761573048139,&
			  0.467913934572691,0.467913934572691 /
!
	data g5 / 0.906179845938664,-0.906179845938664,&
			  0.538469310105683,-0.538469310105683,&
			  0.d0 /
	data w5 / 0.236926885056189,0.236926885056189,&
			  0.478628670499366,0.478628670499366,&
			  0.568888888888889 /
!
	data g4 / 0.861136311594953,-0.861136311594953,&
			  0.339981043584856,-0.339981043584856 /
	data w4 / 0.347854845137454,0.347854845137454,&
			  0.652145154862546,0.652145154862546 /
!
	data g3 / 0.774596669241483,-0.774596669241483,&
			  0.d0 /
	data w3 / 0.555555555555556,0.555555555555556,&
			  0.888888888888889 /
!
	data gi / 10*0.d0 /,wi / 10*0.d0 /
!
!
	goto (3,4,5,6,7,8,9,10) (nptg-2)
3	do i=1,nptg
		gi(i)=g3(i)
		wi(i)=w3(i)
	enddo
	goto 20
4	do i=1,nptg
		gi(i)=g4(i)
		wi(i)=w4(i)
	enddo
	goto 20
5	do i=1,nptg
		gi(i)=g5(i)
		wi(i)=w5(i)
	enddo
	goto 20
6	do i=1,nptg
		gi(i)=g6(i)
		wi(i)=w6(i)
	enddo
	goto 20
7	do i=1,nptg
		gi(i)=g7(i)
		wi(i)=w7(i)
	enddo
	goto 20
8	do i=1,nptg
		gi(i)=g8(i)
		wi(i)=w8(i)
	enddo
	goto 20
9	do i=1,nptg
		gi(i)=g9(i)
		wi(i)=w9(i)
	enddo
	goto 20
10	do i=1,nptg
		gi(i)=g10(i)
		wi(i)=w10(i)
	enddo
20	continue
!
!
!	------------------- initializing -------------------
!
	nn=3*kbem
	do i=1,kbem
		do j=1,nn
			segstd1(i,j)=0.d0
			sehstd1(i,j)=0.d0
			segd1(i,j)=0.d0
			sehd1(i,j)=0.d0
			sItn(i,j)=0.d0
		enddo
	enddo
	xja2=0.5d0*(ee2-ee1) !deltae/2=1/nsub
!
!
!	------------------- numerical integration -------------------
!	-------------------------------------------------------------
!
	do 600 i=1,nptg
!
!	---------- evaluating coefficients ----------
!
		e=0.5d0*(1-gi(i))*ee1+0.5d0*(1+gi(i))*ee2 ! coordinate of gauss point in subregion
!
!	---------- shape functions at integration points ----------
!
		f(1)=e*(e-1.d0)*0.5d0
		f(2)=1.d0-e**2
		f(3)=e*(e+1.d0)*0.5d0
!
!	------ geometrical properties at integration points -------
!
		xco=xelm(1)*f(1)+xelm(2)*f(2)+xelm(3)*f(3)
		yco=yelm(1)*f(1)+yelm(2)*f(2)+yelm(3)*f(3)
!
		dist(1)=xco-xp
		dist(2)=yco-yp
		xja1=dsqrt ((e*ajab+bjab)**2+(e*cjab+djab)**2)
		eta(1)=(e*cjab+djab)/xja1
		eta(2)=-(e*ajab+bjab)/xja1
		ra=dsqrt (dist(1)**2+dist(2)**2)
		rd(1)=dist(1)/ra	!e1=n1=dr/dx1=r1/r
		rd(2)=dist(2)/ra	!e2=n2=dr/dx2=r2/r
		rdn= rd(1)*eta(1)+rd(2)*eta(2)
		coef=xja1*xja2*wi(i)
!
		ddist(1)=dist(1)/zx
		ddist(2)=dist(2)/zx
		dra=ra/zx
!
!	---- computing the needed parameters for K-correction -----
!
!		if ((itime.eq.2).and.(kcorf.eq.1)) then
!			call saicor(sai,dtime,eta,rd,ra)
!		endif
!
!	--------- compute segstd1 & sehstd1 (static case) ---------
!
		if (itime.gt.1) goto 150
		stfundsol: SELECT CASE (ibem)
			CASE (1)
				call elstfsol(stg,stf)
			CASE (2)
				if (kfsol.eq.0.and.Icqm.eq.1) then
					call satpoelstfsolQST(stg,stf)
				else
					call satpoelstfsolQST(stg,stf) !satpoelstfsol(stg,stf)
				endif
			CASE (3)
				call unsatpoelstfsol1(stg,stf)
		end SELECT stfundsol
!
		do k=1,3
			do j1=1,kbem
				idf(j1)=kbem*(k-1)+j1
				do j2=1,kbem
					segstd1(j2,idf(j1))=segstd1(j2,idf(j1))+stg(j2,j1)*f(k)*coef
					sehstd1(j2,idf(j1))=sehstd1(j2,idf(j1))+stf(j2,j1)*f(k)*coef
				enddo
			enddo
		enddo
150		continue
!
!	----------- compute segd1 & sehd1 (dynamic case) ----------
!
		if (idyn.eq.0) goto 500
		dynfundsol: SELECT CASE (ibem)
			CASE (1)
				if (Icqm.eq.0) then
					call eldynfsol(itm,itime,dtime,dg,df)
				else
					call owcqm(dg,df)
				endif
			CASE (2)
				if (kfsol.eq.0.and.Icqm.eq.1) call owcqm(dg,df)
				if (kfsol.eq.1) then
					call satcoefficientsLap ! Aij,Bij,Cij & dAij,dBij
					call satdynfsolGLap(itm,itime,dtime,dg) ! Gij
					call satdynfsolTLap(itm,itime,dtime,df) ! Hij
				else if (kfsol.eq.2) then
					call satcoefficientsAlt
					call satdynfsolGAlt(itm,itime,dtime,dg)
					call satdynfsolTAlt(itm,itime,dtime,df)
				else if (kfsol.eq.3.and.Icqm.eq.1) then
					call satcoefficientsLap
					call owcqm(dg,df)
				endif
			CASE (3)
					call owcqm(dg,df)
			end SELECT dynfundsol
!
		do k=1,3
			do j1=1,kbem
				idf(j1)=kbem*(k-1)+j1
				do j2=1,kbem
					segd1(j2,idf(j1))=segd1(j2,idf(j1))+dg(j2,j1)*f(k)*coef
					if (itime.gt.1 .OR. k.ne.nodo) goto 175
					if (kbem.eq.2) then
						goto 170
					else if (kbem.eq.3) then
						if (j2.eq.3 .and. j1.ne.3) goto 175
						if (j1.eq.3 .and. j2.ne.3) goto 175
					else
						if (j2.eq.3 .and. j1.ne.3) goto 175
						if (j2.eq.4 .and. j1.ne.4) goto 175
						if (j2.ne.3 .and. j1.eq.3) goto 175
						if (j2.ne.4 .and. j1.eq.4) goto 175
					endif
170			        sehd1(j2,idf(j1))= sehd1(j2,idf(j1))+(df(j2,j1)-stf(j2,j1))*f(k)*coef
					goto 200
175			        sehd1(j2,idf(j1))= sehd1(j2,idf(j1))+df(j2,j1)*f(k)*coef
200					continue
!					if (itime.eq.2) then
!						sItn(j2,idf(j1))=sItn(j2,idf(j1))+sai(j2,j1)*f(k)*coef
!					endif
				enddo
			enddo
		enddo
500		continue
600	continue
!
	return
	end

!
	subroutine FS(RTSD,DRTSD,K0,K1)
!
	implicit double precision (a-h,o-z)
!
	common /mate/zlamda,zmuy,znuu,znu,zalpha,zrho,zrhof,zm,zk,zr,zn
!
	COMPLEX(8) W1,W,Z,GammaZ,SLAP,RTSD(4),DRTSD(2),K0(3),K1(3),&
			   UY(200),FUYANAL(200),FNOMEGAUY(200)
	dimension ABSUY(200),DUY(200),yUY(200),f(200),H1(200),R(3),rd(2),eta(2),dist(2),dlt(2,2)
!
	xco=1.d0
	yco=0.d0
	xp=0.d0
	yp=0.d0
!
	dist(1)=xco-xp
	dist(2)=yco-yp
	eta(1)=0
	eta(2)=1
	ra=dsqrt (dist(1)**2+dist(2)**2)
	rd(1)=dist(1)/ra
	rd(2)=dist(2)/ra
	rdn= rd(1)*eta(1)+rd(2)*eta(2)
	dlt(1,1)=1.d0
	dlt(2,2)=1.d0
	dlt(1,2)=0.d0
	dlt(2,1)=0.d0

!
	pi=DACOS (-1.d0)
	ntime=100
	dtime=0.0002
!
	epsilon=10.d0**(-10.d0)
	taw=epsilon**(1.d0/(2.d0*ntime))
	W1=CMPLX(0.d0,2.d0*pi/ntime)
	W=CDexp(W1)
!
10	do l=1,ntime
		Z=taw*W**(l-1)
		GammaZ=0.d0
		do i=1,2
			GammaZ=GammaZ+((1.d0-Z)**i)/i
		enddo
		SLAP=GammaZ/dtime
		K=1
		J=1
		do i=1,3
			R(i)=(2.d0*rd(k)*rd(j)-dlt(k,j))*DRTSD(i)*K1(i)/ra+rd(k)*rd(j)*RTSD(i)*K0(i)
		enddo
		UY(L)=( ( (RTSD(4)-RTSD(2))/(RTSD(1)-RTSD(2)) )*R(1) -&
			( (RTSD(4)-RTSD(1))/(RTSD(1)-RTSD(2)) )*R(2) +&
				dlt(k,j)*RTSD(3)*K0(3)-R(3) )/(2.D0*pi*zrho*SLAP**2.d0)
		absUY(L)=Cdabs (UY(L))
	enddo
!
	call  DFTF (ntime,UY,FUYANAL)

	do l=1,200
		FNOMEGAUY(l)=taw**(-l+1)*FUYANAL(l)/ntime
		DUY(l)=REAL(FNOMEGAUY(l))
	enddo
!
	do n=1,200
		do k=1,n
			if ((k-1)*dtime.ge.0.d0) H1(k)=1.d0
			f(k)=H1(k)
			yUY(n)=yUY(n)+DUY(n-k+1)*f(k)
		enddo
	enddo
!
	return
	end

!****************************************************************************************
!
!							           FSTRESS subroutine
!
!	This subroutine calculates stresses of the elements.
!
!	this subroutine is called in: CONTROL
!
!	variables used are:
!
!		sig8d	:
!		sig8c	:
!		sig8u	:
!		sig8di	:
!		sig8ci	:
!		sig8ui	:
!		sigmad	:
!		sigmac	:
!		sigmau	:
!		sig8do	:
!		sig8co	:
!		sig8uo	:
!		di		:
!		disp	:
!		dt		:
!		dpad	:
!		dpac	:
!		dpau	:
!		alfad	:
!		alfac	:
!		alfau	:
!		iysd	:
!		iysc	:
!		iysu	:
!		dt		:
!		dti		:
!		xe(8)	:
!		ye(8)	:
!		npw		: numbers of nodes with given (temporary) p-water values
!		nnpw	: node number of nodes with temporary imposed pw  (vector of dim = npw)
!		vnpw	: values of imposed p-water  (vector of dim = npw)

!		npa		: numbers of nodes with given (temporary) p-air values
!		nnpa	: node number of nodes with temporary imposed pa  (vector of dim = npa)
!		vnpa	: values of imposed p-air  (vector of dim = npa)
!
!
!	INPUT	:
!
!	OUTPUT	:
!
!
!****************************************************************************************
!
	subroutine fstress(ie8d,ie8c,ie8u,ie3d1,id,x,sm8d,sm8c,sm8u,sig3u,&
					   rbed1,gihd1,gid1,sig8d,sig8c,sig8u,sig8di,sig8ci,sig8ui,sigmad,sigmac,sigmau,sig8do,&
					   sig8co,sig8uo,di,disp,dt,dpad,dpac,dpau,alfad,alfac,alfau,iysd,iysc,&
					   iysu)
!
 	implicit double precision (a-h,o-z)
!
	common /a/ nnp,ifem,ibem,idyn,ieaq,nr
	common /a1/ ne8d,ne8c,ne8u
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /a3/ ne8d1,ne8c1,ne8u1,ne3u1,nload1,nbcx1,nbcy1,nbcw1,nbca1,m8d1,m8c1,m8u1,&
				nbed11,nbed21,nbed31,nbeint1,ne8dd,ne8cc,ne8uu
	common /c/ nload,ntime,itmax,itime,iter,iload,time,dtime,dtimei,dtimem,xinc
   	common /h/ mdof,mdofn,ls,ls1,isolv
	common /g/ init,iprint
	common /chrg1/ nstep,istep
	common /chrg2/ npw,npa,nnpw(200),nnpa(200),vnpw(200),vnpa(200)
	common /ini1/ vnpwi(400),vnpai(400),nnpwi(400),nnpai(400)
	common /ini2/ npwi,npai,ndispi,nveli,nacci
	common /unbem/ nelmu,ie3u(3,500,5)
!
	dimension ie8d(9,ne8d1),ie8c(9,ne8c1),ie8u(9,ne8u1),ie3d1(npbmax,nbed11),xe(8),ye(8),&
			  id(4,nnp),&
			  x(2,nnp),sm8d(60,m8d1),sm8c(60,m8c1),sm8u(60,m8u1),sig8d(9,5,ne8d1),&
			  sig8c(9,6,ne8c1),sig8u(9,7,ne8u1),sig8di(9,5,ne8d1),sig8ci(9,6,ne8c1),&
			  sig8ui(9,7,ne8u1),sigmad(9,5,ne8d1),sigmac(9,6,ne8c1),sigmau(9,7,ne8u1),&
			  sig8do(9,5,ne8d1),sig8co(9,6,ne8c1),sig8uo(9,7,ne8u1),sige(9,6),sigei(9,6),&
			  sigma(9,6),sigo(9,6),dpad(3,3,ne8dd),dpac(3,3,ne8cc),dpau(3,3,ne8uu),&
			  alfad(4,10,ne8dd),alfac(4,10,ne8cc),alfau(4,10,ne8uu),iysd(2,ne8dd),&
			  iysc(2,ne8cc),iysu(2,ne8uu),di(mdofn),disp(mdofn),dt(mdofn),pw(3),dpw(3),&
			  pa(3),dpa(3),du(6),u(6),xd(npbmax),xelm(3),yelm(3),sig3u(5,500,nbed11),sigeu(5),&
			  ud1(n12gh,nbed11),td1(n12gh,nbed11),rbed1(ntime,n12gh,nbed11),&
			  gihd1(n12gh,n12gh,nbed11),gid1(n12gh,n12gh,nbed11),trac(6)
!
!
!	-------------- calculate drained element stresses --------------
!
	if (ne8d.eq.0) goto 10
	do m=1,ne8d
		ield=m
        mt=ie8d(9,m)
        do i=1,8
			j=ie8d(i,m)
			xe(i)=x(1,j)
			ye(i)=x(2,j)
			k1=id(1,j)
			k2=id(2,j)
			du(2*i-1)=dt(k1)-disp(k1)
			du(2*i)  =dt(k2)-disp(k2)
			u(2*i-1) =dt(k1)
			u(2*i)   =dt(k2)
		enddo
		do j=1,9
			do i=1,5
				sig8di(j,i,m)=sig8d(j,i,m)
				sige(j,i)=sig8d(j,i,m)
				sigma(j,i)=sigmad(j,i,m)
				sigo(j,i)=sig8do(j,i,m)
			enddo
		enddo
		call str8nd(ield,mt,xe,ye,du,u,sige,sigei,sigma,sigo,sm8d,sm8c,dpad,alfad,iysd)
		do j=1,9
			do i=1,5
				sig8d(j,i,m)=sige(j,i)
			enddo
		enddo
	enddo
!
!
!	-------------- calculate consolidation element stresses --------------
!
10	if (ne8c.eq.0) goto 15
	do m=1,ne8c
		ielc=m
        mt=ie8c(9,m)
        do i=1,8
			j=ie8c(i,m)
			xe(i)=x(1,j)
			ye(i)=x(2,j)
			k1=id(1,j)
			k2=id(2,j)
			k3=id(3,j)
			du(2*i-1)=dt(k1)-disp(k1)
			du(2*i)  =dt(k2)-disp(k2)
			u(2*i-1)=dt(k1)
			u(2*i)  =dt(k2)
			pw(i)=dt(k3)
			dpw(i)=di(k3)
		enddo
		do j=1,9
			do i=1,6
				sig8ci(j,i,m)=sig8c(j,i,m)
				sige(j,i)=sig8c(j,i,m)
				sigma(j,i)=sigmac(j,i,m)
				sigo(j,i)=sig8co(j,i,m)
			enddo
		enddo
		call str8nc(ielc,mt,xe,ye,du,u,pw,dpw,sige,sigei,sigma,sigo,sm8d,sm8c,dpac,alfac,iysc)
		do j=1,9
			do i=1,6
				sig8c(j,i,m)=sige(j,i)
			enddo
		enddo
	enddo
!
!
15	if (ibem.ne.3) goto 20
	if (nbed1.eq.0) goto 16
	do k1=1,nbed1
		do k=1,nnpbed1(k1)
			j=ie3d1(k,k1)
			do i=1,kbem
				ud1(kbem*k-kbem+i,k1)=dt(id(i,j))
			enddo
		enddo
		do i=1,nnpbed1(k1)*kbem
			td1(i,k1)=0.d0
			do j=1,nnpbed1(k1)*kbem
				td1(i,k1) = td1(i,k1)+gihd1(i,j,k1)*ud1(j,k1)-&
								    gid1(i,j,k1)*rbed1(itime,j,k1)
			enddo
		enddo
!
		n=nnpbed1(k1)
		mt=ie3d1(n+1,k1)
		do m=1,nelmu
			do i=1,3
				j=ie3u(i,m,k1)
				trac(2*i-1)=td1(j*kbem-3,k1)
				trac(2*i)=td1(j*kbem-2,k1)
				xelm(i)=x(1,j)
				yelm(i)=x(2,j)
				nu1=id(1,j)
				nu2=id(2,j)
				nw =id(3,j)
				na =id(4,j)
				u(2*i-1)=dt(nu1)
				u(2*i)  =dt(nu2)
				pw(i)   =dt(nw)
				pa(i)   =dt(na)
				dpw(i)  =di(nw)
				dpa(i)  =di(na)
			enddo
!			call str3nu(m,mt,trac,xelm,yelm,u,pw,pa,dpw,dpa,sigeu,sm8u)
!			do i=1,5
!				sig3u(i,m,k1)=sigeu(i)
!			enddo
		enddo
	enddo
!
16	if (nbed2.eq.0) goto 17
17	if (nbed3.eq.0) goto 20


20	if (npw.eq.0) goto 25
	do i=1,npw
		j=nnpw(i)
		i3=id(3,j)
		if (i3.le.mdofn) dt(i3)=vnpw(i)
	enddo
!
25	if (npa.eq.0) goto 30
	do i=1,npa
		j=nnpa(i)
		i4=id(4,j)
		if (i4.le.mdofn) dt(i4)=vnpa(i)
	enddo
!
30	if (init.eq.0) goto 50
	do i=1,mdof
		dt(i)=0.d0
	enddo
	if (npwi.eq.0) goto 50
	do i=1,npwi
		j=nnpwi(i)
		i3=id(3,j)
		if (i3.le.mdof) dt(i3)=vnpwi(i)
		if (i3.le.mdof) disp(i3)=vnpwi(i)
	enddo
	if (npai.eq.0) goto 50
	do i=1,npai
		j=nnpai(i)
		i4=id(4,j)
		if (i4.le.mdof) dt(i4)=vnpwi(i)
		if (i4.le.mdof) disp(i4)=vnpwi(i)
	enddo
!
50	continue
!
	return
	end

!****************************************************************************************
!
!							            ELSTFSOL subroutine
!
!	This subroutine computes the elasto(EL)-static(ST) fundamental(F) solution(SOL)
!
!	This subroutine is called in these subroutines: EXTINeq1
!
!	variables used are:
!
!		stg		: fundamental solutions for the displacements
!		   stg(1,1): solution for displacement in x-dir due to a unit load in x-dir
!		   stg(1,2): solution for displacement in y-dir due to a unit load in x-dir; G21=Gyx
!		   stg(2,2): solution for displacement in y-dir due to a unit load in y-dir
!		   stg(2,1): solution for displacement in x-dir due to a unit load in y-dir

!		stf		: fundamental solutions for the tractions
!
!	INPUT	: zlamda,zmuy,ra,rd(1),rd(2),rdn,eta(1),eta(2)
!	OUTPUT	: stg,stf
!
!****************************************************************************************
!
	subroutine elstfsol(stg,stf)
!
	implicit double precision (a-h,o-z)
!
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /point/dist(2),ra,rd(2),rdn,eta(2)
	common /mate/zlamda,zmuy,znuu,znu,zalpha,zrho,zrhof,zm,zk,zr,zn
!
	dimension stg(kbem,kbem),stf(kbem,kbem)
!
	pi=DACOS(-1.d0)
	znu=zlamda/(2.d0*(zlamda+zmuy))
	de=4*pi*(1.d0-znu)
	stg(1,1)=((3.d0-4.d0*znu)*LOG(1.d0/ra)+rd(1)**2.d0)/(2.d0*de*zmuy)	! Uxx
	stg(1,2)=rd(1)*rd(2)/(2.d0*de*zmuy)									! Uxy
	stg(2,1)=stg(1,2)													! Uyx
	stg(2,2)=((3.d0-4.d0*znu)*LOG(1.d0/ra)+rd(2)**2.d0)/(2.d0*de*zmuy)	! Uyy
	stf(1,1)=-rdn*((1.d0-2.d0*znu)+2.d0*rd(1)**2.d0)/(ra*de)			! Txx
	stf(1,2)=-(rdn*2.d0*rd(1)*rd(2)+(1.d0-2.d0*znu)*(eta(1)*rd(2)-eta(2)*rd(1)))/(ra*de)
																		! Txy
	stf(2,1)=-(rdn*2.d0*rd(1)*rd(2)+(1.d0-2.d0*znu)*(eta(2)*rd(1)-eta(1)*rd(2)))/(ra*de)
																		! Tyx
	stf(2,2)=-rdn*((1.d0-2.d0*znu)+2.d0*rd(2)**2.d0)/(ra*de)			! Tyy
!
	return
	end
!
!
!****************************************************************************************
!
!****************************************************************************************
!
!							       SATPOELSTFSOL subroutine
!
!	This subroutine calculates SATurated-POro-ELasto-STatique Fundamental SOLutions.
!	This subroutine is called in: EXTINeq1
!	This subroutine calls: -
!
!	Gij-> i:response; j=loading direction
!						G=	[ G11   G21   G31 ]
!							[ G12   G22   G32 ]
!							[ G13   G23   G33 ]
!
!	variables used are:
!
!		dist(1)	: distance from collocation point to the gauss integration
!				  points on the x-plane; = xco-xp
!		dist(2)	: distance from collocation point to the gauss integration
!				  points on the y-plane; =yco-yp
!		ra		: distance from collocation point to the gauss integration points on the
!				  boundary elements
!		rd(1)	: radius derivative dr/dx1
!		rd(2)	: radius derivative dr/dx2
!		rdn		: radius derivative dr/dn
!		eta(1)	: component of the unit normal to the element in direction 1
!				  V3=Vxsi*Vz=[d y/d xsi, -d x/d xsi, 0] then
!				  eta(1)=Vxsi(y)=(e*cjab+djab)/xja1
!		eta(2)	: component of the unit normal to the element in direction 2
!				  V3=Vxsi*Vz=[d y/d xsi, -d x/d xsi, 0] then
!				  eta(2)= -Vxsi(x)=(e*ajab+bjab)/xja1
!		zlamda	: Lam�'s coefficient which is equal to K-2MU/3
!		zmuy	: Lam�'s coefficient which is equal to 3EK/(9K-E)
!		zalpha	: alpha, in saturated formulation
!		zrho	: density of body
!		zrhof	: density of fluid (water) in a saturated volume
!		zm		: 1/Q, in saturated formulation
!		zk		: water permeability in z-direction
!
!	INPUT		: stg,stf
!	OUTPUT		: solPEstG,solPEstT
!
!****************************************************************************************
!
	subroutine satpoelstfsol(solPEstG,solPEstT)
!
	implicit double precision (a-h,o-z)
!
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /point/ dist(2),ra,rd(2),rdn,eta(2)
	common /Dpoint/ ddist(2),dra
	common /mate/zlamda,zmuy,znuu,znu,zalpha,zrho,zrhof,zm,zk,zr,zn
	common /Dmate/cp,zt,zx,dzrho,dzrhof,dzlamda,dzmuy,dzk
	common /ocqm/ EPScqm,taw,Icqm,Lcqm,NPcqm
!
	dimension solPEstG(kbem,kbem),solPEstT(kbem,kbem),dist1(2)
!
!
	pi=DACOS(-1.d0)
	znu=zlamda/(2*(zlamda+zmuy))
	ZE=zmuy*(3.d0*zlamda+2.d0*zmuy)/(zlamda+zmuy)
	de=4*pi*(1.d0-znu)
!
!
	if (Icqm.eq.1.and.kfsol.eq.3) then
		solPEstG(1,1)=((3.d0-4.d0*znu)*dlog(1.d0/ra)+rd(1)**2.d0)/(2.d0*de*zmuy)
		solPEstG(1,2)=rd(1)*rd(2)/(2.d0*de*zmuy)
		solPEstG(2,1)=solPEstG(1,2)
		solPEstG(2,2)=((3.d0-4.d0*znu)*dlog(1.d0/ra)+rd(2)**2.d0)/(2.d0*de*zmuy)
		solPEstG(1,3)=0.d0
		solPEstG(2,3)=0.d0
		solPEstG(3,1)=0.d0
		solPEstG(3,2)=0.d0
		solPEstG(3,3)=dlog(ra)/(2.d0*pi*zk)
!
		solPEstT(1,1)=-rdn*( (1.d0-2.d0*znu)+2.d0*rd(1)**2.d0 )/(ra*de)
		solPEstT(1,2)=( -2.d0*rd(1)*rd(2)*rdn+(1.d0-2.d0*znu)*(rd(1)*eta(2)-rd(2)*eta(1)) )/&
					   (ra*de)
		solPEstT(2,1)=( -2.d0*rd(1)*rd(2)*rdn+(1.d0-2.d0*znu)*(rd(2)*eta(1)-rd(1)*eta(2)) )/&
					   (ra*de)
		solPEstT(2,2)=-rdn*( (1.d0-2.d0*znu)+2.d0*rd(2)**2.d0 )/(ra*de)
		solPEstT(1,3)=eta(1)*zalpha*(1.d0-2*znu)/(8*pi*zmuy*(1.d0-znu))*dlog(ra)
		solPEstT(2,3)=eta(2)*zalpha*(1.d0-2*znu)/(8*pi*zmuy*(1.d0-znu))*dlog(ra)
		solPEstT(3,1)=0.d0
		solPEstT(3,2)=0.d0
		solPEstT(3,3)=-rdn/(2*pi*ra)
		goto 60
	endif
!
!
	do i=1,3	! direction of load [Ux, Uy or Pw] in source point
		do j=1,3	! response direction of field point
			if (i.eq.3 .and. j.ne.3) goto 20
			if (j.eq.3 .and. i.ne.3) goto 30
			if (i.eq.3 .and. j.eq.3) goto 40
!
!
!	---------------- stg(1,2) & stf(1,2) ----------------
!
			x=rd(i)*rd(j)
			if (i.eq.j) x=x-(3.d0-4.d0*znu)*dlog(ra)
			solPEstG(i,j)=x/(8.d0*pi*zmuy*(1.d0-znu))
			x=2*rd(i)*rd(j)
			if (i.eq.j) x=x+(1.d0-2*znu)
			x=-x*rdn + (1.d0-2*znu)*(rd(i)*eta(j)-rd(j)*eta(i))
			solPEstT(i,j)=x/(4*pi*ra*(1.d0-znu))
			goto 50
!
!
!	---------------- stg(3,j) & stf(3,j) ----------------
!
20			solPEstG(i,j)=0.D0
			solPEstT(i,j)=0.D0
			goto 50
!
!
!	---------------- stg(i,3) & stf(i,3) ----------------
!
30			solPEstG(i,j)=zalpha/(8.d0*pi*zk*(zlamda+2*zmuy))*dist(i)*(1.d0-2.d0*dlog(ra))
			solPEstT(i,j)=zalpha*zmuy/(4*pi*zk*(zlamda+2*zmuy))*&
						              ((1.d0+2*dlog(ra))*eta(i)-2*rd(i)*rdn)
			goto 50
!
!
!	---------------- stg(3,3) & stf(3,3) ----------------
!
40			solPEstG(i,j)=-dlog(ra)/(2*pi*zk)
			solPEstT(i,j)=rdn/(2*pi*ra)
!
50			continue
		enddo
	enddo
!
	G12=solPEstG(1,2)
	T12=solPEstT(1,2)
	G21=solPEstG(2,1)
	T21=solPEstT(2,1)
	G13=solPEstG(1,3)
	T13=solPEstT(1,3)
	G23=solPEstG(2,3)
	T23=solPEstT(2,3)
	G31=solPEstG(3,1)
	T31=solPEstT(3,1)
	G32=solPEstG(3,2)
	T32=solPEstT(3,2)
!
	solPEstG(1,2)=G21
	solPEstG(2,1)=G12
	solPEstT(1,2)=T21
	solPEstT(2,1)=T12
!
	solPEstG(1,3)=G31
	solPEstG(2,3)=G32
	solPEstG(3,1)=G13
	solPEstG(3,2)=G23
!
	solPEstT(1,3)=T31
	solPEstT(2,3)=T32
	solPEstT(3,1)=T13
	solPEstT(3,2)=T23
!
60	continue
	return
	end
!
!
!****************************************************************************************
!
!****************************************************************************************
!
!							       SATPOELSTFSOLQST subroutine
!
!	This subroutine calculates SATurated-POro-ELasto-STatique Fundamental SOLutions for
!	Quasi-Static problem by CONVOLUTION QUADRATURE METHOD.
!	This subroutine is called in: EXTINeq1
!	This subroutine calls: -
!
!	Gij-> i:response; j=loading direction
!						G=	[ G11   G21   G31 ]
!							[ G12   G22   G32 ]
!							[ G13   G23   G33 ]
!
!	INPUT		: stg,stf
!	OUTPUT		: solPEstG,solPEstT
!
!****************************************************************************************
!
	subroutine satpoelstfsolQST(solPEstG,solPEstT)
!
	implicit double precision (a-h,o-z)
!
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /point/ dist(2),ra,rd(2),rdn,eta(2)
	common /mate/zlamda,zmuy,znuu,znu,zalpha,zrho,zrhof,zm,zk,zr,zn
!
	dimension solPEstG(kbem,kbem),solPEstT(kbem,kbem)
!
!
	pi=DACOS(-1.d0)
!
	do i=1,3	! direction of load [Ux, Uy or Pw] in source point
		do j=1,3	! response direction of field point
			if (i.eq.3 .and. j.ne.3) goto 20
			if (j.eq.3 .and. i.ne.3) goto 30
			if (i.eq.3 .and. j.eq.3) goto 40
!
!
!	---------------- stg(1,2)=G21 & stf(1,2)=F21 ----------------
!
			x=rd(i)*rd(j)
			if (i.eq.j) x=x-(3.d0-4.d0*znu)*dlog(ra)
			solPEstG(i,j)=x/(8.d0*pi*zmuy*(1.d0-znu))
			x=2*rd(i)*rd(j)
			if (i.eq.j) x=x+(1.d0-2*znu)
			x=-x*rdn + (1.d0-2*znu)*(rd(i)*eta(j)-rd(j)*eta(i))
			solPEstT(i,j)=x/(4*pi*ra*(1.d0-znu))
			goto 50
!
!
!	---------------- stg(3,j)=Gi3 & stf(3,j)=Ti3 ----------------
!
20			solPEstT(i,j)=0.d0
			solPEstG(i,j)=0.d0
			goto 50
!
!
!	---------------- stg(i,3)=G3j & stf(i,3)=T3j ----------------
!
30			solPEstG(i,j)=0.D0
			solPEstT(i,j)=zalpha*(1.D0-2.D0*znu)/( 8*pi*zmuy*(1.D0-znu) )*dlog(ra)*eta(i)
			goto 50
!
!
!	---------------- stg(3,3)=G33 & stf(3,3)=T33 ----------------
!
40			solPEstG(i,j)=dlog(ra)/(2*pi*zk)
			solPEstT(i,j)=-rdn/(2*pi*ra)
!
50			continue
		enddo
	enddo
!
	return
	end
!
!
!****************************************************************************************
!
!****************************************************************************************

	subroutine unsatpoelstfsol(unstG,unstF)
!
	implicit double precision (a-h,o-z)
!
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /e/ gammaw,gammas,gammaa,grav,atmp
	common /umate/ d(3,3),dfs(3),xn,xkw,cww,cwg,xka,cgg,cgw,ev,sat,g1,rs,rw,raa,rmix,xm1,&
				   bt,et,xlambda,xmu,henry
	common /point/ dist(2),ra,rd(2),rdn,eta(2)
!
	dimension unstG(kbem,kbem),unstF(kbem,kbem),dlt(2,2)
!
!
	do i=1,kbem
		do j=1,kbem
			unstG(i,j)=0.d0
			unstF(i,j)=0.d0
		enddo
	enddo
!
	dlt(1,1)=1.d0
	dlt(2,2)=1.d0
	dlt(1,2)=0.d0
	dlt(2,1)=0.d0
!
	pi=DACOS (-1.d0)
	d1=3.d0*bt/(9.d0*bt-et)
	d3=3.d0*bt-et
	zlamda=d1*d3
	zmu=d1*et
	znu=zlamda/(2*(zlamda+zmu))
	fs=dfs(1)
	zkw=xkw*rw/gammaw
	zka=xka*raa/gammaa
!
!
	do i=1,4		! direction of load [Ux,Uy,Pw,Pa] in source point
		do j=1,4	! response direction of field point
			if (i.lt.3 .and. j.eq.3) goto 30 !if (i.lt.3 .and. j.ge.3) goto 20
			if (i.lt.3 .and. j.eq.4) goto 60 !if (i.eq.3 .and. j.lt.3) goto 30
			if (i.eq.3 .and. j.eq.3) goto 40
			if (i.eq.3 .and. j.eq.4) goto 50
			if (i.ge.3 .and. j.lt.3) goto 20 !if (i.eq.4 .and. j.lt.3) goto 20
			if (i.eq.4 .and. j.eq.3) goto 55
			if (i.eq.4 .and. j.eq.4) goto 70
!
!	---------------- Gij & Fij ----------------
!
			x=rd(i)*rd(j)
			if (i.eq.j) x=x-(3.d0-4.d0*znu)*dlog(ra)
			unstG(i,j)=x/(8.d0*pi*zmu*(1.d0-znu))
			x=2*rd(i)*rd(j)
			if (i.eq.j) x=x+(1.d0-2*znu)
			x=-x*rdn + (1.d0-2*znu)*(rd(i)*eta(j)-rd(j)*eta(i))
			unstF(i,j)=x/(4*pi*ra*(1.d0-znu))
			goto 80
!
!	------------ G3j-G4j & F3j-F4j ------------
!
20			unstG(i,j)=0.d0
			unstF(i,j)=0.d0
			goto 80
!
!	---------------- Gi3 & Fi3 ----------------
!
30			if (xkw.eq.0.d0) then
				unstG(i,j)=0.d0
				unstF(i,j)=0.d0
				goto 80
			endif
			unstG(i,j)=0.d0 !-fs*dist(i)*(1.D0-2.D0*dlog(ra))/&
					   !(8.d0*pi*(zlamda+2*zmu)*xkw)
!			unstG(i,j)=-(fs*gammaw*xka-henry*gammaa*(1.d0-fs)*xkw)*dist(i)*&
!					   (1.D0-2.D0*dlog(ra))/&
!					   (8.d0*pi*(zlamda+2*zmu)*xkw*xka*rw)
!
			unstF(i,j)=fs/(8.d0*pi*(zlamda+2*zmu))*&
					   ((1.d0-2.d0*dlog(ra))*eta(i)-2*rd(i)*rdn) !((1.d0+2.d0*dlog(ra))*eta(i)+2*rd(i)*rdn)
			goto 80
!
!	---------------- G33 & T33 ----------------
!
40			if (xkw.eq.0.d0) then
				unstG(i,j)=0.d0
				unstF(i,j)=0.d0
				goto 80
			endif
			unstG(i,j)=-dlog(ra)/(2.D0*pi*xkw) !+
!			unstG(i,j)=-dlog(ra)/(2.D0*pi*zkw)
!
			unstF(i,j)=rdn/(2*pi*ra) !-
			goto 80
!
!	------------ G34 & T34 ------------
!
50			unstG(i,j)=0.d0
			unstF(i,j)=0.d0
			goto 80
!
!	------------ G43 & T43 ------------
!
55			unstG(i,j)=0.D0 !henry*gammaa*dlog(ra)/(2.D0*pi*rw*xka) !XK73*PH11
			unstF(i,j)=0.d0
			goto 80
!
!	---------------- Gi4 & Fi4 ----------------
!
60			if (xka.eq.0.d0) then
				unstG(i,j)=0.d0
				unstF(i,j)=0.d0
				goto 80
			endif
			unstG(i,j)=0.d0 !-(1.D0-fs)*dist(i)*(1.D0-2.D0*dlog(ra))/&
					   !(8.d0*pi*(zlamda+2*zmu)*xka)
!			unstG(i,j)=-(1.D0-fs)*dist(i)*(1.D0-2.D0*dlog(ra))/&
!					   (8.d0*pi*(zlamda+2*zmu)*zka)  !0.d0
!
			unstF(i,j)=(1.d0-fs)/(8.d0*pi*(zlamda+2*zmu))*&
					   ((1.d0-2.d0*dlog(ra))*eta(i)-2*rd(i)*rdn) !((1.d0+2.d0*dlog(ra))*eta(i)+2*rd(i)*rdn)
			goto 80
!
!	---------------- G44 & T44 ----------------
!
70			if (xka.eq.0.d0) then
				unstG(i,j)=0.d0
				unstF(i,j)=0.d0
				goto 80
			endif
			unstG(i,j)=-dlog(ra)/(2.D0*pi*xka) !+
!			unstG(i,j)=-dlog(ra)/(2.D0*pi*zka)
!
			unstF(i,j)=rdn/(2*pi*ra) !-
			goto 80
!
80			continue
		enddo
	enddo
!
	return
	end
!
!
!*******************************************************************************************
!
!****************************************************************************************
!
!****************************************************************************************

	subroutine unsatpoelstfsol1(unstG,unstF)
!
	implicit double precision (a-h,o-z)
!
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /e/ gammaw,gammas,gammaa,grav,atmp
	common /umate/ d(3,3),dfs(3),xn,xkw,cww,cwg,xka,cgg,cgw,ev,sat,g1,rs,rw,raa,rmix,xm1,&
				   bt,et,xlambda,xmu,henry
	common /point/ dist(2),ra,rd(2),rdn,eta(2)
!
	dimension unstG(kbem,kbem),unstF(kbem,kbem),dlt(2,2)
!
!
	do i=1,kbem
		do j=1,kbem
			unstG(i,j)=0.d0
			unstF(i,j)=0.d0
		enddo
	enddo
!
	dlt(1,1)=1.d0
	dlt(2,2)=1.d0
	dlt(1,2)=0.d0
	dlt(2,1)=0.d0
!
	pi=DACOS (-1.d0)
	xE=xmu*(3.d0*xlambda+2.d0*xmu)/(xlambda+xmu)
	xnu=xlambda/(2.d0*(xlambda+xmu))
	fs=dfs(1)
!
!
	do i=1,4		! direction of load [Ux,Uy,Pw,Pa] in source point
		do j=1,4	! response direction of field point
			if (i.lt.3 .and. j.eq.3) goto 30 !if (i.lt.3 .and. j.ge.3) goto 20
			if (i.lt.3 .and. j.eq.4) goto 60 !if (i.eq.3 .and. j.lt.3) goto 30
			if (i.eq.3 .and. j.eq.3) goto 40
			if (i.eq.3 .and. j.eq.4) goto 50
			if (i.ge.3 .and. j.lt.3) goto 20 !if (i.eq.4 .and. j.lt.3) goto 20
			if (i.eq.4 .and. j.eq.3) goto 55
			if (i.eq.4 .and. j.eq.4) goto 70
!
!	---------------- Gij & Fij ----------------
!
			x=rd(i)*rd(j)
			if (i.eq.j) x=x-(3.d0-4.d0*xnu)*dlog(ra)
			unstG(i,j)=x/(8.d0*pi*xmu*(1.d0-xnu))
			x=2*rd(i)*rd(j)
			if (i.eq.j) x=x+(1.d0-2*xnu)
			x=-x*rdn + (1.d0-2*xnu)*(rd(i)*eta(j)-rd(j)*eta(i))
			unstF(i,j)=x/(4*pi*ra*(1.d0-xnu))
			goto 80
!
!	------------ G3j-G4j & F3j-F4j ------------
!
20			unstG(i,j)=0.d0
			unstF(i,j)=0.d0
			goto 80
!
!	---------------- Gi3 & Fi3 ----------------
!
30			if (xkw.eq.0.d0) then
				unstG(i,j)=0.d0
				unstF(i,j)=0.d0
				goto 80
			endif
			unstG(i,j)=0.D0 !fs*dist(i)*(1.D0-2.D0*dlog(ra))/&
					   !(8.d0*pi*(xlambda+2*xmu)*rw*xkw)
			unstF(i,j)=0.D0 !-fs/(8.d0*pi*(xlambda+2*xmu))*&
					   !((1.d0-2.d0*dlog(ra))*eta(i)-2*rd(i)*rdn)
			goto 80
!
!	---------------- G33 & T33 ----------------
!
40			if (xkw.eq.0.d0) then
				unstG(i,j)=0.d0
				unstF(i,j)=0.d0
				goto 80
			endif
			unstG(i,j)=-dlog(ra)/(2.D0*pi*rw*xkw)
			unstF(i,j)=rdn/(2*pi*ra)
			goto 80
!
!	------------ G34 & T34 ------------
!
50			unstG(i,j)=0.d0
			unstF(i,j)=0.d0
			goto 80
!
!	------------ G43 & T43 ------------
!
55			unstG(i,j)=0.D0
			unstF(i,j)=0.d0
			goto 80
!
!	---------------- Gi4 & Fi4 ----------------
!
60			if (xka.eq.0.d0) then
				unstG(i,j)=0.d0
				unstF(i,j)=0.d0
				goto 80
			endif
			unstG(i,j)=0.D0 !(1.D0-fs)*dist(i)*(1.D0-2.D0*dlog(ra))/&
					   !(8.d0*pi*(xlambda+2*xmu)*raa*xka)
			unstF(i,j)=0.D0 !-(1.d0-fs)/(8.d0*pi*(xlambda+2*xmu))*&
					   !((1.d0-2.d0*dlog(ra))*eta(i)-2*rd(i)*rdn)
			goto 80
!
!	---------------- G44 & T44 ----------------
!
70			if (xka.eq.0.d0) then
				unstG(i,j)=0.d0
				unstF(i,j)=0.d0
				goto 80
			endif
			unstG(i,j)=-dlog(ra)/(2.D0*pi*raa*xka)
			unstF(i,j)=rdn/(2*pi*ra)
			goto 80
!
80			continue
		enddo
	enddo
!
	return
	end
!
!
!*******************************************************************************************

!****************************************************************************************
!
!								   SATDYNFSOLGALT subroutine
!
!	This subroutine calculates the SATurated-DYNamique Fundamental SOLutions of G AnaLyTic
!	ot dynamic solutions of displacement in the case of incompressible composents.
!	This subroutine is called in: EXTINeq1
!	This subroutine calls: -
!
!	Gij-> i:response; j=loading direction
!					constant interpolation: G=[ G11       G21       int(G31) ]
!											  [ G12       G22       int(G32) ]
!											  [ dt(G13)   dt(G23)   G33      ]
!					linear   interpolation: G=[ int(G11)  int(G21)  int2(G31)]
!											  [ int(G12)  int(G22)  int2(G32)]
!											  [ G13       G23       int(G33) ]
!					mix      interpolation: G=[ G11       G21       int(G31) ]
!											  [ G12       G22       int(G32) ]
!											  [ G13       G23       int(G33) ]
!	variables used are:
!
!		itime	: time increment from 1 to "ntime"
!							[1] : first time increment,
!							[>1]: higher time increment
!		itm		:
!		dtime	: time step; in seismic loading for all "1:ntime", it is constant
!		solPEdyG: dg, displacement fundamental solution
!		isol	: first index in Gij which is the response of medium, displc. or pressures
!		jsol	: second index in Gij which is the direction of applied force
!		msol	:
!		nh		:
!		nd		: [-2,-1,0,1]
!							[-2]: two times of integration
!							[-1]: one time of integration
!							[0]	: normal
!							[1]	: one time of derivation
!
!****************************************************************************************
!
	subroutine satdynfsolGAlt(itm,itime,dtime,solPEdyG)
!
	implicit double precision (a-h,o-z)
!
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /Dmate/cp,zt,zx,dzrho,dzrhof,dzlamda,dzmuy,dzk
!
	dimension solPEdyG(kbem,kbem)
!
	common /indicesol/ isol,jsol,msol,nh,nd
	common /stability/ anew1,anew2,wil,intb
!
!
	t3=(itime-itm)*dtime+wil*dtime
	if (itime.gt.1) then
		t2=(itime-itm-1)*dtime+wil*dtime
        t1=(itime-itm-2)*dtime+wil*dtime
	else
		t2=0.d0
		t1=0.d0
	endif
!
!	t3=t3/zt
!	t2=t2/zt
!	t1=t1/zt
!	dtime1=dtime
!	dtime=dtime/zt
!
	do isol=1,3
		do jsol=1,3
			SolPEdyG(isol,jsol)=0.d0
			if (jsol.ne.3) goto 10
			if (jsol.eq.3) goto 20
!
!
!	----------------- Gij & G3j -----------------
!
10			if (intb.eq.1) then		! constant temporal interpollation
				nd=0
				if (isol.eq.3) nd=-1
				solPEdyG(isol,jsol)= soldyGA(t2,t3)
			endif
!
			if (intb.eq.2) then		! linear temporal interpollation
				nd=-1
				if (isol.eq.3) nd=-2
				if (itime.gt.1) then
					solPEdyG(isol,jsol)= (soldyGA(t2,t3)-soldyGA(t1,t2))/dtime
				else
					solPEdyG(isol,jsol)= soldyGA(t2,t3)/dtime
				endif
			endif
!
			if (intb.eq.3) then		! mixte temporal interpollation
				nd=0
				if (isol.eq.3) nd=-1
				if (itime.eq.1) then
					solPEdyG(isol,jsol)= soldyGA(t1,t3)
				else
					solPEdyG(isol,jsol)= soldyGA(t1,t3)/2.d0
				endif
			endif
!
			goto 30
!
!
!	----------------- Gi3 & G33 -----------------
!
20			if (intb.eq.1) then		! constant temporal interpollation
				nd=1
				if (isol.eq.3) nd=0
				solPEdyG(isol,jsol)= soldyGA(t2,t3)
			endif
!
			if (intb.gt.1) then		! linear & mixture temporal interpollation
				nd=0
				if (isol.eq.3) nd=-1
				if (itime.gt.1) then
					solPEdyG(isol,jsol)= (soldyGA(t2,t3)-soldyGA(t1,t2))/dtime
				else
					solPEdyG(isol,jsol)= soldyGA(t2,t3)/dtime
				endif
			endif
!
30			continue
!
		enddo
	enddo
!
	G11=solPEdyG(1,1)
	G12=solPEdyG(1,2)
	G13=solPEdyG(1,3)
	G21=solPEdyG(2,1)
	G22=solPEdyG(2,2)
	G23=solPEdyG(2,3)
	G31=solPEdyG(3,1)
	G32=solPEdyG(3,2)
	G33=solPEdyG(3,3)
!
	solPEdyG(1,2)=G21
	solPEdyG(1,3)=G31
	solPEdyG(2,1)=G12
	solPEdyG(2,3)=G32
	solPEdyG(3,1)=G13
	solPEdyG(3,2)=G23
!
!	dtime=dtime1
!
	return
	end
!
!
!****************************************************************************************
!
!								   SATDYNFSOLTALT subroutine
!
!	This subroutine calculates the SATurated-DYNamique Fundamental SOLutions of T AnaLyTic
!	ot dynamic solutions of traction in the case of incompressible composents.

!****************************************************************************************
!
	subroutine satdynfsolTAlt(itm,itime,dtime,solPEdyT)
	implicit double precision (a-h,o-z)
!
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /point/ dist(2),ra,rd(2),rdn,eta(2)
	common /Dpoint/ ddist(2),dra
	common /mate/zlamda,zmuy,znuu,znu,zalpha,zrho,zrhof,zm,zk,zr,zn
	common /Dmate/cp,zt,zx,dzrho,dzrhof,dzlamda,dzmuy,dzk
	common /indicesol/ isol,jsol,msol,nh,nd
	common /stability/ anew1,anew2,wil,intb
!
	dimension dpsolG1(kbem,kbem,2),dpsolG2(kbem,kbem,2),dpsolG3(kbem,kbem,2),&
			  solPEdyT(kbem,kbem),dist1(2)
!
!
!	zlamda1=zlamda
!	zmuy1=zmuy
!	zrho1=zrho
!!	zrhof1=zrhof
!	zk1=zk
!	dist1(1)=dist(1)
!	dist1(2)=dist(2)
!	ra1=ra
!
!	zlamda=dzlamda
!	zmuy=dzmuy
!	zrho=dzrho
!	zrhof=dzrhof
!	zk=dzk
!	dist(1)=ddist(1)
!	dist(2)=ddist(2)
!	ra=dra
!
	t3=(itime-itm)*dtime+wil*dtime
	if (itime.gt.1) then
		t2=(itime-itm-1)*dtime+wil*dtime
        t1=(itime-itm-2)*dtime+wil*dtime
	else
		t2=0.d0
		t1=0.d0
	endif
!
!	t3=t3/zt
!	t2=t2/zt
!	t1=t1/zt
!	dtime1=dtime
!	dtime=dtime/zt
!
	if (intb.eq.1 .OR.intb.eq.2) then
		k=0
		call dpsoldyG(t2,t3,dpsolG1,k)
		if (intb.eq.2 .and. itime.gt.1) call dpsoldyG(t1,t2,dpsolG2,k)
	else
		k=1
		call dpsoldyG(t2,t3,dpsolG1,k)
		if (itime.gt.1) call dpsoldyG(t1,t2,dpsolG2,k)
		k=2
		call dpsoldyG(t1,t3,dpsolG3,k)
	endif
!
	do i=1,3
		do j=1,3
			temp1=0.d0
			temp2=0.d0
			isol=i
			jsol=j
			if (isol.ne.3 .and. jsol.eq.3) goto 20
			if (isol.eq.3 .and. jsol.ne.3) goto 30
			if (isol.eq.3 .and. jsol.eq.3) goto 40
!
!	Tij
			nd=0
			if (intb.gt.1) nd=-1
			do m=1,2
				temp1=temp1 + (dpsolG1(i,j,m)+dpsolG1(m,j,i))*eta(m)
				if (intb.gt.1 .and. itime.gt.1)&
					temp2=temp2 + (dpsolG2(i,j,m)+dpsolG2(m,j,i))*eta(m)
			enddo
			isol=3
			pres1=soldyGA(t2,t3)
			if (intb.ne.1 .and. itime.gt.1) pres2=soldyGA(t1,t2)
			isol=i
			dila1=dpsolG1(1,j,1)+dpsolG1(2,j,2)
			if (intb.ne.1 .and. itime.gt.1)&
				dila2=dpsolG2(1,j,1)+dpsolG2(2,j,2)
			temp1=temp1*zmuy+(dila1*zlamda-pres1*zalpha)*eta(i)
			if (intb.ne.1 .and. itime.gt.1)&
				temp2=temp2*zmuy+(dila2*zlamda-pres2*zalpha)*eta(i)
			if (intb.eq.1) then
				solPEdyT(i,j)=temp1
			else
				solPEdyT(i,j)=(temp1-temp2)/dtime
			endif
			goto 50
!
!	Ti3
!
20			if (intb.eq.1) then
				nd=1
				do m=1,2
					temp1=temp1 + (dpsolG1(i,j,m)+dpsolG1(m,j,i))*eta(m)
				enddo
				isol=3
				pres1=soldyGA(t2,t3)
				isol=i
				dila1=dpsolG1(1,j,1)+dpsolG1(2,j,2)
				temp1=temp1*zmuy+(dila1*zlamda-pres1*zalpha)*eta(i)
				solPEdyT(i,j)=temp1
			endif
			if (intb.eq.3) then
				nd=1
				do m=1,2
					temp1=temp1 + (dpsolG3(i,j,m)+dpsolG3(m,j,i))*eta(m)
				enddo
				isol=3
				pres1=soldyGA(t1,t3)
				isol=i
				dila1=dpsolG1(1,j,1)+dpsolG1(2,j,2)
				temp1=temp1*zmuy+(dila1*zlamda-pres1*zalpha)*eta(i)
				if (itime.eq.1) then
					solPEdyT(i,j)=temp1
				else
					solPEdyT(i,j)=temp1/2.d0
				endif
			endif
			if (intb.eq.2) then
				nd=0
				do m=1,2
					temp1=temp1 + (dpsolG1(i,j,m)+dpsolG1(m,j,i))*eta(m)
					if (itime.gt.1)&
						temp2=temp2 + (dpsolG2(i,j,m)+dpsolG2(m,j,i))*eta(m)
				enddo
				isol=3
				pres1=soldyGA(t2,t3)
				if (itime.gt.1) pres2=soldyGA(t1,t2)
				isol=i
				dila1=dpsolG1(1,j,1)+dpsolG1(2,j,2)
				if (itime.gt.1)&
					dila2=dpsolG2(1,j,1)+dpsolG2(2,j,2)
				temp1=temp1*zmuy+(dila1*zlamda-pres1*zalpha)*eta(i)
				if (itime.gt.1)&
					temp2=temp2*zmuy+(dila2*zlamda-pres2*zalpha)*eta(i)
				solPEdyT(i,j)=(temp1-temp2)/dtime
			endif
			goto 50
!
!	T3j
!
30			nd=-1
			if (intb.gt.2) nd=-2 !****
			do m=1,2
				isol=m
				nd=nd+2
				dep1=soldyGA(t2,t3)
				if (intb.gt.1 .and. itime.gt.1) dep2=soldyGA(t1,t2)
				isol=i
				nd=nd-2
				temp1=temp1+(dpsolG1(i,j,m)+zrhof*dep1)*eta(m)
				if (intb.gt.1 .and. itime.gt.1)&
					temp2=temp2+(dpsolG2(i,j,m)+zrhof*dep2)*eta(m)
			enddo
			if (intb.eq.1) then
				solPEdyT(i,j)=-temp1*zk
			else
				solPEdyT(i,j)=-(temp1-temp2)*zk/dtime
			endif
			goto 50
!
!	T33
!
40			if (intb.eq.1) then
				nd=0
				do m=1,2
					isol=m
					nd=nd+2
					dep1=soldyGA(t2,t3)
					isol=i
					nd=nd-2
					temp1=temp1+(dpsolG1(i,j,m)+zrhof*dep1)*eta(m)
				enddo
				solPEdyT(i,j)=-temp1*zk
			endif
			if (intb.eq.3) then
				nd=0
				do m=1,2
					isol=m
					nd=nd+2
					dep1=soldyGA(t1,t3)
					isol=i
					nd=nd-2
					temp1=temp1+(dpsolG3(i,j,m)+zrhof*dep1)*eta(m)
				enddo
				if (itime.eq.1) then
					solPEdyT(i,j)=-temp1*zk
				else
					solPEdyT(i,j)=-temp1*zk/2.d0
				endif
			endif
			if (intb.eq.2) then
				do m=1,2
					isol=m
					nd=nd+2
					dep1=soldyGA(t2,t3)
					if (itime.gt.1) dep2=soldyGA(t1,t2)
					isol=i
					nd=nd-2
					temp1=temp1+(dpsolG1(i,j,m)+zrhof*dep1)*eta(m)
					if (itime.gt.1)&
						temp2=temp2+(dpsolG2(i,j,m)+zrhof*dep2)*eta(m)
				enddo
				solPEdyT(i,j)=-(temp1-temp2)*zk/dtime
			endif
!
50			continue
		enddo
	enddo
!
	T11=solPEdyT(1,1)
	T12=solPEdyT(1,2)
	T13=solPEdyT(1,3)
	T21=solPEdyT(2,1)
	T22=solPEdyT(2,2)
	T23=solPEdyT(2,3)
	T31=solPEdyT(3,1)
	T32=solPEdyT(3,2)
	T33=solPEdyT(3,3)
!
	solPEdyT(1,2)=T21
	solPEdyT(1,3)=T31
	solPEdyT(2,1)=T12
	solPEdyT(2,3)=T32
	solPEdyT(3,1)=T13
	solPEdyT(3,2)=T23
!
!	zlamda=zlamda1
!	zmuy=zmuy1
!	zrho=zrho1
!	zrhof=zrhof1
!	zk=zk1
!	dist(1)=dist1(1)
!	dist(2)=dist1(2)
!	ra=ra1
!
!	dtime=dtime1
!
	return
	end
!
!
!****************************************************************************************
!
!					          SATCOEFFICIENTSALT subroutine
!
!	This subroutine calculates the coefficients of analytical SATurated-ELasto-DYNamique
!	Fundamental SOLutions & their spatial derivatives in the case of incompressible
!	composents.
!	This subroutine is called in: EXTINeq1
!	This subroutine calls: -
!
!	variables used are:
!
!
!	INPUT	:
!
!	OUTPUT	:
!
!
!****************************************************************************************
!
	subroutine satcoefficientsAlt
!
	implicit double precision (a-h,o-z)
!
	dimension a(2,2),b(2,2),c(2,2),d(2,2),dpa(2,2,2),dpb(2,2,2),dpd(2,2,2),dist1(2)
!
	common /point/ dist(2),ra,rd(2),rdn,eta(2)
	common /Dpoint/ ddist(2),dra
	common /mate/zlamda,zmuy,znuu,znu,zalpha,zrho,zrhof,zm,zk,zr,zn
	common /Dmate/cp,zt,zx,dzrho,dzrhof,dzlamda,dzmuy,dzk
	common /coefA/ ee(9,2,2),dpee(8,2,2,2),dd(6),ff(4),gg(6),v(2),bt
!
!
!	------------------- constant -------------------
!
	pi=DACOS(-1.d0)
!
!	zlamda1=zlamda
!	zmuy1=zmuy
!	zrho1=zrho
!	zrhof1=zrhof
!	zk1=zk
!	dist1(1)=dist(1)
!	dist1(2)=dist(2)
!	ra1=ra
!
!	zlamda=dzlamda
!	zmuy=dzmuy
!	zrho=dzrho
!	zrhof=dzrhof
!	zk=dzk
!	dist(1)=ddist(1)
!	dist(2)=ddist(2)
!	ra=dra
!
	v(1)=dsqrt((zlamda+2*zmuy)/(zrho-zalpha*zrhof))
    v(2)=dsqrt(zmuy/zrho)
	bt=zalpha**2.d0/(2.d0*zk*(zrho-zalpha*zrhof))
!
!
!	----------------- Aij & dpAijm -----------------
!
	temp=2.d0*pi*ra**3.d0
	a(1,1)=(2.d0*dist(1)**2.d0 - ra**2.d0)/temp
	a(1,2)=2.d0*dist(1)*dist(2)/temp
	a(2,1)=a(1,2)
	a(2,2)=(2.d0*dist(2)**2.d0 - ra**2.d0)/temp
!
	temp=temp*ra
	dpa(1,1,1)=(4.d0*ra*dist(1)-6.d0*rd(1)*dist(1)**2.d0+ra**2.d0*rd(1))/temp
	dpa(1,2,1)=(2.d0*ra*dist(2)-6.d0*rd(1)*dist(1)*dist(2))/temp
	dpa(2,1,1)=dpa(1,2,1)
	dpa(2,2,1)=(-6.d0*rd(1)*dist(2)**2.d0+ra**2.d0*rd(1))/temp
	dpa(2,2,2)=(4.d0*ra*dist(2)-6.d0*rd(2)*dist(2)**2.d0+ra**2.d0*rd(2))/temp
	dpa(1,2,2)=(2.d0*ra*dist(1)-6.d0*rd(2)*dist(1)*dist(2))/temp
	dpa(2,1,2)=dpa(1,2,2)
	dpa(1,1,2)=(-6.d0*rd(2)*dist(1)**2.d0+ra**2.d0*rd(2))/temp
!
!
!	----------------- Bij & dpBijm -----------------
!
	temp=2.d0*pi*ra**2.d0
	b(1,1)=dist(1)**2.d0/temp
	b(2,2)=dist(2)**2.d0/temp
	b(1,2)=dist(1)*dist(2)/temp
	b(2,1)=b(1,2)
!
	temp=temp*ra
	dpb(1,1,1)=(2.d0*ra*dist(1)-2.d0*rd(1)*dist(1)**2.d0)/temp
	dpb(1,2,1)=(ra*dist(2)-2.d0*rd(1)*dist(1)*dist(2))/temp
	dpb(2,1,1)=dpb(1,2,1)
	dpb(2,2,1)=-2.d0*rd(1)*dist(2)**2.d0/temp
	dpb(2,2,2)=(2.d0*ra*dist(2)-2.d0*rd(2)*dist(2)**2.d0)/temp
	dpb(1,2,2)=(ra*dist(1)-2.d0*rd(2)*dist(1)*dist(2))/temp
	dpb(2,1,2)=dpb(1,2,2)
	dpb(1,1,2)=-2.d0*rd(2)*dist(1)**2.d0/temp
!
!
!	----------------- Cij (dpCijm=0) -----------------
!
	temp=1.d0/(2.d0*pi*zmuy)
	c(1,1)=temp
	c(2,2)=temp
	c(1,2)=0.d0
	c(2,1)=0.d0
!
!
!	----------------- Dij (dpDij) -----------------
!
	do i=1,2
		do j=1,2
			d(i,j)=a(i,j)/ra
		enddo
	enddo
	temp=pi*ra**5.d0
	dpd(1,1,1)=(2.d0*ra*dist(1)+ra**2.d0*rd(1)-4.d0*rd(1)*dist(1)**2.d0)/temp
	dpd(1,2,1)=(ra*dist(2)-4.d0*rd(1)*dist(1)*dist(2))/temp
	dpd(2,1,1)=dpd(1,2,1)
	dpd(2,2,1)=(ra**2.d0*rd(1)-4.d0*rd(1)*dist(2)**2.d0)/temp
	dpd(2,2,2)=(2.d0*ra*dist(2)+ra**2.d0*rd(2)-4.d0*rd(2)*dist(2)**2.d0)/temp
	dpd(1,2,2)=(ra*dist(1)-4.d0*rd(2)*dist(1)*dist(2))/temp
	dpd(2,1,2)=dpd(1,2,2)
	dpd(1,1,2)=(ra**2.d0*rd(2)-4.d0*rd(2)*dist(1)**2.d0)/temp
!
!
!	----------------- EEij (dpEEijm) -----------------
!
	do i=1,2
		do j=1,2
			temp=zlamda+2.d0*zmuy
			ee(1,i,j)=a(i,j)/temp
			ee(2,i,j)=b(i,j)/temp
			ee(3,i,j)=d(i,j)/zrho
!
			temp=zrho-zalpha*zrhof
			ee(4,i,j)=-d(i,j)/(temp*2.d0*bt)
			ee(5,i,j)=d(i,j)/(temp*4.d0*bt**2.d0)
			ee(6,i,j)=-a(i,j)/zmuy
			ee(7,i,j)=-b(i,j)/zmuy+c(i,j)
			ee(8,i,j)=ee(7,i,j)-ee(6,i,j)*ra/2.d0
			ee(9,i,j)=ee(5,i,j)/(2.d0*bt)
		enddo
	enddo
	do i=1,2
		do j=1,2
			do m=1,2
				temp=zlamda+2.d0*zmuy
				dpee(1,i,j,m)=dpa(i,j,m)/temp
				dpee(2,i,j,m)=dpb(i,j,m)/temp
				dpee(3,i,j,m)=dpd(i,j,m)/zrho
				temp=zrho-zalpha*zrhof
				dpee(4,i,j,m)=-dpd(i,j,m)/(temp*2.d0*bt)
				dpee(5,i,j,m)=dpd(i,j,m)/(temp*4.d0*bt**2.d0)
				dpee(6,i,j,m)=-dpa(i,j,m)/zmuy
				dpee(7,i,j,m)=-dpb(i,j,m)/zmuy
			enddo
		enddo
	enddo
!
!
!	----------------- DD -----------------
!
	temp=2.d0*pi*(zlamda+2.d0*zmuy)
	dd(1)=-zalpha/(temp*zk)
	dd(2)=zrhof/temp
	dd(3)=1.d0/(2.d0*pi*zalpha)
	temp=zrho-zalpha*zrhof
	dd(4)=-zrho/(2.d0*pi*zalpha*temp)
	dd(5)=dd(4)/(2.d0*bt)
	dd(6)=dd(5)/(2.d0*bt)
!
!
!	----------------- FF -----------------
!
	ff(1)=dd(1)
	ff(2)=dd(3)
	ff(3)=-zk*temp/(2.d0*pi*zalpha**3.d0)
	ff(4)=ff(2)*2.d0*bt
!
	gg(1)=1.d0/(2.d0*pi*zk)
	gg(2)=-zrho/(2.d0*pi*zk*temp)
	gg(3)=gg(1)+gg(2)
	gg(4)=-gg(2)*2.d0*bt
	gg(5)=gg(1)*2.d0*bt
	gg(6)=gg(2)/(2.d0*bt)
!
!	zlamda=zlamda1
!	zmuy=zmuy1
!	zrho=zrho1
!	zrhof=zrhof1
!	zk=zk1
!	dist(1)=dist1(1)
!	dist(2)=dist1(2)
!	ra=ra1
!
	return
	end
!
!
!	--------------- the h funtions & their spatial, temporal derivatives ---------------
!
	double precision function hklA(t)
!
	implicit double precision (a-h,o-z)
!
	common /point/ dist(2),ra,rd(2),rdn,eta(2)
	common /Dpoint/ ddist(2),dra
	common /coefA/ ee(9,2,2),dpee(8,2,2,2),dd(6),ff(4),gg(6),v(2),bt
	common /indicesol/ isol,jsol,msol,nh,nd
!
	dimension dist1(2)
!
!	dist1(1)=dist(1)
!	dist1(2)=dist(2)
!	ra1=ra
!
!	dist(1)=ddist(1)
!	dist(2)=ddist(2)
!	ra=dra
!
	hklA=0.d0
	goto (1,2,3,4,5,6) nh
!
!	----------------- h11 -----------------
!
1	if (t.le.(ra/v(1))) goto 100
	x=dsqrt(t**2-(ra/v(1))**2)
	hklA=Dexp(-bt*t)*DCOSH(bt*x)/x
	goto 100
!
!
!	----------------- h12 -----------------
!
2	if (t.le.(ra/v(1))) goto 100
	x=dsqrt(t**2-(ra/v(1))**2)
	hklA=v(1)**2*Dexp(-bt*t)*DSINH(bt*x)/(ra*bt)
	goto 100
!
!
!	----------------- dth12 -----------------
!
3	if (t.le.(ra/v(1))) goto 100
	x=dsqrt(t**2-(ra/v(1))**2)
	hklA=v(1)**2*Dexp(-bt*t)/ra*(DCOSH(bt*x)*t/x - DSINH(bt*x))
	goto 100
!
!
!	----------------- dph12 -----------------
!
4	if (t.le.(ra/v(1))) goto 100
	x=dsqrt(t**2-(ra/v(1))**2)
	hklA=-Dexp(-bt*t)*(DCOSH(bt*x)/x+v(1)**2*DSINH(bt*x)/(bt*ra**2))
	goto 100
!
!
!	----------------- h21 -----------------
!
5	if (t.le.(ra/v(2))) goto 100
	x=dsqrt(t**2-(ra/v(2))**2)
	hklA=1.d0/x
	goto 100
!
!
!	----------------- h22 -----------------
!
6	if (t.le.(ra/v(2))) goto 100
	x=dsqrt(t**2-(ra/v(2))**2)
	hklA=v(2)**2*x/ra
!
100 continue
!
!	dist(1)=dist1(1)
!	dist(2)=dist1(2)
!	ra=ra1
!
	return
	end
!
!
!****************************************************************************************
!
!				     fonction locale 1 pour int�gration num�rique
!
!****************************************************************************************
!
	 double precision function flocal1(tau,t)
!
	implicit double precision (a-h,o-z)
!
	common /point/ dist(2),ra,rd(2),rdn,eta(2)
	common /Dpoint/ ddist(2),dra
	common /coefA/ ee(9,2,2),dpee(8,2,2,2),dd(6),ff(4),gg(6),v(2),bt
	common /indicesol/ isol,jsol,msol,nh,nd
!
	dimension dist1(2)
!
!	dist1(1)=dist(1)
!	dist1(2)=dist(2)
!	ra1=ra
!
!	dist(1)=ddist(1)
!	dist(2)=ddist(2)
!	ra=dra
!
	flocal1=0.d0
    if (nh.ge.7) goto 10
!	hklA(tau)
	flocal1=hklA(tau)
    goto 200
10  goto (20,30,40,50,60,70,80,90,100,110,120) (nh-6)
!
!	pour int Gij
!
20	nh=2 ; temp=(t-tau)*hklA(tau)*ee(1,isol,jsol)
    nh=1 ; flocal1=temp + (t-tau)*hklA(tau)*ee(2,isol,jsol)
    nh=7 ; goto 200
!
!	pour Gij
!
30	nh=2 ; temp=hklA(tau)*ee(1,isol,jsol)
    nh=1 ; flocal1=temp + hklA(tau)*ee(2,isol,jsol)
    nh=8 ; goto 200
!
!	pour int int G3j
!
40	nh=2 ; flocal1=(dd(1)*(t-tau)+dd(2))*hklA(tau)
    nh=9 ; goto 200
!
!	pour int G33
!
50	nh=1 ; flocal1=(gg(1)*(t-tau)+gg(6)*(1.d0-dexp(-2*bt*(t-tau))))*hklA(tau)
	nh=10 ;goto 200
!
!	pour G33
!
60	nh=1 ; flocal1=(gg(1)+gg(2)*dexp(-2*bt*(t-tau)))*hklA(tau)
	nh=11 ; goto 200
!
!	pour dt G33
!
70	nh=1 ; flocal1=dexp(-2*bt*(t-tau))*hklA(tau)
    nh=12 ; goto 200
!
!	pour dp int Gij
!
80	nh=4 ; temp=ee(1,isol,jsol)*rd(msol)*(t-tau)*hklA(tau)
	nh=2 ; temp=temp+dpee(1,isol,jsol,msol)*(t-tau)*hklA(tau)
    nh=1 ; flocal1=temp+dpee(2,isol,jsol,msol)*(t-tau)*hklA(tau)
    nh=2 ; temp=ee(2,isol,jsol)*rd(msol)*2*bt*hklA(tau)/(v(1)**2)
    flocal1=flocal1-temp
    nh=13 ; goto 200
!
!	pour dp Gij
!
90	nh=4 ; temp=ee(1,isol,jsol)*rd(msol)*hklA(tau)
    nh=2 ; temp=temp+dpee(1,isol,jsol,msol)*hklA(tau)
    nh=1 ; flocal1=temp+dpee(2,isol,jsol,msol)*hklA(tau)
    nh=14 ; goto 200
!
!	pour dp int int G3j
!
100	if (jsol.eq.msol) rjm=(1.d0-rd(jsol)*rd(msol))/ra
    if (jsol.ne.msol) rjm=-rd(jsol)*rd(msol)/ra
    nh=4 ; temp=rd(jsol)*(dd(1)*(t-tau)+dd(2))*hklA(tau)*rd(msol)
    nh=2 ; flocal1=temp + rjm*(dd(1)*(t-tau)+dd(2))*hklA(tau)
    nh=15 ; goto 200
!
!	pour dp int G3j
!
110	if (jsol.eq.msol) rjm=(1.d0-rd(jsol)*rd(msol))/ra
    if (jsol.ne.msol) rjm=-rd(jsol)*rd(msol)/ra
    nh=4 ; temp=rd(jsol)*hklA(tau)*rd(msol)
    nh=2 ; flocal1=temp + rjm*hklA(tau)
    nh=16 ; goto 200
!
!	pour dp Gi3
!
120	nh=4 ; temp=rd(isol)*hklA(tau)*rd(msol)
    if (isol.eq.msol) rim=(1.d0-rd(isol)*rd(msol))/ra
    if (isol.ne.msol) rim=-rd(isol)*rd(msol)/ra
    nh=2 ; flocal1=temp + rim*hklA(tau)
    nh=17 ; goto 200
!
200	continue
!
!	dist(1)=dist1(1)
!	dist(2)=dist1(2)
!	ra=ra1
!
	return
    end
!
!****************************************************************************************
!****************************************************************************************
!
    double precision function  soldyGA(t1,t2)
!

    implicit double precision (a-h,o-z)
	EXTERNAL flocal1
!
	common /point/dist(2),ra,rd(2),rdn,eta(2)
	common /Dpoint/ ddist(2),dra
	common /coefA/ee(9,2,2),dpee(8,2,2,2),dd(6),ff(4),gg(6),v(2),bt
	common /indicesol/ isol,jsol,msol,nh,nd
!
	dimension dist1(2)
!
!	dist1(1)=dist(1)
!	dist1(2)=dist(2)
!	ra1=ra
!
!	dist(1)=ddist(1)
!	dist(2)=ddist(2)
!	ra=dra
!
	soldyGA=0.d0
	i=isol
	j=jsol
	zmax=t2
	zmin=ra/v(1)
	if (t1.gt.zmin) zmin=t1
	if (i.eq.3 .and. j.ne.3) goto 19 !G3j
	if (i.ne.3 .and. j.eq.3) goto 29 !Gi3
	if (i.eq.3 .and. j.eq.3) goto 39 !G33
!
9	goto (12,10,11) (nd+2) !nd:-1,0,1
!
!	------------- int Gij -------------
!
12	if (t2.le.(ra/v(1))) goto 17
	if (t1.le.(ra/v(1))) then
		nh=7
		call dqag(flocal1,zmin,t2,t2,res)
        soldyGA=soldyGA + res
	else
        nh=7 ; call dqag(flocal1,ra/v(1),t2,t2,res)
        soldyGA=soldyGA + res
        nh=7 ; call dqag(flocal1,ra/v(1),t1,t1,res)
        soldyGA=soldyGA - res
	endif
17  soldyGA=soldyGA + ee(3,i,j)*(t2**3.D0-t1**3.D0)/6.D0+ee(4,i,j)*(t2**2.D0-t1**2.D0)/2.D0+&
			ee(5,i,j)*(t2-t1)+ ee(9,i,j)*(Dexp(-2.D0*bt*t2)-Dexp(-2.D0*bt*t1))
	if (t2.le.(ra/v(2))) goto 18
	if (t1.le.(ra/v(2))) then
		x2=dsqrt(t2**2.D0-(ra/v(2))**2.D0)
		soldyGA=soldyGA + ee(6,i,j)*(v(2)**2.D0)*(x2**3.D0)/(6.D0*ra)+ee(8,i,j)*&
				(t2*dlog(t2+x2)-x2-t2*dlog(ra/v(2)))
	else
		x2=dsqrt(t2**2.D0-(ra/v(2))**2.D0)
		x1=dsqrt(t1**2.D0-(ra/v(2))**2.D0)
		soldyGA=soldyGA +ee(8,i,j)*(t2*dlog(t2+x2)-x2-t2*dlog(ra/v(2))) &
                        -ee(8,i,j)*(t1*dlog(t1+x1)-x1-t1*dlog(ra/v(2))) &
						+ee(6,i,j)*(v(2)**2.D0)*(x2**3.D0-x1**3.D0)/(6.D0*ra)
	endif
18	goto 100
!
!	------------- Gij -------------
!
10	if (t2.le.(ra/v(1))) goto 13
	nh=8
	call dqag(flocal1,zmin,zmax,0.d0,res)
	soldyGA=soldyGA + res
13	soldyGA=soldyGA + ee(3,I,J)*(t2**2.D0-t1**2.D0)/2.D0 + ee(4,I,J)*(t2-t1)-&
					  ee(5,I,J)*(Dexp(-2.D0*bt*t2)-Dexp(-2.D0*bt*t1))
	if (t2.le.(ra/v(2))) goto 14
	if (t1.le.(ra/v(2))) then
		x2=dsqrt (t2**2.D0-(ra/v(2))**2.D0)
		soldyGA=soldyGA + ee(8,I,J)*dlog(t2+x2)+ee(6,I,J)*v(2)**2.D0*t2*x2/(2.D0*ra)-&
				ee(8,I,J)*dlog(ra/v(2))
	else
		x2=dsqrt(t2**2.D0-(ra/v(2))**2.D0)
		x1=dsqrt(t1**2.D0-(ra/v(2))**2.D0)
		soldyGA=soldyGA + ee(8,I,J)*dlog((t2+x2)/(t1+x1))+ee(6,I,J)*v(2)**2.D0*(t2*x2-t1*x1)/(2.D0*ra)
	endif
14	goto 100
!
!	------------- dt Gij -------------
!
11	if (t2.le.(ra/v(1))) goto 15
	nh=2 ; soldyGA=soldyGA + ee(1,i,j)*(hklA(t2)-hklA(t1))
	nh=1 ; soldyGA=soldyGA + ee(2,i,j)*(hklA(t2)-hklA(t1))
15	soldyGA=soldyGA + ee(3,i,j)*(t2-t1)-ee(4,i,j)*(Dexp(-2.D0*bt*t2)-Dexp(-2.D0*bt*t1))
	if (t2.le.(ra/v(2))) goto 16
	nh=6 ; soldyGA=soldyGA + ee(6,i,j)*(hklA(t2)-hklA(t1))
	nh=5 ; soldyGA=soldyGA + ee(7,i,j)*(hklA(t2)-hklA(t1))
16	goto 100
!
19	goto (26,24,20) (nd+3)
!
!	------------- int int G3j -------------
!
26	if (t2.le.(ra/v(1))) goto 27
	if (t1.le.(ra/v(1))) then
		nh=9 ; call dqag(flocal1,zmin,t2,t2,res)
        soldyGA=soldyGA + res*rd(j)
	else
        nh=9 ; call dqag(flocal1,ra/v(1),t2,t2,res)
        soldyGA=soldyGA + res*rd(j)
        nh=9 ; call dqag(flocal1,ra/v(1),t1,t1,res)
        soldyGA=soldyGA - res*rd(j)
	endif
27  soldyGA=soldyGA + rd(j)/ra*(dd(3)*(t2**2.D0-t1**2.D0)/2.D0+dd(5)*(t2-t1)+dd(6)*&
			(Dexp(-2.D0*bt*t2)-Dexp(-2.D0*bt*t1)))
	goto 100
!
!	------------- int G3j -------------
!
24	if (t2.le.(ra/v(1))) goto 25
	nh=2 ; call dqag(flocal1,zmin,zmax,0.d0,res)
	soldyGA=soldyGA + dd(1)*res
	nh=2 ; soldyGA=soldyGA + dd(2)*(hklA(t2)-hklA(t1))
	soldyGA=soldyGA*rd(j)
25  soldyGA=soldyGA + rd(j)/ra*(dd(3)*(t2-t1)-dd(5)*(Dexp(-2.D0*bt*t2)-Dexp(-2.D0*bt*t1)))
    goto 100
!
!	------------- G3j -------------
!
20	if (t2.le.(ra/v(1))) goto 22
	nh=2 ; soldyGA=soldyGA + dd(1)*(hklA(t2)-hklA(t1))
	nh=3 ; soldyGA=soldyGA + dd(2)*(hklA(t2)-hklA(t1))
	soldyGA=soldyGA*rd(j)
22  soldyGA=soldyGA + rd(j)/ra*dd(4)*(Dexp(-2.D0*bt*t2)-Dexp(-2.D0*bt*t1))
    if (t1.eq.0.0) soldyGA=soldyGA+rd(j)/ra*(dd(3)+dd(4))
    goto 100
!
29	goto (30,31,32) (nd+1)
!
!	------------- Gi3 -------------
!
30	if (t2.le.(ra/v(1))) goto 34
	nh=2 ; call dqag(flocal1,zmin,zmax,0.d0,res)
	soldyGA=soldyGA + rd(i)*ff(1)*res
34  soldyGA=soldyGA + rd(i)/ra*(ff(2)*(t2-t1)-ff(3)*(Dexp(-2.D0*bt*t2)-Dexp(-2.D0*bt*t1)))
    goto 100
!
!	------------- dt Gi3 -------------
!
31	if (t2.le.(ra/v(1))) goto 37
    nh=2 ; soldyGA=soldyGA + rd(i)*ff(1)*(hklA(t2)-hklA(t1))
37  soldyGA=soldyGA - rd(i)/ra*ff(2)*(Dexp(-2.D0*bt*t2)-Dexp(-2.D0*bt*t1))
    goto 100
!
!	------------- ddt Gi3 -------------
!
32	if (t2.le.(ra/v(1))) goto 35
	nh=3 ; soldyGA=soldyGA + rd(i)*ff(1)*(hklA(t2)-hklA(t1))
35  soldyGA=soldyGA + rd(i)/ra*ff(4)*(Dexp(-2.D0*bt*t2)-Dexp(-2.D0*bt*t1))
    if (t1.eq.0.d0) soldyGA=soldyGA + rd(i)/ra*ff(4)
    goto 100
!
!
39	goto (45,40,41) (nd+2)
!
!	------------- int G33 -------------
!
45  if (t2.le.(ra/v(1))) goto 46
	if (t1.le.(ra/v(1))) then
		nh=10 ; call dqag(flocal1,zmin,t2,t2,res)
        soldyGA=soldyGA+res
	else
        nh=10 ; call dqag(flocal1,ra/v(1),t2,t2,res)
        soldyGA=soldyGA + res
        nh=10 ; call dqag(flocal1,ra/v(1),t1,t1,res)
        soldyGA=soldyGA - res
	endif
46  soldyGA=soldyGA - gg(6)*dlog(ra)*(Dexp(-2.D0*bt*t2)-Dexp(-2.D0*bt*t1))
    goto 100
!
!	------------- G33 -------------
!
40	if (t2.le.(ra/v(1))) goto 43
	if (t1.le.(ra/v(1))) then
		nh=11 ; call dqag(flocal1,zmin,t2,t2,res)
        soldyGA=soldyGA + res
    else
        nh=11 ; call dqag(flocal1,ra/v(1),t2,t2,res)
        soldyGA=soldyGA + res
        nh=11 ; call dqag(flocal1,ra/v(1),t1,t1,res)
        soldyGA=soldyGA - res
    endif
43  soldyGA=soldyGA + gg(2)*dlog(ra)*(Dexp(-2.D0*bt*t2)-Dexp(-2.D0*bt*t1))
    if (t1.eq.0.d0) soldyGA=soldyGA + gg(2)*dlog(ra)
    goto 100
!
!	------------- dt G33 -------------
!
41  if (t2.le.(ra/v(1))) goto 44
	nh=1 ; soldyGA=soldyGA + gg(3)*(hklA(t2)-hklA(t1))
    if (t1.le.(ra/v(1))) then
		nh=12 ; call dqag(flocal1,zmin,t2,t2,res)
        soldyGA=soldyGA + gg(4)*res
	else
        nh=12 ; call dqag(flocal1,ra/v(1),t2,t2,res)
        soldyGA=soldyGA + gg(4)*res
        nh=12 ; call dqag(flocal1,ra/v(1),t1,t1,res)
        soldyGA=soldyGA - gg(4)*res
	endif
44  soldyGA=soldyGA + gg(4)*dlog(ra)*(Dexp(-2.D0*bt*t2)-Dexp(-2.D0*bt*t1))
    if (t1.eq.0.d0) soldyGA=soldyGA + gg(4)*dlog(ra)
!
100	continue
!
!	dist(1)=dist1(1)
!	dist(2)=dist1(2)
!	ra=ra1
!
	return
	end
!
!
!****************************************************************************************
!
!              d�riv�e spatiales des solutions dynamiques en d�placement
!
!****************************************************************************************
!
	subroutine dpsoldyG(t1,t2,dpsolG,k)
	implicit double precision (a-h,o-z)
	dimension dpsolG(3,3,2),dist1(2)
	EXTERNAL flocal1
!
	common /point/ dist(2),ra,rd(2),rdn,eta(2)
	common /Dpoint/ ddist(2),dra
	common /coefA/ ee(9,2,2),dpee(8,2,2,2),dd(6),ff(4),gg(6),v(2),bt
	common /indicesol/ isol,jsol,msol,nh,nd
	common /stability/ anew1,anew2,wil,intb
!
	do i=1,3
		do j=1,3
			do m=1,2
				dpsolG(i,j,m)=0.d0
			enddo
		enddo
	enddo
!
!	dist1(1)=dist(1)
!	dist1(2)=dist(2)
!	ra1=ra
!
!	dist(1)=ddist(1)
!	dist(2)=ddist(2)
!	ra=dra
!
	zmax=t2
	zmin=ra/v(1)
	if (t1.gt.zmin) zmin=t1
!
	do i=1,3
		isol=i
		do j=1,3
			jsol=j
			if (i.eq.3 .and. j.ne.3) goto 140
			if (i.ne.3 .and. j.eq.3) goto 160
			if (i.eq.3 .and. j.eq.3) goto 180
!
!	dp int Gij
!
			if (k.eq.2) goto 900
			if (intb.eq.1) goto 240
			do m=1,2
				msol=m
				temp=0.d0
				if (t2.le.(ra/v(1))) goto 220
				if (t1.le.(ra/v(1))) then
					nh=13
					call dqag(flocal1,zmin,t2,t2,res)
					temp=res
				else
					nh=13
					call dqag(flocal1,ra/v(1),t2,t2,res)
					temp=res
					nh=13
					call dqag(flocal1,ra/v(1),t1,t1,res)
					temp=temp-res
				endif
				nh=2
				x=hklA(t2)-hklA(t1)
!         nh=2 ; call dqag(flocal1,zmin,zmax,0.d0,res) !� d�activer
!         x=x+2*bt*res !� d�activer
				temp=temp - ee(2,i,j)*rd(m)*x/(v(1)**2)
220				temp=temp+ dpee(3,i,j,m)*(t2**3-t1**3)/6&
                 + dpee(4,i,j,m)*(t2**2-t1**2)/2&
                 + dpee(5,i,j,m)*(t2-t1)&
                 + dpee(5,i,j,m)/(2*bt)*(Dexp(-2*bt*t2)-Dexp(-2*bt*t1))
				 if (t2.le.(ra/v(2))) goto 230
				 x2=dsqrt(t2**2-(ra/v(2))**2)
				 if (t1.le.(ra/v(2))) then
					temp=temp&
					     +(dpee(7,i,j,m)-dpee(6,i,j,m)*ra/2-ee(6,i,j)*rd(m)/2)&
					     *(t2*dlog(t2+x2)-x2-t2*dlog(ra/v(2)))&
!     $     *(t2*dlog(t2+x2)-x2-(t2+ra/v(2))*dlog(ra/v(2)))
					     +(dpee(6,i,j,m)-ee(6,i,j)*rd(m)/ra)*(v(2)**2)*(x2**3)/(6*ra)&
					     -ee(7,i,j)*rd(m)*x2/ra
				else
					x1=dsqrt(t1**2-(ra/v(2))**2)
					temp=temp&
						+(dpee(7,i,j,m)-dpee(6,i,j,m)*ra/2-ee(6,i,j)*rd(m)/2)&
						*(t2*dlog(t2+x2)-t1*dlog(t1+x1)-x2+x1-(t2-t1)*dlog(ra/v(2)))&
						+(dpee(6,i,j,m)-ee(6,i,j)*rd(m)/ra)&
						*v(2)**2*(x2**3-x1**3)/(6*ra)&
						-ee(7,i,j)*rd(m)/ra*(x2-x1)
				endif
230				dpsolG(i,j,m)=temp
			enddo
			goto 900
!
!     dp Gij
!
240			do m=1,2
				msol=m
				temp=0.d0
				if (t2.le.(ra/v(1))) goto 120
				nh=14
				call dqag(flocal1,zmin,zmax,0.d0,res)
				temp=res
				nh=3
				x=hklA(t2)-hklA(t1)
				nh=2
				x=x+2*bt*(hklA(t2)-hklA(t1))
				temp=temp - ee(2,i,j)*rd(m)*x/(v(1)**2)
120				temp=temp + dpee(3,i,j,m)*(t2-t1)*(t2+t1)/2&
		             + dpee(4,i,j,m)*(t2-t1)&
		             - dpee(5,i,j,m)*(dexp(-2*bt*t2)-dexp(-2*bt*t1))
				if (t2.le.(ra/v(2))) goto 130
				x2=dsqrt(t2**2-(ra/v(2))**2)
				if (t1.le.(ra/v(2))) then
					temp=temp&
					     +(dpee(7,i,j,m)-dpee(6,i,j,m)*ra/2-ee(6,i,j)*rd(m)/2)&
					     *dlog((t2+x2)*v(2)/ra)&
					     +(dpee(6,i,j,m)-ee(6,i,j)*rd(m)/ra)*v(2)**2*t2*x2/(2*ra)&
					     -ee(7,i,j)*rd(m)*t2/(ra*x2)
				else
					x1=dsqrt(t1**2-(ra/v(2))**2)
					temp=temp&
					     +(dpee(7,i,j,m)-dpee(6,i,j,m)*ra/2-ee(6,i,j)*rd(m)/2)&
						 *dlog((t2+x2)/(t1+x1))&
						 +(dpee(6,i,j,m)-ee(6,i,j)*rd(m)/ra)&
						 *v(2)**2*(t2*x2-t1*x1)/(2*ra)&
						 -ee(7,i,j)*rd(m)/ra*(t2/x2-t1/x1)
				endif
130				dpsolG(i,j,m)=temp
			enddo
			goto 900
!
!    dp int int G3j
!
140			if (k.eq.2) goto 900
			if (intb.eq.1) goto 260
			do m=1,2
				msol=m
				temp=0.d0
				if (t2.le.(ra/v(1))) goto 250
				if (t1.le.(ra/v(1))) then
					nh=15
					call dqag(flocal1,zmin,t2,t2,res)
					temp=res
				else
					nh=15
					call dqag(flocal1,ra/v(1),t2,t2,res)
					temp=res
					nh=15
					call dqag(flocal1,ra/v(1),t1,t1,res)
					temp=temp-res
				endif
250				x =dd(3)*(t2**2-t1**2)/2 + dd(5)*(t2-t1)&
				   +dd(6)*(dexp(-2*bt*t2)-dexp(-2*bt*t1))
				if (j.eq.m) temp=temp + x*(1.d0-2*rd(j)*rd(m))/(ra**2)
				if (j.ne.m) temp=temp - x*2*rd(j)*rd(m)/(ra**2)
				dpsolG(i,j,m)=temp
			enddo
			goto 900
!
!	dp int G3j
!
260			do m=1,2
				msol=m
				temp=0.d0
				if (t2.le.(ra/v(1))) goto 150
				nh=16
				call dqag(flocal1,zmin,zmax,0.d0,res)
				temp=res*dd(1)
				nh=4
				temp=temp + rd(j)*dd(2)*(hklA(t2)-hklA(t1))*rd(m)
				nh=2
				x=dd(2)*(hklA(t2)-hklA(t1))
				if (j.eq.m) temp=temp + x*(1.d0-rd(j)*rd(m))/ra
				if (j.ne.m) temp=temp - x*rd(j)*rd(m)/ra
150				x=dd(3)*(t2-t1)-dd(5)*(dexp(-2*bt*t2)-dexp(-2*bt*t1))
				if (j.eq.m) temp=temp + x*(1.d0-2*rd(j)*rd(m))/(ra**2)
				if (j.ne.m) temp=temp - x*2*rd(j)*rd(m)/(ra**2)
				dpsolG(i,j,m)=temp
			enddo
			goto 900
!
!	dp Gi3
!
160			if (k.eq.1) goto 180
			if (intb.eq.1 .or. intb.eq.3) goto 280
			do m=1,2
				msol=m
				temp=0.d0
				if (t2.le.(ra/v(1))) goto 270
				nh=17
				call dqag(flocal1,zmin,zmax,0.d0,res)
				temp=res*ff(1)
270				x=ff(2)*(t2-t1)-ff(3)*(dexp(-2*bt*t2)-dexp(-2*bt*t1))
				if (i.eq.m) temp=temp + x*(1.d0-2*rd(i)*rd(m))/(ra**2)
				if (i.ne.m) temp=temp - x*2*rd(i)*rd(m)/(ra**2)
				dpsolG(i,j,m)=temp
			enddo
			goto 900
!
!     dp dt Gi3
!
280			do m=1,2
				msol=m
				temp=0.d0
				if (t2.le.(ra/v(1))) goto 170
				nh=4
				temp=rd(i)*ff(1)*(hklA(t2)-hklA(t1))*rd(m)
				nh=2
				x = ff(1)*(hklA(t2)-hklA(t1))
				if (i.eq.m) temp=temp + x*(1.d0-rd(i)*rd(m))/ra
				if (i.ne.m) temp=temp - x*rd(i)*rd(m)/ra
170				x=-ff(2)*(dexp(-2*bt*t2)-dexp(-2*bt*t1))
				if (i.eq.m) temp=temp + x*(1.d0-2*rd(i)*rd(m))/(ra**2)
				if (i.ne.m) temp=temp - x*2*rd(i)*rd(m)/(ra**2)
				dpsolG(i,j,m)=temp
			enddo
			goto 900
!
!     dp int G33
!
180			if (k.eq.1) goto 900
			if (intb.eq.1 .or. intb.eq.3) goto 300
			do m=1,2
				msol=m
				temp=0.d0
				if (t2.le.(ra/v(1))) goto 290
				nh=2
				call dqag(flocal1,zmin,zmax,0.d0,res)
				nh=2
				x=gg(3)*(hklA(t2)-hklA(t1)) + gg(5)*res
				temp=-rd(m)*x/v(1)**2
290				temp=temp - gg(6)*(dexp(-2*bt*t2)-dexp(-2*bt*t1))*rd(m)/ra
				dpsolG(i,j,m)=temp
			enddo
			goto 900
!
!     dp G33
!
300			do m=1,2
				msol=m
				temp=0.d0
				if (t2.le.(ra/v(1))) goto 190
				nh=3
				x=gg(3)*(hklA(t2)-hklA(t1))
				nh=2
				x=x+gg(5)*(hklA(t2)-hklA(t1))
				temp=-rd(m)*x/v(1)**2
190				temp=temp + gg(2)*(dexp(-2*bt*t2)-dexp(-2*bt*t1))*rd(m)/ra
				if (t1.eq.0.d0) temp=temp+gg(2)*rd(m)/ra
				dpsolG(i,j,m)=temp
			enddo
			goto 900
!
900			continue
		enddo
	enddo
!
!	dist(1)=dist1(1)
!	dist(2)=dist1(2)
!	ra=ra1
!
	return
	end
!
!
!****************************************************************************************
!
!								   SATDYNFSOLGLAP subroutine
!
!	This subroutine calculates the SATurated-DYNamique Fundamental SOLutions of G
!	numerically (by inverse of Laplace)
!	ot dynamic solutions of displacement in the case of incompressible composents.

!****************************************************************************************
!
	subroutine satdynfsolGLap(itm,itime,dtime,solPEdyG)
	implicit double precision (a-h,o-z)
!
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /indicesol/ isol,jsol,msol,nh,nd
	common /stability/ anew1,anew2,wil,intb
!
	dimension solPEdyG(kbem,kbem)
	doUBLE COMPLEX solLapG
	EXTERNAL solLapG
!
!
	t3=(itime-itm)*dtime+wil*dtime
	if (itime.gt.1) then
		t2=(itime-itm-1)*dtime+wil*dtime
		if (intb.eq.3) t2=(itime-itm-2)*dtime+wil*dtime
	else
        t2=0.d0
    endif
!
	if (itime.gt.1 .and. intb.eq.2) then
		t1=(itime-itm-2)*dtime+wil*dtime
	else
		t1=0.d0
	endif
!
	do isol=1,3
		do jsol=1,3
!	constant (mixte) approximation
			nd=0
			if (isol.eq.3 .and. jsol.ne.3) nd=-1 !G3j
			if (isol.ne.3 .and. jsol.eq.3) nd=1  !Gi3
!	lin�aire approximation
			if (intb.eq.2) nd=nd-1
!
			call dlainv(solLapG,t3,0.d0,5.d-8,5.d-16,500,rslt)
			temp1=rslt
			if (t2.ne.0.d0) then
				call dlainv(solLapG,t2,0.d0,5.d-8,5.d-16,500,rslt)
				temp1=temp1-rslt
			endif
			temp2=0.d0
			if (intb.eq.2 .and. itime.gt.1) then
				call dlainv(solLapG,t2,0.d0,5.d-8,5.d-16,500,rslt)
				temp2=rslt
				if (t1.ne.0.d0) then
					call dlainv(solLapG,t1,0.d0,5.d-8,5.d-16,500,rslt)
					temp2=temp2-rslt
				endif
			endif
!
			if (intb.eq.1 .OR. intb.eq.3) then !?????
				solPEdyG(isol,jsol)=temp1
			else
				solPEdyG(isol,jsol)=(temp1-temp2)/dtime
			endif
		enddo
	enddo
!
	return
	end
!
!
!****************************************************************************************
!
!								   SATDYNFSOLTLAP subroutine
!
!	This subroutine calculates the SATurated-DYNamique Fundamental SOLutions of G
!	numerically (by inverse of Laplace)
!	ot dynamic solutions of traction in the case of incompressible composents.
!
!****************************************************************************************
!
	subroutine satdynfsolTLap(itm,itime,dtime,solPEdyT)
	implicit double precision (a-h,o-z)
!
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /indicesol/ isol,jsol,msol,nh,nd
	common /stability/ anew1,anew2,wil,intb
!
	dimension solPEdyT(kbem,kbem)
	doUBLE COMPLEX solLapT
	EXTERNAL solLapT
!
!
	t3=(itime-itm)*dtime+wil*dtime
	t2=(itime-itm-1)*dtime+wil*dtime
	if (itime.eq.1) t2=0.d0
	if (itime.gt.1 .and. intb.ne.1) then
		t1=(itime-itm-2)*dtime+wil*dtime
	else
		t1=0.d0
	endif
	nd=0
!
	do isol=1,3
		do jsol=1,3
			call dlainv(solLapT,t3,0.d0,5.d-8,5.d-16,500,rslt)
			temp1=rslt
			if (t2.ne.0.d0) then
				call dlainv(solLapT,t2,0.d0,5.d-8,5.d-16,500,rslt)
				temp1=temp1-rslt
			endif
			temp2=0.d0
			if (intb.ne.1 .and. itime.gt.1) then
				call dlainv(solLapT,t2,0.d0,5.d-8,5.d-16,500,rslt)
				temp2=rslt
				if (t1.ne.0.d0) then
					call dlainv(solLapT,t1,0.d0,5.d-8,5.d-16,500,rslt)
					temp2=temp2-rslt
				endif
			endif
			if (intb.eq.1) then
				solPEdyT(isol,jsol)=temp1
			else
				solPEdyT(isol,jsol)=(temp1-temp2)/dtime
			endif
		enddo
	enddo
!
	return
	end
!
!
!****************************************************************************************
!
!     SOLUTIONS FONDAMENTALES NUMERIQUE PORO-ELASTO-DYNAMIQUE
!     par transformation Laplace inverse
!
!****************************************************************************************
!
!     calcul des coefficients et leurs d�riv�es spaciales
!
	subroutine satcoefficientsLap
	implicit double precision (a-h,o-z)
!
	common /point/dist(2),ra,rd(2),rdn,eta(2)
	common /mate/zlamda,zmuy,znuu,znu,zalpha,zrho,zrhof,zm,zk,zr,zn
	common /coefL/aa(2,2),dpaa(2,2,2),bb(2,2),dpbb(2,2,2),cc(2,2)
!
	pi=DACOS(-1.d0)
!
!	Aij et dpAijm
!
	temp=2*pi*ra**3
    aa(1,1)=(2*dist(1)**2 - ra**2)/temp
    aa(1,2)=2*dist(1)*dist(2)/temp
    aa(2,1)=aa(1,2)
    aa(2,2)=(2*dist(2)**2 - ra**2)/temp
    temp=temp*ra
    dpaa(1,1,1)=(4*ra*dist(1)-6*rd(1)*dist(1)**2+ra**2*rd(1))/temp
    dpaa(1,2,1)=(2*ra*dist(2)-6*rd(1)*dist(1)*dist(2))/temp
    dpaa(2,1,1)=dpaa(1,2,1)
    dpaa(2,2,1)=(-6*rd(1)*dist(2)**2+ra**2*rd(1))/temp
    dpaa(2,2,2)=(4*ra*dist(2)-6*rd(2)*dist(2)**2+ra**2*rd(2))/temp
    dpaa(1,2,2)=(2*ra*dist(1)-6*rd(2)*dist(1)*dist(2))/temp
    dpaa(2,1,2)=dpaa(1,2,2)
    dpaa(1,1,2)=(-6*rd(2)*dist(1)**2+ra**2*rd(2))/temp
!
!	Bij et dpBijm
!
	temp=2*pi*ra**2
    bb(1,1)=dist(1)**2/temp
    bb(2,2)=dist(2)**2/temp
    bb(1,2)=dist(1)*dist(2)/temp
    bb(2,1)=bb(1,2)
    temp=temp*ra
    dpbb(1,1,1)=(2*ra*dist(1)-2*rd(1)*dist(1)**2)/temp
    dpbb(1,2,1)=(ra*dist(2)-2*rd(1)*dist(1)*dist(2))/temp
    dpbb(2,1,1)=dpbb(1,2,1)
    dpbb(2,2,1)=-2*rd(1)*dist(2)**2/temp
    dpbb(2,2,2)=(2*ra*dist(2)-2*rd(2)*dist(2)**2)/temp
    dpbb(1,2,2)=(ra*dist(1)-2*rd(2)*dist(1)*dist(2))/temp
    dpbb(2,1,2)=dpbb(1,2,2)
    dpbb(1,1,2)=-2*rd(2)*dist(1)**2/temp
!
!	Cij (dpCijm=0)
!
	temp=1./(2*pi*zmuy)
    cc(1,1)=temp
    cc(2,2)=temp
    cc(1,2)=0.d0
    cc(2,1)=0.d0
!
	return
	end
!
!------------------------------------------------------------------
!     d�riv�e spatiales des solutions dynamiques en d�placement
!------------------------------------------------------------------
!
	subroutine dpsolLapG(s,dpsolG)
	implicit double precision (a-h,o-z)
	doUBLE COMPLEX s,dpsolG(3,3,2),lbd(3),tmp,lbdH,alp2,hklL
!
	common /point/ dist(2),ra,rd(2),rdn,eta(2)
	common /mate/zlamda,zmuy,znuu,znu,zalpha,zrho,zrhof,zm,zk,zr,zn
	common /coefL/ aa(2,2),dpaa(2,2,2),bb(2,2),dpbb(2,2,2),cc(2,2)
	common /indicesol/ isol,jsol,msol,nh,nd
!
	pi=DACOS(-1.d0)
	alp2=zalpha-zrhof*zk*s
	lbdH=zrho*s**2/(zlamda+2*zmuy)
!
	lbd(3)=dsqrt(zrho/zmuy)*s
	tmp=Cdsqrt((s/zk*(1.d0/zm+zalpha*alp2/(zlamda+2*zmuy))+lbdH)**2-4*s*lbdH/(zk*zm))
	lbd(1)=Cdsqrt(0.5d0*(s/zk*(1.d0/zm+zalpha*alp2/(zlamda+2*zmuy))+lbdH+tmp))
	lbd(2)=Cdsqrt(0.5d0*(s/zk*(1.d0/zm+zalpha*alp2/(zlamda+2*zmuy))+lbdH-tmp))
!
	i=isol
	j=jsol
	if (isol.eq.3 .and. jsol.ne.3) goto 30
	if (isol.ne.3 .and. jsol.eq.3) goto 50
	if (isol.eq.3 .and. jsol.eq.3) goto 70
!
!	dp Gij
!
	do i=1,2
		do m=1,2
			dpsolG(i,j,m) = &
            ( aa(i,j)*hklL(lbd(1),m,4)+bb(i,j)*hklL(lbd(1),m,2)&
             +dpaa(i,j,m)*hklL(lbd(1),0,3)+dpbb(i,j,m)*hklL(lbd(1),0,1))&
             *lbd(1)**2*(lbdH-lbd(2)**2)/tmp&
             -(aa(i,j)*hklL(lbd(2),m,4)+bb(i,j)*hklL(lbd(2),m,2)&
             +dpaa(i,j,m)*hklL(lbd(2),0,3)+dpbb(i,j,m)*hklL(lbd(2),0,1))&
             *lbd(2)**2*(lbdH-lbd(1)**2)/tmp&
             -(aa(i,j)*hklL(lbd(3),m,4)+bb(i,j)*hklL(lbd(3),m,2)&
             +dpaa(i,j,m)*hklL(lbd(3),0,3)+dpbb(i,j,m)*hklL(lbd(3),0,1))&
             *lbd(3)**2
         dpsolG(i,j,m)= dpsolG(i,j,m)/(zrho*s**3)&
                        + cc(i,j)*hklL(lbd(3),m,2)/s
		enddo
	enddo
	goto 90
!
!	dp G3j
!
30	do m=1,2
		if (j.eq.m) rim=(1.d0-rd(j)*rd(m))/ra
		if (j.ne.m) rim=- rd(j)*rd(m)/ra
		dpsolG(i,j,m)=&
		      -( lbd(1)**2*(rim*hklL(lbd(1),0,3)+rd(j)*hklL(lbd(1),m,4))&
			  -lbd(2)**2*(rim*hklL(lbd(2),0,3)+rd(j)*hklL(lbd(2),m,4)))&
			  *alp2/(2*pi*(zlamda+2*zmuy)*zk*tmp)
	enddo
	goto 90
!
!	dp Gi3
!
50	do i=1,2
		do m=1,2
			if (i.eq.m) rim=(1.d0-rd(i)*rd(m))/ra
			if (i.ne.m) rim= -rd(i)*rd(m)/ra
			dpsolG(i,j,m)=&
			      -( lbd(1)**2*(rim*hklL(lbd(1),0,3)+rd(i)*hklL(lbd(1),m,4))&
				  -lbd(2)**2*(rim*hklL(lbd(2),0,3)+rd(i)*hklL(lbd(2),m,4)))&
				  *zalpha/(2*pi*(zlamda+2*zmuy)*zk*tmp*s)
		enddo
	enddo
	goto 90
!
!	dp G33
!
70	do m=1,2
		dpsolG(i,j,m)=&
		       ((lbd(1)**2-lbdH)*hklL(lbd(1),m,2)&
			   -(lbd(2)**2-lbdH)*hklL(lbd(2),m,2))&
			   /(2*pi*zk*tmp*s)
	enddo
!
90	return
	end
!
!--------------------------------------
!     solution Laplacienne en d�placement
!--------------------------------------
!
	doUBLE COMPLEX function solLapG(s)
	implicit double precision (a-h,o-z)
	doUBLE COMPLEX s,lbd(3),tmp,alp2,lbdH,hklL
!
	common /point/ dist(2),ra,rd(2),rdn,eta(2)
	common /mate/zlamda,zmuy,znuu,znu,zalpha,zrho,zrhof,zm,zk,zr,zn
	common /coefL/ aa(2,2),dpaa(2,2,2),bb(2,2),dpbb(2,2,2),cc(2,2)
	common /indicesol/ isol,jsol,msol,nh,nd
!
	solLapG=(0.d0,0.d0)
	pi=DACOS(-1.d0)
	alp2=zalpha-zrhof*zk*s
	lbdH=zrho*s**2/(zlamda+2*zmuy)
!
	lbd(3)=dsqrt(zrho/zmuy)*s
	tmp=Cdsqrt((s/zk*(1.d0/zm+zalpha*alp2/(zlamda+2*zmuy))+lbdH)**2- 4*s*lbdH/(zk*zm))
	lbd(1)=Cdsqrt(0.5d0*(s/zk*(1.d0/zm+zalpha*alp2/(zlamda+2*zmuy))+lbdH+tmp))
	lbd(2)=Cdsqrt(0.5d0*(s/zk*(1.d0/zm+zalpha*alp2/(zlamda+2*zmuy))+lbdH-tmp))
!
	i=isol
	j=jsol
	if (i.eq.3 .and. j.ne.3) goto 30
	if (i.ne.3 .and. j.eq.3) goto 50
	if (i.eq.3 .and. j.eq.3) goto 70
!
!	Gij
!
	solLapG= (aa(i,j)*hklL(lbd(1),0,3)+bb(i,j)*hklL(lbd(1),0,1))&
              *(lbd(1)**2)*(lbdH-lbd(2)**2)/tmp&
			  -(aa(i,j)*hklL(lbd(2),0,3)+bb(i,j)*hklL(lbd(2),0,1))&
			  *(lbd(2)**2)*(lbdH-lbd(1)**2)/tmp&
			  -(aa(i,j)*hklL(lbd(3),0,3)+bb(i,j)*hklL(lbd(3),0,1))&
			  *lbd(3)**2
	 solLapG=solLapG/(zrho*s**3)+cc(i,j)*hklL(lbd(3),0,1)/s
	 if (nd.eq.-1) solLapG=solLapG/s
	 goto 100
!
!	G3j
!
30	solLapG=-(lbd(1)**2*hklL(lbd(1),0,3)-lbd(2)**2*hklL(lbd(2),0,3))&
              *rd(j)*alp2/(2*pi*(zlamda+2*zmuy)*zk*tmp)
	if (nd.eq.-1) solLapG=solLapG/s
	if (nd.eq.-2) solLapG=solLapG/s**2
	goto 100
!
!	Gi3
!
50	solLapG=-(lbd(1)**2*hklL(lbd(1),0,3)-lbd(2)**2*hklL(lbd(2),0,3))&
              *rd(i)*zalpha/(2*pi*(zlamda+2*zmuy)*zk*tmp*s)
	if (nd.eq.1) solLapG=solLapG*s
	goto 100
!
!	G33
!
70	solLapG=( (lbd(1)**2-lbdH)*hklL(lbd(1),0,1)&
              -(lbd(2)**2-lbdH)*hklL(lbd(2),0,1))/(2*pi*zk*tmp*s)
	if (nd.eq.-1) solLapG=solLapG/s
	if (nd.eq.1) solLapG=solLapG*s
!
100	return
	end
!
!--------------------------------------------
!     solutions Laplaciennes en traction
!--------------------------------------------
!
	doUBLE COMPLEX function solLapT(s)
	implicit double precision (a-h,o-z)
	doUBLE COMPLEX dpsolG(3,3,2),s,pres,dep,solLapG
!
	common /point/ dist(2),ra,rd(2),rdn,eta(2)
	common /mate/zlamda,zmuy,znuu,znu,zalpha,zrho,zrhof,zm,zk,zr,zn
	common /indicesol/ isol,jsol,msol,nh,nd
	common /stability/ anew1,anew2,wil,intb
!
	solLapT=(0.d0,0.d0)
	call dpsolLapG(s,dpsolG)
	i=isol
	j=jsol
	if (isol.eq.3) goto 20
!
!	Tij et dt Ti3
!
	do m=1,2
		solLapT=solLapT+(dpsolG(i,j,m)+dpsolG(m,j,i))*eta(m)
	enddo
	isol=3
	pres=solLapG(s)
	isol=i
	solLapT=solLapT*zmuy+((dpsolG(1,j,1)+dpsolG(2,j,2))*zlamda&
						   -pres*zalpha)*eta(i)
	if (j.eq.3) solLapT=solLapT*s
	goto 30
!
!	int T3j et T33
!
20	do m=1,2
		isol=m
		dep=solLapG(s)
		isol=i
		solLapT=solLapT+(dpsolG(i,j,m)+zrhof*dep*s**2)*eta(m)
	enddo
	solLapT=-solLapT*zk
	if (j.ne.3) solLapT=solLapT/s
!
!	lin�aire (mixte) approximation
!
30	if (intb.ne.1) solLapT=solLapT/s
!
	return
	end

!****************************************************************************************
!
!							           ELDYNFSOL subroutine
!
!	This subroutine calculates ELasto-DYNamique Fundamental SOLutions.
!	This subroutine is called in: EXTINeq1
!	This subroutine calls: -
!
!	variables used are:
!
!		itime	: time increment from 1 to "ntime"
!							[1] : first time increment,
!							[>1]: higher time increment
!		itm		:
!		dtime	: time step; in seismic loading for all "1:ntime", it is constant
!		dg		: displacement fundamental solution
!		df		: traction fundamental solution
!		pa		:
!		qa		:
!		ca		:
!		sa		:
!		dlt		: kronecker delta
!		vel(1)	: Vp, speed of primary(P) wave
!		vel(2)	: Vs, speed of secondary(S) wave
!		coef1	: costant coefficient used in G formulation: 1/(2*zrho*pi)
!		coef2	: costant coefficient used in F formulation: zmuy/(2*zrho*pi*ra)
!		dist(1)	: distance from collocation point to the gauss integration
!				  points on the x-plane; = xco-xp
!		dist(2)	: distance from collocation point to the gauss integration
!				  points on the y-plane; =yco-yp
!		ra		: distance from collocation point to the gauss integration points on the
!				  boundary elements
!		rd(1)	: radius derivative dr/dx1
!		rd(2)	: radius derivative dr/dx2
!		rdn		: radius derivative dr/dn
!		eta(1)	: component of the unit normal to the element in direction 1
!				  V3=Vxsi*Vz=[d y/d xsi, -d x/d xsi, 0] then
!				  eta(1)=Vxsi(y)=(e*cjab+djab)/xja1
!		eta(2)	: component of the unit normal to the element in direction 2
!				  V3=Vxsi*Vz=[d y/d xsi, -d x/d xsi, 0] then
!				  eta(2)= -Vxsi(x)=(e*ajab+bjab)/xja1
!		zlamda	: Lam�'s coefficient which is equal to K-2MU/3
!		zmuy	: Lam�'s coefficient which is equal to 3EK/(9K-E)
!		zalpha	: alpha, in saturated formulation
!		zrho	: density of body
!		zrhof	: density of fluid (water) in a saturated volume
!		zm		: 1/Q, in saturated formulation
!		zk		: water permeability in z-direction
!		anew1	: Gamma coefficient in Newmark direct integration method. this coef must be
!				  [>=0.5] because that the convergence be inconditionally stable.
!				  In this program, it is chosen 0.5. this choice corresponds to the rule of
!				  the trapezoid
!		anew2	: Beta coefficient in Newmark direct integration method. this coef must be
!				  [>=0.25] because that the convergence be inconditionally stable.
!				  In this program, it is chosen 0.25. this choice corresponds to the rule of
!				  the trapezoid
!		wil		: TETA coefficient in TETA-Wilson method for the stability of numerical
!				  solution [>1]. the coefficient plays a numerical damping role of artificial
!				  wave propagation due to the perturbations. When TETA=[1] it is meaning that
!				  TETA-Wilson method is not considered
!		intb	: interpolation function of time [1,2,3];
!								[1]: constant: two fields of displacement and stress (u,t)
!											   are supposed to be constant in each time
!											   interval
!								[2]: linear	 : two fields of displacement and stress (u,t)
!											   vary linearly in each time interval
!								[3]: mixte	 : displacement field (u) vary linearly while
!											   stress field (t) is constant in each time
!											   interval
!		fca1	: Acosh (Vm*itime*dtime/r)
!		fca2	: Acosh (Vm*(itime-1)*dtime/r)
!		fea		: fca1-fca2: [Acosh (Vm*itime*dtime/r)-Acosh (Vm*(itime-1)*dtime/r)]
!		fta		: [ itime*SQRT[itime**2-(r/Vm/dtime)**2]-
!					(itime-1)*SQRT[(itime-1)**2-(r/Vm/dtime)**2] ]
!		fva		: [ itime/SQRT[itime**2-(r/Vm/dtime)**2]-
!					(itime-1)/SQRT[(itime-1)**2-(r/Vm/dtime)**2] ]
!
!	INPUT	:	itm,itime,dtime
!
!
!	OUTPUT	:	dg,df
!
!
!****************************************************************************************
!
	subroutine eldynfsol(itm,itime,dtime,dg,df)
!
	implicit double precision (a-h,o-z)
!
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /point/ dist(2),ra,rd(2),rdn,eta(2)
	common /mate/zlamda,zmuy,znuu,znu,zalpha,zrho,zrhof,zm,zk,zr,zn
	common /stability/ anew1,anew2,wil,intb
!
	dimension dg(kbem,kbem),df(kbem,kbem),pa(2),qa(2),ca(2),sa(2),dlt(2,2),vel(2)
!
!
	vel(1)=dsqrt ((zlamda+2*zmuy)/zrho)
	vel(2)=dsqrt (zmuy/zrho)
	pi=DACOS (-1.d0)
	dlt(1,1)=1.d0
	dlt(2,2)=1.d0
	dlt(1,2)=0.d0
	dlt(2,1)=0.d0
	coef1=1.d0/(2*zrho*pi)
	coef2=zmuy/(2*zrho*pi*ra)
!
!
!	------------------- initializing -------------------
!
	do i=1,2
		do j=1,2
			dg(i,j)=0.d0
			df(i,j)=0.d0
		enddo
	enddo
!
!
!	----------- evaluating traction solution df,linear approximation -----------
!
	if (intb.eq.2 .OR. intb.eq.3) then
		do ia=1,2
			pa(ia)=fpa(itime,itm,dtime,ra,vel,ia)
			qa(ia)=fqa(itime,itm,dtime,ra,vel,ia)
			do i=1,2
				do j=1,2
					aa1=fa1(i,j,zmuy,zlamda,eta,rd,rdn)
					aa2=fa2(i,j,eta,rd,rdn,dlt)
					aa3=fa3(i,j,eta,rd,rdn,dlt)
					s1=(-aa1*dlt(1,ia)+aa3*dlt(2,ia))/vel(ia)**2
					s2=((-1)**ia)*aa2*(dtime/ra)**2
					df(i,j)=df(i,j)+(s1*pa(ia)-(2.d0/3.d0)*s2*qa(ia))*coef2
				enddo
			enddo
		enddo
	endif
!
!
!	----------- evaluating displacement solution dg, linear approximation -----------
!
	if (intb.eq.2) then
		do ia=1,2
			ca(ia)=fca(itime,itm,dtime,ra,vel,ia)
			pa(ia)=fpa(itime,itm,dtime,ra,vel,ia)
			qa(ia)=fqa(itime,itm,dtime,ra,vel,ia)
			sa(ia)=fsa(itime,itm,dtime,ra,vel,ia)
			do i=1,2
				do j=1,2
					aa1=dlt(i,j)/2.d0
					aa2=(dlt(i,j)-2*rd(i)*rd(j))
					s1=aa1
					s2=((-1)**ia)/2.d0*aa2*(vel(ia)*dtime/ra)**2
					s3=((-1)**ia)/3.d0*aa2*(vel(ia)*dtime/ra)**2
					s4=((-1)**ia)*(dlt(i,j)*dlt(2,ia)-rd(i)*rd(j))
					dg(i,j)=dg(i,j)+(s1*ca(ia)+s2*sa(ia)-s3*qa(ia)-s4*pa(ia))*coef1/vel(ia)**2
				enddo
			enddo
		enddo
	endif
!
!
!	----------- evaluating displacement solution dg, mixte approximation -----------
!
	if (intb.eq.3) then
		do ia=1,2
			ca(ia)=fda(itime,itm,dtime,ra,vel,ia)
			sa(ia)=fra(itime,itm,dtime,ra,vel,ia)
			do i=1,2
				do j=1,2
					aa1=dlt(i,j)/2.d0
					aa2=(dlt(i,j)-2*rd(i)*rd(j))
					s1=aa1
					s2=((-1)**ia)/2.d0*aa2*(vel(ia)*dtime/ra)**2
					dg(i,j)=dg(i,j)+(s1*ca(ia)+s2*sa(ia))*coef1/2.d0/vel(ia)**2
				enddo
			enddo
		enddo
	endif
!
!
!	----------- evaluating traction solution df,constante approximation -----------
!
	if (intb.eq.1) then
		do ia=1,2
			pa(ia)=fva(itime,itm,dtime,ra,vel,ia)
			qa(ia)=fta(itime,itm,dtime,ra,vel,ia)
			do i=1,2
				do j=1,2
					aa1=fa1(i,j,zmuy,zlamda,eta,rd,rdn)
					aa2=fa2(i,j,eta,rd,rdn,dlt)
					aa3=fa3(i,j,eta,rd,rdn,dlt)
					s1=(-aa1*dlt(1,ia)+aa3*dlt(2,ia))/vel(ia)**2
					s2=((-1)**ia)*2.d0*aa2*(dtime/ra)**2
					df(i,j)=df(i,j)+(s1*pa(ia)-s2*qa(ia))*coef2
				enddo
			enddo
		enddo
	endif
!
!
!	----------- evaluating displacement solution dg, constante approximation -----------
!
	if (intb.eq.1) then
		do ia=1,2
			ca(ia)=fea(itime,itm,dtime,ra,vel,ia)
			sa(ia)=fta(itime,itm,dtime,ra,vel,ia)
			do i=1,2
				do j=1,2
					aa1=dlt(i,j)/2.d0
					aa2=(dlt(i,j)-2*rd(i)*rd(j))
					s1=aa1
					s2=((-1)**ia)/2.d0*aa2*(vel(ia)*dtime/ra)**2
					dg(i,j)=dg(i,j)+(s1*ca(ia)+s2*sa(ia))*coef1/vel(ia)**2
				enddo
			enddo
		enddo
	endif
!
	return
	end
!
!
!	----------------------------------
!	coefficient A1
!
	double precision function fa1(i,j,zmuy,zlamda,eta,rd,rdn)
!
	implicit double precision (a-h,o-z)
!
	dimension eta(2),rd(2)
!
	fa1=(zlamda/zmuy)*eta(j)*rd(i)+2*rd(i)*rd(j)*rdn
!
	return
	end
!
!
!	----------------------------------
!	coefficient A2
!
	double precision function fa2(i,j,eta,rd,rdn,dlt)
!
	implicit double precision (a-h,o-z)
!
	dimension eta(2),rd(2),dlt(2,2)
!
	fa2=eta(i)*rd(j)+eta(j)*rd(i)+rdn*(dlt(i,j)-4*rd(i)*rd(j))
!
	return
	end
!
!
!	----------------------------------
!	coefficient A3
!
	double precision function fa3(i,j,eta,rd,rdn,dlt)
!
	implicit double precision (a-h,o-z)
!
	dimension eta(2),rd(2),dlt(2,2)
!
	fa3=rdn*(2*rd(i)*rd(j)-dlt(i,j))-eta(i)*rd(j)
!
	return
	end
!
!
!	----------------------------------
!	coefficient Pk
!
	double precision function fpa(itime,itm,dtime,ra,vel,ia)
!
	implicit double precision (a-h,o-z)
	double precision n1,n2,n3
!
	common /stability/ anew1,anew2,wil,intb
	dimension vel(2)
!
	x=vel(ia)*dtime
	xx=ra/(vel(ia)*dtime)
	n1=(itime-itm+wil)
	n2=(itime-itm-1+wil)
	if (itime.eq.itm) n2=0.d0
	n3=(itime-itm-2+wil)
!
	pa1=0.d0
	if (ra.lt.x*n1) pa1=dsqrt (n1**2-xx**2)
	pa2=0.d0
	if (ra.lt.x*n2) pa2=dsqrt (n2**2-xx**2)
	pa3=0.d0
	if (ra.lt.x*n3) pa3=dsqrt (n3**2-xx**2)
!
	fpa=pa1-2*pa2+pa3
!
	return
	end
!
!
!	----------------------------------
!	coefficient Qk
!
	double precision function fqa(itime,itm,dtime,ra,vel,ia)
!
	implicit double precision (a-h,o-z)
	double precision n1,n2,n3
!
	common /stability/ anew1,anew2,wil,intb
	dimension vel(2)
!
	x=vel(ia)*dtime
	xx=ra/(vel(ia)*dtime)
	n1=(itime-itm+wil)
	n2=(itime-itm-1+wil)
	if (itime.eq.itm) n2=0.d0
	n3=(itime-itm-2+wil)
!
	qa1=0.d0
	if (ra.lt.x*n1) qa1=(n1**2-xx**2)**1.5
	qa2=0.d0
	if (ra.lt.x*n2) qa2=(n2**2-xx**2)**1.5
	qa3=0.d0
	if (ra.lt.x*n3) qa3=(n3**2-xx**2)**1.5
!
	fqa=qa1-2*qa2+qa3
!
	return
	end
!
!
!	----------------------------------
!	coefficient Sk
!
	double precision function fsa(itime,itm,dtime,ra,vel,ia)
!
	implicit double precision (a-h,o-z)
	double precision n1,n2,n3
!
	common /stability/ anew1,anew2,wil,intb
	dimension vel(2)
!
	x=vel(ia)*dtime
	xx=ra/(vel(ia)*dtime)
	n1=(itime-itm+wil)
	n2=(itime-itm-1+wil)
	if (itime.eq.itm) n2=0.d0
	n3=(itime-itm-2+wil)
!
	sa1=0.0
	if (ra.lt.x*n1) sa1=(n1**2)*dsqrt (n1**2-xx**2)
	sa2=0.0
	if (ra.lt.x*n2) sa2=(n2**2)*dsqrt (n2**2-xx**2)
	sa3=0.0
	if (ra.lt.x*n3) sa3=(n3**2)*dsqrt (n3**2-xx**2)
!
	fsa=sa1-2*sa2+sa3
!
	return
	end
!
!
!	----------------------------------
!	coefficient Rk
!
!
	double precision function fra(itime,itm,dtime,ra,vel,ia)
!
	implicit double precision (a-h,o-z)
	double precision n1,n2
!
	common /stability/ anew1,anew2,wil,intb
	dimension vel(2)
!
	x=vel(ia)*dtime
	xx=ra/(vel(ia)*dtime)
	n1=(itime-itm+wil)
	n2=(itime-itm-2+wil)
	if (itime.eq.itm) n2=0.d0
!
	sa1=0.d0
	if (ra.lt.x*n1) sa1=n1*dsqrt (n1**2-xx**2)
	sa2=0.d0
	if (ra.lt.x*n2) sa2=n2*dsqrt (n2**2-xx**2)
!
	fra=sa1-sa2
!
	return
	end
!
!
!	----------------------------------
!	coefficient Tk
!
	double precision function fta(itime,itm,dtime,ra,vel,ia)
!
	implicit double precision (a-h,o-z)
	double precision n1,n2
!
	common /stability/ anew1,anew2,wil,intb
	dimension vel(2)
!
	x=vel(ia)*dtime
	xx=ra/(vel(ia)*dtime)
	n1=(itime-itm+wil)
	n2=(itime-itm-1+wil)
	if (itime.eq.itm) n2=0.d0
	sa1=0.d0
	if (ra.lt.x*n1) sa1=n1*dsqrt (n1**2-xx**2)
	sa2=0.d0
	if (ra.lt.x*n2) sa2=n2*dsqrt (n2**2-xx**2)
!
	fta=sa1-sa2
!
	return
	end
!
!
!	----------------------------------
!	coefficient Vk
!
	double precision function fva(itime,itm,dtime,ra,vel,ia)
!
	implicit double precision (a-h,o-z)
	double precision n1,n2
!
	common /stability/ anew1,anew2,wil,intb
	dimension vel(2)
!
	x=vel(ia)*dtime
	xx=ra/(vel(ia)*dtime)
	n1=(itime-itm+wil)
	n2=(itime-itm-1+wil)
	if (itime.eq.itm) n2=0.d0
	sa1=0.d0
	if (ra.lt.x*n1) sa1=n1/dsqrt (n1**2-xx**2)
	sa2=0.d0
	if (ra.lt.x*n2) sa2=n2/dsqrt (n2**2-xx**2)
!
	fva=sa1-sa2
!
	return
	end
!
!
!	----------------------------------
!
	double precision function fca(itime,itm,dtime,ra,vel,ia)
!
	implicit double precision (a-h,o-z)
	double precision n1,n2,n3
!
	common /stability/ anew1,anew2,wil,intb
	dimension vel(2)
!
	x=vel(ia)*dtime
	xx=ra/(vel(ia)*dtime)
	n1=(itime-itm+wil)
	n2=(itime-itm-1+wil)
	if (itime.eq.itm) n2=0.d0
	n3=(itime-itm-2+wil)
!
	fca1=0.d0
	if (ra.lt.x*n1) fca1=n1*dacosh(n1/xx)
	fca2=0.d0
	if (ra.lt.x*n2) fca2=n2*dacosh(n2/xx)
	fca3=0.d0
	if (ra.lt.x*n3) fca3=n3*dacosh(n3/xx)
!
	fca=fca1-2*fca2+fca3
!
	return
	end
!
!
!	----------------------------------
!
	double precision function fda(itime,itm,dtime,ra,vel,ia)
!
	implicit double precision (a-h,o-z)
	double precision n1,n2
!
	common /stability/ anew1,anew2,wil,intb
	dimension vel(2)
!
	x=vel(ia)*dtime
	xx=ra/(vel(ia)*dtime)
	n1=(itime-itm+wil)
	n2=(itime-itm-2+wil)
	if (itime.eq.itm) n2=0.d0
!
	fca1=0.d0
	if (ra.lt.x*n1) fca1=dacosh(n1/xx)
	fca2=0.d0
	if (ra.lt.x*n2) fca2=dacosh(n2/xx)
!
	fda=fca1-fca2
!
	return
	end
!
!
!	----------------------------------
!
	double precision function fea(itime,itm,dtime,ra,vel,ia)
!
	implicit double precision (a-h,o-z)
	double precision n1,n2
!
	common /stability/ anew1,anew2,wil,intb
	dimension vel(2)
!
	x=vel(ia)*dtime
	xx=ra/(vel(ia)*dtime)
	n1=(itime-itm+wil)
	n2=(itime-itm-1+wil)
	if (itime.eq.itm) n2=0.d0
!
	fca1=0.d0
	if (ra.lt.x*n1) fca1=DACOSH (n1/xx)
	fca2=0.d0
	if (ra.lt.x*n2) fca2=DACOSH (n2/xx)
!
	fea=fca1-fca2
!
	return
	end
!**************************************************************************************************
!
!								           Version HYBRID 2009
!
!**************************************************************************************************
!
!
!               Gatmiri B. & Kamalian M. & Nguyen K.V. & Dehghan K. & Maghoul P.
!
!     two dimensional plane strain or axisymetric finite element / boundary element program for:
!     dynamic analysis of (poro)elastic dry / saturated / unsaturated media.
!
!     eight noded drained / consolidated finite elements and three noded boundary elements.
!
!     linear, hyperbolic (el.) and prevost (ep.) models are used for drained and consolidated
!	  finite elements. linear (poro)elastic model is used for boundary elements.
!
!     program for hydraulic anisotropic soils
!     program for different boundary conditions
!
!
!**************************************************************************************************

	PROGRAM Hybrid
!
	double precision a(30000000)
	integer lmax
	integer int_time1,int_time2,int_time
!
	int_time1 = TIME()
	call control(a,lmax)
!
	write (*,*)
	write (*,*)
	write (*,*)' Program has run succesfully'
!
	int_time2 = TIME()
	int_time=int_time2-int_time1
!
	write (*,*)
	write (*,*)
	write (*,*)'  Duration of Analysis (sec)  : ',int_time
	write (7,*)'  Duration of Analysis (sec)  : ',int_time
	close (7)
	write (*,*)
	write (*,*)
!
	end


!----------------------------------------------------
!     les fontions h et leurs d�riv�es spaciales
!----------------------------------------------------
!
	doUBLE COMPLEX function hklL(z,msol,nh)
	implicit double precision (a-h,o-z)
!
	doUBLE COMPLEX z,besk0,besk1
	common /point/ dist(2),ra,rd(2),rdn,eta(2)
!
	goto (1,2,3,4) nh
!	hk1
1	call cbesk01(z*ra,besk0,besk1)
	hklL=besk0
	goto 5
!
!	hk1,m
2	call cbesk01(z*ra,besk0,besk1)
	hklL=-besk1*z*rd(msol)
	goto 5
!
!	hk2
3	call cbesk01(z*ra,besk0,besk1)
	hklL=besk1/z
	goto 5
!
!	hk2,m
4	call cbesk01(z*ra,besk0,besk1)
	hklL=-(besk0+besk1/(z*ra))*rd(msol)
!
5	return
	end

!****************************************************************************************
!
	subroutine ghmatd1INT(itm,k1,x,xint,sm8d,sm8c,sm8u,ie3d1,hd1int,gd1int,sig3u)
!
	implicit double precision (a-h,o-z)
!
	common /a/ nnp,ifem,ibem,idyn,ieaq,nr
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /a3/ ne8d1,ne8c1,ne8u1,ne3u1,nload1,nbcx1,nbcy1,nbcw1,nbca1,m8d1,m8c1,m8u1,&
				nbed11,nbed21,nbed31,nbeint1,ne8dd,ne8cc,ne8uu
	common /b/ m8d,m8c,m8u,m8dep(10),m8cep(10),m8uep(10)
	common /c/ nload,ntime,itmax,itime,iter,iload,time,dtime,dtimei,dtimem,xinc
	common /e/ gammaw,gammas,gammaa,grav,atmp
	common /h/ mdof,mdofn,ls,ls1,isolv
	common /mate/zlamda,zmuy,znuu,znu,zalpha,zrho,zrhof,zm,zk,zr,zn
	common /umate/ d(3,3),dfs(3),xn,xkw,cww,cwg,xka,cgg,cgw,ev,sat,g1,rs,rw,raa,rmix,xm1,&
				   bt,et,xlambda,xmu,henry
	common /jab/ajab,bjab,cjab,djab,xelm(3),yelm(3),xp,yp
	common /unbem/ nelmu,ie3u(3,50,5)
!
	dimension x(2,nnp),xint(2,nbeint),sm8d(60,m8d1),sm8c(60,m8c1),sm8u(60,m8u1),&
			  ie3d1(npbmax,nbed11),hd1int(kbem*nbeint,n12gh),gd1int(kbem*nbeint,n12gh),&
			  ehd1int(kbem,3*kbem),egd1int(kbem,3*kbem),xd(npbmax),yd(npbmax),&
			  xmtf(n12gh,n12gh),id(4,nnp),sige(5),sig3u(5,500,nbed11)
!
!
!	------------------- initializing -------------------
!
	n=nnpbed1(k1)
	n1=kbem*n
	ndim=kbem*nbeint
	Iint=1
!
	do i=1,ndim
		do j=1,n1
			hd1int(i,j)=0.d0
			gd1int(i,j)=0.d0
		enddo
	enddo
!
	if (itime.eq.1.OR.itm.eq.1) goto 55
	call INTghsaved(itm,ndim,n1,gd1int,hd1int,'ghd1int')
	goto 2000
55	continue

!
!	--------- material & geometrical properties ---------
!
	mtype=ie3d1(n+1,k1)
	matparam: SELECT CASE (kbem)
				CASE (2)	! dry mat
					zrho=sm8d(11,mtype)/grav
					zmuy=3*sm8d(3,mtype)*sm8d(7,mtype)/(9*sm8d(7,mtype)-sm8d(3,mtype))
					zlamda=sm8d(7,mtype)-2.d0*zmuy/3.d0
				CASE (3) 	! saturated mat
!					zrho=sm8c(11,mtype)/grav
					zmuy=sm8c(3,mtype) !3*sm8c(3,mtype)*sm8c(7,mtype)/(9*sm8c(7,mtype)-sm8c(3,mtype))
					ZNU=sm8c(7,mtype)
					ZNUU=sm8c(14,mtype)
!					zlamda=sm8c(7,mtype)-2.d0*zmuy/3.d0
!					zrhof=gammaw/grav
					zk=sm8c(15,mtype) !sm8c(13,mtype)/gammaw
!					zm=1.d0/sm8c(15,mtype)
					zalpha=sm8c(13,mtype) !sm8c(14,mtype)
			  end SELECT matparam
71	continue
!
!	coordinate point of 3 noded boundary elements
!
	do ll=1,n
		ii=ie3d1(ll,k1)
		xd(ll)=x(1,ii)
        yd(ll)=x(2,ii)
	enddo
	xd(n+1)=xd(1)  ! in the last elm, the 3d node corresponds to the 1st node of the 1st elm
    yd(n+1)=yd(1)
!
!
!	--- selecting collocation points & boundary elements ---
!
	do 400 ll=1,nbeint
		nelm=0
		do 400 i=1,n-1,2
			nelm=nelm+1
			xp=xint(1,ll)
			yp=xint(2,ll)
			xelm(1)=xd(i)
			yelm(1)=yd(i)
			xelm(2)=xd(i+1)
			yelm(2)=yd(i+1)
			xelm(3)=xd(i+2)
			yelm(3)=yd(i+2)
!
			if (ibem.eq.3) then
				do kk=1,5
					sige(kk)=sig3u(kk,nelm,k1)
				enddo
				call dmatUelas (mtype,sige,sm8u)
			endif
!
!
!	------ investigating the occurance of real and/or apparent kind of singularity ------
!	----------------- in the element and also the number of subelements -----------------
!
			call sngtp(n,i,ll,Iint,nodo,nsub,nptg)
!
!
!	------ calculating egd1 & ehd1 in each of the elements ------
!
			call INTextineq(itm,itime,dtime,idyn,nodo,nsub,nptg,ehd1int,egd1int)
!
!
!	------ assembling the ehd1 & egd1 matrices into the hd1 & gd1 ones ------
!
			ii=3*kbem
			do 390 k=1,kbem
				do 380 j=1,ii
					if (i.ne.(n-1))    goto 370
					if (j.le.(2*kbem)) goto 370 !last node  last elm. (N1)
					hd1int(kbem*(ll-1)+k,j-2*kbem)=hd1int(kbem*(ll-1)+k,j-2*kbem)+ehd1int(k,j)
					gd1int(kbem*(ll-1)+k,j-2*kbem)=gd1int(kbem*(ll-1)+k,j-2*kbem)+egd1int(k,j)
					goto 380
!
370					hd1int(kbem*(ll-1)+k,kbem*(i-1)+j)=hd1int(kbem*(ll-1)+k,kbem*(i-1)+j)+&
													   ehd1int(k,j)
					gd1int(kbem*(ll-1)+k,kbem*(i-1)+j)=gd1int(kbem*(ll-1)+k,kbem*(i-1)+j)+&
													   egd1int(k,j)
!
380				continue
390			continue
400	continue
!
!	----------------------------------------------------------------
!
	call INTghsaved(itm,ndim,n1,gd1int,hd1int,'ghd1int')
!
2000 return
	end


!
!****************************************************************************************
!
	subroutine INTghsaved(itm,ndim,n1,gd1,hd1,text)
!
!
	implicit double precision (a-h,o-z)
	character text*7, nts*5
!
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /c/ nload,ntime,itmax,itime,iter,iload,time,dtime,dtimei,dtimem,xinc
!
	dimension hd1(kbem*nbeint,n12gh),gd1(kbem*nbeint,n12gh)
!
	numb=itime-itm+1
	call NumToStr(nts,numb)
	open (2,file=text//nts)
	do i=1,ndim
		do j=1,n1
			if (itm.eq.1) then
				write (2,*) gd1(i,j),hd1(i,j)
			else
				read (2,*) gd1(i,j),hd1(i,j)
			endif
		enddo
	enddo
	close (2)
!
	return
	end
!
!****************************************************************************************
!
!****************************************************************************************
!
	subroutine INTextineq(itm,itime,dtime,idyn,nodo,nsub,nptg,ehd1,egd1)
!
	implicit double precision (a-h,o-z)
!
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /mate/zlamda,zmuy,znuu,znu,zalpha,zrho,zrhof,zm,zk,zr,zn
	common /umate/ d(3,3),dfs(3),xn,xkw,cww,cwg,xka,cgg,cgw,ev,sat,g1,rs,rw,raa,rmix,xm1,&
				   bt,et,henry
	common /jab/ajab,bjab,cjab,djab,xelm(3),yelm(3),xp,yp
!
	dimension ehd1(kbem,3*kbem),egd1(kbem,3*kbem),sehd1(kbem,3*kbem),segd1(kbem,3*kbem)
!
!
!	------------------- initializing -------------------
!
	nn=3*kbem
	do i=1,kbem
		do j=1,nn
			ehd1(i,j)=0.d0
			egd1(i,j)=0.d0
		enddo
	enddo
!
!	-------------- geometrical parameters --------------
!
	e11=-1.d0
	e22=1.d0
	ajab=xelm(3)-2.d0*xelm(2)+xelm(1)
	bjab=(xelm(3)-xelm(1))/2.d0
	cjab=yelm(3)-2.d0*yelm(2)+yelm(1)
	djab=(yelm(3)-yelm(1))/2.d0
	deltae=(e22-e11)/nsub
	ee1=e11
!
!	-- (intehd1 & integd1) in each element --
!
	do k=1,nsub
		ee2=ee1+deltae
!
!	-- (sehstd1 & segstd1) and (sehd1 & segd1) in each sub-region --
!
		call INTextineq1(itm,itime,dtime,nodo,nptg,ee1,ee2,sehd1,segd1)
!
!	-- assembling (ehstd1 & egstd1) and (ehd1 & egd1) in each element --
!
		do 50 i=1,kbem
			do 40 j=1,nn
				if (idyn.eq.0) goto 40
35				ehd1(i,j)=ehd1(i,j)+sehd1(i,j)
				egd1(i,j)=egd1(i,j)+segd1(i,j)
40			continue
50		continue
		ee1=ee2
	enddo
!
200	return
	end
!
!****************************************************************************************
!
!****************************************************************************************
!
	subroutine INTextineq1(itm,itime,dtime,nodo,nptg,ee1,ee2,sehd1,segd1)
!
	implicit double precision (a-h,o-z)
!
	common /a/ nnp,ifem,ibem,idyn,ieaq,nr
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /mate/zlamda,zmuy,znuu,znu,zalpha,zrho,zrhof,zm,zk,zr,zn
	common /umate/ d(3,3),dfs(3),xn,xkw,cww,cwg,xka,cgg,cgw,ev,sat,g1,rs,rw,raa,rmix,xm1,&
				   bt,et,henry
	common /jab/ajab,bjab,cjab,djab,xelm(3),yelm(3),xp,yp
	common /point/dist(2),ra,rd(2),rdn,eta(2)
	common /Kcor/kcorf,mtime,rmtol,c1,c2
	common /ocqm/ EPScqm,taw,Icqm,Lcqm,NPcqm
!
	dimension sehd1(kbem,3*kbem),segd1(kbem,3*kbem),idf(kbem),stg(kbem,kbem),stf(kbem,kbem),&
			  dg(kbem,kbem),df(kbem,kbem),f(3),gi(10),wi(10)
!
!	determine of factors integral
!
	dimension g10(10),w10(10),g9(9),w9(9),g8(8),w8(8),g7(7),w7(7),g6(6),w6(6),g5(5),w5(5),&
			  g4(4),w4(4),g3(3),w3(3)
	dimension sai(4,4),sItn(kbem,3*kbem) !?????
!
!
	save g10,w10,g9,w9,g8,w8,g7,w7,g6,w6,g5,w5,g4,w4,g3,w3
!
	data g10 / 0.973906528517172,-0.973906528517172,&
			   0.865063366688985,-0.865063366688985,&
			   0.679409568299024,-0.679409568299024,&
			   0.433395394129247,-0.433395394129247,&
			   0.148874338981631,-0.148874338981631 /
	data w10 / 0.066671344308688,0.066673344308688,&
			   0.149451349150581,0.149451349150581,&
			   0.219086362515982,0.219086362515982,&
			   0.269266719309996,0.269266719309996,&
			   0.295524224714753,0.295524224714753 /
!
	data g9 / 0.968160239507626,-0.968160239507626,&
			  0.836031107326636,-0.836031107326636,&
			  0.613371432700590,-0.613371432700590,&
			  0.324253423403809,-0.324253423403809,&
			  0.d0 /
	data w9 / 0.081274388361574,0.081274388361574,&
			  0.180648160694857,0.180648160694857,&
			  0.260610696402935,0.260610696402935,&
			  0.312347077040003,0.312347077040003,&
			  0.330239355001260 /
!
	data g8 / 0.960289856497536,-0.960289856497536,&
			  0.796666477413627,-0.796666477413627,&
			  0.525532409916329,-0.525532409916329,&
			  0.183434642495650,-0.183434642495650 /
	data w8 / 0.101228536290376,0.101228536290376,&
			  0.222381034453374,0.222381034453374,&
			  0.313706645877887,0.313706645877887,&
			  0.362683783378362,0.362683783378362 /
!
	data g7 / 0.949107912342759,-0.949107912342759,&
			  0.741531185599394,-0.741531185599394,&
			  0.405845151377397,-0.405845151377397,&
			  0.d0 /
	data w7 / 0.129484966168870,0.129484966168870,&
			  0.279705391489277,0.279705391489277,&
			  0.381830050505119,0.381830050505119,&
			  0.417959183673469 /
!
	data g6 / 0.932469514203152,-0.932469514203152,&
			  0.661209386466265,-0.661209386466265,&
			  0.238619186083197,-0.238619186083197 /
	data w6 / 0.171324492379170,0.171324492379170,&
			  0.360761573048139,0.360761573048139,&
			  0.467913934572691,0.467913934572691 /
!
	data g5 / 0.906179845938664,-0.906179845938664,&
			  0.538469310105683,-0.538469310105683,&
			  0.d0 /
	data w5 / 0.236926885056189,0.236926885056189,&
			  0.478628670499366,0.478628670499366,&
			  0.568888888888889 /
!
	data g4 / 0.861136311594953,-0.861136311594953,&
			  0.339981043584856,-0.339981043584856 /
	data w4 / 0.347854845137454,0.347854845137454,&
			  0.652145154862546,0.652145154862546 /
!
	data g3 / 0.774596669241483,-0.774596669241483,&
			  0.d0 /
	data w3 / 0.555555555555556,0.555555555555556,&
			  0.888888888888889 /
!
	data gi / 10*0.d0 /,wi / 10*0.d0 /
!
!
	goto (3,4,5,6,7,8,9,10) (nptg-2)
3	do i=1,nptg
		gi(i)=g3(i)
		wi(i)=w3(i)
	enddo
	goto 20
4	do i=1,nptg
		gi(i)=g4(i)
		wi(i)=w4(i)
	enddo
	goto 20
5	do i=1,nptg
		gi(i)=g5(i)
		wi(i)=w5(i)
	enddo
	goto 20
6	do i=1,nptg
		gi(i)=g6(i)
		wi(i)=w6(i)
	enddo
	goto 20
7	do i=1,nptg
		gi(i)=g7(i)
		wi(i)=w7(i)
	enddo
	goto 20
8	do i=1,nptg
		gi(i)=g8(i)
		wi(i)=w8(i)
	enddo
	goto 20
9	do i=1,nptg
		gi(i)=g9(i)
		wi(i)=w9(i)
	enddo
	goto 20
10	do i=1,nptg
		gi(i)=g10(i)
		wi(i)=w10(i)
	enddo
20	continue
!
!
!	------------------- initializing -------------------
!
	nn=3*kbem
	do i=1,kbem
		do j=1,nn
			segd1(i,j)=0.d0
			sehd1(i,j)=0.d0
		enddo
	enddo
	xja2=0.5d0*(ee2-ee1) !deltae/2=1/nsub
!
!
!	------------------- numerical integration -------------------
!	-------------------------------------------------------------
!
	do 600 i=1,nptg
!
!	---------- evaluating coefficients ----------
!
		e=0.5d0*(1-gi(i))*ee1+0.5d0*(1+gi(i))*ee2 ! coordinate of gauss point in subregion
!
!	---------- shape functions at integration points ----------
!
		f(1)=e*(e-1.d0)*0.5d0
		f(2)=1.d0-e**2
		f(3)=e*(e+1.d0)*0.5d0
!
!	------ geometrical properties at integration points -------
!
		xco=xelm(1)*f(1)+xelm(2)*f(2)+xelm(3)*f(3)
		yco=yelm(1)*f(1)+yelm(2)*f(2)+yelm(3)*f(3)
!
		dist(1)=xco-xp
		dist(2)=yco-yp
		xja1=dsqrt ((e*ajab+bjab)**2+(e*cjab+djab)**2)
		eta(1)=(e*cjab+djab)/xja1
		eta(2)=-(e*ajab+bjab)/xja1
		ra=dsqrt (dist(1)**2+dist(2)**2)
		rd(1)=dist(1)/ra	!e1=n1=dr/dx1=r1/r
		rd(2)=dist(2)/ra	!e2=n2=dr/dx2=r2/r
		rdn= rd(1)*eta(1)+rd(2)*eta(2)
		coef=xja1*xja2*wi(i)
!
!	----------- compute segd1 & sehd1 (dynamic case) ----------
!
		if (idyn.eq.0) goto 500
		dynfundsol: SELECT CASE (ibem)
			CASE (1)
				if (Icqm.eq.0) then
					call eldynfsol(itm,itime,dtime,dg,df)
				else
					call owcqm(dg,df)
				endif
			CASE (2)
				if (kfsol.eq.0.and.Icqm.eq.1) call owcqm(dg,df)
				if (kfsol.eq.1) then
					call satcoefficientsLap ! Aij,Bij,Cij & dAij,dBij
					call satdynfsolGLap(itm,itime,dtime,dg) ! Gij
					call satdynfsolTLap(itm,itime,dtime,df) ! Hij
				else if (kfsol.eq.2) then
					call satcoefficientsAlt
					call satdynfsolGAlt(itm,itime,dtime,dg)
					call satdynfsolTAlt(itm,itime,dtime,df)
				endif
			CASE (3)
					call owcqm(dg,df)
			end SELECT dynfundsol
!
		do k=1,3
			do j1=1,kbem
				idf(j1)=kbem*(k-1)+j1
				do j2=1,kbem
					segd1(j2,idf(j1))= segd1(j2,idf(j1))+dg(j2,j1)*f(k)*coef
			        sehd1(j2,idf(j1))= sehd1(j2,idf(j1))+df(j2,j1)*f(k)*coef
				enddo
			enddo
		enddo
500		continue
600	continue
!
	return
	end
!
!****************************************************************************************

!****************************************************************************************
!
!							           GHMATD1 subroutine
!
!	This subroutine calculates the g & h matrices of area D1 (finite zone)
!	This subroutine is calles in BEMT1.
!
!	These subroutines are called in this section:
!		ghsaved		:
!		sngtp		: it investigates the element's kind of singularity in the element &
!					  also the number of sub-elements
!		extineq		: calculates the ehd1 & egd1 matrices in ordinary elements
!		dinverse	:
!		mtrfrd		:
!		othermat1	:
!
!	variables used are:
!
!		k1			: counter which shows the number of each D1 (finite) zone (1:nbed1)
!		nnpbed1		: total number of nodes of D1 (finite) BEM areas
!		n			: total number of nodes of each D1 (finite) sub-zone
!		n1			: total number of doF of each D1 (finite) sub-zone; N1 = KBEM * N
!		itime		: time increment; [1]: first time increment, [>1]: higher time increment
!		itm			: n=1:N-1
!
!		hd1			: total H matrix in D1 (finite) area concerning all of the h matrices
!					  of each element (ehd1)
!		gd1			: total G matrix in D1 (finite) area concerning all of the g matrices
!					  of each element (egd1)
!		eItn		:
!		hstd1		: total H matrix in D1 (finite) area for static case concerning all of
!					  the h matrices of each element (ehstd1)
!		gstd1		: total G matrix in D1 (finite) area for static case concerning all of
!					  the g matrices of each element (egstd1)
!
!		kbem		: number of degree of freedom per node in BE zone;
!								if IBEM=[1] -> KBEM=[2] [X,Y]
!								if IBEM=[2] -> KBEM=[3] [X,Y,Pw]
!								if IBEM=[3] -> KBEM=[4] [X,Y,Pw,Pa]
!		mtype		: type (behaviour) of materials. in this study, only the linear behaviour
!					  is considered then: [1]linear,
!		zrho		: density of body
!		zmuy		: Lam�'s coefficient which is equal to 3EK/(9K-E)
!		zlamda		: Lam�'s coefficient which is equal to K-2MU/3
!		zrhof		: density of fluid (water) in a saturated volume
!		zk			: water permeability in z-direction
!		zm			: 1/Q, in saturated formulation
!		zalpha		: alpha, in saturated formulation
!
!		xd			: x-component of nodes forming the each D1 (finite) zone
!		yd			: y-component of nodes forming the each D1 (finite) zone
!		xp			: x-component of collocation point (loading point)
!		yp			: y-component of collocation point (loading point)
!		xelm(1:3)	: x-component of nodes forming the each element in each D1 (finite)
!					  sub-zone
!		yelm(1:3)	: y-component of nodes forming the each element in each D1 (finite)
!					  sub-zone
!
!		gdinv		: inverse matrice of gd1
!		gid1		: inverse matrice of gd1
!		gihd1		: gid1*hd1
!		xgid1		: M*gid1
!		xgihd1		: M*gid1*hd1
!
!
!	INPUT	:
!
!	OUTPUT	:
!
!
!****************************************************************************************
!
	subroutine ghmatd1(itm,k1,x,sm8d,sm8c,sm8u,ie3d1,hd1,gd1,gid1,gihd1,xgid1,xgihd1,&
					   sig3u)
!
	implicit double precision (a-h,o-z)
!
	common /a/ nnp,ifem,ibem,idyn,ieaq,nr
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /a3/ ne8d1,ne8c1,ne8u1,ne3u1,nload1,nbcx1,nbcy1,nbcw1,nbca1,m8d1,m8c1,m8u1,&
				nbed11,nbed21,nbed31,nbeint1,ne8dd,ne8cc,ne8uu
	common /b/ m8d,m8c,m8u,m8dep(10),m8cep(10),m8uep(10)
	common /c/ nload,ntime,itmax,itime,iter,iload,time,dtime,dtimei,dtimem,xinc
	common /d/ nbcx,nbcy,nbcw,nbca,nbc
	common /e/ gammaw,gammas,gammaa,grav,atmp
	common /h/ mdof,mdofn,ls,ls1,isolv
	common /mate/zlamda,zmuy,znuu,znu,zalpha,zrho,zrhof,zm,zk,zr,zn
	common /Dmate/cp,zt,zx,dzrho,dzrhof,dzlamda,dzmuy,dzk
	common /umate/ d(3,3),dfs(3),xn,xkw,cww,cwg,xka,cgg,cgw,ev,sat,g1,rs,rw,raa,rmix,xm1,&
				   bt,et,xlambda,xmu,henry
	common /jab/ajab,bjab,cjab,djab,xelm(3),yelm(3),xp,yp
	common /ocqm/ EPScqm,taw,Icqm,Lcqm,NPcqm
	common /unbem/ nelmu,ie3u(3,50,5)
!
	dimension x(2,nnp),sm8d(60,m8d1),sm8c(60,m8c1),sm8u(60,m8u1),ie3d1(npbmax,nbed11),&
			  hd1(n12gh,n12gh),gd1(n12gh,n12gh),&
			  xgihd1(n12gh,n12gh,nbed11),&
			  hstd1(n12gh,n12gh),gstd1(n12gh,n12gh),ehd1(kbem,3*kbem),egd1(kbem,3*kbem),&
			  ehstd1(kbem,3*kbem),egstd1(kbem,3*kbem),xd(npbmax),yd(npbmax),&
			  sige(5),sig3u(5,500,nbed11),&
			  xmtf(n12gh,n12gh),gdinv(n12gh,n12gh),&
			  gid1(n12gh,n12gh,nbed11),&
			  hd11(nnp,nnp),gihd1(n12gh,n12gh,nbed11),&
			  xgid1(n12gh,n12gh,nbed11)
	dimension eItn(kbem,3*kbem)
!
!
!	------------------- initializing -------------------
!
	n=nnpbed1(k1)
	n1=kbem*n
	Iint=0
	do i=1,n1
		do j=1,n1
			hd1(i,j)=0.d0
			gd1(i,j)=0.d0
!			eItn(i,j)=0.d0
		enddo
	enddo
!
	if (itime.eq.1.OR.itm.eq.1) goto 55
1	call ghsaved(itm,n1,gd1,hd1,'ghd1')
	goto 2000
!
55	if (itime.gt.1) goto 70
	do i=1,n1
		do j=1,n1
			hstd1(i,j)=0.d0
			gstd1(i,j)=0.d0
		enddo
	enddo
70	continue
!
!
!	--------- material & geometrical properties ---------
!
	mtype=ie3d1(n+1,k1)
	matparam: SELECT CASE (kbem)
				CASE (2)	! dry mat
					zrho=sm8d(11,mtype)/grav
					zmuy=3*sm8d(3,mtype)*sm8d(7,mtype)/(9*sm8d(7,mtype)-sm8d(3,mtype))
					zlamda=sm8d(7,mtype)-2.d0*zmuy/3.d0
				CASE (3) 	! saturated mat
					if (kfsol.eq.0.and.Icqm.eq.1) then !QUASI-STATIC PROBLEM BY CQM
						zmuy=sm8c(3,mtype)
						ZNU=sm8c(7,mtype)
						ZNUU=sm8c(14,mtype)
						zk=sm8c(15,mtype)
						zalpha=sm8c(13,mtype)
					else if (kfsol.eq.3.and.Icqm.eq.1) then !DYNAMIC PROBLEM BY CQM
						zmuy=sm8c(7,mtype)
						zlamda=(3.d0*sm8c(3,mtype)-2.d0*zmuy)/3.d0
						zrho=sm8c(11,mtype)
						zrhof=sm8c(9,mtype) !gammaw/grav
						zk=sm8c(15,mtype)
						zr=sm8c(12,mtype)
						zn=sm8c(10,mtype)
						zalpha=sm8c(13,mtype)
					else
						zrho=sm8c(11,mtype)/grav
						zrhof=gammaw/grav
						zmuy=3*sm8c(3,mtype)*sm8c(7,mtype)/(9*sm8c(7,mtype)-sm8c(3,mtype))
						zlamda=sm8c(7,mtype)-2.d0*zmuy/3.d0
						zk=sm8c(13,mtype)/gammaw
						zm=1.d0/sm8c(15,mtype)
						zalpha=sm8c(14,mtype)
					endif
			  end SELECT matparam
71	continue
!
!	cp=dsqrt ((zlamda+2*zmuy)/(zrho-zalpha*zrhof))
!
!	coordinate point of 3 noded boundary elements
!
	do ll=1,n
		ii=ie3d1(ll,k1)
		xd(ll)=x(1,ii)
        yd(ll)=x(2,ii)
	enddo
	xd(n+1)=xd(1)  ! in the last elm, the 3d node corresponds to the 1st node of the 1st elm
    yd(n+1)=yd(1)
!
!
!	--- selecting collocation points & boundary elements ---
!
	do 400 ll=1,n
		nelm=0
		do 400 i=1,n-1,2
			nelm=nelm+1
			xp=xd(ll)
			yp=yd(ll)
			xelm(1)=xd(i)
			yelm(1)=yd(i)
			xelm(2)=xd(i+1)
			yelm(2)=yd(i+1)
			xelm(3)=xd(i+2)
			yelm(3)=yd(i+2)
!
			if (ibem.eq.3) then
				do kk=1,5
					sige(kk)=sig3u(kk,nelm,k1)
				enddo
				call dmatUelas (mtype,sige,sm8u)
			endif
!
!
!	------ investigating the occurance of real and/or apparent kind of singularity ------
!	----------------- in the element and also the number of subelements -----------------
!
			call sngtp(n,ll,i,Iint,nodo,nsub,nptg)
!
!
!	------ calculating egd1 & ehd1 in each of the elements ------
!
			call extineq(itm,itime,dtime,idyn,nodo,nsub,nptg,ehd1,egd1,ehstd1,egstd1,&
						 eItn)
!
!
!	------ assembling the ehd1 & egd1 matrices into the hd1 & gd1 ones ------
!
			ii=3*kbem
			do 390 k=1,kbem
				do 380 j=1,ii
					if (i.ne.(n-1))    goto 370
					if (j.le.(2*kbem)) goto 370 !last node  last elm. (N1)
					if (idyn.eq.0)     goto 365
					hd1(kbem*(ll-1)+k,j-2*kbem)=hd1(kbem*(ll-1)+k,j-2*kbem)+ehd1(k,j)
					gd1(kbem*(ll-1)+k,j-2*kbem)=gd1(kbem*(ll-1)+k,j-2*kbem)+egd1(k,j)
!
365					if (itime.gt.1) goto 380
					hstd1(kbem*(ll-1)+k,j-2*kbem)=hstd1(kbem*(ll-1)+k,j-2*kbem)+ ehstd1(k,j)
					gstd1(kbem*(ll-1)+k,j-2*kbem)=gstd1(kbem*(ll-1)+k,j-2*kbem)+ egstd1(k,j)
					goto 380
!
370					if (idyn.eq.0) goto 375
					hd1(kbem*(ll-1)+k,kbem*(i-1)+j)=hd1(kbem*(ll-1)+k,kbem*(i-1)+j)+ehd1(k,j)
					gd1(kbem*(ll-1)+k,kbem*(i-1)+j)=gd1(kbem*(ll-1)+k,kbem*(i-1)+j)+egd1(k,j)
!
375					if (itime.gt.1) goto 380
					hstd1(kbem*(ll-1)+k,kbem*(i-1)+j)=hstd1(kbem*(ll-1)+k,kbem*(i-1)+j)+ehstd1(k,j)
					gstd1(kbem*(ll-1)+k,kbem*(i-1)+j)=gstd1(kbem*(ll-1)+k,kbem*(i-1)+j)+egstd1(k,j)
!
380				continue
390			continue
400	continue
!
!	----------------------------------------------------------------
!
	if (itime.gt.1) then
		call ghsaved(itm,n1,gd1,hd1,'ghd1')
		goto 2000
	endif
!
!
!	------ calculating diagonal coefficients of hstd1 (static case) matrice ------
!	------------------------ by RIGID BODY MOTION method -------------------------
!
	do i=1,n
		if (kbem.eq.2) then
			i1=2*i-1
			i2=2*i
		else if (kbem.eq.3) then
			i1=3*i-2
			i2=3*i-1
			i3=3*i
		else
			i1=4*i-3
			i2=4*i-2
			i3=4*i-1
			i4=4*i
		endif
!
!	the component Hij (for displacement) are strongly singular in all of the cases (dry,
!	saturated & unsaturated cases). then:
	    hstd1(i1,i1)=0.d0
	    hstd1(i2,i1)=0.d0
		hstd1(i1,i2)=0.d0
		hstd1(i2,i2)=0.d0
!
!	in saturated case, the one strongly singular fundamental solution is H33. H3j is regular
!	while Hi3 is weakly singular (which can be calculated by Gauss quadrature). then:
		if (kbem.eq.3) then
			hstd1(i3,i3)=0.d0
		endif
!
!	in unsaturated case, the one strongly singular fundamental solution are H33 & H44.
!	H3j & Hi3 - H4j & Hi4 are weakly singular while H34 & H43 are regular
!	(which can be calculated by Gauss quadrature). then:
		if (kbem.eq.4) then
			hstd1(i3,i3)=0.d0
			hstd1(i4,i4)=0.d0
		endif
!
!
		do 1600 j=1,n
			if (i.eq.j) goto 1600
			if (kbem.eq.2) then
				j1=2*j-1
				j2=2*j
			else if (kbem.eq.3) then
				j1=3*j-2
				j2=3*j-1
				j3=3*j
			else
				j1=4*j-3
				j2=4*j-2
				j3=4*j-1
				j4=4*j
			endif
		    hstd1(i1,i1)=hstd1(i1,i1)-hstd1(i1,j1)
			hstd1(i2,i1)=hstd1(i2,i1)-hstd1(i2,j1)
			hstd1(i1,i2)=hstd1(i1,i2)-hstd1(i1,j2)
			hstd1(i2,i2)=hstd1(i2,i2)-hstd1(i2,j2)
			if (kbem.eq.3) then
				hstd1(i3,i3)=hstd1(i3,i3)-hstd1(i3,j3)
			endif
			if (kbem.eq.4) then
				hstd1(i3,i3)=hstd1(i3,i3)-hstd1(i3,j3)
				hstd1(i4,i4)=hstd1(i4,i4)-hstd1(i4,j4)
			endif
1600	continue
	enddo
!
!
!	------ calculating diagonal coefficients of hd1 (dynamic case) matrice ------
!	------------------------ by RIGID BODY MOTION method -------------------------
!
	if (idyn.eq.0) goto 1750
	do i=1,n
		if (kbem.eq.2) then
			i1=2*i-1
			i2=2*i
		else if (kbem.eq.3) then
			i1=3*i-2
			i2=3*i-1
			i3=3*i
		else
			i1=4*i-3
			i2=4*i-2
			i3=4*i-1
			i4=4*i
		endif
		hd1(i1,i1)=hd1(i1,i1)+hstd1(i1,i1)
		hd1(i2,i1)=hd1(i2,i1)+hstd1(i2,i1)
		hd1(i1,i2)=hd1(i1,i2)+hstd1(i1,i2)
		hd1(i2,i2)=hd1(i2,i2)+hstd1(i2,i2)
		if (kbem.eq.3) then
			hd1(i3,i3)=hd1(i3,i3)+hstd1(i3,i3)
		endif
		if (kbem.eq.4) then
			hd1(i3,i3)=hd1(i3,i3)+hstd1(i3,i3)
			hd1(i4,i4)=hd1(i4,i4)+hstd1(i4,i4)
		endif
	enddo
1750 continue
!
!
!	------ equalizing hd1 & gd1 to hstd1 & gstd1 in static case ------
!
	if (idyn.gt.0) goto 1850
	do i=1,n1
		do j=1,n1
			hd1(i,j)=hstd1(i,j)
			gd1(i,j)=gstd1(i,j)
		enddo
	enddo
!
1850 continue
!
    call ghsaved(itm,n1,gd1,hd1,'ghd1')
!
!
!	------ calculating inverse matrix of gd1 (gid1) for transforming nodal tractions to
!	----------------------------------- nodal forces -----------------------------------
!
	do i=1,n1
		do j=1,n1
			gdinv(i,j)=gd1(i,j)
		enddo
	enddo
!
	call dinverse (gdinv,n12gh,n1,info)
!
	if (info.gt.0) then
		write (*,*) '**GHMATD1** Matrix gd1 is singular at ',info,' row'
		stop
	endif
!
	if (info.eq.-1) write (*,*) '**GHMATD1** Attention: bad condition',' of matrix gd1'
	write (*,*)'		dinverse concluded'
!
	do i=1,n1
		do j=1,n1
			gid1(i,j,k1)=gdinv(i,j)
		enddo
	enddo
!
!
!	---- calculating global M matrice for transforming nodal tractions to nodal forces ----
!
	call mtrfrd(n,xmtf,xd,yd)
!
!
!	---------------- calculating gihd1 & xgihd1 matrices ----------------
!
!
	call othermat1(k1,n,hd1,xmtf,gid1,gihd1,xgid1,xgihd1)
!
	write (*,*)'		othermat1 concluded'
!
!
2000 return
	end


!****************************************************************************************
!
!							           INIT1 subroutine
!
!	This subroutine initializes the stresses tensors, load, displacement, velocity and
!	acceleration vectors to zero.
!	This subroutine is called in: CONTROL
!	This subroutine calls: -
!
!	variables used are:
!
!		nnp		: total number of nodes in problem (FE &BE zones)
!		ne8d	: number of drained (dry) elements in FE zone
!		ne8c	: number of consolidated (saturated) elements in FE zone
!		ne8u	: number of unsaturated elements in FE zone
!		nbed1	: number of D1 BEM areas (finite areas), it can be at maximum 5 zones
!		nbed2	: number of D2 BEM areas (infinite areas), it can be at maximum 5 zones
!		nbed3	: number of D3 BEM areas (semi-infinite area), it can be at maximum 1 zone
!		nnpbed1	: number of nodes of each D1 BEM area; NNPBED1(1 to NBED1)
!		nnpbed2	: number of nodes of each D2 BEM area; NNPBED2(1 to NBED2)
!		nnpbed3	: number of nodes of one D3 BEM area
!		ie8d	:
!		ie8c	:
!		ie8u	:
!		ie3d1	: ie3d1(1:nnpbed1,1:nbed): nodes' number which form the current D1 BEM area
!										   Numbering direction for external surface close
!										   (finite) domain is ANTICLOCKWISE
!				  ie3d1(nnpbed1(i)+1,1:nbed): it shows the type of just one material which
!											  is used:
!							[1]: dry boundary elm
!							[2]: saturated boundary elm
!							[3]: unsaturated boundary elm
!		ie3d2	: ie3d2(1:nnpbed2,1:nbed): nodes' number which form the current D2 BEM area
!										   Numbering direction for internal surface open
!										   (infinite) domain is CLOCKWISE
!				  ie3d2(nnpbed2(i)+1,1:nbed): it shows the type of just one material which
!											  is used:
!							[1]: dry boundary elm
!							[2]: saturated boundary elm
!							[3]: unsaturated boundary elm
!		ie3d3	: ie3d3(1:nnpbed3): nodes' number which form the current D3 BEM area
!										   Numbering direction for semi-infinite domain is
!										   ANTICLOCKWISE
!				  ie3d3(nnpbed3+1): it shows the type of just one material which
!											  is used:
!							[1]: dry boundary elm
!							[2]: saturated boundary elm
!							[3]: unsaturated boundary elm
!				  ie3d3(nnpbed3+2):
!		init	: code describing the type of initial conditions [0,1,2,3]. when we want to
!				  do a topographical effect (ifem=0) we continue by INIT=[0]
!							[0] if the results are taken from a previous analysis
!							[1] if initial normal stresses are considered equal to the
!								weight of the upper soil layers
!							[2] if initial stresses are computed for the steady state of
!							    the system
!							[3] if initial stresses are calculated as the weight of the
!								upper soil layers, considering that the system has a mere
!							    parallelogram shape, and that the x-direction of the
!							    space-frame corresponds to the physical horizontal
!		mdof	: maximum degree of freedom by considering the boundary conditions
!		mdofn	: maximum degree of freedom without considering the boundary conditions
!
!**		sig8d	: for each node in each dry elm in FE, the stress tensors,
!				  sigX;sigY;sigXY;sigXYM;e, are stocked in this variable.
!**		sig8c	: for each node in each saturated elm in FE, the stress tensors,
!				  sigX;sigY;sigXY;sigXYM;Pw;evs, are stocked in this variable.
!**		sig8u	: for each node in each unsaturated elm in FE, the stress tensors,
!				  sigX;sigY;sigXY;sigXYM;suction;Pg_init;e, are stocked in this variable.
!**		sig8di	: in FE
!**		sig8ci	: in FE
!**		sig8ui	: in FE
!**		sigmad	: in FE
!**		sig8mac	: in FE
!**		sig8mau	: in FE
!**		sig8do	: in FE
!**		sig8co	: in FE
!**		sig8uo	: in FE
!		r		: all nodal forces (mechanical TOTRR & seismic RRR) are saved in this
!				  vector; after that in SOLVE subroutine the displacements & pressures
!				  obtained by solve the equations are saved in this vector in each itime,
!				  this vector is saved in DT vector in UPDATE1 subroutine & initialized in it
!		rr		: it is a vector which saves the mechanical loadings
!		totrr	: it is a vector which saves the sum of all mechanical loads (rr vector)
!		rrr		:
!		rnn		:
!		rnn1	:
!		dt		: R vector in this time step "N"
!		di		: it shows the difference between two time steps "DT-DTI"
!		dti		: R vector in previous time step "N-1"
!		disp	: displacement/pressures vector which consists of initial displacement value
!				  of corresponding nodes in which their identities code of Ux & Uy are defined
!				  in ID and initial water and air values of corresponding nodes in which their
!				  identities code of Pw & Pa are defined in ID. these displacement initial
!				  values can be defined for a free node or for a node with boundary conditions
!				  and the water and air initial values must be defined for a free node;
!				  in each "itime" the displacements/pressures vector (R->DT) in PREVIOUS
!				  TIME STEP (N-1) is transfered in DISP vector
!		vel		: velocity vector which consists of initial velocity value of corresponding
!				  nodes in which their identities code of Ux & Uy are defined in ID. these
!				  initial values can be defined for a free node or for a node with b.c.;
!				  in FEM the velocity vector is obtained by:
!				  V(N+1)=V(N)+[(1-anew1)*A(N)+anew1*A(N+1)]*dtime
!		acc		: acceleration vector which consists of initial acceleration value of
!				  corresponding nodes in which their identities code of Ux & Uy are defined
!				  in ID. these initial values can be defined for a free node or for a node
!				  with boundary conditions;
!				  in FEM the velocity vector is obtained by:
!				  A(N+1)=(1/anew2/dtime**2)*[U(N+1)-U(N)]-(1/anew2/dtime)*V(N)
!						 -(1/2/anew2-1)*A(N)
!		iact	:
!
!	INPUT	: ie8d,ie8c,ie8u,ie3d1,ie3d2,ie3d3,
!
!	OUTPUT	: sig8d=0,sig8c=0,sig8u=0,sig8di=0,sig8ci=0,sig8ui=0,sigmad=0,sigmac=0,
!			  sigmau=0,sig8do=0,sig8co=0,sig8uo=0,r=0,rr=0,dt=0,di=0,dti=0,disp=0,
!			  vel=0,acc=0,totrr=0,rrr=0,rnn=0,rnn1=0,iact
!
!****************************************************************************************
!
	subroutine init1(ie8d,ie8c,ie8u,sig8d,sig8c,sig8u,sig8di,sig8ci,sig8ui,sigmad,&
	                 sigmac,sigmau,sig8do,sig8co,sig8uo,sig3u,r,rr,dt,iact,ie3d1,ie3d2,&
					 ie3d3,di,dti,disp,vel,acc,totrr,rrr,rnn,rnn1)
!
	implicit double precision (a-h,o-z)
!
	common /a/ nnp,ifem,ibem,idyn,ieaq,nr
	common /a1/ ne8d,ne8c,ne8u
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /a3/ ne8d1,ne8c1,ne8u1,ne3u1,nload1,nbcx1,nbcy1,nbcw1,nbca1,m8d1,m8c1,m8u1,&
				nbed11,nbed21,nbed31,nbeint1,ne8dd,ne8cc,ne8uu
	common /g/ init,iprint
	common /h/ mdof,mdofn,ls,ls1,isolv
!
	dimension ie8d(9,ne8d1),ie8c(9,ne8c1),ie8u(9,ne8u1),&
			  sig8d(9,5,ne8d1),sig8c(9,6,ne8c1),sig8u(9,7,ne8u1),&
			  r(mdof),rr(mdof),dt(mdofn),iact(nnp),sig8di(9,5,ne8d1),sig8ci(9,6,ne8c1),&
			  sig8ui(9,7,ne8u1),sig8do(9,5,ne8d1),sig8co(9,6,ne8c1),sig8uo(9,7,ne8u1),&
			  ie3d1(npbmax,nbed11),ie3d2(npbmax,nbed21),ie3d3(npbmax),di(mdofn),&
			  dti(mdofn),disp(mdofn),vel(mdofn),acc(mdofn),totrr(mdof),rrr(mdof),&
			  sigmad(9,5,ne8d1),sigmac(9,6,ne8c1),sigmau(9,7,ne8u1),sig3u(5,500,nbed11),&
			  rnn1(mdof),rnn(mdof)
!
!
!     initilize iact
!
	do i=1,nnp
		iact(i)=0
	enddo
!
!     initialize stress arrays and find active nodes
!
	if (ne8d.eq.0) goto 110
	do i=1,ne8d
		do j=1,8
			k=ie8d(j,i)
			iact(k)=1	! it activates all of the nodes of the dry elms
		enddo
		do k=1,9
			do j=1,5
				sig8d(k,j,i)=0.0d0
				sigmad(k,j,i)=0.0d0
				sig8di(k,j,i)=0.0d0
				sig8do(k,j,i)=0.0d0
			enddo
		enddo
	enddo
!
110	if (ne8c.eq.0) goto 120
	do i=1,ne8c
		do j=1,8
			k=ie8c(j,i)
			iact(k)=1
		enddo
		do k=1,9
			do j=1,6
				sig8c(k,j,i)=0.0d0
				sigmac(k,j,i)=0.0d0
				sig8ci(k,j,i)=0.0d0
				sig8co(k,j,i)=0.0d0
			enddo
		enddo
	enddo
!
120	if (ne8u.eq.0) goto 125
	do i=1,ne8u
		do j=1,8
			k=ie8u(j,i)
			iact(k)=1
		enddo
		do k=1,9
			do j=1,7
				sig8u(k,j,i)=0.0d0
				sigmau(k,j,i)=0.0d0
				sig8ui(k,j,i)=0.0d0
				sig8uo(k,j,i)=0.0d0
			enddo
		enddo
	enddo
!
125	if (ibem.ne.3) goto 130
	do i=1,nbed11
		do k=1,5
			do j=1,50
				sig3u(k,j,i)=0.d0
			enddo
		enddo
	enddo
!
130	if (nbed1.eq.0) goto 230
	do 220 i=1,nbed1
		do 220 j=1,nnpbed1(i)
			k=ie3d1(j,i)
			iact(k)=1
220		continue
230	continue
!
	if (nbed2.eq.0) goto 330
	do 320 i=1,nbed2
		do 320 j=1,nnpbed2(i)
			k=ie3d2(j,i)
			iact(k)=1
320		continue
330	continue
!
	if (nbed3.eq.0) goto 430
	do 420 j=1,nnpbed3
		k=ie3d3(j)
		iact(k)=1
420	continue
430	continue
!
!     initialize displacement, velocity, accelerations & load arrays
!
	do i=1,mdof
		r(i)=0.0d0
		rr(i)=0.0d0
		rrr(i)=0.0d0
		totrr(i)=0.0d0
		rnn1(i)=0.d0
		rnn(i)=0.d0
	enddo
!
	do i=1,mdofn
		di(i)=0.0d0
        dt(i)=0.0d0
        dti(i)=0.0d0
        disp(i)=0.0d0
        vel(i)=0.0d0
        acc(i)=0.0d0
	enddo
!
	return
	end



!****************************************************************************************
!
!							           INIT3 subroutine
!
!	This subroutine inputs initial values of displacememts,velocities,accelerations and
!	pore pressures.
!
!	variables used are:
!
!		id		: degrees of freedom of each node from "1" to "mdof"; for the nodes with
!				  the boundary conditions, each node takes from "mdof+1" to "mdof+nbc"
!				  [X,Y]; which are free for the dry elm &
!				  [X,Y,Pw]; which are free for the saturated elm
!				  [X,Y,Pw,Pa]; which are free for the unsaturated elm
!		vel		: velocity vector which consists of initial velocity value of corresponding
!				  nodes in which their identities code of Ux & Uy are defined in ID. these
!				  initial values can be defined for a free node or for a node with b.c.;
!				  in FEM the velocity vector is obtained by:
!				  V(N+1)=V(N)+[(1-anew1)*A(N)+anew1*A(N+1)]*dtime
!		acc		: acceleration vector which consists of initial acceleration value of
!				  corresponding nodes in which their identities code of Ux & Uy are defined
!				  in ID. these initial values can be defined for a free node or for a node
!				  with boundary conditions;
!				  in FEM the velocity vector is obtained by:
!				  A(N+1)=(1/anew2/dtime**2)*[U(N+1)-U(N)]-(1/anew2/dtime)*V(N)
!						 -(1/2/anew2-1)*A(N)
!		disp	: displacement/pressures vector which consists of initial displacement value
!				  of corresponding nodes in which their identities code of Ux & Uy are defined
!				  in ID and initial water and air values of corresponding nodes in which their
!				  identities code of Pw & Pa are defined in ID. these displacement initial
!				  values can be defined for a free node or for a node with boundary conditions
!				  and the water and air initial values must be defined for a free node;
!				  in each "itime" the displacements/pressures vector (R->DT) in PREVIOUS
!				  TIME STEP (N-1) is transfered in DISP vector
!		nacci	: number of nodes on with an initial acceleration is imposed which can be
!				  at max. 400 nodes
!		nnacci	: nodes' numbers on which an initial acc. is imposed
!		vxacci	: initial acc. values imposed on the preceding nodes (X-dir)
!		vyacci	: initial acc. values imposed on the preceding nodes (Y-dir)
!		nveli	: number of nodes on with an initial velocity is imposed which can be at
!				  max. 400 nodes
!		nnveli	: nodes' numbers on which an initial vel. is imposed
!		vxveli	: initial velocity values imposed on the preceding nodes (X-dir)
!		vyveli	: initial velocity values imposed on the preceding nodes (Y-dir)
!		ndispi	: number of nodes on which an initial displacement is imposed which can be
!				  at max. 400 nodes
!		nndispi	: nodes' numbers on which an initial disp. is imposed
!		vxdispi	: initial displacement values imposed on the preceding nodes (X-dir)
!		vydispi	: initial displacement values imposed on the preceding nodes (Y-dir)
!		npwi	: number of nodes on which an initial water pressure is imposed which can
!				  be at max. 400 nodes
!		nnpwi	: nodes' numbers on which an initial water-p is imposed
!		vnpwi	: initial water pressure values imposed on the preceding nodes
!		npai	: number of nodes on which an initial gas pressure is imposed which can be
!				  at max. 400 nodes
!		nnpai	: nodes' numbers on which an initial air-p is imposed
!		vnpai	: initial air pressure values imposed on the preceding nodes
!		mdofn	: maximum degree of freedom without considering the boundary conditions
!		mdof	: maximum degree of freedom by considering the boundary conditions

!
!	INPUT	: nacci,nveli,ndispi,npwi,npai,id
!
!	OUTPUT	: disp,vel,acc
!
!****************************************************************************************
!
	subroutine init3(id,disp,vel,acc,sig3u,sm8u,ie3d1,ie3d2,ie3d3,x,dti,dt,di)
!
	implicit double precision (a-h,o-z)
!
	common /a/ nnp,ifem,ibem,idyn,ieaq,nr
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /a1/ ne8d,ne8c,ne8u
	common /a3/ ne8d1,ne8c1,ne8u1,ne3u1,nload1,nbcx1,nbcy1,nbcw1,nbca1,m8d1,m8c1,m8u1,&
				nbed11,nbed21,nbed31,nbeint1,ne8dd,ne8cc,ne8uu
	common /g/ init,iprint
	common /h/ mdof,mdofn,ls,ls1,isolv
	common /ini1/ vnpwi(400),vnpai(400),nnpwi(400),nnpai(1000)
	common /ini2/ npwi,npai,ndispi,nveli,nacci
	common /ini3/ nndispi(400),nnveli(400),nnacci(400)
	common /ini4/ vxdispi(400),vydispi(400)
	common /ini5/ vxveli(400),vyveli(400)
	common /ini6/ vxacci(400),vyacci(400)
	common /unbem/ nelmu,ie3u(3,500,5)
!
	dimension id(4,nnp),disp(mdofn),vel(mdofn),acc(mdofn),dti(mdofn),dt(mdofn),di(mdofn),&
			  sig3u(5,500,nbed11),sm8u(60,m8u1),ie3d1(npbmax,nbed11),ie3d2(npbmax,nbed21),&
	          ie3d3(npbmax),pa(3),pw(3),xelm(3),yelm(3),x(2,nnp)
!
	if (nacci.eq.0) goto 250
	do i=1,nacci
		kk=nnacci(i)
		i1=id(1,kk)		! code identity of kk concerning to Ux which can be 1<i1<mdof or
!						  mdof+1<i1<mdof+nbc???
		i2=id(2,kk)		! code identity of kk concerning to Uy which can be 1<i1<mdof or
!						  mdof+1<i1<mdof+nbc???
!		if (i1.le.mdof) acc(i1)=vxacci(i)
!		if (i2.le.mdof) acc(i2)=vyacci(i)
		acc(i1)=vxacci(i)
		acc(i2)=vyacci(i)
	enddo
250	continue
!
	if (nveli.eq.0) goto 350
	do i=1,nveli
		kk=nnveli(i)
		i1=id(1,kk)
		i2=id(2,kk)
!		if (i1.le.mdof) vel(i1)=vxveli(i)
!		if (i2.le.mdof) vel(i2)=vyveli(i)
		vel(i1)=vxveli(i)
		vel(i2)=vyveli(i)
	enddo
350	continue
!
	if (ndispi.eq.0) goto 450
	do i=1,ndispi
		kk=nndispi(i)
		i1=id(1,kk)
		i2=id(2,kk)
!		if (i1.le.mdof) disp(i1)=vxdispi(i)
!		if (i2.le.mdof) disp(i2)=vydispi(i)
		disp(i1)=vxdispi(i)
		disp(i2)=vydispi(i)
	enddo
450	continue
!
	if (npwi.eq.0) goto 550
	do i=1,npwi
		j=nnpwi(i)
		i3=id(3,j)		! code identity of kk concerning to Pw which must be 1<i1<mdof
!		if (i3.le.mdof) then	! we affect the imposed value of Pw if and only if it is
!								  a degree of freedom, i.e. id(i,j)<=mdof
			disp(i3)=vnpwi(i)
!		endif
	enddo
550	continue
!
	if (npai.eq.0) goto 650
	do i=1,npai
		j=nnpai(i)
		i4=id(4,j)		! code identity of kk concerning to Pa which must be 1<i1<mdof
!		if (i4.le.mdof) then	! we affect the imposed value of Pa if and only if it is
!								  a degree of freedom, i.e. id(i,j)<=mdof
			disp(i4)=vnpai(i)
!		endif
	enddo
650	continue
!
!1000  continue
!
!     initializing dti (according to first increment & first iteration)
!
!	do i=1,mdofn
!		di(i)=disp(i)
!		dt(i)=disp(i)
!		dti(i)=disp(i)
!	enddo
!
	if (ibem.ne.3) goto 1
	if (nbed1.eq.0) goto 2
	do k1=1,nbed1
		n=nnpbed1(k1)
		mt=ie3d1(n+1,k1)
		nelm=0
		do i=1,n-1,2
			nelm=nelm+1
			ie3u(1,nelm,k1)=ie3d1(i,k1)
			ie3u(2,nelm,k1)=ie3d1(i+1,k1)
			if (i.eq.n-1) then
				ie3u(3,nelm,k1)=ie3d1(1,k1)
				goto 660
			endif
			ie3u(3,nelm,k1)=ie3d1(i+2,k1)
660			continue
		enddo
		nelmu=nelm
!
		do m=1,nelmu
			if (npwi.eq.0.and.npai.eq.0) then
!				sig3u(5,m,k1)=0.d0
!				sati=sm8u(14,mt)
!				sig3u(4,m,k1)=sati
			else
				do i=1,3
					j=ie3u(i,m,k1)
					xelm(i)=x(1,j)
					yelm(i)=x(2,j)
					nw =id(3,j)
					na =id(4,j)
					pw(i)=disp(nw)
					pa(i)=disp(na)
				enddo
				pair=0.d0
				pwat=0.d0
				do j=1,3
					pair=pair+pa(j)
					pwat=pwat+pw(j)
				enddo
				sig3u(5,m,k1)=pair/3.
				sig3u(4,m,k1)=sig3u(5,m,k1)-pwat/3.d0
			endif
		enddo
	enddo
	goto 1
2	if (nbed2.eq.0) goto 3
3	if (nbed3.eq.0) goto 1
!
1	continue
!
	return
	end

!****************************************************************************************
!
!							           INPUT subroutine
!
!	AIM			: reads mesh,material,boundary & initial conditions,print & loading data
!	callED IN	: CONTROL
!	callS		: COORC,COORP,ELM8D,ELM8C,ELM8U,BLOCD,BLOCC,BLOCU,BEMD1,BEMD2,BEMD3
!				  BEMEN,MATED,MATEC,MATEU,BOUND,PRTOUT,PHASE,INITI
!
!	VARIABLES	:
!
!		x		: coordinates of nodes (nodal data) in FE & BE areas;
!							x(1,1:nnp) : X-coord
!							x(2,1:nnp) : Y-coord
!		ie3d1	: ie3d1(1:nnpbed1,1:nbed): nodes' number which form the current D1 BEM area
!										   Numbering direction for external surface close
!										   (finite) domain is ANTICLOCKWISE
!				  ie3d1(nnpbed1(i)+1,1:nbed): it shows the type of just one material which
!											  is used:
!							[1]: dry boundary elm
!							[2]: saturated boundary elm
!							[3]: unsaturated boundary elm
!		ie3d2	: ie3d2(1:nnpbed2,1:nbed): nodes' number which form the current D2 BEM area
!										   Numbering direction for internal surface open
!										   (infinite) domain is CLOCKWISE
!				  ie3d2(nnpbed2(i)+1,1:nbed): it shows the type of just one material which
!											  is used:
!							[1]: dry boundary elm
!							[2]: saturated boundary elm
!							[3]: unsaturated boundary elm
!		ie3d3	: ie3d3(1:nnpbed3): nodes' number which form the current D3 BEM area
!										   Numbering direction for semi-infinite domain is
!										   ANTICLOCKWISE
!				  ie3d3(nnpbed3+1): it shows the type of just one material which
!											  is used:
!							[1]: dry boundary elm
!							[2]: saturated boundary elm
!							[3]: unsaturated boundary elm
!		ie3d3int	: ie3d3int(1:nnpbed3int): internal nodes' number in D3 area which the
!											  numbering direction is ANTICLOCKWISE
!		xencl: X-coordinates of enclosing elements' nodes (nodal data)
!		yencl: Y-coordinates of enclosing elements' nodes (nodal data)
!		sm8d : all of the parameters corresponding to each behaviour of dry soil are
!			   defined in this matrix.
!						  icpt:[1]-> sm8d (3,i) = E
!									 sm8d (7,i) = K
!									 sm8d (10,i)= Ko
!									 sm8d (11,i)= gammas
!						  icpt:[2]-> sm8d (1,i) = friction angle (�)
!									 sm8d (2,i) = cohesion (Pa)
!									 sm8d (3,i) = E-load, loading modulus used in non-linear
!												  elasticity (adimensionel)
!									 sm8d (4,i) = E-unld, unloading modulus used in
!												  non-linear elasticity (adim)
!									 sm8d (5,i) = n, exponent used to compute the Young
!												  loading and unloading moduli in non-linear
!											      elasticity (adim)
!									 sm8d (6,i) = Rf, coefficient multiplying the ratio
!												  (sigm1-sigm3)/(sigm1-sigm3)ult
!												  in nlin els (adim)
!                                    sm8d (7,i) = Kb, bulk modulus used in non-linear
!											      elasticity (adim)
!	                                 sm8d (8,i) = m, exponent used to compute the bulk
!										          modulus B in nlin els (adim)
!									 sm8d (9,i) = e-f, Young modulus minimal value (Pa)
!									 sm8d (10,i)=Ko, steady state coefficient (adim)
!								     sm8d (11,i)= gammas, soil volumetric weight (N,m-3)
!								     sm8d (12,i)= ten-st, traction resistency (Pa)
!						  icpt:[3]-> sm8d (1,i) = friction angle
!									 sm8d (2,i) = cohesion
!								     sm8d (3,i) = OCR
!								     sm8d (4,i) = k
!                                    sm8d (5,i) = n, exponent used to compute the Young
!												  loading and unloading moduli in non-linear
!												  elasticity (adim)
!									 sm8d (6,i) = e
!									 sm8d (7,i) = Kb, bulk modulus used in non-linear
!												  elasticity (adim)
!									 sm8d (8,i) = m, exponent used to compute the bulk
!												  modulus B in nlin els (adim)
!									 sm8d (10,i) =Ko
!									 sm8d (11,i) =gamma
!						  icpt:[4]-> ??????
!		sm8c	: all of the parameters corresponding to each behaviour of saturated soil
!				  are defined in this matrix.
!						  icpt:[1]-> sm8c (3,i) = E
!									 sm8c (7,i) = K
!									 sm8c (10,i)= Ko
!									 sm8c (11,i)= gamma
!									 sm8c (13,i)= kz, vertical water permeability (m.s^-1)
!									 sm8c (14,i) = alpha, 1-Kd/Ks
!								     sm8c (15,i) = compm, Q
!									 sm8c (16,i) = kx/kz, ratio of the horizontal water
!												   permeability to the vertical water
!												   permeability.
!						  icpt:[2]-> sm8c (1,i) = friction angle (�)
!									 sm8c (2,i) = cohesion (Pa)
!									 sm8c (3,i) = E-load, loading modulus used in non-linear
!												  elasticity (adimensionel)
!									 sm8c (4,i) = E-unld, unloading modulus used in
!												  non-linear elasticity (adim)
!									 sm8c (5,i) = n, exponent used to compute the Young
!												  loading and unloading moduli in non-linear
!												  elasticity (adim)
!									 sm8c (6,i) = Rf, coefficient multiplying the ratio
!												  (sigm1-sigm3)/(sigm1-sigm3)ult in nlin els
!												  (adim)
!                                    sm8c (7,i) = Kb, bulk modulus used in non-linear
!											      elasticity (adim)
!	                                 sm8c (8,i) = m, exponent used to compute the bulk
!												  modulus B in nlin els (adim)
!									 sm8c (9,i) = e-f, Young modulus minimal value (Pa)
!									 sm8c (10,i)=Ko, steady state coefficient (adim)
!								     sm8c (11,i)= gamma, soil volumetric weight (N,m-3)
!								     sm8c (12,i)= ten-st, traction resistency (Pa)
!									 sm8c (13,i)= perm-z
!									 sm8c (14,i)= sat
!									 sm8c (15,i)= compm
!									 sm8c (16,i)= kx/kz
!						  icpt:[3]-> sm8c (1,i) = friction angle
!									 sm8c (2,i) = cohesion
!								     sm8c (3,i) = OCR
!								     sm8c (4,i) = k
!                                    sm8c (5,i) = n, exponent used to compute the Young
!												  loading and unloading moduli in non-linear
!												  elasticity (adim)
!									 sm8c (6,i) = e
!									 sm8c (7,i) = Kb, bulk modulus used in non-linear
!												  elasticity (adim)
!									 sm8c (8,i) = m, exponent used to compute the bulk
!												  modulus B in nlin els (adim)
!									 sm8c (10,i) =Ko
!									 sm8d (11,i) =gamma
!									 sm8c (13,i) = perm-z
!									 sm8c (14,i) = sat
!									 sm8c (15,i) = compm
!									 sm8c (16,i) = kx/kz
!						  icpt:[4]-> ??????
!		sm8u	: all of the parameters corresponding to each behaviour of unsat. soil are
!				  defined in this matrix.
!						  icpt:[2]->
!									sm8u(1,i)= friction angle (�)
!									sm8u(2,i)= cohesion (Pa)
!									sm8u(3,i)= Kl,loading modulus used in non-linear
!											   elasticity(adim.)
!									sm8u(4,i)= Ku,unloading modulus used in non-linear
!											   elasticity (adim.)
!									sm8u(5,i)= n,exponent used to compute the Young loading
!											   and unloading moduli in non-linear elasticity
!									sm8u(6,i)= Rf,coefficient multiplying the ratio
!											   (SIGM1-SIGM3)/(SIGM1-SIGM3)ult in non-linear
!											   elasticity
!									sm8u(7,i)= Kb,bulk modulus used in non-linear elasticity
!									sm8u(8,i)= m,exponent used to compute the bulk modulus B
!											   in non-linear elasticity (adimensional)
!									sm8u(9,i)= k0, steady state coefficient (adimensional)
!									sm8u(10,i)= traction resistance (Pa ?)
!									sm8u(11,i)= Young modulus minimal value (Pa)
!									sm8u(12,i)= OCR, over-consolidation ratio (adimensional)
!									sm8u(13,i)= -
!									sm8u(14,i)= Sw0,initial saturation degree(adimensional)
!									sm8u(15,i)= a, constant (m.s-1) used in the formula of
!						        				water permeability
!									sm8u(16,i)= alpha_w, constant (adimensional) used in
!												the formula of water permeability
!									sm8u(17,i)= Swr,residual saturation degree (adim.),
!												used in the formula of water permeability
!									sm8u(19,i)= b,constant (m2) used in the formula of air
!											    permeability
!									sm8u(20,i)= a,constant (adim.) used in the formula of air
!												permeability
!									sm8u(21,i)= "viscosity",dynamic air viscosity(N.s.m-2),
!												used in the formula of air permeability.
!												Generally,Visc_a=1,846.10-5 N.s.m-2.
!									sm8u(22,i)= Kw_max, maximal authorized water permeability
!												value(m.s-1)
!									sm8u(23,i)= Ka_max, maximal authorized air permeability
!												value (m.s-1)
!									sm8u(24,i)= ae, constant (adimensional) used in the
!												formula of the void ratio state surface
!									sm8u(25,i)= be, constant (adimensional) used in the
!												formula of the void ratio state surface
!									sm8u(27,i)= de, constant (adimensional) used in the
!												formula of the void ratio state surface
!									sm8u(28,i)= SIGMe,maximal traction resistance (Pa)
!												used in the formula of the void ratio
!												state surface
!									sm8u(29,i)= as,constant (adimensional) used in the
!												formula of the saturation degree state
!												surface
!									sm8u (30,i)= bs,constant (adimensional) used in the
!												 formula of the saturation degree state
!												 surface
!									sm8u (33,i)= msuc-e
!									sm8u (34,i)= msuc-c
!
!		kodex	: numbers of nodes in which zero displacement in X-dir is imposed
!		kodey	: numbers of nodes in which zero displacement in Y-dir is imposed
!		kodew	: numbers of nodes in which zero pore water pressure is imposed
!		kodea	: numbers of nodes in which zero pore gas pressure is imposed
!		ie8dout	: numbers of the nodes forming the current drained element that we need
!				  their outputs, IE8doUT (1:NE8do)
!		ie8cout	: numbers of the nodes forming the current saturated element that we need
!				  their outputs, IE8COUT (1:NE8CO)
!		ie8uout	: numbers of the nodes forming the current unsaturated element that we need
!				  their outputs, IE8UOUT (1:NE8UO)
!		iconst	: iconst(1,nload1)= Initial load step number (time step)
!				  iconst(2,nload1)= Number of load steps (nstep)
!				  SUM[iconst(2,1:nload)]=ntime
!				  iconst(3,nload1)= Maximum iterations per load step (itmax)
!
!		ie8d	:
!		ie8c	:
!		ie8u	:
!
!
!	INPUT	:
!
!	OUTPUT	:
!
!****************************************************************************************
!
	subroutine input(ie8d,ie8c,ie8u,iconst,kodex,kodey,kodew,kodea,x,xint,sm8d,sm8c,sm8u,&
					 ie3d1,ie3d2,ie3d3,infile)
!
	implicit double precision (a-h,o-z)
	character (5) :: infile*72, card, wd(21)
	logical :: 	strcomp,nmiss
!
	common /ipt/ istop,list,wd
	common /a/ nnp,ifem,ibem,idyn,ieaq,nr
	common /a1/ ne8d,ne8c,ne8u
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /a3/ ne8d1,ne8c1,ne8u1,ne3u1,nload1,nbcx1,nbcy1,nbcw1,nbca1,m8d1,m8c1,m8u1,&
				nbed11,nbed21,nbed31,nbeint1,ne8dd,ne8cc,ne8uu
	common /a4/ xencl(500),yencl(500)
	common /b/ m8d,m8c,m8u,m8dep(10),m8cep(10),m8uep(10)
	common /c/ nload,ntime,itmax,itime,iter,iload,time,dtime,dtimei,dtimem,xinc
	common /d/ nbcx,nbcy,nbcw,nbca,nbc
	common /e/ gammaw,gammas,gammaa,grav,atmp
	common /g/ init,iprint
	common /h/ mdof,mdofn,ls,ls1,isolv
	common /h1/ ptoll1,ptoll3
	common /i/iterp,itprt(100),ne8do,ne8co,ne8uo,ie8dout(10),ie8cout(10),ie8uout(10),nnpo,&
			  inpout(100)
	common /ini1/ vnpwi(400),vnpai(400),nnpwi(400),nnpai(400)
	common /ini2/ npwi,npai,ndispi,nveli,nacci
	common /ini3/ nndispi(400),nnveli(400),nnacci(400)
	common /ini4/ vxdispi(400),vydispi(400)
	common /ini5/ vxveli(400),vyveli(400)
	common /ini6/ vxacci(400),vyacci(400)
	common /Kcor/ kcorf,mtime,rmtol,c1,c2
	common /ocqm/ EPScqm,taw,Icqm,Lcqm,NPcqm
!
	dimension ie8d(9,ne8d1),ie8c(9,ne8c1),ie8u(9,ne8u1),iconst(3,nload1),kodex(nbcx1),&
			  kodey(nbcy1),kodew(nbcw1),kodea(nbca1),x(2,nnp),xint(2,nbeint1),&
			  sm8d(60,m8d1),sm8c(60,m8c1),sm8u(60,m8u1),ie3d1(npbmax,nbed11),&
			  ie3d2(npbmax,nbed21),ie3d3(npbmax)
!
!
	list=21
	istop=0
	nmiss=.false.
 	wd(1)='coorc'		! the coordinates of nodes in cartesian system
	wd(2)='coorp'		! the coordinates of nodes in polar system
 	wd(3)='coint'		! the coordinates of nodes in cartesian system
	wd(4)='elm8d'		!
	wd(5)='elm8c'		!
	wd(6)='elm8u'		!
    wd(7)='blocd'		!
	wd(8)='blocc'		!
	wd(9)='blocu'		!
	wd(10)='bemd1'		!
	wd(11)='bemd2'		!
    wd(12)='bemd3'		!
	wd(13)='bemen'		!
	wd(14)='mated'		!
	wd(15)='matec'		!
	wd(16)='mateu'		!
    wd(17)='bound'		!
	wd(18)='print'		!
	wd(19)='phase'		!
	wd(20)='initi'		!
    wd(21)='force'		!
!
!
!-------------------------------
!	    to control input
!-------------------------------
!
10	read (4,*) card
	if (card(1:1).eq.'#') goto 10
	do i=1,list
		if (strcomp(card,wd(i))) goto 15
	enddo
	write (*,*) 'KEY IS not CORRECT'
	istop = istop + 1
	goto 10
!
15	if (i.eq.1)  call coorc(nnp,x,nmiss)
	if (i.eq.2)  call coorp(nnp,x,nmiss)
	if (i.eq.3)  call coint(nnp,xint,nmiss)
	if (i.eq.4)  call elm8d(ne8d,ie8d)
	if (i.eq.5)  call elm8c(ne8c,ie8c)
!	if (i.eq.6)  call elm8u(ne8u,ie8u)
	if (i.eq.7)  call blocd(nnp,x,ne8d,ie8d,nmiss)
	if (i.eq.8)  call blocc(nnp,x,ne8c,ie8c,nmiss)
!	if (i.eq.9)  call blocu(nnp,x,ne8u,ie8u,nmiss)
	if (i.eq.10) call bemd1(ibem,nnp,nbed1,nnpbed1,ie3d1,npbmax)
	if (i.eq.11) call bemd2(ibem,nnp,nbed2,nnpbed2,ie3d2)
	if (i.eq.12) call bemd3(ibem,nnp,nbed3,nnpbed3,ie3d3)
	if (i.eq.13) call bemen(nnpen,xencl,yencl)
	if (i.eq.14) call mated(m8d,sm8d,m8dep,ibem,ne8d)
	if (i.eq.15) call matec(m8c,sm8c,m8cep,ibem,ne8c)
	if (i.eq.16) call mateu(m8u,sm8u,m8uep,ibem,ne8u,idyn)
	if (i.eq.17) call bound(nnp,nbcx,nbcy,nbcw,nbca,kodex,kodey,kodew,kodea)
	if (i.eq.18) call prtout
	if (i.eq.19) call phase(nload,iconst)
	if (i.eq.20) call initi(idyn,ifem,ibem)
	if (i.eq.21) goto 20
	goto 10
!
20 continue
!
!
!---------------------------------------
!     to verify the number of errors
!---------------------------------------
!
!	the number of registered nodes
	n=0
	do i=1,nnp
		if (x(1,i).eq.0.d0 .and. x(2,i).eq.0.d0) n=n+1
	enddo
	if (nmiss) n=n-1
	if (n.gt.0) then
		istop=istop+1
        write (*,51) n
	endif
!
!	number of registered dry elements
	n=0
	do i=1,ne8d
		if (ie8d(1,i).eq.0) n=n+1
	enddo
	if (n.gt.0) then
		istop=istop+1
        write (*,52) n
	endif
!
!	number of registered saturated elements
	n=0
	do i=1,ne8c
		if (ie8c(1,i).eq.0) n=n+1
	enddo
	if (n.gt.0) then
		istop=istop+1
        write (*,53) n
	endif
!
	write (*,54) istop
	if (istop.eq.0) goto 25
	stop
!
!
!     --------------------------------
!                write data
!     --------------------------------
!
!     basic data
25	if (iprint.ne.1) goto 900
	write (7,150)
	write (7,155) nnp,ifem,ibem
!
	if (ifem.eq.0) goto 30
	write (7,160) ne8d,ne8c,ne8u
!
30	if (ibem.eq.0) goto 35
	write (7,165) nbed1,nbed2,nbed3,nbeint
	if (nbed1.gt.0) then
		write (7,170)
		write (7,* ) (nnpbed1(i),i=1,nbed1)
	endif
	if (nbed2.gt.0) then
		write (7,175)
		write (7,*) (nnpbed2(i),i=1,nbed2)
	endif
	if (nbed3.gt.0) write (7,180) nnpbed3,nnpen
!
35	write (7,185) m8d,m8c,m8u
	write (7,190)
	if (m8d.eq.0) goto 40
	write (7,195)
	write (7,200)
	write (7,205) (i,m8dep(i),i=1,m8d)
40	if	(m8c.eq.0) goto 45
	write (7,210)
	write (7,200)
	write (7,205) (i,m8cep(i),i=1,m8c)
45	if	(m8u.eq.0) goto 50
	write (7,215)
	write (7,200)
	write (7,205) (i,m8uep(i),i=1,m8u)
!
50	write (7,220) idyn,ieaq,nload,ntime,isolv,nr,dtimei,dtimem,ptoll1,ptoll3
!
	write (7,225) kcorf,mtime,rmtol
	write (7,230) gammaw,gammas,gammaa,grav
!
!
!     nodal coordinates
	write (7,235)
	do i=1,nnp
		write (7,240) i,x(1,i),x(2,i)
	enddo
!
!
!	  drained elements in FE zone
	if (ne8d.eq.0) goto 55
	write (7,245)
	do i=1,ne8d
		write (7,250) i,(ie8d(l,i),l=1,9)
	enddo
!
!
!	  saturated elements in FE zone
55	if (ne8c.eq.0) goto 60
	write (7,255)
	do i=1,ne8c
		write (7,250) i,(ie8c(l,i),l=1,9)
	enddo
!
!
!	  unsaturated elements in FE zone
60	if (ne8u.eq.0) goto 65
	write (7,260)
	do i=1,ne8u
		write (7,250) i,(ie8u(l,i),l=1,9)
	enddo
!
!
!     boundary finite zones d1
65	if (nbed1.eq.0) goto 70
	write (7,265)
	do i=1,nbed1
		write (7,270) i
		write (7,250) (ie3d1(j,i),j=1,nnpbed1(i))
	enddo
!
!
!     boundary infinite zones d2
70	if (nbed2.eq.0) goto 75
	write (7,275)
	do i=1,nbed2
		write (7,270) i
		write (7,250) (ie3d2(j,i),j=1,nnpbed2(i))
	enddo
!
!
!	  boundary semi-infinite zone d3 & 'enclosing element'
75	if (nbed3.eq.0) goto 76
	write (7,280)
	write (7,250) (ie3d3(j),j=1,nnpbed3)
!
!
76	if (nbeint.eq.0) goto 80
	write (7,281)
	do i=1,nbeint
		write (7,240) i,xint(1,i),xint(2,i)
	enddo
!
!	  it gives the BE nodal carachteristics in another file 'noded3.txt'
	if (ieaq.eq.3) then		! it is valid just for ieaq=3 (incident wave)?
!							  for the other code like 1 &2, how is it?????
		open (1,file='noded3.txt')
		do j=1,nnpbed3
			write (1,*) x(1,ie3d3(j)),x(2,ie3d3(j))
		enddo
		write (1,281)
		do j=1,nbeint
			write (1,*) xint(1,j),xint(2,j)
		enddo
		close (1)
	endif
!	  it gives the data of enclosing element
	write (7,285)
	do i=1,nnpen
		write (7,240) i,xencl(i),yencl(i)
	enddo
!
!
!	  mesh file & result readen by GID7.1
	call maillage(x,ie8d,ie8c,ie8u,ie3d1,ie3d2,ie3d3,infile)
!
!     dry materials
80	if (m8d.eq.0) goto 82
	write (7,290)
	do i=1,m8d
		if (m8dep(i).eq.1) then
			write (7,295)
			write (7,300) sm8d(3,i),sm8d(7,i),sm8d(10,i),sm8d(11,i)
		endif
		if (m8dep(i).eq.2) then
			write (7,305)
			write (7,310) (sm8d(k,i),k=1,8)
			write (7,315)
			write (7,310) (sm8d(k,i),k=9,12)
		endif
		if (m8dep(i).eq.3) then
			write (7,320)
			write (7,310) (sm8d(k,i),k=1,8)
            write (7,325)
            write (7,310) sm8d(10,i),sm8d(11,i)
		endif
		if (m8dep(i).eq.4) then
            write (7,330)
            write (7,310) (sm8d(k,i),k=1,7)
            write (7,335)
            write (7,310) (sm8d(k,i),k=9,12)
            write (7,340)
			do nn=1,sm8d(9,i)
				j=4*nn
				write (7,310) sm8d(17+j,i),sm8d(18+j,i),sm8d(19+j,i),sm8d(20+j,i)
			enddo
		endif
	enddo
!
!
!	saturated materials
82	if (m8c.eq.0) goto 83
	write (7,480)
	do i=1,m8c
		if (m8cep(i).eq.1) then
			write (7,485)
			write (7,310) sm8c(3,i),sm8c(7,i),sm8c(10,i),sm8c(11,i),sm8c(13,i),sm8c(14,i),&
						  sm8c(15,i),sm8c(16,i)
		endif
		if (m8cep(i).eq.2) then
			write (7,305)
			write (7,310) (sm8c(k,i),k=1,8)
			write (7,490)
			write (7,310) (sm8c(k,i),k=9,16)
		endif
		if (m8cep(i).eq.3) then
			write (7,320)
			write (7,310) (sm8c(k,i),k=1,8)
			write (7,495)
			write (7,310) sm8c(10,i),sm8c(11,i),(sm8c(k,i),k=13,16)
		endif
		if (m8cep(i).eq.4) then
			write (7,330)
			write (7,310) (sm8c(k,i),k=1,7)
			write (7,500)
			write (7,310) (sm8c(k,i),k=9,16)
			write (7,340)
			do nn=1,sm8c(9,i)
				j=4*nn
				write (7,310) sm8c(17+j,i),sm8c(18+j,i),sm8c(19+j,i),sm8c(20+j,i)
			enddo
		endif
	enddo
!
!
!	unsaturated materials
83	if (m8u.eq.0) goto 85
	write (7,505)
	if (idyn.eq.2) then
		do i=1,m8u
			if (m8uep(i).eq.1) then
				write (7,296)
				write (7,310) sm8u(3,i),sm8u(5,i),sm8u(6,i),sm8u(7,i)
				write (7,511)
				write (7,515) sm8u(11,i),sm8u(14,i)
				write (7,520)
				write (7,516) sm8u(15,i),sm8u(16,i),sm8u(17,i),sm8u(22,i)
				write (7,521)
				write (7,516) sm8u(19,i),sm8u(20,i),sm8u(21,i),sm8u(23,i)
				write (7,525)
				write (7,516) sm8u(24,i),sm8u(25,i),sm8u(27,i),sm8u(28,i)
				write (7,530)
				write (7,517) sm8u(31,i)
				write (7,535)
				write (7,517) sm8u(33,i)
			endif
		enddo
	else
		do i=1,m8u
			if (m8uep(i).eq.1) then
				write (7,2960)
				write (7,310) sm8u(3,i),sm8u(7,i)
				write (7,511)
				write (7,515) sm8u(11,i),sm8u(14,i)
				write (7,520)
				write (7,516) sm8u(15,i),sm8u(16,i),sm8u(17,i),sm8u(22,i)
				write (7,521)
				write (7,516) sm8u(19,i),sm8u(20,i),sm8u(21,i),sm8u(23,i)
				write (7,525)
				write (7,516) sm8u(24,i),sm8u(25,i),sm8u(27,i),sm8u(28,i)
				write (7,530)
				write (7,517) sm8u(31,i)
				write (7,535)
				write (7,517) sm8u(33,i)
				write (7,5350)
				write (7,517) sm8u(34,i)
			endif
		enddo
	endif
!
!
!     fixed boundary conditions
85	write (7,345) nbcx,nbcy,nbcw,nbca
	if (nbcx.eq.0) goto 90
	write (7,350)
	write (7,355) (kodex(i),i=1,nbcx)
90  if (nbcy.eq.0) goto 95
	write (7,360)
	write (7,355) (kodey(i),i=1,nbcy)
95  if (nbcw.eq.0) goto 100
    write (7,365)
    write (7,355) (kodew(i),i=1,nbcw)
100 if (nbca.eq.0) goto 105
    write (7,370)
    write (7,355) (kodea(i),i=1,nbca)
!
!
!     loading phases
105	write (7,375)
	do i=1,nload
		write (7,380) (iconst(j,i),j=1,3)
	enddo
!
!
!     initial values
	write (7,385) npwi,npai,ndispi,nveli,nacci
	if (nacci.eq.0) goto 110
	write (7,390)
	write (7,395) (nnacci(i),i=1,nacci)
	write (7,400)
	write (7,405) (vxacci(i),i=1,nacci)
	write (7,410)
	write (7,405) (vyacci(i),i=1,nacci)
!
110	if (nveli.eq.0) goto 115
	write (7,415)
	write (7,395) (nnveli(i),i=1,nveli)
	write (7,420)
	write (7,405) (vxveli(i),i=1,nveli)
	write (7,425)
	write (7,405) (vyveli(i),i=1,nveli)
!
115	if (ndispi.eq.0) goto 120
	write (7,430)
	write (7,395) (nndispi(i),i=1,ndispi)
	write (7,435)
	write (7,405) (vxdispi(i),i=1,ndispi)
	write (7,440)
	write (7,405) (vydispi(i),i=1,ndispi)
!
120	if (npwi.eq.0) goto 125
	write (7,445)
	write (7,395) (nnpwi(i),i=1,npwi)
	write (7,450)
	write (7,405) (vnpwi(i),i=1,npwi)
!
125	if (npai.eq.0) goto 130
	write (7,455)
	write (7,395) (nnpwi(i),i=1,npwi)
	write (7,460)
	write (7,405) (vnpwi(i),i=1,npwi)
!
!
!     result informations to print
130	write (7,465) init, iprint,iterp
	write (7,470)
	write (7,475) (itprt(i),i=1,iterp)
!
!
!     --------------------------------
!                FORMATS
!     --------------------------------
!
51	FORMAT (5x,'**Error**',i4,' missing nodes ')
52	FORMAT (5x,'**Error**',i4,' missing drained elements ')
53	FORMAT (5x,'**Error**',i4,' missing consolidation elements ')
54	FORMAT (//10x,'Total Number Of Data Errors = ',i5)
!
5	FORMAT (6i5)
8	FORMAT (10i5)
23	FORMAT (3i5)
27	FORMAT (//10x,'Node Numbers & Given Values For Drain Interfaces'//)
28	FORMAT (6(2x,e10.3))
!
!
150	FORMAT (//10x,'basic data'//)
155	FORMAT (5x,'number of total nodes             = ',i5/&
		    5x,'FEM code                          = ',i5/&
            5x,'BEM code                          = ',i5//)
!
160	FORMAT (5x,'number of drained elements        = ',i5/&
            5x,'number of consolidated elements   = ',i5/&
            5x,'number of unsaturated elements	  = ',i5//)

!
165	FORMAT (5x,'number of D1 BEM areas            = ',i5/&
            5x,'number of D2 BEM areas            = ',i5/&
			5x,'number of D3 BEM areas            = ',i5/&
            5x,'number of internal nodes in BEM ar= ',i5//)
170	FORMAT (/5x,'number of node of D1 BEM areas   = ',/5x)
175	FORMAT (/5x,'number of node of D2 BEM areas   = ',/5x)
180	FORMAT (5x,'number of nodes of D3 BEM area    = ',i5/&
            5x,'number of nodes of en.el. area    = ',i5//)
!
185	FORMAT (//10x,'material data'//&
             5x,'number of drained materials       = ',i5/&
             5x,'number of consolidated materials  = ',i5/&
             5x,'number of unsaturated materials   = ',i5/)
!
190	FORMAT (10x,'behavior type:      '/)
195	FORMAT (10x,'drained elements')
200	FORMAT (5x,'material',3x,'type')
205	FORMAT (4x,i5,4x,i5/)
210	FORMAT (10x,'consolidated elements')
215	FORMAT (10x,'unsaturated elements')
!
220	FORMAT (//10x,'loading and construction data'//&
             5x,'dynamic code                      = ',i5/&
             5x,'earthquake code                   = ',i5/&
             5x,'number of load steps              = ',i5/&
             5x,'number of time steps              = ',i5/&
             5x,'solution code                     = ',i5/&
             5x,'modified newton raphson code      = ',i5/&
             5x,'minimum time step                 = ',e10.3/&
             5x,'maximum time step                 = ',e10.3/&
             5x,'integration constant-beta         = ',e10.3/&
             5x,'permissible tollerance            = ',e10.3//)
!
225 FORMAT (//10x,'K-Correction parameters'//&
     		 5x,'K-Correction code                 = ',i5/&
     	     5x,'Cutting Point of near history     = ',i5/&
     	     5x,'K-Correction Tolerance            = ',e10.3//)
!
230 FORMAT (//10x,'physical constants'//&
             5x,'unit weight of water              = ',e10.3/&
             5x,'unit weight of solid              = ',e10.3/&
             5x,'unit weight of AIR                = ',e10.3/&
             5x,'acc. due to gravity               = ',e10.3//)
!
235	FORMAT (//10x,'MESH and MATERIAL data'//&
           10x,'Nodal Data'//&
           1x,'No.',13x,'X-coord',11x,'Y-coord'/)
240	FORMAT (i5,4x,e15.4,3x,e15.4)
!
245	FORMAT (//10x,'Drained Element Data'//&
		   '  No.',4x,'i',4x,'j',4x,'k',4x,'l',&
		   4x,'m',4x,'n',4x,'o',4x,'p',2x,'mtype'/)
250	FORMAT (10i5)
!
255	FORMAT (//10x,'Saturated Element Data'//&
           '  No.',4x,'i',4x,'j',4x,'k',4x,'l',&
           4x,'m',4x,'n',4x,'o',4x,'p',2x,'mtype'/)
!
260	FORMAT (//10x,'Unsaturated Element Data'//&
           '  No.',4x,'i',4x,'j',4x,'k',4x,'l',&
           4x,'m',4x,'n',4x,'o',4x,'p',2x,'mtype'/)
265	FORMAT (//10x,'Nodes number of BEM D1 finite domains'/)
270	FORMAT (/5x,'Domain',i5/)
275	FORMAT (//10x,'Nodes number of BEM D2 infinite domains'/)
280	FORMAT (//10x,'Nodes number of BEM D3 domain'/)
281	FORMAT (//10x,'Internal Nodes number in BEM domain'/)
285	FORMAT (//10x,'Data Of Enclosing Element'//&
              1x,'No.',13x,'X-coord',11x,'Y-coord'/)
290	FORMAT (//5x,'Drained material prepeties')
295	FORMAT (//10x,'Linear Elastic Material'//&
			  3x,'E',8x,'K',8x,'Ko',8x,'gamma'/)
296	FORMAT (//10x,'Linear Elastic Material'//&
			  3x,'Kl',8x,'Ka',8x,'Kw',8x,'Kb'/)
2960 FORMAT (//10x,'Linear Elastic Material'//&
			  3x,'Kl',8x,'Kb'/)
300	FORMAT (4e10.3)
305 FORMAT (//10x,'Non Linear Elastic Hyperbolic Material for Statics'&
			//5x,'phi',4x,'cohesion',3x,'E-load',4x,'E-unld',7x,'n',8x,'Rf',&
			  9x,'K-b',9x,'m'/)
310	FORMAT (8e10.3)
315	FORMAT (//5x,'e-f',8x,'Ko',6x,'gamma',5x,'ten-st'//)
320	FORMAT (//10x,'Non Linear Elastic Hyperbolic Material for Dynamics'&
			//5x,'phi',4x,'cohesion',3x,'OCR',4x,'k',7x,'n',8x,'e',&
			  9x,'K-b',9x,'m'/)
325	FORMAT (//5x,'Ko',6x,'gamma'//)
330	FORMAT (//10x,'Elastoplastic Material'//&
			  5x,'phi',4x,'cohesion',3x,'g1/pref',3x,'b1/pref',5x,'pow',6x,&
			  'etac',6x,'etae'/)
335	FORMAT (//5x,'nys',8x,'ko',6x,'gamma',5x,'ten-st'//)
340	FORMAT (//5x,'alfa0',8x,'hc1/pref',6x,'he1/pref',5x,'radius'//)
!
345	FORMAT (//10x,'boundary condition data'//&
			  5x,'number of b.c. in x-dir           = ',i5/&
			  5x,'number of b.c. in y-dir           = ',i5/&
			  5x,'number of b.c. in w-dir           = ',i5/&
			  5x,'number of b.c. in a-dir           = ',i5//)
350	FORMAT (//10x,'B.C. Data In X-dir'//&
			  5x,'Node Numbers With Zero Displacement'/)
355	FORMAT (10i5)
360	FORMAT (//10x,'B.C. Data In Y-dir'//&
			  5x,'Node Numbers With Zero Displacements'/)
365	FORMAT (//10x,'B.C. Data In Pressure'//&
			  5x,'Node Numbers With Zero Pore Water Pressure'//)
370	FORMAT (//10x,'B.C. Data In Pressure'//&
			  5x,'Node Numbers With Zero Pore Air Pressure'//)
!
375	FORMAT (//10x,'Loading Data'//&
			  5x,'time step',5x,'nstep',5x,'iteration'//)
380	FORMAT (5x,3(i5,5x))
!
385	FORMAT (//10x,'number of nodes with initials values'//&
			  5x,'number of nodes with initial p-water values       = ',i5/&
			  5x,'number of nodes with initial p-air values		    = ',i5/&
			  5x,'number of nodes with initial displacement values  = ',i5/&
			  5x,'number of nodes with initial velocity values      = ',i5/&
			  5x,'number of nodes with initial acceleration values  = ',&
			  i5/)
390	FORMAT (//10x,'node numbers of nodes with initial',&
				  ' accelerations'//)
395	FORMAT (8i5)
400	FORMAT (//10x,'initial accelerations for the above nodes',&
				  ' (x dir)'//)
405	FORMAT (6(2x,e10.3))
410	FORMAT (//10x,'initial accelerations for the above nodes',&
				  ' (y dir)'//)
415	FORMAT (//10x,'node numbers of nodes with initial velocities'//)
420	FORMAT (//10x,'initial velocities for the above nodes (x dir)'//)
425	FORMAT (//10x,'initial velocities for the above nodes (y dir)'//)
430	FORMAT (//10x,'node numbers of nodes with initial',&
				  ' displacements '//)
435	FORMAT (//10x,'initial displacements',&
				  ' for the above nodes (x dir)'//)
440	FORMAT (//10x,'initial displacements',&
				  ' for the above nodes (y dir)'//)
445	FORMAT (//10x,'node numbers of nodes with initial p-water '//)
450	FORMAT (//10x,'initial p-water values for the above nodes'//)
455	FORMAT (//10x,'node numbers of nodes with initial a-water '//)
460	FORMAT (//10x,'initial a-water values for the above nodes'//)
465	FORMAT (//10x,'initial condition data'//&
			  5x,'initial stress code for soil      = ',i5/&
			  5x,'initial data output print code    = ',i5/&
			  5x,'number of printouts               = ',i5//)
!
470	FORMAT (//10x,'time steps for printouts'//)
475	FORMAT (5x,8i5/)
!
480	FORMAT (//5x,'Saturated material prepeties')
485	FORMAT (//10x,'Linear Elastic Material'//&
			  3x,'E',8x,'K',8x,'Ko',8x,'gamma',&
			  5x,'perm-z',5x,'sat',5x,'compm',5x,'kx/kz'/)
490 FORMAT (//5x,'e-f',8x,'Ko',6x,'gamma',5x,'ten-st',5x,&
			  'perm-z',5x,'sat',5x,'compm',5x,'kx/kz'/)
495	FORMAT (//5x,'Ko',6x,'gamma',5x,&
			  'perm-z',5x,'sat',5x,'compm',5x,'kx/kz'/)
500	FORMAT (//5x,'nys',8x,'ko',6x,'gamma',5x,'ten-st',5x,&
			'perm-z',5x,'sat',5x,'compm',5x,'kx/kz'/)
!
505	FORMAT (//5x,'Unsaturated material prepeties')
510 FORMAT (//5x,'Ko',6x,'ten-st',5x,'e-f',8x,'OCR',5x,&
			  'ini-sat'/)
511 FORMAT (//5x,'Emin[Pa]',6x,'Sw0'/)
515	FORMAT (2(4x,e10.3))
516	FORMAT (4(4x,e10.3))
517	FORMAT (8x,e10.3)
520	FORMAT (//5x,'Water flow data'//&
			  8x,'aw',10x,'alphaw',10x,'sat-ru',10x,'perwmax'//)
521	FORMAT (//5x,'Air flow data'//&
			  8x,'b',10x,'alphaa',10x,'viscosity',10x,'peramax'//)
525	FORMAT (//5x,'state surfaces data'//&
			  8x,'a-e',10x,'b-e',10x,'e0',10x,'sig-e'//)
530	FORMAT (//5x,'Sr data'//&
			  7x,'beta-w'//)
535	FORMAT (//5x,'suction data'//&
			  3x,'msuc-e'//)
5350 FORMAT (//5x,'Henry Coefficient'//)

!
!
900	return
	end
!
!****************************************************************************************
!							           COORC subroutine
!	This subroutine reads the coordinates of nodes in cartesian system.
!	This subroutine is called in: INPUT
!	This subroutine calls: -
!	variables used are:
!		xinc	: the increment step which is added to the reference x(n1) to create new
!				  node
!		yinc	: the increment step which is added to the reference y(n1) to create new
!				  node
!		n		: node number
!		ng		: increment code
!	INPUT	: nnp, x(2,nnp) = 0, nmiss = .false.
!	OUTPUT	: x(2,nnp), nmiss, istop
!****************************************************************************************
!
	subroutine coorc(nnp,x,nmiss)
	implicit double precision (a-h,o-z)
	character (5) :: card, wd(21)
	logical :: nmiss
!
	dimension x(2,1)	! here, the dimension
!
	common /ipt/ istop,list,wd
!
	n=0		! n is the number of node at each iteration
	ng=0	! ng is the value of increment which can be equal to 0 (off) or 1 (on).
!
1	read (4,*) card
	if (card(1:1).eq.'#') goto 1
	do i=1,list
		if (card.eq.wd(i)) goto 8
	enddo
	BACKSPACE (4)
!
	l=n
	lg=ng
!
2	read (4,*,err=3) n,ng,x(1,n),x(2,n)
	goto 4
3	write (*,*) '**Error** in reading in COORC'
	stop
!
4	if ((n.gt.nnp).OR.(n.le.0)) then
		write (*,10) n
		istop = istop +1
		goto 1
	endif
!
	if (x(1,n).eq.0.d0 .and. x(2,n).eq.0.d0) nmiss=.true.
!
	if (lg) 5,1,5
5	lg = SIGN (lg,n-l)
	li = (ABS (n-l+lg)-1)/ABS (lg)
	xinc = (x(1,n)-x(1,l))/li
	yinc = (x(2,n)-x(2,l))/li
6	l = l +lg
	if ((n-l)*lg .le. 0) goto 1
	if (l.le.0.OR.l.gt.nnp) then
		write (*,20) l
		istop = istop + 1
		goto 1
	endif
!
	x(1,l) = x(1,l-lg) + xinc
	x(2,l) = x(2,l-lg) + yinc
!
	if (x(1,l).eq.0.d0 .and. x(2,l).eq.0.d0) nmiss=.true.
	goto 6
!
8	BACKSPACE (4)
!
10	FORMAT (5x,'**Error** invalid node number:',i5,' in COORC')
20	FORMAT (5x,'**Error** attemp to generate node:',i5,' in COORC')
!
	return
	end
!
!
!****************************************************************************************
!							           COORP subroutine
!	This subroutine reads the coordinates of nodes in polar system and converts them
!	into the cartesian coordinate system.
!	VARIABLES	:
!					n1		: the first node (you must know the coordinate)
!					n2		: the last node (you must know the coordinate)
!					x(1,n1) : the x-coordinate of the first node
!					x(2,n1) : the y-coordinate of the first node
!					x(1,n2)	: the x-coordinate of the last node
!					x(2,n2) : the y-coordinate of the last node
!					lg		: the increment, in this subroutine, it must be equal to [1]
!					idir	: it denotes the direction required to reach the point from &
!							  the 0� ray or polar axis, [1,-1];
!							  [1] : positive or anticlockwise (counterclockwise) angle &
!							  [-1]: negative or clockwise angle
!					xo,yo	: the coordinate of central point known as the pole
!					r		: radial coordinate
!					p		: angular coordinate
!	INPUT		: nnp, x(2,nnp) = 0, nmiss = .false.
!	OUTPUT		: x(2,nnp), nmiss, istop
!****************************************************************************************
!
	subroutine coorp(nnp,x,nmiss)
	implicit double precision (a-h,o-z)
	character (5) :: card, wd(21)
	logical :: nmiss
!
	dimension x(2,1)	! here, the dimension
!
	common /ipt/ istop,list,wd
!
1	read (4,*) card
	if (card(1:1).eq.'#') goto 1
	do i=1,list
		if (card.eq.wd(i)) goto 8
	enddo
	BACKSPACE (4)
!
2	read (4,*,err=3) n1,n2,x(1,n1),x(2,n1),x(1,n2),x(2,n2),lg,idir,xo,yo
	goto 4
3	write (*,*) '**Error** in reading in COORP'
	stop
!
!	VERifYING THE data
4	if ((n1.gt.nnp).OR.(n1.le.0).OR.(n2.gt.nnp).OR.(n2.le.0)) then
		write (*,10) n1,n2
        istop = istop +1
		goto 1
	endif
!
	if ((x(1,n1).eq.xo.and.x(2,n1).eq.yo).OR.(x(1,n2).eq.xo.and.x(2,n2).eq.yo)) then
		write (*,20) n1,n2
		istop = istop +1
		goto 1
	endif
!
	if (lg.eq.0) then
		write (*,30) n1,n2,lg
		istop = istop +1
		goto 1
	endif
	if (idir.ne.-1 .and. idir.ne.1) then
		write (*,50) idir
		istop=istop+1
		goto 1
	endif
!
!	PREPARATION
	r1=dsqrt ((x(1,n1)-xo)**2+(x(2,n1)-yo)**2)
	p1=DACOS ((x(1,n1)-xo)/r1)
	if (x(2,n1).lt.yo) p1=2*DACOS(-1.d0)-p1
!
	r2=dsqrt((x(1,n2)-xo)**2+(x(2,n2)-yo)**2)
	p2=DACOS ((x(1,n2)-xo)/r2)
	if (x(2,n2).lt.yo) p2=2*DACOS (-1.d0)-p2
!
	lg = SIGN (lg,n2-n1)
	li = (ABS (n2-n1+lg)-1)/ABS(lg)
	rinc = (r2-r1)/li
	pinc = (p2-p1)/li
	if (idir.eq.-1) pinc=(p2-p1-2*DACOS(-1.d0))/li
!
	l= n1
	r= r1
	p= p1
!
!	LOOP
6	l = l +lg
	if ((n2-l)*lg .le. 0) goto 1
	if (l.le.0.OR.l.gt.nnp) then
		write (*,40) l
		istop = istop + 1
		goto 1
	endif
	r= r+rinc
	p= p+pinc
	x(1,l)=xo + r*DCOS(p)
	x(2,l)=yo + r*DSIN(p)
	if (x(1,l).eq.0.d0 .and. x(2,l).eq.0.d0) nmiss=.true.
	goto 6
8	BACKSPACE (4)
!
10	FORMAT (5x,'**Error** invalid node number:',2i5,' in COORP')
20	FORMAT (5x,'**Error** position of node:',2i5,' in COORP')
30	FORMAT (5x,'**Error** invalid step increment',3i5,' in COORP')
40	FORMAT (5x,'**Error** attemp to generate node:',i5,' in COORP')
50	FORMAT (5x,'**Error** invalid direction of generation',i5)
!
	return
	end
!
!
!****************************************************************************************
!							           COINT subroutine
!	This subroutine reads the coordinates of internal nodes in BE area in cartesian system.
!	This subroutine is called in: INPUT
!	This subroutine calls: -
!	variables used are:
!		nbeint	: number of internal nodes in BE zone
!		xint		: coordinates of internal nodes (nodal data) in BE D3 area;
!							xint(1,1:nnpbed3int) : X-coord
!							xint(2,1:nnpbed3int) : Y-coord
!		xinc	: the increment step which is added to the reference xint(n1) to create new
!				  node
!		yinc	: the increment step which is added to the reference yint(n1) to create new
!				  node
!		n		: node number
!		ng		: increment code
!	INPUT	: nbeint, xint(2,nbeint) = 0, nmiss = .false.
!	OUTPUT	: xint(2,nbeint), nmiss, istop
!****************************************************************************************
!
	subroutine coint(nbeint,xint,nmiss)
	implicit double precision (a-h,o-z)
	character (5) :: card, wd(21)
	logical :: nmiss
!
	dimension xint(2,1)	! here, the dimension
!
	common /ipt/ istop,list,wd
!
	n=0		! n is the number of node at each iteration
	ng=0	! ng is the value of increment which can be equal to 0 (off) or 1 (on).
!
1	read (4,*) card
	if (card(1:1).eq.'#') goto 1
	do i=1,list
		if (card.eq.wd(i)) goto 8
	enddo
	BACKSPACE (4)
!
	l=n
	lg=ng
!
2	read (4,*,err=3) n,ng,xint(1,n),xint(2,n)
	goto 4
3	write (*,*) '**Error** in reading in COORCINT'
	stop
!
4	if ((n.gt.nbeint).OR.(n.le.0)) then
		write (*,10) n
		istop = istop +1
		goto 1
	endif
!
	if (xint(1,n).eq.0.d0 .and. xint(2,n).eq.0.d0) nmiss=.true.
!
	if (lg) 5,1,5
5	lg = SIGN (lg,n-l)
	li = (ABS (n-l+lg)-1)/ABS (lg)
	xinc = (xint(1,n)-xint(1,l))/li
	yinc = (xint(2,n)-xint(2,l))/li
6	l = l +lg
	if ((n-l)*lg .le. 0) goto 1
	if (l.le.0.OR.l.gt.nbeint) then
		write (*,20) l
		istop = istop + 1
		goto 1
	endif
!
	xint(1,l) = xint(1,l-lg) + xinc
	xint(2,l) = xint(2,l-lg) + yinc
!
	if (xint(1,l).eq.0.d0 .and. xint(2,l).eq.0.d0) nmiss=.true.
	goto 6
!
8	BACKSPACE (4)
!
10	FORMAT (5x,'**Error** invalid internal node number:',i5,' in COORCINT')
20	FORMAT (5x,'**Error** attemp to generate internal node:',i5,' in COORCINT')
!
	return
	end
!
!
!****************************************************************************************
!							           ELM8D subroutine
!	This subroutine reads the drained element data.
!	VARIABLES	:
!					n: the number of element at each iteration
!					ng: the value of increment which can be equal to 0 (off) or 1 (on).
!					ie8d(1:8,n): the number of global nodes corresponding with the local nodes in counter-clockwise direction and in increasing order
!					ie8d(9,n): the material behaviour type (linear, hyperbolic, elastoplastic)
!	INPUT		: ne8d, ie8d(9,ne8d1) = 0
!	OUTPUT		: ie8d(9,ne8d1), istop
!****************************************************************************************
!
	subroutine elm8d(ne8d,ie8d)
	implicit double precision (a-h,o-z)
	character (5) :: card, wd(21)
!
	dimension ie8d(9,1)	! here, the dimension????????
!
	common /ipt/ istop,list,wd
!
	if (ne8d.eq.0) then
		write (*,*) '**Error** in ELM8D: THERE IS not DRAINED ELEMENT'
		istop=istop+1
		goto 6		! if there is not drained element, go to the consolidated elements
	endif
!
	n=0
	ng=0
!
1	read (4,*) card
	if (card(1:1).eq.'#') goto 1
	do i=1,list
		if (card.eq.wd(i)) goto 5
	enddo
	BACKSPACE (4)
!
	l=n
	lg=ng
!
	read (4,*,err=2) n,ng,(ie8d(i,n),i=1,9)
	goto 3
2	write (*,*) '**Error** in reading in ELM8D'
	stop
!
!	verifying data
3	if ((n.gt.ne8d).OR.(n.le.0)) then
		write (*,20) n
		istop = istop +1
		goto 1
	endif
!
	if (lg.eq.0) goto 1
	lg = SIGN (lg,n-l)
4	l=l+lg
	if (l.eq.n) goto 1
!
	do j=1,8
		ie8d(j,l)=ie8d(j,l-lg)+2
	enddo
	ie8d(4,l)=ie8d(4,l-lg)+1
	ie8d(8,l)=ie8d(8,l-lg)+1
	ie8d(9,l)=ie8d(9,l-lg)
	goto 4
!
5	BACKSPACE (4)
!
20	FORMAT (5x,'**Error** invalid element number',i5,' in ELM8D')
!
6	return
	end
!
!
!****************************************************************************************
!							           ELM8C subroutine
!	This subroutine reads the saturated element data.
!	VARIABLES	:
!					n: the number of element at each iteration
!					ng: the value of increment which can be equal to 0 (off) or 1 (on).
!					ie8c(1:8,n): the number of global nodes corresponding with the local nodes in counter-clockwise direction and in increasing order
!					ie8c(9,n): the material behaviour type (linear, hyperbolic, elastoplastic)
!	INPUT		: ne8c, ie8c(9,ne8c1) = 0
!	OUTPUT		: ie8c(1:9,ne8c1), istop
!****************************************************************************************
!
	subroutine elm8c(ne8c,ie8c)
	implicit double precision (a-h,o-z)
	character (5) :: card, wd(21)
	dimension ie8c(9,1)	! here, the dimension????????
!
	common /ipt/ istop,list,wd
!
	if (ne8c.eq.0) then
		write (*,*) '**Error** in ELM8C: THERE IS not SATURATED ELEMENT'
        istop=istop+1
		goto 6		! if there is not saturated element, go to the unsaturated elements
	endif
!
	n=0
	ng=0
1	read (4,*) card
	if (card(1:1).eq.'#') goto 1
	do i=1,list
		if (card.eq.wd(i)) goto 5
	enddo
	BACKSPACE (4)
!
	l=n
	lg=ng
!
	read (4,*,err=2) n,ng,(ie8c(i,n),i=1,9)
	goto 3
2	write (*,*) '**Error** in reading in ELM8C'
	stop
!
!	verifying data
3	if ((n.gt.ne8c).OR.(n.le.0)) then
		write (*,20) n
		istop = istop +1
		goto 1
	endif
!
	if (lg.eq.0) goto 1
	lg = SIGN (lg,n-l)
4	l=l+lg
	if (l.eq.n) goto 1
!
	do j=1,8
		ie8c(j,l)=ie8c(j,l-lg)+2	!ie8c(j,l)=ie8c(j,l-1)+2 ??
	enddo
	ie8c(4,l)=ie8c(4,l-lg)+1
	ie8c(8,l)=ie8c(8,l-lg)+1
	ie8c(9,l)=ie8c(9,l-lg)
	goto 4
!
5	BACKSPACE (4)
!
20	FORMAT (5x,'**Error** invalid element number',i5,' in ELM8C')
!
6	return
	end
!
!
!****************************************************************************************
!							           BLOCD subroutine
!	This subroutine reads & generates the 8 nodal drained super-element data.
!	avoid to define the nodes presented in BLOCD in COORC, COORP & ELM8D, in the other words
!	The difference between ELM8D & BLOCD is that:
!			ELM8D: the coordinates of all of the nodes are given in "coorc" in input file
!				   & then the vector "X" is made in COORC, completely.
!				   in "elm8d", we can see the carachtetistics of drained elements, for example
!				   the first element is composed of which nodes & ...
!			BLOCD: the coordinates of BE zone are given in "coorc" in input file
!				   & then some of the constituents of "X" is made in COORC.
!				   in "BLOCD", we can see the numbers of nodes of drained elements in FE zone,
!				   the number of elements in each side (x->1;y->2) and the material types of each elements.
!				   it generates all of the nodes, elements & material types in FE zone.
!	The relative positions of the reference nodes must be satisfied:
!								(x1,y1)-------------(x2,y2)
!									  |-------------|
!									  |-------------|
!									  |-------------|
!									  |-------------|
!								(x4,y4)-------------(x3,y3)
!	In these two cases, the order of numeration of nodes (mesh) is linear and layer
!								(n1) ---->----
!								     ---->----
!							    (n3) ---->----
!	VARIABLES	:
!					n1: number of first drained node in FE
!					e1: number of first drained element in FE
!					ne1: number of drained elements in direction x (to generate)
!					ne2: number of drained elements in direction y (to generate)
!					mat: the material behaviour type (linear, hyperbolic, elastoplastic) in drained elements
!					x1,y1: coordinates of node N1 in left-high extreme
!					x2,y2: coordinates of node in right-high extreme
!					x3,y3: coordinates of node in right-low extreme in second line
!					x4,y4: coordinates of node in left-low extreme in second line
!					ie8d(1:8,n): the number of global nodes corresponding with the local nodes in counter-clockwise direction and in increasing order
!					ie8d(9,n): the material behaviour type (linear, hyperbolic, elastoplastic)
!	INPUT		: nnp, x, ne8d, nmiss, ie8d(9,ne8d1) = 0
!	OUTPUT		: ie8d(9,ne8d1),x, istop
!****************************************************************************************
!
	subroutine blocd(nnp,x,ne8d,ie8d,nmiss)
	implicit double precision (a-h,o-z)
	character (5) :: card, wd(21)
	logical :: nmiss
!
	dimension x(2,1),ie8d(9,1)
!
	common /ipt/ istop,list,wd
!
	if (ne8d.eq.0) then
		write (*,*) '**Error** in BLOCD: THERE IS not DRAINED ELEMENT'
		istop=istop+1
		return
	endif
!
3	read (4,*) card
	if (card(1:1).eq.'#') goto 3
	do i=1,list
		if (card.eq.wd(i)) goto 4
	enddo
	BACKSPACE (4)
!
	read (4,*,err=1) n1,e1,ne1,ne2,mat,x1,y1,x2,y2,x3,y3,x4,y4
	goto 2
1	write (*,*) '**Error** in reading in BLOCD'
	stop
!
!	verifying data
2	if (ne1.lt.1 .or. ne2.lt.1) then
		write (*,10) ne1,ne2
		istop=istop + 1
		goto 3
	endif
!
	n3=3*ne1*ne2 + 2*ne2 + 2*ne1 +n1
!
!	then n3 in a 8 nodal element is the first node (in left side of elm) of third line
!
	if (n1.lt.1 .OR. n3.gt.nnp) then
		write (*,20) n1
		istop=istop + 1
		goto 3
	endif
	if (e1.lt.1 .OR. (e1+ne1*ne2-1).gt.ne8d) then
		write (*,30) e1
		istop=istop + 1
		goto 3
	endif
!
!	generation in according to side X1-X4 (vertical-left)
	n= n1
	l= n1
	x(1,n)=x1
	x(2,n)=y1
	if (x(1,n).eq.0.d0 .and. x(2,n).eq.0.d0) nmiss=.true.
	xinc = (x4-x1)/ne2
	yinc = (y4-y1)/ne2
	do i=1,ne2
		l=l+3*ne1+2
		x(1,l) = x(1,n) + xinc
        x(2,l) = x(2,n) + yinc
		if (x(1,l).eq.0.d0 .and. x(2,l).eq.0.d0) nmiss=.true.
		n=l
	enddo
	n= n1
	l= n+2*ne1+1
	xinc= xinc/2.
	yinc = yinc/2.
	x(1,l) = x(1,n) + xinc
    x(2,l) = x(2,n) + yinc
	do i=1, ne2-1
		l=l+3*ne1+2
		n=n+3*ne1+2
        x(1,l) = x(1,n) + xinc
        x(2,l) = x(2,n) + yinc
        if (x(1,l).eq.0.d0 .and. x(2,l).eq.0.d0) nmiss=.true.
	enddo
!
!	generation in according to side X2-X3 (vertical-right)
	n= n1+2*ne1
	l= n
    x(1,n)=x2
	x(2,n)=y2
    xinc = (x3-x2)/ne2
	yinc = (y3-y2)/ne2
	do i=1,ne2
		l=l+3*ne1+2
        x(1,l) = x(1,n) + xinc
        x(2,l) = x(2,n) + yinc
        if (x(1,l).eq.0.d0 .and. x(2,l).eq.0.d0) nmiss=.true.
		n=l
	enddo
	n= n1+2*ne1
	l= n+ne1+1
    xinc= xinc/2.
	yinc = yinc/2.
    x(1,l) = x(1,n) + xinc
    x(2,l) = x(2,n) + yinc
    do i=1, ne2-1
		l=l+3*ne1+2
        n=n+3*ne1+2
        x(1,l) = x(1,n) + xinc
        x(2,l) = x(2,n) + yinc
        if (x(1,l).eq.0.d0 .and. x(2,l).eq.0.d0) nmiss=.true.
	enddo
!
!	generation in according to side X1-X2 (horizontal-high) & X4-X3 (horizontal-low) & middle
	n= n1
	l= n1
    xinc = (x2-x1)/ne1/2.
	yinc = (y2-y1)/ne1/2.
	do i=1,2*ne1-1
		l=l+1
        x(1,l) = x(1,n) + xinc
        x(2,l) = x(2,n) + yinc
        if (x(1,l).eq.0.d0 .and. x(2,l).eq.0.d0) nmiss=.true.
		n=l
	enddo
	do i=1,ne2
		l=l+2
		n=l
        xinc = (x(1,l+ne1)-x(1,l))/ne1
        yinc = (x(2,l+ne1)-x(2,l))/ne1
        do j=1,ne1-1
			l=l+1
            x(1,l) = x(1,n) + xinc
            x(2,l) = x(2,n) + yinc
            if (x(1,l).eq.0.d0 .and. x(2,l).eq.0.d0) nmiss=.true.
            n=l
		enddo
	l=l+2
	n=l
    xinc = (x(1,l+2*ne1)-x(1,l))/ne1/2.
    yinc = (x(2,l+2*ne1)-x(2,l))/ne1/2.
		do j=1,2*ne1-1
			l=l+1
			x(1,l) = x(1,n) + xinc
            x(2,l) = x(2,n) + yinc
            if (x(1,l).eq.0.d0 .and. x(2,l).eq.0.d0) nmiss=.true.
            n=l
		enddo
	enddo
!
!	generation of elements
	l= e1
	n= n1
	do i=1,ne2
		k=2*ne1+1
		do j=1,ne1
			if (x2.gt.x1 .and. y3.lt.y2) then
				ie8d(5,l)=n+2
				ie8d(6,l)=n+1
				ie8d(7,l)=n
				ie8d(8,l)=n+k
				ie8d(4,l)=ie8d(8,l)+1
				ie8d(1,l)=n+3*ne1+2
				ie8d(2,l)=ie8d(1,l)+1
				ie8d(3,l)=ie8d(1,l)+2
			endif
			if (x2.lt.x1 .and. y3.lt.y2) then	! I AM not SURE, I doNT TRY IT!!
				ie8d(7,l)=n+2
				ie8d(6,l)=n+1
				ie8d(5,l)=n
				ie8d(4,l)=n+k
				ie8d(8,l)=ie8d(4,l)+1
				ie8d(3,l)=n+3*ne1+2
				ie8d(2,l)=ie8d(3,l)+1
				ie8d(1,l)=ie8d(3,l)+2
			endif
			if (x2.gt.x1 .and. y3.gt.y2) then   ! I AM not SURE, I doNT TRY IT!!
				ie8d(3,l)=n+2
				ie8d(2,l)=n+1
				ie8d(1,l)=n
				ie8d(8,l)=n+k
				ie8d(4,l)=ie8d(8,l)+1
				ie8d(7,l)=n+3*ne1+2
				ie8d(6,l)=ie8d(7,l)+1
				ie8d(5,l)=ie8d(7,l)+2
			endif
			if (x2.lt.x1 .and. y3.gt.y2) then   ! I AM not SURE, I doNT TRY IT!!
				ie8d(1,l)=n+2
				ie8d(2,l)=n+1
				ie8d(3,l)=n
				ie8d(4,l)=n+k
				ie8d(8,l)=ie8d(4,l)+1
				ie8d(5,l)=n+3*ne1+2
				ie8d(6,l)=ie8d(5,l)+1
				ie8d(7,l)=ie8d(5,l)+2
			endif
			ie8d(9,l)=mat
            k=k-1
			l=l+1
			n=n+2
		enddo
		n=n1+i*(3*ne1+2)
	enddo
	goto 3
!
10	FORMAT (5x,'**Error** invalid number of element of each cote',&
				i5,i5,' in BLOCD')
20	FORMAT (5x,'**Error** invalid node number',&
                i5,' in BLOCD')
30	FORMAT (5x,'**Error** invalid element number',&
                i5,' in BLOCD')
!
4	BACKSPACE (4)
	return
	end
!
!
!****************************************************************************************
!							           BLOCC subroutine
!	This subroutine reads & generates the 8 nodal saturated super-element data.
!	avoid to define the nodes presented in BLOCC in COORC, COORP & ELM8C, in the other words
!	The difference between ELM8C & BLOCC is that:
!			ELM8C: the coordinates of all of the nodes are given in "coorc" in input file
!				   & then the vector "X" is made in COORC, completely.
!				   in "elm8c", we can see the carachtetistics of saturated elements, for example
!				   the first element is composed of which nodes & ...
!			BLOCC: the coordinates of BE zone are given in "coorc" in input file
!				   & then some of the constituents of "X" is made in COORC.
!				   in "BLOCC", we can see the numbers of nodes of saturated elements in FE zone,
!				   the number of elements in each side (x->1;y->2) and the material types of each elements.
!				   it generates all of the nodes, elements & material types in FE zone.
!	The relative positions of the reference nodes must be satisfied:
!								(x1,y1)-------------(x2,y2)
!									  |-------------|
!									  |-------------|
!									  |-------------|
!									  |-------------|
!								(x4,y4)-------------(x3,y3)
!	In these two cases, the order of numeration of nodes (mesh) is linear and layer
!								(n1) ---->----
!								     ---->----
!							    (n3) ---->----
!	VARIABLES	:
!					n1: number of first saturated node in FE
!					e1: number of first saturated element in FE
!					ne1: number of saturated elements in direction x (to generate)
!					ne2: number of saturated elements in direction y (to generate)
!					mat: the material behaviour type (linear, hyperbolic, elastoplastic) in saturated elements
!					x1,y1: coordinates of node N1 in left-high extreme
!					x2,y2: coordinates of node in right-high extreme
!					x3,y3: coordinates of node in right-low extreme in second line
!					x4,y4: coordinates of node in left-low extreme in second line
!					ie8c(1:8,n): the number of global nodes corresponding with the local nodes in counter-clockwise direction and in increasing order
!					ie8c(9,n): the material behaviour type (linear, hyperbolic, elastoplastic)
!	INPUT		: nnp, x, ne8c, nmiss, ie8c(9,ne8c1) = 0
!	OUTPUT		: ie8c(9,ne8c1),x, istop
!****************************************************************************************
!
	subroutine blocc(nnp,x,ne8c,ie8c,nmiss)
	implicit double precision (a-h,o-z)
	character (5) :: card, wd(21)
	logical :: nmiss
!
	dimension x(2,1),ie8c(9,1)
!
	common /ipt/ istop,list,wd
!
	if (ne8c.eq.0) then
		write (*,*) '**Error** in BLOCC: THERE IS not SATURATED ELEMENT'
		istop=istop+1
		return
	endif
!
3	read (4,*) card
	if (card(1:1).eq.'#') goto 3
	do i=1,list
		if (card.eq.wd(i)) goto 4
	enddo
	BACKSPACE (4)
!
	read (4,*,err=1) n1,e1,ne1,ne2,mat,x1,y1,x2,y2,x3,y3,x4,y4
	goto 2
1	write (*,*) '**Error** in reading in BLOCC'
	stop
!
!	verifying data
2	if (ne1.lt.1 .or. ne2.lt.1) then
		write (*,10) ne1,ne2
		istop=istop + 1
		goto 3
	endif
!
	n3=3*ne1*ne2 + 2*ne2 + 2*ne1 +n1
!
!	then n3 in a 8 nodal element is the first node (in left side of elm) of third line
!
	if (n1.lt.1 .OR. n3.gt.nnp) then
		write (*,20) n1
		istop=istop + 1
		goto 3
	endif
	if (e1.lt.1 .OR. (e1+ne1*ne2-1).gt.ne8c) then
		write (*,30) e1
		istop=istop + 1
		goto 3
	endif
!
!	generation in according to side X1-X4 (vertical-left)
	n= n1
	l= n1
	x(1,n)=x1
	x(2,n)=y1
	if (x(1,n).eq.0.d0 .and. x(2,n).eq.0.d0) nmiss=.true.
	xinc = (x4-x1)/ne2
	yinc = (y4-y1)/ne2
	do i=1,ne2
		l=l+3*ne1+2
		x(1,l) = x(1,n) + xinc
        x(2,l) = x(2,n) + yinc
		if (x(1,l).eq.0.d0 .and. x(2,l).eq.0.d0) nmiss=.true.
		n=l
	enddo
	n= n1
	l= n+2*ne1+1
	xinc= xinc/2.
	yinc = yinc/2.
	x(1,l) = x(1,n) + xinc
    x(2,l) = x(2,n) + yinc
	do i=1, ne2-1
		l=l+3*ne1+2
		n=n+3*ne1+2
        x(1,l) = x(1,n) + xinc
        x(2,l) = x(2,n) + yinc
        if (x(1,l).eq.0.d0 .and. x(2,l).eq.0.d0) nmiss=.true.
	enddo
!
!	generation in according to side X2-X3 (vertical-right)
	n= n1+2*ne1
	l= n
    x(1,n)=x2
	x(2,n)=y2
    xinc = (x3-x2)/ne2
	yinc = (y3-y2)/ne2
	do i=1,ne2
		l=l+3*ne1+2
        x(1,l) = x(1,n) + xinc
        x(2,l) = x(2,n) + yinc
        if (x(1,l).eq.0.d0 .and. x(2,l).eq.0.d0) nmiss=.true.
		n=l
	enddo
	n= n1+2*ne1
	l= n+ne1+1
    xinc= xinc/2.
	yinc = yinc/2.
    x(1,l) = x(1,n) + xinc
    x(2,l) = x(2,n) + yinc
    do i=1, ne2-1
		l=l+3*ne1+2
        n=n+3*ne1+2
        x(1,l) = x(1,n) + xinc
        x(2,l) = x(2,n) + yinc
        if (x(1,l).eq.0.d0 .and. x(2,l).eq.0.d0) nmiss=.true.
	enddo
!
!	generation in according to side X1-X2 (horizontal-high) & X4-X3 (horizontal-low) & middle
	n= n1
	l= n1
    xinc = (x2-x1)/ne1/2.
	yinc = (y2-y1)/ne1/2.
	do i=1,2*ne1-1
		l=l+1
        x(1,l) = x(1,n) + xinc
        x(2,l) = x(2,n) + yinc
        if (x(1,l).eq.0.d0 .and. x(2,l).eq.0.d0) nmiss=.true.
		n=l
	enddo
	do i=1,ne2
		l=l+2
		n=l
        xinc = (x(1,l+ne1)-x(1,l))/ne1
        yinc = (x(2,l+ne1)-x(2,l))/ne1
        do j=1,ne1-1
			l=l+1
            x(1,l) = x(1,n) + xinc
            x(2,l) = x(2,n) + yinc
            if (x(1,l).eq.0.d0 .and. x(2,l).eq.0.d0) nmiss=.true.
            n=l
		enddo
	l=l+2
	n=l
    xinc = (x(1,l+2*ne1)-x(1,l))/ne1/2.
    yinc = (x(2,l+2*ne1)-x(2,l))/ne1/2.
		do j=1,2*ne1-1
			l=l+1
			x(1,l) = x(1,n) + xinc
            x(2,l) = x(2,n) + yinc
            if (x(1,l).eq.0.d0 .and. x(2,l).eq.0.d0) nmiss=.true.
            n=l
		enddo
	enddo
!
!	generation of elements
	l= e1
	n= n1
	do i=1,ne2
		k=2*ne1+1
		do j=1,ne1
			if (x2.gt.x1 .and. y3.lt.y2) then
				ie8c(5,l)=n+2
				ie8c(6,l)=n+1
				ie8c(7,l)=n
				ie8c(8,l)=n+k
				ie8c(4,l)=ie8c(8,l)+1
				ie8c(1,l)=n+3*ne1+2
				ie8c(2,l)=ie8c(1,l)+1
				ie8c(3,l)=ie8c(1,l)+2
			endif
			if (x2.lt.x1 .and. y3.lt.y2) then	! I AM not SURE, I doNT TRY IT!!
				ie8c(7,l)=n+2
				ie8c(6,l)=n+1
				ie8c(5,l)=n
				ie8c(4,l)=n+k
				ie8c(8,l)=ie8c(4,l)+1
				ie8c(3,l)=n+3*ne1+2
				ie8c(2,l)=ie8c(3,l)+1
				ie8c(1,l)=ie8c(3,l)+2
			endif
			if (x2.gt.x1 .and. y3.gt.y2) then   ! I AM not SURE, I doNT TRY IT!!
				ie8c(3,l)=n+2
				ie8c(2,l)=n+1
				ie8c(1,l)=n
				ie8c(8,l)=n+k
				ie8c(4,l)=ie8c(8,l)+1
				ie8c(7,l)=n+3*ne1+2
				ie8c(6,l)=ie8c(7,l)+1
				ie8c(5,l)=ie8c(7,l)+2
			endif
			if (x2.lt.x1 .and. y3.gt.y2) then   ! I AM not SURE, I doNT TRY IT!!
				ie8c(1,l)=n+2
				ie8c(2,l)=n+1
				ie8c(3,l)=n
				ie8c(4,l)=n+k
				ie8c(8,l)=ie8c(4,l)+1
				ie8c(5,l)=n+3*ne1+2
				ie8c(6,l)=ie8c(5,l)+1
				ie8c(7,l)=ie8c(5,l)+2
			endif
			ie8c(9,l)=mat
            k=k-1
			l=l+1
			n=n+2
		enddo
		n=n1+i*(3*ne1+2)
	enddo
	goto 3
!
10	FORMAT (5x,'**Error** invalid number of element of each cote',&
				i5,i5,' in BLOCC')
20	FORMAT (5x,'**Error** invalid node number',&
                i5,' in BLOCC')
30	FORMAT (5x,'**Error** invalid element number',&
                i5,' in BLOCC')
!
4	BACKSPACE (4)
	return
	end
!
!
!****************************************************************************************
!									BEMD1 subroutine(PADSA@GARD)
!	This subroutine reads the nodes of finite zones (D1) of BE.
!	To obtaine the results of the interior nodes (in domain), we cut into several zones
!	in which the nodes in boundary elms are the interior nodes.
!	This subroutine is called in: INPUT
!	This subroutine calls: -
!	variables used are:
!		ie3d1	: ie3d1(1:nnpbed1,1:nbed): nodes' number which form the current D1 BEM area
!										   Numbering direction for external surface close
!										   (finite) domain is ANTICLOCKWISE
!				  ie3d1(nnpbed1(i)+1,1:nbed): it shows the number of material which is used
!											  is used: 1:m8d;1:m8c;1:m8u
!		npbmax	: maximum nodes in all of the BE zone + 1; for example if we consider D3 BE
!				  area NPBMAX = nnpbed3 + 1
!		n: number of node
!		ng: increment code
!		izone: number of BEM D1 zone
!	INPUT		: ibem, nnp, nbed1, nnpbed1, npbmax
!	OUTPUT		: ie3d1(1,nbed11)
!****************************************************************************************
!
	subroutine bemd1(ibem,nnp,nbed1,nnpbed1,ie3d1,npbmax)
	implicit double precision (a-h,o-z)
	character (5) :: card, wd(21)
!
	dimension nnpbed1(5),ie3d1(npbmax,nbed1)
	dimension l(nbed1),lg(nbed1),j(nbed1)
!
	common /ipt/ istop,list,wd
!
	if (nbed1.eq.0 .OR. ibem.eq.0) then
		write (*,*) '**Error** in BEMD1: THERE IS not BEM D1 (FINITE) doMAIN'
		istop=istop+1
		return
	endif
!
	read (4,*,err=30) (ie3d1(nnpbed1(i)+1,i),i=1,nbed1)
!
	do i=1,nbed1
		l(i)=0
		lg(i)=0
		j(i)=0
	enddo
!
	n=0
	ng=0
	izone=1
!
10	read (4,*) card
	if (card(1:1).eq.'#') goto 10
	do i=1,list
		if (card.eq.wd(i)) goto 90
	enddo
	BACKSPACE (4)
!
	l(izone)=n
	lg(izone)=ng
!
	read (4,*,err=30) n,ng,izone
	goto 40
30	write (*,*) '**Error** in reading in BEMD1'
	stop
40	if ((n.gt.nnp).OR.(n.le.0)) then
		write (*,1) n
		istop = istop +1
		goto 10
	endif
	if (izone.lt.1 .OR. izone.gt.nbed1) then
		write (*,3) izone
		istop=istop+1
		goto 10
	endif
!
	if (lg(izone)) 50,80,50
!
50	lg(izone) = SIGN (lg(izone),n-l(izone))
60	if ((n-l(izone))*lg(izone) .le. 0) goto 10
    if (l(izone).le.0 .OR. l(izone).gt.nnp) then
		write (*,2) l(izone)
		istop = istop + 1
		goto 10
	endif
80	l(izone)=l(izone)+lg(izone)
	j(izone)=j(izone)+1
	if (j(izone).gt.nnpbed1(izone)) then
		write (*,4) izone,n
		istop=istop+1
		goto 10
	endif
	if (lg(izone).eq.0) then
		ie3d1(j(izone),izone)=n
	else
        ie3d1(j(izone),izone)=l(izone)
	endif
	goto 60
90	BACKSPACE (4)
!
	do izone=1,nbed1
	if (j(izone).lt.nnpbed1(izone)) then
		write (*,5) nnpbed1(izone)-j(izone),izone
		istop=istop+1
	endif
	enddo
!
1	FORMAT (5x,'**Error** invalid node number:',i5,' in BEMD1')
2	FORMAT (5x,'**Error** attemp to generate node:',i5,' in BEMD1')
3	FORMAT (5x,'**Error** invalid zone BEM d1 number',i5)
4	FORMAT (5x,'**Error** exceed the number of',i3,'th BEM d1 node:',i5)
5	FORMAT (5x,'**Error**',i4,' missing nodes of',i3,'th BEM d1 domain')
!
	return
	end
!
!
!****************************************************************************************
!								   BEMD2 subroutine(SA@GARD)
!	This subroutine reads the nodes of infinite zones (D2) of BE.
!	This subroutine is called in: INPUT
!	This subroutine calls: -
!	variables used are:
!		ie3d2	: ie3d2(1:nnpbed2,1:nbed): nodes' number which form the current D2 BEM area
!										   Numbering direction for internal surface open
!										   (infinite) domain is CLOCKWISE
!				  ie3d2(nnpbed2(i)+1,1:nbed): it shows the type of just one material which
!											  is used:
!							[1]: dry boundary elm
!							[2]: saturated boundary elm
!							[3]: unsaturated boundary elm
!		npbmax	: maximum nodes in all of the BE zone + 1; for example if we consider D3 BE
!				  area NPBMAX = nnpbed3 + 1
!		n: number of node
!		ng: increment code
!		izone: number of BEM D2 zone
!	INPUT		: ibem, nnp, nbed2, nnpbed2, npbmax
!	OUTPUT		: ie3d2(1,nbed21)
!****************************************************************************************
!
	subroutine bemd2(ibem,nnp,nbed2,nnpbed2,ie3d2)
	implicit double precision (a-h,o-z)
	character (5) :: card, wd(21)
	logical :: 	strcomp
!
	dimension nnpbed2(5),ie3d2(1,1)
	dimension l(nbed2),lg(nbed2),j(nbed2)
!
	common /ipt/ istop,list,wd
!
	if (nbed2.eq.0 .OR. ibem.eq.0) then
		write (*,*) '**Error** in BEMD2: THERE IS not BEM D2 (INFINITE) doMAIN'
		istop=istop+1
		return
	endif
!
	read (4,*,err=30) (ie3d2(nnpbed2(i)+1,i),i=1,nbed2)
!
	do i=1,nbed2
		l(i)=0
		lg(i)=0
		j(i)=0
	enddo
!
	n=0
	ng=0
	izone=1
!
10	read (4,*) card
	if (card(1:1).eq.'#') goto 10
	do i=1,list
		if (strcomp(card,wd(i))) goto 90
	enddo
	BACKSPACE (4)
!
	l(izone)=n
	lg(izone)=ng
!
	read (4,*,err=30) n,ng,izone
	goto 40
30	write (*,*) '**Error** in reading in BEMD2'
	stop
40	if ((n.gt.nnp).OR.(n.le.0)) then
		write (*,1) n
		istop = istop +1
		goto 10
	endif
	if (izone.lt.1 .OR. izone.gt.nbed2) then
		write (*,3) izone
		istop=istop+1
		goto 10
	endif
!
	if (lg(izone)) 50,80,50
!
50	lg(izone) = SIGN (lg(izone),n-l(izone))
60	if ((n-l(izone))*lg(izone) .le. 0) goto 10
	if (l(izone).le.0 .OR. l(izone).gt.nnp) then
		write (*,2) l(izone)
		istop = istop + 1
		goto 10
	endif
80	l(izone)=l(izone)+lg(izone)
	j(izone)=j(izone)+1
	if (j(izone).gt.nnpbed2(izone)) then
		write (*,4) izone,n
		istop=istop+1
		goto 10
	endif
	if (lg(izone).eq.0) then
		ie3d2(j(izone),izone)=n
	else
		ie3d2(j(izone),izone)=l(izone)
	endif
	goto 60
90	BACKSPACE (4)
!
	do izone=1,nbed2
		if (j(izone).lt.nnpbed2(izone)) then
			write (*,5) nnpbed2(izone)-j(izone),izone
            istop=istop+1
		endif
	enddo
!
1	FORMAT (5x,'**Error** invalid node number:',i5,' in BEMD2')
2	FORMAT (5x,'**Error** attemp to generate node:',i5,' in BEMD2')
3   FORMAT (5x,'**Error** invalid zone BEM d2 number',i5)
4   FORMAT (5x,'**Error** exceed the number of',i3,'th BEM d2 node:',i5)
5   FORMAT (5x,'**Error**',i4,' missing nodes of',i3,'th BEM d2 domain')
!
	return
	end
!
!
!****************************************************************************************
!			                       BEMD3 subroutine(PADSA@GARD)
!	This subroutine reads the nodes of semi-infinite zones (D3) of BE.
!	This subroutine is called in: INPUT
!	This subroutine calls: -
!	variables used are:
!		ie3d3	: ie3d3(1:nnpbed3): nodes' number which form the current D3 BEM area
!										   Numbering direction for semi-infinite domain is
!										   ANTICLOCKWISE
!				  ie3d3(nnpbed3+1): it shows the type of just one material which
!											  is used:
!							[1]: dry boundary elm
!							[2]: saturated boundary elm
!							[3]: unsaturated boundary elm
!		n: number of node
!		ng: increment code
!	INPUT		: ibem, nnp, nbed3, nnpbed3
!	OUTPUT		: ie3d3(1:nnpbed3+2)
!****************************************************************************************
!
	subroutine bemd3(ibem,nnp,nbed3,nnpbed3,ie3d3)
	implicit double precision (a-h,o-z)
	character (5) :: card, wd(21)
!
	dimension ie3d3(1)
!
	common /ipt/ istop,list,wd
!
	if (nbed3.eq.0 .OR. ibem.eq.0) then
		write (*,*) '**Error** in BEMD3: THERE IS not BEM D3 (SEMI-INFINITE) doMAIN'
        istop=istop+1
		return
	endif
!
	read (4,*,err=30) ie3d3(nnpbed3+1)	! when ie3d3(nnpbed3+1).eq.1 -> dry boundary elm
!										  when ie3d3(nnpbed3+1).eq.2 -> saturated boundary elm
!										  when ie3d3(nnpbed3+1).eq.3 -> unsaturated boundary elm
!	read (4,*,err=30) ie3d3(nnpbed3+2)	! it shows the number (num�ro) of (dry,sat,unsat)
!										  material which is used [m8d,m8c,m8u]
!
	n=0
	ng=0
	j=0
!
10	read (4,*) card
	if (card(1:1).eq.'#') goto 10
	do i=1,list
		if (card.eq.wd(i)) goto 90
	enddo
	BACKSPACE (4)
!
	l=n
	lg=ng
	read (4,*,err=30) n,ng
	goto 40
!
30	write (*,*) '**Error** in reading in BEMD3'
	stop
!
40	if ((n.gt.nnp).OR.(n.le.0)) then
		write (*,1) n
		istop = istop +1
		goto 10
	endif
!
	if (lg) 50,80,50
50	lg = SIGN (lg,n-l)
60	if ((n-l)*lg .le. 0) goto 10
	if (l.le.0 .OR. l.gt.nnp) then
		write (*,2) l
		istop = istop + 1
		goto 10
	endif
80	l=l+lg
    j=j+1
	if (j.gt.nnpbed3) then
		write (*,4) n
		istop=istop+1
		goto 10
	endif
	if (lg.eq.0) then
		ie3d3(j)=n
	else
		ie3d3(j)=l
	endif
	goto 60
90	BACKSPACE (4)
	if (j.lt.nnpbed3) then
		write (*,5) nnpbed3-j
		istop=istop+1
	endif
!
1	FORMAT (5x,'**Error** invalid node number:',i5,' in BEMD3')
2	FORMAT (5x,'**Error** attemp to generate node:',i5,' in BEMD3')
4	FORMAT (5x,'**Error** exceed the number of BEM d3 node:',i5)
5   FORMAT (5x,'**Error**',i4,' missing nodes of BEM d3 domain')
!
	return
	end
!
!
!****************************************************************************************
!			                            BEMEN subroutine
!	This subroutine reads the nodes of enclosing elements of BE.
!	The mesh of enclosing elements is different from the BEM mesh & general nodal coordinates
!	in x & elm3d.
!	This subroutine is called in: INPUT
!	This subroutine calls: -
!	variables used are:
!		nnpen	: number of nodes of fictitious enclosing elements in D3 BEM area, it can be
!				  at maximum 500 nodes
!		xencl: X-coordinates of enclosing elements' nodes (nodal data)
!		yencl: Y-coordinates of enclosing elements' nodes (nodal data)
!		n: number of node
!		ng: increment code
!	INPUT		: nnpen
!	OUTPUT		: xencl, yencl
!****************************************************************************************
!
	subroutine bemen(nnpen,xencl,yencl)
	implicit double precision (a-h,o-z)
	character (5) :: card, wd(21)
!
	dimension xencl(1),yencl(1)
!
	common /ipt/ istop,list,wd
!
!	read (4,*,err=3) nnpen
!
	n=0
	ng=0
	j=0
!
1	read (4,*) card
	if (card(1:1).eq.'#') goto 1
	do i=1,list
		if (card.eq.wd(i)) goto 8
	enddo
	BACKSPACE (4)
!
	l=n
	lg=ng
2	read (4,*,err=3) n,ng,xencl(n),yencl(n)
	j=j+1
	goto 4
3	write (*,*) '**Error** in reading in BEMEN'
	stop
!
4	if ((n.gt.nnpen).OR.(n.le.0)) then
		write (*,10) n
		istop = istop +1
        goto 1
	endif
!
	if (lg) 5,1,5
5	lg = SIGN (lg,n-l)
	li = (ABS (n-l+lg)-1)/ABS (lg)
	xinc = (xencl(n)-xencl(l))/li
	yinc = (yencl(n)-yencl(l))/li
!
6	l = l +lg
	if ((n-l)*lg .le. 0) goto 1
	if (l.le.0.OR.l.gt.nnpen) then
		write (*,20) l
		istop = istop + 1
		goto 1
	endif
	xencl(l)=xencl(l-lg)+xinc
	yencl(l)=yencl(l-lg)+yinc
	j=j+1
	goto 6
8	BACKSPACE (4)
	if (j.lt.nnpen) then
		write (*,30) nnpen-j
		istop=istop+1
	endif
!
10	FORMAT (5x,'**Error** invalid node number:',i5,' in BEMEN')
20	FORMAT (5x,'**Error** attemp to generate node:',i5,' in BEMEN')
30	FORMAT (5x,'**Error**',i4,' missing nodes in BEMEN')
!
	return
	end
!
!
!****************************************************************************************
!							           MATED subroutine
!	This subroutine reads the properties of drained materials.
!	This subroutine is called in: INPUT
!	This subroutine calls: -
!	variables used are:
!		m8d		: number of drained (dry) materials in FE & BE zones; in BE area the program
!				  is written for just one material (homogeneous space) then if we use a dry
!				  material in BE zone (by considering the dry materials in FE zone) the
!				  number of material of BE zone is always equal to M8D because the materials
!			      from "1" to "M8D-1" are used in FE zone
!		m8dep	: behaviour type of drained (dry) materials [1,2,3,4];
!							[1]: linear,
!							[2]: hyperbolic for statics,
!							[3]: hyperbolic for dynamics,
!							[4]: elastoplastic
!				  if we use a dry material in BE zone then the behaviour type of BE zone's
!				  material is: M8DEP(m8d)
!		icpt	: behaviour type code [1,2,3,4]
!							[1]linear
!							[2]hyperbolic for statics (Duncan & Chan)
!							[3]hyperbolic for dynamics (Hardin & Drnewich),
!							[4]elastoplastic (Prevost/1985) [THIS SECTION IS UNDER
!							   CONSTRUCTION BY MR.GHANEFAR AT UNIVERSITY OF TEHRAN]
!		imat	: number (num�ro) of drained materials
!		sm8d	: all of the parameters corresponding to each behaviour of dry soil are
!				  defined in this matrix.
!						  icpt:[1]-> sm8d (3,i) = E
!									 sm8d (7,i) = K
!									 sm8d (10,i)= Ko
!									 sm8d (11,i)= gammas
!						  icpt:[2]-> sm8d (1,i) = friction angle (�)
!									 sm8d (2,i) = cohesion (Pa)
!									 sm8d (3,i) = E-load, loading modulus used in non-linear
!												  elasticity (adimensionel)
!									 sm8d (4,i) = E-unld, unloading modulus used in
!												  non-linear elasticity (adim)
!									 sm8d (5,i) = n, exponent used to compute the Young
!												  loading and unloading moduli in non-linear
!											      elasticity (adim)
!									 sm8d (6,i) = Rf, coefficient multiplying the ratio
!												  (sigm1-sigm3)/(sigm1-sigm3)ult
!												  in nlin els (adim)
!                                    sm8d (7,i) = Kb, bulk modulus used in non-linear
!											      elasticity (adim)
!	                                 sm8d (8,i) = m, exponent used to compute the bulk
!										          modulus B in nlin els (adim)
!									 sm8d (9,i) = e-f, Young modulus minimal value (Pa)
!									 sm8d (10,i)=Ko, steady state coefficient (adim)
!								     sm8d (11,i)= gammas, soil volumetric weight (N,m-3)
!								     sm8d (12,i)= ten-st, traction resistency (Pa)
!						  icpt:[3]-> sm8d (1,i) = friction angle
!									 sm8d (2,i) = cohesion
!								     sm8d (3,i) = OCR
!								     sm8d (4,i) = k
!                                    sm8d (5,i) = n, exponent used to compute the Young
!												  loading and unloading moduli in non-linear
!												  elasticity (adim)
!									 sm8d (6,i) = e
!									 sm8d (7,i) = Kb, bulk modulus used in non-linear
!												  elasticity (adim)
!									 sm8d (8,i) = m, exponent used to compute the bulk
!												  modulus B in nlin els (adim)
!									 sm8d (10,i) =Ko
!									 sm8d (11,i) =gamma
!						  icpt:[4]-> ??????
!	INPUT		: ibem, ne8d, m8d, m8dep
!	OUTPUT		: sm8d
!****************************************************************************************
!
	subroutine mated(m8d,sm8d,m8dep,ibem,ne8d)
	implicit double precision (a-h,o-z)
	character (5) :: card, wd(21)
!
	dimension sm8d(60,1),m8dep(1)
!
	common /ipt/ istop,list,wd
!
!
	if (ibem.ne.1 .and. ne8d.eq.0) then
		write (*,*) '**Error** in MATED: THERE IS not DRAINED MATERIAL'
		istop=istop+1
		return
	endif
!
	do i=1,m8d
		do j=1,60
			sm8d(j,i)=0.d0
		enddo
	enddo
	do 70 i=1,m8d
		read (4,*,err=10) imat,icpt
		goto 20
10		write (*,*) '**Error** in reading in MATED'
		stop
20		if (imat.lt.1.OR.imat.gt.m8d.OR.icpt.lt.1.OR.icpt.gt.4) then
			write (*,1) imat,icpt
			istop=istop+1
			goto 70
		endif
!
		m8dep(imat)=icpt
		goto (30,40,50,60) icpt
!
!	  linear elastic: icpt=[1]
30		read (4,*,err=10) sm8d(3,imat),sm8d(7,imat),sm8d(10,imat),sm8d(11,imat)
		goto 70
!
!	  hyperbolic non-linear elastic for statics: icpt=[2]
40		read (4,*,err=10) (sm8d(j,imat),j=1,8)
		sm8d(1,imat)=sm8d(1,imat)*1.74532925199432d-2 !??????????????
		read (4,*,err=10) (sm8d(j,imat),j=9,12)
		goto 70
!
!	  hyperbolic non-linear elastic for for dynamics: icpt=[3]
50		read (4,*,err=10) (sm8d(j,imat),j=1,8)
		sm8d(1,imat)=sm8d(1,imat)*1.74532925199432d-2 !??????????????
		read (4,*,err=10) sm8d(10,imat),sm8d(11,imat)
		goto 70
!
!	  Prevost elasto-plastic 1985: icpt=[4]
60		read (4,*,err=10) (sm8d(j,imat),j=1,7)
		sm8d(1,imat)=sm8d(1,imat)*1.74532925199432d-2  !??????????????
		read (4,*,err=10) (sm8d(j,imat),j=9,12)
		nys=sm8d(9,imat)
		do k=1,nys
			j=4*k
			read (4,*,err=10) sm8d(17+j,imat),sm8d(18+j,imat),sm8d(19+j,imat),sm8d(20+j,imat)
		enddo
70	continue
!
1	FORMAT ('**Error** invalid materiel number',2i5,' in MATED')
	return
	end
!
!
!****************************************************************************************
!							           MATEC subroutine
!	This subroutine reads the properties of saturated materials.
!	This subroutine is called in: INPUT
!	This subroutine calls: -
!	variables used are:
!		m8c		: number of consolidated (saturated) materials in FE & BE zones; in BE area
!				  the program is written for just one material (homogeneous space) then if
!				  we use a saturated material in BE zone (by considering the saturated
!				  materials in FE zone) the number of material of BE zone is always equal to
!				  M8C because the materials from "1" to "M8C-1" are used in FE zone
!		m8cep	: behaviour type of saturated materials [1,2,3,4];
!							[1]: linear,
!							[2]: hyperbolic for statics,
!							[3]: hyperbolic for dynamics,
!							[4]: elastoplastic
!				  if we use a saturated material in BE zone then the behaviour type of BE
!				  zone's material is: M8CEP(m8c)
!		ne8c	: number of consolidated (saturated) elements in FE zone
!		icpt	: behaviour type code [1,2,3,4]
!							[1]linear
!							[2]hyperbolic for statics (Duncan & Chan)
!							[3]hyperbolic for dynamics (Hardin & Drnewich),
!							[4]elastoplastic (Prevost/1985) [THIS SECTION IS UNDER
!							   CONSTRUCTION BY MR.GHANEFAR AT UNIVERSITY OF TEHRAN]
!		imat	: number (num�ro) of saturated materials
!		sm8c	: all of the parameters corresponding to each behaviour of saturated soil
!				  are defined in this matrix.
!						  icpt:[1]-> sm8c (3,i) = G  !E
!									 sm8c (7,i) = NU !K
!									 sm8c (10,i)= PHI !Ko
!									 sm8c (11,i)= - !gamma
!									 sm8c (13,i) = alpha, 1-Kd/Ks
!									 sm8c (14,i)= NUU
!								     sm8c (15,i) = kz, vertical water permeability (m.s^-1)
!									 sm8c (16,i) = ! kx/kz, ratio of the horizontal water
!												   permeability to the vertical water
!												   permeability.
!					DYNAMIC PROBLEM BY CQM
!
!						  icpt:[1]-> sm8c (3,i) = E
!									 sm8c (7,i) = NU
!									 sm8c (10,i)= PHI (POROSITY)
!									 sm8c (11,i)= DENSITY OF SOIL
!									 sm8c (12,i)= R
!									 sm8c (13,i)= alpha, 1-Kd/Ks
!								     sm8c (15,i)= kz, vertical water permeability (m^4/N.s)
!
!						  icpt:[2]-> sm8c (1,i) = friction angle (�)
!									 sm8c (2,i) = cohesion (Pa)
!									 sm8c (3,i) = E-load, loading modulus used in non-linear
!												  elasticity (adimensionel)
!									 sm8c (4,i) = E-unld, unloading modulus used in
!												  non-linear elasticity (adim)
!									 sm8c (5,i) = n, exponent used to compute the Young
!												  loading and unloading moduli in non-linear
!												  elasticity (adim)
!									 sm8c (6,i) = Rf, coefficient multiplying the ratio
!												  (sigm1-sigm3)/(sigm1-sigm3)ult in nlin els
!												  (adim)
!                                    sm8c (7,i) = Kb, bulk modulus used in non-linear
!											      elasticity (adim)
!	                                 sm8c (8,i) = m, exponent used to compute the bulk
!												  modulus B in nlin els (adim)
!									 sm8c (9,i) = e-f, Young modulus minimal value (Pa)
!									 sm8c (10,i)=Ko, steady state coefficient (adim)
!								     sm8c (11,i)= gamma, soil volumetric weight (N,m-3)
!								     sm8c (12,i)= ten-st, traction resistency (Pa)
!									 sm8c (13,i)= perm-z
!									 sm8c (14,i)= sat
!									 sm8c (15,i)= compm
!									 sm8c (16,i)= kx/kz
!						  icpt:[3]-> sm8c (1,i) = friction angle
!									 sm8c (2,i) = cohesion
!								     sm8c (3,i) = OCR
!								     sm8c (4,i) = k
!                                    sm8c (5,i) = n, exponent used to compute the Young
!												  loading and unloading moduli in non-linear
!												  elasticity (adim)
!									 sm8c (6,i) = e
!									 sm8c (7,i) = Kb, bulk modulus used in non-linear
!												  elasticity (adim)
!									 sm8c (8,i) = m, exponent used to compute the bulk
!												  modulus B in nlin els (adim)
!									 sm8c (10,i) =Ko
!									 sm8d (11,i) =gamma
!									 sm8c (13,i) = perm-z
!									 sm8c (14,i) = sat
!									 sm8c (15,i) = compm
!									 sm8c (16,i) = kx/kz
!						  icpt:[4]-> ??????
!	INPUT		: ibem, ne8c, m8c, sm8c, m8cep
!	OUTPUT		: sm8c
!****************************************************************************************
!
	subroutine matec(m8c,sm8c,m8cep,ibem,ne8c)
	implicit double precision (a-h,o-z)
	character (5) :: card, wd(21)
!
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /ocqm/ EPScqm,taw,Icqm,Lcqm,NPcqm
!
	dimension sm8c(60,1),m8cep(1)
!
	common /ipt/ istop,list,wd
!
	if (ibem.ne.2 .and. ne8c.eq.0) then
		write (*,*) '**Error** in MATEC: THERE IS not SATURATED MATERIAL'
		istop=istop+1
		return
	endif
!
	do i=1,m8c
		do j=1,60
			sm8c(j,i)=0.d0
		enddo
	enddo
	do 70 i=1,m8c
		read (4,*,err=10) imat,icpt
		goto 20
10		write (*,*) '**Error** in reading in MATEC'
		stop
20		if (imat.lt.1.OR.imat.gt.m8c.OR.icpt.lt.1.OR.icpt.gt.4) then
			write (*,1) imat,icpt
			istop=istop+1
			goto 70
		endif
!
		m8cep(imat)=icpt
		goto (30,40,50,60) icpt
!
!	  linear elastic: icpt=[1]
30		if (kfsol.eq.0.and.Icqm.eq.1) then !QUASI-STATIC PROBLEM BY CQM
			read (4,*,err=10) sm8c(3,imat),sm8c(7,imat),sm8c(10,imat),sm8c(13,imat),&
							  sm8c(14,imat),sm8c(15,imat)
		else if (kfsol.eq.3.and.Icqm.eq.1) then !DYNAMIC PROBLEM BY CQM
			read (4,*,err=10) sm8c(3,imat),sm8c(7,imat),sm8c(9,imat),&
							  sm8c(10,imat),sm8c(11,imat),&
							  sm8c(12,imat),sm8c(13,imat),sm8c(15,imat)
		else
			read (4,*,err=10) sm8c(3,imat),sm8c(7,imat),sm8c(10,imat),sm8c(11,imat),&
							  sm8c(13,imat),sm8c(14,imat),sm8c(15,imat),sm8c(16,imat)
		endif
		goto 70
!
!	  hyperbolic non-linear elastic for statics: icpt=[2]
40		read (4,*,err=10) (sm8c(j,imat),j=1,8)
		sm8c(1,imat)=sm8c(1,imat)*1.74532925199432d-2 !??????????????
		read (4,*,err=10) (sm8c(j,imat),j=9,16)
		goto 70
!
!	  hyperbolic non-linear elastic for for dynamics: icpt=[3]
50		read (4,*,err=10) (sm8c(j,imat),j=1,8)
		sm8c(1,imat)=sm8c(1,imat)*1.74532925199432d-2 !??????????????
		read (4,*,err=10) sm8c(10,imat),sm8c(11,imat),(sm8c(j,imat),j=13,16)
		goto 70
!
!	  Prevost elasto-plastic 1985: icpt=[4]
60		read (4,*,err=10) (sm8c(j,imat),j=1,7)
		sm8c(1,imat)=sm8c(1,imat)*1.74532925199432d-2  !??????????????
		read (4,*,err=10) (sm8c(j,imat),j=9,16)
		nys=sm8c(9,imat)
		do k=1,nys
			j=4*k
			read (4,*,err=10) sm8c(17+j,imat),sm8c(18+j,imat),sm8c(19+j,imat),sm8c(20+j,imat)
		enddo
70	continue
!
1	FORMAT ('**Error** invalid materiel number',2i5,' in MATEC')
	return
	end
!
!
!****************************************************************************************
!							           MATEU subroutine
!	This subroutine reads the properties of unsaturated materials.
!	This subroutine is called in: INPUT
!	This subroutine calls: -
!	variables used are:
!		m8u		: number of unsaturated materials in FE & BE zones; in BE area the program
!				  is written for just one material (homogeneous space) then if we use an
!				  unsaturated material in BE zone (by considering the unsaturated materials
!				  in FE zone) the number of material of BE zone is always equal to M8U
!				  because the materials from "1" to "M8U-1" are used in FE zone
!		m8uep	: behaviour type of unsaturated materials [1,2,3,4];
!								[1]: linear,
!				  if we use an unsaturated material in BE zone then the behaviour type of BE
!				  zone's material is: M8UEP(m8u)
!		icpt	: behaviour type code [1,2,3,4]
!							[1]: linear
!							[2]: non-linear
!		imat	: number (num�ro) of unsaturated materials
!		sm8u	: all of the parameters corresponding to each behaviour of unsat. soil are
!				  defined in this matrix.
!					  icpt:[1]->
!						 sm8u(3,i) = Kl,loading modulus used in non-linear elasticity(adim.)(E)
!						 sm8u(5,i) = Ca, Compressibility of air   (=droa/(roa*dpa))
!						 sm8u(6,i) = Cw, Compressibility of water (=drow/(row*dpw))
!						 sm8u(7,i) = Kb,bulk modulus used in non-linear elasticity (adim.)(B)
!						 sm8u(8,i) = m, exponent used to compute the bulk modulus B in non-linear
!									 elasticity (adimensional)
!						 sm8u(11,i)= Young modulus minimal value (Pa)
!						 sm8u(14,i)= Sw0,initial saturation degree(adimensional)
!						 sm8u(15,i)= aw, constant (m.s-1) used in the formula of
!						        	 water permeability
!						 sm8u(16,i)= alpha_w, constant (adimensional) used in
!									 the formula of water permeability
!						 sm8u(17,i)= Swr,residual saturation degree (adim.),
!									 used in the formula of water permeability
!						 sm8u(19,i)= b,constant (m2) used in the formula of air permeability
!						 sm8u(20,i)= alpha_,constant (adim.) used in the formula of air
!									 permeability
!						 sm8u(21,i)= "viscosity",dynamic air viscosity(N.s.m-2),used in the
!									 formula of air permeability.
!									 Generally,Visc_a=1,846.10-5 N.s.m-2.
!						 sm8u(22,i)= Kw_max, maximal authorized water permeability value(m.s-1)
!						 sm8u(23,i)= Ka_max, maximal authorized air permeability value (m.s-1)
!						 sm8u(24,i)= ae, constant (adimensional) used in the formula of the
!									 void ratio state surface
!						 sm8u(25,i)= be, constant (adimensional) used in the formula of the
!									 void ratio state surface
!						 sm8u(27,i)= e0,initial void ration used in the formula of the
!									 void ratio state surface
!						 sm8u(28,i)= SIGMe,maximal traction resistance (Pa) used in the
!									 formula of the void ratio state surface
!						 sm8u(31,i)= beta_w,constant (1/Pa) used in the formula of
!									 the saturation degree state surface
!						 sm8u(33,i)= msuc-e
!						 sm8u(34,i)= Henry coefficient
!	INPUT		: ibem, ne8u, m8u, m8uep
!	OUTPUT		: sm8u
!****************************************************************************************
!
	subroutine mateu(m8u,sm8u,m8uep,ibem,ne8u,idyn)
	implicit double precision (a-h,o-z)
	character (5) :: card, wd(21)
!
	dimension sm8u(60,1),m8uep(1)
!
	common /ipt/ istop,list,wd
!
	if (ibem.ne.3 .and. ne8u.eq.0) then
		write (*,*) '**Error** in MATEU: THERE IS not UNSATURATED MATERIAL'
		istop=istop+1
		return
	endif
!
	do i=1,m8u
		do j=1,60
			sm8u(j,i)=0.d0
		enddo
	enddo
	do 70 i=1,m8u
		read (4,*,err=10) imat,icpt
		goto 20
10		write (*,*) '**Error** in reading in MATEU'
		stop
20		if (imat.lt.1.OR.imat.gt.m8u.OR.icpt.lt.1.OR.icpt.gt.4) then
			write (*,1) imat,icpt
			istop=istop+1
			goto 70
		endif
!
		m8uep(imat)=icpt
		goto (30,40,50,60) icpt
!
!	  linear elastic: icpt=[1]
30		if (idyn.eq.2) then
			read (4,*,err=10) sm8u(3,imat),sm8u(5,imat),sm8u(6,imat),sm8u(7,imat),&
				sm8u(8,imat),sm8u(11,imat),sm8u(14,imat),sm8u(15,imat),sm8u(16,imat),&
				sm8u(17,imat),sm8u(19,imat),sm8u(20,imat),sm8u(21,imat),&
				sm8u(22,imat),sm8u(23,imat),sm8u(24,imat),sm8u(25,imat),&
				sm8u(27,imat),sm8u(28,imat),sm8u(31,imat),sm8u(33,imat)
		else
			read (4,*,err=10) sm8u(3,imat),sm8u(7,imat),sm8u(8,imat),sm8u(11,imat),&
				sm8u(14,imat),sm8u(15,imat),sm8u(16,imat),sm8u(17,imat),&
				sm8u(19,imat),sm8u(20,imat),sm8u(21,imat),sm8u(22,imat),&
				sm8u(23,imat),sm8u(24,imat),sm8u(25,imat),sm8u(27,imat),&
				sm8u(28,imat),sm8u(31,imat),sm8u(33,imat),sm8u(34,imat)
		endif
		goto 70
!
!	  non-linear elastic: icpt=[2]
40		write (*,*) 'UNDER CONSTRUCTION'
		stop
		goto 70
!
50		write (*,*) 'UNDER CONSTRUCTION'
		stop
		goto 70
!
60		write (*,*) 'UNDER CONSTRUCTION'
		stop
		goto 70
!
70	continue
!
1	FORMAT ('**Error** invalid materiel number',2i5,' in MATEU')
!
	return
	end
!
!
!****************************************************************************************
!							           BOUND subroutine
!	This subroutine reads & generates the fixed degree of freedom.
!	This subroutine is called in: INPUT
!	This subroutine calls: -
!	variables used are:
!
!		nbcx	: number of (null) boundary conditions introduced for x-displacements
!		nbcy	: number of (null) boundary conditions introduced for y-displacements
!		nbcw	: number of (null) boundary conditions introduced for nodal water pressure
!		nbca	: number of (null) boundary conditions introduced for nodal air pressure
!		nbc		: total number of null boundary conditions; nbc = nbcx+nbcy+nbcw+nbca
!		kodex	: numbers of nodes in which zero displacement in X-dir is imposed
!		kodey	: numbers of nodes in which zero displacement in Y-dir is imposed
!		kodew	: numbers of nodes in which zero pore water pressure is imposed
!		kodea	: numbers of nodes in which zero pore gas pressure is imposed
!		kode	: the code of boundary conditions [1,2,3,4];
!							[1]: b.c for X-direction (Ux=0)
!							[2]: b.c for Y-direction (Uy=0)
!							[3]: b.c for Pw-direction (Px=0)
!							[4]: b.c for Pg-direction (Pg=0)
!	INPUT		: nnp, nbcx, nbcy, nbcw, nbca
!	OUTPUT		: kodex, kodey, kodew, kodea
!****************************************************************************************
!
	subroutine bound(nnp,nbcx,nbcy,nbcw,nbca,kodex,kodey,kodew,kodea)
	implicit double precision (a-h,o-z)
	character*5 card,wd(21)
!
	dimension kodex(1),kodey(1),kodew(1),kodea(1)
!
	common /ipt/ istop,list,wd
!
	ix=0
	iy=0
	iw=0
	ia=0
!
10	read (4,*) card
	if (card(1:1).eq.'#') goto 10
	do i=1,list
		if (card.eq.wd(i)) goto 90
	enddo
	BACKSPACE (4)
!
20	read (4,*,err=30) n1,n2,lg,kode
	goto 40
!
30	write (*,*) '**Error** in reading in BOUND'
	stop
!
!     verification of data
40	if ((n1.gt.nnp).OR.(n1.le.0).OR.(n2.gt.nnp).OR.(n2.le.0)) then
		write (*,1) n1,n2
		istop = istop +1
		goto 10
	endif
	if ((kode.gt.4).OR.(kode.le.0)) then
		write (*,3) kode
		istop = istop +1
		goto 10
	endif
!
	lg = SIGN (lg,n2-n1)
!	lg=sign(max(abs(lg),1),n2-n1)
	l= n1
!
60	if ((n2-l)*lg .lt. 0) goto 10
	if (l.le.0 .OR. l.gt.nnp) then
		write (*,2) l
		istop = istop + 1
		goto 10
	endif
!
!     read & generate b.c. data for x-dir
	if (kode.eq.1) then
		ix=ix+1
		if (ix.gt.nbcx) goto 80
		kodex(ix)=l
	endif
!
!     read & generate b.c. data for y-dir
	if (kode.eq.2) then
		iy=iy+1
		if (iy.gt.nbcy) goto 80
		kodey(iy)=l
	endif
!
!     read & generate b.c. data for Pw-dir
	if (kode.eq.3) then
		iw=iw+1
		if (iw.gt.nbcw) goto 80
		kodew(iw)=l
	endif
!
!     read & generate b.c. data for Pg-dir
	if (kode.eq.4) then
		ia=ia+1
		if (ia.gt.nbca) goto 80
		kodea(ia)=l
	endif
!
	if (lg.eq.0) goto 10
	l=l+lg
	goto 60
!
80	write (*,4)
	istop=istop+1
90	BACKSPACE (4)
!
	if (ix.lt.nbcx .OR. iy.lt.nbcy .OR. iw.lt.nbcw .OR. ia.lt.nbca) then
		write (*,*) '   **Error** number of b.c is not sufficient'
		istop=istop+1
		return
	endif
!
1	FORMAT (5x,'**Error** invalid node number:',i5,' in BOUND')
2	FORMAT (5x,'**Error** attemp to generate node:',i5,' in BOUND')
3	FORMAT (5x,'**Error** invalid b.c code',i5,' in BOUND')
4	FORMAT (5x,'**Error** attemp to number of b.c in BOUND')
!
	return
	end
!
!
!****************************************************************************************
!			                            PRTOUT subroutine
!	This subroutine reads number of nodes which we want to print their results.
!	This subroutine is called in: INPUT
!	This subroutine calls: -
!	variables used are:
!		init	: code describing the type of initial conditions [0,1,2,3]. when we want to
!				  do a topographical effect (ifem=0) we continue by INIT=[0]
!							[0] if the results are taken from a previous analysis
!							[1] if initial normal stresses are considered equal to the
!								weight of the upper soil layers
!							[2] if initial stresses are computed for the steady state of
!							    the system
!							[3] if initial stresses are calculated as the weight of the
!								upper soil layers, considering that the system has a mere
!							    parallelogram shape, and that the x-direction of the
!							    space-frame corresponds to the physical horizontal
!		iprint	: initial data output print code (control code for the input and output files)
!				  [0,1];
!							[0] to avoid the outprint of the mesh file
!							[1] to print out the mesh file
!		iterp	: number of iloads for which the results are printed out in the output file
!		itprt	: iload's numbers which are selected to print out the results;
!							[0] if no printout is expected
!		nnpo	: number of nodes that we need their outputs, it must be at maximum 80 nodes.
!		ne8do	: number of drained elements that we need their outputs, it must be at maximum
!				  10 elms.
!		ne8co	: number of saturated elements that we need their outputs
!		ne8uo	: number of unsaturated elements that we need their outputs
!		inpout	: nodes' numbers which we need their outputs, INPOUT (1:nnpo)
!		ie8dout	: numbers of the nodes forming the current drained element that we need
!				  their outputs, IE8doUT (1:NE8do)
!		ie8cout	: numbers of the nodes forming the current saturated element that we need
!				  their outputs, IE8COUT (1:NE8CO)
!		ie8uout	: numbers of the nodes forming the current unsaturated element that we need
!				  their outputs, IE8UOUT (1:NE8UO)
!	INPUT		:
!	OUTPUT		: init, iprint, iterp, nnpo, ne8do, ne8co, ne8uo, itprt, inpout, ie8dout, ie8cout, ie8uout
!****************************************************************************************
!
	subroutine prtout
	implicit double precision (a-h,o-z)
!
	common /g/init,iprint
	common /i/iterp,itprt(100),ne8do,ne8co,ne8uo,ie8dout(10),ie8cout(10),ie8uout(10),nnpo,&
			  inpout(100)
!
!
	read (4,*,err=40) init,iprint,iterp,nnpo,ne8do,ne8co,ne8uo
	read (4,*) (itprt(i),i=1,iterp)
!
	if (nnpo.eq.0) goto 10
	read (4,*,err=40) (inpout(i),i=1,nnpo)
!
10	if (ne8do.eq.0) goto 20
	read (4,*,err=40) (ie8dout(i),i=1,ne8do)
!
20	if (ne8co.eq.0) goto 30
	read (4,*,err=40) (ie8cout(i),i=1,ne8co)
!
30	if (ne8uo.eq.0) goto 50
	read (4,*,err=40) (ie8uout(i),i=1,ne8uo)
!
	goto 50
40	write (*,*) '**Error** in reading in PRINT'
	stop
50	return
	end
!
!
!****************************************************************************************
!			                            PHASE subroutine
!	This subroutine reads the construction phases.
!	This subroutine is called in: INPUT
!	This subroutine calls: -
!	variables used are:
!		nload	: number of load steps; load is applied in "nload" steps, in each
!				  "iload=1:nload" the characteristics of loading such as 'intensity-nstep-...
!				  are differents
!		iconst	: iconst(1,nload1)= Initial load step number (time step)
!				  iconst(2,nload1)= Number of load steps (nstep)
!				  SUM[iconst(2,1:nload)]=ntime
!				  iconst(3,nload1)= Maximum iterations per load step (itmax)
!	INPUT		: nload
!	OUTPUT		: iconst
!****************************************************************************************
!
	subroutine phase(nload,iconst)
	integer i,nload,iconst
!
	dimension iconst(3,1)
!
	do i=1,nload
		read (4,*,err=10) (iconst(j,i),j=1,3)
	enddo
!
	goto 20
!
10	write (*,*) '**Error** in reading in PHASE'
	stop
20  return
	end
!
!
!****************************************************************************************
!			                            INITI subroutine
!	This subroutine reads the initial values.
!	This subroutine is called in: INPUT
!	This subroutine calls: -
!	variables used are:
!		idyn	: dynamic code or type of loading [0,1,2];
!								[0]: static
!								[1]: quasi-static
!								[2]: dynamic problem
!		npwi	: number of nodes on which an initial water pressure is imposed which can
!				  be at max. 400 nodes
!		npai	: number of nodes on which an initial gas pressure is imposed which can be
!				  at max. 400 nodes
!		ndispi	: number of nodes on which an initial displacement is imposed which can be
!				  at max. 400 nodes
!		nveli	: number of nodes on with an initial velocity is imposed which can be at
!				  max. 400 nodes
!		nacci	: number of nodes on with an initial acceleration is imposed which can be
!				  at max. 400 nodes
!		nnpwi	: nodes' numbers on which an initial water-p is imposed
!		vnpwi	: initial water pressure values imposed on the preceding nodes
!		nnpai	: nodes' numbers on which an initial air-p is imposed
!		vnpai	: initial air pressure values imposed on the preceding nodes
!		nndispi	: nodes' numbers on which an initial disp. is imposed
!		vxdispi	: initial displacement values imposed on the preceding nodes (X-dir)
!		vydispi	: initial displacement values imposed on the preceding nodes (Y-dir)
!		nnveli	: nodes' numbers on which an initial vel. is imposed
!		vxveli	: initial velocity values imposed on the preceding nodes (X-dir)
!		vyveli	: initial velocity values imposed on the preceding nodes (Y-dir)
!		nnacci	: nodes' numbers on which an initial acc. is imposed
!		vxacci	: initial acc. values imposed on the preceding nodes (X-dir)
!		vyacci	: initial acc. values imposed on the preceding nodes (Y-dir)
!	INPUT		: idyn
!	OUTPUT		: npwi,nnpwi,vnpwi; npai,nnpai,vnpai; ndispi,nndispi,vxdispi,vydispi;
!				  nveli,nnveli,vxveli,vyveli; nacci,nnacci,vxacci,vyacci
!****************************************************************************************
!
	subroutine initi(idyn,ifem,ibem)
	implicit double precision (a-h,o-z)
!
	common /ini1/ vnpwi(400),vnpai(400),nnpwi(400),nnpai(400)
	common /ini2/ npwi,npai,ndispi,nveli,nacci
	common /ini3/ nndispi(400),nnveli(400),nnacci(400)
	common /ini4/ vxdispi(400),vydispi(400)
	common /ini5/ vxveli(400),vyveli(400)
	common /ini6/ vxacci(400),vyacci(400)
!
!	if (idyn.eq.0) then
!		write (*,*) '**Error** STATIC PROBLEM HAVE not INITIAL VALUES'
!		istop=istop+1
!		return
!	endif
!
	read (4,*,err=50) nacci,nveli,ndispi,npwi,npai
!
	if (npai.gt.400 .OR. npwi.gt.400 .OR. ndispi.gt.400 .OR. nveli.gt.400 .OR.&
	    nacci.gt.400) then
			write (*,*) '**Error** number of initial values is greater than 400'
			istop=istop+1
			return
	endif
!
	if (nacci.eq.0) goto 10
	read (4,*,err=50) (nnacci(i),i=1,nacci)
	read (4,*,err=50) (vxacci(i),i=1,nacci)
	read (4,*,err=50) (vyacci(i),i=1,nacci)
!
10	if (nveli.eq.0) goto 20
	read (4,*,err=50) (nnveli(i),i=1,nveli)
	read (4,*,err=50) (vxveli(i),i=1,nveli)
	read (4,*,err=50) (vyveli(i),i=1,nveli)
!
20	if (ndispi.eq.0) goto 30
	read (4,*,err=50) (nndispi(i),i=1,ndispi)
	read (4,*,err=50) (vxdispi(i),i=1,ndispi)
	read (4,*,err=50) (vydispi(i),i=1,ndispi)
!
30	if (npwi.eq.0) goto 40
	read (4,*,err=50) (nnpwi(i),i=1,npwi)
	read (4,*,err=50) (vnpwi(i),i=1,npwi)
!
40	if (npai.eq.0) goto 60
	read (4,*,err=50) (nnpai(i),i=1,npai)
	read (4,*,err=50) (vnpai(i),i=1,npai)
!
	goto 60
!
50	write (*,*) '**Error** in reading in INITI'
	stop
1	FORMAT (//2x,'number of initial values is greater than 400')
!
60	return
	end
!
!
!****************************************************************************************
!			                            MAILLAGE subroutine
!	This subroutine: maillage et resultt a lire par Gid
!****************************************************************************************
!
	subroutine maillage(x,ie8d,ie8c,ie8u,ie3d1,ie3d2,ie3d3,infile)
!
	implicit double precision (a-h,o-z)
	character*72 mesh,infile,res
!
	common /a/ nnp,ifem,ibem,idyn,ieaq,nr
	common /a1/ ne8d,ne8c,ne8u
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /a3/ ne8d1,ne8c1,ne8u1,ne3u1,nload1,nbcx1,nbcy1,nbcw1,nbca1,m8d1,m8c1,m8u1,&
				nbed11,nbed21,nbed31,nbeint1,ne8dd,ne8cc,ne8uu
	common /a4/ xencl(500),yencl(500)
	dimension x(2,nnp),ie8d(9,ne8d1),ie8c(9,ne8c1),ie8u(9,ne8u1),ie3d1(npbmax,nbed11),&
			  ie3d2(npbmax,nbed21),ie3d3(npbmax)
!
	l=LEN_TRIM(infile)-4
	mesh=infile(1:l)//'.flavia.msh'
	res=infile(1:l)//'.flavia.res'
	open (8,file=mesh)
	open (9,file=res)
	write (9,*) 'GiD Post Results File 1.0'
!
	ielm=0
	if (ifem.eq.0) goto 30
!	elements finis secs
	if (ne8d.eq.0) goto 10
	write (8,3)
	do i=1,nnp
		write (8,1) i,x(1,i),x(2,i)
	enddo
	write (8,*) 'end coordinates'
	write (8,*) 'Elements'
	do i=1,ne8d
		write (8,2) i+ielm,ie8d(1,i),ie8d(3,i),ie8d(5,i),ie8d(7,i),&
                    ie8d(2,i),ie8d(4,i),ie8d(6,i),ie8d(8,i),ie8d(9,i)
	enddo
	write (8,*) 'end elements'
	ielm=ielm+ne8d
!
!     elements finis saturated
10	if (ne8c.eq.0) goto 30
	write (8,4)
	if (ne8c.eq.0) then
		do i=1,nnp
			write (8,1) i,x(1,i),x(2,i)
		enddo
	endif
	write (8,*) 'end coordinates'
	write (8,*) 'Elements'
	do i=1,ne8c
		write (8,2) i+ielm,ie8c(1,i),ie8c(3,i),ie8c(5,i),ie8c(7,i),&
                    ie8c(2,i),ie8c(4,i),ie8c(6,i),ie8c(8,i),ie8c(9,i)
	enddo
	ielm=ielm+ne8c
	write (8,*) 'end elements'
!
!     �l�ments fronti�res zone 1
30	if (nbed1.eq.0) goto 50
	write (8,5)
	if (ifem.ne.0) goto 40
	do i=1,nnp
		write (8,1) i,x(1,i),x(2,i)
	enddo
	write (8,*) 'end coordinates'
40	continue
	write (8,*) 'Elements'
	n=nnpbed1(1)/2
	do i=1,n
		n1=ie3d1(2*i-1,1)
		n2=ie3d1(2*i+1,1)
		n3=ie3d1(2*i,1)
		if (i.eq.n) n2=ie3d1(1,1)
		write (8,9) i+ielm,n1,n2,n3
	enddo
	ielm=ielm+n
	write (8,*) 'end elements'
!
!	�l�ments fronti�res zone 2
50	if (nbed2.eq.0) goto 70
	write (8,6)
	if (ifem.ne.0 .or. nbed1.ne.0) goto 60
	do i=1,nnp
		write (8,1) i,x(1,i),x(2,i)
	enddo
	write (8,*) 'end coordinates'
60	continue
	write (8,*) 'Elements'
	n=nnpbed2(1)/2
	do i=1,n
		n1=ie3d2(2*i-1,1)
		n2=ie3d2(2*i+1,1)
		n3=ie3d2(2*i,1)
		if (i.eq.n) n2=ie3d2(1,1)
		write (8,9) i+ielm,n1,n2,n3
	enddo
	ielm=ielm+n
	write (8,*) 'end elements'
!
!	�l�ments fronti�res zone 3 and enclosing element
70	if (nbed3.eq.0) goto 90
	write (8,7)
	if (ifem.ne.0 .OR. nbed1.ne.0 .OR. nbed2.ne.0) goto 80
	do i=1,nnp
		write (8,1) i,x(1,i),x(2,i)
	enddo
80	do i=nnp+1,nnp+nnpen
		write (8,1) i,xencl(i-nnp),yencl(i-nnp)
	enddo
	write (8,*) 'end coordinates'
	write (8,*) 'Elements'
	n=nnpbed3/2
	do i=1,n
		n1=ie3d3(2*i-1)
		n2=ie3d3(2*i+1)
		n3=ie3d3(2*i)
		write (8,9) i+ielm,n1,n2,n3
	enddo
	ielm=ielm+n+1
	write (8,9) ielm,ie3d3(nnpbed3),nnp+2,nnp+1
	nn=nnpen/2
	do i=1,nn
		n1=nnp+2*i
		n2=nnp+2*i+2
		n3=nnp+2*i+1
		if (i.eq.nn) n2=ie3d3(1)
		write (8,9) i+ielm,n1,n2,n3
	enddo
	write (8,*) 'end elements'
!
90	close (8)
!
1	FORMAT (i5,4x,e15.4,3x,e15.4)
2	FORMAT (10i5)
3	FORMAT ('MESH "FEM dry element" dimension 2 ElemType ',&
			'Quadrilateral Nnode 8'/'Coordinates')
4	FORMAT ('MESH "FEM satured element" dimension 2 ElemType ',&
			'Quadrilateral Nnode 8'/'Coordinates')
5	FORMAT ('MESH "BEM element, zone 1" dimension 2 ElemType ',&
			'Linear Nnode 3'/'Coordinates')
6	FORMAT ('MESH "BEM element, zone 2" dimension 2 ElemType ',&
			'Linear Nnode 3'/'Coordinates')
7	FORMAT ('MESH "BEM element, zone 3 and enclosing element" ',&
			'dimension 2 ElemType Linear Nnode 3'/'Coordinates')
9	FORMAT (4i5)
!
	return
	end
!
!
!****************************************************************************************
!			                            OUTGID subroutine
!	This subroutine: resultat en deplacement � lire par Gid
!****************************************************************************************
!
	subroutine outgid(itime,nnp,mdofn,iact,id,dt)
!
	implicit double precision (a-h,o-z)
	dimension iact(1),id(3,1),dt(1)
!
	write (9,*) 'Result "Displacment" "Load" ',itime,' Vector OnNodes'
	write (9,*) 'Values'
	do 10 i=1,nnp
		j=iact(i)
        if (j.eq.0) goto 10
		i1=id(1,i)
		i2=id(2,i)
		u1=0.0d0
		u2=0.0d0
		u3=0.d0
		if (i1.le.mdofn) u1=dt(i1)
		if (i2.le.mdofn) u2=dt(i2)
		write (9,20) i,u1,u2,u3
10		continue
	write (9,*) 'End Values'
!
20	FORMAT (i5,3e15.4)
!
	return
	end
!****************************************************************************************
!
	subroutine intdisp(id,x,xint,sm8d,sm8c,sm8u,ie3d1,ie3d2,ie3d3,gihd1,gihd2,gihd3,gid1,&
					   gid2,gid3,rbed1,rbed2,rbed3,disp,sig3u,tbed1,tbed2,tbed3,ubed1,ubed2,&
					   ubed3)
!
	implicit double precision (a-h,o-z)
!
	common /a/ nnp,ifem,ibem,idyn,ieaq,nr
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /a3/ ne8d1,ne8c1,ne8u1,ne3u1,nload1,nbcx1,nbcy1,nbcw1,nbca1,m8d1,m8c1,m8u1,&
				nbed11,nbed21,nbed31,nbeint1,ne8dd,ne8cc,ne8uu
	common /b/ m8d,m8c,m8u,m8dep(10),m8cep(10),m8uep(10)
	common /c/ nload,ntime,itmax,itime,iter,iload,time,dtime,dtimei,dtimem,xinc
	common /h/ mdof,mdofn,ls,ls1,isolv
	common /i/iterp,itprt(100),ne8do,ne8co,ne8uo,ie8dout(10),ie8cout(10),ie8uout(10),nnpo,&
			  inpout(100)
!
	dimension id(4,nnp),x(2,nnp),xint(2,nbeint),sm8d(60,m8d1),sm8c(60,m8c1),sm8u(60,m8u1),&
			  ie3d1(npbmax,nbed11),ie3d2(npbmax,nbed21),ie3d3(npbmax),sig3u(5,500,nbed11),&
			  disp(mdofn),gihd1(n12gh,n12gh,nbed11),gihd2(n12gh,n12gh,nbed21),&
			  gihd3(n12gh,n12gh),gid1(n12gh,n12gh,nbed11),gid2(n12gh,n12gh,nbed21),&
			  gid3(n12gh,n12gh),tbed1(ntime,n12gh,nbed11),tbed2(ntime,n12gh,nbed21),&
			  tbed3(ntime,n12gh),ubed1(ntime,n12gh,nbed11),ubed2(ntime,n12gh,nbed21),&
			  ubed3(ntime,n12gh),rbed1(ntime,n12gh,nbed11),rbed2(ntime,n12gh,nbed21),&
			  rbed3(ntime,n12gh),uint1(kbem*nbeint,nbed11),uint2(kbem*nbeint,nbed21),&
			  uint3(kbem*nbeint),ubdn1(n12gh,nbed11),tbdn1(n12gh,nbed11)
!
!
!	---------------------- BEM area D1 ----------------------
!
	if (nbed1.eq.0) goto 40
	call intdispd1(id,x,xint,sm8d,sm8c,sm8u,ie3d1,tbed1,ubed1,rbed1,gihd1,gid1,disp,sig3u,&
				   uint1)
	write (*,*)'			intdispd1 concluded'
40	continue
!
!	---------------------- BEM area D2 ----------------------
!
	if (nbed2.eq.0) goto 60
!	call intdispd2(id,x,xint,sm8d,sm8c,sm8u,ie3d2,rbed2,gihd2,gid2,disp,sig3u,uint2)
	write (*,*)'			intdispd2 concluded'
60	continue
!
!	---------------------- BEM area D3 ----------------------
!
	if (nbed3.eq.0) goto 80
!	call intdispd3(id,x,xint,sm8d,sm8c,sm8u,ie3d3,rbed3,gihd3,gid3,disp,sig3u,uint3)
	write (*,*)'			intdispd3 concluded'
80	continue
!
!
90	continue
!
240	FORMAT (i5,4x,e15.4,3x,e15.4)
281	FORMAT (//10x,'Internal Nodes number in BEM domain'/)
!
	return
	end



!****************************************************************************************
!
	subroutine intdispd1(id,x,xint,sm8d,sm8c,sm8u,ie3d1,tbed1,ubed1,rbed1,gihd1,gid1,disp,&
						 sig3u,uint1)
!
	implicit double precision (a-h,o-z)
!
	common /a/ nnp,ifem,ibem,idyn,ieaq,nr
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /a3/ ne8d1,ne8c1,ne8u1,ne3u1,nload1,nbcx1,nbcy1,nbcw1,nbca1,m8d1,m8c1,m8u1,&
				nbed11,nbed21,nbed31,nbeint1,ne8dd,ne8cc,ne8uu
	common /b/ m8d,m8c,m8u,m8dep(10),m8cep(10),m8uep(10)
	common /c/ nload,ntime,itmax,itime,iter,iload,time,dtime,dtimei,dtimem,xinc
	common /h/ mdof,mdofn,ls,ls1,isolv
	common /i/iterp,itprt(100),ne8do,ne8co,ne8uo,ie8dout(10),ie8cout(10),ie8uout(10),nnpo,&
			  inpout(100)
	common /stability/anew1,anew2,wil,intb
!
	dimension disp(mdofn),id(4,nnp),x(2,nnp),xint(2,nbeint),sm8d(60,m8d1),sm8c(60,m8c1),&
			  sm8u(60,m8u1),ie3d1(npbmax,nbed11),ubed1(ntime,n12gh,nbed11),&
			  tbed1(ntime,n12gh,nbed11),rbed1(ntime,n12gh,nbed11),uint1(kbem*nbeint,nbed11),&
			  gihd1(n12gh,n12gh,nbed11),gid1(n12gh,n12gh,nbed11),sig3u(5,500,nbed11),&
			  rbdn1(kbem*nbeint,nbed11),hd1int(kbem*nbeint,n12gh),gd1int(kbem*nbeint,n12gh),&
			  ubdn1(n12gh,nbed11),tbdn1(n12gh,nbed11)
!
!
	ndim=kbem*nbeint
!
	do k=1,nbed1
		do i=1,ndim
			uint1(i,k)=0.d0
		enddo
	enddo
!
	do 1000 k1=1,nbed1
!
		n=nnpbed1(k1)*kbem
!
!	******** calculate displacement & traction vector of time step itime (i.e. n) ********
!
!		displacement/pressure vector ubed1 at time: itime
!
		do k=1,nnpbed1(k1)
			j=ie3d1(k,k1)
			do i=1,kbem
				ubdn1(kbem*k-kbem+i,k1)=disp(id(i,j))
			enddo
		enddo
!
!		traction/flow vector tbed1 at time: itime
!
		do i=1,n
			tbdn1(i,k1)=0.d0
			do j=1,n
				tbdn1(i,k1) = tbdn1(i,k1)+gihd1(i,j,k1)*ubdn1(j,k1)-&
										  gid1(i,j,k1) *rbed1(itime,j,k1)
			enddo
		enddo
!
!	**************************************************************************************
!
!	initializing vector rbdn
!
		do i=1,ndim
			rbdn1(i,k1)=0.d0
		enddo
!
!	this loop is dedicated to calculate Z(N)=SUM(n=1:N-1) [G(N+1-n)*T(n)-H(N+1-n)*U(n)]
!
10		do 500 itm=1,itime-1

!	calculating G & H matices at itime > 1

			call ghmatd1INT(itm,k1,x,xint,sm8d,sm8c,sm8u,ie3d1,hd1int,gd1int,sig3u)

			do i=1,ndim
				do j=1,n
					rbdn1(i,k1)=rbdn1(i,k1)-( gd1int(i,j)*tbed1(itm,j,k1)-&
										      hd1int(i,j)*ubed1(itm,j,k1) )/wil
				enddo
			enddo
!
500		continue
!
!	----- last moment itime -----
!
		itm=itime
		call ghmatd1INT(itm,k1,x,xint,sm8d,sm8c,sm8u,ie3d1,hd1int,gd1int,sig3u)
		do i=1,ndim
			do j=1,n
				rbdn1(i,k1)= rbdn1(i,k1)-( gd1int(i,j)*tbdn1(j,k1)-&
								                       hd1int(i,j)*ubdn1(j,k1) )
			enddo
		enddo
		if (itime.eq.1) GOTO 100

		itm=itime-1
		do i=1,ndim
			do j=1,n
				rbdn1(i,k1)= rbdn1(i,k1)-( gd1int(i,j)*tbed1(itm,j,k1)-&
								        hd1int(i,j)*ubed1(itm,j,k1) )*(1.d0-wil)/wil
			enddo
		enddo
!
100		continue
!
		do i=1,ndim
			uint1(i,k1)=rbdn1(i,k1)
		enddo
!
		j=9
		if (nnpo.eq.0) goto 8
		do i=1,nnpo
			j=j+i
		enddo
		j=j+1
8		do i=1,nbeint
			j=j+i
			if (ibem.eq.1) then
				write (j,15) time,uint1(i*kbem-1,k1),uint1(i*kbem,k1)
			else if (ibem.eq.2) then
				write (j,15) time,uint1(i*kbem-2,k1),uint1(i*kbem-1,k1),&
								   uint1(i*kbem,k1)
			else
				write (j,15) time,uint1(i*kbem-3,k1),uint1(i*kbem-2,k1),uint1(i*kbem-1,k1),&
							                          uint1(i*kbem,k1)
			endif
		enddo
!
1000 continue
!
15	FORMAT (10(1x,e15.6))
!
	return
	end


!****************************************************************************************
!
!							           MMATRICE subroutine?????? KODE NABAYAD = 2 bashad?
!				chon mikhahim F ra tolid konim pas bayad Ux & Uy ra dar nazar begirim???
!
!
!	This subroutine calculates local matrice Xmm in order to transform (line) element nodal
!	tractions to nodal forces.
!	This subroutine is called in: MTRFRD
!
!	variables used are:
!
!		xmm			: M (distributed matric) of each elm which is equal to
!					  Me=[N11,N12,N13;N21,N22,N23;N31,N32,N33] where
!					  Nlj=SUM[1:gi]{ Nl*Nj*jacobian*W] }
!					  this SUM for curved form element is obtained by 10 Gauss points.
!					  but for linear form element ????????
!		xelm(1:3)	: x-component of nodes forming the each element in one D3
!					  (semi-infinite) zone. in GHMATD3 subroutine in the first party to
!					  calculate the diagonal composents of fundamnetal solution with rigid
!					  body motion we consider NNPBED3 & NNPEN. but after obtaining final GD3
!					  & HD3 we consider only NNBED3 nodes
!		yelm(1:3)	: y-component of nodes forming the each element in one D3
!					  (semi-infinite) zone. in GHMATD3 subroutine in the first party to
!					  calculate the diagonal composents of fundamnetal solution with rigid
!					  body motion we consider NNPBED3 & NNPEN. but after obtaining final GD3
!					  & HD3 we consider only NNBED3 nodes
!		kode		: when this subroutine is called from MTRFRD it is "kbem" number of
!					  degree of freedom per node in BE zone;
!								if IBEM=[1] -> KBEM=[2] [X,Y]
!								if IBEM=[2] -> KBEM=[3] [X,Y,Pw]
!								if IBEM=[3] -> KBEM=[4] [X,Y,Pw,Pa]
!					  but when this subroutine is called from .... it is "2" ....
!		gi			: gauss points
!		wi			: weights of gauss points
!		e			: ???????????????
!		f			: shape functions
!		a			: coefficient in tangential vector to 1D element (Vxsi) in x-direction
!					  Vxsi(x)=(a) xsi + b
!		b			: coefficient in tangential vector to 1D element (Vxsi) in x-direction
!					  Vxsi(x)=(a) xsi + b
!		c			: coefficient in tangential vector to 1D element (Vxsi) in y-direction
!					  Vxsi(y)=(c) xsi + d
!		d			: coefficient in tangential vector to 1D element (Vxsi) in y-direction
!					  Vxsi(y)=(c) xsi + d
!
!
!	INPUT	:
!
!	OUTPUT	:
!
!
!****************************************************************************************
!
	subroutine mmatrice(xmm,xelm,yelm,kode)
!
	implicit double precision (a-h,o-z)
!
	dimension gi(10),wi(10),xelm(3),yelm(3),e(3,3),f(3),xmm(12,12)
!
	data gi/0.9739065285,-0.9739065285,&
			0.8650633666,-0.8650633666,&
			0.6794095683,-0.6794095683,&
			0.4333953941,-0.4333953941,&
			0.1488743389,-0.1488743389/
	data wi/0.0666713443,0.0666733443,&
			0.1494513491,0.1494513491,&
			0.2190863625,0.2190863625,&
			0.2692667193,0.2692667193,&
			0.2955242247,0.2955242247/
!
!
!	-------------------------- initializing --------------------------
!
	n1=3*kode
	do i=1,n1
		do j=1,n1
         xmm(i,j)=0.d0
		enddo
	enddo
!
!	if (xelm(2).eq.(xelm(1)+xelm(3))/2 .and. yelm(2).eq.(yelm(1)+yelm(3))/2) goto 200
!
!	geometrical parameters
	a=xelm(3)-2.d0*xelm(2)+xelm(1)
	b=(xelm(3)-xelm(1))/2.d0
	c=yelm(3)-2.d0*yelm(2)+yelm(1)
	d=(yelm(3)-yelm(1))/2.d0
!
!	evaluation of xmm matrice ,curved element
	do i=1,10
		 f(1)=gi(i)*(gi(i)-1)*0.5d0
         f(2)=1.d0 - gi(i)**2
         f(3)=gi(i)*(gi(i)+1)*0.5d0
         coef=dsqrt ((gi(i)*a+b)**2+(gi(i)*c+d)**2)*wi(i)
		 do ii=1,3
			do j=1,3
				do k=1,kode
					xmm((ii-1)*kode+k,(j-1)*kode+k)=&
					xmm((ii-1)*kode+k,(j-1)*kode+k)+f(ii)*f(j)*coef
				enddo
			enddo
		enddo
	enddo
	goto 300
!
!	evaluation of xmm matrice ,linear element ???????????????????
!
200	demi=dsqrt ((xelm(3)-xelm(1))**2+(yelm(3)-yelm(1))**2)/2
	e(1,1)=demi*4.d0/15.d0
	e(1,2)=demi*2.d0/15.d0
	e(1,3)=-demi/15.d0
	e(2,1)=e(1,2)
	e(2,2)=demi*16.d0/15.d0
	e(2,3)=demi*2.d0/15.d0
	e(3,1)=e(1,3)
	e(3,2)=e(2,3)
	e(3,3)=e(1,1)
	do i=1,3
		do j=1,3
			do k=1,kode
				xmm((i-1)*kode+k,(j-1)*kode+k)=e(i,j)
			enddo
		enddo
	enddo
!
300	return
	end

!****************************************************************************************
!
!							           MTRFRD subroutine
!
!	This subroutine calculates the global M matrice in order to transform nodal tractions
!	to nodal forces.
!	This subroutine is called in: GHMATD3
!	This subroutine calls:		  MMATRICE
!
!	variables used are:
!
!		n			: total number of nodes of one D3 (semi-infinite) area NNPBED3
!		xmtf		: M (distributed matric) of all of elements of one D3 area
!		xmm			:
!		xd			: x-component of nodes forming the D3 zone from line 1:nnpbed3
!					  & nodes forming the enclosed zone from line nnpbed3+1:nnpen
!		yd			: y-component of nodes forming the D3 zone from line 1:nnpbed3
!					  & nodes forming the enclosed zone from line nnpbed3+1:nnpen
!		kbem		: number of degree of freedom per node in BE zone;
!								if IBEM=[1] -> KBEM=[2] [X,Y]
!								if IBEM=[2] -> KBEM=[3] [X,Y,Pw]
!								if IBEM=[3] -> KBEM=[4] [X,Y,Pw,Pa]
!		xelm(1:3)	: x-component of nodes forming the each element in one D3
!					  (semi-infinite) zone. in GHMATD3 subroutine in the first party to
!					  calculate the diagonal composents of fundamnetal solution with rigid
!					  body motion we consider NNPBED3 & NNPEN. but after obtaining final GD3
!					  & HD3 we consider only NNBED3 nodes
!		yelm(1:3)	: y-component of nodes forming the each element in one D3
!					  (semi-infinite) zone. in GHMATD3 subroutine in the first party to
!					  calculate the diagonal composents of fundamnetal solution with rigid
!					  body motion we consider NNPBED3 & NNPEN. but after obtaining final GD3
!					  & HD3 we consider only NNBED3 nodes
!
!
!	INPUT	:
!
!	OUTPUT	:
!
!
!****************************************************************************************
!
	subroutine mtrfrd(n,xmtf,xd,yd)
!
	implicit double precision (a-h,o-z)
!
	common /a/ nnp,ifem,ibem,idyn,ieaq,nr
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
!
	dimension xd(npbmax),yd(npbmax),xelm(3),yelm(3),&
			  xmm(12,12),xmtf(n12gh,n12gh)
!
!
!	-------------------------- initializing --------------------------
!
	n1=kbem*n
	do i=1,n1
		do j=1,n1
			xmtf(i,j)=0.d0
		enddo
	enddo
!
!	------------------- evaluating the xmtf matrice -------------------
!
!	geometrical parameters
	do 100 i=1,n-1,2
		xelm(1)=xd(i)
		yelm(1)=yd(i)
		xelm(2)=xd(i+1)
		yelm(2)=yd(i+1)
		xelm(3)=xd(i+2)
		yelm(3)=yd(i+2)
!
!	call subroutine in order to evaluate elements participation in xmtf
!
		call mmatrice(xmm,xelm,yelm,kbem)
!
!	assemble
!
		nn=3*kbem
		do k=1,nn
			do l=1,nn
				ik=kbem*(i-1)+k
				il=kbem*(i-1)+l
				if (i.lt.(n-1)) goto 40
				if (k.gt.2*kbem) ik=k-2*kbem	! for the last elm
				if (l.gt.2*kbem) il=l-2*kbem	! for the last elm
40				xmtf(ik,il)=xmtf(ik,il)+xmm(k,l)
			enddo
		enddo
100	continue
!
	return
	end


!****************************************************************************************
!
!							           OTHERMAT1 subroutine
!
!	This subroutine calculates the gihd1 & xgihd1 matrices.
!	This subroutine is called in these subroutines: GHMATD1
!
!	variables used are:
!
!		k1			: counter which shows the number of each D1 (finite) zone (1:nbed1)
!		n			: total number of nodes of each D1 (finite) sub-zone
!		kbem		: number of degree of freedom per node in BE zone;
!								if IBEM=[1] -> KBEM=[2] [X,Y]
!								if IBEM=[2] -> KBEM=[3] [X,Y,Pw]
!								if IBEM=[3] -> KBEM=[4] [X,Y,Pw,Pa]
!		hd1			: total H matrix in D1 (finite) area concerning all of the h matrices
!					  of each element (ehd1)
!		gid1		: inverse matrice of gd1
!		xmtf		: M (distributed matric) of all of elements
!		gihd1		: gid1*hd1
!		xgid1		: M*gid1
!		xgihd1		: M*gid1*hd1
!
!
!	INPUT	:
!
!	OUTPUT	:
!
!****************************************************************************************
!
	subroutine othermat1(k1,n,hd1,xmtf,gid1,gihd1,xgid1,xgihd1)
!
	implicit double precision (a-h,o-z)
!
	common /a/ nnp,ifem,ibem,idyn,ieaq,nr
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /a3/ ne8d1,ne8c1,ne8u1,ne3u1,nload1,nbcx1,nbcy1,nbcw1,nbca1,m8d1,m8c1,m8u1,&
				nbed11,nbed21,nbed31,nbeint1,ne8dd,ne8cc,ne8uu
!
	dimension hd1(n12gh,n12gh),xmtf(n12gh,n12gh),gid1(n12gh,n12gh,nbed11),&
			  gihd1(n12gh,n12gh,nbed11),xgid1(n12gh,n12gh,nbed11),xgihd1(n12gh,n12gh,nbed11)
!
!
!	------------------- initializing -------------------
!
	n1=kbem*n
	do iy=1,n1
		do ix=1,n1
			gihd1(iy,ix,k1)=0.d0
			xgid1(iy,ix,k1)=0.d0
			xgihd1(iy,ix,k1)=0.d0
		enddo
	enddo
!
!
!	------------------- gihd1 matrice -------------------
!
	do iy=1,n1
		do ix=1,n1
			do i=1,n1
				gihd1(iy,ix,k1)=gihd1(iy,ix,k1)+gid1(iy,i,k1)*hd1(i,ix)
			enddo
		enddo
	enddo
!
!	------------------- xgid1 matrice -------------------
!
	do iy=1,n1
		do ix=1,n1
			do i=1,n1
				xgid1(iy,ix,k1)=xgid1(iy,ix,k1)+xmtf(iy,i)*gid1(i,ix,k1)
			enddo
		enddo
	enddo
!
!	------------------- xgihd1 matrice -------------------
!
	do iy=1,n1
		do ix=1,n1
			do i=1,n1
				xgihd1(iy,ix,k1)=xgihd1(iy,ix,k1)+xmtf(iy,i)*gihd1(i,ix,k1)
			enddo
		enddo
	enddo
!
	return
	end


!****************************************************************************************
!
!							           OUTPUT subroutine
!
!	This subroutine prints out results of each load step.
!	This subroutine is called in: CONTROL
!	This subroutine calls:		  SHAP8N,BMAT8N
!
!	variables used are:
!
!		ie8d	:
!		ie8c	:
!		ie8u	:
!		x		: coordinates of nodes (nodal data) in FE & BE areas;
!							x(1,1:nnp) : X-coord
!							x(2,1:nnp) : Y-coord
!		id		: degrees of freedom of each node from "1" to "mdof"; for the nodes with
!				  the boundary conditions, each node takes from "mdof+1" to "mdof+nbc"
!				  [X,Y]; which are free for the dry elm &
!				  [X,Y,Pw]; which are free for the saturated elm
!				  [X,Y,Pw,Pa]; which are free for the unsaturated elm
!		sig8d	:
!		sig8c	:
!		sig8u	:
!		dt		: R vector in this time step "N"
!		iact	:
!		iterp	: number of iloads for which the results are printed out in the output file
!		itprt	: iload's numbers which are selected to print out the results;
!							[0] if no printout is expected
!		itime	: time increment from 1 to "ntime"
!							[1] : first time increment,
!							[>1]: higher time increment
!		dtime	: time step; in seismic loading for all "1:ntime", it is constant
!		dt		: R vector in this time step "N"
!		epi		:
!		xe		:
!		ye		:
!		unode	:
!
!	INPUT	:
!
!
!	OUTPUT	:
!
!
!****************************************************************************************
!
	subroutine output(ie8d,ie8c,ie8u,id,x,sig8d,sig8c,sig8u,dt,iact,sig3u)
!
	implicit double precision (a-h,o-z)
!
	common /a/ nnp,ifem,ibem,idyn,ieaq,nr
	common /a1/ ne8d,ne8c,ne8u
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /a3/ ne8d1,ne8c1,ne8u1,ne3u1,nload1,nbcx1,nbcy1,nbcw1,nbca1,m8d1,m8c1,m8u1,&
				nbed11,nbed21,nbed31,nbeint1,ne8dd,ne8cc,ne8uu
	common /b/ m8d,m8c,m8u,m8dep(10),m8cep(10),m8uep(10)
	common /c/ nload,ntime,itmax,itime,iter,iload,time,dtime,dtimei,dtimem,xinc
	common /h/ mdof,mdofn,ls,ls1,isolv
	common /chrg1/ nstep,istep
	common /i/iterp,itprt(100),ne8do,ne8co,ne8uo,ie8dout(10),ie8cout(10),ie8uout(10),nnpo,&
			  inpout(100)
	common /unbem/ nelmu,ie3u(3,500,5)
!
	dimension ie8d(9,ne8d1),ie8c(9,ne8c1),ie8u(9,ne8u1),id(4,nnp),x(2,nnp),epi(3),&
			  sig8d(9,5,ne8d1),sig8c(9,6,ne8c1),sig8u(9,7,ne8u1),dt(mdofn),iact(nnp),&
			  xe(8),ye(8),xelm(3),yelm(3),unode(200),f(8),pfs(8),pft(8),b(3,16),u(16),&
			  sig3u(5,500,nbed11)

!
!
	call outgid(itime,nnp,mdofn,iact,id,dt)
!
!	------- deciding whether to print all element stresses & nodal variables -------
!
	l=0
	do 10 i=1,iterp
		j=itprt(i)
		if (j.eq.itime) l=1
10	continue
	if (l.ne.1) goto 40
!
!
!	----- print nodal variables and set displacements of builtup nodes to zero -----
!
	write (7,1) itime,time,dtime
	write (7,2)
	do 20 i=1,nnp
		j=iact(i)
		if (j.eq.0) goto 20
		i1=id(1,i)
		i2=id(2,i)
		i3=id(3,i)
		i4=id(4,i)
		u1=0.0d0
		u2=0.0d0
		pw=0.0d0
		pa=0.0d0
		if (i1.le.mdofn) u1=dt(i1)
        if (i2.le.mdofn) u2=dt(i2)
        if (i3.le.mdofn) pw=dt(i3)
		if (i4.le.mdofn) pa=dt(i4)
		write (7,3) i,x(1,i),x(2,i),u1,u2,pw,pa
20	continue
!
!
	if (ne8d.eq.0) goto 30
	write (7,4)
	do i=1,ne8d
		do j=1,8
			k=ie8d(j,i)
			xe(j)=x(1,k)
			ye(j)=x(2,k)
		enddo
		xm=0.125*(xe(1)+xe(2)+xe(3)+xe(4)+xe(5)+xe(6)+xe(7)+xe(8))
		ym=0.125*(ye(1)+ye(2)+ye(3)+ye(4)+ye(5)+ye(6)+ye(7)+ye(8))
		sx=-sig8d(5,1,i)
		sy=-sig8d(5,2,i)
        sxy=sig8d(5,3,i)
!		sm=sig8d(5,4,i)
		write (7,5) i,xm,ym,sx,sy,sxy !,sm
	enddo
30	continue
!
!
	if (ne8c.eq.0) goto 35
	write (7,6)
	do i=1,ne8c
		do j=1,8
			k=ie8c(j,i)
			xe(j)=x(1,k)
			ye(j)=x(2,k)
		enddo
		xm=0.125*(xe(1)+xe(2)+xe(3)+xe(4)+xe(5)+xe(6)+xe(7)+xe(8))
        ym=0.125*(ye(1)+ye(2)+ye(3)+ye(4)+ye(5)+ye(6)+ye(7)+ye(8))
        sx=-sig8c(5,1,i)
		sy=-sig8c(5,2,i)
        sxy=sig8c(5,3,i)
		pp=sig8c(5,6,i)
		write (7,5) i,xm,ym,sx,sy,sxy,pp
	enddo
!
!
35	continue
!
!
	if (ne8u.eq.0) goto 36
	write (7,7)
	do i=1,ne8u
		do j=1,8
			k=ie8u(j,i)
			xe(j)=x(1,k)
			ye(j)=x(2,k)
		enddo
		xm=0.125*(xe(1)+xe(2)+xe(3)+xe(4)+xe(5)+xe(6)+xe(7)+xe(8))
        ym=0.125*(ye(1)+ye(2)+ye(3)+ye(4)+ye(5)+ye(6)+ye(7)+ye(8))
        sx=-sig8u(5,1,i)
		sy=-sig8u(5,2,i)
        sxy=sig8u(5,3,i)
		ppw=sig8u(5,6,i)
		ppa=sig8u(5,7,i)
		write (7,5) i,xm,ym,sx,sy,sxy,ppw,ppa
	enddo
!
36	continue
!
!
	if (ibem.ne.3) goto 40
	if (nbed1.eq.0) goto 40
	write (7,8)
	do k1=1,nbed1
		do i=1,nelmu
			do j=1,3
				k=ie3u(j,i,k1)
				xelm(j)=x(1,k)
				yelm(j)=x(2,k)
			enddo
			xm=(xelm(1)+xelm(2)+xelm(3))/3.
			ym=(yelm(1)+yelm(2)+yelm(3))/3.
			sx=-sig3u(1,i,k1)
			sy=-sig3u(2,i,k1)
			sxy=sig3u(3,i,k1)
			ppw=-sig3u(4,i,k1)+sig3u(5,i,k1)
			ppa=sig3u(5,i,k1)
			write (7,5) i,xm,ym,sx,sy,sxy,ppw,ppa
		enddo
	enddo
!
!	----- printing results of selected elements -----
!
40	continue
!
	ii=9
!	nodes
	if (nnpo.eq.0) goto 50
	do i=1,nnpo
		ii=ii+i
        j=inpout(i)
        i1=id(1,j)
		i2=id(2,j)
		i3=id(3,j)
		i4=id(4,j)
        u1=0.0d0
		u2=0.0d0
		pw=0.0d0
		pa=0.0d0
		if (i1.le.mdofn) u1=dt(i1)
		if (i2.le.mdofn) u2=dt(i2)
		if (i3.le.mdofn) pw=dt(i3)
		if (i4.le.mdofn) pa=dt(i4)
		write (ii,15) time,u1,u2,pw,pa
!
		unode(2*i-1)=0.0d0
		unode(2*i)=0.0d0
		if (i1.le.mdofn) unode(2*i-1)=dt(i1)
		if (i2.le.mdofn) unode(2*i)=dt(i2)
	enddo
	ii=ii+1
	write (ii,17) time,(unode(i),i=1,200)
50	continue
!
!
!	----- dry elements at center -----
!
	if (ne8do.eq.0) goto 60
	do i=1,ne8do
		ii=ii+i
		j=ie8dout(i)
		sx=-sig8d(5,1,j)
		sy=-sig8d(5,2,j)
		sxy=sig8d(5,3,j)
		sm=(sx+sy)/2.d0
		tau=dsqrt ((sx-sy)**2+sxy**2)
		p=(2*sx+sy)/3
		q=(sy-sx)
!
!	deformation
		do jj=1,8
			k=ie8d(jj,j)
			xe(jj)=x(1,k)
			ye(jj)=x(2,k)
			k1=id(1,k)
			k2=id(2,k)
			u(2*jj-1)=dt(k1)
			u(2*jj)=dt(k2)
		enddo
		call shap8n(0.d0,0.d0,f,pfs,pft)
		call bmat8n(xe,ye,pfs,pft,b,detj)
		do jj=1,3
			epi(jj)=0.0d0
			do k=1,16
				epi(jj)=epi(jj)+b(jj,k)*u(k)
			enddo
		enddo
		gamma=dsqrt ((epi(1)-epi(2))**2+epi(3)**2)
		write (ii,15) time,sx,sy,sxy,sm,tau,-epi(1),-epi(2),epi(3),gamma
	enddo
60	continue
!
!
!	----- consolidation elements at center -----
!
	if (ne8co.eq.0) goto 65
	do i=1,ne8co
		ii=ii+i
		j=ie8cout(i)
		sx=-sig8c(5,1,j)
		sy=-sig8c(5,2,j)
		pp=sig8c(5,6,j)
		sxy=sig8c(5,3,j)
		sm=(sx+sy)/2.d0
		tau=dsqrt ((sx-sy)**2+sxy**2)
		p=(2*sx+sy)/3
		q=(sy-sx)
!
!	deformation
		do jj=1,8
			k=ie8c(jj,j)
			xe(jj)=x(1,k)
			ye(jj)=x(2,k)
			k1=id(1,k)
			k2=id(2,k)
			u(2*jj-1)=dt(k1)
			u(2*jj)=dt(k2)
		enddo
		call shap8n(0.d0,0.d0,f,pfs,pft)
		call bmat8n(xe,ye,pfs,pft,b,detj)
		do jj=1,3
			epi(jj)=0.0d0
			do k=1,16
				epi(jj)=epi(jj)+b(jj,k)*u(k)
			enddo
		enddo
		gamma=dsqrt ((epi(1)-epi(2))**2+epi(3)**2)
		write (ii,16) time,sx,sy,sxy,sm,tau,pp,-epi(1),-epi(2),epi(3),gamma
	enddo
!
!
!	----- unsaturated elements at center -----
!
65	if (ne8uo.eq.0) goto 70
	do i=1,ne8uo
		ii=ii+i
		j=ie8uout(i)
		sx=-sig8u(5,1,j)
		sy=-sig8u(5,2,j)
		suc=sig8u(5,6,j)
		ppa=sig8u(5,7,j)
		sxy=sig8u(5,3,j)
		sm=(sx+sy)/2.d0
		tau=dsqrt((sx-sy)**2+sxy**2)
		p=(2*sx+sy)/3
		q=(sy-sx)
!
!	deformation
		do jj=1,8
			k=ie8u(jj,j)
			xe(jj)=x(1,k)
			ye(jj)=x(2,k)
			k1=id(1,k)
			k2=id(2,k)
			u(2*jj-1)=dt(k1)
			u(2*jj)=dt(k2)
		enddo
		call shap8n(0.d0,0.d0,f,pfs,pft)
		call bmat8n(xe,ye,pfs,pft,b,detj) !for unsat is not yet ready
		do jj=1,3
			epi(jj)=0.0d0
			do k=1,16
				epi(jj)=epi(jj)+b(jj,k)*u(k)
			enddo
		enddo
		gamma=dsqrt ((epi(1)-epi(2))**2+epi(3)**2)
		write (ii,16) time,sx,sy,sxy,sm,tau,suc,ppa,-epi(1),-epi(2),epi(3),gamma
	enddo
70	continue
!
!
!	------------- formats -------------
!
1	FORMAT (//10x,'solution data'//&
			  5x,'time step number              = ',i5/&
			  5x,'time                          = ',e10.3/&
			  5x,'time step size                = ',e10.3/)
2	FORMAT (//10x,'nodal displacements and excess pore pressures'//&
			      'node no.',4x,'x-coord',4x,'y-coord',6x,'u-disp',6x,&
				  'v-disp',6x,'p-water',6x,'p-air'//)
3	FORMAT (2x,i5,6(2x,e10.3))
4	FORMAT (//10x,'drained element stresses'//&
			      'e.no.',5x,'x',7x,'y',3x,'sig-x',5x,'sig-y',6x,'sig-xy',6x,'sig-xym'//)
5	FORMAT (i4,2(1x,f6.2),4(1x,e15.6))
6	FORMAT (//10x,'consolidation element effective stresses'//&
				  'e.no.',5x,'x',7x,'y',3x,'sig-x',5x,'sig-y',4x,'sig-xy',3x,&
				  'sig-xym',4x,'pore-pp'//)
7	FORMAT (//10x,'unsaturated element effective stresses'//&
				  'e.no.',5x,'x',7x,'y',3x,'sig-x',5x,'sig-y',4x,'sig-xy',3x,&
				  'sig-xym',4x,'Wpore-ppw',4x,'Apore-ppa'//)
8	FORMAT (//10x,'unsaturated element in BE zone'//&
				  'e.no.',5x,'x',7x,'y',3x,'sig-x',5x,'sig-y',4x,'sig-xy',3x,&
				  'sig-xym',4x,'Wpore-ppw',4x,'Apore-ppa'//)
9	FORMAT (i4,2(1x,f6.2),5(1x,e10.3))
15	FORMAT (10(1x,e15.6))
16	FORMAT (12(1x,e10.4))
17	FORMAT (e12.4,200e12.4)
!
	return
	end

  SUBROUTINE DFTF(N, DAT,FSUM)
        IMPLICIT REAL*8 (A-H,O-Z)
        COMPLEX*16 DAT(*),FSUM(*),TEMP
        DATA PI/3.141592653589793D0/
        IM=1
        IB=-N/2
        IE=IB+N-1
        DO M=1,N
          FSUM(IM)=0
          IDAT=1
          DO I=IB,IE
            ARG=(2*PI*I*M)/N
            TEMP=EXP((0,1)*ARG)   
            FSUM(IM)=FSUM(IM)+TEMP*DAT(IDAT)
            IDAT=IDAT+1
          ENDDO
          IM=IM+1
        ENDDO
        RETURN
        END
	

		
	subroutine owcqm(dg,df)
	
	implicit double precision (a-h,o-z)
	common /a/ nnp,ifem,ibem,idyn,ieaq,nr
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /c/ nload,ntime,itmax,itime,iter,iload,time,dtime,dtimei,dtimem,xinc
	common /point/ dist(2),ra,rd(2),rdn,eta(2)
	common /coefL/aa(2,2),dpaa(2,2,2),bb(2,2),dpbb(2,2,2),cc(2,2)
	common /mate/zlamda,zmuy,znuu,znu,zalpha,zrho,zrhof,zm,zk,zr,zn
	common /umate/ d(3,3),dfs(3),xn,xkw,cww,cwg,xka,cgg,cgw,ev,sat,g1,rs,rw,raa,rmix,xm1,&
				   bt,et,xlambda,xmu,henry
	common /ocqm/ EPScqm,taw,Icqm,Lcqm,NPcqm
	common /comp/ icomp
	COMPLEX(8) W1,W,Z,GammaZ,SLAP,cqmG(kbem,kbem),cqmF(kbem,kbem),RTD(4),RTQ(2),RTSD(4),&
			   RTSDI(2),RTSQ,&
			   zz,CBS(2),CBSS(4),K0(4),K1(4),KQ0(2),KQ1(2),KQ2(2),KSD0(3),KSD1(3),DKSD0(3),DKSD1(3),&
			   KSDI0(2),KSDI1(2),DKSDI0(2),DKSDI1(2),&
			   KS0,KS1,KS2,KS3,DKQ0(2),DKQ1(2),DK0(4),DK1(4),&
			   GAM11,GAM12,GAM13,GAM21,GAM22,GAM31,GAM32,GAM33,GAM1,GAM2,GAM3,&
			   cqmG11(ntime),cqmG12(ntime),cqmG13(ntime),cqmG14(ntime),&
			   cqmG21(ntime),cqmG22(ntime),cqmG23(ntime),cqmG24(ntime),&
			   cqmG31(ntime),cqmG32(ntime),cqmG33(ntime),cqmG34(ntime),&
			   cqmG41(ntime),cqmG42(ntime),cqmG43(ntime),cqmG44(ntime),&
			   cqmF11(ntime),cqmF12(ntime),cqmF13(ntime),cqmF14(ntime),&
			   cqmF21(ntime),cqmF22(ntime),cqmF23(ntime),cqmF24(ntime),&
			   cqmF31(ntime),cqmF32(ntime),cqmF33(ntime),cqmF34(ntime),&
			   cqmF41(ntime),cqmF42(ntime),cqmF43(ntime),cqmF44(ntime),&
			   FcqmG11(ntime),FcqmG12(ntime),FcqmG13(ntime),FcqmG14(ntime),&
			   FcqmG21(ntime),FcqmG22(ntime),FcqmG23(ntime),FcqmG24(ntime),&
			   FcqmG31(ntime),FcqmG32(ntime),FcqmG33(ntime),FcqmG34(ntime),&
			   FcqmG41(ntime),FcqmG42(ntime),FcqmG43(ntime),FcqmG44(ntime),&
			   FcqmF11(ntime),FcqmF12(ntime),FcqmF13(ntime),FcqmF14(ntime),&
			   FcqmF21(ntime),FcqmF22(ntime),FcqmF23(ntime),FcqmF24(ntime),&
			   FcqmF31(ntime),FcqmF32(ntime),FcqmF33(ntime),FcqmF34(ntime),&
			   FcqmF41(ntime),FcqmF42(ntime),FcqmF43(ntime),FcqmF44(ntime),&
			   FcqmG(kbem,kbem),FcqmF(kbem,kbem),FNOMEGAG(kbem,kbem),FNOMEGAF(kbem,kbem)

	dimension  dg(kbem,kbem),df(kbem,kbem),vel(2),uvel(2) !,cqmID(100000,3) 
	
!
!
	pi=DACOS (-1.d0)
	if (ibem.eq.1) then
		vel(1)=dsqrt ((zlamda+2*zmuy)/zrho)
		vel(2)=dsqrt (zmuy/zrho)
	endif
	if (idyn.eq.2.and.ibem.eq.2.and.kfsol.eq.3) then
		vel(1)=dsqrt ( (zlamda+2*zmuy)/(zrho-zalpha*zrhof) ) !dsqrt( (zlamda+2*zmuy+zalpha**2.d0*zr/zn)/zrho)
		vel(2)=dsqrt (zmuy/zrho)
	endif
	if (idyn.eq.2.and.ibem.eq.3) then
!		uvel(1)=dsqrt ()
!		uvel(2)=dsqrt ()
	endif
!
!
	do i=1,kbem
		do j=1,kbem
			FNOMEGAG(i,j)=0.d0
			FNOMEGAF(i,j)=0.d0
			FcqmG(i,j)=0.d0
			FcqmF(i,j)=0.d0
			cqmG(i,j)=0.d0
			cqmF(i,j)=0.d0
			dg(i,j)=0.d0
			df(i,j)=0.d0
		enddo
	enddo
!
	do i=1,ntime
		cqmG11(i)=0.d0
		cqmG12(i)=0.d0
		cqmG21(i)=0.d0
		cqmG22(i)=0.d0
		cqmG13(i)=0.d0
		cqmG23(i)=0.d0
		cqmG33(i)=0.d0
		cqmG14(i)=0.d0
		cqmG24(i)=0.d0
		cqmG34(i)=0.d0
		cqmG43(i)=0.d0
		cqmG44(i)=0.d0
		cqmF11(i)=0.d0
		cqmF12(i)=0.d0
		cqmF21(i)=0.d0
		cqmF22(i)=0.d0
		cqmF13(i)=0.d0
		cqmF23(i)=0.d0
		cqmF33(i)=0.d0
		cqmF14(i)=0.d0
		cqmF24(i)=0.d0
		cqmF34(i)=0.d0
		cqmF43(i)=0.d0
		cqmF44(i)=0.d0
	enddo
!
	if (ibem.eq.1) then
		x=vel(1)*dtime
		if (ra.gt.x*itime) goto 100
	endif
!	if (idyn.eq.2.and.ibem.eq.2.and.kfsol.eq.3) then
!		x=vel(1)*dtime
!		if (ra.gt.x*itime) goto 100
!	endif
!	if (idyn.eq.2.and.ibem.eq.3) then
!		x=uvel(1)*dtime
!		if (ra.gt.x*itime) goto 100
!	endif
!
	Lcqm=ntime
	taw=EPScqm**(1.d0/(2.d0*ntime))
	W1=CMPLX(0.d0,2.d0*pi/Lcqm)
	W=CDexp(W1)
!
10	do l=1,Lcqm
		Z=taw*W**(l-1)
		GammaZ=0.d0
		do i=1,NPcqm
			GammaZ=GammaZ+((1.d0-Z)**i)/i
		enddo
		SLAP=GammaZ/dtime
!
		FundamentalSolution: SELECT CASE (ibem)
			CASE (1)
				call weldynfsol(W,l,SLAP,cqmG,cqmF)
!
!		SATURATED SOIL
!
			CASE (2)
				call ROOTS(SLAP,RTD,RTQ,RTSD,RTSDI,RTSQ)
				if (kfsol.eq.0.and.idyn.eq.1) then        !QUASI-STATIC
					zz=Cdsqrt(RTSQ)*ra
					call cbesk01 (0.d0,zz,4,CBSS)
					KS0=CBSS(1)
					KS1=CBSS(2)
					KS2=CBSS(3)
					KS3=CBSS(4)
					call wsatQSTfsolG(SLAP,RTSQ,KS0,KS1,KS2,cqmG)
					call wsatQSTfsolF(SLAP,RTSQ,KS0,KS1,KS2,KS3,cqmF)
!
				else if (kfsol.eq.3.and.idyn.eq.2.and.icomp.eq.1) then   !DYNAMIC COMPRESSIBLE
					do j=1,3
						zz=Cdsqrt(RTSD(j))*ra
						CBS(1)=(0.0d0,0.0d0)
						CBS(2)=(0.0d0,0.0d0)
						if (Cdabs(zz).gt.150.d0) then
							KSD0(j)=0.d0
							KSD1(j)=0.d0
							goto 12
						endif
						call cbesk01 (0.d0,zz,2,CBS)
						KSD0(j)=CBS(1)
						KSD1(j)=CBS(2)
12						continue
						DKSD0(j)=-Cdsqrt(RTSD(j))*KSD1(j)
						DKSD1(j)=-Cdsqrt(RTSD(j))*(KSD0(j)+1.d0*KSD1(j)/zz)
					enddo
					call wsatDYNfsolG(SLAP,RTSD,RTSDI,KSD0,KSD1,KSDI0,KSDI1,cqmG)
					call wsatDYNfsolF(SLAP,RTSD,RTSDI,KSD0,KSD1,KSDI0,KSDI1,DKSD0,DKSD1,&
									  DKSDI0,DKSDI1,cqmG,cqmF)
!
				else if (kfsol.eq.3.and.idyn.eq.2.and.icomp.eq.0) then   !DYNAMIC INCOMPRESSIBLE
					do j=1,2
						zz=Cdsqrt(RTSDI(j))*ra*SLAP !SLAP
						CBS(1)=(0.0d0,0.0d0)
						CBS(2)=(0.0d0,0.0d0)
						if (Cdabs(zz).gt.150.d0) then
							KSDI0(j)=0.d0
							KSDI1(j)=0.d0
							goto 13
						endif
						call cbesk01 (0.d0,zz,2,CBS)
						KSDI0(j)=CBS(1)
						KSDI1(j)=CBS(2)
13						continue
						DKSDI0(j)=-Cdsqrt(RTSDI(j))*SLAP*KSDI1(j)
						DKSDI1(j)=-Cdsqrt(RTSDI(j))*SLAP*(KSDI0(j)+1.d0*KSDI1(j)/zz)
					enddo
					call wsatDYNfsolG(SLAP,RTSD,RTSDI,KSD0,KSD1,KSDI0,KSDI1,cqmG)
					call wsatDYNfsolF(SLAP,RTSD,RTSDI,KSD0,KSD1,KSDI0,KSDI1,DKSD0,DKSD1,&
									  DKSDI0,DKSDI1,cqmG,cqmF)
!
				endif
!	calcule analytique saturated soil
!
					call SATANAL
!
!		UNSATURATED SOIL
!
			CASE (3)
				call ROOTS(SLAP,RTD,RTQ,RTSD,RTSDI,RTSQ)
				if (idyn.eq.2) then
					do j=1,4
						zz=Cdsqrt(RTD(j))*ra
						CBS(1)=(0.0d0,0.0d0)
						CBS(2)=(0.0d0,0.0d0)
						if (Cdabs(zz).gt.150.d0) then
							K0(j)=0.d0
							K1(j)=0.d0
							goto 15
						endif
						call cbesk01 (0.d0,zz,2,CBS)
						K0(j)=CBS(1)
						K1(j)=CBS(2)
15						continue
						DK0(j)=-Cdsqrt(RTD(j))*K1(j)
						DK1(j)=-Cdsqrt(RTD(j))*(K0(j)+1.d0*K1(j)/zz)
					enddo
					call wunsatdynfsolG(SLAP,RTD,K0,K1,cqmG)
					call wunsatdynfsolF(SLAP,RTD,K0,K1,DK0,DK1,cqmG,cqmF)
!
				else if (idyn.eq.1) then
					do j=1,2
						zz=Cdsqrt(RTQ(j))*ra
						CBS(1)=(0.0d0,0.0d0)
						CBS(2)=(0.0d0,0.0d0)
						if (Cdabs(zz).gt.150.d0.OR.Cdabs(zz).eq.0.d0) then
							KQ0(j)=0.d0
							KQ1(j)=0.d0
							goto 20
						endif
						CBSS = Bessel_jn(3,real(zz))  !call DCBKS (0.d0,zz,3,CBSS)
						KQ0(j)=CBSS(1)
						KQ1(j)=CBSS(2)
						KQ2(j)=CBSS(3)
20						continue
						DKQ0(j)=-Cdsqrt(RTQ(j))*KQ1(j)
						DKQ1(j)=-Cdsqrt(RTQ(j))*(KQ0(j)+1.d0*KQ1(j)/zz)
					enddo
					call wunsatqstfsolG1(SLAP,RTQ,KQ0,KQ1,KQ2,cqmG)
					call wunsatqstfsolF1(SLAP,RTQ,KQ0,KQ1,KQ2,DKQ0,DKQ1,cqmG,cqmF)
				endif
!
		end SELECT FundamentalSolution
!
!
		cqmG11(l)=cqmG(1,1)
		cqmG22(l)=cqmG(2,2)
		cqmG12(l)=cqmG(1,2)
		cqmG21(l)=cqmG(2,1)
		cqmF11(l)=cqmF(1,1)
		cqmF22(l)=cqmF(2,2)
		cqmF12(l)=cqmF(1,2)
		cqmF21(l)=cqmF(2,1)
		if (kbem.ge.3) then
			cqmG13(l)=cqmG(1,3)
			cqmG23(l)=cqmG(2,3)
			cqmG31(l)=cqmG(3,1)
			cqmG32(l)=cqmG(3,2)
			cqmG33(l)=cqmG(3,3)
			cqmF13(l)=cqmF(1,3)
			cqmF23(l)=cqmF(2,3)
			cqmF31(l)=cqmF(3,1)
			cqmF32(l)=cqmF(3,2)
			cqmF33(l)=cqmF(3,3)
		endif
		if (kbem.eq.4) then
			cqmG14(l)=cqmG(1,4)
			cqmG24(l)=cqmG(2,4)
			cqmG34(l)=cqmG(3,4)
			cqmG41(l)=cqmG(4,1)
			cqmG42(l)=cqmG(4,2)
			cqmG43(l)=cqmG(4,3)
			cqmG44(l)=cqmG(4,4)
			cqmF14(l)=cqmF(1,4)
			cqmF24(l)=cqmF(2,4)
			cqmF34(l)=cqmF(3,4)
			cqmF41(l)=cqmF(4,1)
			cqmF42(l)=cqmF(4,2)
			cqmF43(l)=cqmF(4,3)
			cqmF44(l)=cqmF(4,4)
		endif
	enddo
!
!
	call  DFTF (ntime,cqmG11,FcqmG11)
	call  DFTF (ntime,cqmG22,FcqmG22)
	call  DFTF (ntime,cqmG12,FcqmG12)
	call  DFTF (ntime,cqmG21,FcqmG21)
	call  DFTF (ntime,cqmF11,FcqmF11)
	call  DFTF (ntime,cqmF22,FcqmF22)
	call  DFTF (ntime,cqmF12,FcqmF12)
	call  DFTF (ntime,cqmF21,FcqmF21)
	if (kbem.ge.3) then
		call  DFTF (ntime,cqmG13,FcqmG13)
		call  DFTF (ntime,cqmG23,FcqmG23)
		call  DFTF (ntime,cqmG31,FcqmG31)
		call  DFTF (ntime,cqmG32,FcqmG32)
		call  DFTF (ntime,cqmG33,FcqmG33)
		call  DFTF (ntime,cqmF13,FcqmF13)
		call  DFTF (ntime,cqmF23,FcqmF23)
		call  DFTF (ntime,cqmF31,FcqmF31)
		call  DFTF (ntime,cqmF32,FcqmF32)
		call  DFTF (ntime,cqmF33,FcqmF33)
	endif
	if (kbem.eq.4) then
		call  DFTF (ntime,cqmG14,FcqmG14)
		call  DFTF (ntime,cqmG24,FcqmG24)
		call  DFTF (ntime,cqmG34,FcqmG34)
		call  DFTF (ntime,cqmG41,FcqmG41)
		call  DFTF (ntime,cqmG42,FcqmG42)
		call  DFTF (ntime,cqmG43,FcqmG43)
		call  DFTF (ntime,cqmG44,FcqmG44)
		call  DFTF (ntime,cqmF14,FcqmF14)
		call  DFTF (ntime,cqmF24,FcqmF24)
		call  DFTF (ntime,cqmF34,FcqmF34)
		call  DFTF (ntime,cqmF41,FcqmF41)
		call  DFTF (ntime,cqmF42,FcqmF42)
		call  DFTF (ntime,cqmF43,FcqmF43)
		call  DFTF (ntime,cqmF44,FcqmF44)
	endif
!
!	do i=1,ntime
		FcqmG(1,1)=FcqmG11(itime)
	 	FcqmG(2,2)=FcqmG22(itime)
		FcqmG(1,2)=FcqmG12(itime)
		FcqmG(2,1)=FcqmG21(itime)
		FcqmF(1,1)=FcqmF11(itime)
		FcqmF(2,2)=FcqmF22(itime)
		FcqmF(1,2)=FcqmF12(itime)
		FcqmF(2,1)=FcqmF21(itime)
		if (kbem.ge.3) then
			FcqmG(1,3)=FcqmG13(itime)
			FcqmG(2,3)=FcqmG23(itime)
			FcqmG(3,1)=FcqmG31(itime)
			FcqmG(3,2)=FcqmG32(itime)
			FcqmG(3,3)=FcqmG33(itime)
			FcqmF(1,3)=FcqmF13(itime)
			FcqmF(2,3)=FcqmF23(itime)
			FcqmF(3,1)=FcqmF31(itime)
			FcqmF(3,2)=FcqmF32(itime)
			FcqmF(3,3)=FcqmF33(itime)
		endif
		if (kbem.eq.4) then
			FcqmG(1,4)=FcqmG14(itime)
			FcqmG(2,4)=FcqmG24(itime)
			FcqmG(3,4)=FcqmG34(itime)
			FcqmG(4,1)=FcqmG41(itime)
			FcqmG(4,2)=FcqmG42(itime)
			FcqmG(4,3)=FcqmG43(itime)
			FcqmG(4,4)=FcqmG44(itime)
			FcqmF(1,4)=FcqmF14(itime)
			FcqmF(2,4)=FcqmF24(itime)
			FcqmF(3,4)=FcqmF34(itime)
			FcqmF(4,1)=FcqmF41(itime)
			FcqmF(4,2)=FcqmF42(itime)
			FcqmF(4,3)=FcqmF43(itime)
			FcqmF(4,4)=FcqmF44(itime)
		endif
!	enddo
!
!
	do i=1,kbem
		do j=1,kbem
			FNOMEGAG(i,j)=taw**(-itime+1)*FcqmG(i,j)/ntime
			FNOMEGAF(i,j)=taw**(-itime+1)*FcqmF(i,j)/ntime
			dg(i,j)=REAL(FNOMEGAG(i,j))
			df(i,j)=REAL(FNOMEGAF(i,j))
		enddo
	enddo
!
!
100	continue
!
	return
	end


!*******************************************************************************************
!
	subroutine ROOTDYN(SLAP,RTD)
!
	implicit double precision (a-h,o-z)
	COMPLEX(8) RTD(4),A(3),D1,D3,c3p,c2m,c3m,p,q,rr,qq,cc,ww,SLAP
	common /e/ gammaw,gammas,gammaa,grav,atmp
	common /umate/ d(3,3),dfs(3),xn,xkw,cww,cwg,xka,cgg,cgw,ev,sat,g1,rs,rw,raa,rmix,xm1,&
				   bt,et,xlambda,xmu,henry
!
!
	do i=1,4
		RTD(i)=0.d0
	enddo
!
	RTD(1)=(rmix/xmu)*SLAP**2.d0
!
!	calculate rtd2^2+rtd3^2+rtd4^2
!
	c1=xlambda+2.d0*xmu
	fs=dfs(1)
	D1=DCMPLX (0.d0,0.d0)
	D1=(rmix/c1+rw*fs/c1+raa*(1.d0-fs)/c1)*SLAP**2.d0
	if (xka.ne.0.d0) D1=D1-(cgg/xka+(1.d0-sat)*(1.d0-fs)/(c1*xka))*SLAP
	if (xkw.ne.0.d0) D1=D1-(cww/xkw+sat*fs/(c1*xkw))*SLAP
	c3p=D1
	A(3)=-c3p
!
!	calculate rtd2^2*rtd3^2+rtd2^2*rtd4^2+rtd3^2*rtd4^2
!
	c2=fs*cgg-(1.d0-fs)*cwg
	c3=-fs*cwg+cww*(1.d0-fs)
	D3=DCMPLX (0.d0,0.d0)
	if (xka.ne.0.d0) D3=-(rmix*cgg/(c1*xka)+rw*c2/(c1*xka))*SLAP**3.d0
	if (xkw.ne.0.d0) D3=D3-(rmix*cww/(c1*xkw)+raa*c3/(c1*xkw))*SLAP**3.d0
	if (xka.ne.0.d0.and.xkw.ne.0.d0) D3=D3+((cww*cgg-cwg**2.d0)/(xkw*xka)+&
							sat*c2/(c1*xka*xkw)+(1.d0-sat)*c3/(c1*xkw*xka))*SLAP**2.d0
	c2m=D3
	A(2)=c2m
!
!	calculate rtd2^2*rtd3^2*rtd4^2
!
	if (xka.ne.0.d0.and.xkw.ne.0.d0) c3m=rmix*SLAP**4.d0*(cww*cgg-cwg**2.d0)/(c1*xkw*xka)
	A(1)=-c3m
!
!	solve the cubic polynomial with complex coefficients
!
	pi=DACOS (-1.d0)
	onethird=1.d0/3.d0
	twothird=2.d0/3.d0
	p=A(2)-onethird*(A(3)**2.d0)
	q=2.d0*(onethird*A(3))**3.d0-onethird*A(2)*A(3)+A(1)
	rr=-0.5d0*q
	qq=onethird*p
!
!	The following construction is made in order to avoid numeric overflow
!	in evaluation of the discriminant, in a a manner similar to that used
!	in standard routines for evaluation of sqrt(x^2+y^2), for example.
!
	absrr=Cdabs(rr)
	absqq=Cdabs(qq)
	absrr23=absrr**twothird
	if ((absrr23).lt.(absqq)) then
		absqq32=absqq**1.5d0
		cc=rr+absqq32*Cdsqrt((rr/absqq32)**2.d0+(qq/absqq)**3.d0)
	else
		cc=rr+absrr*Cdsqrt((rr/absrr)**2.d0+(qq/absrr23)**3.d0)
	endif
!
!	Evaluate the solutions for ww from the binomial equation ww^3=c, and
!	form the solutions z(k) from the three obtained roots.
!
	cc=cc**onethird
	do 10 k=0,2
		phi=k*twothird*pi
		ww=cc*CMPLX (COS(phi),SIN(phi))
		RTD(k+2)=ww-onethird*(a(3)+p/ww)
10	continue
!
	return
	end subroutine
!
!*******************************************************************************************
!
!*******************************************************************************************
!
	subroutine ROOTSATDYN(SLAP,RTSD)
!
	implicit double precision (a-h,o-z)
	COMPLEX(8) RTSD(4),D1,D3,SLAP,zalphaT
	common /e/ gammaw,gammas,gammaa,grav,atmp
	common /mate/zlamda,zmuy,znuu,znu,zalpha,zrho,zrhof,zm,zk,zr,zn
!
!
	do i=1,4
		RTSD(i)=0.d0
	enddo
!
	RTSD(3)=(zrho/zmuy)*SLAP**2.d0
	RTSD(4)=(zrho/(zlamda+2.d0*zmuy))*SLAP**2.d0
!
	zalphaT=zalpha-SLAP*zrhof*zk
!
	D1=(zn**2.d0)*SLAP/(zk*zr)+zalpha*SLAP*(zalpha-SLAP*zrhof*zk)/((zlamda+2.d0*zmuy)*zk)+&
	   (SLAP**2.d0)*zrho/(zlamda+2.d0*zmuy)
	D3=D1**2.d0-4.d0*(SLAP**2.d0)*(zn**2.d0)*zrho*SLAP/(zr*(zlamda+2.d0*zmuy)*zk)
!
	RTSD(1)=0.5d0*(D1+Cdsqrt(D3))
	RTSD(2)=0.5d0*(D1-Cdsqrt(D3))
!
	return
	end
!
!*******************************************************************************************
!
!*******************************************************************************************
!
	subroutine ROOTSATDYNI(SLAP,RTSDI)
!
	implicit double precision (a-h,o-z)
	COMPLEX(8) RTSDI(2),D1,D3,SLAP,zalphaT,beta
	common /e/ gammaw,gammas,gammaa,grav,atmp
	common /mate/zlamda,zmuy,znuu,znu,zalpha,zrho,zrhof,zm,zk,zr,zn
!
!
	do i=1,2
		RTSDI(i)=0.d0
	enddo
!
	beta=zk*zrhof*(zn*SLAP)**2.d0/(SLAP*zn**2.d0+zk*zn*zrhof*SLAP**2.d0)
!
!	RTSDI(2)=(zrho/zmuy)*SLAP**2.d0
	RTSDI(2)=(zrho-beta*zrhof)*SLAP**2.d0/zmuy
!
	zalphaT=zalpha-SLAP*zrhof*zk
	D1=SLAP*zalpha*zalphaT/( (zlamda+2.d0*zmuy)*zk )
	D3=zrho*SLAP**2.d0/(zlamda+2.d0*zmuy)
!
!	RTSDI(1)=D1+D3
	RTSDI(1)=(zrho+zrhof*(1.d0/beta-2.d0))*SLAP**2.d0/(zlamda+2.d0*zmuy)
!
	return
	end


!*******************************************************************************************
!
	subroutine ROOTQST(SLAP,RTQ)
!
	implicit double precision (a-h,o-z)
	COMPLEX(8) RTQ(2),SLAP,c11,c12,c13,c14,c21,c22,c23,c24,c25,c31,c32,c33,c34,D1,D2,D3,&
			   RAD,SLAP1
	common /e/ gammaw,gammas,gammaa,grav,atmp
	common /umate/ d(3,3),dfs(3),xn,xkw,cww,cwg,xka,cgg,cgw,ev,sat,g1,rs,rw,raa,rmix,xm1,&
				   bt,et,xlambda,xmu,henry
!
!
	do i=1,2
		RTQ(i)=0.d0
	enddo
!
!	xlambda1=xlambda/(xlambda+2.d0*xmu)
!	xmu1=xmu/(xlambda+2.d0*xmu)
!	SLAP1=SLAP*rw*xkw
!	xkw1=xkw/xkw
!	xka1=xka/xkw
!	g11=g1*xlambda
!
	D1=- (	dfs(1)*sat*SLAP/(xkw*(xlambda+2.d0*xmu))+&
			(1.d0-dfs(1))*(1.d0-sat)*SLAP/(xka*(xlambda+2.d0*xmu))+&
			g1*xn*SLAP*(xkw+xka)/(xkw*xka) )
	D3=D1**2.d0-4.d0*g1*xn*SLAP*SLAP/(xkw*xka*(xlambda+2.d0*xmu))
!
!	D1=( dfs(1)*sat*SLAP/(xkw*(zlamda+2.d0*zmu))+&
!	     (1.d0-dfs(1))*(1.d0-sat)*SLAP/(xka*(zlamda+2.d0*zmu))-&
!	     g1*xn*SLAP*(xkw+xka)/(xkw*xka) )
!	D3=D1**2.d0+4.d0*g1*xn*SLAP*SLAP/(xkw*xka*(zlamda+2.d0*zmu))
!
	RTQ(1)=0.5d0*(D1+Cdsqrt(D3))
	RTQ(2)=0.5d0*(D1-Cdsqrt(D3))
!
	return
	end


!*******************************************************************************************
!
	subroutine ROOTS(SLAP,RTD,RTQ,RTSD,RTSDI,RTSQ)
!
	implicit double precision (a-h,o-z)
!
	COMPLEX(8) SLAP,RTD(4),RTQ(2),RTSD(4),RTSDI(2),RTSQ
	common /a/ nnp,ifem,ibem,idyn,ieaq,nr
	common /e/ gammaw,gammas,gammaa,grav,atmp
	common /mate/zlamda,zmuy,znuu,znu,zalpha,zrho,zrhof,zm,zk,zr,zn
	common /umate/ d(3,3),dfs(3),xn,xkw,cww,cwg,xka,cgg,cgw,ev,sat,g1,rs,rw,raa,rmix,xm1,&
				   bt,et,xlambda,xmu,henry
	common /comp/ icomp
!
!
	if (ibem.eq.3) then
		Xnu=xlambda/(2.*(xlambda+xmu))
		if (idyn.eq.2) call ROOTDYN(SLAP,RTD) !dynamic
		if (idyn.eq.1) call ROOTQST(SLAP,RTQ) !quasi-static
	else if (ibem.eq.2) then
		if (idyn.eq.2.and.icomp.eq.1) call ROOTSATDYN(SLAP,RTSD)
		if (idyn.eq.2.and.icomp.eq.0) call ROOTSATDYNI(SLAP,RTSDI)
		if (idyn.eq.1) RTSQ=SLAP*(zalpha**2.D0)*(1.D0-znuu)*(1.D0-2.D0*znu)**2.D0/&
							 ( 2.D0*zk*zmuy*(1.D0-znu)*(znuu-znu) )
	endif
!
	return
	end subroutine
!
!*******************************************************************************************



!
	subroutine SATANAL
!
	implicit double precision (a-h,o-z)
!
	common /mate/zlamda,zmuy,znuu,znu,zalpha,zrho,zrhof,zm,zk,zr,zn
!
	COMPLEX(8) W1,W,Z,GammaZ,SLAP,UY(200),P(200),RTSD(4),DRTSD(2),DI(2),ZCOSH,ZSINH,&
			   UYANAL(200),PANAL(200),FUYANAL(200),FPANAL(200),FNOMEGAUY(200),FNOMEGAP(200)
	dimension ABSUY(200),ABSP(200),DUY(200),DP(200),yUY(200),yP(200),f(200),H1(200)
	! EXTERNAL CCOSH,ZSINH
!
	Y=0.d0
	LL=3.d0
!
	pi=DACOS (-1.d0)
	ntime=200
	dtime=0.0001
!
	epsilon=10.d0**(-10.d0)
	taw=epsilon**(1.d0/(2.d0*ntime))
	W1=CMPLX(0.d0,2.d0*pi/ntime)
	W=CDexp(W1)
!
10	do l=1,ntime
		Z=taw*W**(l-1)
		GammaZ=0.d0
		do i=1,2
			GammaZ=GammaZ+((1.d0-Z)**i)/i
		enddo
		SLAP=GammaZ/dtime
!	do i=1,200
!		W=50.d0*i
!		SLAP=DCMPLX (0.d0,W)
		call ROOTSATDYN (SLAP,RTSD)
		do N=1,2
			DRTSD(N)=Cdsqrt( RTSD(N) )
		enddo
		do k=1,2
			DI(k)=( (zlamda+2.d0*zmuy)*RTSD(k)-zrho*SLAP**2.d0 )/( zalpha*DRTSD(k) )
		enddo
		UY(L)=( DI(2)*SINH(DRTSD(1)*Y)/COSH(DRTSD(1)*LL)-DI(1)*SINH(DRTSD(2)*Y)/COSH(DRTSD(2)*LL) )/&
			  ((zlamda+2.d0*zmuy)*(DI(1)*DRTSD(2)-DI(2)*DRTSD(1)))
		P(L) =DI(1)*DI(2)*(COSH(DRTSD(1)*Y)/COSH(DRTSD(1)*LL) -COSH(DRTSD(2)*Y)/COSH(DRTSD(2)*LL) )/&
			  ((zlamda+2.d0*zmuy)*(DI(1)*DRTSD(2)-DI(2)*DRTSD(1)))
		absUY(L)=Cdabs (UY(L))
		absP(L)=Cdabs (P(L)) 
		
	enddo
!
	call  DFTF (ntime,UY,FUYANAL)
	call  DFTF (ntime,P,FPANAL)

	do l=1,200
		FNOMEGAUY(l)=taw**(-l+1)*FUYANAL(l)/ntime
		FNOMEGAP (l)=taw**(-l+1)*FPANAL(l)/ntime
		DUY(l)=REAL(FNOMEGAUY(l))
		DP(l)=REAL(FNOMEGAP(l))
	enddo
!
	do n=1,200
		do k=1,n
			if ((k-1)*dtime.ge.0.d0) H1(k)=1.d0
			f(k)=H1(k)
			yUY(n)=yUY(n)+DUY(n-k+1)*f(k)
			yP(n)=yP(n)+DP(n-k+1)*f(k)
		enddo
	enddo
!
	return
	end

!****************************************************************************************
!
!					                  SEISMIC subroutine
!
!	This subroutine takes into account imposed displacement by seismic.
!
!	this subroutine is called in this section: CONTROL
!
!
!	variables used are:
!
!
!		r		:
!		neaq	: number of nodes on bedrock
!		nneaq	: node numbers of nodes on bedrock
!		uxeaq	: bedrock motion in x direction
!		uyeaq	: bedrock motion in y direction
!		id		:
!		na		:
!		s1		:
!		s2		:
!		s11		:
!		s22		:
!		iter	:
!
!
!	INPUT	:
!
!
!	OUTPUT	:
!
!
!****************************************************************************************
!
	subroutine seismic(id,na,s1,s2,s11,s22,r,iter)
!
	implicit double precision (a-h,o-z)
!
	common /a/ nnp,ifem,ibem,idyn,ieaq,nr
	common /chrg3/ nstime,neaq,nneaq(500),uxeaq,uyeaq,axeaq,ayeaq,uxeq(3,500),uyeq(3,500)
	common /h/ mdof,mdofn,ls,ls1,isolv
!
	dimension s1(ls),s2(ls1),id(4,nnp),r(mdof),na(mdof),s11(ls),s22(ls1)
!
!
	if (iter.eq.1) then
		do i=1,ls
			s11(i)=s1(i)
		enddo
		do i=1,ls1
			s22(i)=s2(i)
		enddo
	else
		do i=1,ls
			s1(i)=s11(i)
		enddo
		do i=1,ls1
			s2(i)=s22(i)
		enddo
	endif
!
	do 30 l=1,neaq
		i=nneaq(l)
		do 20 ii=1,2
			j=id(ii,i)
			if (j.gt.mdof) goto 20
			nrowj=1
			if (j.gt.1) nrowj=j+1-(na(j)-na(j-1))
        do 10 k=nrowj,mdof
          if (k.lt.j) then
            lpl=k-j+na(j)
            if (ii.eq.1) r(k)=r(k)-s1(lpl)*uxeaq
            if (ii.eq.2) r(k)=r(k)-s1(lpl)*uyeaq
            s1(lpl)=0.d0
            if (isolv.eq.3) s2(lpl)=0.0d0
          elseif (k.eq.j) then
            s1(na(k))=1.d0
            if (isolv.eq.3) s2(na(k))=1.d0
            if (ii.eq.1) r(k)=uxeaq
            if (ii.eq.2) r(k)=uyeaq
          else
            nrowk=k+1-(na(k)-na(k-1))
            if (nrowk.le.j) then
              lpl=j-k+na(k)
              if (isolv.eq.1) then
                if (ii.eq.1) r(k)=r(k)-s1(lpl)*uxeaq
                if (ii.eq.2) r(k)=r(k)-s1(lpl)*uyeaq
                s1(lpl)=0.d0
              else
                if (ii.eq.1) r(k)=r(k)-s2(lpl)*uxeaq
                if (ii.eq.2) r(k)=r(k)-s2(lpl)*uyeaq
                s1(lpl)=0.d0 ; s2(lpl)=0.d0
              endif
            endif
          endif
10      continue
20    continue
30    continue
!
      return
      end


!****************************************************************************************
!
!							           SHAP8N subroutine
!
!	This subroutine calculates the shape functions and their derivatives.
!									3---2---1
!
!	this subroutine is called in: STR8ND
!
!	variables used are:
!
!	INPUT	:
!
!	OUTPUT	:
!
!
!****************************************************************************************
!
	subroutine shap3n(s,f,pfs)
!
 	implicit double precision (a-h,o-z)
!
	dimension f(3),pfs(3)
!
	sp=1.-s
	sm=1.+s
!
!	h functions
!
	h1=0.5*sm
	h2=sp*sm
	h3=0.5*sp
!
!	h derivatives to s
!
	h1s=0.5
	h2s=-2.*s
	h3s=-0.5
!
!	shape functions and their derivatives
!
!	f(1) to f(3)
!
	f(1)=h1-0.5*h2
	f(2)=h2
	f(3)=h3-0.5*h2
!
!	pfs(1) to pfs(3)
!
	pfs(1)=h1s-0.5*h2s
	pfs(2)=h2s
	pfs(3)=h3s-0.5*h2s
!
	return
	end


!****************************************************************************************
!
!							           SHAP8N subroutine
!
!	This subroutine calculates the shape functions and their derivatives.
!									3---2---1
!									---------
!									4-------8
!									---------
!									5---6---7
!
!	this subroutine is called in: STR8ND
!
!	variables used are:
!
!	INPUT	:
!
!	OUTPUT	:
!
!
!****************************************************************************************
!
	subroutine shap8n(s,t,f,pfs,pft)
!
 	implicit double precision (a-h,o-z)
!
	dimension f(8),pfs(8),pft(8)
!
	sp=1.-s
	sm=1.+s
	tp=1.-t
	tm=1.+t
!
!	h functions
!
	h1=0.25*sm*tm
	h2=0.5*sp*sm*tm
	h3=0.25*sp*tm
	h4=0.5*tp*tm*sp
	h5=0.25*sp*tp
	h6=0.5*sp*sm*tp
	h7=0.25*sm*tp
	h8=0.5*tp*tm*sm
!
!	h derivatives to s
!
	h1s=0.25*tm
	h2s=-s*tm
	h3s=-0.25*tm
	h4s=-0.5*tp*tm
	h5s=-0.25*tp
	h6s=-s*tp
	h7s=0.25*tp
	h8s=0.5*tp*tm
!
!	h derivatives to t
!
	h1t=0.25*sm
	h2t=0.5*sp*sm
	h3t=0.25*sp
	h4t=-t*sp
	h5t=-0.25*sp
	h6t=-0.5*sp*sm
	h7t=-0.25*sm
	h8t=-t*sm
!
!	shape functions and their derivatives
!
!	f(1) to f(8)
!
	f(1)=h1-0.5*h8-0.5*h2
	f(2)=h2
	f(3)=h3-0.5*h2-0.5*h4
	f(4)=h4
	f(5)=h5-0.5*h4-0.5*h6
	f(6)=h6
	f(7)=h7-0.5*h6-0.5*h8
	f(8)=h8
!
!	pfs(1) to pfs(8)
!
	pfs(1)=h1s-0.5*h8s-0.5*h2s
	pfs(2)=h2s
	pfs(3)=h3s-0.5*h2s-0.5*h4s
	pfs(4)=h4s
	pfs(5)=h5s-0.5*h4s-0.5*h6s
	pfs(6)=h6s
	pfs(7)=h7s-0.5*h6s-0.5*h8s
	pfs(8)=h8s
!
!	pft(1) to pft(8)
!
	pft(1)=h1t-0.5*h8t-0.5*h2t
	pft(2)=h2t
	pft(3)=h3t-0.5*h2t-0.5*h4t
    pft(4)=h4t
    pft(5)=h5t-0.5*h4t-0.5*h6t
    pft(6)=h6t
    pft(7)=h7t-0.5*h6t-0.5*h8t
    pft(8)=h8t
!
	return
	end

!****************************************************************************************
!
!							            SNGTP subroutine
!
!	This subroutine evaluates the element's kinds of singularity & the number of sub-elm
!	the error estimate that is introduced in this subroutine g(1:8) allows us to ensure
!	that the error made by numerical integration is nearly constant, regardless of the
!	proximity of point Pi.
!   In some cases, when point Pi is near the element, the number of Gauss points required
!	will exceed 4 in table 6.1. In this case it is necessary to subdivide the element into
!	sub regions of integration. A simple approach is to subdivide the element into equal
!	subdivisions depending on the value of R/L. If according to the R/L value the maximum
!	number of Gauss points available is exceeded, the element is subdivided into K regions
!	where K = INT [(R/L)min / (R/L) ]
!
!	This subroutine is called in: GHMATD1,GHMATD2,GHMATD3
!	This subroutine calls: -
!
!	variables used are:
!
!		eps		: upper bound of error in the Stroud and Secrest method = 1.d-16
!		j		: counter which shows the number of Gauss points between 3 & 10
!		g		: Stroud and Secrest method; R/L for 3-10 Gauss points
!		itime	: time increment; [1]: first time increment, [>1]: higher time increment
!		xp		: x-component of collocation point (loading point)
!		yp		: y-component of collocation point (loading point)
!		xelm	: x-component of nodes forming the each element in each D1 (finite) sub-zone
!		yelm	: y-component of nodes forming the each element in each D1 (finite) sub-zone
!		r1		: distance between the collocation point & the first node of element
!		r2		: distance between the collocation point & the second node of element
!		r3		: distance between the collocation point & the third node of element
!		rmin	: minimum distance between the collocation point & the nodes of element
!		rl		: length of element (L)
!		n		: nnpbed1
!		ll		: counter which shows the coolocation point (loading point) ll=1,nnpbed1
!		i		: counter which shows the element node (studied point), i=1,nnpbed1-1,2
!		nodo	: code of singularity of element [0,1,2,3]
!									[0]: no singularity, loading point is outside of the
!										 boundary element studied
!									[1]: singularity in the first node of element
!									[2]: singularity in the second node of element
!									[3]: singularity in the third node of element
!		nsub	: number of subelements
!		nptg	: number of Gauss points in each sub-region which is determined by the
!				  Stroud and Secrest method
!
!
!	INPUT	: n,ll,i
!
!	OUTPUT	: nodo,nsub,nptg
!
!
!****************************************************************************************
!
	subroutine sngtp(n,ll,i,Iint,nodo,nsub,nptg)
!
	implicit double precision (a-h,o-z)
!
	common /a/ nnp,ifem,ibem,idyn,ieaq,nr
	common /c/ nload,ntime,itmax,itime,iter,iload,time,dtime,dtimei,dtimem,xinc
	common /jab/ajab,bjab,cjab,djab,xelm(3),yelm(3),xp,yp
!
	dimension g(8)
	save g,eps
!
	if (eps.ne.0.d0) goto 10
	eps=1.d-16
	do j =3,10
		g(j-2)=4*((eps/2)**(0.5d0/j))
	enddo
10	continue
!
	if (itime.eq.1) goto 40
	nodo=0
	nsub=1
	nptg=6
	goto 50
!
40	r1=dsqrt ((xp-xelm(1))**2+(yp-yelm(1))**2)
	r2=dsqrt ((xp-xelm(2))**2+(yp-yelm(2))**2)
	r3=dsqrt ((xp-xelm(3))**2+(yp-yelm(3))**2)
	rmin=DMIN1 (r1,r2,r3)
	rl=dsqrt ((xelm(3)-xelm(1))**2+(yelm(3)-yelm(1))**2)
!
!
!	------------- singularity of element -------------
!
	nodo = 0					! loading point is outside of the boundary element studied
!
!	loading point is inside of the boundary element studied:
	if (ll.eq.i)     nodo=1		! singularity in the first node of element
	if (ll.eq.(i+1)) nodo=2		! singularity in the second node of element
	if (ll.eq.(i+2)) nodo=3		! singularity in the third node of element
	if (ll.eq.1 .and. i.eq.n-1) nodo=3	! in the last element singularity occures in the
!										  third node
	if (Iint.eq.1) nodo = 0
!
!
!	-------- evaluating automatic number of subelement --------
!
	if (nodo.eq.0) then ! no singularity, when point Pi is not one of the elm. nodes
		ratio=rl/rmin
	else if (nodo.eq.2) then	! singularity in the midside node of element
		ratio=2.d0
	else	! singularity in the side nodes of element
        ratio=1.d0
	endif
!
	nsub=INT (ratio/g(8))+1
!
!
!	-------- evaluating automatic number of points Gauss --------
!
	ratio=ratio/nsub
	do j=3,10
		if (ratio.lt.g(j-2)) then
			nptg=j
			goto 50
		endif
	enddo
50	continue
!
	return
	end




!****************************************************************************************
!
!							           SOLVE subroutine
!
!	This subroutine solves for the global displacements and pore pressures:
!		skyline solver for symmetric matrix by MONDKAR & POWEL:	(isolv=1)
!		skyline solver for unsymmetric matrix by TAYLOR:		(isolv=3)
!	This subroutine is called in: CONTROL
!	This subroutine calls: -
!
!	variables used are:
!
!		als		: symmetric bem stiffness matrix (S1)
!		als1	: unsymmetric bem stiffness matrix (S2)
!		isolv	: solution type indicator [1,3]
!							[1]: for a symmetric computation
!							[3]: for a non-symmetric computation
!		b		: R, all nodal forces (mechanical TOTRR & seismic RRR) are saved in this
!				  vector; after that in SOLVE subroutine the displacements & pressures
!				  obtained by solve the equations are saved in this vector in each itime
!
!
!	INPUT	:
!
!	OUTPUT	: b
!
!
!****************************************************************************************
!
	subroutine solve(als,als1,b,na)
!
	implicit double precision (a-h,o-z)
!
	common /h/ mdof,mdofn,ls,ls1,isolv
!
	dimension a(ls),b(1),c(ls1),na(1),als(1),als1(1)
!
!
	neq=mdof
	do i=1,ls
		a(i)=als(i)
	enddo
	do i=1,ls1
		c(i)=als1(i)
	enddo
!
	goto (1000,2000,3000) isolv
!
1000 continue
!
!
!	---------------------------- reduce coefficient matrix a ----------------------------
!
      neqq=neq-1
      il=1
      najp=na(1)
      do 7 j=2,neq
      naj=na(j)
      jk=naj-j
      if=1-jk+najp
      if (if.ge.j) goto 6
      if1=if+1
      if(if1.gt.il) goto 4
      jia=2+najp
      naip=na(if)
      kl=if+jk
      do 3 i=if1,il
      nai=na(i)
      ik=nai-i
      ii=i+1-nai+naip
      if(ii.ge.i) goto 2
      kf=max0(ii,if)+jk
      jj=ik-jk
      aa=0.0
      do 1 k=kf,kl
1     aa=aa+a(k)*a(jj+k)
      a(jia)=a(jia)-aa
2     jia=jia+1
      kl=kl+1
3     naip=nai
4     kf=jk+if
      kl=naj-1
      aa=0.0
      do 5 k=kf,kl
      nai=na(if)
      cc=a(k)/a(nai)
      aa=aa+a(k)*cc
      a(k)=cc
5     if=if+1
      a(naj)=a(naj)-aa
6     il=il+1
7     najp=naj
!
!     reduce vector b and back substitute
!
      do 8 n=1,neqq
      if(b(n).ne.0.0) goto 9
8     continue
      n=neqq
9     n1=n+1
      i1=n1+1
      kl=n
      naip=na(n)
      do 12 i=n1,neq
      nai=na(i)
      ii=i1-nai+naip
      if(ii.ge.i) goto 11
      kf=max0(ii,n)
      ik=nai-i
      ika=ik+kf
      bb=0.0
      do 10 k=kf,kl
      bb=bb+a(ika)*b(k)
10    ika=ika+1
      b(i)=b(i)-bb
11    i1=i1+1
      kl=kl+1
12    naip=nai
      do 13 i=n,neq
      nai=na(i)
13    b(i)=b(i)/a(nai)
      j=neq
      naj=na(neq)
      do 16 i=1,neqq
      najp=na(j-1)
      jka=najp+1
      ii=j-naj+jka
      if(ii.ge.j) goto 15
      kl=j-1
      bb=b(j)
      do 14 k=ii,kl
      b(k)=b(k)-a(jka)*bb
14    jka=jka+1
15    j=j-1
16    naj=najp
      return
!
!
!
2000 continue
	goto 4000
3000 continue
!
!
!	---------------------------- unsymmetric matrix ----------------------------
!
	jr=0
	do 300 j=1,neq
		jd=na(j)
		jh=jd-jr
		if (jh.le.1) goto 300
		is=j+1-jh
		ie=j-1
		k=jr+1
		id=0
!
!	reduce all equations except diagonal
!
		do 200 i=is,ie
			ir=id
			id=na(i)
			ih= MIN0 (id-ir-1,i-is)
			if (ih.eq.0) goto 150
			a(k)=a(k)-dot(a(k-ih),c(id-ih),ih)
			c(k)=c(k)-dot(c(k-ih),a(id-ih),ih)
150			if (a(id).ne.0.0) c(k)=c(k)/a(id)
200			k=k+1
!
!	reduce diagonal term
!
		a(jd)=a(jd)-dot(a(jr+1),c(jr+1),jh-1)
!
!	forward reduce the r.h.s
!
		b(j)=b(j)-dot(c(jr+1),b(is),jh-1)
300	jr=jd
!
!
!	back substitution
!
	j=neq
	jd=na(j)
500	if (a(jd).ne.0.0) b(j)=b(j)/a(jd)
	d=b(j)
	j=j-1
	if (j.le.0) return
	jr=na(j)
	if (jd-jr.le.1) goto 700
	is=j-jd+jr+2
	k=jr-is+1
	do 600 i=is,j
600		b(i)=b(i)-a(i+k)*d
!
700	jd=jr
	goto 500
4000 continue
!
	return
	end
!
!
!
!	************************ vector dot product ************************
!
	function dot(a,b,n)
!
	implicit double precision (a-h,o-z)
!
	dimension a(1),b(1)
!
	dot=0.0d0
	do i=1,n
		dot=dot+a(i)*b(i)
	enddo
!
	return
	end


!****************************************************************************************
!
!							           STR3NU subroutine
!
!	This subroutine calculates stresses for unsaturated elements at the centroid.
!
!	this subroutine is called in: FSTRESS
!
!	variables used are:
!
!	INPUT	:
!
!	OUTPUT	:
!
!
!****************************************************************************************
!
	subroutine str3nu(m,mt,trac,xe,ye,u,pw,pa,dpw,dpa,sige,sm8u)
!
 	implicit double precision (a-h,o-z)
!
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /a3/ ne8d1,ne8c1,ne8u1,ne3u1,nload1,nbcx1,nbcy1,nbcw1,nbca1,m8d1,m8c1,m8u1,&
				nbed11,nbed21,nbed31,nbeint1,ne8dd,ne8cc,ne8uu
	common /b/ m8d,m8c,m8u,m8dep(10),m8cep(10),m8uep(10)
	common /c/ nload,ntime,itmax,itime,iter,iload,time,dtime,dtimei,dtimem,xinc
	common /chrg2/ npw,npa,nnpw(200),nnpa(200),vnpw(200),vnpa(200)
	common /g/ init,iprint
	common /ini2/ npwi,npai,ndispi,nveli,nacci
	common /umate/ d(3,3),dfs(3),xn,xkw,cww,cwg,xka,cgg,cgw,ev,sat,g1,rs,rw,raa,rmix,xm1,&
				   bt,et,xlambda,xmu,henry
!
	dimension xe(3),ye(3),u(6),pw(3),pa(3),dpw(3),dpa(3),sige(5),sm8u(60,m8u1),&
			  f(3),pfs(3),b(3,6),epi(3),trac(6)
!
!
	imat=m8uep(mt)
	zmu=3.d0*bt*et/(9.d0*bt-et)
	zlamda=bt-2.d0*zmu/3.d0
	znu=zlamda/(2*(zlamda+zmu))
!
	s=0.0
	call shap3n(s,f,pfs)
	call bmat3n(xe,ye,pfs,b,detj)
!
!	calculate (incremental) strains
!
	epix=0.0d0
	do i=1,3
		epix=epix+( pfs(i)*u(2*i-1) )/detj
	enddo
!
!	calculate stresses
!
!	linear elastic material
!
	if (imat.eq.1) then
		temp0=0.d0
		temp1=0.d0
		do i=1,3
			temp0=temp0+f(i)*trac(2*i-1)
			temp1=temp1+f(i)*trac(2*i)
		enddo
		sige(2)=temp1 !sigy
		sige(3)=temp0 !sigxy
		sige(1)=(2.d0*zmu*epix+znu*sige(2))/(1.d0-znu)
	endif
!
!	calculate pore pressure and suction at centroid
!
	temp2=0.0d0
	temp3=0.0d0
	do i=1,3
		temp2=temp2+f(i)*pw(i)
		temp3=temp3+f(i)*pa(i)
	enddo
	suc=temp3-temp2
	sige(4)=suc
	sige(5)=temp3
!
	return
	end

!****************************************************************************************
!
!							           STR8ND subroutine
!
!	This subroutine calculates stresses for consolidated elements at point Gauss.
!
!	this subroutine is called in: FSTRESS
!
!	variables used are:
!
!	INPUT	:
!
!	OUTPUT	:
!
!
!****************************************************************************************
!
	subroutine str8nc(ielc,mt,xe,ye,du,u,pw,dpw,sige,sigei,sigma,sigo,sm8d,sm8c,dpac,&
					  alfac,iysc)
!
 	implicit double precision (a-h,o-z)
!
	common /a3/ ne8d1,ne8c1,ne8u1,ne3u1,nload1,nbcx1,nbcy1,nbcw1,nbca1,m8d1,m8c1,m8u1,&
				nbed11,nbed21,nbed31,nbeint1,ne8dd,ne8cc,ne8uu
	common /b/ m8d,m8c,m8u,m8dep(10),m8cep(10),m8uep(10)
	common /c/ nload,ntime,itmax,itime,iter,iload,time,dtime,dtimei,dtimem,xinc
	common /chrg2/ npw,npa,nnpw(200),nnpa(200),vnpw(200),vnpa(200)
	common /g/ init,iprint
	common /ini2/ npwi,npai,ndispi,nveli,nacci
!
	dimension xe(8),ye(8),du(16),u(16),pw(8),dpw(8),ss(9),tt(9),sige(9,6),sigei(9,6),&
			  sigma(9,6),ssige(6),ssigma(6),sm8d(60,m8d1),sm8c(60,m8c1),sigo(9,6),&
			  ssigo(6),dpac(3,3,ne8cc),alfac(4,10,ne8cc),iysc(2,ne8cc),f(8),pfs(8),&
			  pft(8),b(3,16),d(3,3),depi(3),epi(3)
!
	data ss/0.77459,0.0,-0.77459,0.77459,0.0,-0.77459,0.77459,0.0,-0.77459/
	data tt/0.77459,0.77459,0.77459,0.0,0.0,0.0,-0.77459,-0.77459,-0.77459/
!
!
	imat=m8cep(mt)
	if (imat.eq.2 .OR. imat.eq.3) then
		phi=sm8c(1,mt)
		c=sm8c(2,mt)
	endif
!
	do 100 ii=1,9
		s=ss(ii)
		t=tt(ii)
		call shap8n(s,t,f,pfs,pft)
		call bmat8n(xe,ye,pfs,pft,b,detj)
!
!	calculate (incremental) strains
!
		do i=1,3
			depi(i)=0.0d0
			epi(i)=0.0d0
			do j=1,16
				depi(i)=depi(i)+b(i,j)*du(j)
				epi(i)=epi(i)+b(i,j)*u(j)
			enddo
		enddo
!
!	calculate stresses
!
!	linear elastic material
!
	if (imat.eq.1) then
		call dmatl(2,mt,d,sm8d,sm8c)
		do i=1,3
			temp=0.d0
			do j=1,3
				temp=temp+d(i,j)*depi(j)
			enddo
			sige(ii,i)=sigma(ii,i)+temp
		enddo
	endif
!
!	hyperbolic elastic material for statics
!
	if (imat.eq.2) then
		do i=1,6
			ssige(i)=sige(ii,i)
			ssigma(i)=sigma(ii,i)
		enddo
        if (init.gt.0) then
			call dmath1(2,mt,d,sm8d,sm8c,ssige)
			do i=1,3
				temp=0.d0
				do j=1,3
					temp=temp+d(i,j)*depi(j)
				enddo
				sige(ii,i)=sigma(ii,i)+temp
			enddo
		else
			call stresshyp1(2,mt,sm8d,sm8c,phi,c,ssigma,ssige,depi)
			do i=1,6
				sige(ii,i)=ssige(i)
			enddo
		endif
	endif
!
!	hyperbolic elastic material for dynamics
!
	if (imat.eq.3) then
		do i=1,6
			ssige(i)=sige(ii,i)
			ssigma(i)=sigma(ii,i)
			ssigo(i)=sigo(ii,i)
		enddo
		if (init.gt.0) then
			call dmath2(2,0,mt,d,sm8d,sm8c,ssige,epi)
			do i=1,3
				temp=0.d0
				do j=1,3
					temp=temp+d(i,j)*epi(j)
				enddo
				sige(ii,i)=temp
			enddo
		else
			call stresshyp3(2,mt,sm8d,sm8c,ssige,ssigma,ssigo,epi,depi)
			do i=1,6
				sige(ii,i)=ssige(i)
			enddo
		endif
	endif
!
!	if (imat.eq.4) then
!        idmat=0
!        call dmatcp(ielc,mt,idmat,d,sige,dsige,sigei,sm8c,depi,
!     $    dpac,alfac,iysc,xkw,comp,sat,coef)
!      endif
!
!	calculate pore pressure
!
	if (npwi.eq.0) goto 10
	if (itime.gt.1) goto 10
	temp=0.0d0
	do i=1,8
		temp=temp+f(i)*(pw(i)+dpw(i))
	enddo
	sige(ii,6)=temp
	goto 100
!
10	if (npw.eq.0) goto 20
	if (npwi.ne.0) goto 20
	if (itime.eq.1) then
		temp=0.0d0
          do  i=1,8
            temp=temp+f(i)*(pw(i)+dpw(i))
          enddo
          sige(ii,6)=temp
          goto 100
        endif
20    temp=0.0d0
      do i=1,8
        temp=temp+f(i)*pw(i)
      enddo
      sige(ii,6)=temp
100   continue
!
	return
	end

!****************************************************************************************
!
!							           STR8ND subroutine
!
!	This subroutine calculates stresses for drained elements at point Gauss.
!
!	this subroutine is called in: FSTRESS
!
!	variables used are:
!
!	INPUT	:
!
!	OUTPUT	:
!
!
!****************************************************************************************
!
	subroutine str8nd(ield,mt,xe,ye,du,u,sige,sigei,sigma,sigo,sm8d,sm8c,dpad,alfad,iysd)
!
 	implicit double precision (a-h,o-z)
!
	common /a3/ ne8d1,ne8c1,ne8u1,ne3u1,nload1,nbcx1,nbcy1,nbcw1,nbca1,m8d1,m8c1,m8u1,&
				nbed11,nbed21,nbed31,nbeint1,ne8dd,ne8cc,ne8uu
	common /b/ m8d,m8c,m8u,m8dep(10),m8cep(10),m8uep(10)
	common /c/ nload,ntime,itmax,itime,iter,iload,time,dtime,dtimei,dtimem,xinc
	common /g/ init,iprint
!
    dimension xe(8),ye(8),du(16),u(16),ss(9),tt(9),&
			  sige(9,6),sigei(9,6),sigma(9,6),ssige(6),ssigma(6),&
			  sm8d(60,m8d1),sm8c(60,m8c1),sigo(9,6),ssigo(6),&
			  dpad(3,3,ne8dd),alfad(4,10,ne8dd),iysd(2,ne8dd),&
			  f(8),pfs(8),pft(8),b(3,16),d(3,3),depi(3),epi(3)
!
	data ss/0.77459,0.0,-0.77459,0.77459,0.0,-0.77459,0.77459,0.0,-0.77459/
	data tt/0.77459,0.77459,0.77459,0.0,0.0,0.0,-0.77459,-0.77459,-0.77459/
!
!
	imat=m8dep(mt)
	if (imat.eq.2 .OR. imat.eq.3) then
		phi=sm8d(1,mt)
		c=sm8d(2,mt)
	endif
!
	do 550 ii=1,9
		s=ss(ii)
		t=tt(ii)
		call shap8n(s,t,f,pfs,pft)
		call bmat8n(xe,ye,pfs,pft,b,detj)
!
!     calculate (incremental) strains
!
	do i=1,3
		depi(i)=0.d0
		epi(i)=0.d0
		do j=1,16
			depi(i)=depi(i)+b(i,j)*du(j)
			epi(i)=epi(i)+b(i,j)*u(j)
		enddo
	enddo
!
!     calculate stresses
!
!     linear elastic material
!
	if (imat.eq.1) then
		call dmatl(1,mt,d,sm8d,sm8c)
		do i=1,3
			temp=0.d0
			do j=1,3
				temp=temp+d(i,j)*depi(j)
			enddo
			sige(ii,i)=sigma(ii,i)+temp
		enddo
	endif
!
!    hyperbolic elastic material for statics	????
!
	if (imat.eq.2) then
		do i=1,5
			ssige(i)=sige(ii,i)
			ssigma(i)=sigma(ii,i)
		enddo
		if (init.gt.0) then
			call dmath1(1,mt,d,sm8d,sm8c,ssige)
			do i=1,3
				temp=0.d0
				do j=1,3
					temp=temp+d(i,j)*depi(j)
				enddo
				sige(ii,i)=sigma(ii,i)+temp
			enddo
		else
			call stresshyp1(1,mt,sm8d,sm8c,phi,c,ssigma,ssige,depi)
			do i=1,5
				sige(ii,i)=ssige(i)
			enddo
		endif
	endif
!
!     hyperbolic elastic material for dynamics	????
!
	if (imat.eq.3) then
		do i=1,5
			ssige(i)=sige(ii,i)
			ssigma(i)=sigma(ii,i)
			ssigo(i)=sigo(ii,i)
		enddo
		if (init.gt.0) then
			call dmath2(1,0,mt,d,sm8d,sm8c,ssige,epi)
			do i=1,3
				temp=0.d0
				do j=1,3
					temp=temp+d(i,j)*depi(j)
				enddo
				sige(ii,i)=sigma(ii,i)+temp
			enddo
		else
			call stresshyp2(1,mt,sm8d,sm8c,ssige,ssigma,ssigo,epi,depi)
			do i=1,5
				sige(ii,i)=ssige(i)
			enddo
		endif
	endif
!
!     elastoplastic material 	????
!
!	if (imat.eq.4) then
!		idmat=0
!		call dmatdp(ield,mt,idmat,d,sige,dsige,sigei,sm8d,depi,dpad,alfad,iysd)
!	endif
!
550	continue
!
	return
	end

!!!!!!c*************************************!!!!!!c*********************************

!!!c INTEGRATION ADAPTIVE DE LA HYPERBOLIC CONSTITUTIVE RELATION
!!!!!c              POUR CALCULER LA CONTRAITE

!!!!!c  c'est le sous-programme de 'str8nd' et 'str8nc'

!!!!!!c*************************************!!!!!!c*********************************

      subroutine stresshyp1(kode,mt,sm8d,sm8c,phi,c,sigma,sige,depi)


      implicit double precision (a-h,o-z)
      common /a3/ne8d1,ne8c1,nload1,nbcx1,nbcy1,nbcw1,m8d1,m8c1,nbed11,nbed21,nbed31,ne8dd,ne8cc
      dimension sigma(6),sige(6),depi(3),d(3,3),dsig1(3),dsig2(3),sig1(6),sig2(3),ddep(3),erreur(3),sm8d(60,m8d1),sm8c(60,m8c1)

      toll=1.d-3
      do j=1,3
        sige(j)=sigma(j)
      enddo
      t=0.d0 ; dt=1.d0

10    if (dt.gt.(1.d0-t)) dt=1.d0-t

!!c pr�diction sigma1 par Euler forward
      call dmath1(kode,mt,d,sm8d,sm8c,sige)
      do j=1,3
        ddep(j)=dt*depi(j)
      enddo
      do j=1,3
        dsig1(j)=0.d0
        do k=1,3
          dsig1(j)=dsig1(j)+d(j,k)*ddep(k)
        enddo
        sig1(j)=sige(j)+dsig1(j)
      enddo
      sig1(4)=sige(4)
!!!csig1(5)=sige(5)
!!!ccall princp(-sig1(1),-sig1(2),-sig1(3),s1,s3)
!!!csr=(1.d0-dsin(phi))*(s1-s3)/2.d0
!!c $    /(c*dcos(phi)+s3*dsin(phi))
!!!cif (sr.gt.sige(4)) sig1(4)=sr

!!c correction par r�gle de trap�ze
      rerr=0.d0 ; rsig=0.d0
      call dmath1(kode,mt,d,sm8d,sm8c,sig1)
      do j=1,3
        dsig2(j)=0.d0
        do k=1,3
          dsig2(j)=dsig2(j)+d(j,k)*ddep(k)
        enddo
        sig2(j)=sige(j)+(dsig1(j)+dsig2(j))/2.d0
        erreur(j)=(dsig2(j)-dsig1(j))/2.d0
        rsig=rsig+sig2(j)**2
        rerr=rerr+erreur(j)**2
      enddo

!!c verification
      r=dsqrt(rerr)/dsqrt(rsig)
      if (r.le.toll) goto 20
        bt=0.8*dsqrt(toll/r)
        dt=dt*bt
        goto 10
20    do j=1,3
        sige(j)=sig2(j)
      enddo
!!!ccall princp(-sige(1),-sige(2),-sige(3),s1,s3)
!!!csr=(1.d0-dsin(phi))*(s1-s3)/2.d0
!!c $    /(c*dcos(phi)+s3*dsin(phi))
!!!cif (sr.gt.sige(4)) sige(4)=sr
      t=t+dt
      if (t.lt.1.d0) goto 10

30    call princp(-sige(1),-sige(2),-sige(3),s1,s3)
      sr=(1.d0-dsin(phi))*(s1-s3)/2.d0/(c*dcos(phi)+s3*dsin(phi))
!!!cif (     (sr.ge.sige(4)*0.99 .and. sige(5).eq.0.d0)
!!c $    .or. (sr.ge.sige(4)*1.01 .and. sige(5).eq.1.d0)) then
!!!!c sige(5)=0.d0
!!!!c sige(4)=sr
         if (sr.gt.sige(4)) sige(4)=sr
!!!celse
!!!!c sige(5)=1.d0
!!!cendif

      return
      end


!!!!!!c*************************************!!!!!!c*********************************

!!!!!c                CALCULER LA CONTRAITE
!!!!!c         c'est le sous-programme de 'str8nd'

!!!!!!c*************************************!!!!!!c*********************************

      subroutine stresshyp2(kode,mt,sm8d,sm8c,sige,sigma,sigo,epi,depi)

      implicit double precision (a-h,o-z)
      common /a3/ne8d1,ne8c1,nload1,nbcx1,nbcy1,nbcw1,m8d1,m8c1,nbed11,nbed21,nbed31,ne8dd,ne8cc
      common /e/gammaw,grav,atmp
	  dimension sige(6),sigma(6),sigo(6),epi(3),depi(3),d(3,3),sm8d(60,m8d1),sm8c(60,m8c1),epi1(3),dsig1(3),dsig2(3), &
	  			sig1(6),sig2(3),ddep(3),erreur(3)

!!!ccall dmath2(kode,1,mt,d,sm8d,sm8c,sigo,epi)
      call dmath2(kode,3,mt,d,sm8d,sm8c,sigo,epi)
      do i=1,3
        temp=0.d0
        do j=1,3
          temp=temp+d(i,j)*epi(j)
        enddo
        sige(i)=sigo(i)+temp
      enddo
      goto 100

      toll=1.d-4
      do j=1,3
        sige(j)=sigma(j)
        epi(j)=epi(j)-depi(j)
      enddo
      t=0.d0 ; dt=1.d0

10    if (dt.gt.(1.d0-t)) dt=1.d0-t

!!c pr�diction sigma1 par Euler forward
      call dmath2(kode,2,mt,d,sm8d,sm8c,sige,epi)
      do j=1,3
        ddep(j)=dt*depi(j)
      enddo
      do j=1,3
        dsig1(j)=0.d0
        do k=1,3
          dsig1(j)=dsig1(j)+d(j,k)*ddep(k)
        enddo
        sig1(j)=sige(j)+dsig1(j)
        epi1(j)=epi(j)+ddep(j)
      enddo

!!c correction par r�gle de trap�ze
      rerr=0.d0 ; rsig=0.d0
      call dmath2(kode,2,mt,d,sm8d,sm8c,sig1,epi1)
      do j=1,3
        dsig2(j)=0.d0
        do k=1,3
          dsig2(j)=dsig2(j)+d(j,k)*ddep(k)
        enddo
        sig2(j)=sige(j)+(dsig1(j)+dsig2(j))/2.d0
        erreur(j)=(dsig2(j)-dsig1(j))/2.d0
        rsig=rsig+sig2(j)**2
        rerr=rerr+erreur(j)**2
      enddo

!!c verification
      r=dsqrt(rerr)/dsqrt(rsig)
      if (r.le.toll) goto 20
        bt=0.8*dsqrt(toll/r)
        dt=dt*bt
        goto 10
20    do j=1,3
        sige(j)=sig2(j)
        epi(j)=epi1(j)
      enddo
      t=t+dt
      if (t.lt.1.d0) goto 10

100   return
      end


!!!!!!c*************************************!!!!!!c*********************************

!!!c INTEGRATION DE LA HYPERBOLIC CONSTITUTIVE RELATION
!!!!!c              POUR CALCULER LA CONTRAITE

!!!!!c        c'est le sous-programme de 'str8nc'

!!!!!!c*************************************!!!!!!c*********************************

      subroutine stresshyp3(kode,mt,sm8d,sm8c,sige,sigma,sigo,epi,depi)


      implicit double precision (a-h,o-z)
      common /a3/ne8d1,ne8c1,nload1,nbcx1,nbcy1,nbcw1,m8d1,m8c1,nbed11,nbed21,nbed31,ne8dd,ne8cc
      dimension sigma(6),sige(6),sigo(6),epi(3),depi(3),d(3,3),sm8d(60,m8d1),sm8c(60,m8c1)

      call dmath2(kode,3,mt,d,sm8d,sm8c,sige,epi)
      do i=1,3
      sige(i)=sigo(i)
      do j=1,3
        sige(i)=sige(i)+d(i,j)*epi(j)
      enddo
      enddo
      goto 10

      n=50
      do i=1,3
        sige(i)=sigma(i)
        epi(i)=epi(i)-depi(i)
      enddo
      do i=1,3
        depi(i)=depi(i)/n
      enddo

      do i=1,n
        call dmath2(kode,2,mt,d,sm8d,sm8c,sige,epi)
        do j=1,3
        epi(j)=epi(j)+depi(j)
        do k=1,3
          sige(j)=sige(j)+d(j,k)*depi(k)
        enddo
        enddo
      enddo

10    return
      end


!****************************************************************************************
!
!							           TOLL subroutine
!
!	This subroutine calculates square of the error.
!	This subroutine is called in: CONTROL
!	This subroutine calls: -
!
!	variables used are:
!
!		dt		: R vector in this time step "N"
!		di		: it shows the difference between two time steps "DT-DTI"
!
!	INPUT	:
!	OUTPUT	:
!
!
!****************************************************************************************
!
	subroutine toll(di,dt,ie8d,ie8c,ie8u,id,rnn,rnn1,tol1,tol2,tol3,tol4,tol5)
!
	implicit double precision (a-h,o-z)
	common /a/ nnp,ifem,ibem,idyn,ieaq,nr
	common /a1/ ne8d,ne8c,ne8u
	common /a3/ ne8d1,ne8c1,ne8u1,ne3u1,nload1,nbcx1,nbcy1,nbcw1,nbca1,m8d1,m8c1,m8u1,&
				nbed11,nbed21,nbed31,nbeint1,ne8dd,ne8cc,ne8uu
	common /h/ mdof,mdofn,ls,ls1,isolv
!
	dimension di(mdofn),dt(mdofn),ie8d(9,ne8d1),ie8c(9,ne8c1),ie8u(9,ne8u1),id(4,nnp),&
			  rnn(mdof),rnn1(mdof)
!
!
!	norm 1a en d�placements g�n�ralis�s
!	tol1 = xnorm(di,mdofn)/xnorm(dt,mdofn)
!	tol2 = 1.d0
!
!	norm 1b en d�placements et en pression
!
	ts1=0.d0
	xms1=0.d0
	ts2=0.d0
	xms2=0.d0
	ts3=0.d0
	xms3=0.d0
	if (ne8d.eq.0) goto 10
	do ield=1,ne8d
		do i=1,8
			j=ie8d(i,ield)
			ts1=ts1+di(id(1,j))**2+di(id(2,j))**2
			xms1=xms1+dt(id(1,j))**2+dt(id(2,j))**2
		enddo
	enddo
!
10	if (ne8c.eq.0) goto 15
	do ield=1,ne8c
		do i=1,8
			j=ie8c(i,ield)
			ts1=ts1+di(id(1,j))**2+di(id(2,j))**2
			xms1=xms1+dt(id(1,j))**2+dt(id(2,j))**2
			ts2=ts2+di(id(3,j))**2
			xms2=xms2+dt(id(3,j))**2
		enddo
	enddo
!
15	if (ne8u.eq.0) goto 20
	do ield=1,ne8u
		do i=1,8
			j=ie8u(i,ield)
			ts1=ts1+di(id(1,j))**2+di(id(2,j))**2
			xms1=xms1+dt(id(1,j))**2+dt(id(2,j))**2
			ts2=ts2+di(id(3,j))**2
			xms2=xms2+dt(id(3,j))**2
			ts3=ts3+di(id(4,j))**2
			xms3=xms3+dt(id(4,j))**2
		enddo
	enddo
!
20	tol1=dsqrt(ts1)/dsqrt(xms1)
	tol2=0.d0
	if (xms2.ne.0.d0) tol2=dsqrt(ts2)/dsqrt(xms2)
	tol3=0.d0
	if (xms3.ne.0.d0) tol3=dsqrt(ts2)/dsqrt(xms2)
!
!	norm 2 en r�sidu
!
	ts1=0.d0
	xms1=0.d0
	do i=1,mdof
		ts1=ts1+rnn1(i)**2
        xms1=xms1+rnn(i)**2
	enddo
	tol4=0.d0
	if (xms1.ne.0.d0) tol4=dsqrt(ts1)/dsqrt(xms1)
!
!	norm 3 en �nergie
!
	tol5=1.d0
!
	return
	end


!****************************************************************************************
!								COMPARISON OF TWO character STRINGS
!****************************************************************************************
!
	logical function strcomp(a,b)
	character*5 a,b
	integer i,inc
!
	data inc/32/
!
	strcomp = .false.
	do i = 1,5
		ia= ichar (a(i:i))
		ib= ichar (b(i:i))
		if (ia.ne.ib .and. (ia+inc).ne.ib .and. ia.ne.(ib+inc)) return
	enddo
	strcomp = .true.
	return
	end
!
!
!****************************************************************************************
!							        NUMTOSTRING subroutine
!	This subroutine converse a integer number to a string.
!	This subroutine is called in these subroutines: GHsaveD
!****************************************************************************************
!
	subroutine NumToStr(a,i)
!
	integer i,j,k,n
	character*(*) a
!
	j=i
	n=0
!
!	n is number of digit de i
10	if (j.ne.0) then
		n=n+1
		j=j/10
		goto 10
	endif
	j=i
	a=''
!
20	if (j.ne.0) then
		k=j-(j/10)*10
		j=j/10
		if (k.eq.0) a(n:n) =  '0'
		if (k.eq.1) a(n:n) =  '1'
		if (k.eq.2) a(n:n) =  '2'
		if (k.eq.3) a(n:n) =  '3'
		if (k.eq.4) a(n:n) =  '4'
		if (k.eq.5) a(n:n) =  '5'
		if (k.eq.6) a(n:n) =  '6'
		if (k.eq.7) a(n:n) =  '7'
		if (k.eq.8) a(n:n) =  '8'
		if (k.eq.9) a(n:n) =  '9'
		n=n-1
		goto 20
	endif
	return
	end
!
!****************************************************************************************
!								ARC (inverse) HYPERBOLIC COSINE
!	The result is based on the relation arccosh x = ln(x+(x**2-1)**0.5)
!	where the square root is taken to be positive.
!	This form is used directly for 1 < x < 10^k, where k = n/2 + 1, and the
!	machine uses approximately n decimal place arithmetic.
!	For x>10^k, (x**2-1)**0.5 is equal to x**0.5 to within the accuracy of
!	the machine and hence we can guard against premature overflow and,
!	without loss of accuracy, calculate arccosh x = ln2+lnx.
!****************************************************************************************
!
	double precision function dacosh(x)
!
	double precision x, dln2, xmax
	save dln2, xmax
	data dln2 / 0.69314718055994530941723212145818d0 /
	data xmax / 0.d0 /
!
	if (xmax.eq.0.d0) xmax = 2.0d0/2.2204404925031d-18
!
	if (x.lt.1.d0) goto 10
	if (x.lt.xmax) dacosh = dlog (x+dsqrt (x*x-1.0d0))
	if (x.ge.xmax) dacosh = dln2 + dlog (x)
	return
!
10	write (*,*) 'DACOSH:  x less than 1'
	return
	end
!
!****************************************************************************************
!								INTEGRATION NUMERIQUE
!            f      - double precision
!                     function subprogam defining the integrand
!            a      - double precision
!                     lower limit of integration (in x)
!            b      - double precision (in x)
!                     upper limit of integration
!            y      - second variable of f(x,y)
!            epsabs - double precision
!                     absolute accoracy requested
!            epsrel - double precision
!                     relative accuracy requested
!            key    - integer
!                     key for choice of local integration rule
!                     a gauss-kronrod pair is used with
!                       7 - 15 points if key.lt.2,
!                      10 - 21 points if key = 2,
!                      15 - 31 points if key = 3,
!                      20 - 41 points if key = 4,
!                      25 - 51 points if key = 5,
!                      30 - 61 points if key.gt.5.
!            res -    double precision
!                     approximation to the integral
!            abserr - double precision
!                     estimate of the modulus of the absolute error,
!                     which should equal or exceed abs(i-res)
!            neval  - integer
!                     number of integrand evaluations
!            ier    - integer
!                     ier = 0 normal and reliable termination of the
!                             routine. it is assumed that the requested
!                             accuracy has been achieved.
!                     ier.gt.0 abnormal termination of the routine
!                             the estimates for result and error are
!                             less reliable. it is assumed that the
!                             requested accuracy has not been achieved.
!                      error messages
!                     ier = 1 maximum number of subdivisions allowed
!                             has been achieved. one can allow more
!                             subdivisions by increasing the value of
!                             limit (and taking the according dimension
!                             adjustments into account). however, if
!                             this yield no improvement it is advised
!                             to analyze the integrand in order to
!                             determine the integration difficulaties.
!                             if the position of a local difficulty can
!                             be determined (i.e.singularity,
!                             discontinuity within the interval) one
!                             will probably gain from splitting up the
!                             interval at this point and calling the
!                             integrator on the subranges. if possible,
!                             an appropriate special-purpose integrator
!                             should be used which is designed for
!                             handling the type of difficulty involved.
!                         = 2 the occurrence of roundoff error is
!                             detected, which prevents the requested
!                             tolerance from being achieved.
!                         = 3 extremely bad integrand behaviour occurs
!                             at some points of the integration
!                             interval.
!                         = 6 the input is invalid, because
!                             (epsabs.le.0 and
!                              epsrel.lt.dmax1(50*rel.mach.acc.,0.5d-28))
!                             or limit.lt.1 or lenw.lt.limit*4.
!                             res, abserr, neval, last are set
!                             to zero.
!                             except when lenw is invalid, iwork(1),
!                             work(limit*2+1) and work(limit*3+1) are
!                             set to zero, work(1) is set to a and
!                             work(limit+1) to b.
!
      subroutine dqpsrt(limit,last,maxerr,ermax,elist,iord,nrmax)
!
!     this routine maintains the descending ordering in the
!     list of the local error estimated resulting from the
!     interval subdivision process. at each call two error
!     estimates are inserted using the sequential search
!     method, top-down for the largest error estimate and
!     bottom-up for the smallest error estimate.
!
      double precision elist,ermax,errmax,errmin
      integer i,ibeg,ido,iord,isucc,j,jbnd,jupbn,k,last,limit,maxerr,nrmax
      dimension elist(last),iord(last)
!
!     check whether the list contains more than
!     two error estimates.
!
      if(last.gt.2) go to 10
      iord(1) = 1
      iord(2) = 2
      go to 90
!
!     this part of the routine is only executed if, due to a
!     difficult integrand, subdivision increased the error
!     estimate. in the normal case the insert procedure should
!     start after the nrmax-th largest error estimate.
!
   10 errmax = elist(maxerr)
      if(nrmax.eq.1) go to 30
      ido = nrmax-1
      do 20 i = 1,ido
        isucc = iord(nrmax-1)
!     jump out of do-loop
        if(errmax.le.elist(isucc)) go to 30
        iord(nrmax) = isucc
        nrmax = nrmax-1
   20 continue
!
!     compute the number of elements in the list to be maintained
!     in descending order. this number depends on the number of
!     subdivisions still allowed.
!
   30 jupbn = last
      if(last.gt.(limit/2+2)) jupbn = limit+3-last
      errmin = elist(last)
!
!     insert errmax by traversing the list top-down,
!     starting comparison from the element elist(iord(nrmax+1)).
!
      jbnd = jupbn-1
      ibeg = nrmax+1
      if(ibeg.gt.jbnd) go to 50
      do 40 i=ibeg,jbnd
        isucc = iord(i)
!     jump out of do-loop
        if(errmax.ge.elist(isucc)) go to 60
        iord(i-1) = isucc
   40 continue
   50 iord(jbnd) = maxerr
      iord(jupbn) = last
      go to 90
!
!     insert errmin by traversing the list bottom-up.
!
   60 iord(i-1) = maxerr
      k = jbnd
      do 70 j=i,jbnd
        isucc = iord(k)
!     jump out of do-loop
        if(errmin.lt.elist(isucc)) go to 80
        iord(k+1) = isucc
        k = k-1
   70 continue
      iord(i) = last
      go to 90
   80 iord(k+1) = last
!
!     set maxerr and ermax.
!
   90 maxerr = iord(nrmax)
      ermax = elist(maxerr)
      return
      end
!
!****************************************************************************************
!							            GHsaveD subroutine
!	This subroutine prints the Gd1&Hd1 matrices. it reads & writes the g,h for increments
!	This subroutine is called in these subroutines: GHMATD1
!	This subroutine calles this subroutine:
!	NumToStr	:
!****************************************************************************************
!
	subroutine ghsaved(itm,n1,gd1,hd1,text)
!
!
	implicit double precision (a-h,o-z)
	character text*4, nts*5
!
	dimension hd1(n12gh,n12gh),gd1(n12gh,n12gh)
!
    common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /c/ nload,ntime,itmax,itime,iter,iload,time,dtime,dtimei,dtimem,xinc
!
	numb=itime-itm+1
	call NumToStr(nts,numb)
	open (2,file=text//nts)
	do i=1,n1
		do j=1,n1
			if (itm.eq.1) then
				write (2,*) gd1(i,j),hd1(i,j)
			else
				read (2,*) gd1(i,j),hd1(i,j)
			endif
		enddo
	enddo
	close (2)
!
	return
	end
!
!****************************************************************************************
!								INTEGRATION NUMERIQUE
!            f      - double precision
!                     function subprogam defining the integrand
!            a      - double precision
!                     lower limit of integration (in x)
!            b      - double precision (in x)
!                     upper limit of integration
!            y      - second variable of f(x,y)
!            epsabs - double precision
!                     absolute accoracy requested
!            epsrel - double precision
!                     relative accuracy requested
!            key    - integer
!                     key for choice of local integration rule
!                     a gauss-kronrod pair is used with
!                       7 - 15 points if key.lt.2,
!                      10 - 21 points if key = 2,
!                      15 - 31 points if key = 3,
!                      20 - 41 points if key = 4,
!                      25 - 51 points if key = 5,
!                      30 - 61 points if key.gt.5.
!            res -    double precision
!                     approximation to the integral
!            abserr - double precision
!                     estimate of the modulus of the absolute error,
!                     which should equal or exceed abs(i-res)
!            neval  - integer
!                     number of integrand evaluations
!            ier    - integer
!                     ier = 0 normal and reliable termination of the
!                             routine. it is assumed that the requested
!                             accuracy has been achieved.
!                     ier.gt.0 abnormal termination of the routine
!                             the estimates for result and error are
!                             less reliable. it is assumed that the
!                             requested accuracy has not been achieved.
!                      error messages
!                     ier = 1 maximum number of subdivisions allowed
!                             has been achieved. one can allow more
!                             subdivisions by increasing the value of
!                             limit (and taking the according dimension
!                             adjustments into account). however, if
!                             this yield no improvement it is advised
!                             to analyze the integrand in order to
!                             determine the integration difficulaties.
!                             if the position of a local difficulty can
!                             be determined (i.e.singularity,
!                             discontinuity within the interval) one
!                             will probably gain from splitting up the
!                             interval at this point and calling the
!                             integrator on the subranges. if possible,
!                             an appropriate special-purpose integrator
!                             should be used which is designed for
!                             handling the type of difficulty involved.
!                         = 2 the occurrence of roundoff error is
!                             detected, which prevents the requested
!                             tolerance from being achieved.
!                         = 3 extremely bad integrand behaviour occurs
!                             at some points of the integration
!                             interval.
!                         = 6 the input is invalid, because
!                             (epsabs.le.0 and
!                              epsrel.lt.dmax1(50*rel.mach.acc.,0.5d-28))
!                             or limit.lt.1 or lenw.lt.limit*4.
!                             res, abserr, neval, last are set
!                             to zero.
!                             except when lenw is invalid, iwork(1),
!                             work(limit*2+1) and work(limit*3+1) are
!                             set to zero, work(1) is set to a and
!                             work(limit+1) to b.
!
!****************************************************************************************
!
      subroutine dqag(f,a,b,y,res)
!
      double precision a,abserr,b,epsabs,epsrel,f,res,work,y
      integer ier,iwork,key,last,lenw,limit,l1,l2,l3,neval
      parameter (limit=10, epsabs = 0.5d-16,epsrel=0.5d-4,lenw=41,key= 1)
!     limit > 0  lenw > 4*limit   key = 1 (2,3,4,5,6)
!     epsabs > 0 epsrel > 1.d-32
      dimension iwork(limit),work(lenw)
      external f
!
      ier = 6
      neval = 0
      last = 0
      res = 0.0d+00
      abserr = 0.0d+00
!
      l1 = limit+1
      l2 = limit+l1
      l3 = limit+l2
!
      call dqage(f,a,b,y,epsabs,epsrel,key,limit,res,abserr,neval,ier,work(1),work(l1),work(l2),work(l3),iwork,last)
!
      return
      end
!
!****************************************************************************************

      subroutine dqage(f,a,b,y,epsabs,epsrel,key,limit,res,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
!
	  double precision a,abserr,alist,area,area1,area12,area2,a1,a2,b,blist,b1,b2,d1mach,dabs,defabs,defab1, &
	  				   defab2,dmax1,elist,epmach,epsabs,epsrel,errbnd,errmax,error1,error2,erro12,errsum,f,resabs,res,rlist,y,uflow
      integer ier,iord,iroff1,iroff2,k,key,keyf,last,limit,maxerr,neval,nrmax
!
      dimension alist(limit),blist(limit),elist(limit),iord(limit),rlist(limit)
      common /machine/ epmach,uflow
      external f
!
      epmach=d1mach(4)
      ufow=d1mach(1)
!      epmach = 2.2204404925032D-18
!      uflow = 2.0D0*2.225073858507201D-308
!
!     valeurs initiales des parametres
      ier = 0
      neval = 0
      last = 0
      res = 0.0d+00
      abserr = 0.0d+00
      alist(1) = a
      blist(1) = b
      rlist(1) = 0.0d+00
      elist(1) = 0.0d+00
      iord(1) = 0
      if(epsabs.le.0.0d+00.and.epsrel.lt.dmax1(0.5d+02*epmach,0.5d-28))ier = 6
      if(ier.eq.6) go to 999
!
!     first approximation to the integral
!     -----------------------------------
      keyf = key
      if(key.le.0) keyf = 1
      if(key.ge.7) keyf = 6
      neval = 0
      if(keyf.eq.1) call dqk15(f,a,b,y,res,abserr,defabs,resabs)
!      if(keyf.eq.2) call dqk21(f,a,b,y,res,abserr,defabs,resabs)
!      if(keyf.eq.3) call dqk31(f,a,b,y,res,abserr,defabs,resabs)
!      if(keyf.eq.4) call dqk41(f,a,b,y,res,abserr,defabs,resabs)
!      if(keyf.eq.5) call dqk51(f,a,b,y,res,abserr,defabs,resabs)
!      if(keyf.eq.6) call dqk61(f,a,b,y,res,abserr,defabs,resabs)
      last = 1
      rlist(1) = res
      elist(1) = abserr
      iord(1) = 1
!
!     test on accuracy
      errbnd = dmax1(epsabs,epsrel*dabs(res))
      if(abserr.le.0.5d+02*epmach*defabs.and.abserr.gt.errbnd) ier = 2
      if(limit.eq.1) ier = 1
      if(ier.ne.0.or.(abserr.le.errbnd.and.abserr.ne.resabs).or.abserr.eq.0.0d+00) go to 60
!
!     initialization
      errmax = abserr
      maxerr = 1
      area = res
      errsum = abserr
      nrmax = 1
      iroff1 = 0
      iroff2 = 0
!
!     main do-loop
!     ------------
      do 30 last = 2,limit
!
!     bisect the subinterval with the largest error estimate.
        a1 = alist(maxerr)
        b1 = 0.5d+00*(alist(maxerr)+blist(maxerr))
        a2 = b1
        b2 = blist(maxerr)
        if(keyf.eq.1) call dqk15(f,a1,b1,y,area1,error1,resabs,defab1)
!        if(keyf.eq.2) call dqk21(f,a1,b1,y,area1,error1,resabs,defab1)
!        if(keyf.eq.3) call dqk31(f,a1,b1,y,area1,error1,resabs,defab1)
!        if(keyf.eq.4) call dqk41(f,a1,b1,y,area1,error1,resabs,defab1)
!        if(keyf.eq.5) call dqk51(f,a1,b1,y,area1,error1,resabs,defab1)
!        if(keyf.eq.6) call dqk61(f,a1,b1,y,area1,error1,resabs,defab1)
        if(keyf.eq.1) call dqk15(f,a2,b2,y,area2,error2,resabs,defab2)
!        if(keyf.eq.2) call dqk21(f,a2,b2,y,area2,error2,resabs,defab2)
!        if(keyf.eq.3) call dqk31(f,a2,b2,y,area2,error2,resabs,defab2)
!        if(keyf.eq.4) call dqk41(f,a2,b2,y,area2,error2,resabs,defab2)
!        if(keyf.eq.5) call dqk51(f,a2,b2,y,area2,error2,resabs,defab2)
!        if(keyf.eq.6) call dqk61(f,a2,b2,y,area2,error2,resabs,defab2)
!
!     improve previous approximations to integral
!     and error and test for accuracy.
!
        neval = neval+1
        area12 = area1+area2
        erro12 = error1+error2
        errsum = errsum+erro12-errmax
        area = area+area12-rlist(maxerr)
        if(defab1.eq.error1.or.defab2.eq.error2) go to 5
        if(dabs(rlist(maxerr)-area12).le.0.1d-04*dabs(area12).and.erro12.ge.0.99d+00*errmax) iroff1 =iroff1+1
        if(last.gt.10.and.erro12.gt.errmax) iroff2 = iroff2+1
5       rlist(maxerr) = area1
        rlist(last) = area2
        errbnd = max(epsabs,epsrel*dabs(area))
        if(errsum.le.errbnd) go to 8
!
!     test for roundoff error and eventually set error flag.
        if(iroff1.ge.6.or.iroff2.ge.20) ier = 2
!
!     set error flag in the case that the number of subintervals equals limit.
        if(last.eq.limit) ier = 1
!
!     set error flag in the case of bad integrand behaviour
!     at a point of the integration range.
        if(dmax1(dabs(a1),dabs(b2)).le.(0.1d+01+0.1d+03*epmach)*(dabs(a2)+0.1d+04*uflow)) ier = 3
!
!     append the newly-created intervals to the list.
!
8       if(error2.gt.error1) go to 10
        alist(last) = a2
        blist(maxerr) = b1
        blist(last) = b2
        elist(maxerr) = error1
        elist(last) = error2
        go to 20
   10   alist(maxerr) = a2
        alist(last) = a1
        blist(last) = b1
        rlist(maxerr) = area2
        rlist(last) = area1
        elist(maxerr) = error2
        elist(last) = error1
!
!     call subroutine dqpsrt to maintain the descending ordering
!     in the list of error estimates and select the subinterval
!     with the largest error estimate (to be bisected next).
!
   20   call dqpsrt(limit,last,maxerr,errmax,elist,iord,nrmax)
!     jump out of do-loop
      if(ier.ne.0.or.errsum.le.errbnd) go to 40
   30 continue
!
!     compute final result.
!
   40 res = 0.0d+00
      do 50 k=1,last
        res = res+rlist(k)
   50 continue
      abserr = errsum
   60 if(keyf.ne.1) neval = (10*keyf+1)*(2*neval+1)
      if(keyf.eq.1) neval = 30*neval+15
!
  999 return
      end
!
!c----------------------------------------------------------------------
      subroutine dqk15(f,a,b,y,res,abserr,resabs,resasc)
!c
	  double precision a,absc,abserr,b,centr,dabs,dhlgth,dmax1,dmin1,f,fc,fsum,fval1,fval2,fv1,fv2, & 
	  				   hlgth,resabs,resasc,resg,resk,reskh,res,wg,wgk,xgk,y,epmach,uflow
      integer j,jtw,jtwm1
      common /machine/ epmach,uflow
      external f
      dimension fv1(7),fv2(7),wg(4),wgk(8),xgk(8)
      save wg,xgk,wgk
!c
      data wg  (  1) /0.129484966168869693270611432679082d0/
      data wg  (  2) /0.279705391489276667901467771423780d0/
      data wg  (  3) /0.381830050505118944950369775488975d0/
      data wg  (  4) /0.417959183673469387755102040816327d0/
!c
      data xgk (  1) /0.991455371120812639206854697526329d0/
      data xgk (  2) /0.949107912342758524526189684047851d0/
      data xgk (  3) /0.864864423359769072789712788640926d0/
      data xgk (  4) /0.741531185599394439863864773280788d0/
      data xgk (  5) /0.586087235467691130294144838258730d0/
      data xgk (  6) /0.405845151377397166906606412076961d0/
      data xgk (  7) /0.207784955007898467600689403773245d0/
      data xgk (  8) /0.000000000000000000000000000000000d0/
!c
      data wgk (  1) /0.022935322010529224963732008058970d0/
      data wgk (  2) /0.063092092629978553290700663189204d0/
      data wgk (  3) /0.104790010322250183839876322541518d0/
      data wgk (  4) /0.140653259715525918745189590510238d0/
      data wgk (  5) /0.169004726639267902826583426598550d0/
      data wgk (  6) /0.190350578064785409913256402421014d0/
      data wgk (  7) /0.204432940075298892414161999234649d0/
      data wgk (  8) /0.209482141084727828012999174891714d0/
!c
      centr = 0.5d+00*(a+b)
      hlgth = 0.5d+00*(b-a)
      dhlgth = dabs(hlgth)
!c
!c     compute the 15-point kronrod approximation to
!c     the integral, and estimate the absolute error.
!c
      fc = f(centr,y)
      resg = fc*wg(4)
      resk = fc*wgk(8)
      resabs = dabs(resk)
      do 10 j=1,3
        jtw = j*2
        absc = hlgth*xgk(jtw)
        fval1 = f(centr-absc,y)
        fval2 = f(centr+absc,y)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(dabs(fval1)+dabs(fval2))
   10 continue
      do 15 j = 1,4
        jtwm1 = j*2-1
        absc = hlgth*xgk(jtwm1)
        fval1 = f(centr-absc,y)
        fval2 = f(centr+absc,y)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(dabs(fval1)+dabs(fval2))
   15 continue
      reskh = resk*0.5d+00
      resasc = wgk(8)*dabs(fc-reskh)
      do 20 j=1,7
        resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))
   20 continue
      res = resk*hlgth
      resabs = resabs*dhlgth
      resasc = resasc*dhlgth
      abserr = dabs((resk-resg)*hlgth)
      if(resasc.ne.0.0d+00.and.abserr.ne.0.0d+00) abserr = resasc*dmin1(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00)
      if(resabs.gt.uflow/(0.5d+02*epmach)) abserr = dmax1 ((epmach*0.5d+02)*resabs,abserr)
      return
      end
!****************************************************************************************
!
!							           UPDATE1 subroutine
!
!	This subroutine updates displacement arrays after each iteration.
!	This subroutine is called in: CONTROL
!	This subroutine calls: -
!
!	variables used are:
!
!		r		: all nodal forces (mechanical TOTRR & seismic RRR) are saved in this
!				  vector; after that in SOLVE subroutine the displacements & pressures
!				  obtained by solve the equations are saved in this vector in each itime,
!				  this vector is saved in DT vector in UPDATE1 subroutine & initialized in it
!		dt		: R vector in this time step "N"
!		dti		: R vector in previous time step "N-1"
!		di		: it shows the difference between two time steps "DT-DTI"
!
!
!	INPUT	:
!
!
!	OUTPUT	:
!
!
!****************************************************************************************
!
	subroutine update1(r,di,dt,dti)
!
 	implicit double precision (a-h,o-z)
!
	common /g/ init,iprint
	common /h/ mdof,mdofn,ls,ls1,isolv
	common /a/ nnp,ifem,ibem,idyn,ieaq,nr
	common /chrg3/ nstime,neaq,nneaq(500),uxeaq,uyeaq,axeaq,ayeaq,uxeq(3,500),uyeq(3,500)
!
	dimension r(mdof),di(mdofn),dt(mdofn),dti(mdofn)
!
!
!	---------- init=1,2 ----------
!
	if (init.eq.0) goto 100
	do i=1,mdof
		dt(i)=r(i)
		r(i)=0.0d0
	enddo
	goto 200
!
!
!	----------- init=0 -----------
!
100	do i=1,mdof
		dti(i)=dt(i)
        dt(i)=r(i)
		r(i)=0.0d0
        di(i)=dt(i)-dti(i)
	enddo
!
200	return
	end

!****************************************************************************************
!
!							           UPDATE3 subroutine
!
!	This subroutine updates displacement, velocity, acceleration and stress arrays after
!	each increment.
!	This subroutine is called in: CONTROL
!	This subroutine calls: -
!
!	variables used are:
!
!		disp	: displacement/pressures vector which consists of initial displacement value
!				  of corresponding nodes in which their identities code of Ux & Uy are defined
!				  in ID and initial water and air values of corresponding nodes in which their
!				  identities code of Pw & Pa are defined in ID. these displacement initial
!				  values can be defined for a free node or for a node with boundary conditions
!				  and the water and air initial values must be defined for a free node;
!				  in each "itime" the displacements/pressures vector (R->DT) in PREVIOUS
!				  TIME STEP (N-1) is transfered in DISP vector
!		vel		: velocity vector which consists of initial velocity value of corresponding
!				  nodes in which their identities code of Ux & Uy are defined in ID. these
!				  initial values can be defined for a free node or for a node with b.c.;
!				  in FEM the velocity vector is obtained by:
!				  V(N)=V(N-1)+[(1-anew1)*A(N-1)+anew1*A(N)]*dtime
!		acc		: acceleration vector which consists of initial acceleration value of
!				  corresponding nodes in which their identities code of Ux & Uy are defined
!				  in ID. these initial values can be defined for a free node or for a node
!				  with boundary conditions;
!				  in FEM the velocity vector is obtained by:
!				  A(N)=(1/anew2/dtime**2)*[U(N)-U(N-1)]-(1/anew2/dtime)*V(N-1)
!						 -(1/2/anew2-1)*A(N-1)
!		temp	: acceleration vector in previous time step A(N-1)
!		sig8d	:
!		sig8c	:
!		sig8u	:
!		sigmad	:
!		sigmac	:
!		sigmau	:
!		dt		: R vector in this time step "N"
!		r		: all nodal forces (mechanical TOTRR & seismic RRR) are saved in this
!				  vector; after that in SOLVE subroutine the displacements & pressures
!				  obtained by solve the equations are saved in this vector in each itime,
!				  this vector is saved in DT vector in UPDATE1 subroutine & initialized in it
!		dti		: R vector in previous time step "N-1"
!		di		: it shows the difference between two time steps "DT-DTI"
!
!
!	INPUT	:
!
!
!	OUTPUT	:
!
!
!****************************************************************************************
!
	subroutine update3(disp,vel,acc,sig8d,sig8c,sig8u,sigmad,sigmac,sigmau,dt)
!
 	implicit double precision (a-h,o-z)
!
	common /a/ nnp,ifem,ibem,idyn,ieaq,nr
	common /a1/ ne8d,ne8c,ne8u
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /a3/ ne8d1,ne8c1,ne8u1,ne3u1,nload1,nbcx1,nbcy1,nbcw1,nbca1,m8d1,m8c1,m8u1,&
				nbed11,nbed21,nbed31,nbeint1,ne8dd,ne8cc,ne8uu
	common /c/ nload,ntime,itmax,itime,iter,iload,time,dtime,dtimei,dtimem,xinc
	common /h/ mdof,mdofn,ls,ls1,isolv
	common /stability/anew1,anew2,wil,intb
!
	dimension disp(mdofn),vel(mdofn),acc(mdofn),sig8d(9,5,ne8d1),sig8c(9,6,ne8c1),&
			  sig8u(9,7,ne8u1),sigmad(9,5,ne8d1),sigmac(9,6,ne8c1),sigmau(9,7,ne8u1),&
			  dt(mdofn),anew(8),temp(mdofn) !,&
			  !sig8di(9,5,ne8d1),sig8ci(9,6,ne8c1),sig8ui(9,6,ne8u1)
!
!
!	----------------------- Newmark's coefficients -----------------------
!
	anew(1)=1.d0/(anew2*dtime**2)
!	anew(2)=anew1/(anew2*dtime)
	anew(3)=1.d0/(anew2*dtime)
	anew(4)=1.d0/(2*anew2)-1.d0
!	anew(5)=anew1/anew2-1.d0
!	anew(6)=dtime/2.d0*(anew1/anew2-2.d0)
	anew(7)=dtime*(1.d0-anew1)
	anew(8)=anew1*dtime
!
!
!	------------------------ update accelerations ------------------------
!
	if (idyn.ne.2) goto 10
	do i=1,mdofn
		temp(i)=acc(i)
        acc(i)=(dt(i)-disp(i))*anew(1)-vel(i)*anew(3)-acc(i)*anew(4)
	enddo
!
!
!	------------------------- update velocities --------------------------
!
10	if (idyn.eq.0) goto 20
	do i=1,mdofn
!		vel(i)=(dt(i)-disp(i))*(2/dtime)-vel(i)
        vel(i)=vel(i)+temp(i)*anew(7)+acc(i)*anew(8)
      enddo
!
!
!	------------------------ update displacements ------------------------
!
20	do i=1,mdofn
		disp(i)=dt(i)
	enddo
!
!
!	------------------- update drained element stresses ------------------
!
	if (ne8d.eq.0) goto 30
	do i=1,ne8d
		do k=1,9
			do j=1,5
!				sig8di(k,j,i)=sig8d(k,j,i)
				sigmad(k,j,i)=sig8d(k,j,i)
			enddo
		enddo
	enddo
!
!
!	---------------- update consolidation element stresses ---------------
!
30	if (ne8c.eq.0) goto 35
	do i=1,ne8c
		do k=1,9
			do j=1,6
!				sig8ci(k,j,i)=sig8c(k,j,i)
				sigmac(k,j,i)=sig8c(k,j,i)
			enddo
		enddo
	enddo
!
!
!	----------------- update unsaturated element stresses ----------------
!
35	if (ne8u.eq.0) goto 40
	do i=1,ne8u
		do k=1,9
			do j=1,7
!				sig8ui(k,j,i)=sig8u(k,j,i)
				sigmau(k,j,i)=sig8u(k,j,i)
			enddo
		enddo
	enddo
!
!
40	return
	end



	subroutine wdpsoldyG(SLAP,DRT,RT,K0,K1,DK0,DK1,CsS,CsW,CsA,CwS,CaS,CCwW,CwA,&
						  CaW,CaA,dpsolG)
!
	implicit double precision (a-h,o-z)
!
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /point/ dist(2),ra,rd(2),rdn,eta(2)
	common /umate/ d(3,3),dfs(3),xn,xkw,cww,cwg,xka,cgg,cgw,ev,sat,g1,rs,rw,raa,rmix,xm1,&
				   bt,et,xlambda,xmu,henry
	common /ocqm/ EPScqm,taw,Icqm,Lcqm,NPcqm
	common /UcoefR/ r1(2,2),dpr1(2,2,2),r2(2,2),dpr2(2,2,2),c(2,2)
!
	COMPLEX(8) SLAP,RT(4),DRT(4),K0(4),K1(4),DK0(4),DK1(4),dpsolG(kbem,kbem,2),temp,&
			   CsS(4),CsW(3),CsA(3),CwS(3),CaS(3),CCwW(3),CwA(3),CaW(3),CaA(3)
	dimension dlt(2,2)
!
!
	do i=1,kbem
		do j=1,kbem
			do m=1,2
				dpsolG(i,j,m)=0.d0
			enddo
		enddo
	enddo
!
	pi=DACOS (-1.d0)
	zmu=3.d0*bt*et/(9.d0*bt-et)
	zlamda=bt-2.d0*zmu/3.d0
	c1=1./(2.*pi*zmu)
	sata=1.-sat
!
	dlt(1,1)=1.d0
	dlt(2,2)=1.d0
	dlt(1,2)=0.d0
	dlt(2,1)=0.d0
!
	do i=1,4
		do j=1,4
			if (i.eq.3.and.j.lt.3) GOTO 10
			if (i.eq.3.and.j.eq.3) GOTO 20
			if (i.eq.3.and.j.eq.4) GOTO 30
			if (i.eq.4.and.j.lt.3) GOTO 40
			if (i.eq.4.and.j.eq.3) GOTO 50
			if (i.eq.4.and.j.eq.4) GOTO 60
			if (i.lt.3.and.j.eq.3) GOTO 70
			if (i.lt.3.and.j.eq.4) GOTO 80
!
!	--- dpGij ---
!
			do m=1,2
				temp=DCMPLX (0.d0,0.d0)
				do k=1,4
					temp=temp+CsS(k)*( DRT(k)*(dpr1(i,j,m)*K1(k)+rd(m)*r1(i,j)*DK1(k))+&
									   RT(k)* (dpr2(i,j,m)*K0(k)+rd(m)*r2(i,j)*DK0(k)) )
				enddo
				temp=temp+c(i,j)*rd(m)*DK0(1)
				dpsolG(i,j,m)=temp
			enddo
			goto 100

!	--- dpG3j ---
!
10			do m=1,2
				temp=DCMPLX (0.d0,0.d0)
				do k=2,4
					temp=temp+CwS(k-1)*( (dlt(j,m)-rd(j)*rd(m))*K1(k)/ra+&
										 rd(j)*rd(m)*DK1(k) )
				enddo
				dpsolG(3,j,m)=temp
			enddo
			goto 100
!
!	--- dpG33 ---
!
20			do m=1,2
				temp=DCMPLX (0.d0,0.d0)
				do k=2,4
					temp=temp+CCwW(k-1)*rd(m)*DK0(k)
				enddo
				dpsolG(3,3,m)=temp
			enddo
			goto 100
!
!	--- dpG34 ---
!
30			do m=1,2
				temp=DCMPLX (0.d0,0.d0)
				do k=2,4
					temp=temp+CwA(k-1)*rd(m)*DK0(k)
				enddo
				dpsolG(3,4,m)=temp
			enddo
			goto 100
!
!	--- dpG4j ---
!
40			do m=1,2
				temp=DCMPLX (0.d0,0.d0)
				do k=2,4
					temp=temp+CaS(k-1)*( (dlt(j,m)-rd(j)*rd(m))*K1(k)/ra+&
										 rd(j)*rd(m)*DK1(k) )
				enddo
				dpsolG(4,j,m)=temp
			enddo
			goto 100
!
!	--- dpG43 ---
!
50			do m=1,2
				temp=DCMPLX (0.d0,0.d0)
				do k=2,4
					temp=temp+CaW(k-1)*rd(m)*DK0(k)
				enddo
				dpsolG(4,3,m)=temp
			enddo
			goto 100
!
!	--- dpG44 ---
!
60			do m=1,2
				temp=DCMPLX (0.d0,0.d0)
				do k=2,4
					temp=temp+CaA(k-1)*rd(m)*DK0(k)
				enddo
				dpsolG(4,4,m)=temp
			enddo
			goto 100
!
!	--- dpGi3 ---
!
70			do m=1,2
				temp=DCMPLX (0.d0,0.d0)
				do k=2,4
					temp=temp+CsW(k-1)*( (dlt(i,m)-rd(i)*rd(m))*K1(k)/ra+&
										 rd(i)*rd(m)*DK1(k) )
				enddo
				dpsolG(i,3,m)=temp
			enddo
			goto 100
!
!	--- dpGi4 ---
!
80			do m=1,2
				temp=DCMPLX (0.d0,0.d0)
				do k=2,4
					temp=temp+CsA(k-1)*( (dlt(i,m)-rd(i)*rd(m))*K1(k)/ra+&
										 rd(i)*rd(m)*DK1(k) )
				enddo
				dpsolG(i,4,m)=temp
			enddo
			goto 100
100		continue
		enddo
	enddo
!
	return
	end



!*******************************************************************************************
!
	subroutine wdpsolqsG(SLAP,DRTQ,RTQ,KQ0,KQ1,DKQ0,DKQ1,GAM11,GAM12,GAM13,GAM21,GAM22,&
						  GAM31,GAM32,GAM33,GAM1,GAM2,GAM3,dpsolG)
!
	implicit double precision (a-h,o-z)
!
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /point/ dist(2),ra,rd(2),rdn,eta(2)
	common /umate/ d(3,3),dfs(3),xn,xkw,cww,cwg,xka,cgg,cgw,ev,sat,g1,rs,rw,raa,rmix,xm1,&
				   bt,et,xlambda,xmu,henry
	common /ocqm/ EPScqm,taw,Icqm,Lcqm,NPcqm
	common /UcoefR/ r1(2,2),dpr1(2,2,2),r2(2,2),dpr2(2,2,2),c(2,2)
!
	COMPLEX(8) SLAP,RTQ(2),DRTQ(2),KQ0(2),KQ1(2),DKQ0(2),DKQ1(2),GAM11,GAM12,GAM13,GAM21,GAM22,&
			   GAM31,GAM32,GAM33,GAM1,GAM2,GAM3,dpsolG(kbem,kbem,2),temp,&
			   K11,K12,K13,K21,K22,K23,K31,K32,K41,K42,K51,K52,K61,K62,K71,K72,K73,K74,&
			   K75,K76,K77,dpGAM11(2),dpGAM12(2),dpGAM13(2),dpGAM21(2),dpGAM22(2),&
			   dpGAM31(2),dpGAM32(2),dpGAM33(2),dpGAM1(2),dpGAM2(2),dpGAM3(2)
	dimension dlt(2,2)
!
!
	do i=1,kbem
		do j=1,kbem
			do m=1,2
				dpsolG(i,j,m)=0.d0
			enddo
		enddo
	enddo
!
	call WFSOLGCOEFQST(SLAP,K11,K12,K13,K21,K22,K23,K31,K32,K41,K42,K51,K52,K61,K62,K71,K72,&
					   K73,K74,K75,K76,K77)
!
	dlt(1,1)=1.d0
	dlt(2,2)=1.d0
	dlt(1,2)=0.d0
	dlt(2,1)=0.d0
!
	do i=1,2
		DRTQ(i)=Cdsqrt(RTQ(i))
	enddo
!
!	--- dpGAM ---
!
	do m=1,2
		dpGAM11(m)=-rd(m)*( KQ1(2)*RTQ(2)**1.5d0-KQ1(1)*RTQ(1)**1.5d0 )/( (RTQ(2)-RTQ(1)) )
		dpGAM12(m)=-rd(m)*SLAP*( KQ1(2)*DRTQ(2)-KQ1(1)*DRTQ(1) )/(RTQ(2)-RTQ(1))
		dpGAM13(m)=-rd(m)*(SLAP**2.d0)*( KQ1(2)/DRTQ(2)-KQ1(1)/DRTQ(1) )/(RTQ(2)-RTQ(1))
!
		dpGAM21(m)=-rd(m)*SLAP*( (RTQ(2)*KQ0(2)-RTQ(1)*KQ0(1))+&
						    (DRTQ(2)*KQ1(2)/ra-DRTQ(1)*KQ1(1)/ra) )/(RTQ(2)-RTQ(1))
		dpGAM22(m)=-rd(m)*(SLAP**2.d0)*( (KQ0(2)-KQ0(1))+&
							(KQ1(2)/(DRTQ(2)*ra)-KQ1(2)/(DRTQ(1)*ra)) )/(RTQ(2)-RTQ(1))
!
		dpGAM31(m)=-rd(m)*( (RTQ(2)*KQ0(2)-RTQ(1)*KQ0(1))+&
						    (DRTQ(2)*KQ1(2)/ra-DRTQ(1)*KQ1(1)/ra) )/(RTQ(2)-RTQ(1))
		dpGAM32(m)=-rd(m)*SLAP*( (KQ0(2)-KQ0(1))+&
							(KQ1(2)/(DRTQ(2)*ra)-KQ1(2)/(DRTQ(1)*ra)) )/(RTQ(2)-RTQ(1))
		dpGAM33(m)=-rd(m)*(SLAP**2.d0)*( (KQ0(2)/RTQ(2)-KQ0(1)/RTQ(1))+&
							(KQ1(2)/(ra*RTQ(2)**1.5d0)-KQ1(1)/(ra*RTQ(1)**1.5d0)) )/&
							(RTQ(2)-RTQ(1))
!
		dpGAM1(m)=K11*dpGAM11(m)+K12*dpGAM12(m)+K13*dpGAM13(m)
		dpGAM2(m)=K21*dpGAM31(m)+K22*dpGAM32(m)+K23*dpGAM33(m)
		dpGAM3(m)=K21*dpGAM11(m)+K22*dpGAM12(m)+K23*dpGAM13(m)
	enddo

!
	do i=1,4
		do j=1,4
			if (i.eq.3.and.j.lt.3) GOTO 10
			if (i.eq.3.and.j.eq.3) GOTO 20
			if (i.eq.3.and.j.eq.4) GOTO 30
			if (i.eq.4.and.j.lt.3) GOTO 40
			if (i.eq.4.and.j.eq.3) GOTO 50
			if (i.eq.4.and.j.eq.4) GOTO 60
			if (i.lt.3.and.j.eq.3) GOTO 70
			if (i.lt.3.and.j.eq.4) GOTO 80
!
!	--- dpGij ---
!
			do m=1,2
				temp=DCMPLX (0.d0,0.d0)
				temp=temp+dlt(i,j)*dpGAM1(m)+dpr1(i,j,m)*GAM2+r1(i,j)*dpGAM2(m)+&
					 dpr2(i,j,m)*GAM3+r2(i,j)*dpGAM3(m)
				dpsolG(i,j,m)=temp
			enddo
			goto 100
!
!	--- dpG3j ---
!
10			do m=1,2
				temp=DCMPLX (0.d0,0.d0)
				temp=temp+( (rd(j)*rd(m)-dlt(j,m))/ra )*(K41*GAM21+K42*GAM22)-&
					         rd(j)*( K41*dpGAM21(m)+K42*dpGAM22(m) )
				dpsolG(3,j,m)=temp
			enddo
			goto 100
!
!	--- dpG33 ---
!
20			do m=1,2
				temp=DCMPLX (0.d0,0.d0)
				temp=temp+K76*dpGAM11(m)+K77*dpGAM12(m)
				dpsolG(3,3,m)=temp
			enddo
			goto 100
!
!	--- dpG34 ---
!
30			do m=1,2
				temp=DCMPLX (0.d0,0.d0)
				temp=temp+K73*dpGAM11(m)+K74*dpGAM12(m)
				dpsolG(3,4,m)=temp
			enddo
			goto 100
!
!	--- dpG4j ---
!
40			do m=1,2
				temp=DCMPLX (0.d0,0.d0)
				temp=temp+( (rd(j)*rd(m)-dlt(j,m))/ra )*(K31*GAM21+K32*GAM22)-&
					         rd(j)*( K31*dpGAM21(m)+K32*dpGAM22(m) )
				dpsolG(4,j,m)=temp
			enddo
			goto 100
!
!	--- dpG43 ---
!
50			do m=1,2
				temp=DCMPLX (0.d0,0.d0)
				temp=temp+K75*dpGAM12(m)
				dpsolG(4,3,m)=temp
			enddo
			goto 100
!
!	--- dpG44 ---
!
60			do m=1,2
				temp=DCMPLX (0.d0,0.d0)
				temp=temp+K71*dpGAM11(m)+K72*dpGAM12(m)
				dpsolG(4,4,m)=temp
			enddo
			goto 100
!
!	--- dpGi3 ---
!
70			do m=1,2
				temp=DCMPLX (0.d0,0.d0)
				temp=temp+( (rd(i)*rd(m)-dlt(i,m))/ra )*(K61*GAM21+K62*GAM22)-&
					         rd(i)*( K61*dpGAM21(m)+K62*dpGAM22(m) )
				dpsolG(i,3,m)=temp
			enddo
			goto 100
!
!	--- dpGi4 ---
!
80			do m=1,2
				temp=DCMPLX (0.d0,0.d0)
				temp=temp+( (rd(i)*rd(m)-dlt(i,m))/ra )*(K51*GAM21+K52*GAM22)-&
					         rd(i)*( K51*dpGAM21(m)+K52*dpGAM22(m) )
				dpsolG(i,4,m)=temp
			enddo
			goto 100
100		continue
		enddo
	enddo
!
	return
	end



!
	subroutine weldynfsol(W,l,SLAP,cqmG,cqmF)
!
	implicit double precision (a-h,o-z)
!
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /c/ nload,ntime,itmax,itime,iter,iload,time,dtime,dtimei,dtimem,xinc
	common /point/ dist(2),ra,rd(2),rdn,eta(2)
	common /mate/zlamda,zmuy,znuu,znu,zalpha,zrho,zrhof,zm,zk,zr,zn
	common /ocqm/ EPScqm,taw,Icqm,Lcqm,NPcqm
!
	COMPLEX(8) W,SLAP,z1,z2,cbk0,cdk0,cbk1,cdk1,K0z1,K0z2,K1z1,K1z2,K2z2,K2z1,rr1,rr2,rr3,&
			   rr4,cqmG(kbem,kbem),cqmF(kbem,kbem),A,B,dAdr,P,Q,dBdr,R,S
	dimension dlt(2,2),vel(2)
!
!
	do i=1,kbem
		do j=1,kbem
			cqmG(i,j)=0.d0
			cqmF(i,j)=0.d0
		enddo
	enddo
!
	pi=DACOS (-1.d0)
!
	vel(1)=dsqrt ((zlamda+2*zmuy)/zrho)
	vel(2)=dsqrt (zmuy/zrho)
	rap=vel(2)/vel(1)
!
	dlt(1,1)=1.d0
	dlt(2,2)=1.d0
	dlt(1,2)=0.d0
	dlt(2,1)=0.d0
!
	coef1=1.d0/(2*zmuy*pi)
	coef2=1.d0/(2*pi)
!
	z1=ra*SLAP/vel(1)
	z2=ra*SLAP/vel(2)
!
	cbk0=(0.0D0,0.0D0)
	cbk1=(0.0D0,0.0D0)
!	call MBK01 (z1,cbk0,cbk1)
	K0z1=cbk0
	K1z1=cbk1
!
	cbk0=(0.0D0,0.0D0)
	cbk1=(0.0D0,0.0D0)
!	call MBK01 (z2,cbk0,cbk1)
	K0z2=cbk0
	K1z2=cbk1
!
	K2z2=K0z2+(2.d0/z2)*K1z2
	K2z1=K0z1+(2.d0/z1)*K1z1
!
	A=K0z2+(K1z2-rap*K1z1)/z2
	B=K2z2-rap**2.d0*K2z1
!
	rr1=(K0z1*z1+2.d0*K1z1)
	rr2=(K0z2*z2+2.d0*K1z2)
	dAdr=(rap*rr1-rr2-K1z2*z2**2.d0)/(ra*z2)
!
	P=dAdr-B/ra
	Q=-2.d0*B/ra
!
	rr3=K1z1*z1+2.d0*(K0z1+2.d0*K1z1/z1)
	rr4=K1z2*z2+2.d0*(K0z2+2.d0*K1z2/z2)
	dBdr=(rap**2.d0*rr3-rr4)/ra
!
	R=-2.d0*dBdr
	S=((1.d0/rap)**2.d0-2.d0)*(dAdr-dBdr-B/ra) !***
!
	do i=1,2
		do j=1,2
			cqmG(i,j)=coef1*(A*dlt(i,j)-B*rd(i)*rd(j))
			cqmF(i,j)=coef2*( P*(rdn*dlt(i,j)+rd(j)*eta(i))+Q*(eta(j)*rd(i)-2.d0*rd(i)*rd(j)*rdn)+&
						      R*rd(i)*rd(j)*rdn+S*rd(i)*eta(j) ) !***
		enddo
	enddo
!
	return
	end





!*******************************************************************************************
!
	subroutine wfsolFcoef(SLAP,DRT,RT,CsS,CsW,CsA,CwS,CaS,CCwW,CwA,CaW,CaA)
!
	implicit double precision (a-h,o-z)
!
	common /e/ gammaw,gammas,gammaa,grav,atmp
	common /umate/ d(3,3),dfs(3),xn,xkw,cww,cwg,xka,cgg,cgw,ev,sat,g1,rs,rw,raa,rmix,xm1,&
				   bt,et,xlambda,xmu,henry
!
	COMPLEX(8) SLAP,DRT(4),RT(4),xKss1,xKss2,CsS(4),CsW(3),CsA(3),CwS(3),CaS(3),CCwW(3),&
			   CwA(3),CaW(3),CaA(3),cc3,cc4,cc5,cc6,cc7,cc8,cc9,cc10,cc11,cc12,cc13,cc14,cc15,&
			   cc16,cc17,cc18,cc19,cc20
!
!
	pi=DACOS (-1.d0)
!
	c1=xlambda+2.d0*xmu
	sata=1.D0-sat
	fs=dfs(1)
!
	cc1=1.d0/(2.d0*pi*xmu)
	cc2=-(xlambda+xmu)/c1
!
	cc3=DCMPLX (0.d0,0.d0)
	if (xkw.ne.0.d0) cc3=(sat-rw*xkw*SLAP)*SLAP/(2.d0*pi*c1*xkw)
!
	cc4=DCMPLX (0.d0,0.d0)
	if (xka.ne.0.d0) cc4=( cwg*SLAP*(sata-raa*xka*SLAP)-&
										   cgg*SLAP*(sat-rw*xkw*SLAP) )/&
										   (xka*(sat-rw*xkw*SLAP))
!
	cc5=DCMPLX (0.d0,0.d0)
	if (xka.ne.0.d0) cc5=(sata-raa*xka*SLAP)*SLAP/(2.d0*pi*c1*xka)
!
	cc6=DCMPLX (0.d0,0.d0)
	if (xkw.ne.0.d0) cc6=( cwg*SLAP*(sat-rw*xkw*SLAP)-&
										   cww*SLAP*(sata-raa*xka*SLAP) )/&
										   (xkw*(sata-raa*xka*SLAP))
!
	cc7=DCMPLX (0.d0,0.d0)
	cc7=-fs/(2.d0*pi*c1)
!
	cc8=DCMPLX (0.d0,0.d0)
	if (xka.ne.0.d0) cc8=( cwg*(1.d0-fs)-fs*cgg )*SLAP/(xka*fs)
!
	cc9=DCMPLX (0.d0,0.d0)
	cc9=-(1.d0-fs)/(2.d0*pi*c1)
!
	cc10=DCMPLX (0.d0,0.d0)
	if (xkw.ne.0.d0) cc10=( cww*(1.d0-fs)-fs*cwg )*SLAP/(-xkw*(1.d0-fs))
!
	cc11=DCMPLX (0.d0,0.d0)
	cc11=(raa*(1.d0-fs)+rmix)*SLAP**2.d0/c1
	if (xka.ne.0.d0) cc11=cc11-( sata*(1.d0-fs)/(c1*xka)+cgg/xka )*SLAP
!
	cc12=DCMPLX (0.d0,0.d0)
	if (xka.ne.0.d0) cc12=cc12-rmix*cgg*SLAP**3.d0/(c1*xka)
!
	cc21=1.d0/(2.d0*pi)
!
	cc13=DCMPLX (0.d0,0.d0)
	if (xka.ne.0.d0) cc13=SLAP/(2.d0*pi*xka*c1)
!
	cc14=DCMPLX (0.d0,0.d0)
	cc14=-c1*cwg+(raa*xka*SLAP-sata)*fs
!
	cc15=DCMPLX (0.d0,0.d0)
	cc15=rmix*cwg*SLAP**2.d0
!
	cc16=DCMPLX (0.d0,0.d0)
	cc16=( rw*fs/c1+rmix/c1 )*SLAP**2.d0
	if (xkw.ne.0.d0) cc16=cc16-( sat*fs/(c1*xkw)+cww/xkw )*SLAP
!
	cc17=DCMPLX (0.d0,0.d0)
	if (xkw.ne.0.d0) cc17=-rmix*cww*SLAP**3.d0/(c1*xkw)
!
	cc22=1.d0/(2.d0*pi)
!
	cc18=DCMPLX (0.d0,0.d0)
	if (xkw.ne.0.d0.and.xka.ne.0.d0) cc18=SLAP/(2.d0*pi*xkw*c1)
!
	cc19=DCMPLX (0.d0,0.d0)
	cc19=-c1*cwg+(rw*xkw*SLAP-sat)*(1.d0-fs)
!
	cc20=DCMPLX (0.d0,0.d0)
	cc20=rmix*cwg*SLAP**2.d0
!

	xKss1=DCMPLX (0.d0,0.d0)
	xKss1=(rw*fs+raa*(1.d0-fs))*SLAP**2.d0/(xlambda+xmu)
	if (xka.ne.0.d0) xKss1=xKss1-sata*(1.d0-fs)*SLAP/(xka*(xlambda+xmu))
	if (xkw.ne.0.d0) xKss1=xKss1-sat*fs*SLAP/(xkw*(xlambda+xmu))
	if (xka.ne.0.d0.and.xkw.ne.0.d0) xKss1=xKss1-(cgg*xkw+cww*xka)*SLAP/(xkw*xka)
!
	xKss2=DCMPLX (0.d0,0.d0)
	if (xka.ne.0.d0) xKss2=xKss2-rw*(cgg*fs-cwg*(1.-fs))*SLAP**3.d0/(xka*(xlambda+xmu))
	if (xkw.ne.0.d0) xKss2=xKss2-raa*(-cwg*fs+cww*(1.d0-fs))*SLAP**3.d0/(xkw*(xlambda+xmu))
	if (xka.ne.0.d0.and.xkw.ne.0.d0) xKss2=xKss2+( (-cwg**2.d0+cww*cgg)+&
										  sat*(-cwg*(1.d0-fs)+cgg*fs)/(xlambda+xmu)+&
										  sata*(-cwg*fs+cww*(1.d0-fs))/(xlambda+xmu) )*&
										  SLAP**2.d0/(xkw*xka)
!
	do i=1,4
		CsS(i)=0.d0
		i1=i
		i2=i+1
		i3=i+2
		i4=i+3
		if (i2.gt.4) i2=i2-4
		if (i3.gt.4) i3=i3-4
		if (i4.gt.4) i4=i4-4
		CsS(i)=cc1*cc2*(RT(i)**2.d0-xKss1*RT(i)+xKss2)/((RT(i1)-RT(i2))*(RT(i1)-RT(i3))*&
												  (RT(i1)-RT(i4)))
	enddo
!
	do i=1,3
		CsW(i)=0.d0;CsA(i)=0.d0;CwS(i)=0.d0;CaS(i)=0.d0;CCwW(i)=0.d0
		CwA(i)=0.d0;CaA(i)=0.d0;CaW(i)=0.d0
		i1=i+1
		i2=i+2
		i3=i+3
		if (i2.gt.4) i2=i2-3
		if (i3.gt.4) i3=i3-3
		CsW(i)=cc3*(RT(i1)-cc4)/( (RT(i3)-RT(i1))*(RT(i2)-RT(i1)) )
		CsA(i)=cc5*(RT(i1)-cc6)/( (RT(i3)-RT(i1))*(RT(i2)-RT(i1)) )
		CwS(i)=cc7*(RT(i1)-cc8)/( (RT(i3)-RT(i1))*(RT(i2)-RT(i1)) )
		CaS(i)=cc9*(RT(i1)-cc10)/( (RT(i3)-RT(i1))*(RT(i2)-RT(i1)) )
		CCwW(i)=cc21*(RT(i1)**2.d0-cc11*RT(i1)+cc12)/( (RT(i3)-RT(i1))*(RT(i2)-RT(i1)) )
		CwA(i)=cc13*(cc14*RT(i1)+cc15)/( (RT(i3)-RT(i1))*(RT(i2)-RT(i1)) )
		CaA(i)=cc22*(RT(i1)**2.d0-cc16*RT(i1)+cc17)/( (RT(i3)-RT(i1))*(RT(i2)-RT(i1)) )
		CaW(i)=cc18*(cc19*RT(i1)+cc20)/( (RT(i3)-RT(i1))*(RT(i2)-RT(i1)) )
	enddo
!
	return
	end


!*******************************************************************************************
!
	subroutine WFSOLGCOEFQST(SLAP,CK11,CK12,CK13,CK21,CK22,CK23,CK31,CK32,CK41,CK42,CK51,CK52,&
							 CK61,CK62,CK71,CK72,CK73,CK74,CK75,CK76,CK77)
!
	implicit double precision (a-h,o-z)
!
	COMPLEX(8) SLAP,c11,c12,c13,c14,c21,c22,c23,c24,c25,c31,c32,c33,c34,D1,D11,F11,F12,F13,&
			   F21,F22,F23,F31,F32,F41,F42,F51,F52,F61,F62,F71,F72,F73,F74,F75,F76,F77,&
			   CK11,CK12,CK13,CK21,CK22,CK23,&
			   CK31,CK32,CK41,CK42,CK51,CK52,CK61,CK62,CK71,CK72,&
			   CK73,CK74,CK75,CK76,CK77

	common /e/ gammaw,gammas,gammaa,grav,atmp
	common /point/ dist(2),ra,rd(2),rdn,eta(2)
	common /umate/ d(3,3),dfs(3),xn,xkw,cww,cwg,xka,cgg,cgw,ev,sat,g1,rs,rw,raa,rmix,xm1,&
				   bt,et,xlambda,xmu,henry
	common /UcoefR/ r1(2,2),dpr1(2,2,2),r2(2,2),dpr2(2,2,2),c(2,2)
!
!
	CK11=DCMPLX (0.d0,0.d0)
	CK12=DCMPLX (0.d0,0.d0)
	CK13=DCMPLX (0.d0,0.d0)
	CK21=DCMPLX (0.d0,0.d0)
	CK22=DCMPLX (0.d0,0.d0)
	CK23=DCMPLX (0.d0,0.d0)
	CK31=DCMPLX (0.d0,0.d0)
	CK32=DCMPLX (0.d0,0.d0)
	CK41=DCMPLX (0.d0,0.d0)
	CK42=DCMPLX (0.d0,0.d0)
	CK51=DCMPLX (0.d0,0.d0)
	CK52=DCMPLX (0.d0,0.d0)
	CK61=DCMPLX (0.d0,0.d0)
	CK62=DCMPLX (0.d0,0.d0)
	CK71=DCMPLX (0.d0,0.d0)
	CK72=DCMPLX (0.d0,0.d0)
	CK73=DCMPLX (0.d0,0.d0)
	CK74=DCMPLX (0.d0,0.d0)
	CK75=DCMPLX (0.d0,0.d0)
	CK76=DCMPLX (0.d0,0.d0)
	CK77=DCMPLX (0.d0,0.d0)
!
	pi=DACOS (-1.d0)
	zmu=3.d0*bt*et/(9.d0*bt-et)
	zlamda=bt-2.d0*zmu/3.d0
	znu=zlamda/(2.d0*(zlamda+zmu))
!
	c11=zlamda+zmu
	c12=zmu
	c13=1.d0-dfs(1) !****
	c14=dfs(1)      !****
	c21=SLAP*raa*( 1.d0-sat*(1.d0-henry) )
	c22=-SLAP*raa*g1*xn*(1.d0-henry)
	c23=-raa*xka/gammaa
	c24=SLAP*raa*g1*xn*(1.d0-henry)
	c25=-henry*raa*xkw/gammaw
	c31=SLAP*rw*sat
	c32=SLAP*rw*g1*xn
	c33=-rw*xkw/gammaw
	c34=-SLAP*rw*g1*xn
	D1=c12*(c11+c12)*c23*c33
	D11=1.d0/(2.d0*pi*D1)
!
!	--- Fij coefficients ---
!
!	F11=(c11+c12)*c23*c33
!	F12=( -c14*c23*c31+c13*(c25*c31-c21*c33)-(c11+c12)*(c25*c32-c22*c33-c23*c34) )
!	F13=( c14*(c21*c32-c22*c31)+c13*(c24*c31-c21*c34)-(c11+c12)*(c24*c32-c22*c34) )
!	F21=-c11*c23*c33
!	F22=( c14*c23*c31+c13*(c21*c33-c25*c31)+c11*(c25*c32-c22*c33-c23*c34) )
!	F23=( c14*(c22*c31-c21*c32)+c13*(c21*c34-c24*c31)+c11*(c24*c32-c22*c34) )
	F31=c12*c13*c33            !-
	F32=-c12*(c14*c32-c13*c34) !+
	F41=-c12*(c13*c25-c14*c23) !+
	F42=-c12*(c13*c24-c14*c22) !+
	F51=-c12*(c25*c31-c21*c33) !+
	F52=-c12*(c24*c31-c21*c34) !+
	F61=c12*c23*c31            !-
	F62=-c12*(c21*c32-c22*c31) !+
	F71=c12*(c11+c12)*c33
	F72=c12*(-c14*c31+(c11+c12)*c34)
	F73=-c12*(c11+c12)*c25
	F74=-c12*(-c14*c21+(c11+c12)*c24)
	F75=-c12*(-c13*c31+(c11+c12)*c32)
	F76=c12*(c11+c12)*c23
	F77=c12*(-c13*c21+(c11+c12)*c22)
!
!	--- Kij coefficients ---
!
!	K11=F11*D11
!	K12=F12*D11/SLAP
!	K13=F13*D11/SLAP**2.d0
!	K21=F21*D11
!	K22=F22*D11/SLAP
!	K23=F23*D11/SLAP**2.d0
	CK31=F31*D11
	CK32=F32*D11/SLAP
	CK41=F41*D11
	CK42=F42*D11/SLAP
	CK51=F51*D11/SLAP
	CK52=F52*D11/SLAP**2.d0
	CK61=F61*D11/SLAP
	CK62=F62*D11/SLAP**2.d0
	CK71=F71*D11
	CK72=F72*D11/SLAP
	CK73=F73*D11
	CK74=F74*D11/SLAP
	CK75=F75*D11/SLAP
	CK76=F76*D11
	CK77=F77*D11/SLAP
!
!	CK11=1.D0/(2.D0*pi*zmu)
!	CK12=( xn*g1*(zlamda+2.D0*zmu)*(xka*gammaw+xkw*gammaa)+xkw*gammaa*(1.d0-sat)*(1.d0-dfs(1))+&
!		 xka*gammaw*sat*dfs(1) )/&
!		(2.D0*pi*zmu*(zlamda+2.D0*zmu)*xkw*xka)
!!	CK13=xn*g1*gammaw*gammaa/(2.D0*pi*zmu*(zlamda+2.D0*zmu)*xkw*xka)
!	CK21=-(zlamda+zmu)/(2.D0*pi*zmu*(zlamda+2.D0*zmu))
!	CK22=-( xn*g1*(zlamda+zmu)*(xka*gammaw+xkw*gammaa)+xkw*gammaa*(1.d0-sat)*(1.d0-dfs(1))+&
!		 xka*gammaw*sat*dfs(1) )/&
!		(2.D0*pi*zmu*(zlamda+2.D0*zmu)*xkw*xka)
!	CK23=-CK13
!



!	R1ij & dpR1ijm
!
	temp=ra**3.d0
	r1(1,1)=(2.d0*dist(1)**2.d0-ra**2.d0)/temp
	r1(1,2)=(2.d0*dist(1)*dist(2))/temp
	r1(2,1)=r1(1,2)
	r1(2,2)=(2.d0*dist(2)**2.d0-ra**2.d0)/temp
	temp=temp*ra
	dpr1(1,1,1)=(4.d0*ra*dist(1)-6.d0*rd(1)*dist(1)**2.d0+ra**2.d0*rd(1))/temp
	dpr1(1,2,1)=(2.d0*ra*dist(2)-6.d0*rd(1)*dist(1)*dist(2))/temp
	dpr1(2,1,1)=dpr1(1,2,1)
	dpr1(2,2,1)=(-6.d0*rd(1)*dist(2)**2.d0+ra**2.d0*rd(1))/temp
	dpr1(2,2,2)=(4.d0*ra*dist(2)-6.d0*rd(2)*dist(2)**2.d0+ra**2.d0*rd(2))/temp
	dpr1(1,2,2)=(2.d0*ra*dist(1)-6.d0*rd(2)*dist(1)*dist(2))/temp
	dpr1(2,1,2)=dpr1(1,2,2)
	dpr1(1,1,2)=(-6.d0*rd(2)*dist(1)**2.d0+ra**2.d0*rd(2))/temp
!
!	R2ij & dpR2ijm
	temp=ra**2.d0
	r2(1,1)=dist(1)**2.d0/temp
	r2(2,2)=dist(2)**2.d0/temp
	r2(1,2)=dist(1)*dist(2)/temp
	r2(2,1)=r1(1,2)
	temp=temp*ra
	dpr2(1,1,1)=(2.d0*ra*dist(1)-2.d0*rd(1)*dist(1)**2.d0)/temp
	dpr2(1,2,1)=(ra*dist(2)-2.d0*rd(1)*dist(1)*dist(2))/temp
	dpr2(2,1,1)=dpr2(1,2,1)
	dpr2(2,2,1)=(-2.d0*rd(1)*dist(2)**2.d0)/temp
	dpr2(2,2,2)=(2.d0*ra*dist(2)-2.d0*rd(2)*dist(2)**2.d0)/temp
	dpr2(1,2,2)=(ra*dist(1)-2.d0*rd(2)*dist(1)*dist(2))/temp
	dpr2(2,1,2)=dpr2(1,2,2)
	dpr2(1,1,2)=(-2.d0*rd(2)*dist(1)**2.d0)/temp
!
	return
	end


!*******************************************************************************************
!
	subroutine wsatDYNfsolF(SLAP,RTSD,RTSDI,K0,K1,KI0,KI1,DK0,DK1,DKI0,DKI1,cqmG,cqmF)
!
	implicit double precision (a-h,o-z)
!
	COMPLEX(8) SLAP,RTSD(4),RTSDI(2),DRTSD(4),DRTSDI(2),K0(3),K1(3),KI0(2),KI1(2),&
			   DK0(3),DK1(3),DKI0(2),DKI1(2),R7,temp,temp1,temp2,&
			   cqmG(kbem,kbem),cqmF(kbem,kbem),F12,F21,F13,F23,F31,F32,XKSI,zalphaT,&
			   dpsolG(2,2,2),beta,RN(2),DRN(2),RL4
!
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /e/ gammaw,gammas,gammaa,grav,atmp
	common /point/ dist(2),ra,rd(2),rdn,eta(2)
	common /coefL/aa(2,2),dpaa(2,2,2),bb(2,2),dpbb(2,2,2),cc(2,2)
	common /mate/zlamda,zmuy,znuu,znu,zalpha,zrho,zrhof,zm,zk,zr,zn
	common /ocqm/ EPScqm,taw,Icqm,Lcqm,NPcqm
	common /comp/ icomp
!
	dimension dlt(2,2),abscqmF(kbem,kbem),dpdd(2,2,2)
!
!
	do i=1,kbem
		do j=1,kbem
			cqmF(i,j)=0.d0
		enddo
	enddo
!
	pi=DACOS (-1.d0)
	XKSI=zrho*SLAP*SLAP/(zlamda+2.d0*zmuy)
	zalphaT=zalpha-SLAP*zrhof*zk
	beta=zk*zrhof*(zn**2.d0)*SLAP*SLAP/(SLAP*zn**2.d0+zk*zn*zrhof*SLAP*SLAP)
	RL4=(zrho-beta*zrhof)/(zlamda+2.d0*zmuy)
!
	dlt(1,1)=1.d0
	dlt(2,2)=1.d0
	dlt(1,2)=0.d0
	dlt(2,1)=0.d0
!
	do i=1,4
		DRTSD(i)=Cdsqrt( RTSD(i) )
	enddo
	do i=1,2
		DRTSDI(i)=Cdsqrt( RTSDI(i) )
	enddo
!
	if (icomp.eq.0) goto 1000 !INCOMPRESSIBLE F.S
!
!	COMPRESSIBLE F.S
!
	do k=1,3
		do j=1,3
			temp=DCMPLX (0.d0,0.d0)
			if (k.lt.3.and.j.eq.3) GOTO 10
			if (k.eq.3.and.j.lt.3) GOTO 30
			if (k.eq.3.and.j.eq.3) GOTO 40
!
!	Tij
!
			temp1=0.d0
			temp2=0.d0
			temp1=rd(K)*eta(J)*( RTSD(1)*DRTSD(1)*K1(1)*(RTSD(2)-RTSD(4))-&
								 RTSD(2)*DRTSD(2)*K1(2)*(RTSD(1)-RTSD(4)) )/&
								 ( 2.D0*pi*zrho*(SLAP*SLAP)*(RTSD(1)-RTSD(2)) )
			temp1=temp1*zlamda+zalpha*SLAP*zalpha*rd(K)*eta(J)*( DRTSD(1)*K1(1)/(RTSD(1)-RTSD(2))+&
							   DRTSD(2)*K1(2)/(RTSD(2)-RTSD(1)) )/(2.D0*pi*zk*(zlamda+2.d0*zmuy))
!
			R7=rdn*( dlt(k,j)-4.d0*rd(j)*rd(k) )+rd(K)*eta(J)+rd(J)*eta(K)
!
			temp2=(1.d0/pi)*(&
					  ( (RTSD(4)-RTSD(2))/( RTSD(3)*(RTSD(1)-RTSD(2)) ))*&
							( R7*DRTSD(1)*(DRTSD(1)*K0(1)+2.d0*K1(1)/ra)/ra-rd(j)*rd(k)*rdn*&
							  RTSD(1)*DRTSD(1)*K1(1) )-&
					  ( (RTSD(4)-RTSD(1))/( RTSD(3)*(RTSD(1)-RTSD(2)) ))*&
							( R7*DRTSD(2)*(DRTSD(2)*K0(2)+2.d0*K1(2)/ra)/ra-rd(j)*rd(k)*rdn*&
							  RTSD(2)*DRTSD(2)*K1(2) )-&
					   R7*(DRTSD(3)*K0(3)+2.d0*K1(3)/ra)/( DRTSD(3)*ra )-&
					   0.5d0*( rdn*(dlt(k,j)-2.d0*rd(k)*rd(j))+rd(J)*eta(K) )*DRTSD(3)*K1(3)  )
!
			temp=temp2+temp1
			cqmF(k,j)=temp
			goto 100
!
!	Ti3
!
10			temp=( rd(k)*rdn*zalpha*RTSD(4)*(RTSD(2)*K0(2)-RTSD(1)*K0(1))/(RTSD(1)-RTSD(2))+&
				   (2.d0*rd(k)*rdn-eta(k))*zalpha*RTSD(4)*(DRTSD(2)*K1(2)-DRTSD(1)*K1(1))/&
				   (ra*(RTSD(1)-RTSD(2))) )/(2.D0*pi*zrho*SLAP*SLAP)
!
			cqmF(k,j)=temp
			goto 100
!
!	T3j
!
30			temp=( 2.d0*rd(j)*rdn*SLAP*(zalpha-SLAP*zrhof*zk)*zmuy*&
							( (RTSD(2)*K0(2)+DRTSD(2)*K1(2)/ra)-(RTSD(1)*K0(1)+DRTSD(1)*K1(1)/ra) )-&
				   2.d0*(eta(j)-rd(j)*rdn)*SLAP*(zalpha-SLAP*zrhof*zk)*zmuy*&
											   (DRTSD(2)*K1(2)/ra-DRTSD(1)*K1(1)/ra)+&
				   eta(j)*( zlamda*SLAP*(zalpha-SLAP*zrhof*zk)*RTSD(2)-&
							zalpha*SLAP*(zlamda+2.d0*zmuy)*(RTSD(2)-RTSD(4)) )*K0(2)-&
				   eta(j)*( zlamda*SLAP*(zalpha-SLAP*zrhof*zk)*RTSD(1)-&
							zalpha*SLAP*(zlamda+2.d0*zmuy)*(RTSD(1)-RTSD(4)) )*K0(1) )/&
				   (2.D0*pi*zk*(RTSD(1)-RTSD(2))*(zlamda+2.d0*zmuy))
!
			cqmF(k,j)=temp
			goto 100
!
!	T33
!
40			temp=rdn*( (RTSD(2)-RTSD(4))*DRTSD(2)*K1(2)-(RTSD(1)-RTSD(4))*DRTSD(1)*K1(1) )/&
					   ( 2.D0*pi*(RTSD(1)-RTSD(2)) )
!
			cqmF(k,j)=temp
			goto 100
!
100			continue
!
		enddo
	enddo
	goto 1
!
!	INCOMPRESSIBLE F.S
!
1000 continue
!
!	--- dpGij ---
!
	do i=1,2
		do j=1,2
			do m=1,2
				dpdd(i,j,m)=(-4.D0*rd(i)*rd(j)*rd(m)+rd(i)*dlt(j,m)+rd(j)*dlt(i,m)+rd(m)*dlt(i,j))/ra
				DRN(1)=2.D0*dpdd(i,j,m)*DRTSDI(1)*SLAP*KI1(1)/ra**2.D0-&
					   (2.D0*rd(i)*rd(j)*rd(m)-rd(m)*dlt(i,j))*DRTSDI(1)*SLAP*KI1(1)/ra**2.D0+&
					   dpdd(i,j,m)*SLAP*SLAP*RTSDI(1)*KI0(1)-&
					   rd(i)*rd(j)*rd(m)*SLAP*SLAP*SLAP*RTSDI(1)*DRTSDI(1)*KI1(1)
!
				DRN(2)=2.D0*dpdd(i,j,m)*DRTSDI(2)*SLAP*KI1(2)/ra**2.D0-&
					   (2.D0*rd(i)*rd(j)*rd(m)-rd(m)*dlt(i,j))*DRTSDI(2)*SLAP*KI1(2)/ra**2.D0+&
					   dpdd(i,j,m)*SLAP*SLAP*RTSDI(2)*KI0(2)-&
					   rd(i)*rd(j)*rd(m)*SLAP*SLAP*SLAP*RTSDI(2)*DRTSDI(2)*KI1(2)
!
				temp=DCMPLX (0.d0,0.d0)
				temp=temp+( RL4*DRN(1)/RTSDI(1)-&
				            2.D0*(RL4-RTSDI(1))*dpdd(i,j,m)/(RTSDI(1)*ra**2.D0)-DRN(2)-&
							dlt(i,j)*rd(m)*SLAP*SLAP*SLAP*RTSDI(2)*DRTSDI(2)*KI1(2) )/&
						  (2.D0*pi*SLAP*SLAP*(zrho-beta*zrhof))
!
!				temp=temp+( dpaa(i,j,m)*KI1(1)/DRTSDI(1)+aa(i,j)*rd(m)*DKI1(1)/DRTSDI(1)+&
!					dpbb(i,j,m)*KI0(1)+bb(i,j)*rd(m)*DKI0(1) )*XKSI/(zrho*SLAP**2.d0)-&
!					dpaa(i,j,m)*( (XKSI-RTSDI(1))/(ra*RTSDI(1)*zrho*SLAP**2.d0) )+&
!					aa(i,j)*rd(m)*( (XKSI-RTSDI(1))/(ra**2.d0*RTSDI(1)*zrho*SLAP**2.d0) )-&
!					( dpaa(i,j,m)*KI1(2)/DRTSDI(2)+aa(i,j)*rd(m)*DKI1(2)/DRTSDI(2)+&
!					dpbb(i,j,m)*KI0(2)+bb(i,j)*rd(m)*DKI0(2) )*RTSDI(2)/(zrho*SLAP**2.d0)+&
!					cc(i,j)*rd(m)*DKI0(2)
!
				dpsolG(i,j,m)=temp
			enddo
		enddo
	enddo
!

!
	do i=1,3
		do j=1,3
			temp=DCMPLX (0.d0,0.d0)
			if (i.lt.3.and.j.eq.3) GOTO 50
			if (i.eq.3.and.j.lt.3) GOTO 60
			if (i.eq.3.and.j.eq.3) GOTO 70
!
!	Tij
!
			do m=1,2
				temp=temp+(dpsolG(i,j,m)+dpsolG(m,j,i))*eta(m)
			enddo
			dila1=dpsolG(1,j,1)+dpsolG(2,j,2)
			temp=temp*zmuy+(dila1*zlamda-SLAP*zalpha*cqmG(i,3))*eta(i)
!
			cqmF(i,j)=temp
			goto 200
!
!	Ti3
!
50			RN(1)=(2.d0*rdn*rd(i)-eta(i))*DRTSDI(1)*SLAP*KI1(1)/ra+rdn*rd(i)*RTSDI(1)*KI0(1)*SLAP*SLAP
			RN(2)=(2.d0*rdn*rd(i)-eta(i))*DRTSDI(2)*SLAP*KI1(2)/ra+rdn*rd(i)*RTSDI(2)*KI0(2)*SLAP*SLAP
			temp=beta*(1.d0-beta)*( ((eta(i)-2.D0*rdn*rd(i))/ra)*(DRTSDI(1)*SLAP*KI1(1)-1.d0/ra)-&
								      RTSDI(1)*SLAP*SLAP*rdn*rd(i)*KI0(1) )/&
								  (2.D0*pi*(zlamda+2.d0*zmuy)*RTSDI(1)*SLAP*SLAP)-&
				 beta*(RL4*RN(1)/RTSDI(1)-(RL4-RTSDI(1))*(2.d0*rdn*rd(i)-eta(i))/(RTSDI(1)*ra**2.d0)-&
				        RN(2)+eta(i)*RTSDI(2)*KI0(2)*SLAP*SLAP)/(2.D0*pi*SLAP*SLAP*(zrho-beta*zrhof))
!
!			temp=zalpha*(&
!						((eta(i)-2.d0*rdn*rd(i))/ra)*(KI1(1)/DRTSDI(1)-1.d0/(RTSDI(1)*ra))-&
!						  rdn*rd(i)*KI0(1) )/(2.D0*pi*(zlamda+2.d0*zmuy))
!
			cqmF(i,j)=temp
			goto 200
!
!	T3j
!
60			temp=zrhof*SLAP*SLAP*eta(j)*(-(1.d0-beta)*zlamda*KI0(1)/(zlamda+2.d0*zmuy)+&
				                          zalpha*((RTSDI(1)-RL4)*KI0(1)-RL4*dlog(ra))/&
				 RTSDI(1))/(2.D0*pi*beta)+zmuy*zrhof*(1.d0-beta)*&
				 ( ((eta(j)-2.D0*rdn*rd(j))/ra)*(DRTSDI(1)*SLAP*KI1(1)-1.d0/ra)-&
				   RTSDI(1)*SLAP*SLAP*rd(j)*rdn*KI0(1) )/(pi*beta*RTSDI(1)*(zlamda+2.d0*zmuy))
!
!			temp=zmuy*zalphaT*SLAP*(&
!						 ((eta(j)-2.d0*rd(j)*rdn)/ra)*( 2.d0*KI1(1)/DRTSDI(1)-2.d0/(RTSDI(1)*ra) )-&
!						   2.d0*rd(j)*rdn*KI0(1) )/(2.D0*pi*zk*(zlamda+2.d0*zmuy))+&
!				 zalpha*SLAP*eta(j)*(&
!				        zalphaT*SLAP*zalpha*KI0(1)/(zk*(zlamda+2.d0*zmuy))-XKSI*dlog(ra) )-&
!				 zlamda*zalphaT*SLAP*eta(j)*KI0(1)/(2.D0*pi*zk*(zlamda+2.d0*zmuy))
!
			cqmF(i,j)=temp
			goto 200
!
!	T33
!
70			temp=-rdn*( -(RTSDI(1)-RL4)*DKI1(1)+RL4/ra+zrhof**2.d0*(1.d0-beta)*SLAP*&
				         (DRTSDI(1)*SLAP*KI1(1)-1.d0/ra)/(zlamda+2.d0*zmuy) )/(2.D0*pi*RTSDI(1))
!			temp=-rdn*( DRTSDI(1)*(RTSDI(1)-XKSI)*KI1(1)+XKSI/ra )/(2.D0*pi*RTSDI(1))
!
			cqmF(i,j)=temp
			goto 200
!
200			continue
!
		enddo
	enddo
!
!
!
1	return
	end subroutine
!
!*******************************************************************************************

!*******************************************************************************************
!
	subroutine wsatDYNfsolG(SLAP,RTSD,RTSDI,K0,K1,KI0,KI1,cqmG)
!
	implicit double precision (a-h,o-z)
!
	COMPLEX(8) SLAP,RTSD(4),DRTSD(4),RTSDI(2),DRTSDI(2),K0(3),K1(3),KI0(2),KI1(2),&
			   R(3),temp,cqmG(kbem,kbem),G13,G23,G31,G32,G33,XKSI,zalphaT,beta,RL4
!
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /e/ gammaw,gammas,gammaa,grav,atmp
	common /point/ dist(2),ra,rd(2),rdn,eta(2)
	common /coefL/aa(2,2),dpaa(2,2,2),bb(2,2),dpbb(2,2,2),cc(2,2)
	common /mate/zlamda,zmuy,znuu,znu,zalpha,zrho,zrhof,zm,zk,zr,zn
	common /ocqm/ EPScqm,taw,Icqm,Lcqm,NPcqm
	common /comp/ icomp
!
	dimension dlt(2,2),abscqmG(kbem,kbem)
!
!
	do i=1,kbem
		do j=1,kbem
			cqmG(i,j)=0.d0
		enddo
	enddo
!
	pi=DACOS (-1.d0)
	XKSI=zrho*(SLAP*SLAP)/(zlamda+2.d0*zmuy)
	zalphaT=zalpha-SLAP*zrhof*zk
	beta=zk*zrhof*(zn**2.d0)*SLAP*SLAP/(SLAP*zn**2.d0+zk*zn*zrhof*SLAP*SLAP)
	RL4=(zrho-beta*zrhof)/(zlamda+2.d0*zmuy)
!
	dlt(1,1)=1.d0
	dlt(2,2)=1.d0
	dlt(1,2)=0.d0
	dlt(2,1)=0.d0
!
	do i=1,4
		DRTSD(i)=Cdsqrt( RTSD(i) )
	enddo
	do i=1,2
		DRTSDI(i)=Cdsqrt( RTSDI(i) )
	enddo
!
	if (icomp.eq.0) goto 1000 !INCOMPRESSIBLE F.S
!
!	call FS(RTSD,DRTSD,K0,K1)
!
!	COMPRESSIBLE F.S
!
	do k=1,3
		do j=1,3
			temp=DCMPLX (0.d0,0.d0)
			if (k.lt.3.and.j.eq.3) GOTO 10
			if (k.eq.3.and.j.lt.3) GOTO 30
			if (k.eq.3.and.j.eq.3) GOTO 40
!
!	Gij
!
			do i=1,3
				R(i)=(2.d0*rd(k)*rd(j)-dlt(k,j))*DRTSD(i)*K1(i)/ra+rd(k)*rd(j)*RTSD(i)*K0(i)
			enddo
			temp=( ( (RTSD(4)-RTSD(2))/(RTSD(1)-RTSD(2)) )*R(1) -&
				   ( (RTSD(4)-RTSD(1))/(RTSD(1)-RTSD(2)) )*R(2) +&
				   dlt(k,j)*RTSD(3)*K0(3)-R(3) )/(2.D0*pi*zrho*SLAP*SLAP)
			cqmG(k,j)=temp
			goto 100
!
!	Gi3
!
10			temp=-zalpha*rd(k)*( DRTSD(1)*K1(1)/(RTSD(1)-RTSD(2))+&
								 DRTSD(2)*K1(2)/(RTSD(2)-RTSD(1)) )/(2.D0*pi*zk*(zlamda+2.d0*zmuy))
			cqmG(k,j)=temp
			goto 100
!
!	G3j
!
30			temp=SLAP*zalphaT*rd(j)*( DRTSD(1)*K1(1)/(RTSD(1)-RTSD(2))+&
									  DRTSD(2)*K1(2)/(RTSD(2)-RTSD(1)) )/(2.D0*pi*zk*(zlamda+2.d0*zmuy))
			cqmG(k,j)=temp
			goto 100
!
!	G33
!
40			temp= -( K0(1)*(RTSD(1)-RTSD(4))-K0(2)*(RTSD(2)-RTSD(4)) )/(2.D0*pi*zk*(RTSD(1)-RTSD(2)))
			cqmG(k,j)=temp
			goto 100
!
100			continue
!
		enddo
	enddo
	goto 1
!
!	INCOMPRESSIBLE F.S
!
1000 continue
!
	do i=1,3
		do j=1,3
			temp=DCMPLX (0.d0,0.d0)
			if (i.lt.3.and.j.eq.3) GOTO 50
			if (i.eq.3.and.j.lt.3) GOTO 60
			if (i.eq.3.and.j.eq.3) GOTO 70
!
!	Gij
!
			do k=1,2
				R(k)=(2.d0*rd(i)*rd(j)-dlt(i,j))*DRTSDI(k)*SLAP*KI1(k)/ra+rd(i)*rd(j)*RTSDI(k)*KI0(k)*SLAP*SLAP
			enddo
			temp=(RL4*R(1)/RTSDI(1)-(RL4-RTSDI(1))*(2.d0*rd(i)*rd(j)-dlt(i,j))/(RTSDI(1)*ra**2.d0)-&
				  R(2)+dlt(i,j)*RTSDI(2)*KI0(2)*SLAP*SLAP)/(2.D0*pi*SLAP*SLAP*(zrho-beta*zrhof))
!
!			temp=( aa(i,j)*KI1(1)/DRTSDI(1)+bb(i,j)*KI0(1) )*XKSI/(zrho*SLAP**2.d0)-&
!				 ( aa(i,j)*(XKSI-RTSDI(1))/(ra*RTSDI(1)*zrho*SLAP**2.d0) )-&
!				 ( aa(i,j)*KI1(2)/DRTSDI(2)+bb(i,j)*KI0(2) )*RTSDI(2)/(zrho*SLAP**2.d0)+&
!				   cc(i,j)*KI0(2)
!
			cqmG(i,j)=temp
			goto 200
!
!	Gi3
!
50			temp=-rd(i)*zrhof*((1.d0-beta)*(DRTSDI(1)*SLAP*ra*KI1(1)-1.d0))/(2.D0*pi*SLAP*beta*&
																			 RTSDI(1)*ra*(zlamda+2.d0*zmuy))
!			temp=-zalpha*rd(i)*( KI1(1)/DRTSDI(1)-1.d0/(RTSDI(1)*ra) )/&
!				 (2.D0*pi*zk*(zlamda+2.d0*zmuy))
!
			cqmG(i,j)=temp
			goto 200
!
!	G3j
!
60			temp=-SLAP*cqmG(j,3)
!			temp=(zalpha-SLAP*zrhof*zk)*SLAP*rd(j)*( KI1(1)/DRTSDI(1)-1.d0/(RTSDI(1)*ra) )/&
!								                   (2.D0*pi*zk*(zlamda+2.d0*zmuy))
!
			cqmG(i,j)=temp
			goto 200
!
!	G33
!
70			temp=-SLAP*zrhof*((RTSDI(1)-RL4)*KI0(1)-RL4*dlog(ra))/(2.D0*pi*beta*RTSDI(1))
!			temp=-( (RTSDI(1)-XKSI)*KI0(1)-XKSI*dlog(ra) )/(2.D0*pi*zk*RTSDI(1))
!
			cqmG(i,j)=temp
			goto 200
!
200			continue
!
		enddo
	enddo
!
1	return
	end subroutine
!
!*******************************************************************************************

!*******************************************************************************************
!
	subroutine wsatqstfsolF(SLAP,RTSQ,KS0,KS1,KS2,KS3,cqmF)
!
	implicit double precision (a-h,o-z)
!
	COMPLEX(8) SLAP,RTSQ,DRTSQ,KS0,KS1,KS2,KS3,XKSI,temp,cqmF(kbem,kbem)
!
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /e/ gammaw,gammas,gammaa,grav,atmp
	common /point/ dist(2),ra,rd(2),rdn,eta(2)
	common /UcoefR/ r1(2,2),dpr1(2,2,2),r2(2,2),dpr2(2,2,2),c(2,2)
	common /mate/zlamda,zmuy,znuu,znu,zalpha,zrho,zrhof,zm,zk,zr,zn
	common /ocqm/ EPScqm,taw,Icqm,Lcqm,NPcqm
!
	dimension dlt(2,2),abscqmF(kbem,kbem)
!
!
	do i=1,kbem
		do j=1,kbem
			cqmF(i,j)=0.d0
		enddo
	enddo
!
	pi=DACOS (-1.d0)
!
	dlt(1,1)=1.d0
	dlt(2,2)=1.d0
	dlt(1,2)=0.d0
	dlt(2,1)=0.d0
!
	DRTSQ=Cdsqrt( RTSQ )
	XKSI=DRTSQ*ra
!
	do k=1,3
		do j=1,3
			temp=DCMPLX (0.d0,0.d0)
			if (k.lt.3.and.j.eq.3) GOTO 10
			if (k.eq.3.and.j.lt.3) GOTO 30
			if (k.eq.3.and.j.eq.3) GOTO 40
!
!	Fij
!
			temp=( (1.D0-2.D0*ZNUU)/(2.D0*(1.D0-ZNUU))*( eta(j)*rd(k)-eta(k)*rd(j)-dlt(k,j)*rdn )/ra-&
				             rd(k)*rd(j)*rdn/((1.D0-ZNUU)*ra)+&
							 (ZNUU-ZNU)*DRTSQ/((1.D0-ZNUU)*(1.D0-ZNU))*( eta(j)*rd(k)*(KS3-3.d0*KS2/XKSI-2.D0/XKSI**(3.d0))+&
																		(eta(k)*rd(j)+dlt(k,j)*rdn)*(KS2/XKSI-2.D0/XKSI**(3.d0))+&
																		 rd(k)*rd(j)*rdn*(8.D0/XKSI**(3.d0)-KS3) ))/(2.D0*pi)
			cqmF(k,j)=temp
			goto 100
!
!	Fi3
!
10			temp=zalpha*(1.D0-2.D0*ZNU)/(4.d0*pi*ZMUY*(1.D0-ZNU))*( rd(k)*rdn*(2.D0/XKSI**(2.d0)-KS2)+eta(k)*(KS1/XKSI-1.D0/XKSI**(2.d0)) )
			cqmF(k,j)=temp
			goto 100
!
!	F3j
!
30			temp=SLAP*zalpha*(1.D0-2.D0*ZNU)/(4.d0*pi*ZK*(1.D0-ZNU))*( eta(j)*(KS2+KS0-2.D0/XKSI**(2.d0))+&
																	   rd(j)*rdn*(4.D0/XKSI**(2.d0)-2.d0*KS2) )

			cqmF(k,j)=temp
			goto 100
!
!	F33
!
40			temp=-DRTSQ*rdn*KS1/(2.D0*pi)
			cqmF(k,j)=temp
			goto 100
!
100			continue
!
		enddo
	enddo
!
!
	return
	end subroutine
!
!*******************************************************************************************


!*******************************************************************************************
!
	subroutine wsatqstfsolG(SLAP,RTSQ,KS0,KS1,KS2,cqmG)
!
	implicit double precision (a-h,o-z)
!
	COMPLEX(8) SLAP,RTSQ,DRTSQ,KS0,KS1,KS2,XKSI,temp,cqmG(kbem,kbem)
!
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /e/ gammaw,gammas,gammaa,grav,atmp
	common /point/ dist(2),ra,rd(2),rdn,eta(2)
	common /UcoefR/ r1(2,2),dpr1(2,2,2),r2(2,2),dpr2(2,2,2),c(2,2)
	common /mate/zlamda,zmuy,znuu,znu,zalpha,zrho,zrhof,zm,zk,zr,zn
	common /ocqm/ EPScqm,taw,Icqm,Lcqm,NPcqm
!
	dimension dlt(2,2),abscqmG(kbem,kbem)
!
!
	do i=1,kbem
		do j=1,kbem
			cqmG(i,j)=0.d0
		enddo
	enddo
!
	pi=DACOS (-1.d0)
!
	dlt(1,1)=1.d0
	dlt(2,2)=1.d0
	dlt(1,2)=0.d0
	dlt(2,1)=0.d0
!
	DRTSQ=Cdsqrt( RTSQ )
	XKSI=DRTSQ*ra
!
	do k=1,3
		do j=1,3
			temp=DCMPLX (0.d0,0.d0)
			if (k.lt.3.and.j.eq.3) GOTO 10
			if (k.eq.3.and.j.lt.3) GOTO 30
			if (k.eq.3.and.j.eq.3) GOTO 40
!
!	Gij
!
			temp=( -(3.D0-4.D0*ZNUU)/(4.D0*ZMUY*(1.D0-ZNUU))*dlt(k,j)*dlog(ra)+&
				             rd(k)*rd(j)/(4.D0*ZMUY*(1.D0-ZNUU))+&
							 (ZNUU-ZNU)/(2.D0*ZMUY*(1.D0-ZNUU)*(1.D0-ZNU))*( dlt(k,j)*&
							    (XKSI**(-2.D0)-KS1/XKSI)+rd(k)*rd(j)*(KS2-2.D0*XKSI**(-2.D0)) ))/(2.D0*pi)
			cqmG(k,j)=temp
			goto 100
!
!	Gi3
!
10			temp=-(ZNUU-ZNU)/(SLAP*zalpha*(1.D0-2.D0*ZNU)*(1.D0-ZNUU))*DRTSQ*rd(k)*(KS1-1.D0/XKSI)/(2.D0*pi)
			cqmG(k,j)=temp
			goto 100
!
!	G3j
!
30			temp=(ZNUU-ZNU)/(zalpha*(1.D0-2.D0*ZNU)*(1.D0-ZNUU))*rd(j)*DRTSQ*(KS1-1.D0/XKSI)/(2.D0*pi)
			cqmG(k,j)=temp
			goto 100
!
!	G33
!
40			temp=-KS0/(2.D0*pi*ZK)
			cqmG(k,j)=temp
			goto 100
!
100			continue
!
		enddo
	enddo
!
!
	return
	end subroutine
!
!*******************************************************************************************




!*******************************************************************************************
!
	subroutine WFSOLGCOEF(SLAP,rab1,rab2,r331,r332,r441,r442,r341,r342,r431,r432,ra41,ra31,r4b1,&
					r3b1,rab3,rab4)
!
	implicit double precision (a-h,o-z)
!
	COMPLEX(8) SLAP,r331,r332,r441,r442,r341,r342,r431,r432,ra41,&
			   ra31,r4b1,r3b1,rab1,rab2,rab3,rab4
!
	common /e/ gammaw,gammas,gammaa,grav,atmp
	common /point/ dist(2),ra,rd(2),rdn,eta(2)
	common /umate/ d(3,3),dfs(3),xn,xkw,cww,cwg,xka,cgg,cgw,ev,sat,g1,rs,rw,raa,rmix,xm1,&
				   bt,et,xlambda,xmu,henry
!
!
	rab1=DCMPLX (0.d0,0.d0)
	rab2=DCMPLX (0.d0,0.d0)
	rab3=DCMPLX (0.d0,0.d0)
	rab4=DCMPLX (0.d0,0.d0)
	r331=DCMPLX (0.d0,0.d0)
	r332=DCMPLX (0.d0,0.d0)
	r441=DCMPLX (0.d0,0.d0)
	r442=DCMPLX (0.d0,0.d0)
	r341=DCMPLX (0.d0,0.d0)
	r342=DCMPLX (0.d0,0.d0)
	r442=DCMPLX (0.d0,0.d0)
	r431=DCMPLX (0.d0,0.d0)
	r432=DCMPLX (0.d0,0.d0)
	ra31=DCMPLX (0.d0,0.d0)
	ra41=DCMPLX (0.d0,0.d0)
	r4b1=DCMPLX (0.d0,0.d0)
	r3b1=DCMPLX (0.d0,0.d0)
!
	c1=xlambda+2.d0*xmu
	sata=1.d0-sat
	fs=dfs(1)
!
	pi=DACOS (-1.d0)
!
!	--- Gab coefficients ---
!
	rab1=(rw*fs+raa*(1.d0-fs))*SLAP**2.d0/(xlambda+xmu)
	if (xka.ne.0.d0) rab1=rab1-sata*(1.d0-fs)*SLAP/(xka*(xlambda+xmu))
	if (xkw.ne.0.d0) rab1=rab1-sat*fs*SLAP/(xkw*(xlambda+xmu))
	if (xka.ne.0.d0.and.xkw.ne.0.d0) rab1=rab1-(cgg*xkw+cww*xka)*SLAP/(xkw*xka)
!
	if (xka.ne.0.d0) rab2=rab2-rw*(cgg*fs-cwg*(1.d0-fs))*SLAP**3.d0/(xka*(xlambda+xmu))
	if (xkw.ne.0.d0) rab2=rab2-raa*(-cwg*fs+cww*(1.d0-fs))*SLAP**3.d0/(xkw*(xlambda+xmu))
	if (xka.ne.0.d0.and.xkw.ne.0.d0) rab2=rab2+( (-cwg**2.d0+cww*cgg)+&
										  sat*(-cwg*(1.d0-fs)+cgg*fs)/(xlambda+xmu)+&
										  sata*(-cwg*fs+cww*(1.d0-fs))/(xlambda+xmu) )*&
										  SLAP**2.d0/(xkw*xka)
!
	rab3=rmix*SLAP**2.d0/c1
	rab4=-(xlambda+xmu)/(rmix*SLAP**2.d0)
!
!	--- G33 coefficients ---
!
	r331=(raa*(1.d0-fs)+rmix)*SLAP**2.d0/c1
	if (xka.ne.0.d0) r331=r331-( sata*(1.d0-fs)/(c1*xka)+cgg/xka )*SLAP
!
	if (xka.ne.0.d0) r332=r332-rmix*cgg*SLAP**3.d0/(c1*xka)
!
!	--- G44 coefficients ---
!
	r441=( rw*fs/c1+rmix/c1 )*SLAP**2.d0
	if (xkw.ne.0.d0) r441=r441-( sat*fs/(c1*xkw)+cww/xkw )*SLAP
!
	if (xkw.ne.0.d0) r442=r442-rmix*cww*SLAP**3.d0/(c1*xkw)
!
!	--- G34 coefficients ---
!
	r341=-c1*cwg+(raa*xka*SLAP-sata)*fs
	r342=rmix*cwg*SLAP**2.d0
!
!	--- G43 coefficients ---
!
	r431=-c1*cwg+(rw*xkw*SLAP-sat)*(1.d0-fs)
	r432=rmix*cwg*SLAP**2.d0
!
!	--- Ga4 coefficients ---
!
	if (xka.ne.0.d0.and.xkw.ne.0.d0) ra41=(cwg*(sat-rw*xkw*SLAP)-cww*(sata-raa*xka*SLAP))*&
										   SLAP/(xkw*(sata-raa*xka*SLAP))
!
!	--- Ga3 coefficients ---
!
	if (xka.ne.0.d0.and.xkw.ne.0.d0) ra31=(cwg*(sata-raa*xka*SLAP)-cgg*(sat-rw*xkw*SLAP))*&
										   SLAP/(xka*(sat-rw*xkw*SLAP))
!
!	--- G4b coefficients ---
!
	if (xkw.ne.0.d0) r4b1=(-cwg*fs+cww*(1.d0-fs))*SLAP/(-xkw*(1-fs))
!
!	--- G3b coefficients ---
!
	if (xka.ne.0.d0) r3b1=(cwg*(1.d0-fs)-cgg*fs)*SLAP/(xka*fs)
!
	return
	end



!*******************************************************************************************
!
	subroutine wunsatdynfsolF(SLAP,RT,K0,K1,DK0,DK1,cqmG,cqmF)
!
	implicit double precision (a-h,o-z)
!
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /e/ gammaw,gammas,gammaa,grav,atmp
	common /point/ dist(2),ra,rd(2),rdn,eta(2)
	common /umate/ d(3,3),dfs(3),xn,xkw,cww,cwg,xka,cgg,cgw,ev,sat,g1,rs,rw,raa,rmix,xm1,&
				   bt,et,xlambda,xmu,henry
	common /ocqm/ EPScqm,taw,Icqm,Lcqm,NPcqm
	common /UcoefR/ r1(2,2),dpr1(2,2,2),r2(2,2),dpr2(2,2,2),c(2,2)
!
	COMPLEX(8) SLAP,RT(4),DRT(4),K0(4),K1(4),DK0(4),DK1(4),CsS(4),CsW(3),CsA(3),CwS(3),&
			   CaS(3),CCwW(3),CwA(3),CaW(3),CaA(3),cqmG(kbem,kbem),cqmF(kbem,kbem),&
			   temp,temp1,temp2
!
	dimension dlt(2,2),abscqmF(kbem,kbem)
!
!
	do i=1,kbem
		do j=1,kbem
			cqmF(i,j)=0.d0
		enddo
	enddo
!
	pi=DACOS (-1.d0)
	cc1=1.D0/(2.D0*pi*xmu)
!
	dlt(1,1)=1.d0
	dlt(2,2)=1.d0
	dlt(1,2)=0.d0
	dlt(2,1)=0.d0
!
	do i=1,4
		DRT(i)=Cdsqrt(RT(i))
	enddo
!
	call wfsolFcoef(SLAP,DRT,RT,CsS,CsW,CsA,CwS,CaS,CCwW,CwA,CaW,CaA)
!
	do i=1,4
		do j=1,4
			temp=DCMPLX (0.d0,0.d0)
			temp1=DCMPLX (0.d0,0.d0)
			temp2=DCMPLX (0.d0,0.d0)
			if (i.lt.3.and.j.eq.3) GOTO 10
			if (i.lt.3.and.j.eq.4) GOTO 20
			if (i.eq.3.and.j.lt.3) GOTO 30
			if (i.eq.3.and.j.eq.3) GOTO 40
			if (i.eq.3.and.j.eq.4) GOTO 50
			if (i.eq.4.and.j.lt.3) GOTO 60
			if (i.eq.4.and.j.eq.3) GOTO 70
			if (i.eq.4.and.j.eq.4) GOTO 80
!
!	Tab
!
			R3=rdn*( dlt(i,j)-4.d0*rd(I)*rd(J) )+rd(i)*eta(J)+rd(J)*eta(i)
			R4=rd(I)*rd(J)*rdn
!
			do k=1,4
				temp1=temp1-rd(I)*eta(J)*CsS(k)*RT(k)*DRT(k)*K1(k)
			enddo
			temp1=temp1-(rd(I)*eta(J)/(2.d0*pi*xmu))*DRT(1)*K1(1)
!
			do k=1,4
				temp2=temp2+2.d0*CsS(k)*( R3*RT(k)*K0(k)+(2.d0*R3/ra-R4*RT(k))*DRT(k)*K1(k) )
			enddo
			temp2=temp2-DRT(1)*K1(1)*(dlt(i,j)*rdn+rd(J)*eta(i))/(2.d0*pi*xmu)
!
			temp=temp2*xmu+(temp1*xlambda-SLAP*sat*cqmG(i,3)-SLAP*(1.d0-sat)*cqmG(i,4) )*eta(J)
			cqmF(i,j)=temp
			GOTO 100
!
!	Ta3
!
10			do m=1,3
				temp=CwS(m)*( (eta(i)-2.d0*rdn*rd(I))*DRT(m+1)*K1(m+1)/ra-rdn*rd(I)*RT(m+1)*K0(m+1) )
			enddo
			cqmF(i,j)=temp
			GOTO 100
!
!	Ta4
!
20			do m=1,3
				temp=CaS(m)*( (eta(i)-2.d0*rdn*rd(I))*DRT(m+1)*K1(m+1)/ra-rdn*rd(I)*RT(m+1)*K0(m+1) )
			enddo
			cqmF(i,j)=temp
			GOTO 100
!
!	T3b
!
30			do m=1,3
				temp1=temp1+2.D0*CsW(m)*( (eta(j)-2.d0*rdn*rd(j))*DRT(m+1)*K1(m+1)/ra-rdn*rd(j)*RT(m+1)*K0(m+1) )
				temp2=temp2-CsW(m)*RT(m+1)*K0(m+1)
			enddo
			temp=temp1*xmu+(temp2*xlambda-SLAP*sat*cqmG(3,3)-SLAP*(1.d0-sat)*cqmG(3,4))*eta(j)
			cqmF(i,j)=temp
			GOTO 100
!
!	T33
!
40			do m=1,3
				temp=temp-CCwW(m)*rdn*DRT(m+1)*K1(m+1)
			enddo
			cqmF(i,j)=temp
			GOTO 100
!
!	T34
!
50			do m=1,3
				temp=temp-CaW(m)*rdn*DRT(m+1)*K1(m+1)
			enddo
			cqmF(i,j)=temp
			GOTO 100
!
!	T4b
!
60			do m=1,3
				temp1=temp1+2.D0*CsA(m)*( (eta(j)-2.d0*rdn*rd(j))*DRT(m+1)*K1(m+1)/ra-rdn*rd(j)*RT(m+1)*K0(m+1) )
				temp2=temp2-CsA(m)*RT(m+1)*K0(m+1)
			enddo
			temp=temp1*xmu+(temp2*xlambda-SLAP*sat*cqmG(4,3)-SLAP*(1.d0-sat)*cqmG(4,4))*eta(j)
			cqmF(i,j)=temp
			GOTO 100
!
!	T43
!
70			do m=1,3
				temp=temp-CwA(m)*rdn*DRT(m+1)*K1(m+1)
			enddo
			cqmF(i,j)=temp
			GOTO 100
!
!	T44
!
80			do m=1,3
				temp=temp-CaA(m)*rdn*DRT(m+1)*K1(m+1)
			enddo
			cqmF(i,j)=temp
			GOTO 100
!
100		continue
		enddo
	enddo
!
	return
	end


!*******************************************************************************************
!
	subroutine wunsatdynfsolG(SLAP,RT,K0,K1,cqmG)
!
	implicit double precision (a-h,o-z)
!
	COMPLEX(8) SLAP,RT(4),DRT(4),K0(4),K1(4),r331,r332,r441,r442,r341,r342,r431,r432,ra41,&
			   ra31,r4b1,r3b1,rab1,rab2,rab3,rab4,cqmG(kbem,kbem)
!
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /e/ gammaw,gammas,gammaa,grav,atmp
	common /point/ dist(2),ra,rd(2),rdn,eta(2)
	common /umate/ d(3,3),dfs(3),xn,xkw,cww,cwg,xka,cgg,cgw,ev,sat,g1,rs,rw,raa,rmix,xm1,&
				   bt,et,xlambda,xmu,henry
	common /ocqm/ EPScqm,taw,Icqm,Lcqm,NPcqm
!
	dimension dlt(2,2),abscqmG(kbem,kbem)
!
!
	do i=1,kbem
		do j=1,kbem
			cqmG(i,j)=0.d0
		enddo
	enddo
!
	znu=xlambda/(2.D0*(xlambda+xmu))
!
	c1=xlambda+2.D0*xmu
	sata=1.D0-sat
	fs=dfs(1)
!
	pi=DACOS (-1.d0)
	dlt(1,1)=1.d0
	dlt(2,2)=1.d0
	dlt(1,2)=0.d0
	dlt(2,1)=0.d0
!
	call WFSOLGCOEF(SLAP,rab1,rab2,r331,r332,r441,r442,r341,r342,r431,r432,ra41,ra31,r4b1,&
					r3b1,rab3,rab4)
!
	do i=1,4
		DRT(i)=Cdsqrt( RT(i) )
	enddo
!
!	Gab
!
	do k=1,2
		do j=1,2
			R1=(2.D0*rd(k)*rd(j)-dlt(k,j))/ra
			R2=rd(k)*rd(j)
			do i=1,4
				i1=i
				i2=i+1
				i3=i+2
				i4=i+3
				if (i2.gt.4) i2=i2-4
				if (i3.gt.4) i3=i3-4
				if (i4.gt.4) i4=i4-4
				cqmG(k,j)=cqmG(k,j)+(RT(i1)**2.D0-rab1*RT(i1)+rab2)*&
									(R1*DRT(i1)*K1(i1)+R2*RT(i1)*K0(i1))/&
									((RT(i1)-RT(i2))*(RT(i1)-RT(i3))*(RT(i1)-RT(i4)))*&
									rab3*rab4/(2.D0*pi*xmu)
			enddo
			cqmG(k,j)=cqmG(k,j)+dlt(k,j)*K0(1)/(2.D0*pi*xmu)
		enddo
	enddo
!
	do i=2,4
		i1=i
		i2=i+1
		i3=i+2
		if (i2.gt.4) i2=i2-3
		if (i3.gt.4) i3=i3-3
!
!	G33
		if (xkw.ne.0.d0) cqmG(3,3)=cqmG(3,3)-((RT(i1)**2.D0-r331*RT(i1)+r332)/&
							((RT(i3)-RT(i1))*(RT(i2)-RT(i1))))*K0(i1)/(2.D0*pi*xkw)
!
!	G44
!
		if (xka.ne.0.d0) cqmG(4,4)=cqmG(4,4)-((RT(i1)**2.D0-r441*RT(i1)+r442)/&
							((RT(i3)-RT(i1))*(RT(i2)-RT(i1))))*K0(i1)/(2.D0*pi*xka)
!
!	G34
!
		if (xka.ne.0.d0.and.xkw.ne.0.d0) cqmG(3,4)=cqmG(3,4)-(r431*RT(i1)+r432)/&
											((RT(i3)-RT(i1))*(RT(i2)-RT(i1)))*&
											  K0(i1)*SLAP/(2.D0*c1*pi*xka*xkw)
!
!	G43
!
		if (xka.ne.0.d0.and.xkw.ne.0.d0) cqmG(4,3)=cqmG(4,3)-(r341*RT(i1)+r342)/&
											((RT(i3)-RT(i1))*(RT(i2)-RT(i1)))*&
											  K0(i1)*SLAP/(2.D0*c1*pi*xka*xkw)
!
!	Ga3 & Ga4 & G3b & G4b
!
		do k=1,2
			if (xkw.ne.0.d0) cqmG(k,3)=cqmG(k,3)-(-DRT(i1)*(RT(i1)-r3b1)/( (RT(i3)-RT(i1))*&
												(RT(i2)-RT(i1))))*K1(i1)*rd(k)*&
												fs/(2.d0*c1*pi*xkw)
			if (xka.ne.0.d0) cqmG(k,4)=cqmG(k,4)-(-DRT(i1)*(RT(i1)-r4b1)/((RT(i3)-RT(i1))*&
												(RT(i2)-RT(i1))))*K1(i1)*rd(k)*&
												(1.d0-fs)/(2.d0*c1*pi*xka)
			if (xkw.ne.0.d0) cqmG(3,k)=cqmG(3,k)+( -DRT(i1)*(RT(i1)-ra31)/( (RT(i3)-RT(i1))*&
												(RT(i2)-RT(i1))))*K1(i1)*rd(k)*&
												(rw*xkw*SLAP-sat)*SLAP/(2.d0*c1*pi*xkw)
			if (xka.ne.0.d0) cqmG(4,k)=cqmG(4,k)+(-DRT(i1)*(RT(i1)-ra41)/((RT(i3)-RT(i1))*&
												(RT(i2)-RT(i1))))*K1(i1)*rd(k)*&
												(raa*xka*SLAP-sata)*SLAP/(2.d0*c1*pi*xka)
		enddo
!
	enddo
!
	do i=1,4
		do j=1,4
			abscqmG(i,j)=Cdabs (cqmG(i,j))
		enddo
	enddo
!
	return
	end subroutine
!
!*******************************************************************************************
!*******************************************************************************************
!
	subroutine wunsatqstfsolf(SLAP,RTQ,KQ0,KQ1,DKQ0,DKQ1,cqmG,cqmF)
!
	implicit double precision (a-h,o-z)
!
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /e/ gammaw,gammas,gammaa,grav,atmp
	common /point/ dist(2),ra,rd(2),rdn,eta(2)
	common /umate/ d(3,3),dfs(3),xn,xkw,cww,cwg,xka,cgg,cgw,ev,sat,g1,rs,rw,raa,rmix,xm1,bt,et,xlambda,xmu,henry
	common /ocqm/ EPScqm,taw,Icqm,Lcqm,NPcqm
	common /UcoefR/ r1(2,2),dpr1(2,2,2),r2(2,2),dpr2(2,2,2),c(2,2)
!
	complex(8) SLAP,RTQ(2),DRTQ(2),KQ0(2),KQ1(2),DKQ0(2),DKQ1(2),cqmG(kbem,kbem),GAM11,GAM12, & 
			   GAM13,GAM21,GAM22,GAM31,GAM32,GAM33,GAM1,GAM2,GAM3,cqmF(kbem,kbem),dpsolG(kbem,kbem,2), &
			   temp,dila1,T12,T21,temp1,temp2,R8,R9,K11,K12,K13,K21,K22,K23,K31,K32,K41,K42,K51,K52,K61,K62, &
			   K71,K72,K73,K74,K75,K76,K77
!
	dimension dlt(2,2),abscqmF(kbem,kbem)
!
!
	do i=1,kbem
		do j=1,kbem
			cqmF(i,j)=0.d0
		enddo
	enddo
!
	pi=DACOS (-1.d0)
	zmu=3.d0*bt*et/(9.d0*bt-et)
	zlamda=bt-2.d0*zmu/3.d0
	znu=zlamda/(2.d0*(zlamda+zmu))
!
	dlt(1,1)=1.d0
	dlt(2,2)=1.d0
	dlt(1,2)=0.d0
	dlt(2,1)=0.d0
!
	do i=1,2
		DRTQ(i)=Cdsqrt(RTQ(i))
	enddo
!
!	call WFSOLGCOEFQST(SLAP,K11,K12,K13,K21,K22,K23,K31,K32,K41,K42,K51,K52,K61,K62,K71,K72,&
!					   K73,K74,K75,K76,K77)
!
!	call wdpsolqsG (SLAP,DRTQ,RTQ,KQ0,KQ1,DKQ0,DKQ1,GAM11,GAM12,GAM13,GAM21,GAM22,&
!					GAM31,GAM32,GAM33,GAM1,GAM2,GAM3,dpsolG)
!
	do i=1,4
		do j=1,4
			temp=DCMPLX (0.d0,0.d0)
			temp1=DCMPLX (0.d0,0.d0)
			temp2=DCMPLX (0.d0,0.d0)
			if (i.lt.3.and.j.eq.3) GOTO 10
			if (i.lt.3.and.j.eq.4) GOTO 20
			if (i.eq.3.and.j.lt.3) GOTO 30
			if (i.eq.3.and.j.eq.3) GOTO 40
			if (i.eq.3.and.j.eq.4) GOTO 50
			if (i.eq.4.and.j.lt.3) GOTO 60
			if (i.eq.4.and.j.eq.3) GOTO 70
			if (i.eq.4.and.j.eq.4) GOTO 80
!
!	Tab
!
			CK11=1.D0/(2.D0*pi*zmu)
			CK12=( xn*g1*(zlamda+2.D0*zmu)*(xka+xkw)+xkw*(1.d0-sat)*(-1.d0+dfs(1))-&
				   xka*sat*dfs(1) )/&
				 (2.D0*pi*zmu*(zlamda+2.D0*zmu)*xkw*xka)
			CK13=-xn*g1/(2.D0*pi*zmu*(zlamda+2.D0*zmu)*xkw*xka)
			CK21=-(zlamda+zmu)/(2.D0*pi*zmu*(zlamda+2.D0*zmu))
			CK22=-( xn*g1*(zlamda+zmu)*(xka+xkw)+xkw*(1.d0-sat)*(-1.d0+dfs(1))-&
					xka*sat*dfs(1) )/&
				  (2.D0*pi*zmu*(zlamda+2.D0*zmu)*xkw*xka)
			CK23=-CK13
!
			temp1=rd(I)*eta(J)*&
				( (RTQ(1)*RTQ(1)*(CK11+CK21)+SLAP*RTQ(1)*(CK12+CK22)+SLAP*SLAP*(CK13+CK23))*KQ1(1)/&
				   DRTQ(1)-&
				  (RTQ(2)*RTQ(2)*(CK11+CK21)+SLAP*RTQ(2)*(CK12+CK22)+SLAP*SLAP*(CK13+CK23))*KQ1(2)/&
				   DRTQ(2) )/&
				(RTQ(2)-RTQ(1))
!
!			temp1=temp1*zlamda+rd(I)*eta(J)*SLAP*&
!								( (sat*dfs(1)*gammaw*xka+(1.d0-sat)*(1.d0-dfs(1))*gammaa*xkw)*&
!										(KQ1(2)*DRTQ(2)-KQ1(1)*DRTQ(1))+&
!								  (g1*gammaw*gammaa*xn*SLAP)*(KQ1(2)/DRTQ(2)-KQ1(1)/DRTQ(1)) )/&
!							    (2.D0*pi*(zlamda+2.d0*zmu)*xkw*xka*(RTQ(1)-RTQ(2)))
!
			temp1=temp1*zlamda+rd(J)*eta(I)*SLAP*&
								( (sat*dfs(1)*xka+(1.d0-sat)*(1.d0-dfs(1))*xkw)*&
										(KQ1(2)*DRTQ(2)-KQ1(1)*DRTQ(1))+&
								  (g1*xn*SLAP)*(KQ1(2)/DRTQ(2)-KQ1(1)/DRTQ(1)) )/&
							    (2.D0*pi*(zlamda+2.d0*zmu)*xkw*xka*(RTQ(1)-RTQ(2)))
!
			R7=rdn*( dlt(i,j)-4.d0*rd(I)*rd(J) )+rd(i)*eta(J)+rd(J)*eta(i)
			R8=R7*DRTQ(1)*( DRTQ(1)*KQ0(1)+2.D0*KQ1(1)/ra )/ra-&
			   rd(i)*rd(j)*rdn*RTQ(1)*DRTQ(1)*KQ1(1)
			R9=R7*DRTQ(2)*( DRTQ(2)*KQ0(2)+2.D0*KQ1(2)/ra )/ra-&
			   rd(i)*rd(j)*rdn*RTQ(2)*DRTQ(2)*KQ1(2)
!
			temp2=(rdn*dlt(i,j)+rd(I)*eta(J))*&
				  ( (RTQ(1)*RTQ(1)*CK11+SLAP*RTQ(1)*CK12+SLAP*SLAP*CK13)*KQ1(1)/DRTQ(1)-&
				    (RTQ(2)*RTQ(2)*CK11+SLAP*RTQ(2)*CK12+SLAP*SLAP*CK13)*KQ1(2)/DRTQ(2)  )/&
				  (RTQ(2)-RTQ(1))+&
				  2.d0*(RTQ(1)*RTQ(1)*CK21+SLAP*RTQ(1)*CK22+SLAP*SLAP*CK23)*R8/&
				       (RTQ(1)*RTQ(1)*(RTQ(1)-RTQ(2)))-&
				  2.d0*(RTQ(2)*RTQ(2)*CK21+SLAP*RTQ(2)*CK22+SLAP*SLAP*CK23)*R9/&
				       (RTQ(2)*RTQ(2)*(RTQ(1)-RTQ(2)))
!
			temp=temp2*zmu+temp1
!
			cqmF(i,j)=temp
			GOTO 100
!
!	Ta3
!
10			temp=-(1.d0/(2.d0*pi*(zlamda+2.D0*zmu)*xka*(RTQ(1)-RTQ(2))))*&
				 ( (eta(i)-2.d0*rdn*rd(i))*( dfs(1)*xka*&
										     (DRTQ(2)*KQ1(2)/ra-DRTQ(1)*KQ1(1)/ra)+&
											  g1*xn*SLAP*&
											  (KQ1(2)/(ra*DRTQ(2))-KQ1(1)/(ra*DRTQ(1))) )-&
				   rd(i)*rdn*( dfs(1)*xka*(RTQ(2)*KQ0(2)-RTQ(1)*KQ0(1))+&
							   g1*xn*SLAP*(KQ0(2)-KQ0(1)) )	)
!
			cqmF(i,j)=temp
			GOTO 100
!
!	Ta4
!
20			temp=-(1.d0/(2.d0*pi*(zlamda+2.D0*zmu)*xkw*(RTQ(1)-RTQ(2))))*&
				 ( (eta(i)-2.d0*rdn*rd(i))*( (1.d0-dfs(1))*xkw*&
										     (DRTQ(2)*KQ1(2)/ra-DRTQ(1)*KQ1(1)/ra)+&
											  g1*xn*SLAP*&
											  (KQ1(2)/(ra*DRTQ(2))-KQ1(1)/(ra*DRTQ(1))) )-&
				   rd(i)*rdn*( (1.d0-dfs(1))*xkw*(RTQ(2)*KQ0(2)-RTQ(1)*KQ0(1))+&
								g1*xn*SLAP*(KQ0(2)-KQ0(1)) )	)
!
			cqmF(i,j)=temp
			GOTO 100
!
!	T3b
!
30			temp1=eta(j)*(xka*sat*SLAP*(RTQ(2)*KQ0(2)-RTQ(1)*KQ0(1))+&
						  g1*xn*SLAP*SLAP*(KQ0(2)-KQ0(1)))/&
						 (2.d0*pi*(zlamda+2.D0*zmu)*xkw*xka*(RTQ(1)-RTQ(2)))
!
			temp1=temp1*zlamda+eta(j)*SLAP*( sat*xka*(RTQ(2)*KQ0(2)-RTQ(1)*KQ0(1))+&
											 g1*xn*SLAP*(KQ0(2)-KQ0(1)) ) /&
							   (2.d0*pi*xkw*xka*(RTQ(1)-RTQ(2)))
!
			temp2=(1.d0/( pi*(zlamda+2.D0*zmu)*xka*xkw*(RTQ(2)-RTQ(1))))*&
				  ( ((eta(j)-2.d0*rdn*rd(j))/ra)*&
				    (xka*sat*SLAP*(DRTQ(2)*KQ1(2)-DRTQ(1)*KQ1(1))+&
					 g1*xn*SLAP*SLAP*(KQ1(2)/DRTQ(2)-KQ1(1)/DRTQ(1)))-&
					 rdn*rd(j)*(xka*sat*SLAP*(RTQ(2)*KQ0(2)-RTQ(1)*KQ0(1))+&
					            g1*xn*SLAP*SLAP*(KQ0(2)-KQ0(1)))    )
!
			temp=temp2*zmu+temp1
!
			cqmF(i,j)=temp
			GOTO 100
!
!	T33
!
40			temp=rdn*(1.d0/(2.d0*pi*(zlamda+2.D0*zmu)*xka*(RTQ(2)-RTQ(1))))*&
				 ( xka*(zlamda+2.D0*zmu)*(KQ1(2)*RTQ(2)*DRTQ(2)-KQ1(1)*RTQ(1)*DRTQ(1))+&
				   ( (zlamda+2.D0*zmu)*g1*xn-(1.d0-sat)*(1.d0-dfs(1)) )*SLAP*&
				          (KQ1(2)*DRTQ(2)-KQ1(1)*DRTQ(1)) )
!
			cqmF(i,j)=temp
			GOTO 100
!
!	T34
!
50			temp=rdn*( (sat*(1.d0-dfs(1))+g1*xn*(zlamda+2.D0*zmu))/&
						(2.d0*pi*(zlamda+2.D0*zmu)*xkw) )*SLAP*&
					 (KQ1(2)*DRTQ(2)-KQ1(1)*DRTQ(1))/(RTQ(2)-RTQ(1))
!
			cqmF(i,j)=temp
			GOTO 100
!
!	T4b
!
60			temp1=eta(j)*(xkw*(1.d0-sat)*SLAP*(RTQ(2)*KQ0(2)-RTQ(1)*KQ0(1))+&
						  g1*xn*SLAP*SLAP*(KQ0(2)-KQ0(1)))/&
						 (2.d0*pi*(zlamda+2.D0*zmu)*xkw*xka*(RTQ(1)-RTQ(2)))
!
			temp1=temp1*zlamda+eta(j)*SLAP*( (1.d0-sat)*xkw*(RTQ(2)*KQ0(2)-RTQ(1)*KQ0(1))+&
											 g1*xn*SLAP*(KQ0(2)-KQ0(1)) ) /&
							   (2.d0*pi*xkw*xka*(RTQ(1)-RTQ(2)))
!
			temp2=(1.d0/( pi*(zlamda+2.D0*zmu)*xka*xkw*(RTQ(2)-RTQ(1))))*&
				  ( ((eta(j)-2.d0*rdn*rd(j))/ra)*&
				    (xkw*(1.d0-sat)*SLAP*(DRTQ(2)*KQ1(2)-DRTQ(1)*KQ1(1))+&
					 g1*xn*SLAP*SLAP*(KQ1(2)/DRTQ(2)-KQ1(1)/DRTQ(1)))-&
					 rdn*rd(j)*(xkw*(1.d0-sat)*SLAP*(RTQ(2)*KQ0(2)-RTQ(1)*KQ0(1))+&
					            g1*xn*SLAP*SLAP*(KQ0(2)-KQ0(1)))    )
!
			temp=temp2*zmu+temp1
!
			cqmF(i,j)=temp
			GOTO 100
!
!	T43
!
70			temp=rdn*(dfs(1)*(1.d0-sat)+(zlamda+2.D0*zmu)*g1*xn)*SLAP*&
                 (KQ1(2)*DRTQ(2)-KQ1(1)*DRTQ(1))/(2.d0*pi*(zlamda+2.D0*zmu)*xka*(RTQ(2)-RTQ(1)))
!
			cqmF(i,j)=temp
			GOTO 100
!
!	T44
!
80			temp=rdn*(xkw*(zlamda+2.D0*zmu)*(KQ1(2)*RTQ(2)*DRTQ(2)-KQ1(1)*RTQ(1)*DRTQ(1))-&
                      (sat*dfs(1)-(zlamda+2.D0*zmu)*g1*xn)*SLAP*&
					         (KQ1(2)*DRTQ(2)-KQ1(1)*DRTQ(1)))/&
				 (2.d0*pi*(zlamda+2.D0*zmu)*xkw*(RTQ(2)-RTQ(1)))
!
			cqmF(i,j)=temp
			GOTO 100
!
100		continue
		enddo
	enddo
!
	return
	end

!*******************************************************************************************
!
	subroutine wunsatqstfsolF1(SLAP,RTQ,KQ0,KQ1,KQ2,DKQ0,DKQ1,cqmG,cqmF)
!
	implicit double precision (a-h,o-z)
!
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /e/ gammaw,gammas,gammaa,grav,atmp
	common /point/ dist(2),ra,rd(2),rdn,eta(2)
	common /umate/ d(3,3),dfs(3),xn,xkw,cww,cwg,xka,cgg,cgw,ev,sat,g1,rs,rw,raa,rmix,xm1,&
				   bt,et,xlambda,xmu,henry
	common /ocqm/ EPScqm,taw,Icqm,Lcqm,NPcqm
	common /UcoefR/ r1(2,2),dpr1(2,2,2),r2(2,2),dpr2(2,2,2),c(2,2)
!
	COMPLEX(8) SLAP,RTQ(2),DRTQ(2),KQ0(2),KQ1(2),KQ2(2),DKQ0(2),DKQ1(2),cqmG(kbem,kbem),&
			   cqmF(kbem,kbem),dpsolG(kbem,kbem,2),temp,temp1,temp2,temp3,temp4,&
			   CK11,CK22,CK33,CK44
!
	dimension dlt(2,2),abscqmF(kbem,kbem)
!
!
	do i=1,kbem
		do j=1,kbem
			cqmF(i,j)=0.d0
		enddo
	enddo
!
	pi=DACOS (-1.d0)
	dlt(1,1)=1.d0
	dlt(2,2)=1.d0
	dlt(1,2)=0.d0
	dlt(2,1)=0.d0
!
	do i=1,2
		DRTQ(i)=Cdsqrt(RTQ(i))
	enddo
!
!
	do i=1,4
		do j=1,4
			temp=DCMPLX (0.d0,0.d0)
			temp1=DCMPLX (0.d0,0.d0)
			temp2=DCMPLX (0.d0,0.d0)
			if (i.lt.3.and.j.eq.3) GOTO 10
			if (i.lt.3.and.j.eq.4) GOTO 20
			if (i.eq.3.and.j.lt.3) GOTO 30
			if (i.eq.3.and.j.eq.3) GOTO 40
			if (i.eq.3.and.j.eq.4) GOTO 50
			if (i.eq.4.and.j.lt.3) GOTO 60
			if (i.eq.4.and.j.eq.3) GOTO 70
			if (i.eq.4.and.j.eq.4) GOTO 80
!
!	Tab
!
			CK11=xn*g1*SLAP*(xka+xkw)/(xka*xkw)
			CK22=CK11/(2.d0*pi*(xlambda+2.D0*xmu))
!
			R7=( rdn*( dlt(i,j)-4.d0*rd(I)*rd(J) )+rd(i)*eta(J)+rd(J)*eta(i) )/ra
			R8=rd(I)*(eta(J)-2.d0*rd(J)*rdn)/ra
!
			temp1=R8/(2.d0*pi*xmu)+(1.d0/(pi*(xlambda+2.D0*xmu)*(RTQ(1)-RTQ(2))))*&
			      ( R7*(KQ2(1)*RTQ(1)-KQ2(2)*RTQ(2))-rd(I)*rd(J)*rdn*(KQ1(1)*DRTQ(1)*RTQ(1)-KQ1(2)*DRTQ(2)*RTQ(2))+&
				    CK11*(R7*(KQ2(1)-KQ2(2))-rd(I)*rd(J)*rdn*(KQ1(1)*DRTQ(1)-KQ1(2)*DRTQ(2))) )+&
				  4.d0*CK22*R7/(RTQ(1)*RTQ(2)*ra**2.d0)
!
			temp2=-( rd(I)*eta(J)/(2.d0*pi*(xlambda+2.D0*xmu)*(RTQ(1)-RTQ(2))) )*&
			       ( (DRTQ(1)*RTQ(1)*KQ1(1)-DRTQ(2)*RTQ(2)*KQ1(2))+&
					  CK11*(DRTQ(1)*KQ1(1)-DRTQ(2)*KQ1(2))  )
!
			temp3=(SLAP*rw*sat*(-cqmG(i,3))+SLAP*raa*(1.d0-sat)*(-cqmG(i,4)))*eta(J)

!
			temp=xlambda*temp2+xmu*temp1+temp3
!
			cqmF(i,j)=temp
			GOTO 100
!
!	Ta3
!
10			temp=(xkw/(2.d0*pi*SLAP*(RTQ(1)-RTQ(2))))*&
			      ( ((eta(i)-2.d0*rdn*rd(i))/ra)*( (RTQ(2)+dfs(1)*SLAP/((xlambda+2.D0*xmu)*xkw))*DRTQ(1)*KQ1(1)+&
				                                   (RTQ(1)-RTQ(2))/ra-&
												   (RTQ(1)+dfs(1)*SLAP/((xlambda+2.D0*xmu)*xkw))*DRTQ(2)*KQ1(2) )-&
				     rd(i)*rdn*( (RTQ(2)+dfs(1)*SLAP/((xlambda+2.D0*xmu)*xkw))*RTQ(1)*KQ0(1)-&
					             (RTQ(1)+dfs(1)*SLAP/((xlambda+2.D0*xmu)*xkw))*RTQ(2)*KQ0(2)  )   )
!
			cqmF(i,j)=temp
			GOTO 100
!
!	Ta4
!
20			temp=(xka/(2.d0*pi*SLAP*(RTQ(1)-RTQ(2))))*&
			      ( ((eta(i)-2.d0*rdn*rd(i))/ra)*( (RTQ(2)+(1.d0-dfs(1))*SLAP/((xlambda+2.D0*xmu)*xka))*DRTQ(1)*KQ1(1)+&
				                                   (RTQ(1)-RTQ(2))/ra-&
												   (RTQ(1)+(1.d0-dfs(1))*SLAP/((xlambda+2.D0*xmu)*xka))*DRTQ(2)*KQ1(2) )-&
				     rd(i)*rdn*( (RTQ(2)+(1.d0-dfs(1))*SLAP/((xlambda+2.D0*xmu)*xka))*RTQ(1)*KQ0(1)-&
					             (RTQ(1)+(1.d0-dfs(1))*SLAP/((xlambda+2.D0*xmu)*xka))*RTQ(2)*KQ0(2)  )   )
!
			cqmF(i,j)=temp
			GOTO 100
!
!	T3b
!
30			temp1=-( eta(J)/(2.d0*pi*(RTQ(1)-RTQ(2))) )*( (RTQ(2)+sat*SLAP/((xlambda+2.D0*xmu)*xkw))*RTQ(1)*KQ0(1)-&
                                                          (RTQ(1)+sat*SLAP/((xlambda+2.D0*xmu)*xkw))*RTQ(2)*KQ0(2)  )
			temp2=(1.d0/(pi*(RTQ(1)-RTQ(2))))*&
			      ( ((eta(J)-2.d0*rdn*rd(J))/ra)*( (RTQ(2)+sat*SLAP/((xlambda+2.D0*xmu)*xkw))*DRTQ(1)*KQ1(1)+&
				                                   (RTQ(1)-RTQ(2))/ra-&
												   (RTQ(1)+sat*SLAP/((xlambda+2.D0*xmu)*xkw))*DRTQ(2)*KQ1(2) )-&
				     rd(J)*rdn*( (RTQ(2)+sat*SLAP/((xlambda+2.D0*xmu)*xkw))*RTQ(1)*KQ0(1)-&
					             (RTQ(1)+sat*SLAP/((xlambda+2.D0*xmu)*xkw))*RTQ(2)*KQ0(2)  )   )

			temp=xlambda*temp1+xmu*temp2+(SLAP*rw*sat*(-cqmG(3,3))+SLAP*raa*(1.d0-sat)*(-cqmG(3,4)))
!
			cqmF(i,j)=temp
			GOTO 100
!
!	T33
!
40			temp=(rdn/(2.d0*pi*(RTQ(1)-RTQ(2))))*&
			      ( (RTQ(1)+(1.d0-dfs(1))*SLAP*(1.d0-sat)/((xlambda+2.D0*xmu)*xka)+g1*xn*SLAP/xka)*DRTQ(1)*KQ1(1)-&
					(RTQ(2)+(1.d0-dfs(1))*SLAP*(1.d0-sat)/((xlambda+2.D0*xmu)*xka)+g1*xn*SLAP/xka)*DRTQ(2)*KQ1(2) )
!
			cqmF(i,j)=temp
			GOTO 100
!
!	T34
!
50			temp=-(rdn/(2.d0*pi*xkw*(RTQ(2)-RTQ(1))))*&
				  ( -g1*xn+sat*(1.d0-dfs(1))/(xlambda+2.D0*xmu) )*SLAP*&
				          (KQ1(2)*DRTQ(2)-KQ1(1)*DRTQ(1))
!
			cqmF(i,j)=temp
			GOTO 100
!
!	T4b
!
60			temp1=-( eta(J)/(2.d0*pi*(RTQ(1)-RTQ(2))) )*( (RTQ(2)+(1.d0-sat)*SLAP/((xlambda+2.D0*xmu)*xka))*RTQ(1)*KQ0(1)-&
                                                          (RTQ(1)+(1.d0-sat)*SLAP/((xlambda+2.D0*xmu)*xka))*RTQ(2)*KQ0(2)  )
			temp2=(1.d0/(pi*(RTQ(1)-RTQ(2))))*&
			      ( ((eta(J)-2.d0*rdn*rd(J))/ra)*( (RTQ(2)+(1.d0-sat)*SLAP/((xlambda+2.D0*xmu)*xka))*DRTQ(1)*KQ1(1)+&
				                                   (RTQ(1)-RTQ(2))/ra-&
												   (RTQ(1)+(1.d0-sat)*SLAP/((xlambda+2.D0*xmu)*xka))*DRTQ(2)*KQ1(2) )-&
				     rd(J)*rdn*( (RTQ(2)+(1.d0-sat)*SLAP/((xlambda+2.D0*xmu)*xka))*RTQ(1)*KQ0(1)-&
					             (RTQ(1)+(1.d0-sat)*SLAP/((xlambda+2.D0*xmu)*xka))*RTQ(2)*KQ0(2)  )   )

			temp=xlambda*temp1+xmu*temp2+(SLAP*rw*sat*(-cqmG(4,3))+SLAP*raa*(1.d0-sat)*(-cqmG(4,4)))
!
			cqmF(i,j)=temp
			GOTO 100
!
!	T43
!
70			temp=-(rdn/(2.d0*pi*xka*(RTQ(1)-RTQ(2))))*&
				  ( -g1*xn+(1.d0-sat)*dfs(1)/(xlambda+2.D0*xmu) )*SLAP*&
				          (KQ1(1)*DRTQ(1)-KQ1(2)*DRTQ(2))
!
			cqmF(i,j)=temp
			GOTO 100
!
!	T44
!
80			temp=(rdn/(2.d0*pi*(RTQ(2)-RTQ(1))))*&
				 ( (KQ1(2)*RTQ(2)*DRTQ(2)-KQ1(1)*RTQ(1)*DRTQ(1))+&
				   ( (g1*xn/xkw)+sat*dfs(1)/((xlambda+2.D0*xmu)*xkw) )*SLAP*&
				          (KQ1(2)*DRTQ(2)-KQ1(1)*DRTQ(1)) )
!
			cqmF(i,j)=temp
			GOTO 100
!
100		continue
		enddo
	enddo
!
	return
	end


!*******************************************************************************************
!
	subroutine wunsatqstfsolF2(SLAP,RTQ,KQ0,KQ1,DKQ0,DKQ1,cqmG,cqmF)
!
	implicit double precision (a-h,o-z)
!
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /e/ gammaw,gammas,gammaa,grav,atmp
	common /point/ dist(2),ra,rd(2),rdn,eta(2)
	common /umate/ d(3,3),dfs(3),xn,xkw,cww,cwg,xka,cgg,cgw,ev,sat,g1,rs,rw,raa,rmix,xm1,&
				   bt,et,xlambda,xmu,henry
	common /ocqm/ EPScqm,taw,Icqm,Lcqm,NPcqm
	common /UcoefR/ r1(2,2),dpr1(2,2,2),r2(2,2),dpr2(2,2,2),c(2,2)
!
	COMPLEX(8) SLAP,RTQ(2),DRTQ(2),KQ0(2),KQ1(2),DKQ0(2),DKQ1(2),cqmG(kbem,kbem),&
			   cqmF(kbem,kbem),dpsolG(kbem,kbem,2),temp,temp1,temp2,temp3,temp4,&
			   CK11,CK22,CK33,CK44,R8,R9
!
	dimension dlt(2,2),abscqmF(kbem,kbem)
!
!
	do i=1,kbem
		do j=1,kbem
			cqmF(i,j)=0.d0
		enddo
	enddo
!
	pi=DACOS (-1.d0)
	dlt(1,1)=1.d0
	dlt(2,2)=1.d0
	dlt(1,2)=0.d0
	dlt(2,1)=0.d0
!
	do i=1,2
		DRTQ(i)=Cdsqrt(RTQ(i))
	enddo
!
!
	do i=1,4
		do j=1,4
			temp=DCMPLX (0.d0,0.d0)
			temp1=DCMPLX (0.d0,0.d0)
			temp2=DCMPLX (0.d0,0.d0)
			if (i.lt.3.and.j.eq.3) GOTO 10
			if (i.lt.3.and.j.eq.4) GOTO 20
			if (i.eq.3.and.j.lt.3) GOTO 30
			if (i.eq.3.and.j.eq.3) GOTO 40
			if (i.eq.3.and.j.eq.4) GOTO 50
			if (i.eq.4.and.j.lt.3) GOTO 60
			if (i.eq.4.and.j.eq.3) GOTO 70
			if (i.eq.4.and.j.eq.4) GOTO 80
!
!	Tab
!
			CK11=( xn*g1*(xka+xkw)/(xka*xkw) )*SLAP
			CK44=xn*g1*SLAP*SLAP/(xkw*xka)
!
			R6=rdn*( dlt(i,j)-4.d0*rd(I)*rd(J) )+rd(i)*eta(J)+rd(J)*eta(i)
			R7=rd(I)*rd(J)*rdn
			R8=R6*DRTQ(2)*(KQ0(2)*DRTQ(2)+2.D0*KQ1(2)/ra)/ra-RTQ(2)*DRTQ(2)*rd(I)*rd(J)*rdn*KQ1(2)
			R9=R6*DRTQ(1)*(KQ0(1)*DRTQ(1)+2.D0*KQ1(1)/ra)/ra-RTQ(1)*DRTQ(1)*rd(I)*rd(J)*rdn*KQ1(1)

!
			temp1=( -2.d0*R6/(pi*ra**3.d0) )+&
			      ( xmu/(pi*(xlambda+2.D0*xmu)*(RTQ(1)-RTQ(2))) )*&
				  ( (RTQ(1)+CK11)*R9-(RTQ(2)+CK11)*R8    )
!
			temp2=( xlambda*rd(I)*eta(J)/(2.d0*pi*(xlambda+2.D0*xmu)*(RTQ(2)-RTQ(1))) )*&
				   ( (RTQ(1)+CK11)*KQ1(1)*RTQ(1)*DRTQ(1)-(RTQ(2)+CK11)*KQ1(2)*RTQ(2)*DRTQ(2) )+&
!
				  ( rd(I)*eta(J)/(2.d0*pi*(xlambda+2.D0*xmu)*(RTQ(1)-RTQ(2))) )*&
				  ( CK44*(KQ1(1)*DRTQ(1)-KQ1(2)*DRTQ(2))+&
				    (dfs(1)*sat*SLAP/xkw)*(KQ1(1)*DRTQ(1)*RTQ(1)-KQ1(2)*DRTQ(2)*RTQ(2))+&
					((1.d0-dfs(1))*(1.d0-sat)*SLAP/xka)*(KQ1(1)*DRTQ(1)*RTQ(1)-KQ1(2)*DRTQ(2)*RTQ(2)) )
!
			temp=temp1+temp2
!
			cqmF(i,j)=temp
			GOTO 100
!
!	Ta3
!
10			temp=(1.d0/(2.d0*pi*(xlambda+2.D0*xmu)*(RTQ(1)-RTQ(2))))*&
				 ( ((eta(i)-2.d0*rdn*rd(i))/ra)*&
						( dfs(1)*(DRTQ(1)*RTQ(1)*KQ1(1)-DRTQ(2)*RTQ(2)*KQ1(2))+&
						  (g1*xn*SLAP/xka)*(KQ1(1)*DRTQ(1)-KQ1(2)*DRTQ(2)) )-&
				   rd(i)*rdn*( dfs(1)*(RTQ(1)*RTQ(1)*KQ0(1)-RTQ(2)*RTQ(2)*KQ0(2))+&
							   (g1*xn*SLAP/xka)*(RTQ(1)*KQ0(1)-RTQ(2)*KQ0(2)) )	)
!
			cqmF(i,j)=temp
			GOTO 100
!
!	Ta4
!
20			temp=(1.d0/(2.d0*pi*(xlambda+2.D0*xmu)*(RTQ(1)-RTQ(2))))*&
				 ( ((eta(i)-2.d0*rdn*rd(i))/ra)*&
						( (1.d0-dfs(1))*(DRTQ(1)*RTQ(1)*KQ1(1)-DRTQ(2)*RTQ(2)*KQ1(2))+&
						  (g1*xn*SLAP/xkw)*(KQ1(1)*DRTQ(1)-KQ1(2)*DRTQ(2)) )-&
				   rd(i)*rdn*( (1.d0-dfs(1))*(RTQ(1)*RTQ(1)*KQ0(1)-RTQ(2)*RTQ(2)*KQ0(2))+&
							   (g1*xn*SLAP/xkw)*(RTQ(1)*KQ0(1)-RTQ(2)*KQ0(2)) )	)
!
			cqmF(i,j)=temp
			GOTO 100
!
!	T3b
!
30			temp=(1.d0/(pi*(xlambda+2.D0*xmu)*xkw*(RTQ(1)-RTQ(2))))*&
				 (	-(xlambda+xmu)*eta(j)*( (g1*xn*SLAP*SLAP/xka)*(KQ0(1)*RTQ(1)-KQ0(2)*RTQ(2))+&
											 sat*SLAP*(RTQ(1)*RTQ(1)*KQ0(1)-RTQ(2)*RTQ(2)*KQ0(2)) )+&
					xmu*( ((eta(j)-2.d0*rdn*rd(j))/ra)*&
								( (g1*xn*SLAP*SLAP/xka)*(KQ1(1)*DRTQ(1)-KQ1(2)*DRTQ(2))+&
								   sat*SLAP*(DRTQ(1)*RTQ(1)*KQ1(1)-DRTQ(2)*RTQ(2)*KQ1(2)) )-&
						    rd(j)*rdn*( (g1*xn*SLAP*SLAP/xka)*(RTQ(1)*KQ0(1)-RTQ(2)*KQ0(2))+&
						                 sat*SLAP*(RTQ(1)*RTQ(1)*KQ0(1)-RTQ(2)*RTQ(2)*KQ0(2)) )  )  )
!
			cqmF(i,j)=temp
			GOTO 100
!
!	T33
!
40			temp=(rdn/(2.d0*pi*(RTQ(1)-RTQ(2))))*&
				 ( KQ1(1)*RTQ(1)*DRTQ(1)*( RTQ(1)+(g1*xn/xka+(1.d0-sat)*(1.d0-dfs(1))/&
				                                    ((xlambda+2.D0*xmu)*xka))*SLAP )-&
				   KQ1(2)*RTQ(2)*DRTQ(2)*( RTQ(2)+(g1*xn/xka+(1.d0-sat)*(1.d0-dfs(1))/&
				                                    ((xlambda+2.D0*xmu)*xka))*SLAP )  )
!
			cqmF(i,j)=temp
			GOTO 100
!
!	T34
!
50			temp=-(rdn/(2.d0*pi*xkw*(RTQ(2)-RTQ(1))))*&
				  ( -g1*xn+sat*(1.d0-dfs(1))/(xlambda+2.D0*xmu) )*SLAP*&
				          (KQ1(2)*DRTQ(2)*RTQ(2)-KQ1(1)*DRTQ(1)*RTQ(1))
!
			cqmF(i,j)=temp
			GOTO 100
!
!	T4b
!
60			temp=(1.d0/(pi*(xlambda+2.D0*xmu)*xka*(RTQ(1)-RTQ(2))))*&
				 (	-(xlambda+xmu)*eta(j)*( (g1*xn*SLAP*SLAP/xkw)*(KQ0(1)*RTQ(1)-KQ0(2)*RTQ(2))+&
											(1.d0-sat)*SLAP*(RTQ(1)*RTQ(1)*KQ0(1)-RTQ(2)*RTQ(2)*KQ0(2)) )+&
					xmu*( ((eta(j)-2.d0*rdn*rd(J))/ra)*&
								( (g1*xn*SLAP*SLAP/xkw)*(KQ1(1)*DRTQ(1)-KQ1(2)*DRTQ(2))+&
								  (1.d0-sat)*SLAP*(DRTQ(1)*RTQ(1)*KQ1(1)-DRTQ(2)*RTQ(2)*KQ1(2)) )-&
						   rd(j)*rdn*( (g1*xn*SLAP*SLAP/xkw)*(RTQ(1)*KQ0(1)-RTQ(2)*KQ0(2))+&
						               (1.d0-sat)*SLAP*(RTQ(1)*RTQ(1)*KQ0(1)-RTQ(2)*RTQ(2)*KQ0(2)) )  )  )
!
			cqmF(i,j)=temp
			GOTO 100
!
!	T43
!
70			temp=-(rdn/(2.d0*pi*xka*(RTQ(1)-RTQ(2))))*&
				  ( -g1*xn+(1.d0-sat)*dfs(1)/(xlambda+2.D0*xmu) )*SLAP*&
				          (KQ1(1)*DRTQ(1)*RTQ(1)-KQ1(2)*DRTQ(2)*RTQ(2))
!
			cqmF(i,j)=temp
			GOTO 100
!
!	T44
!
80			temp=(rdn/(2.d0*pi*(RTQ(1)-RTQ(2))))*&
				 ( KQ1(1)*RTQ(1)*DRTQ(1)*( RTQ(1)+(g1*xn/xkw+sat*dfs(1)/&
				                                    ((xlambda+2.D0*xmu)*xkw))*SLAP )-&
				   KQ1(2)*RTQ(2)*DRTQ(2)*( RTQ(2)+(g1*xn/xka+sat*dfs(1)/&
				                                    ((xlambda+2.D0*xmu)*xkw))*SLAP )  )
!
			cqmF(i,j)=temp
			GOTO 100
!
100		continue
		enddo
	enddo
!
	return
	end

!*******************************************************************************************
!
	subroutine wunsatqstfsolG(SLAP,RTQ,KQ0,KQ1,cqmG)
!
	implicit double precision (a-h,o-z)
!
	COMPLEX(8) SLAP,RTQ(2),DRTQ(2),KQ0(2),KQ1(2),&
			   GAM11,GAM12,GAM13,GAM21,GAM22,GAM31,GAM32,GAM33,GAM1,GAM2,GAM3,temp,&
			   cqmG(kbem,kbem),R7(2)
!
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /e/ gammaw,gammas,gammaa,grav,atmp
	common /point/ dist(2),ra,rd(2),rdn,eta(2)
	common /UcoefR/ r1(2,2),dpr1(2,2,2),r2(2,2),dpr2(2,2,2),c(2,2)
	common /umate/ d(3,3),dfs(3),xn,xkw,cww,cwg,xka,cgg,cgw,ev,sat,g1,rs,rw,raa,rmix,xm1,&
				   bt,et,xlambda,xmu,henry
!
	dimension dlt(2,2),abscqmG(kbem,kbem)
!
!
	do i=1,kbem
		do j=1,kbem
			cqmG(i,j)=0.d0
		enddo
	enddo
!
	pi=DACOS (-1.d0)
	zmu=3.d0*bt*et/(9.d0*bt-et)
	zlamda=bt-2.d0*zmu/3.d0
	znu=zlamda/(2.d0*(zlamda+zmu))
!	dfs(1)=-dfs(1)
!
	dlt(1,1)=1.d0
	dlt(2,2)=1.d0
	dlt(1,2)=0.d0
	dlt(2,1)=0.d0
!
	do i=1,2
		DRTQ(i)=Cdsqrt( RTQ(i) )
	enddo
!
!	call WFSOLGCOEFQST(SLAP,CK11,CK12,CK13,CK21,CK22,CK23,CK31,CK32,CK41,CK42,CK51,CK52,CK61,&
!							 CK62,CK71,CK72,CK73,CK74,CK75,CK76,CK77)
!
!	GAM11=( KQ0(2)*RTQ(2)-KQ0(1)*RTQ(1) )/( (RTQ(2)-RTQ(1)) )
!	GAM12=SLAP*( KQ0(2)-KQ0(1) )/(RTQ(2)-RTQ(1))
!	GAM13=(SLAP**2.d0)*( KQ0(2)/RTQ(2)-KQ0(1)/RTQ(1) )/(RTQ(2)-RTQ(1))
!
!	GAM21=SLAP*( KQ1(2)*DRTQ(2)-KQ1(1)*DRTQ(1) )/(RTQ(2)-RTQ(1))
!	GAM22=(SLAP**2.d0)*( KQ1(2)/DRTQ(2)-KQ1(1)/DRTQ(1) )/(RTQ(2)-RTQ(1))
!
!	GAM31=( KQ1(2)*DRTQ(2)-KQ1(1)*DRTQ(1) )/(RTQ(2)-RTQ(1))
!	GAM32=SLAP*( KQ1(2)/DRTQ(2)-KQ1(1)/DRTQ(1) )/(RTQ(2)-RTQ(1))
!	GAM33=(SLAP**2.d0)*( KQ1(2)/DRTQ(2)**3.d0-KQ1(1)/DRTQ(1)**3.d0 )/(RTQ(2)-RTQ(1))
!
!	GAM1=CK11*GAM11+CK12*GAM12+CK13*GAM13
!	GAM2=CK21*GAM31+CK22*GAM32+CK23*GAM33
!	GAM3=CK21*GAM11+CK22*GAM12+CK23*GAM13
!
	do k=1,4
		do j=1,4
			temp=DCMPLX (0.d0,0.d0)
			if (k.lt.3.and.j.eq.3) GOTO 10
			if (k.lt.3.and.j.eq.4) GOTO 20
			if (k.eq.3.and.j.lt.3) GOTO 30
			if (k.eq.3.and.j.eq.3) GOTO 40
			if (k.eq.3.and.j.eq.4) GOTO 50
			if (k.eq.4.and.j.lt.3) GOTO 60
			if (k.eq.4.and.j.eq.3) GOTO 70
			if (k.eq.4.and.j.eq.4) GOTO 80
!
!	Gij
!
!			CK11=1.D0/(2.D0*pi*zmu)
!			CK12=( xn*g1*(zlamda+2.D0*zmu)*(xka*gammaw+xkw*gammaa)+xkw*gammaa*(1.d0-sat)*(-1.d0+dfs(1))-&
!				   xka*gammaw*sat*dfs(1) )/&
!				 (2.D0*pi*zmu*(zlamda+2.D0*zmu)*xkw*xka)
!			CK13=-xn*g1*gammaw*gammaa/(2.D0*pi*zmu*(zlamda+2.D0*zmu)*xkw*xka)
!			CK21=-(zlamda+zmu)/(2.D0*pi*zmu*(zlamda+2.D0*zmu))
!			CK22=-( xn*g1*(zlamda+zmu)*(xka*gammaw+xkw*gammaa)+xkw*gammaa*(1.d0-sat)*(-1.d0+dfs(1))-&
!					xka*gammaw*sat*dfs(1) )/&
!				  (2.D0*pi*zmu*(zlamda+2.D0*zmu)*xkw*xka)
!
!			CK11=1.D0/(2.D0*pi*zmu)
!			CK12=( xn*g1*(zlamda+2.D0*zmu)*(xka+xkw)+xkw*(1.d0-sat)*(-1.d0+dfs(1))-&
!				   xka*sat*dfs(1) )/&
!				 (2.D0*pi*zmu*(zlamda+2.D0*zmu)*xkw*xka)
!			CK13=-xn*g1/(2.D0*pi*zmu*(zlamda+2.D0*zmu)*xkw*xka)
!			CK21=-(zlamda+zmu)/(2.D0*pi*zmu*(zlamda+2.D0*zmu))
!			CK22=-( xn*g1*(zlamda+zmu)*(xka+xkw)+xkw*(1.d0-sat)*(-1.d0+dfs(1))-&
!					xka*sat*dfs(1) )/&
!				  (2.D0*pi*zmu*(zlamda+2.D0*zmu)*xkw*xka)
!			CK23=-CK13
!
!			temp=dlt(k,j)*(CK11*(KQ0(2)*RTQ(2)-KQ0(1)*RTQ(1))+CK12*SLAP*(KQ0(2)-KQ0(1))+&
!						   CK13*SLAP*SLAP*(KQ0(2)/RTQ(2)-KQ0(1)/RTQ(1)))/(RTQ(2)-RTQ(1))+&
!				 CK21*( (2.D0*rd(k)*rd(j)-dlt(k,j))*(KQ1(2)*DRTQ(2)-KQ1(1)*DRTQ(1))+&
!						 rd(k)*rd(j)*ra*(KQ0(2)*RTQ(2)-KQ0(1)*RTQ(1)) )/(ra*(RTQ(2)-RTQ(1)))+&
!				 CK22*SLAP*( (2.D0*rd(k)*rd(j)-dlt(k,j))*(KQ1(2)/DRTQ(2)-KQ1(1)/DRTQ(1))+&
!							  rd(k)*rd(j)*ra*(KQ0(2)-KQ0(1)) )/(ra*(RTQ(2)-RTQ(1)))+&
!				 CK23*SLAP*SLAP*( (2.D0*rd(k)*rd(j)-dlt(k,j))*(KQ1(2)/(DRTQ(2)*RTQ(2))-KQ1(1)/(DRTQ(1)*RTQ(1)))+&
!								   rd(k)*rd(j)*ra*(KQ0(2)/RTQ(2)-KQ0(1)/RTQ(1)) )/(ra*(RTQ(2)-RTQ(1)))
!
!
			temp=(dlt(k,j)/(RTQ(2)-RTQ(1)))*&
				 (  (KQ0(2)*RTQ(2)-KQ0(1)*RTQ(1))/(2.D0*pi*xmu)+&
					(sat*dfs(1)/((xlambda+2.D0*xmu)*xkw)+&
					 (1.d0-sat)*(1.d0-dfs(1))/((xlambda+2.D0*xmu)*xka)+&
					 xn*g1*(xka+xkw)/(xkw*xka))*SLAP*(KQ0(2)-KQ0(1))/(2.D0*pi*xmu)+&
					(g1*xn*SLAP*SLAP/(2.D0*pi*xmu*(xlambda+2.D0*xmu)*xkw*xka))*&
					(KQ0(2)/RTQ(2)-KQ0(1)/RTQ(1))  )   -&
!
				 (1.d0/(2.D0*pi*xmu*(xlambda+2.D0*xmu)*(RTQ(2)-RTQ(1))))*&
				 (  (xlambda+xmu)*((2.D0*rd(k)*rd(j)-dlt(k,j))*(KQ1(2)*DRTQ(2)-KQ1(1)*DRTQ(1))/ra+&
				                    rd(k)*rd(j)*(KQ0(2)*RTQ(2)-KQ0(1)*RTQ(1)))+&
				    (xn*g1*(xlambda+xmu)*(xka+xkw)/(xka*xkw)+(1.d0-sat)*(1.d0-dfs(1))/xka+&
					 sat*dfs(1)/xkw)*SLAP*((2.D0*rd(k)*rd(j)-dlt(k,j))*&
					                       (KQ1(2)/DRTQ(2)-KQ1(1)/DRTQ(1))/ra+&
				                            rd(k)*rd(j)*(KQ0(2)-KQ0(1)))+&
					(g1*xn*SLAP*SLAP/(xkw*xka))*((2.D0*rd(k)*rd(j)-dlt(k,j))*&
					                             (KQ1(2)/(RTQ(2)*DRTQ(2))-KQ1(1)/(RTQ(1)*DRTQ(1)))/ra+&
				                                  rd(k)*rd(j)*(KQ0(2)/RTQ(2)-KQ0(1)/RTQ(1)))  )
!
!
			cqmG(k,j)=temp
			goto 100
!
!	Gi3
!
10			temp=-( rd(k)/(2.d0*pi*(xlambda+2.D0*xmu)*rw*xkw*(RTQ(2)-RTQ(1))) ) *&
				       (   dfs(1)*(KQ1(2)*DRTQ(2)-KQ1(1)*DRTQ(1))+&
						   (g1*xn*SLAP/xka)*(KQ1(2)/DRTQ(2)-KQ1(1)/DRTQ(1))  )
!
!			temp=rd(k)*(1.d0/(2.d0*pi*(zlamda+2.D0*zmu)*xkw*xka*(RTQ(2)-RTQ(1))))*&
!				       (   (dfs(1)*xka)*&
!							      (KQ1(2)*DRTQ(2)-KQ1(1)*DRTQ(1))+&
!							g1*xn*SLAP*&
!								  (KQ1(2)/DRTQ(2)-KQ1(1)/DRTQ(1))  )
!
			cqmG(k,j)=temp
			goto 100
!
!	Gi4
!
20			temp=-( rd(k)/(2.d0*pi*(xlambda+2.D0*xmu)*raa*xka*(RTQ(2)-RTQ(1))) ) *&
				       (   (1.d0-dfs(1))*(KQ1(2)*DRTQ(2)-KQ1(1)*DRTQ(1))+&
						   (g1*xn*SLAP/xkw)*(KQ1(2)/DRTQ(2)-KQ1(1)/DRTQ(1))  )
!
!			temp=rd(k)*(1.d0/(2.d0*pi*(zlamda+2.D0*zmu)*xkw*xka*(RTQ(2)-RTQ(1))))*&
!				       (   ((1.d0-dfs(1))*xkw)*&
!							      (KQ1(2)*DRTQ(2)-KQ1(1)*DRTQ(1))+&
!							g1*xn*SLAP*&
!								  (KQ1(2)/DRTQ(2)-KQ1(1)/DRTQ(1))  )
!
			cqmG(k,j)=temp
			goto 100
!
!	G3j
!
30			temp=( rd(j)/(2.d0*pi*(xlambda+2.D0*xmu)*xkw*(RTQ(2)-RTQ(1))) )*&
					   ( sat*SLAP*(KQ1(2)*DRTQ(2)-KQ1(1)*DRTQ(1))+&
					     (g1*xn/xka)*SLAP*SLAP*(KQ1(2)/DRTQ(2)-KQ1(1)/DRTQ(1)) )
!
!			temp=rd(j)*(1.d0/(2.d0*pi*(zlamda+2.D0*zmu)*xkw*xka*(RTQ(2)-RTQ(1))))*&
!					   ( xka*gammaw*sat*SLAP*(KQ1(2)*DRTQ(2)-KQ1(1)*DRTQ(1))+&
!					     g1*gammaw*gammaa*xn*SLAP*SLAP*(KQ1(2)/DRTQ(2)-KQ1(1)/DRTQ(1)) )
!
			cqmG(k,j)=temp
			goto 100
!
!	G33
!
40			temp=(1.d0/(2.d0*pi*rw*xkw*(RTQ(2)-RTQ(1))))*&
						( (KQ0(2)*RTQ(2)-KQ0(1)*RTQ(1))+&
					      (g1*xn/xka+(1.d0-dfs(1))*(1.d0-sat)/(xka*(xlambda+2.D0*xmu)))*&
						  SLAP*(KQ0(2)-KQ0(1)) )
!
!			temp=gammaw*(1.d0/(2.d0*pi*(zlamda+2.D0*zmu)*rw*xkw*xka*(RTQ(2)-RTQ(1))))*&
!						(xka*(zlamda+2.D0*zmu)*(KQ0(2)*RTQ(2)-KQ0(1)*RTQ(1))+&
!					     gammaa*(g1*xn*(1.d0-henry)*(zlamda+2.D0*zmu)-&
!								 (1.d0-dfs(1))*(1.d0-sat*(1.d0-henry)))*SLAP*&
!								 (KQ0(2)-KQ0(1)) )
!
			cqmG(k,j)=temp
			goto 100
!
!	G34
!
50			temp=-(1.d0/(2.d0*pi*raa*xkw*xka*(RTQ(2)-RTQ(1))))*&
				  (sat*(1.d0-dfs(1))/(xlambda+2.D0*xmu)-g1*xn)*SLAP*(KQ0(2)-KQ0(1))

!			temp=((sat*(1.d0-dfs(1))+g1*xn*(zlamda+2.D0*zmu))/&
!					   (2.d0*pi*(zlamda+2.D0*zmu)*xkw*xka*(RTQ(2)-RTQ(1))))*&
!						SLAP*(KQ0(2)-KQ0(1))
!
			cqmG(k,j)=temp
			goto 100
!
!	G4j
!
60			temp=rd(j)*(1.d0/(2.d0*pi*(xlambda+2.D0*xmu)*xka*(RTQ(2)-RTQ(1))))*&
					   ( (1.d0-sat)*SLAP*(KQ1(2)*DRTQ(2)-KQ1(1)*DRTQ(1))+&
					     (g1*xn/xkw)*SLAP*SLAP*(KQ1(2)/DRTQ(2)-KQ1(1)/DRTQ(1)) )

!			temp=rd(j)*(1.d0/(2.d0*pi*(zlamda+2.D0*zmu)*xkw*xka*(RTQ(2)-RTQ(1))))*&
!					   ( xkw*(1.d0-sat)*SLAP*(KQ1(2)*DRTQ(2)-KQ1(1)*DRTQ(1))+&
!					     g1*xn*SLAP*SLAP*(KQ1(2)/DRTQ(2)-KQ1(1)/DRTQ(1)) )
!
			cqmG(k,j)=temp
			goto 100
!
!	G43
!
70			temp=-(1.d0/(2.d0*pi*rw*xkw*xka*(RTQ(2)-RTQ(1))))*&
				  ( (1.d0-sat)*dfs(1)/(xlambda+2.D0*xmu)-g1*xn)*SLAP*(KQ0(2)-KQ0(1))
!
!			temp=(  (g1*xn*(zlamda+2.D0*zmu)+dfs(1)*(1.d0-sat))/&
!					(2.d0*pi*(zlamda+2.D0*zmu)*xkw*xka*(RTQ(2)-RTQ(1)))  )*&
!				  SLAP*(KQ0(2)-KQ0(1))
!
			cqmG(k,j)=temp
			goto 100
!
!	G44
!
80			temp=(1.d0/(2.d0*pi*raa*xka*(RTQ(2)-RTQ(1))))*&
						( (KQ0(2)*RTQ(2)-KQ0(1)*RTQ(1))+&
					      (g1*xn/xkw+dfs(1)*sat/(xkw*(xlambda+2.D0*xmu)))*&
						  SLAP*(KQ0(2)-KQ0(1)) )
!
!			temp=-(1.d0/(2.d0*pi*(zlamda+2.D0*zmu)*xkw*xka*(RTQ(2)-RTQ(1))))*&
!						( -xkw*(zlamda+2.D0*zmu)*(KQ0(2)*RTQ(2)-KQ0(1)*RTQ(1))+&
!					      (dfs(1)*sat-g1*xn*(zlamda+2.D0*zmu))*SLAP*(KQ0(2)-KQ0(1)) )

!
			cqmG(k,j)=temp
			goto 100
!
100			continue
!
		enddo
	enddo
!
	return
	end subroutine
!
!*******************************************************************************************


!*******************************************************************************************
!
	subroutine wunsatqstfsolG1(SLAP,RTQ,KQ0,KQ1,KQ2,cqmG)
!
	implicit double precision (a-h,o-z)
!
	COMPLEX(8) SLAP,RTQ(2),DRTQ(2),KQ0(2),KQ1(2),KQ2(2),temp,cqmG(kbem,kbem),CK11,CK22,CK33,CK44
!
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /e/ gammaw,gammas,gammaa,grav,atmp
	common /point/ dist(2),ra,rd(2),rdn,eta(2)
	common /UcoefR/ r1(2,2),dpr1(2,2,2),r2(2,2),dpr2(2,2,2),c(2,2)
	common /umate/ d(3,3),dfs(3),xn,xkw,cww,cwg,xka,cgg,cgw,ev,sat,g1,rs,rw,raa,rmix,xm1,&
				   bt,et,xlambda,xmu,henry
!
	dimension dlt(2,2),abscqmG(kbem,kbem)
!
!
	do i=1,kbem
		do j=1,kbem
			cqmG(i,j)=0.d0
		enddo
	enddo
!
	pi=DACOS (-1.d0)
	dlt(1,1)=1.d0
	dlt(2,2)=1.d0
	dlt(1,2)=0.d0
	dlt(2,1)=0.d0
!
	do i=1,2
		DRTQ(i)=Cdsqrt( RTQ(i) )
	enddo
!
	do k=1,4
		do j=1,4
			temp=DCMPLX (0.d0,0.d0)
			if (k.lt.3.and.j.eq.3) GOTO 10
			if (k.lt.3.and.j.eq.4) GOTO 20
			if (k.eq.3.and.j.lt.3) GOTO 30
			if (k.eq.3.and.j.eq.3) GOTO 40
			if (k.eq.3.and.j.eq.4) GOTO 50
			if (k.eq.4.and.j.lt.3) GOTO 60
			if (k.eq.4.and.j.eq.3) GOTO 70
			if (k.eq.4.and.j.eq.4) GOTO 80
!
!	Gij
!
			CK11=( xn*g1*(xka+xkw)/(xka*xkw)+&
				  (1.d0-sat)*(1.d0-dfs(1))/((xlambda+2.D0*xmu)*xka)+&
				   sat*dfs(1)/((xlambda+2.D0*xmu)*xkw) )*SLAP
			CK22=xn*g1*SLAP*SLAP/((xlambda+2.D0*xmu)*xkw*xka)
			CK33=( (xlambda+xmu)*xn*g1*(xka+xkw)/(xka*xkw)+&
				   (1.d0-sat)*(1.d0-dfs(1))/xka+&
				   sat*dfs(1)/xkw )*SLAP
			CK44=xn*g1*SLAP*(xkw+xka)/(xkw*xka)
!
			R7=( (2.d0*rd(K)*rd(J)-dlt(K,J))/ra )
			R8=rd(K)*rd(J)

!
!			temp=( dlt(k,j)/(2.d0*pi*xmu*(RTQ(2)-RTQ(1))) )*&
!				 (  (KQ0(2)*RTQ(2)-KQ0(1)*RTQ(1))+&
!					 CK11*(KQ0(2)-KQ0(1))+&
!					 CK22*(KQ0(2)/RTQ(2)-KQ0(1)/RTQ(1))  )   -&
!
!				 ( 1.d0/(2.D0*pi*xmu*(xlambda+2.D0*xmu)*(RTQ(2)-RTQ(1))) )*&
!				 (  (xlambda+xmu)*( R7*(KQ1(2)*DRTQ(2)-KQ1(1)*DRTQ(1))+&
!				                    R8*(KQ0(2)*RTQ(2)-KQ0(1)*RTQ(1)) )+&
!				    CK33*( R7*(KQ1(2)/DRTQ(2)-KQ1(1)/DRTQ(1))+&
!				           R8*(KQ0(2)-KQ0(1)) )+&
!					CK44*( R7*(KQ1(2)/(RTQ(2)*DRTQ(2))-KQ1(1)/(RTQ(1)*DRTQ(1)))+&
!				           R8*(KQ0(2)/RTQ(2)-KQ0(1)/RTQ(1)) )  )
!
			temp=-dlt(K,J)*dlog(ra)/(4.d0*pi*xmu)+R8/(4.d0*pi*xmu)+&
		         (1.d0/(2.d0*pi*(xlambda+2.D0*xmu)*(RTQ(1)-RTQ(2))))*&
				 ( R8*(RTQ(1)*KQ2(1)-RTQ(2)*KQ2(2))-dlt(K,J)*(KQ1(1)*DRTQ(1)/ra-KQ1(2)*DRTQ(2)/ra)+&
				   CK44*( R8*(KQ2(1)-KQ2(2))-dlt(K,J)*(KQ1(1)/(DRTQ(1)*ra)-KQ1(2)/(DRTQ(2)*ra)) +&
                   R7*(RTQ(1)-RTQ(2))/(ra*RTQ(1)*RTQ(2)) ) )
!
			cqmG(k,j)=temp
			goto 100
!
!	Gi3
!
10			temp=-( rd(k)/(2.d0*pi*rw*SLAP*(RTQ(1)-RTQ(2))) ) *&
				       (   (RTQ(2)+dfs(1)*SLAP/(xkw*(xlambda+2.D0*xmu)))*KQ1(1)*DRTQ(1)-&
					       (RTQ(1)+dfs(1)*SLAP/(xkw*(xlambda+2.D0*xmu)))*KQ1(2)*DRTQ(2)+&
						   (RTQ(1)-RTQ(2))/ra  )
!
			cqmG(k,j)=temp
			goto 100
!
!	Gi4
!
20			temp=-( rd(k)/(2.d0*pi*raa*SLAP*(RTQ(1)-RTQ(2))) ) *&
				       (   (RTQ(2)+(1.d0-dfs(1))*SLAP/(xka*(xlambda+2.D0*xmu)))*KQ1(1)*DRTQ(1)-&
					       (RTQ(1)+(1.d0-dfs(1))*SLAP/(xka*(xlambda+2.D0*xmu)))*KQ1(2)*DRTQ(2)+&
						   (RTQ(1)-RTQ(2))/ra  )
!
			cqmG(k,j)=temp
			goto 100
!
!	G3j
!
30			temp=( rd(j)/(2.d0*pi*(RTQ(1)-RTQ(2))) )*&
					   ( (RTQ(2)+sat*SLAP/(xkw*(xlambda+2.D0*xmu)))*KQ1(1)*DRTQ(1)-&
					     (RTQ(1)+sat*SLAP/(xkw*(xlambda+2.D0*xmu)))*KQ1(2)*DRTQ(2)+&
					     (RTQ(1)-RTQ(2))/ra )
!
			cqmG(k,j)=temp
			goto 100
!
!	G33
!
40			temp=(1.d0/(2.d0*pi*rw*xkw*(RTQ(1)-RTQ(2))))*&
						( (RTQ(1)+(g1*xn/xka+(1.d0-dfs(1))*(1.d0-sat)/(xka*(xlambda+2.D0*xmu)))*SLAP)*KQ0(1)-&
						  (RTQ(2)+(g1*xn/xka+(1.d0-dfs(1))*(1.d0-sat)/(xka*(xlambda+2.D0*xmu)))*SLAP)*KQ0(2) )
!
			cqmG(k,j)=temp
			goto 100
!
!	G34
!
50			temp=(-1.d0/(2.d0*pi*raa*xkw*xka*(RTQ(1)-RTQ(2))))*&
				  (sat*(1.d0-dfs(1))/(xlambda+2.D0*xmu)-g1*xn)*SLAP*(KQ0(1)-KQ0(2))
!
			cqmG(k,j)=temp
			goto 100
!
!	G4j
!
60			temp=( rd(j)/(2.d0*pi*(RTQ(1)-RTQ(2))) )*&
					   ( (RTQ(2)+(1.d0-sat)*SLAP/(xka*(xlambda+2.D0*xmu)))*KQ1(1)*DRTQ(1)-&
					     (RTQ(1)+(1.d0-sat)*SLAP/(xka*(xlambda+2.D0*xmu)))*KQ1(2)*DRTQ(2)+&
					     (RTQ(1)-RTQ(2))/ra )
!
			cqmG(k,j)=temp
			goto 100
!
!	G43
!
70			temp=-(1.d0/(2.d0*pi*rw*xkw*xka*(RTQ(1)-RTQ(2))))*&
				  ( (1.d0-sat)*dfs(1)/(xlambda+2.D0*xmu)-g1*xn)*SLAP*(KQ0(1)-KQ0(2))
!
			cqmG(k,j)=temp
			goto 100
!
!	G44
!
80			temp=(1.d0/(2.d0*pi*raa*xka*(RTQ(1)-RTQ(2))))*&
						( (RTQ(1)+(g1*xn/xkw+dfs(1)*sat/(xkw*(xlambda+2.D0*xmu)))*SLAP)*KQ0(1)-&
						  (RTQ(2)+(g1*xn/xkw+dfs(1)*sat/(xkw*(xlambda+2.D0*xmu)))*SLAP)*KQ0(2) )
!
			cqmG(k,j)=temp
			goto 100
!
100			continue
!
		enddo
	enddo
!
	return
	end subroutine
!
!*******************************************************************************************

!*******************************************************************************************
!
	subroutine wunsatqstfsolG2(SLAP,RTQ,KQ0,KQ1,cqmG)
!
	implicit double precision (a-h,o-z)
!
	COMPLEX(8) SLAP,RTQ(2),DRTQ(2),KQ0(2),KQ1(2),temp,cqmG(kbem,kbem),CK11,CK22,CK33,CK44
!
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /e/ gammaw,gammas,gammaa,grav,atmp
	common /point/ dist(2),ra,rd(2),rdn,eta(2)
	common /UcoefR/ r1(2,2),dpr1(2,2,2),r2(2,2),dpr2(2,2,2),c(2,2)
	common /umate/ d(3,3),dfs(3),xn,xkw,cww,cwg,xka,cgg,cgw,ev,sat,g1,rs,rw,raa,rmix,xm1,&
				   bt,et,xlambda,xmu,henry
!
	dimension dlt(2,2),abscqmG(kbem,kbem)
!
!
	do i=1,kbem
		do j=1,kbem
			cqmG(i,j)=0.d0
		enddo
	enddo
!
	pi=DACOS (-1.d0)
	dlt(1,1)=1.d0
	dlt(2,2)=1.d0
	dlt(1,2)=0.d0
	dlt(2,1)=0.d0
!
	do i=1,2
		DRTQ(i)=Cdsqrt( RTQ(i) )
	enddo
!
	do k=1,4
		do j=1,4
			temp=DCMPLX (0.d0,0.d0)
			if (k.lt.3.and.j.eq.3) GOTO 10
			if (k.lt.3.and.j.eq.4) GOTO 20
			if (k.eq.3.and.j.lt.3) GOTO 30
			if (k.eq.3.and.j.eq.3) GOTO 40
			if (k.eq.3.and.j.eq.4) GOTO 50
			if (k.eq.4.and.j.lt.3) GOTO 60
			if (k.eq.4.and.j.eq.3) GOTO 70
			if (k.eq.4.and.j.eq.4) GOTO 80
!
!	Gij
!
			R7=( (2.d0*rd(K)*rd(J)-dlt(K,J))/ra )
			R8=rd(K)*rd(J)

!
			temp=( (-1.d0/(2.D0*pi*xmu))*R7/ra )+&
				 ( 1.d0/(2.D0*pi*(xlambda+2.D0*xmu)*(RTQ(1)-RTQ(2))) )*&
				 (  ( R7*(KQ1(1)*DRTQ(1)*RTQ(1)-KQ1(2)*DRTQ(2)*RTQ(2))+&
				      R8*(KQ0(1)*RTQ(1)*RTQ(1)-KQ0(2)*RTQ(2)*RTQ(2)) )+&
				    ( xn*g1*(xka+xkw)*SLAP/(xkw*xka) )*&
					( R7*(KQ1(1)*DRTQ(1)-KQ1(2)*DRTQ(2))+R8*(RTQ(1)*KQ0(1)-RTQ(2)*KQ0(2)) )  )
!
			cqmG(k,j)=temp
			goto 100
!
!	Gi3
!
10			temp=-( rd(k)/(2.d0*pi*(xlambda+2.D0*xmu)*rw*xkw*(RTQ(1)-RTQ(2))) ) *&
				       (   dfs(1)*(KQ1(1)*DRTQ(1)*RTQ(1)-KQ1(2)*DRTQ(2)*RTQ(2))+&
						   (g1*xn*SLAP/xka)*(KQ1(1)*DRTQ(1)-KQ1(2)*DRTQ(2))  )
!
			cqmG(k,j)=temp
			goto 100
!
!	Gi4
!
20			temp=-( rd(k)/(2.d0*pi*(xlambda+2.D0*xmu)*raa*xka*(RTQ(1)-RTQ(2))) ) *&
				       (   (1.d0-dfs(1))*(KQ1(1)*DRTQ(1)*RTQ(1)-KQ1(2)*DRTQ(2)*RTQ(2))+&
						   (g1*xn*SLAP/xkw)*(KQ1(1)*DRTQ(1)-KQ1(2)*DRTQ(2))  )
!
			cqmG(k,j)=temp
			goto 100
!
!	G3j
!
30			temp=( rd(j)/(2.d0*pi*(xlambda+2.D0*xmu)*xkw*(RTQ(1)-RTQ(2))) )*&
					   ( sat*SLAP*(KQ1(1)*DRTQ(1)*RTQ(1)-KQ1(2)*DRTQ(2)*RTQ(2))+&
					     (g1*xn/xka)*SLAP*SLAP*(KQ1(1)*DRTQ(1)-KQ1(2)*DRTQ(2)) )
!
			cqmG(k,j)=temp
			goto 100
!
!	G33
!
40			temp=(1.d0/(2.d0*pi*rw*xkw*(RTQ(1)-RTQ(2))))*&
						( (KQ0(1)*RTQ(1)*RTQ(1)-KQ0(2)*RTQ(2)*RTQ(2))+&
					      (g1*xn/xka+(1.d0-dfs(1))*(1.d0-sat)/(xka*(xlambda+2.D0*xmu)))*&
						  SLAP*(RTQ(1)*KQ0(1)-RTQ(2)*KQ0(2)) )
!
			cqmG(k,j)=temp
			goto 100
!
!	G34
!
50			temp=-(1.d0/(2.d0*pi*raa*xkw*xka*(RTQ(1)-RTQ(2))))*&
				  (sat*(1.d0-dfs(1))/(xlambda+2.D0*xmu)-g1*xn)*SLAP*(RTQ(1)*KQ0(1)-RTQ(2)*KQ0(2))
!
			cqmG(k,j)=temp
			goto 100
!
!	G4j
!
60			temp=rd(j)*(1.d0/(2.d0*pi*(xlambda+2.D0*xmu)*xka*(RTQ(1)-RTQ(2))))*&
					   ( (1.d0-sat)*SLAP*(KQ1(1)*DRTQ(1)*RTQ(1)-KQ1(2)*DRTQ(2)*RTQ(2))+&
					     (g1*xn/xkw)*SLAP*SLAP*(KQ1(1)*DRTQ(1)-KQ1(2)*DRTQ(2)) )
!
			cqmG(k,j)=temp
			goto 100
!
!	G43
!
70			temp=-(1.d0/(2.d0*pi*rw*xkw*xka*(RTQ(1)-RTQ(2))))*&
				  ( (1.d0-sat)*dfs(1)/(xlambda+2.D0*xmu)-g1*xn)*SLAP*(RTQ(1)*KQ0(1)-RTQ(2)*KQ0(2))
!
			cqmG(k,j)=temp
			goto 100
!
!	G44
!
80			temp=(1.d0/(2.d0*pi*raa*xka*(RTQ(1)-RTQ(2))))*&
						( (KQ0(1)*RTQ(1)*RTQ(1)-KQ0(2)*RTQ(2)*RTQ(2))+&
					      (g1*xn/xkw+dfs(1)*sat/(xkw*(xlambda+2.D0*xmu)))*&
						  SLAP*(RTQ(1)*KQ0(1)-RTQ(2)*KQ0(2)) )
!
			cqmG(k,j)=temp
			goto 100
!
100			continue
!
		enddo
	enddo
!
	return
	end subroutine
!
!*******************************************************************************************

!****************************************************************************************
!
!							           BEM subroutine
!
!	This subroutine is the main controller of BEM in problem.
!	This subroutine is called in: CONTROL
!	This subroutine calls:		  BEMINIT,BEMT1,BEMTN,CHARGE5,CHARGE6,BEMS1S2
!
!	variables used are:
!
!		id		: degrees of freedom of each node from "1" to "mdof"; for the nodes with
!				  the boundary conditions, each node takes from "mdof+1" to "mdof+nbc"
!				  [X,Y]; which are free for the dry elm &
!				  [X,Y,Pw]; which are free for the saturated elm
!				  [X,Y,Pw,Pa]; which are free for the unsaturated elm
!		mdof	: maximum degree of freedom by considering the boundary conditions
!		mdofn	: maximum degree of freedom without considering the boundary conditions
!		na		: this locates the diagonal coefficient for skyline storage. we economise
!				  the space by neglecting the stockage of certain coefs which are zero
!				  because of the absence of the connectivity of doF. we delete just the zeros
!				  lied to the non connectivity. therefore na(i) shows the i-th coef to stock.
!		npbmax	: maximum nodes in all of the BE zone + 1; for example if we consider D3 BE
!				  area NPBMAX = nnpbed3 + 1
!		n12gh	: in geomechanical problems there is more than one unknown per node, so the
!				  matrices in BEM are assambled according to the doF, rather than node
!				  number. for dry soil each node has [2] doF (ux,uy), for for saturated soil
!				  each node has [3] (ux,uy,pw) doF & for unsaturated soil each node has
!				  [4] doF(ux,uy,pw,pa). then total number of doF in BE=(npbmax-1)*kbem
!		ntime	: number of time steps or total number of "nstep" in all "nload"s
!				  loading is divided into SUM[nstep]s, when the characteristics of "nstep"s
!				  are similar we put these "nsteps" in one group which is named iload
!		rrr		:
!		disp	: displacement/pressures vector which consists of initial displacement value
!				  of corresponding nodes in which their identities code of Ux & Uy are defined
!				  in ID and initial water and air values of corresponding nodes in which their
!				  identities code of Pw & Pa are defined in ID. these displacement initial
!				  values can be defined for a free node or for a node with boundary conditions
!				  and the water and air initial values must be defined for a free node;
!				  in each "itime" the displacements/pressures vector (R->DT) is transfered
!			      in DISP vector
!		gid1	: inverse matrice of gd1, G1^(-1)
!		gid2	: inverse matrice of gd2, G2^(-1)
!		gid3	: inverse matrice of gd3, G3^(-1)
!		gihd1	: gid1*hd1, G1^(-1)*H1
!		gihd2	: gid2*hd2, G2^(-1)*H2
!		gihd3	: gid3*hd3, G3^(-1)*H3
!		xgid1	: M*gid1, M*G1^(-1)
!		xgid2	: M*gid2, M*G2^(-1)
!		xgid3	: M*gid3, M*G3^(-1)
!		xgihd1	: M*gid1*hd1, M*G1^(-1)*H1
!		xgihd2	: M*gid1*hd1, M*G1^(-1)*H2
!		xgihd3	: M*gid1*hd1, M*G1^(-1)*H3
!		frbed1	: nodal forces in D1 BEM area (finite zone)
!		frbed2	: nodal forces in D2 BEM area (infinite zone)
!		frbed3	: nodal forces in D3 BEM area (semi-infinite zone)
!		rbed1	: stress tensor in D1 BEM area which includes JUST effect of the past
!				  dynamic history; Z(N)=SUM(n=1:N-1) [G(N+1-n)*T(n)-H(N+1-n)*U(n)];
!				  in the first time step Z(1)=rbed1(1,ii)=0
!		rbed2	: stress tensor in D2 BEM area which includes JUST effect of the past
!				  dynamic history; Z(N)=SUM(n=1:N-1) [G(N+1-n)*T(n)-H(N+1-n)*U(n)];
!				  in the first time step Z(1)=rbed2(1,ii)=0
!		rbed3	: stress tensor in D3 BEM area which includes JUST effect of the past
!				  dynamic history; Z(N)=SUM(n=1:N-1) [G(N+1-n)*T(n)-H(N+1-n)*U(n)];
!				  in the first time step Z(1)=rbed3(1,ii)=0
!		tbed1	: traction/flow [tx,ty,qw,qa] vector T(N) at time: "itime-1"
!		tbed2	: traction/flow [tx,ty,qw,qa] vector T(N) at time: "itime-1"
!		tbed3	: traction/flow [tx,ty,qw,qa] vector T(N) at time: "itime-1"
!		ubed1	: displacement/pressure [ux,uy,Pw,Pa] vector U(N) at time: "itime-1" which
!				  is calculated in the previous time step in SOLVE subroutine &
!				  replaced in "DISP" vector in UPDATE3 subroutine
!		ubed2	: displacement/pressure [ux,uy,Pw,Pa] vector U(N) at time: "itime-1" which
!				  is calculated in the previous time step in SOLVE subroutine &
!				  replaced in "DISP" vector in UPDATE3 subroutine
!		ubed3	: displacement/pressure [ux,uy,Pw,Pa] vector U(N) at time: "itime-1" which
!				  is calculated in the previous time step in SOLVE subroutine &
!				  replaced in "DISP" vector in UPDATE3 subroutine
!		ls		: lenght of the global rigidity matrix in symmetric case. if the matrix is
!				  not symmetric, 2 triangles must be filled.
!		ls1		: if the problem is non symmetric then this partie must be taken into account
!				  then for this partie there is the same number of coefs to stocke like ls.
!		s1		:
!		s2		:
!		iter	: iteration in each time step (ITIME) which can be varied between 1 & ITMAX,
!				  but in this program it is considered always equal 1;
!							 [1] : first iteration
!		itime	: time increment from 1 to "ntime"
!							[1] : first time increment,
!							[>1]: higher time increment
!		itm		:
!
!
!	INPUT	:	id,na,x,sm8d,sm8c,sm8u,ie3d1,ie3d2,ie3d3,rrr,disp,
!
!
!	OUTPUT	:	s1,s2,xgihd1,xgihd2,xgihd3,gihd1,gihd2,gihd3,xgid1,xgid2,xgid3,gid1,gid2,
!				gid3,tbed1,tbed2,tbed3,ubed1,ubed2,ubed3,rbed1,rbed2,rbed3
!
!
!****************************************************************************************
!
	subroutine bem(id,na,x,xint,sm8d,sm8c,sm8u,s1,s2,ie3d1,ie3d2,ie3d3,sig3u,r,rrr,disp,&
				   xgihd1,xgihd2,xgihd3,gihd1,gihd2,gihd3,xgid1,xgid2,xgid3,gid1,gid2,gid3,&
				   tbed1,tbed2,tbed3,ubed1,ubed2,ubed3,rbed1,rbed2,rbed3)
!
	implicit double precision (a-h,o-z)
!
	common /a/ nnp,ifem,ibem,idyn,ieaq,nr
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /a3/ ne8d1,ne8c1,ne8u1,ne3u1,nload1,nbcx1,nbcy1,nbcw1,nbca1,m8d1,m8c1,m8u1,&
				nbed11,nbed21,nbed31,nbeint1,ne8dd,ne8cc,ne8uu
	common /a4/ xencl(500),yencl(500)
	common /b/ m8d,m8c,m8u,m8dep(10),m8cep(10),m8uep(10)
	common /c/ nload,ntime,itmax,itime,iter,iload,time,dtime,dtimei,dtimem,xinc
	common /h/ mdof,mdofn,ls,ls1,isolv
	common /g/ init,iprint
	common /chrg2/ npw,npa,nnpw(200),nnpa(200),vnpw(200),vnpa(200)
!
	dimension id(4,nnp),na(mdof),x(2,nnp),xint(2,nbeint1),sm8d(60,m8d1),sm8c(60,m8c1),&
			  sm8u(60,m8u1),s1(ls),s2(ls1),ie3d1(npbmax,nbed11),ie3d2(npbmax,nbed21),&
			  ie3d3(npbmax),sig3u(5,500,nbed11),rrr(mdof),disp(mdofn),&
			  xgihd2(n12gh,n12gh,nbed21),xgihd3(n12gh,n12gh),&
			  gihd2(n12gh,n12gh,nbed21),gihd3(n12gh,n12gh),&
			  gid2(n12gh,n12gh,nbed21),gid3(n12gh,n12gh),&
			  xgid2(n12gh,n12gh,nbed21),xgid3(n12gh,n12gh),tbed1(ntime,n12gh,nbed11),&
			  tbed2(ntime,n12gh,nbed21),tbed3(ntime,n12gh),ubed1(ntime,n12gh,nbed11),&
			  ubed2(ntime,n12gh,nbed21),ubed3(ntime,n12gh),rbed1(ntime,n12gh,nbed11),&
			  rbed2(ntime,n12gh,nbed21),rbed3(ntime,n12gh),frbed1(n12gh,nbed11),&
			  frbed2(n12gh,nbed21),frbed3(n12gh),feaq(n12gh),r(mdof),&
			  xgihd1(n12gh,n12gh,nbed11),&
			  gihd1(n12gh,n12gh,nbed11),&
			  gid1(n12gh,n12gh,nbed11),&
			  xgid1(n12gh,n12gh,nbed11)
!
!
	if (iter.gt.1) goto 80
	if (itime.gt.1) goto 50
!
!
!	----------------------- control bem in first time increment -----------------------
!
	call beminit(rbed1,rbed2,rbed3,frbed1,frbed2,frbed3)
	write (*,*)'		beminit concluded'
	itm=1
	call bemt1(itm,x,sm8d,sm8c,sm8u,ie3d1,ie3d2,ie3d3,xgihd1,xgihd2,xgihd3,gihd1,gihd2,gihd3,&
			   gid1,gid2,gid3,xgid1,xgid2,xgid3,disp,id,sig3u)
	write (*,*)'		bemt1 concluded'
!
	if (init.ne.0) goto 80	! initial stresses are unknown, the necessary parameters are
!							  obtained in BEMT1 subroutine
	goto 60
!
!
!	---------------------- control bem in higher time increments -----------------------
!
50	call bemtn(id,x,sm8d,sm8c,sm8u,ie3d1,ie3d2,ie3d3,gihd1,gihd2,gihd3,gid1,gid2,gid3,&
			   xgihd1,xgihd2,xgihd3,xgid1,xgid2,xgid3,tbed1,tbed2,tbed3,ubed1,ubed2,ubed3,&
			   rbed1,rbed2,rbed3,frbed1,frbed2,frbed3,disp,sig3u)
	write (*,*)'     bemtn concluded'
!
!
!	---------------- apply seismic loading at bem halfspace (d3) nodes -----------------
!
60	if (ieaq.ne.3) goto 70
	call charge5(feaq,xgid3,ie3d3)
	write (*,*)'     charge5 concluded'
!
!
!	-- assemble bem nodal (internal & seismic) force vectors and transfer it to force --
!	------------------------------------- vecor rrr ------------------------------------
!
70	if (idyn.eq.0) goto 80
	call charge6(id,ie3d1,ie3d2,ie3d3,rrr,frbed1,frbed2,frbed3,feaq)
	write (*,*)'     charge6 concluded'
!
!
!	----- assemble bem stiffness matrices into global vectors s1 & s2 -----
!
80	call bems1s2(id,na,s1,s2,ie3d1,ie3d2,ie3d3,xgihd1,xgihd2,xgihd3,r)
	write (*,*)'     bems1s2 concluded'
!
100	return
	end

!****************************************************************************************
!
!							           OTHERMAT3 subroutine
!
!	This subroutine calculates the gihd3 & xgihd3 matrices.
!	This subroutine is called in: GHMATD3
!	This subroutine calls: -
!
!	variables used are:
!
!		n			: total number of nodes of one D3 (semi-infinite) area NNPBED3
!		kbem		: number of degree of freedom per node in BE zone;
!								if IBEM=[1] -> KBEM=[2] [X,Y]
!								if IBEM=[2] -> KBEM=[3] [X,Y,Pw]
!								if IBEM=[3] -> KBEM=[4] [X,Y,Pw,Pa]
!		hd3			: fundamental solution matrix of stress/flow;
!					  total H matrix in D3 (semi-infinite) area concerning all of the h
!					  matrices of each element (ehd3). global coefficient h matrices
!					  assembled by gathering element contributions.
!					  in this subroutine, we consider just the nnpbed3 first lines of this
!					  matrix (square matrix).
!		gid3		: inverse matrice of gd3(kbem*nnpbed3,kbem*nnpbed3)
!		xmtf		: M (distributed matric) of all of elements
!		gihd3		: gid3*hd3
!		xgid3		: M*gid3
!		xgihd3		: M*gid3*hd3
!
!
!	INPUT	:
!
!	OUTPUT	:
!
!****************************************************************************************
!
	subroutine othermat3(n,hd3,xmtf,gid3,gihd3,xgid3,xgihd3)
!
	implicit double precision (a-h,o-z)
!
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
!
	dimension hd3(n12gh,n3gh),xmtf(n12gh,n12gh),gid3(n12gh,n12gh),gihd3(n12gh,n12gh),&
			  xgid3(n12gh,n12gh),xgihd3(n12gh,n12gh)
!
!
!	------------------- initializing -------------------
!
	n1=kbem*n
	do iy=1,n1
		do ix=1,n1
			gihd3(iy,ix)=0.d0
			xgid3(iy,ix)=0.d0
			xgihd3(iy,ix)=0.d0
		enddo
	enddo
!
!
!	------------------- gihd3 matrice -------------------
!
	do iy=1,n1
		do ix=1,n1
			do i=1,n1
				gihd3(iy,ix)=gihd3(iy,ix)+gid3(iy,i)*hd3(i,ix)
			enddo
		enddo
	enddo
!
!
!	------------------- xgid3 matrice -------------------
!
	do iy=1,n1
		do ix=1,n1
			do i=1,n1
				xgid3(iy,ix)=xgid3(iy,ix)+xmtf(iy,i)*gid3(i,ix)
			enddo
		enddo
	enddo
!
!
!	------------------- xgihd3 matrice -------------------
!
	do iy=1,n1
		do ix=1,n1
			do i=1,n1
				xgihd3(iy,ix)=xgihd3(iy,ix)+xmtf(iy,i)*gihd3(i,ix)
			enddo
		enddo
	enddo
!
	return
	end
!****************************************************************************************
!
!							           GHMATD3 subroutine
!
!	This subroutine calculates the g & h matrices of area D3 (semi-infinite zone)
!	This subroutine is called in: BEMT1
!	This subroutine calls:		  GHsaveD,SNGTP,EXTINeq,DINVERSE,MTRFRD,OTHERMAT3
!
!	variables used are:
!
!		itm		:
!		x		: coordinates of nodes (nodal data) in FE & BE areas;
!							x(1,1:nnp) : X-coord
!							x(2,1:nnp) : Y-coord
!		sm8d	: all of the parameters corresponding to each behaviour of dry soil are
!				  defined in this matrix.
!						  icpt:[1]-> sm8d (3,i) = E
!									 sm8d (7,i) = K
!									 sm8d (10,i)= Ko
!									 sm8d (11,i)= gammas
!						  icpt:[2]-> sm8d (1,i) = friction angle (�)
!									 sm8d (2,i) = cohesion (Pa)
!									 sm8d (3,i) = E-load, loading modulus used in non-linear
!												  elasticity (adimensionel)
!									 sm8d (4,i) = E-unld, unloading modulus used in
!												  non-linear elasticity (adim)
!									 sm8d (5,i) = n, exponent used to compute the Young
!												  loading and unloading moduli in non-linear
!											      elasticity (adim)
!									 sm8d (6,i) = Rf, coefficient multiplying the ratio
!												  (sigm1-sigm3)/(sigm1-sigm3)ult
!												  in nlin els (adim)
!                                    sm8d (7,i) = Kb, bulk modulus used in non-linear
!											      elasticity (adim)
!	                                 sm8d (8,i) = m, exponent used to compute the bulk
!										          modulus B in nlin els (adim)
!									 sm8d (9,i) = e-f, Young modulus minimal value (Pa)
!									 sm8d (10,i)=Ko, steady state coefficient (adim)
!								     sm8d (11,i)= gammas, soil volumetric weight (N,m-3)
!								     sm8d (12,i)= ten-st, traction resistency (Pa)
!						  icpt:[3]-> sm8d (1,i) = friction angle
!									 sm8d (2,i) = cohesion
!								     sm8d (3,i) = OCR
!								     sm8d (4,i) = k
!                                    sm8d (5,i) = n, exponent used to compute the Young
!												  loading and unloading moduli in non-linear
!												  elasticity (adim)
!									 sm8d (6,i) = e
!									 sm8d (7,i) = Kb, bulk modulus used in non-linear
!												  elasticity (adim)
!									 sm8d (8,i) = m, exponent used to compute the bulk
!												  modulus B in nlin els (adim)
!									 sm8d (10,i) =Ko
!									 sm8d (11,i) =gamma
!						  icpt:[4]-> ??????
!		sm8c	: all of the parameters corresponding to each behaviour of saturated soil
!				  are defined in this matrix.
!						  icpt:[1]-> sm8c (3,i) = E
!									 sm8c (7,i) = K
!									 sm8c (10,i)= Ko
!									 sm8c (11,i)= gamma
!									 sm8c (13,i)= kz, vertical water permeability (m.s^-1)
!									 sm8c (14,i) = alpha, 1-Kd/Ks
!								     sm8c (15,i) = compm, Q
!									 sm8c (16,i) = kx/kz, ratio of the horizontal water
!												   permeability to the vertical water
!												   permeability.
!						  icpt:[2]-> sm8c (1,i) = friction angle (�)
!									 sm8c (2,i) = cohesion (Pa)
!									 sm8c (3,i) = E-load, loading modulus used in non-linear
!												  elasticity (adimensionel)
!									 sm8c (4,i) = E-unld, unloading modulus used in
!												  non-linear elasticity (adim)
!									 sm8c (5,i) = n, exponent used to compute the Young
!												  loading and unloading moduli in non-linear
!												  elasticity (adim)
!									 sm8c (6,i) = Rf, coefficient multiplying the ratio
!												  (sigm1-sigm3)/(sigm1-sigm3)ult in nlin els
!												  (adim)
!                                    sm8c (7,i) = Kb, bulk modulus used in non-linear
!											      elasticity (adim)
!	                                 sm8c (8,i) = m, exponent used to compute the bulk
!												  modulus B in nlin els (adim)
!									 sm8c (9,i) = e-f, Young modulus minimal value (Pa)
!									 sm8c (10,i)=Ko, steady state coefficient (adim)
!								     sm8c (11,i)= gamma, soil volumetric weight (N,m-3)
!								     sm8c (12,i)= ten-st, traction resistency (Pa)
!									 sm8c (13,i)= perm-z
!									 sm8c (14,i)= sat
!									 sm8c (15,i)= compm
!									 sm8c (16,i)= kx/kz
!						  icpt:[3]-> sm8c (1,i) = friction angle
!									 sm8c (2,i) = cohesion
!								     sm8c (3,i) = OCR
!								     sm8c (4,i) = k
!                                    sm8c (5,i) = n, exponent used to compute the Young
!												  loading and unloading moduli in non-linear
!												  elasticity (adim)
!									 sm8c (6,i) = e
!									 sm8c (7,i) = Kb, bulk modulus used in non-linear
!												  elasticity (adim)
!									 sm8c (8,i) = m, exponent used to compute the bulk
!												  modulus B in nlin els (adim)
!									 sm8c (10,i) =Ko
!									 sm8d (11,i) =gamma
!									 sm8c (13,i) = perm-z
!									 sm8c (14,i) = sat
!									 sm8c (15,i) = compm
!									 sm8c (16,i) = kx/kz
!						  icpt:[4]-> ??????
!		sm8u	: all of the parameters corresponding to each behaviour of dry soil are
!				  defined in this matrix.
!						  icpt:[1]-> sm8u (3,i) = E
!									 sm8u (7,i) = Nou
!									 sm8u (10,i)= Ko
!									 sm8u (11,i)= gamma
!									 sm8u (13,i)= vertical water permeability Kw
!									 sm8u (14,i)= vertical gas permeability Kg
!									 sm8u (14,i) = sat
!								     sm8u (15,i) = compm
!									 sm8u (16,i) = kx/kz
!		ie3d3	: ie3d3(1:nnpbed3): nodes' number which form the current D3 BEM area
!										   Numbering direction for semi-infinite domain is
!										   ANTICLOCKWISE
!				  ie3d3(nnpbed3+1): it shows the type of just one material which
!											  is used:
!							[1]: dry boundary elm
!							[2]: saturated boundary elm
!							[3]: unsaturated boundary elm
!				  ie3d3(nnpbed3+2): it shows the number (num�ro) of (dry,sat,unsat)
!										  material which is used [m8d,m8c,m8u]
!		kbem		: number of degree of freedom per node in BE zone;
!								if IBEM=[1] -> KBEM=[2] [X,Y]
!								if IBEM=[2] -> KBEM=[3] [X,Y,Pw]
!								if IBEM=[3] -> KBEM=[4] [X,Y,Pw,Pa]
!		npbmax	: maximum nodes in all of the BE zone + 1; for example if we consider D3 BE
!				  area NPBMAX = nnpbed3 + 1
!		n12gh	: in geomechanical problems there is more than one unknown per node, so the
!				  matrices in BEM are assambled according to the doF, rather than node
!				  number. for dry soil each node has [2] doF (ux,uy), for for saturated soil
!				  each node has [3] (ux,uy,pw) doF & for unsaturated soil each node has
!				  [4] doF(ux,uy,pw,pa). then total number of doF in BE=(npbmax-1)*kbem
!		n3gh	: in semi-infinite problem in BEM we have to consider a Fictitious enclosing
!				  surface to perform point collocation method. N3GH shows the total number
!				  of doF in (D3+enclosing surface) in BE technique=(npbmax-1+nnpen)*kbem
!		nnpbed3		: total number of nodes of D3 (semi-infinite) BEM area
!		n			: total number of nodes of one D3 (semi-infinite) area
!		nn			: total number of nodes of one D3 (semi-infinite) area & enclosed area
!		n1			: total number of doF of one D3 (semi-infinite) area; N1 = KBEM * N
!					  which presents the collocation points
!		n2			: total number of doF of one D3 (semi-infinite) area & enclosed area;
!					  N2 = KBEM * NN which presents the doF of nodes in D3
!		itime		: time increment; [1]: first time increment, [>1]: higher time increment
!		mtype		: type (behaviour) of materials. in this study, only the linear
!					  behaviour is considered then: [1]linear,
!		zrho		: density of body
!		zmuy		: Lam�'s coefficient which is equal to 3EK/(9K-E)
!		zlamda		: Lam�'s coefficient which is equal to K-2MU/3
!		zrhof		: density of fluid (water) in a saturated volume
!		zk			: water permeability in z-direction
!		zm			: 1/Q, in saturated formulation
!		zalpha		: alpha, in saturated formulation
!		init		: code describing the type of initial conditions;
!					            [0] if the results are taken from a previous analysis
!								[1] if initial normal stresses are considered equal to
!									the weight of the upper soil layers
!								[2] if initial stresses are computed for the steady state
!									of the system
!								[3] if initial stresses are calculated as the weight of
!									the upper soil layers, considering that the system has
!									a mere parallelogram shape, and that the x-direction of
!									the space-frame corresponds to the physical horizontal
!				                [4] if initial stresses are given on some elements
!									(or on whole elements)
!		einit		: initial E
!		binit		: initial bulk modulus K
!		c1			: Vp, speed of primary(P) wave
!		c2			: Vs, speed of secondary(S) wave
!		xd			: x-component of nodes forming the D3 zone from line 1:nnpbed3
!					  & nodes forming the enclosed zone from line nnpbed3+1:nnpen
!		yd			: y-component of nodes forming the D3 zone from line 1:nnpbed3
!					  & nodes forming the enclosed zone from line nnpbed3+1:nnpen
!		xp			: x-component of collocation point (loading point)
!		yp			: y-component of collocation point (loading point)
!		xelm(1:3)	: x-component of nodes forming the each element in one D3
!					  (semi-infinite) zone. in GHMATD3 subroutine in the first party to
!					  calculate the diagonal composents of fundamnetal solution with rigid
!					  body motion we consider NNPBED3 & NNPEN. but after obtaining final GD3
!					  & HD3 we consider only NNBED3 nodes
!		yelm(1:3)	: y-component of nodes forming the each element in one D3
!					  (semi-infinite) zone. in GHMATD3 subroutine in the first party to
!					  calculate the diagonal composents of fundamnetal solution with rigid
!					  body motion we consider NNPBED3 & NNPEN. but after obtaining final GD3
!					  & HD3 we consider only NNBED3 nodes
!		hd3			: total H matrix in D3 (semi-infinite) area concerning all of the h
!					  matrices of each element (ehd3). global coefficient h matrices
!					  assembled by gathering element contributions.
!		gd3			: total G matrix in D3 (semi-infinite) area concerning all of the g
!					  matrices of each element (egd3). global coefficient g matrices
!					  assembled by gathering element contributions.
!		gdinv		: inverse matrice of gd3
!		gid3		: inverse matrice of gd3, G3^(-1)
!		gihd3		: gid3*hd3, G3^(-1)*H3
!		xgid3		: M*gid3, M*G3^(-1)
!		xgihd3		: M*gid1*hd1, M*G1^(-1)*H3
!		ed3		: ??????????????????????????????????????????
!		eItn		: ??????????????????????????
!		hstd3		: total H matrix in D3 (semi-infinite) area for static case concerning
!					  all of the h matrices of each element (ehstd3)
!		gstd3		: total G matrix in D3 (semi-infinite) area for static case concerning
!					  all of the g matrices of each element (egstd3)
!
!
!	INPUT	: itm,x,sm8d,sm8c,sm8u,ie3d3,hd3,kbem
!
!	OUTPUT	: hd3,gd3,gid3,gihd3,xgid3,xgihd3,ed3
!
!
!****************************************************************************************
!
	subroutine ghmatd3(itm,x,sm8d,sm8c,sm8u,ie3d3,hd3,gd3,gid3,gihd3,xgid3,xgihd3,ed3)
!
	implicit double precision (a-h,o-z)
!	double precision gdinv(n12gh,n12gh)
!
	common /a/ nnp,ifem,ibem,idyn,ieaq,nr
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /a3/ ne8d1,ne8c1,ne8u1,nload1,nbcx1,nbcy1,nbcw1,nbca1,m8d1,m8c1,m8u1,nbed11,&
				nbed21,nbed31,ne8dd,ne8cc,ne8uu
	common /a4/ xencl(500),yencl(500)
	common /b/ m8d,m8c,m8u,m8dep(10),m8cep(10),m8uep(10)
	common /c/ nload,ntime,itmax,itime,iter,iload,time,dtime,dtimei,dtimem,xinc
	common /e/ gammaw,gammas,gammaa,grav,atmp
	common /g/ init,iprint
	common /mate/zlamda,zmuy,znuu,znu,zalpha,zrho,zrhof,zm,zk,zr,zn
	common /umate/ d(3,3),dfs(3),xn,xkw,cww,cwg,xka,cgg,cgw,ev,sat,g1,rs,rw,raa,rmix,xm1,&
				   bt,et,xlambda,xmu,henry
	common /jab/ajab,bjab,cjab,djab,xelm(3),yelm(3),xp,yp
	common /Kcor/ kcorf,mtime,rmtol,c1,c2
	common /ocqm/ EPScqm,taw,Icqm,Lcqm,NPcqm
	common /unbem/ nelmu,ie3u(3,50,5)
!
	dimension x(2,nnp),sm8d(60,m8d1),sm8c(60,m8c1),sm8u(60,m8u1),ie3d3(npbmax),&
			  hd3(n12gh,n3gh),gd3(n12gh,n3gh),gid3(n12gh,n12gh),gihd3(n12gh,n12gh),&
			  xgid3(n12gh,n12gh),xgihd3(n12gh,n12gh),hstd3(n12gh,n3gh),gstd3(n12gh,n3gh),&
			  ehd3(kbem,3*kbem),egd3(kbem,3*kbem),ehstd3(kbem,3*kbem),egstd3(kbem,3*kbem),&
			  xmtf(n12gh,n12gh),xd(npbmax+nnpen),yd(npbmax+nnpen)
	dimension eItn(kbem,3*kbem),ed3(4,n3gh) ! dimension ed3????????????????
!
!
!	------------------- initializing -------------------
!
	n=nnpbed3
	nn=nnpbed3+nnpen
	n1=kbem*n
	n2=kbem*nn
	nnn=nn-1
	Iint=0
!
	if (itime.eq.1) goto 30
	nnn=n-2		!????????????????
30	continue
!
	do i=1,n1	! i shows the collocation point which can be one of the node in D3 zone
!				  for each iteration
		do j=1,n2	! j shows the field point which can be one of the node in D3 zone or
!					  in enclosed zone for each iteration
			hd3(i,j)=0.d0
			gd3(i,j)=0.d0
			ed3(1,j)=0.d0
			ed3(2,j)=0.d0
			ed3(3,j)=0.d0
			ed3(4,j)=0.d0
		enddo
	enddo
!
	if (itime.eq.1.OR.itm.eq.1) goto 55	! itime -> CONTROL, in first time step itm=1
!																in BEM
	call ghsaved(itm,n1,gd3,hd3,'ghd3')
	goto 2000
!
55	if (itime.gt.1) goto 70
	do i=1,n1
		do j=1,n2
			hstd3(i,j)=0.d0
			gstd3(i,j)=0.d0
		enddo
	enddo
70	continue
!
!
!	--------- material & geometrical properties ---------
!
	mtype=ie3d3(n+1)
	matparam: SELECT CASE (kbem)
				CASE (2)	! dry mat
					if (mtype.eq.1) then ! mtype=1:linear
						zrho=sm8d(11,mtype)/grav
						zmuy=3*sm8d(3,mtype)*sm8d(7,mtype)/(9*sm8d(7,mtype)-sm8d(3,mtype))
						zlamda=sm8d(7,mtype)-2.d0*zmuy/3.d0
					endif
					if (init.ne.0) then	! to calculate initial stresses, we are interested
!									to have the Lam� coefficients matrix. we want to
!									calculate the weight of soil, for that we transform
!									the weight of soil to the overloads F->{u}=inv[D]{F}
!									then {u}->{eps} and {sigma}=[D]*inv[D]*{f}*...
!									then in this stage, it is not very important the value
!									of Lam� coefficients -> [D]*inv[D]=I
						einit=48.d15
						binit=40.d15
						zmuy=3*einit*binit/(9*binit-einit)
						zlamda=binit-2.d0*zmuy/3.d0
					endif
!					if (mtype.eq.2) then  ! mtype=2:non-linear
!						call dmatD( )
!					endif
				CASE (3) 	! saturated mat, ie3d3(n+1)=2 -> mtype=1:linear
					if (mtype.eq.1) then ! mtype=1:linear
						zrho=sm8c(11,mtype)/grav
						zmuy=3*sm8c(3,mtype)*sm8c(7,mtype)/(9*sm8c(7,mtype)-sm8c(3,mtype))
						zlamda=sm8c(7,mtype)-2.d0*zmuy/3.d0
						zrhof=gammaw/grav
						zk=sm8c(13,mtype)/gammaw
						zm=1.d0/sm8c(15,mtype)
						zalpha=sm8c(14,mtype)
					endif
!					if (mtype.eq.2) then  ! mtype=2:non-linear
!						call dmatC( )
!					endif
!				CASE (4) 	! unsaturated mat, ie3d3(n+1)=3 -> mtype=1:linear
!					mtype=m8uep(ie3d3(n+2))
!					if (mtype.eq.1) then ! mtype=1:linear
!						call dmatU (d,xlambda,xmu,dfs,sm8u,xn,xkw,cww,cwg,xka,cgg,cgw,ev,&
!									sat,g1,rs,rw,ra,rmix,xm1,bt,et)
!					end
			  end SELECT matparam
!
!	----------------- for K-correction -----------------
!
	c1=dsqrt ((zlamda+2*zmuy)/zrho)
	c2=dsqrt (zmuy/zrho)
!
!	--- coordinate point of 3 noded boundary elements ---
!
	do ll=1,n
		ii=ie3d3(ll)
		xd(ll)=x(1,ii)
		yd(ll)=x(2,ii)
	enddo
	do ll=1,nnpen
		lll=n+ll
		xd(lll)=xencl(ll)
		yd(lll)=yencl(ll)
	enddo
	xd(nn+1)=xd(1)
	yd(nn+1)=yd(1)
!
!
!	--- selecting collocation points & boundary elements ---
!
	do 400 ll=1,n
		do 400 i=1,nnn,2
			xp=xd(ll)
			yp=yd(ll)
			xelm(1)=xd(i)
			yelm(1)=yd(i)
			xelm(2)=xd(i+1)
			yelm(2)=yd(i+1)
			xelm(3)=xd(i+2)
			yelm(3)=yd(i+2)
!
!
!	------ investigating the occurance of real and/or apparent kind of singularity ------
!	----------------- in the element and also the number of subelements -----------------
!
			call sngtp(nn,ll,i,Iint,nodo,nsub,nptg)
!
!
!	------ calculating egd3 & ehd3 in each of the elements ------
!
			call extineq(itm,itime,dtime,idyn,nodo,nsub,nptg,ehd3,egd3,ehstd3,egstd3,&
						 eItn)
!
!
!	------ assembling the ehd3 & egd3 matrices into the hd3 & gd3 ones ------
!
!
			ii=3*kbem
300			do 390 k=1,kbem
				do 380 j=1,ii
					if (itime.gt.1) goto 370
					if (i.ne.(nn-1)) goto 370	! if i=(nn-1) we consider the last elm.
					if (j.le.(2*kbem)) goto 370 	! if j>(2*kbem) we consider the last node
!												  (N1=1) of the last elm.
360					if (idyn.eq.0) goto 365
					hd3(kbem*(ll-1)+k,j-2*kbem)=hd3(kbem*(ll-1)+k,j-2*kbem)+ehd3(k,j)
					gd3(kbem*(ll-1)+k,j-2*kbem)=gd3(kbem*(ll-1)+k,j-2*kbem)+egd3(k,j)
365					if (itime.gt.1) goto 380
					hstd3(kbem*(ll-1)+k,j-2*kbem)=hstd3(kbem*(ll-1)+k,j-2*kbem)+ehstd3(k,j)
					gstd3(kbem*(ll-1)+k,j-2*kbem)=gstd3(kbem*(ll-1)+k,j-2*kbem)+egstd3(k,j)
					goto 380
!
370					if (idyn.eq.0) goto 375
					hd3(kbem*(ll-1)+k,kbem*(i-1)+j)=hd3(kbem*(ll-1)+k,kbem*(i-1)+j)+ehd3(k,j)
					gd3(kbem*(ll-1)+k,kbem*(i-1)+j)=gd3(kbem*(ll-1)+k,kbem*(i-1)+j)+egd3(k,j)
!					if (itime.eq.2) then
					ed3(k,kbem*(i-1)+j)=ed3(k,kbem*(i-1)+j)+eItn(k,j)
!					endif
!
375					if (itime.gt.1) goto 380
					hstd3(kbem*(ll-1)+k,kbem*(i-1)+j)=hstd3(kbem*(ll-1)+k,kbem*(i-1)+j)+&
													  ehstd3(k,j)
					gstd3(kbem*(ll-1)+k,kbem*(i-1)+j)=gstd3(kbem*(ll-1)+k,kbem*(i-1)+j)+&
													  egstd3(k,j)
380				continue
390			continue
400	continue	! global coefficient matrices G & H are assembled by gathering element
!				  contributions
!
!	----------------------------------------------------------------
!
	if (itime.gt.1) then
		call ghsaved(itm,n1,gd3,hd3,'ghd3')
		goto 2000
	endif
!
!
!	------ calculating diagonal coefficients of hstd3 (static case) matrice ------
!	------- by MODifIED RIGID BODY MOTION method (ahmad & banerjee, 1988) --------
!
	do i=1,n
		if (kbem.eq.2) then
			i1=2*i-1
			i2=2*i
		else if (kbem.eq.3) then
			i1=3*i-2
			i2=3*i-1
			i3=3*i
		else
			i1=4*i-3
			i2=4*i-2
			i3=4*i-1
			i4=4*i
		endif
!
!	the component Hij (for displacement) are strongly singular in all of the cases (dry,
!	saturated & unsaturated cases). then:
		hstd3(i1,i1)=0.d0
        hstd3(i2,i1)=0.d0
        hstd3(i1,i2)=0.d0
        hstd3(i2,i2)=0.d0
!
!	in saturated case, the one strongly singular fundamental solution is H33. H3j is regular
!	while Hi3 is weakly singular (which can be calculated by Gauss quadrature). then:
		if (kbem.eq.3) then
			hstd3(i3,i3)=0.d0
		endif
!
!	in unsaturated case, ...
		if (kbem.eq.4) then 	!???????????
			hstd3(i3,i3)=0.d0
			hstd3(i4,i4)=0.d0
		endif
!
!
		do 1600 j=1,nn
			if (i.eq.j) goto 1600
			if (kbem.eq.2) then
				j1=2*j-1
				j2=2*j
			else if (kbem.eq.3) then
				j1=3*j-2
				j2=3*j-1
				j3=3*j
			else
				j1=4*j-3
				j2=4*j-2
				j3=4*j-1
				j4=4*j
			endif
			hstd3(i1,i1)=hstd3(i1,i1)-hstd3(i1,j1)
			hstd3(i2,i1)=hstd3(i2,i1)-hstd3(i2,j1)
			hstd3(i1,i2)=hstd3(i1,i2)-hstd3(i1,j2)
			hstd3(i2,i2)=hstd3(i2,i2)-hstd3(i2,j2)
			if (kbem.eq.3) then
				hstd3(i3,i3)=hstd3(i3,i3)-hstd3(i3,j3)
			endif
			if (kbem.eq.4) then		!???????????
				hstd3(i3,i3)=hstd3(i3,i3)-hstd3(i3,j3)
				hstd3(i4,i4)=hstd3(i4,i4)-hstd3(i4,j4)
			endif
1600	continue
	enddo
!
!
!	------ calculating diagonal coefficients of hd3 (dynamic case) matrice ------
!	------ by MODifIED RIGID BODY MOTION method (ahmad & banerjee, 1988) --------
!
	if (idyn.eq.0) goto 1750
	do i=1,n
		if (kbem.eq.2) then
			i1=2*i-1
			i2=2*i
		else if (kbem.eq.3) then
			i1=3*i-2
			i2=3*i-1
			i3=3*i
		else
			i1=4*i-3
			i2=4*i-2
			i3=4*i-1
			i4=4*i
		endif
		hd3(i1,i1)=hd3(i1,i1)+hstd3(i1,i1)
        hd3(i2,i1)=hd3(i2,i1)+hstd3(i2,i1)
        hd3(i1,i2)=hd3(i1,i2)+hstd3(i1,i2)
        hd3(i2,i2)=hd3(i2,i2)+hstd3(i2,i2)
		if (kbem.eq.3) then
			hd3(i3,i3)=hd3(i3,i3)+hstd3(i3,i3)
		endif
		if (kbem.eq.4) then		!???????????
			hd3(i3,i3)=hd3(i3,i3)+hstd3(i3,i3)
			hd3(i4,i4)=hd3(i4,i4)+hstd3(i4,i4)
		endif
	enddo
1750 continue
!
	call ghsaved(itm,n1,gd3,hd3,'ghd3')
!
!
!	------ equalizing hd3 & gd3 to hstd3 & gstd3 in static case ------
!
	if (idyn.gt.0) goto 1850
	do i=1,n1
		do j=1,n2
			hd3(i,j)=hstd3(i,j)
			gd3(i,j)=gstd3(i,j)
		enddo
	enddo
1850 continue
!
!
!	------ calculating inverse matrice of gd3 (gid3) for transforming nodal tractions to
!	----------------------------------- nodal forces -----------------------------------
!
!	from now on, we consider just the nnpbed3's nodes (square matrix)
	do i=1,n1
		do j=1,n1
			gid3(i,j)=gd3(i,j)
		enddo
	enddo
	call dinverse(gid3,n12gh,n1,info)
!
	if (info.gt.0) then
		write (*,*) '**GHMATD3** Matrix gd3 is singular at ',info,' row'
		stop
	endif
	if (info.eq.-1) write (*,*) '**GHMATD3** Attention: bad condition',' of matrix gd3'
	write (*,*)'               inverse of gd3 concluded'
!
!
!	---- calculating global M matrice for transforming nodal tractions to nodal forces ----
!
	call mtrfrd(n,xmtf,xd,yd)
!
!
!	---------------- calculating gihd3 & xgihd3 matrices ----------------
!
	call othermat3(n,hd3,xmtf,gid3,gihd3,xgid3,xgihd3)
	write (*,*) '               othermat 3 concluded'
!
2000 continue
	return
	end

!****************************************************************************************
!
!							           BEMTND3 subroutine
!
!	This subroutine controlls BEM of area D3 at higher time increments.
!	This subroutine is called in: BEMTN
!	This subroutine calls:		  BEMTD3,GHMATD3,BEMRD3
!
!	variables used are:
!
!		id		: degrees of freedom of each node from "1" to "mdof"; for the nodes with
!				  the boundary conditions, each node takes from "mdof+1" to "mdof+nbc"
!				  [X,Y]; which are free for the dry elm &
!				  [X,Y,Pw]; which are free for the saturated elm
!				  [X,Y,Pw,Pa]; which are free for the unsaturated elm
!		x		: coordinates of nodes (nodal data) in FE & BE areas;
!							x(1,1:nnp) : X-coord
!							x(2,1:nnp) : Y-coord
!		sm8d	: all of the parameters corresponding to each behaviour of dry soil are
!				  defined in this matrix.
!						  icpt:[1]-> sm8d (3,i) = E
!									 sm8d (7,i) = K
!									 sm8d (10,i)= Ko
!									 sm8d (11,i)= gammas
!						  icpt:[2]-> sm8d (1,i) = friction angle (�)
!									 sm8d (2,i) = cohesion (Pa)
!									 sm8d (3,i) = E-load, loading modulus used in non-linear
!												  elasticity (adimensionel)
!									 sm8d (4,i) = E-unld, unloading modulus used in
!												  non-linear elasticity (adim)
!									 sm8d (5,i) = n, exponent used to compute the Young
!												  loading and unloading moduli in non-linear
!											      elasticity (adim)
!									 sm8d (6,i) = Rf, coefficient multiplying the ratio
!												  (sigm1-sigm3)/(sigm1-sigm3)ult
!												  in nlin els (adim)
!                                    sm8d (7,i) = Kb, bulk modulus used in non-linear
!											      elasticity (adim)
!	                                 sm8d (8,i) = m, exponent used to compute the bulk
!										          modulus B in nlin els (adim)
!									 sm8d (9,i) = e-f, Young modulus minimal value (Pa)
!									 sm8d (10,i)=Ko, steady state coefficient (adim)
!								     sm8d (11,i)= gammas, soil volumetric weight (N,m-3)
!								     sm8d (12,i)= ten-st, traction resistency (Pa)
!						  icpt:[3]-> sm8d (1,i) = friction angle
!									 sm8d (2,i) = cohesion
!								     sm8d (3,i) = OCR
!								     sm8d (4,i) = k
!                                    sm8d (5,i) = n, exponent used to compute the Young
!												  loading and unloading moduli in non-linear
!												  elasticity (adim)
!									 sm8d (6,i) = e
!									 sm8d (7,i) = Kb, bulk modulus used in non-linear
!												  elasticity (adim)
!									 sm8d (8,i) = m, exponent used to compute the bulk
!												  modulus B in nlin els (adim)
!									 sm8d (10,i) =Ko
!									 sm8d (11,i) =gamma
!						  icpt:[4]-> ??????
!		sm8c	: all of the parameters corresponding to each behaviour of saturated soil
!				  are defined in this matrix.
!						  icpt:[1]-> sm8c (3,i) = E
!									 sm8c (7,i) = K
!									 sm8c (10,i)= Ko
!									 sm8c (11,i)= gamma
!									 sm8c (13,i)= kz, vertical water permeability (m.s^-1)
!									 sm8c (14,i) = alpha, 1-Kd/Ks
!								     sm8c (15,i) = compm, Q
!									 sm8c (16,i) = kx/kz, ratio of the horizontal water
!												   permeability to the vertical water
!												   permeability.
!						  icpt:[2]-> sm8c (1,i) = friction angle (�)
!									 sm8c (2,i) = cohesion (Pa)
!									 sm8c (3,i) = E-load, loading modulus used in non-linear
!												  elasticity (adimensionel)
!									 sm8c (4,i) = E-unld, unloading modulus used in
!												  non-linear elasticity (adim)
!									 sm8c (5,i) = n, exponent used to compute the Young
!												  loading and unloading moduli in non-linear
!												  elasticity (adim)
!									 sm8c (6,i) = Rf, coefficient multiplying the ratio
!												  (sigm1-sigm3)/(sigm1-sigm3)ult in nlin els
!												  (adim)
!                                    sm8c (7,i) = Kb, bulk modulus used in non-linear
!											      elasticity (adim)
!	                                 sm8c (8,i) = m, exponent used to compute the bulk
!												  modulus B in nlin els (adim)
!									 sm8c (9,i) = e-f, Young modulus minimal value (Pa)
!									 sm8c (10,i)=Ko, steady state coefficient (adim)
!								     sm8c (11,i)= gamma, soil volumetric weight (N,m-3)
!								     sm8c (12,i)= ten-st, traction resistency (Pa)
!									 sm8c (13,i)= perm-z
!									 sm8c (14,i)= sat
!									 sm8c (15,i)= compm
!									 sm8c (16,i)= kx/kz
!						  icpt:[3]-> sm8c (1,i) = friction angle
!									 sm8c (2,i) = cohesion
!								     sm8c (3,i) = OCR
!								     sm8c (4,i) = k
!                                    sm8c (5,i) = n, exponent used to compute the Young
!												  loading and unloading moduli in non-linear
!												  elasticity (adim)
!									 sm8c (6,i) = e
!									 sm8c (7,i) = Kb, bulk modulus used in non-linear
!												  elasticity (adim)
!									 sm8c (8,i) = m, exponent used to compute the bulk
!												  modulus B in nlin els (adim)
!									 sm8c (10,i) =Ko
!									 sm8d (11,i) =gamma
!									 sm8c (13,i) = perm-z
!									 sm8c (14,i) = sat
!									 sm8c (15,i) = compm
!									 sm8c (16,i) = kx/kz
!						  icpt:[4]-> ??????
!		sm8u	: all of the parameters corresponding to each behaviour of dry soil are
!				  defined in this matrix.
!						  icpt:[1]-> sm8u (3,i) = E
!									 sm8u (7,i) = K
!									 sm8u (10,i)= Ko
!									 sm8u (11,i)= gamma
!									 sm8u (13,i)= vertical water permeability Kw
!									 sm8u (14,i)= vertical gas permeability Kg
!									 sm8u (14,i) = sat
!								     sm8u (15,i) = compm
!									 sm8u (16,i) = kx/kz
!		ie3d3	: ie3d3(1:nnpbed3): nodes' number which form the current D3 BEM area
!										   Numbering direction for semi-infinite domain is
!										   ANTICLOCKWISE
!				  ie3d3(nnpbed3+1): it shows the type of just one material which
!											  is used:
!							[1]: dry boundary elm
!							[2]: saturated boundary elm
!							[3]: unsaturated boundary elm
!				  ie3d3(nnpbed3+2):
!		tbed3		: traction/flow [tx,ty,qw,qa] vector T(N) at time: "itime-1"
!		ubed3		: displacement/pressure [ux,uy,Pw,Pa] vector U(N) at time: "itime-1" which
!					  is calculated in the previous time step in SOLVE subroutine &
!					  replaced in "DISP" vector in UPDATE3 subroutine
!		rbed3		: stress tensor in D3 BEM area which includes JUST effect of the past
!					  dynamic history; Z(N)=SUM(n=1:N-1) [G(N+1-n)*T(n)-H(N+1-n)*U(n)];
!					  in the first time step Z(1)=rbed3(1,ii)=0
!		frbed3		:
!		gihd3		: gid3*hd3, G3^(-1)*H3
!		gid3		: inverse matrice of gd3, G3^(-1)
!		xgihd3		: M*gid1*hd1, M*G1^(-1)*H3
!		xgid3		: M*gid3, M*G3^(-1)
!		disp		: displacement/pressures vector which consists of initial displacement
!					  value of corresponding nodes in which their identities code of Ux & Uy
!					  are defined in ID and initial water and air values of corresponding
!					  nodes in which their identities code of Pw & Pa are defined in ID.
!					  these displacement initial values can be defined for a free node or
!					  for a node with boundary conditions and the water and air initial values
!					  must be defined for a free node; in each "itime" the displacements/
!					  pressures vector (R->DT) in PREVIOUS TIME STEP (N-1) is transfered
!					  in DISP vector
!		kcorf		: K-correction parameters: K-correction code [0,1];
!								[0]: wrobel method
!								[1]: time truncated method
!		ed3			:
!		red3		:
!		kpos		:
!
!
!	INPUT	:
!
!	OUTPUT	:
!
!
!****************************************************************************************
!
	subroutine bemtnd3(id,x,sm8d,sm8c,sm8u,ie3d3,tbed3,ubed3,rbed3,frbed3,gihd3,gid3,&
					   xgihd3,xgid3,disp,sig3u)
!
	implicit double precision (a-h,o-z)
!
	character tx*25
!
	common /a/ nnp,ifem,ibem,idyn,ieaq,nr
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /a3/ ne8d1,ne8c1,ne8u1,ne3u1,nload1,nbcx1,nbcy1,nbcw1,nbca1,m8d1,m8c1,m8u1,&
				nbed11,nbed21,nbed31,nbeint1,ne8dd,ne8cc,ne8uu
	common /a4/ xencl(500),yencl(500)
	common /b/ m8d,m8c,m8u,m8dep(10),m8cep(10),m8uep(10)
	common /c/ nload,ntime,itmax,itime,iter,iload,time,dtime,dtimei,dtimem,xinc
	common /h/ mdof,mdofn,ls,ls1,isolv
	common /stability/anew1,anew2,wil,intb
	common /Kcor/ kcorf,mtime,rmtol,c1,c2
!
	dimension id(4,nnp),x(2,nnp),sm8d(60,m8d1),sm8c(60,m8c1),sm8u(60,m8u1),ie3d3(npbmax+1),&
			  tbed3(ntime,n12gh),ubed3(ntime,n12gh),rbed3(ntime,n12gh),frbed3(n12gh),&
			  gihd3(n12gh,n12gh),gid3(n12gh,n12gh),xgihd3(n12gh,n12gh),xgid3(n12gh,n12gh),&
			  disp(mdofn),gd3(n12gh,n3gh),hd3(n12gh,n3gh),hd31(n12gh,n12gh),gd31(n12gh,n12gh),&
			  sig3u(5,500,nbed11)
	dimension ed3(4,n3gh),red3(4)
!
!
!	------------- calculate traction vector of time step itime-1 (i.e. n-1) -------------
!
	call bemtd3(id,ie3d3,tbed3,ubed3,rbed3,gihd3,gid3,disp)
!
!
!	--------------- determine the effect of previous time increments ---------------
!	------------------------------	on the current one -----------------------------
!
!	initializing vector rbed3
!
	n=nnpbed3*kbem
	do i=1,n
		rbed3(itime,i)=0.d0
	enddo
!
	red3(1)=0.d0
	red3(2)=0.d0
	red3(3)=0.d0
	red3(4)=0.d0
!
	if (kcorf.eq.1) then
!
!	----- calculating traction vectors in 2nd step -----
!
		if (itime.eq.2) then
			kpos=0
			itm=1
!	calculating g & h matices at itime=2
			call ghmatd3(itm,x,sm8d,sm8c,sm8u,ie3d3,hd3,gd3,gid3,gihd3,xgid3,xgihd3,&
						 ed3)
!	increasing rbed3 vector due to itm time increment -> itm=1 (second time increment)
			call bemrd3(itm,tbed3,ubed3,rbed3,hd3,gd3,ed3,R,kpos)
			temp=1
!	finding step limit for K-correction
			do WHILE (temp.ge.rmtol)
				lim=lim+1
				temp=ABS (R**lim)
				PRINT *,temp
			enddo
		endif
!
!	----- calculating traction vectors in other steps -----
!
		if (itime.gt.2) then
!	for steps from limit to cutting point (far history)
			if (itime.gt.(lim+mtime)) then	!LIM??????????????????????
				kpos=1
				init=itime-mtime-lim
			else
!	for steps after cutting point (near history)
				kpos=0
				init=1
			endif
			do 500 itm=init,itime-1
				if (kpos.eq.0) then
					tx='     exact'
				else if ((kpos.eq.1).and.(itm.lt.(itime-mtime))) then
					tx='     aproximated'
				else if ((kpos.eq.1).and.(itm.ge.(itime-mtime))) then
					tx='     exact'
				endif
				PRINT *,itm,tx
!	calculating g & h matices at itm time increment
				call ghmatd3(itm,x,sm8d,sm8c,sm8u,ie3d3,hd3,gd3,gid3,gihd3,xgid3,xgihd3,&
							 ed3)
!	increasing rbed3 vector due to itm time increment
				call bemrd3(itm,tbed3,ubed3,rbed3,hd3,gd3,ed3,R,kpos)
500			continue
		endif
!
	else
		kpos=0
		do itm=1,itime-1
!	calculating g & h matices at itime>1
			call ghmatd3(itm,x,sm8d,sm8c,sm8u,ie3d3,hd3,gd3,gid3,gihd3,xgid3,xgihd3,&
						 ed3)
!	increasing rbed3 vector due to itm time increment
			call bemrd3(itm,tbed3,ubed3,rbed3,hd3,gd3,ed3,R,kpos)
		enddo
	endif
!
!
!	----- last moment itime-1 -----
!
	call ghsaved(itime,n,gd31,hd31,'ghd3') ! this reads G(1) & H(1) in each itime
	do i=1,n
		do j=1,n
			rbed3(itime,i)=rbed3(itime,i)+( gd31(i,j)*tbed3(itm,j)-&
						   hd31(i,j)*ubed3(itm,j) )*(1.d0-wil)/wil
		enddo
	enddo
!
!
!	----- transfer vector rbed3 to its corresponding force vector form -----
!
	call bemfd3(frbed3,rbed3,xgid3)
!
	return
	end
!****************************************************************************************
!
!							           BEMTD3 subroutine
!
!	This subroutine calculates the tractions of area D3 at time increment itime-1.
!	this subroutine is called in this section: BEMTND3
!
!	variables used are:
!
!		id		: degrees of freedom of each node from "1" to "mdof"; for the nodes with
!				  the boundary conditions, each node takes from "mdof+1" to "mdof+nbc"
!				  [X,Y]; which are free for the dry elm &
!				  [X,Y,Pw]; which are free for the saturated elm
!				  [X,Y,Pw,Pa]; which are free for the unsaturated elm
!		ie3d3	: ie3d3(1:nnpbed3): nodes' number which form the current D3 BEM area
!										   Numbering direction for semi-infinite domain is
!										   ANTICLOCKWISE
!				  ie3d3(nnpbed3+1): it shows the type of just one material which
!											  is used:
!							[1]: dry boundary elm
!							[2]: saturated boundary elm
!							[3]: unsaturated boundary elm
!				  ie3d3(nnpbed3+2):
!		tbed3	: traction/flow [tx,ty,qw,qa] vector T(N) at time: "itime-1"
!		ubed3	: displacement/pressure [ux,uy,Pw,Pa] vector U(N) at time: "itime-1" which
!				  is calculated in the previous time step in SOLVE subroutine &
!				  replaced in "DISP" vector in UPDATE3 subroutine
!		rbed3	: stress tensor in D3 BEM area which includes JUST effect of the past
!				  dynamic history; Z(N)=SUM(n=1:N-1) [G(N+1-n)*T(n)-H(N+1-n)*U(n)];
!				  in the first time step Z(1)=rbed3(1,ii)=0
!		gihd3	: gid3*hd3, G3^(-1)*H3
!		gid3	: inverse matrice of gd3, G3^(-1)
!		kbem	: number of degree of freedom per node in BE zone;
!							if IBEM=[1] -> KBEM=[2] [X,Y]
!							if IBEM=[2] -> KBEM=[3] [X,Y,Pw]
!							if IBEM=[3] -> KBEM=[4] [X,Y,Pw,Pa]
!		disp	: displacement/pressures vector which consists of initial displacement
!				  value of corresponding nodes in which their identities code of Ux & Uy
!				  are defined in ID and initial water and air values of corresponding
!				  nodes in which their identities code of Pw & Pa are defined in ID.
!				  these displacement initial values can be defined for a free node or
!				  for a node with boundary conditions and the water and air initial values
!				  must be defined for a free node; in each "itime" the displacements/
!				  pressures vector (R->DT) in PREVIOUS TIME STEP (N-1) is transfered
!				  in DISP vector
!		ieaq	: earthquake code [0,1,2,3];
!							[0]: static case
!							[1]: accelaration
!							[2]: imposed deplacement
!							[3]: incident wave, in most of the time IEAQ=[3]
!		nnpbed3	: number of nodes of one D3 BEM area
!		ueaq	: if TETA-method is used, displacement field due to incident waves transforms
!				  to "Ueq(N)-[(teta-1)/teta]*Uex(N-1)"
!				  in this subroutine, vector of U_inc at time: "itime-1"
!		ubed3	: displacement/pressure [ux,uy,Pw,Pa] vector U(N) at time: "itime-1" which
!				  is calculated in the previous time step in SOLVE subroutine &
!				  replaced in "DISP" vector in UPDATE3 subroutine
!		disp	: displacement/pressures vector which consists of initial displacement
!				  value of corresponding nodes in which their identities code of Ux & Uy
!				  are defined in ID and initial water and air values of corresponding
!				  nodes in which their identities code of Pw & Pa are defined in ID.
!				  these displacement initial values can be defined for a free node or
!				  for a node with boundary conditions and the water and air initial values
!				  must be defined for a free node; in each "itime" the displacements/
!				  pressures vector (R->DT) in PREVIOUS TIME STEP (N-1) is transfered
!				  in DISP vector
!		uxeq	: (incident wave) displacement in x dir
!							uxeq (3,1:nnpbed3): in current time -> itime
!							uxeq (2,1:nnpbed3): in previous time -> itime-1
!							uxeq (1,1:nnpbed3): in two previous time -> itime-2
!		uyeq	: (incident wave) displacement in y dir
!							uyeq (3,1:nnpbed3): in current time -> itime
!							uyeq (2,1:nnpbed3): in previous time -> itime-1
!							uyeq (1,1:nnpbed3): in two previous time -> itime-2
!
!
!
!	INPUT	:
!
!	OUTPUT	: ubed3,tbed3
!
!
!****************************************************************************************
!
	subroutine bemtd3(id,ie3d3,tbed3,ubed3,rbed3,gihd3,gid3,disp)
!
	implicit double precision (a-h,o-z)
!
	common /a/ nnp,ifem,ibem,idyn,ieaq,nr
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /c/ nload,ntime,itmax,itime,iter,iload,time,dtime,dtimei,dtimem,xinc
	common /g/ init,iprint
	common /h/ mdof,mdofn,ls,ls1,isolv
	common /chrg3/ nstime,neaq,nneaq(500),uxeaq,uyeaq,axeaq,ayeaq,uxeq(3,500),uyeq(3,500)
	common /stability/anew1,anew2,wil,intb
!
	dimension id(4,nnp),ie3d3(npbmax),tbed3(ntime,n12gh),ubed3(ntime,n12gh),&
			  rbed3(ntime,n12gh),gihd3(n12gh,n12gh),gid3(n12gh,n12gh),disp(mdofn),ueaq(n12gh)
!
!
	if (ieaq.eq.3) then
		do i=1,nnpbed3
			ueaq(kbem*(i-1)+1)=uxeq(2,i) -(wil-1.d0)/wil*uxeq(1,i)
			ueaq(kbem*(i-1)+2)=uyeq(2,i) -(wil-1.d0)/wil*uyeq(1,i)
			if (kbem.eq.3) ueaq(kbem*i)=0.d0
			if (kbem.eq.4) then
				ueaq(kbem*i-1)=0.d0
				ueaq(kbem*i)  =0.d0
			endif
		enddo
	endif
!
!	displacement/pressure vector ubed3 at time: itime-1
	it=itime-1
	do k=1,nnpbed3
		j=ie3d3(k)
		do i=1,kbem
			ubed3(it,kbem*k-kbem+i)=disp(id(i,j))
		enddo
	enddo
!
!	traction/flow vector tbed3
	n=nnpbed3*kbem
!	if (ifem.eq.0.OR.init.eq.0) then
!		do i=1,n
!			tbed3(it,i)=0.d0
!		enddo
!	else
		do i=1,n
			tbed3(it,i)=0.d0
			do j=1,n
				tbed3(it,i) = tbed3(it,i) + gihd3(i,j)*ubed3(it,j) - gid3(i,j)*rbed3(it,j)
				if (ieaq.eq.3) tbed3(it,i) = tbed3(it,i) - gid3(i,j)*ueaq(j)
			enddo
		enddo
!	endif
!
	return
	end
!****************************************************************************************
!
!							           BEMRD3 subroutine
!
!	This subroutine calculates the rbed3 vector's increment due to time incement itm.
!	this subroutine is called in this section: BEMTND3
!
!	variables used are:
!
!		tbed3		: traction/flow [tx,ty,qw,qa] vector T(N) at time: "itime-1"
!		ubed3		: displacement/pressure [ux,uy,Pw,Pa] vector U(N) at time: "itime-1" which
!					  is calculated in the previous time step in SOLVE subroutine &
!					  replaced in "DISP" vector in UPDATE3 subroutine
!		rbed3		: stress tensor in D3 BEM area which includes JUST effect of the past
!					  dynamic history; Z(N)=SUM(n=1:N-1) [G(N+1-n)*T(n)-H(N+1-n)*U(n)];
!					  in the first time step Z(1)=rbed3(1,ii)=0
!		hd3			: total H matrix in D3 (semi-infinite) area concerning all of the h
!					  matrices of each element (ehd3). global coefficient h matrices
!					  assembled by gathering element contributions.
!		gd3			: total G matrix in D3 (semi-infinite) area concerning all of the g
!					  matrices of each element (egd3). global coefficient g matrices
!					  assembled by gathering element contributions.
!		kbem		: number of degree of freedom per node in BE zone;
!								if IBEM=[1] -> KBEM=[2] [X,Y]
!								if IBEM=[2] -> KBEM=[3] [X,Y,Pw]
!								if IBEM=[3] -> KBEM=[4] [X,Y,Pw,Pa]
!		ed3			: ??????????????????????????????????????????
!		kpos		:
!		R			:
!		red3		:
!		kcorf		: K-correction parameters: K-correction code [0,1];
!								[0]: wrobel method
!								[1]: time truncated method
!
!
!
!	INPUT	:
!
!	OUTPUT	:
!
!
!****************************************************************************************
!
	subroutine bemrd3(itm,tbed3,ubed3,rbed3,hd3,gd3,ed3,R,kpos)
!
	implicit double precision (a-h,o-z)
!
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /c/ nload,ntime,itmax,itime,iter,iload,time,dtime,dtimei,dtimem,xinc
	common /stability/anew1,anew2,wil,intb
	common /Kcor/ kcorf,mtime,rmtol,c1,c2
!
	dimension tbed3(ntime,n12gh),ubed3(ntime,n12gh),rbed3(ntime,n12gh),hd3(n12gh,n3gh),&
			  gd3(n12gh,n3gh)
	dimension ed3(4,n3gh),red3(4)	!???????
!
!
	red3(1)=0.d0
	red3(2)=0.d0
	red3(3)=0.d0
	red3(4)=0.d0
!
!	Solution whitout K-correction
!
	if (kcorf.eq.0) then
		n=kbem*nnpbed3
		do i=1,n
			do j=1,n
				rbed3(itime,i)=rbed3(itime,i)+(gd3(i,j)*tbed3(itm,j)-&
											   hd3(i,j)*ubed3(itm,j))/wil
			enddo
		enddo
	else
!	Solution with K-correction
!	exact solution for the times before reaching to cuuting point
		if (kpos.eq.0) then
			n=kbem*nnpbed3
			do i=1,n
				do j=1,n
					rbed3(itime,i)=rbed3(itime,i)+(gd3(i,j)*tbed3(itm,j)-&
												   hd3(i,j)*ubed3(itm,j))/wil
				enddo
			enddo
		endif
!	aproximated solution for far history
		if ((kpos.eq.1).and.(itm.lt.(itime-mtime))) then
			n=kbem*nnpbed3
			do i=1,n
				do j=1,n
					coef=R**(itime-mtime-itm-1)		! R ????????????????????????
					rbed3(itime,i)=rbed3(itime,i)+(gd3(i,j)*tbed3(itm,j)-&
												   hd3(i,j)*ubed3(itm,j))*coef/wil
				enddo
			enddo
		endif
!	exact solution for near history
		if ((kpos.eq.1).and.(itm.ge.(itime-mtime))) then
			n=kbem*nnpbed3
			do i=1,n
				do j=1,n
					rbed3(itime,i)=rbed3(itime,i)+(gd3(i,j)*tbed3(itm,j)-&
												   hd3(i,j)*ubed3(itm,j))/wil
				enddo
			enddo
		endif
!	calculating the K-correction final factor
		if (itime.eq.2) then	!????????????????????????????????
			do 400 k=1,nnpbed3
				idf1=2*k-1
				idf2=2*k
				red3(1)=red3(1)+ed3(1,idf1)
				red3(2)=red3(2)+ed3(1,idf2)
				red3(3)=red3(3)+ed3(2,idf1)
				red3(4)=red3(4)+ed3(2,idf2)
400			continue
			R1=red3(1)+red3(2)
			R2=red3(3)+red3(4)
			R=R2/R1
		endif
	endif
!
	return
	end
!****************************************************************************************
!
!							           BEMFD3 subroutine
!
!	This subroutine converts the area D3's RBED3 (traction) vector to a (M/G1)*RBED3 (force)
!	one.
!	FRBED3=(M/G1)*{[SUN(n=1:N-1)[G(N-n+1)*T(n)-H(N-n+1)*U(n)]]/wil+
!								(1-wil/wil)*[G(1)*T(N-1)-H(1)*U(N-1)]}
!	This subroutine is called in: BEMTND3
!	This subroutine calls: -
!
!	variables used are:
!
!		rbed3		: stress tensor in D3 BEM area which includes JUST effect of the past
!					  dynamic history; Z(N)=SUM(n=1:N-1) [G(N+1-n)*T(n)-H(N+1-n)*U(n)];
!					  in the first time step Z(1)=rbed3(1,ii)=0
!		frbed3		: [M/G(1)]*RBED3(N,:)
!		kbem		: number of degree of freedom per node in BE zone;
!								if IBEM=[1] -> KBEM=[2] [X,Y]
!								if IBEM=[2] -> KBEM=[3] [X,Y,Pw]
!								if IBEM=[3] -> KBEM=[4] [X,Y,Pw,Pa]
!		nnpbed3		: number of nodes of one D3 BEM area
!
!
!	INPUT	:
!
!	OUTPUT	:
!
!
!****************************************************************************************
!
	subroutine bemfd3(frbed3,rbed3,xgid3)
!
	implicit double precision (a-h,o-z)
!
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /c/ nload,ntime,itmax,itime,iter,iload,time,dtime,dtimei,dtimem,xinc
!
	dimension frbed3(n12gh),rbed3(ntime,n12gh),xgid3(n12gh,n12gh)
!
!
	n=kbem*nnpbed3
	do i=1,n
		frbed3(i)=0.d0
		do j=1,n
			frbed3(i) = frbed3(i) + xgid3(i,j)*rbed3(itime,j)
		enddo
	enddo
!
	return
	end
!****************************************************************************************
!
!							           BEMFD1 subroutine
!
!	This subroutine converts the area D1's Z (traction) vector to a Z' (force) one.
!	Z'=(M/G1)*Z
!	this subroutine is called in this section: BEMTND1
!
!	variables used are:
!
!		rbed1		: stress tensor in D1 BEM area which includes JUST effect of the past
!					  dynamic history; Z(N)=SUM(n=1:N-1) [G(N+1-n)*T(n)-H(N+1-n)*U(n)];
!					  in the first time step Z(1)=rbed1(1,ii)=0
!		frbed1		:
!
!
!
!	INPUT	:
!
!	OUTPUT	:
!
!
!****************************************************************************************
!
	subroutine bemfd1(k1,frbed1,rbed1,xgid1)
!
	implicit double precision (a-h,o-z)
!
	common /a/ nnp,ifem,ibem,idyn,ieaq,nr
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /a3/ ne8d1,ne8c1,ne8u1,ne3u1,nload1,nbcx1,nbcy1,nbcw1,nbca1,m8d1,m8c1,m8u1,&
				nbed11,nbed21,nbed31,nbeint1,ne8dd,ne8cc,ne8uu
	common /c/ nload,ntime,itmax,itime,iter,iload,time,dtime,dtimei,dtimem,xinc
!
	dimension frbed1(n12gh,nbed11),rbed1(ntime,n12gh,nbed11),&
			  xgid1(n12gh,n12gh,nbed11)

!
!
	n=kbem*nnpbed1(k1)
	do i=1,n
		frbed1(i,k1)=0.d0
		do j=1,n
			frbed1(i,k1) = frbed1(i,k1) + xgid1(i,j,k1)*rbed1(itime,j,k1)
		enddo
	enddo
!
!
	return
	end
!****************************************************************************************
!
!							           BEMINIT subroutine
!
!	This subroutine initializes the r (stresses) & f (nodal forces) vectors
!	due to first time increment to zero.
!	This subroutine is called in: BEM
!	This subroutine calls: -
!
!	variables used are:
!
!		frbed1	: nodal forces in D1 BEM area (finite zone)
!		frbed2	: nodal forces in D2 BEM area (infinite zone)
!		frbed3	: nodal forces in D3 BEM area (semi-infinite zone)
!		rbed1	: stress tensor in D1 BEM area which includes JUST effect of the past
!				  dynamic history; Z(N)=SUM(n=1:N-1) [G(N+1-n)*T(n)-H(N+1-n)*U(n)];
!				  in the first time step Z(1)=rbed1(1,ii)=0
!		rbed2	: stress tensor in D2 BEM area which includes JUST effect of the past
!				  dynamic history; Z(N)=SUM(n=1:N-1) [G(N+1-n)*T(n)-H(N+1-n)*U(n)];
!				  in the first time step Z(1)=rbed2(1,ii)=0
!		rbed3	: stress tensor in D3 BEM area which includes JUST effect of the past
!				  dynamic history; Z(N)=SUM(n=1:N-1) [G(N+1-n)*T(n)-H(N+1-n)*U(n)];
!				  in the first time step Z(1)=rbed3(1,ii)=0
!		ibem	: BEM code [0,1,2,3];
!							[0]: program doesn't consider BEM analysis. except some
!								 simple examples IBEM must never be [0]
!							[1]: BEM_dry
!							[2]: BEM_saturate
!							[3]: BEM_unsaturate
!		kbem	: number of degree of freedom per node in BE zone [2,3,4];
!							if IBEM=1 -> KBEM=2 [X,Y]
!							if IBEM=2 -> KBEM=3 [X,Y,Pw]
!							if IBEM=4 -> KBEM=4 [X,Y,Pw,Pa]
!		nbed1	: number of D1 BEM areas (finite areas), it can be at maximum 5 zones
!		nbed2	: number of D2 BEM areas (infinite areas), it can be at maximum 5 zones
!		nbed3	: number of D3 BEM areas (semi-infinite area), it can be at maximum 1 zone
!		nnpbed1	: number of nodes of each D1 BEM area; NNPBED1(1 to NBED1)
!		nnpbed2	: number of nodes of each D2 BEM area; NNPBED2(1 to NBED2)
!		nnpbed3	: number of nodes of one D3 BEM area

!
!	INPUT	: nbed1,nbed2,nbed3,nnpbed1,nnpbed2,nnpbed3,kbem,ibem
!
!	OUTPUT	: rbed1,rbed2,rbed3,frbed1,frbed2,frbed3
!
!
!****************************************************************************************
!
	subroutine beminit(rbed1,rbed2,rbed3,frbed1,frbed2,frbed3)
!
	implicit double precision (a-h,o-z)
!
	common /a/ nnp,ifem,ibem,idyn,ieaq,nr
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /a3/ ne8d1,ne8c1,ne8u1,ne3u1,nload1,nbcx1,nbcy1,nbcw1,nbca1,m8d1,m8c1,m8u1,&
				nbed11,nbed21,nbed31,nbeint1,ne8dd,ne8cc,ne8uu
	common /c/ nload,ntime,itmax,itime,iter,iload,time,dtime,dtimei,dtimem,xinc
!
	dimension rbed1(ntime,n12gh,nbed11),rbed2(ntime,n12gh,nbed21),rbed3(ntime,n12gh),&
			  frbed1(n12gh,nbed11),frbed2(n12gh,nbed21),frbed3(n12gh)
!
!
!	  D1 (finite) areas
	if (nbed1.eq.0) goto 150
	do k1=1,nbed1
		do ii=1,nnpbed1(k1)*kbem
			rbed1(1,ii,k1)=0.d0
			frbed1(ii,k1)=0.d0
		enddo
	enddo
150	continue
!
!     D2 areas; IT IS not VERifIED
	if (nbed2.eq.0) goto 350
	do k1=1,nbed2
		do ii=1,nnpbed2(k1)*kbem
			rbed2(1,ii,k1)=0.d0
			frbed2(ii,k1)=0.d0
		enddo
	enddo
350	continue
!
!     D3 area (semi-infinite zone)
	if (nbed3.eq.0) goto 550
	do ii=1,nnpbed3*kbem
		rbed3(1,ii)=0.d0
		frbed3(ii)=0.d0
	enddo
550	continue
!
	return
	end
!****************************************************************************************
!
!							           BEMRD1 subroutine
!
!	This subroutine calculates the rbed1 vector's increment due to time incement itm.
!	this subroutine is called in this section: BEMTND1
!
!	variables used are:
!
!		ubed1		: displacement/pressure [ux,uy,Pw,Pa] vector at time: "itime-1" which
!					  is calculated in the previous time step
!		tbed1		: traction/flow [tx,ty,qw,qa] vector at time: "itime-1" which is
!					  calculated in the previous time step
!		rbed1		: stress tensor in D3 BEM area which includes JUST effect of the past
!					  dynamic history; Z(N)=SUM(n=1:N-1) [G(N+1-n)*T(n)-H(N+1-n)*U(n)];
!					  in the first time step Z(1)=rbed3(1,ii)=0
!		red1		:
!
!
!	INPUT	:
!
!	OUTPUT	:
!
!
!****************************************************************************************
!
	subroutine bemrd1(itm,k1,tbed1,ubed1,rbed1,hd1,gd1)
!
	implicit double precision (a-h,o-z)
!
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /a3/ ne8d1,ne8c1,ne8u1,ne3u1,nload1,nbcx1,nbcy1,nbcw1,nbca1,m8d1,m8c1,m8u1,&
				nbed11,nbed21,nbed31,nbeint1,ne8dd,ne8cc,ne8uu
	common /c/ nload,ntime,itmax,itime,iter,iload,time,dtime,dtimei,dtimem,xinc
	common /stability/anew1,anew2,wil,intb
!
	dimension tbed1(ntime,n12gh,nbed11),ubed1(ntime,n12gh,nbed11),rbed1(ntime,n12gh,nbed11),&
			  gd1(n12gh,n12gh),hd1(n12gh,n12gh)
!
!
	n=kbem*nnpbed1(k1)
	do i=1,n
		do j=1,n
			rbed1(itime,i,k1)=rbed1(itime,i,k1)+( gd1(i,j)*tbed1(itm,j,k1)-&
												 hd1(i,j)*ubed1(itm,j,k1) )/wil
		enddo
	enddo
!
	return
	end
!****************************************************************************************
!
!							           BEMS1S2 subroutine
!
!	This subroutine assembles bem stiffness (mgih) matrices and transfers them to vectors
!	s1 & s2.
!	This subroutine is called in: BEM
!	This subroutine calls: -
!
!	variables used are:
!
!		id		: degrees of freedom of each node from "1" to "mdof"; for the nodes with
!				  the boundary conditions, each node takes from "mdof+1" to "mdof+nbc"
!				  [X,Y]; which are free for the dry elm &
!				  [X,Y,Pw]; which are free for the saturated elm
!				  [X,Y,Pw,Pa]; which are free for the unsaturated elm
!		na		: this locates the diagonal coefficient for skyline storage. we economise
!				  the space by neglecting the stockage of certain coefs which are zero
!				  because of the absence of the connectivity of doF. we delete just the zeros
!				  lied to the non connectivity. therefore na(i) shows the i-th coef to stock.
!		s1		:
!		s2		:
!		xgihd1	: M*gid1*hd1
!		xgihd2	: M*gid2*hd2
!		xgihd3	: M*gid3*hd3
!		ide		: id matrix for BE zone which put in verctorial format
!		jdf		: global doF number of a node (1:nnpbed1) with the local dof number kbem
!		ndfe	: number of degree of freedom for BE zone
!		ls		: lenght of the global rigidity matrix in symmetric case. if the matrix
!				  is not symmetric, 2 triangles must be filled.
!		ls1		: if the problem is non symmetric then this partie must be taken into
!				  account then for this partie there is the same number of coefs to
!				  stocke like ls.
!		mdof	: maximum degree of freedom by considering the boundary conditions
!		isolve	: solution type indicator [1,3];
!							[1]: for a symmetric computation
!							[3]: for a non-symmetric computation
!
!	INPUT	:
!
!	OUTPUT	:
!
!
!****************************************************************************************
!
	subroutine bems1s2(id,na,s1,s2,ie3d1,ie3d2,ie3d3,xgihd1,xgihd2,xgihd3,r)
!
	implicit double precision (a-h,o-z)
!
	common /a/ nnp,ifem,ibem,idyn,ieaq,nr
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /a3/ ne8d1,ne8c1,ne8u1,ne3u1,nload1,nbcx1,nbcy1,nbcw1,nbca1,m8d1,m8c1,m8u1,&
				nbed11,nbed21,nbed31,nbeint1,ne8dd,ne8cc,ne8uu
	common /b/ m8d,m8c,m8u,m8dep(10),m8cep(10),m8uep(10)
	common /c/ nload,ntime,itmax,itime,iter,iload,time,dtime,dtimei,dtimem,xinc
	common /h/ mdof,mdofn,ls,ls1,isolv
	common /chrg2/ npw,npa,nnpw(200),nnpa(200),vnpw(200),vnpa(200)
!
	dimension id(4,nnp),na(mdof),s1(ls),s2(ls1),ie3d1(npbmax,nbed11),ie3d2(npbmax,nbed21),&
			  ie3d3(npbmax+1),xgihd2(n12gh,n12gh,nbed21),&
			  xgihd3(n12gh,n12gh),ide(n12gh),jdf(4),r(mdof),&
			  xgihd1(n12gh,n12gh,nbed11)
!
!
!	------------------- initializing -------------------
!
	do i=1,ls
		s1(i)=0.d0
	enddo
	if (isolv.ne.3) goto 300 ! symmetric case
	do i=1,ls1
		s2(i)=0.d0
	enddo
300	continue
!
!
!	------------------- bem D1 (finite) area -------------------
!
	if (nbed1.eq.0) goto 500
	do 450 k1=1,nbed1
		do k2=1,nnpbed1(k1)
			ii=ie3d1(k2,k1)
			do k=1,kbem
				jdf(k)=kbem*k2-kbem+k
				ide(jdf(k))=id(k,ii)
			enddo
		enddo
		ndfe=kbem*nnpbed1(k1)
		do 430 j=1,ndfe
			k=ide(j)
			if (k.gt.mdof) goto 430 ! node with boundary conditions
			l=na(k)-k
			do 420 i=1,ndfe
				m1=ide(i)
				if (m1.gt.k.OR.m1.gt.mdof) goto 420
				m1=m1+l
				s1(m1)=s1(m1)+xgihd1(i,j,k1)
				s2(m1)=s2(m1)+xgihd1(j,i,k1)
420			continue
430		continue
450	continue
500	continue
!
!
!	------------------- bem D2 (infinite) area -------------------
!
	if (nbed2.eq.0) goto 600
	do 550 k1=1,nbed2
		do k2=1,nnpbed2(k1)
			ii=ie3d2(k2,k1)
			do k=1,kbem
				jdf(k)=kbem*k2-kbem+k
				ide(jdf(k))=id(k,ii)
			enddo
		enddo
		ndfe=kbem*nnpbed2(k1)
		do 530 j=1,ndfe
			k=ide(j)
			if (k.gt.mdof) goto 530
			l=na(k)-k
			do 520 i=1,ndfe
				m1=ide(i)
				if (m1.gt.k.OR.m1.gt.mdof) goto 520
				m1=m1+l
				s1(m1)=s1(m1)+xgihd2(i,j,k1)
				s2(m1)=s2(m1)+xgihd2(j,i,k1)
520			continue
530		continue
550	continue
600	continue
!
!
!	------------------- bem D3 (semi-infinite) area -------------------
!
	if (nbed3.eq.0) goto 700
	do k2=1,nnpbed3
		ii=ie3d3(k2)
		do k=1,kbem
			jdf(k)=kbem*k2-kbem+k
			ide(jdf(k))=id(k,ii)
		enddo
	enddo
	ndfe=kbem*nnpbed3
	do 630 j=1,ndfe
		k=ide(j)
		if (k.gt.mdof) goto 630
		l=na(k)-k
		do 620 i=1,ndfe
			m1=ide(i)
			if (m1.gt.k.OR.m1.gt.mdof) goto 620
			m1=m1+l
			s1(m1)=s1(m1)+xgihd3(i,j)
			s2(m1)=s2(m1)+xgihd3(j,i)
620		continue
630	continue
700	continue
!
	if (ifem.eq.1) then
		write (*,*) 'ERROR! THE BEM RIGIDITY MATRIX IS MODifIED FOR INITIAL VALUE JUST FOR &
					 BEM ANALYSIS'
		goto 1
	endif
!
!	modify stiffness matrix for imposed pw-value
!
	if (npw.eq.0) goto 6000
	do 6100 l=1,npw
		i=nnpw(l)
		j=id(3,i)
		if (j.gt.mdof) goto 6100
		ndbl=0
		if (j.gt.1) ndbl=na(j-1)
		nrow=j+1-(na(j)-ndbl)
		do 6200 k=nrow,mdof
			if (k.gt.j) then
				lpl=j-k+na(k)
				if (lpl.le.na(k-1)) goto 6200
				r(j)=r(j)-s1(lpl)*(vnpw(l)) !-dt(j))
!
				s1(lpl)=0.D0
				if (isolv.eq.3) then
					r(j)=r(j)-s2(lpl)*(vnpw(l)) !-dt(j))
					s2(lpl)=0.D0
				endif
			else
				lpl=k-j+na(j)
				r(j)=r(j)-s1(lpl)*(vnpw(l)) !-dt(j))
!
				s1(lpl)=0.D0
				if (isolv.eq.3) then
					r(j)=r(j)-s2(lpl)*(vnpw(l)) !-dt(j))
					s2(lpl)=0.D0
				endif
			endif
6200	continue
!
		s1(na(j))=1.D0
		r(j)=(vnpw(l)) !-dt(j))
		if (isolv.eq.3) then
			s2(na(j))=1.D0
			r(j)=(vnpw(l)) !-dt(j))
		endif
6100 continue
6000 continue
!
!	modify stiffness matrix for imposed pa-value
!
	if (npa.eq.0) goto 6500
	do 6600 l=1,npa
		i=nnpa(l)
		j=id(4,i)
		if (j.gt.mdof) goto 6600
		ndbl=0
		if (j.gt.1) ndbl=na(j-1)
		nrow=j+1-(na(j)-ndbl)
		do 6700 k=nrow,mdof
			if (k.gt.j) then
				lpl=j-k+na(k)
				if (lpl.le.na(k-1)) goto 6700
				r(j)=r(j)-s1(lpl)*(vnpa(l)) !-dt(j))
!
				s1(lpl)=0.D0
				if (isolv.eq.3) then
					r(j)=r(j)-s2(lpl)*(vnpa(l)) !-dt(j))
					s2(lpl)=0.D0
				endif
			else
				lpl=k-j+na(j)
				r(j)=r(j)-s1(lpl)*(vnpa(l)) !-dt(j))
				s1(lpl)=0.D0
				if (isolv.eq.3) then
					r(j)=r(j)-s2(lpl)*(vnpa(l)) !-dt(j))
					s2(lpl)=0.D0
				endif
			endif
6700	continue
		s1(na(j))=1.D0
		r(j)=(vnpa(l)) !-dt(j))
!
		if (isolv.eq.3) then
			s2(na(j))=1.D0
			r(j)=(vnpa(l)) !-dt(j))
		endif
6600 continue
6500 continue
!
1	continue
!
	return
	end
!****************************************************************************************
!
!							           BEMT1 subroutine
!
!	This subroutine controlls BEM at first time increment.
!	This subroutine is called in: BEM
!	This subroutine calls:		  GHMATD1,GHMATD2,GHMATD3
!
!	variables used are:
!
!		itm		:
!		x		: coordinates of nodes (nodal data) in FE & BE areas;
!							x(1,1:nnp) : X-coord
!							x(2,1:nnp) : Y-coord
!		kbem	: number of degree of freedom per node in BE zone [2,3,4];
!							if IBEM=1 -> KBEM=2 [X,Y]
!							if IBEM=2 -> KBEM=3 [X,Y,Pw]
!							if IBEM=4 -> KBEM=4 [X,Y,Pw,Pa]
!		npbmax	: maximum nodes in all of the BE zone + 1; for example if we consider D3 BE
!				  area NPBMAX = nnpbed3 + 1
!		n12gh	: in geomechanical problems there is more than one unknown per node, so the
!				  matrices in BEM are assambled according to the doF, rather than node
!				  number. for dry soil each node has [2] doF (ux,uy), for for saturated soil
!				  each node has [3] (ux,uy,pw) doF & for unsaturated soil each node has
!				  [4] doF(ux,uy,pw,pa). then total number of doF in BE=(npbmax-1)*kbem
!		n3gh	: in semi-infinite problem in BEM we have to consider a Fictitious enclosing
!				  surface to perform point collocation method. N3GH shows the total number
!				  of doF in (D3+enclosing surface) in BE technique=(npbmax-1+nnpen)*kbem
!		gd1		: fundamental solution matrix of displacement/pressure
!		gd2		: fundamental solution matrix of displacement/pressure
!		gd3		: fundamental solution matrix of displacement/pressure; total G matrix in D3
!				  (semi-infinite) area concerning all of the g matrices of each element (egd3)
!				  global coefficient g matrices assembled by gathering element contributions
!		hd1		: fundamental solution matrix of stress/flow
!		hd2		: fundamental solution matrix of stress/flow
!		hd3		: fundamental solution matrix of stress/flow; total H matrix in D3
!				  (semi-infinite) area concerning all of the h matrices of each element (ehd3)
!				  global coefficient h matrices assembled by gathering element contributions
!		ed3		: ??????????????????????????????????????????
!		gid1	: inverse matrice of gd1, G1^(-1)
!		gid2	: inverse matrice of gd2, G2^(-1)
!		gid3	: inverse matrice of gd3, G3^(-1)
!		xgihd1	: M*gid1*hd1, M*G1^(-1)*H1
!		xgihd2	: M*gid1*hd1, M*G1^(-1)*H2
!		xgihd3	: M*gid1*hd1, M*G1^(-1)*H3
!		gihd1	: gid1*hd1, G1^(-1)*H1
!		gihd2	: gid2*hd2, G2^(-1)*H2
!		gihd3	: gid3*hd3, G3^(-1)*H3
!		xgid1	: M*gid1, M*G1^(-1)
!		xgid2	: M*gid2, M*G2^(-1)
!		xgid3	: M*gid3, M*G3^(-1)
!
!
!	INPUT	:
!
!	OUTPUT	:
!
!
!****************************************************************************************
!
	subroutine bemt1(itm,x,sm8d,sm8c,sm8u,ie3d1,ie3d2,ie3d3,xgihd1,xgihd2,xgihd3,gihd1,&
					 gihd2,gihd3,gid1,gid2,gid3,xgid1,xgid2,xgid3,disp,id,sig3u)
!
	implicit double precision (a-h,o-z)
!
	common /a/ nnp,ifem,ibem,idyn,ieaq,nr
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /a3/ ne8d1,ne8c1,ne8u1,ne3u1,nload1,nbcx1,nbcy1,nbcw1,nbca1,m8d1,m8c1,m8u1,&
				nbed11,nbed21,nbed31,nbeint1,ne8dd,ne8cc,ne8uu
	common /a4/ xencl(500),yencl(500)
	common /b/ m8d,m8c,m8u,m8dep(10),m8cep(10),m8uep(10)
	common /h/ mdof,mdofn,ls,ls1,isolv
!
	dimension x(2,nnp),sm8d(60,m8d1),sm8c(60,m8c1),sm8u(60,m8u1),ie3d1(npbmax,nbed11),&
			  ie3d2(npbmax,nbed21),ie3d3(npbmax),disp(mdofn),&
			  xgihd2(n12gh,n12gh,nbed21),xgihd3(n12gh,n12gh),&
			  gihd2(n12gh,n12gh,nbed21),gihd3(n12gh,n12gh),&
			  gid2(n12gh,n12gh,nbed21),gid3(n12gh,n12gh),&
			  xgid2(n12gh,n12gh,nbed21),xgid3(n12gh,n12gh),hd1(n12gh,n12gh),&
			  gd1(n12gh,n12gh),hd2(n12gh,n12gh),gd2(n12gh,n12gh),hd3(n12gh,n3gh),&
			  gd3(n12gh,n3gh),id(4,nnp),sig3u(5,500,nbed11),&
			  xgihd1(n12gh,n12gh,nbed11),&
			  gihd1(n12gh,n12gh,nbed11),&
			  gid1(n12gh,n12gh,nbed11),&
			  xgid1(n12gh,n12gh,nbed11)

	dimension ed3(4,n3gh)  ! dimension ed3????????????????
!
!
!	D1 bem area (finite zone)
	if (nbed1.eq.0) goto 40
	do k1=1,nbed1
		call ghmatd1(itm,k1,x,sm8d,sm8c,sm8u,ie3d1,hd1,gd1,gid1,gihd1,xgid1,xgihd1,&
					 sig3u)
	enddo
	write (*,*)'ghmatd1 concluded'
40	continue
!
!	D2 bem area (infinite zone)
	if (nbed2.eq.0) goto 60
	do k1=1,nbed2
!		call ghmatd2(itm,k1,x,sm8d,sm8c,sm8u,ie3d2,hd2,gd2,gid2,gihd2,xgid2,xgihd2)
	enddo
	write (*,*)'ghmatd2 concluded'
60	continue
!
!	D3 bem area (semi-infinite zone)
	if (nbed3.eq.0) goto 80
	call ghmatd3(itm,x,sm8d,sm8c,sm8u,ie3d3,hd3,gd3,gid3,gihd3,xgid3,xgihd3,ed3)
	write (*,*)'ghmatd3 concluded'
80	continue
!
	return
	end
!****************************************************************************************
!
!							           BEMTD1 subroutine
!
!	This subroutine calculates the tractions of area D1 at time increment itime-1.
!	this subroutine is called in this section: BEMTND1
!
!	variables used are:
!
!		tbed1	: traction/flow [tx,ty,qw,qa] vector T(N) at time: "itime-1"
!		ubed1	: displacement/pressure [ux,uy,Pw,Pa] vector U(N) at time: "itime-1" which
!				  is calculated in the previous time step in SOLVE subroutine &
!				  replaced in "DISP" vector in UPDATE3 subroutine
!		rbed1	: stress tensor in D3 BEM area which includes JUST effect of the past
!				  dynamic history; Z(N)=SUM(n=1:N-1) [G(N+1-n)*T(n)-H(N+1-n)*U(n)];
!				  in the first time step Z(1)=rbed1(1,ii)=0
!		gihd1	: gid1*hd1, G1^(-1)*H1
!		gid1	: inverse matrice of gd1, G1^(-1)
!		disp	: displacement/pressures vector which consists of initial displacement
!				  value of corresponding nodes in which their identities code of Ux & Uy
!				  are defined in ID and initial water and air values of corresponding
!				  nodes in which their identities code of Pw & Pa are defined in ID.
!				  these displacement initial values can be defined for a free node or
!				  for a node with boundary conditions and the water and air initial values
!				  must be defined for a free node; in each "itime" the displacements/
!				  pressures vector (R->DT) in PREVIOUS TIME STEP (N-1) is transfered
!				  in DISP vector
!		ieaq	: earthquake code [0,1,2,3];
!							[0]: static case
!							[1]: accelaration
!							[2]: imposed deplacement
!							[3]: incident wave, in most of the time IEAQ=[3]
!		ueaq	: if TETA-method is used, displacement field due to incident waves transforms
!				  to "Ueq(N)-[(teta-1)/teta]*Uex(N-1)"
!				  in this subroutine, vector of U_inc at time: "itime-1"
!		ubed1	: displacement/pressure [ux,uy,Pw,Pa] vector U(N) at time: "itime-1" which
!				  is calculated in the previous time step in SOLVE subroutine &
!				  replaced in "DISP" vector in UPDATE3 subroutine
!		disp	: displacement/pressures vector which consists of initial displacement
!				  value of corresponding nodes in which their identities code of Ux & Uy
!				  are defined in ID and initial water and air values of corresponding
!				  nodes in which their identities code of Pw & Pa are defined in ID.
!				  these displacement initial values can be defined for a free node or
!				  for a node with boundary conditions and the water and air initial values
!				  must be defined for a free node; in each "itime" the displacements/
!				  pressures vector (R->DT) in PREVIOUS TIME STEP (N-1) is transfered
!				  in DISP vector
!		uxeq	: (incident wave) displacement in x dir
!							uxeq (3,1:nnpbed3): in current time -> itime
!							uxeq (2,1:nnpbed3): in previous time -> itime-1
!							uxeq (1,1:nnpbed3): in two previous time -> itime-2
!		uyeq	: (incident wave) displacement in y dir
!							uyeq (3,1:nnpbed3): in current time -> itime
!							uyeq (2,1:nnpbed3): in previous time -> itime-1
!							uyeq (1,1:nnpbed3): in two previous time -> itime-2
!
!	INPUT	:
!
!	OUTPUT	: ubed1,tbed1
!
!
!****************************************************************************************
!
	subroutine bemtd1(k1,id,ie3d1,tbed1,ubed1,rbed1,gihd1,gid1,disp)
!
	implicit double precision (a-h,o-z)
!
	common /a/ nnp,ifem,ibem,idyn,ieaq,nr
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /a3/ ne8d1,ne8c1,ne8u1,ne3u1,nload1,nbcx1,nbcy1,nbcw1,nbca1,m8d1,m8c1,m8u1,&
				nbed11,nbed21,nbed31,nbeint1,ne8dd,ne8cc,ne8uu
	common /c/ nload,ntime,itmax,itime,iter,iload,time,dtime,dtimei,dtimem,xinc
	common /h/ mdof,mdofn,ls,ls1,isolv
!
	dimension ie3d1(npbmax,nbed11),disp(mdofn),id(4,nnp),tbed1(ntime,n12gh,nbed11),&
			  ubed1(ntime,n12gh,nbed11),rbed1(ntime,n12gh,nbed11),&
			  gihd1(n12gh,n12gh,nbed11),gid1(n12gh,n12gh,nbed11)
!
!
!	displacement/pressure vector ubed1 at time: itime-1
!
	it=itime-1
	do k=1,nnpbed1(k1)
		j=ie3d1(k,k1)
		do i=1,kbem
			ubed1(it,kbem*k-kbem+i,k1)=disp(id(i,j))
		enddo
	enddo
!
!	traction/flow vector tbed1
!
	n=nnpbed1(k1)*kbem
	do i=1,n
		tbed1(it,i,k1)=0.d0
		do j=1,n
			tbed1(it,i,k1) = tbed1(it,i,k1)+gihd1(i,j,k1)*ubed1(it,j,k1)-&
							 gid1(i,j,k1)*rbed1(it,j,k1)
		enddo
	enddo
!
	return
	end
!****************************************************************************************
!
!							           BEMTN subroutine
!
!	This subroutine controlls BEM at higher time increments.
!	This subroutine is called in: BEM
!	This subroutine calls:		  BEMTND1,BEMTND2,BEMTND3
!
!	variables used are:
!
!		id		: degrees of freedom of each node from "1" to "mdof"; for the nodes with
!				  the boundary conditions, each node takes from "mdof+1" to "mdof+nbc"
!				  [X,Y]; which are free for the dry elm &
!				  [X,Y,Pw]; which are free for the saturated elm
!				  [X,Y,Pw,Pa]; which are free for the unsaturated elm
!		gid1	: inverse matrice of gd1, G1^(-1)
!		gid2	: inverse matrice of gd2, G2^(-1)
!		gid3	: inverse matrice of gd3, G3^(-1)
!		xgihd1	: M*gid1*hd1, M*G1^(-1)*H1
!		xgihd2	: M*gid1*hd1, M*G1^(-1)*H2
!		xgihd3	: M*gid1*hd1, M*G1^(-1)*H3
!		gihd1	: gid1*hd1, G1^(-1)*H1
!		gihd2	: gid2*hd2, G2^(-1)*H2
!		gihd3	: gid3*hd3, G3^(-1)*H3
!		xgid1	: M*gid1, M*G1^(-1)
!		xgid2	: M*gid2, M*G2^(-1)
!		xgid3	: M*gid3, M*G3^(-1)
!		tbed1	:
!		tbed2	:
!		tbed3	: traction/flow [tx,ty,qw,qa] vector T(N) at time: "itime-1"
!		ubed1	:
!		ubed2	:
!		ubed3	: displacement/pressure [ux,uy,Pw,Pa] vector U(N) at time: "itime-1" which
!				  is calculated in the previous time step in SOLVE subroutine &
!				  replaced in "DISP" vector in UPDATE3 subroutine
!		rbed1	:
!		rbed2	:
!		rbed3	: stress tensor in D3 BEM area which includes JUST effect of the past
!				  dynamic history; Z(N)=SUM(n=1:N-1) [G(N+1-n)*T(n)-H(N+1-n)*U(n)];
!				  in the first time step Z(1)=rbed3(1,ii)=0
!		frbed1	:
!		frbed2	:
!		frbed3	:
!		disp	: displacement/pressures vector which consists of initial displacement value
!				  of corresponding nodes in which their identities code of Ux & Uy are defined
!				  in ID and initial water and air values of corresponding nodes in which their
!				  identities code of Pw & Pa are defined in ID. these displacement initial
!				  values can be defined for a free node or for a node with boundary conditions
!				  and the water and air initial values must be defined for a free node;
!				  in each "itime" the displacements/pressures vector (R->DT) in PREVIOUS
!				  TIME STEP (N-1) is transfered in DISP vector
!
!	INPUT	:
!
!	OUTPUT	:
!
!
!****************************************************************************************
!
	subroutine bemtn(id,x,sm8d,sm8c,sm8u,ie3d1,ie3d2,ie3d3,gihd1,gihd2,gihd3,gid1,gid2,gid3,&
					 xgihd1,xgihd2,xgihd3,xgid1,xgid2,xgid3,tbed1,tbed2,tbed3,ubed1,ubed2,&
					 ubed3,rbed1,rbed2,rbed3,frbed1,frbed2,frbed3,disp,sig3u)
!
	implicit double precision (a-h,o-z)
!
	common /a/ nnp,ifem,ibem,idyn,ieaq,nr
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /a3/ ne8d1,ne8c1,ne8u1,ne3u1,nload1,nbcx1,nbcy1,nbcw1,nbca1,m8d1,m8c1,m8u1,&
				nbed11,nbed21,nbed31,nbeint1,ne8dd,ne8cc,ne8uu
	common /a4/ xencl(500),yencl(500)
	common /c/ nload,ntime,itmax,itime,iter,iload,time,dtime,dtimei,dtimem,xinc
	common /h/ mdof,mdofn,ls,ls1,isolv
!
	dimension id(4,nnp),x(2,nnp),sm8d(60,m8d1),sm8c(60,m8c1),sm8u(60,m8u1),&
			  ie3d1(npbmax,nbed11),ie3d2(npbmax,nbed21),ie3d3(npbmax),&
			  gihd2(n12gh,n12gh,nbed21),gihd3(n12gh,n12gh),&
			  gid2(n12gh,n12gh,nbed21),gid3(n12gh,n12gh),&
			  xgihd2(n12gh,n12gh,nbed21),xgihd3(n12gh,n12gh),&
			  xgid2(n12gh,n12gh,nbed21),xgid3(n12gh,n12gh),&
			  tbed1(ntime,n12gh,nbed11),tbed2(ntime,n12gh,nbed21),tbed3(ntime,n12gh),&
			  ubed1(ntime,n12gh,nbed11),ubed2(ntime,n12gh,nbed21),ubed3(ntime,n12gh),&
			  rbed1(ntime,n12gh,nbed11),rbed2(ntime,n12gh,nbed21),rbed3(ntime,n12gh),&
			  frbed1(n12gh,nbed11),frbed2(n12gh,nbed21),frbed3(n12gh),disp(mdofn),&
			  sig3u(5,500,nbed11),&
			  xgihd1(n12gh,n12gh,nbed11),&
			  gihd1(n12gh,n12gh,nbed11),&
			  gid1(n12gh,n12gh,nbed11),&
			  xgid1(n12gh,n12gh,nbed11)

!
!
!	---------------------- BEM area D1 ----------------------
!
	if (nbed1.eq.0) goto 40
	call bemtnd1(id,x,sm8d,sm8c,sm8u,ie3d1,tbed1,ubed1,rbed1,frbed1,gihd1,gid1,xgihd1,&
				 xgid1,disp,sig3u)
	write (*,*)'          bemtnd1 concluded'
40	continue
!
!
!	---------------------- BEM area D2 ----------------------
!
	if (nbed2.eq.0) goto 60
!	call bemtnd2(id,x,sm8d,sm8c,sm8u,ie3d2,tbed2,ubed2,rbed2,frbed2,gihd2,gid2,xgihd2,&
!				 xgid2,disp,sig3u)
	write (*,*)'          bemtnd2 concluded'
60	continue
!
!
!	---------------------- BEM area D3 ----------------------
!
	if (nbed3.eq.0) goto 80
	call bemtnd3(id,x,sm8d,sm8c,sm8u,ie3d3,tbed3,ubed3,rbed3,frbed3,gihd3,gid3,xgihd3,&
				 xgid3,disp,sig3u)
	write (*,*)'          bemtnd3 concluded'
80	continue
!
	return
	end
!****************************************************************************************
!
!							           BEMTND1 subroutine
!
!	This subroutine controlls BEM of area D1 at higher time increments.
!	these subroutines are called in this section:
!		bemtd1	:
!		ghmatd1	:
!		bemrd1	:
!
!	variables used are:
!
!		ubed1		: displacement/pressure [ux,uy,Pw,Pa] vector U(N) at time: "itime-1" which
!					  is calculated in the previous time step in SOLVE subroutine &
!					  replaced in "DISP" vector in UPDATE3 subroutine
!		tbed1		: traction/flow [tx,ty,qw,qa] vector T(N) at time: "itime-1"
!		rbed1		: stress tensor in D3 BEM area which includes JUST effect of the past
!					  dynamic history; Z(N)=SUM(n=1:N-1) [G(N+1-n)*T(n)-H(N+1-n)*U(n)];
!					  in the first time step Z(1)=rbed3(1,ii)=0
!		red1		:
!		kcorf		:
!
!	INPUT	:
!
!	OUTPUT	:
!
!
!****************************************************************************************
!
	subroutine bemtnd1(id,x,sm8d,sm8c,sm8u,ie3d1,tbed1,ubed1,rbed1,frbed1,gihd1,gid1,&
					   xgihd1,xgid1,disp,sig3u)
!
	implicit double precision (a-h,o-z)
!
	common /a/ nnp,ifem,ibem,idyn,ieaq,nr
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /a3/ ne8d1,ne8c1,ne8u1,ne3u1,nload1,nbcx1,nbcy1,nbcw1,nbca1,m8d1,m8c1,m8u1,&
				nbed11,nbed21,nbed31,nbeint1,ne8dd,ne8cc,ne8uu
	common /b/ m8d,m8c,m8u,m8dep(10),m8cep(10),m8uep(10)
	common /c/ nload,ntime,itmax,itime,iter,iload,time,dtime,dtimei,dtimem,xinc
	common /h/ mdof,mdofn,ls,ls1,isolv
	common /stability/anew1,anew2,wil,intb
!
	dimension id(4,nnp),x(2,nnp),sm8d(60,m8d1),sm8c(60,m8c1),sm8u(60,m8u1),&
			  ie3d1(npbmax,nbed11),tbed1(ntime,n12gh,nbed11),ubed1(ntime,n12gh,nbed11),&
			  rbed1(ntime,n12gh,nbed11),frbed1(n12gh,nbed11),&
			  disp(mdofn),gd1(n12gh,n12gh),hd1(n12gh,n12gh),sig3u(5,500,nbed11),&
			  xgihd1(n12gh,n12gh,nbed11),&
			  gihd1(n12gh,n12gh,nbed11),&
			  gid1(n12gh,n12gh,nbed11),&
			  xgid1(n12gh,n12gh,nbed11)
!
!
	do 1000 k1=1,nbed1
!
!	------------- calculate traction vector of time step itime-1 (i.e. n-1) -------------
!
		call bemtd1(k1,id,ie3d1,tbed1,ubed1,rbed1,gihd1,gid1,disp)
!
!
!	--------------- evaluating effects of previous time increments ---------------
!	------------------------------	on the current one -----------------------------
!
!	initializing vector rbed1
!
		n=nnpbed1(k1)*kbem
		do i=1,n
			rbed1(itime,i,k1)=0.d0
		enddo
!
!	this loop is dedicated to calculate Z(N)=SUM(n=1:N-1) [G(N+1-n)*T(n)-H(N+1-n)*U(n)]
		do 500 itm=1,itime-1
!
!	calculating G & H matices at itime > 1
			call ghmatd1(itm,k1,x,sm8d,sm8c,sm8u,ie3d1,hd1,gd1,gid1,gihd1,xgid1,xgihd1,&
						 sig3u)
!
!	increasing rbed1 vector due to itm time increment
			call bemrd1(itm,k1,tbed1,ubed1,rbed1,hd1,gd1)
!
500		continue
!
!
!	----- last moment itime-1 -----
!
		itm=itime
		call ghmatd1(itm,k1,x,sm8d,sm8c,sm8u,ie3d1,hd1,gd1,gid1,gihd1,xgid1,xgihd1,sig3u)
!
		itm=itime-1
		do i=1,n
			do j=1,n
				rbed1(itime,i,k1)= rbed1(itime,i,k1)+( gd1(i,j)*tbed1(itm,j,k1)-&
								   hd1(i,j)*ubed1(itm,j,k1) )*(1.d0-wil)/wil
			enddo
		enddo
!
!
!	----- convert vector rbed1 to its corresponding force vector form -----
!
		call bemfd1(k1,frbed1,rbed1,xgid1)
!
1000 continue
!
	return
	end
!****************************************************************************************
!
!							           BMAT3N subroutine
!
!	This subroutine calculates the strain displacement matrix for a 3 noded element.
!	{dN/dx;dN/dy}=(1/detJ)[dy/dt,-dy/ds;-dx/dt,dx/ds]*{dN/ds;dN/dt}
!	eps=B*U
!	B=[dN1/dx,0,dN2/dx,0,dN3/dx,0;
!	   0,dN1/dy,0,dN2/dy,0,dN3/dy;
!	   dN1/dy,dN1/dx,dN2/dy,dN2/dx,dN3/dy,dN3/dx]
!	this subroutine is called in: STR3ND
!
!	variables used are:
!
!	INPUT	:
!
!	OUTPUT	:
!
!
!****************************************************************************************
!
	subroutine bmat3n(xe,ye,pfs,b,detj)
!
 	implicit double precision (a-h,o-z)
!
	dimension xe(3),ye(3),pfs(3),b(3,6)
!
	pxs=0.0d0
	pys=0.0d0
!
	do 100 i=1,3
		pxs=pxs+pfs(i)*xe(i)
        pys=pys+pfs(i)*ye(i)
100	continue
!
	detj=dsqrt(pxs**2.+pys**2.)
!
!	ii=1
!	jj=2
!	do 200 i=1,8
!		b(1,ii)=(pfs(i))/detj
!       b(1,jj)=0.0
!       b(2,jj)=(pfs(i))/detj
!       b(3,ii)=b(2,jj)
!       b(3,jj)=b(1,ii)
!       ii=ii+2
!       jj=jj+2
!200	continue
!
	return
	end
!****************************************************************************************
!
!							           BMAT8N subroutine
!
!	This subroutine calculates the strain displacement matrix for a 8 noded element.
!	{dN/dx;dN/dy}=(1/detJ)[dy/dt,-dy/ds;-dx/dt,dx/ds]*{dN/ds;dN/dt}
!	eps=B*U
!	B=[dN1/dx,0,dN2/dx,0,dN3/dx,0,dN4/dx,0,dN5/dx,0,dN6/dx,0,dN7/dx,0,dN8/dx,0;
!	   0,dN1/dy,0,dN2/dy,0,dN3/dy,0,dN4/dy,0,dN5/dy,0,dN6/dy,0,dN7/dx,0,dN8/dx;
!	   dN1/dy,dN1/dx,dN2/dy,dN2/dx,dN3/dy,dN3/dx,dN4/dy,dN4/dx,dN5/dy,dN5/dx,&
!					dN6/dy,dN6/dx,dN7/dy,dN7/dx,dN8/dy,dN8/dx				   ]
!	this subroutine is called in: STR8ND
!
!	variables used are:
!
!	INPUT	:
!
!	OUTPUT	:
!
!
!****************************************************************************************
!
	subroutine bmat8n(xe,ye,pfs,pft,b,detj)
!
 	implicit double precision (a-h,o-z)
!
	dimension xe(8),ye(8),pfs(8),pft(8),b(3,16)
!
	pxs=0.0d0
	pxt=0.0d0
	pys=0.0d0
	pyt=0.0d0
!
	do 100 i=1,8
		pxs=pxs+pfs(i)*xe(i)
        pxt=pxt+pft(i)*xe(i)
        pys=pys+pfs(i)*ye(i)
        pyt=pyt+pft(i)*ye(i)
100	continue
!
	detj=pxs*pyt-pxt*pys
!
	ii=1
	jj=2
	do 200 i=1,8
		b(1,ii)=(pyt*pfs(i)-pys*pft(i))/detj
        b(1,jj)=0.0
        b(2,jj)=(pxs*pft(i)-pxt*pfs(i))/detj
        b(2,ii)=0.0
        b(3,ii)=b(2,jj)
        b(3,jj)=b(1,ii)
        ii=ii+2
        jj=jj+2
200	continue
!
	return
	end
!!c ****************************************************************

!!c FONCTIONS DE BESSEL MODifIEES K0,K1 COMPLEX

!!c ****************************************************************

      subroutine cbesk01(s,cbesk0,cbesk1)
      double precision x,y,re0, im0, re1, im1
      double complex s,cbesk0,cbesk1

      x=dble(s) ; y=dimag(s)
      call kzeone (x, y, re0, im0, re1, im1)
      cbesk0=dcmplx(re0,im0) ; cbesk1=dcmplx(re1,im1)
      return
      end

      subroutine kzeone(x, y, re0, im0, re1, im1)
      double precision x, y, x2, y2, re0, im0, re1, im1,r1, r2, t1, t2, p1, p2, rterm, iterm, exsq(8), tsq(8)

!!c the arrays tsq and exsq contain the square of the
!!c abscissas and the weight factors used in the gauss-
!!c hermite quadrature.

	  data tsq(1) /0.0d0/, tsq(2) /3.19303633920635d-1/,tsq(3) /1.29075862295915d0/, &
		   tsq(4)/2.95837445869665d0/, tsq(5) /5.40903159724444d0/,tsq(6) /8.80407957805676d0/, & 
		   tsq(7)/1.34685357432515d1/, tsq(8) /2.02499163658709d1/,exsq(1) /0.5641003087264d0/, &
		   exsq(2)/0.4120286874989d0/, exsq(3) /0.1584889157959d0/,exsq(4) /0.3078003387255d-1/, &
		   exsq(5)/0.2778068842913d-2/, exsq(6) /0.1000044412325d-3/,exsq(7) /0.1059115547711d-5/, exsq(8)/0.1522475804254d-8/

      r2 = x*x + y*y
      if (x.gt.0.0d0 .or. r2.ne.0.0d0) go to 10
      write (6,100)
      stop
   10 if (r2.ge.1.96d2) go to 50
      if (r2.ge.1.849d1) go to 30

!!c this section calculates the functions using the series expansions
!!c -----------------------------

      x2 = x/2.0d0
      y2 = y/2.0d0
      p1 = x2*x2
      p2 = y2*y2
      t1 = -(dlog(p1+p2)/2.0d0+0.5772156649015329d0)
!!c this constant is euler*s constant
      t2 = -datan2(y,x)
      x2 = p1 - p2
      y2 = x*y2
      rterm = 1.0d0
      iterm = 0.0d0
      re0 = t1
      im0 = t2
      t1 = t1 + 0.5d0
      re1 = t1
      im1 = t2
      p2 = dsqrt(r2)
      l = 2.106d0*p2 + 4.4d0
      if (p2.lt.8.0d-1) l = 2.129d0*p2 + 4.0d0
      do 20 n=1,l
        p1 = n
        p2 = n*n
        r1 = rterm
        rterm = (r1*x2-iterm*y2)/p2
        iterm = (r1*y2+iterm*x2)/p2
        t1 = t1 + 0.5d0/p1
        re0 = re0 + t1*rterm - t2*iterm
        im0 = im0 + t1*iterm + t2*rterm
        p1 = p1 + 1.0d0
        t1 = t1 + 0.5d0/p1
        re1 = re1 + (t1*rterm-t2*iterm)/p1
        im1 = im1 + (t1*iterm+t2*rterm)/p1
   20 continue
      r1 = x/r2 - 0.5d0*(x*re1-y*im1)
      r2 = -y/r2 - 0.5d0*(x*im1+y*re1)
      re1 = r1
      im1 = r2
      return

!!c this section calculates the functions using the integral
!!c representation, eqn 3, evaluated with 15 point gauss-
!!c hermite quadrature
!!c ---------------------------

   30 x2 = 2.0d0*x
      y2 = 2.0d0*y
      r1 = y2*y2
      p1 = dsqrt(x2*x2+r1)
      p2 = dsqrt(p1+x2)
      t1 = exsq(1)/(2.0d0*p1)
      re0 = t1*p2
      im0 = t1/p2
      re1 = 0.0d0
      im1 = 0.0d0
      do 40 n=2,8
        t2 = x2 + tsq(n)
        p1 = dsqrt(t2*t2+r1)
        p2 = dsqrt(p1+t2)
        t1 = exsq(n)/p1
        re0 = re0 + t1*p2
        im0 = im0 + t1/p2
        t1 = exsq(n)*tsq(n)
        re1 = re1 + t1*p2
        im1 = im1 + t1/p2
   40 continue
      t2 = -y2*im0
      re1 = re1/r2
      r2 = y2*im1/r2
      rterm = 1.41421356237309d0*dcos(y)
      iterm = -1.41421356237309d0*dsin(y)
!!c this constant is sqrt(2.0).
      im0 = re0*iterm + t2*rterm
      re0 = re0*rterm - t2*iterm
      t1 = re1*rterm - r2*iterm
      t2 = re1*iterm + r2*rterm
      re1 = (t1*x + t2*y)/dexp(x)
      im1 = (-t1*y + t2*x)/dexp(x)
      re0=re0/dexp(x)
      im0=im0/dexp(x)
      return

!!c this section calculates the functions using the asymptotic expansions
!!c -----------------------------------

   50 rterm = 1.0d0
      iterm = 0.0d0
      re0 = 1.0d0
      im0 = 0.0d0
      re1 = 1.0d0
      im1 = 0.0d0
      p1 = 8.0d0*r2
      p2 = dsqrt(r2)
      l = 3.91d0+8.12d1/p2
      r1 = 1.0d0
      r2 = 1.0d0
      m = -8
      k = 3
      do 60 n=1,l
        m = m + 8
        k = k - m
        r1 = float(k-4)*r1
        r2 = float(k)*r2
        t1 = float(n)*p1
        t2 = rterm
        rterm = (t2*x+iterm*y)/t1
        iterm = (-t2*y+iterm*x)/t1
        re0 = re0 + r1*rterm
        im0 = im0 + r1*iterm
        re1 = re1 + r2*rterm
        im1 = im1 + r2*iterm
   60 continue
      t1 = dsqrt(p2+x)
      t2 = -y/t1
      p1 = 8.86226925452758d-1/p2
!!c this constant is sqrt(pi)/2.0
      rterm = p1*dcos(y)
      iterm = -p1*dsin(y)
      r1 = re0*rterm - im0*iterm
      r2 = re0*iterm + im0*rterm
      re0 = t1*r1 - t2*r2
      im0 = t1*r2 + t2*r1
      r1 = re1*rterm - im1*iterm
      r2 = re1*iterm + im1*rterm
      re1 = (t1*r1 - t2*r2)/dexp(x)
      im1 = (t1*r2 + t2*r1)/dexp(x)
      re0=re0/dexp(x)
      im0=im0/dexp(x)

      return
100   format ('CBESK01: argument is zero, or lies in left halfcomplex plane')
      end
!****************************************************************************************
!
!							           CHARGE1 subroutine
!
!	This subroutine checks for loading step and zeros the rrr vector.
!	This subroutine is called in: CONTROL
!	This subroutine calls: -
!
!	variables used are:
!
!		iconst	: iconst(1,nload1)= Initial load step number (time step)
!				  iconst(2,nload1)= Number of load steps (nstep)
!				  SUM[iconst(2,1:nload)]=ntime
!				  iconst(3,nload1)= Maximum iterations per load step (itmax)
!		nload	: number of load steps; load is applied in "nload" steps, in each
!				  "iload=1:nload" the characteristics of loading such as 'intensity-nstep-...
!				  are differents
!		rrr		: seismic loadings
!		ntime	: number of time steps or total number of "nstep"s in all "nload"s. loading
!				  is divided into SUM[nstep]s, when the characteristics of "nstep"s are
!				  similar we put these "nsteps" in one group which is named iload
!		itime	: time increment from 1 to "ntime"
!							[1] : first time increment,
!							[>1]: higher time increment
!		k		: if in each time step of calculation (each 1 sec) load is applied -> k=1
!		mdof	: maximum degree of freedom by considering the boundary conditions
!
!	INPUT	: iconst,nload,itime
!
!
!	OUTPUT	: k,ii
!
!****************************************************************************************
!
	subroutine charge1(iconst,rrr,nload,itime,k,ii)
!
	implicit double precision (a-h,o-z)
!
	common /h/ mdof,mdofn,ls,ls1,isolv
!
	dimension iconst(3,nload),rrr(mdof)
!
!
!	---------------- check for loading step ----------------
!
	do ii=1,nload
		j=iconst(1,ii)
		if (j.eq.itime) k=1
		if (k.eq.1) goto 100
	enddo
!
!
!	------------------- zero rrr vector --------------------
!
100	do i=1,mdof
		rrr(i)=0.d0
	enddo
!
	return
	end

!****************************************************************************************
!
!							           CHARGE2 subroutine
!
!	This subroutine determines the loading condition and calls the load routines.
!	for the current load step:
!			- data for loading (in this step)
!			- b. loading and flow data
!			- nb of nodes with given values
!			- n� nd with given pw b.c
!			- given pw values
!			- n� nd with given pa b.c
!			- given pa values
!
!	This subroutine is called in: CONTROL
!	This subroutine calls:		  LOADNODE, FLOWNODE, MMATRICE
!
!	variables used are:
!
!		id		: degrees of freedom of each node from "1" to "mdof"; for the nodes with
!				  the boundary conditions, each node takes from "mdof+1" to "mdof+nbc"
!				  [X,Y]; which are free for the dry elm &
!				  [X,Y,Pw]; which are free for the saturated elm
!				  [X,Y,Pw,Pa]; which are free for the unsaturated elm
!		iconst	: iconst(1,nload1)= Initial load step number (time step)
!				  iconst(2,nload1)= Number of load steps (nstep);
!				  SUM[iconst(2,1:nload)]=ntime
!				  iconst(3,nload1)= Maximum iterations per load step (itmax)
!		x		: coordinates of nodes (nodal data) in FE & BE areas;
!							x(1,1:nnp) : X-coord
!							x(2,1:nnp) : Y-coord
!		rr		: it is a vector which saves the mechanical loadings
!		ntime	: number of time steps or total number of "nstep" in all "nload"s
!				  loading is divided into SUM[nstep]s, when the characteristics of "nstep"s
!				  are similar we put these "nsteps" in one group which is named iload
!		nload	: load is applied in "nload" step, in each "iload=1:nload" the character-
!				  istics of loading such as 'intensity - nstep - ...' are differents
!		ic		: time step's number in each iload in which the load is applied;
!				  iconst(1,iload)
!		nstep	: number of substeps in each "iload"
!		itmax	: number of iterations which in this program it is considered [1]
!		mdof	: maximum degree of freedom by considering the boundary conditions
!		lpri	: print code
!		iterp	: number of iloads for which the results are printed out in the output file
!		itprt	: iload's numbers which are selected to print out the results;
!							[0] if no printout is expected
!		ieaq	: earthquake code [0,1,2,3];
!								[0]: static case
!								[1]: accelaration
!								[2]: imposed displacement
!								[3]: incident wave, in most of the time IEAQ=[3]
!		itime	: time increment from 1 to "ntime"
!							[1] : first time increment,
!							[>1]: higher time increment
!		nstime	: number of seismic events. when dynamic load is applied there are not the
!				  other types of loading
!		dtime	: time step; in seismic loading for all "1:ntime", it is constant
!		xinc	: time step increment factor (coef for time step) which is always equal 1
!		neaq	: number of nodes on bedrock when ieaq=[2], it means that dynamic analysis
!				  will be done with the imposed displacement in certain nodes (nneaq)
!		nneaq	: nodes' numbers on bedrock when ieaq=[2], it means that dynamic analysis
!				  will be done with the imposed displacement in certain nodes (nneaq)
!		uxeaq	: bedrock motion in x direction when ieaq=[2], it means that dynamic analysis
!				  will be done with the imposed displacement in certain nodes (nneaq)
!		uyeaq	: bedrock motion in y direction when ieaq=[2], it means that dynamic analysis
!				  will be done with the imposed displacement in certain nodes (nneaq)
!		axeaq	: earthquake acceleration in x direction when ieaq=[1], it means that dynamic
!				  analysis will be done with the accelerogram data
!		ayeaq	: earthquake acceleration in y direction when ieaq=[1], it means that dynamic
!				  analysis will be done with the accelerogram data
!		nnpbed3	: number of nodes of one D3 BEM area
!		uxeq	: (incident wave) displacement in x dir
!							uxeq (3,1:nnpbed3): in current time -> itime
!							uxeq (2,1:nnpbed3): in previous time -> itime-1
!							uxeq (1,1:nnpbed3): in two previous time -> itime-2
!		uyeq	: (incident wave) displacement in y dir
!							uyeq (3,1:nnpbed3): in current time -> itime
!							uyeq (2,1:nnpbed3): in previous time -> itime-1
!							uyeq (1,1:nnpbed3): in two previous time -> itime-2
!		nnw		: number of surfaces with imposed qw
!		infw	: (1,n) number of first node with imposed qw
!				  (2,n) number of second node with imposed qw
!				  (3,n) number of third node with imposed qw
!		aqw		: value of qw
!		nna		: number of surfaces with imposed qa
!		infa	: (1,n) number of first node with imposed qa
!				  (2,n) number of second node with imposed qa
!				  (3,n) number of third node with imposed qa
!		aqa		: value of qa
!
!		npw		: numbers of nodes with given (temporary) p-water values
!		nnpw	: nodes' numbers with temporary imposed pw (vector of dim = npw)
!		vnpw	: values of imposed p-water  (vector of dim = npw)

!		npa		: numbers of nodes with given (temporary) p-air values
!		nnpa	: nodes' numbers with temporary imposed pa (vector of dim = npa)
!		vnpa	: values of imposed p-air  (vector of dim = npa)
!
!		nnload	: number of loaded surfaces (load node)
!		ip		: (1,i) number of first node of ith loaded surface
!				  (2,i) number of second node of ith loaded surface
!				  (3,i) number of third node of ith loaded surface
!		px		: (1,i) Lineic load 1 in x/r direction
!				  (2,i)	Lineic load 2 in x/r direction
!				  (3,i)	Lineic load 3 in x/r direction
!				  (4,i)	Lineic load 1 in y/z direction
!				  (5,i)	Lineic load 2 in y/z direction
!				  (6,i)	Lineic load 3 in y/z direction
!				  (7,i) 1 or 0 ???????????????????????
!
!	INPUT	:
!	OUTPUT	:
!
!****************************************************************************************
!
	subroutine charge2(id,iconst,x,rr,disp)
!
	implicit double precision (a-h,o-z)
!
	common /a/ nnp,ifem,ibem,idyn,ieaq,nr
	common /a1/ ne8d,ne8c,ne8u
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /a3/ ne8d1,ne8c1,ne8u1,ne3u1,nload1,nbcx1,nbcy1,nbcw1,nbca1,m8d1,m8c1,m8u1,&
				nbed11,nbed21,nbed31,nbeint1,ne8dd,ne8cc,ne8uu
	common /c/ nload,ntime,itmax,itime,iter,iload,time,dtime,dtimei,dtimem,xinc
	common /i/ iterp,itprt(100),ne8do,ne8co,ne8uo,ie8dout(10),ie8cout(10),ie8uout(10),nnpo,&
			   inpout(80)
	common /h/ mdof,mdofn,ls,ls1,isolv
	common /chrg1/ nstep,istep
	common /chrg2/ npw,npa,nnpw(200),nnpa(200),vnpw(200),vnpa(200)
	common /chrg3/ nstime,neaq,nneaq(500),uxeaq,uyeaq,axeaq,ayeaq,uxeq(3,500),uyeq(3,500)
	common /sfn/ nnload,nnw,nna,nnw1,nna1
	common /snwm/ ip(3,100),px(7,100),infw(3,100),aqw(100),aqw1(100),infw1(3,100),&
				  infa(3,100),aqa(100),aqa1(100),infa1(3,100)
				  ! infw(3,100) or infw(2,100)??????????
!
	dimension id(4,nnp),iconst(3,nload1),x(2,nnp),rr(mdof),disp(mdofn)
!
!
!	---------------- initilize vector rr ----------------
!
	do i=1,mdof
		rr(i)=0.d0
	enddo
!
!
!	---------------- determine loading condition ----------------
!
	ic=iconst(1,iload)
	nstep=iconst(2,iload)
	itmax=iconst(3,iload)
!
	lpri=0
	k=ic+nstep-1
	do i=1,iterp
		kk=itprt(i)
		if (kk.eq.k) lpri=1
	enddo
	if (lpri.eq.1) write (7,1) ic,nstep,itmax
!
!
	if (ieaq.eq.0) goto 20	! static case
!
!	---------------- seismic loading condition ----------------
!
!	basic data of seismic loading
	if (itime.gt.1) goto 10
	if (lpri.eq.1) write (7,*) 'basic data of seismic loading'
	read (4,*) nstime,dtime,xinc
	write (7,5) nstime,dtime,xinc
!
	if (ieaq.eq.2) then
		read (4,*) neaq
		read (4,*) (nneaq(i),i=1,neaq)
		if (lpri.eq.1) then
			write (7,4) neaq
			write (7,2)
			write (7,17) (nneaq(i),i=1,neaq)
		endif
	endif
!
!	reading acceleration or displacement
10	if (itime.le.nstime) goto 15	! if itime>nstime the seismic loading is finished &
!									  we want to know the other loading conditions like
!									  water or air flow,disp,...
	ieaq=0
	goto 20
!
15	if (ieaq.eq.1) read (4,*) axeaq,ayeaq ! in FEM (dynamic loading:acceleration)
	if (ieaq.eq.2) read (4,*) uxeaq,uyeaq
	if (ieaq.eq.3) then
		do i=1,nnpbed3
			if (itime.gt.1) then
				uxeq(1,i)=uxeq(2,i)
				uxeq(2,i)=uxeq(3,i)
				uyeq(1,i)=uyeq(2,i)
				uyeq(2,i)=uyeq(3,i)
			else
				uxeq(1,i)=0.d0
				uxeq(2,i)=0.d0
				uyeq(1,i)=0.d0
				uyeq(2,i)=0.d0
			endif
		enddo
		read (4,*) (uxeq(3,i),i=1,nnpbed3)
		read (4,*) (uyeq(3,i),i=1,nnpbed3)
	endif
!
	if (lpri.eq.1) then
		if (ieaq.eq.1) then
			write (7,*) 'earthquake acceleration in x direction: ',axeaq
			write (7,*) 'earthquake acdeleration in y direction: ',ayeaq
		endif
		if (ieaq.eq.2) then
			write (7,*) 'bedrock motion in x direction: ',uxeaq
			write (7,*) 'bedrock motion in y direction: ',uyeaq
		endif
		if (ieaq.eq.3) then
			write (7,*) '(incident wave) displacement in x dir: '
			write (7,22) (uxeq(3,i),i=1,nnpbed3)
			write (7,*) '(incident wave) displacement in y dir: '
			write (7,22) (uyeq(3,i),i=1,nnpbed3)
		endif
	endif
	goto 1900
20	continue
!
!
!	---------------- other loading conditions (non-seismic) ----------------
!
!	read loading data
!
	if (nnw.eq.0) goto 250	! it gives us the value of qw[t(n)] in [t(n),t(n+1)]
	nnw1=nnw
	do n=1,nnw1
		do i=1,3
			infw1(i,n)=infw(i,n)
			aqw1(n)=aqw(n)
		enddo
	enddo
!
250	if (nna.eq.0) goto 300
	nna1=nna
	do n=1,nna1
		do i=1,3
			infa1(i,n)=infa(i,n)
			aqa1(n)=aqa(n)
		enddo
	enddo
!
300	continue
!
	read (4,*) dtime,xinc
	if (lpri.eq.1) write (7,3) dtime,xinc
!
!	read boundary loading and flow conditions
!
	read (4,*) npw,npa,nnload,nnw,nna
	if (lpri.eq.1) write (7,23) npw,npa,nnload,nnw,nna
!
	if (npw.eq.0) goto 75
	read (4,*) (nnpw(i),i=1,npw)
	read (4,*) (vnpw(i),i=1,npw)
	do n=1,npw
		l=nnpw(n)
		k=id(3,l)
		if (itime.eq.1) then
			disp(k)=vnpw(n)
		endif
	enddo
	if (lpri.eq.1) then
		write (7,18)
		write (7,11) (nnpw(i),i=1,npw)
		write (7,19)
		write (7,22) (vnpw(i),i=1,npw)
	endif
!
75	if (npa.eq.0) goto 80
	read (4,*) (nnpa(i),i=1,npa)
	read (4,*) (vnpa(i),i=1,npa)
	do n=1,npa
		l=nnpa(n)
		k=id(4,l)
		if (itime.eq.1) then
			disp(k)=vnpa(n)
		endif
	enddo
	if (lpri.eq.1) then
		write (7,38)
		write (7,11) (nnpa(i),i=1,npa)
		write (7,39)
		write (7,22) (vnpa(i),i=1,npa)
	endif
!
80	if (nnload.eq.0) goto 100
	if (lpri.eq.1) write (7,24)
	i=0
50	read (4,*)  i1,i2,i3
	if (i3.ne.0) then
		i=i+1
        ip(1,i)=i1
		ip(2,i)=i2
		ip(3,i)=i3
		read (4,*) px(1,i),px(2,i),px(3,i),px(4,i),px(5,i),px(6,i)
        px(7,i)=1.d0
		if (i.ge.nnload) goto 40
		goto 50
	endif
!
!	if i3=0:
	read (4,*) pxx,pyy,kode
	if (kode.eq.1) then
		if ( MOD(i2-i1,2).ne.0) then
			write (*,*) '**Error** CHARGE2: ',i1,i2
			stop
		endif
		j=i+(i2-i1)/2
		if (j.gt.nnload) then
			write (*,*) '**Error** CHARGE2: exceed of number of ','load surface',i1,i2
			stop
		endif
        j1=i1
60      i=i+1
        ip(1,i)=j1
		ip(2,i)=j1+1
		ip(3,i)=j1+2
		px(1,i)=pxx
		px(2,i)=pxx
		px(3,i)=pxx
		px(4,i)=pyy
		px(5,i)=pyy
		px(6,i)=pyy
		px(7,i)=1.d0
        j1=j1+2
		if (i.lt.j) goto 60
		if (j.ge.nnload) goto 40
		goto 50
	else
		j=i+i2-i1
		if (j.gt.nnload) then
			write (*,*) '**Error** CHARGE2: exceed of number of ','load node',i1,i2
			stop
		endif
        j1=i1
70      i=i+1
        ip(1,i)=j1
		px(1,i)=pxx
		px(4,i)=pyy
		px(7,i)=0.d0
        j1=j1+1
		if (i.lt.j) goto 70
		if (j.ge.nnload) goto 40
		goto 50
	endif
!
40	do i=1, nnload
		if (lpri.eq.1 .and. px(7,i).eq.1.d0) write (7,25) ip(1,i),ip(2,i),ip(3,i),px(1,i),&
													px(2,i),px(3,i),px(4,i),px(5,i),px(6,i)
		if (lpri.eq.1 .and. px(7,i).eq.0.d0) write (7,28) ip(1,i),px(1,i),px(4,i)
	enddo
!
100	continue
!
	if (nnw.eq.0) goto 85
	if (lpri.eq.1) write (7,26)
	do n=1,nnw
		read (4,*) infw(1,n),infw(2,n),infw(3,n),aqw(n)
		if (lpri.eq.1) write (7,27) infw(1,n),infw(2,n),infw(3,n),aqw(n)
	enddo
!
85	if (nna.eq.0) goto 90
	if (lpri.eq.1) write (7,36)
	do n=1,nna
		read (4,*) infa(1,n),infa(2,n),infa(3,n),aqa(n)
		if (lpri.eq.1) write (7,27) infa(1,n),infa(2,n),infa(3,n),aqa(n)
	enddo
!
90	continue
!
	if (nnload.ne.0) call loadnode(x,id,rr)
	if (nnw.ne.0) call flownode(x,id,rr)
	if (nna.ne.0) call floanode(x,id,rr)
!
!	--------------------- formats ---------------------
!
1	FORMAT (//10x,'data for loading'//&
			  5x,'time step number                              = ',i5/&
			  5x,'number of substeps                            = ',i5/&
			  5x,'number of iterations                          = ',i5//)
5	FORMAT (//10x,'time step data'//&
			  5x,'number of seismic events                      = ',i5/&
			  5x,'time step                                     = ',e10.3/&
			  5x,'time step increment factor                    = ',e10.3//)
4	FORMAT (/5x,'number of nodes on bedrock                    = ',i5/)
2	FORMAT (//10x,'node numbers of nodes on bedrock')
3	FORMAT (//10x,'time step data'//&
			  5x,'time step                                     = ',e10.3/&
			  5x,'time step increment factor                    = ',e10.3//)
11	FORMAT (8i5)
18	FORMAT (/10x,'node numbers of nodes with given p-water b.c.'//)
38	FORMAT (/10x,'node numbers of nodes with given p-air b.c.'//)
19	FORMAT (/10x,'given p-water values for the above nodes'//)
39	FORMAT (/10x,'given p-air values for the above nodes'//)
21	FORMAT (/10x,'given incident wave diplacements for be area d3'//)
22	FORMAT (6(2x,d10.3))
17	FORMAT(6(2x,i5))
23	FORMAT (/10x,'numbers of nodes with given p-water values	= ',i5/&
			 10x,'numbers of nodes with given p-air values		= ',i5/&
			 10x,'number of loading surfaces (load node)		= ', i5/&
			 10x,'number of water flow surfaces					= ', i5/&
			 10x,'number of air flow surfaces					= ', i5//)
24	FORMAT (5x,'n1',5x,'n2',5x,'n3',9x,'fx1',9x,'fx2',9x,'fx3',9x,'fy1',9x,'fy2',9x,&
			'fy3',//)
25	FORMAT (3(2x,i5),6(2x,e10.3))
26  FORMAT (//10x,'n1',10x,'n2',10x,'water flow',//)
36  FORMAT (//10x,'n1',10x,'n2',10x,'air flow',//)
27  FORMAT (7x,i5,7x,i5,7x,i5,10x,e10.3)
28  FORMAT (2x,i5,16x,e10.3,26x,e10.3)
!
1900	continue
!
	return
	end
!
!
!****************************************************************************************
!							           LOADNODE subroutine
!	This subroutine calculates nodal loads due to surface loads (stress to force).
!	this subroutine is called in this section: CHARGE2
!	variables used are:
!		rr		:
!		nnload	: number of loaded surfaces (load node)
!		ip		: (1,i) number of first node of ith loaded surface
!				  (2,i) number of second node of ith loaded surface
!				  (3,i) number of third node of ith loaded surface
!		px		: (1,i) Lineic load 1 in x/r direction
!				  (2,i)	Lineic load 2 in x/r direction
!				  (3,i)	Lineic load 3 in x/r direction
!				  (4,i)	Lineic load 1 in y/z direction
!				  (5,i)	Lineic load 2 in y/z direction
!				  (6,i)	Lineic load 3 in y/z direction
!				  (7,i) 1 or 0 ???????????????????????
!****************************************************************************************
!
	subroutine loadnode(x,id,rr)
!
	implicit double precision (a-h,o-z)
!
	common /a/ nnp,ifem,ibem,idyn,ieaq,nr
	common /h/ mdof,mdofn,ls,ls1,isolv
	common /chrg1/ nstep,istep
	common /sfn/ nnload,nnw,nna,nnw1,nna1
	common /snwm/ ip(3,100),px(7,100),infw(3,100),aqw(100),aqw1(100),infw1(3,100),&
				  infa(3,100),aqa(100),aqa1(100),infa1(3,100)
				  ! infw(3,100) or infw(2,100)??????????
!
	dimension x(2,nnp),id(4,nnp),rr(mdof),q(6),qf(6),xelm(3),yelm(3),xmm(12,12)
!
!
	do i=1,nnload
!
!	loading surface
		if (px(7,i).eq.1.d0) then
			i1=ip(1,i)
			i2=ip(2,i)
			i3=ip(3,i)
			xelm(1)=x(1,i1)
			xelm(2)=x(1,i2)
			xelm(3)=x(1,i3)
			yelm(1)=x(2,i1)
			yelm(2)=x(2,i2)
			yelm(3)=x(2,i3)
			a=xelm(3)-2*xelm(2)+xelm(1)
			b=(xelm(3)-xelm(1))/2
			c=yelm(3)-2*yelm(2)+yelm(1)
			d=(yelm(3)-yelm(1))/2
			do k=1,6
				q(k)=0.d0
				qf(k)=0.d0 !qforce?????
				do j=1,6	!3node*2[Ux,Uy]=6
					xmm(k,j)=0.d0
				enddo
			enddo
			qf(1)=px(1,i)
			qf(2)=px(4,i)
			qf(3)=px(2,i)
			qf(4)=px(5,i)
			qf(5)=px(3,i)
			qf(6)=px(6,i)
!
!	  call subroutine in order to evaluate elements participation in xmtf
!	  in order to transform (line) element nodal tractions to nodal forces.
			call mmatrice(xmm,xelm,yelm,2) ! kode=2; [Fx,Fy]
			do k=1,6
				do j=1,6
					q(k)=q(k)+xmm(k,j)*qf(j)
				enddo
			enddo
			j1=id(1,i1)
			j2=id(2,i1)
			j3=id(1,i2)
			j4=id(2,i2)
			j5=id(1,i3)
			j6=id(2,i3)
			if (j1.le.mdof) rr(j1)=rr(j1)+q(1)/nstep
			if (j2.le.mdof) rr(j2)=rr(j2)+q(2)/nstep
			if (j3.le.mdof) rr(j3)=rr(j3)+q(3)/nstep
			if (j4.le.mdof) rr(j4)=rr(j4)+q(4)/nstep
			if (j5.le.mdof) rr(j5)=rr(j5)+q(5)/nstep
			if (j6.le.mdof) rr(j6)=rr(j6)+q(6)/nstep
!
!     load node
		else
			j1=id(1,i1)
			j2=id(2,i1)
			if (j1.le.mdof) rr(j1)=rr(j1)+px(1,i)/nstep
			if (j2.le.mdof) rr(j2)=rr(j2)+px(4,i)/nstep
		endif
	enddo
!
	return
	end
!
!
!****************************************************************************************
!									FLOWNODE subroutine
!	This subroutine calculates nodal water loads due to surface water flow. as qw(m/s)
!	then we have to integrate it in time to abtain in m.
!	qw(t)=qw(t_n)+( (t-t_n)/(t_n+1-t_n) )*( qw(t_n+1)-qw(t_n) ) also we have:
!	INTG[qw(t)dt;tn,t(n+1)]=[(1-beta)*qw(tn)+beta*qw(t_n+1)]*dtime
!	Dans le syst�me lin�aire ainsi obtenu, le param�tre beta peut �tre consid�r� comme
!	degr� d�implicit� qui correspond au type d�approximation choisie pour l�int�gration
!	dans le temps. beta= 0 correspond au sch�ma explicite (m�thode d�Euler progressive).
!	Cette m�thode ne donne pas de r�sultats satisfaisants pour les �quations diff�rentielles
!	mal conditionn�es. En prenant beta = 1, la r�solution sera implicite (m�thode d�Euler
!	r�trograde) et avec beta= 1/2, on obtient l�algorithme de Crank - Nicholson. Ainsi,
!	le choix de beta influence la stabilit� et la pr�cision de l�algorithme num�rique.
!	Pour beta>=1/2 le sch�ma est inconditionnellement stable quelle que soit la condition
!	initiale et la seule contrainte est celle de la pr�cision (Zienkiewicz, 1977).
!
!	this subroutine is called in this section: CHARGE2
!	variables used are:
!
!		rr		:
!		nnw		:
!		ip		:
!		px		:
!****************************************************************************************
!
	subroutine flownode(x,id,rr)
!
	implicit double precision (a-h,o-z)
!
	common /a/ nnp,ifem,ibem,idyn,ieaq,nr
	common /h/ mdof,mdofn,ls,ls1,isolv
	common /c/ nload,ntime,itmax,itime,iter,iload,time,dtime,dtimei,dtimem,xinc
	common /sfn/ nnload,nnw,nna,nnw1,nna1
	common /snwm/ ip(3,100),px(7,100),infw(3,100),aqw(100),aqw1(100),infw1(3,100),&
				  infa(3,100),aqa(100),aqa1(100),infa1(3,100)
!
	dimension x(2,nnp),id(4,nnp),rr(mdof),q(3),qf(3),qw(3,100),xelm(3),yelm(3),xmm(12,12)
!
!
	if (nnw.eq.0) goto 200
!
	beta=0.67
!
	do i=1,3
		do j=1,100
			qw(i,j)=0.d0
		enddo
	enddo
!
	do 10 i=1,3
		do 10 j=1,nnw1
10			qw(i,j)=(1.-beta)*aqw1(j)
!
	do 100 i=1,nnw1
		i1=infw1(1,i)
		i2=infw1(2,i)
		i3=infw1(3,i)
		if (i2.eq.0.and.i3.eq.0) goto 30
		xelm(1)=x(1,i1)
		xelm(2)=x(1,i2)
		xelm(3)=x(1,i3)
		yelm(1)=x(2,i1)
		yelm(2)=x(2,i2)
		yelm(3)=x(2,i3)
		a=xelm(3)-2*xelm(2)+xelm(1)
		b=(xelm(3)-xelm(1))/2
		c=yelm(3)-2*yelm(2)+yelm(1)
		d=(yelm(3)-yelm(1))/2
		do k=1,3
			q(k)=0.d0
			qf(k)=0.d0
			do j=1,3	!3node*1[Pw]=3
				xmm(k,j)=0.d0
			enddo
		enddo
		qf(1)=qw(1,i)
		qf(2)=qw(2,i)
		qf(3)=qw(3,i)
		call mmatrice(xmm,xelm,yelm,1) ! kode=1; [Pw]
		do k=1,3
			do j=1,3
				q(k)=q(k)+xmm(k,j)*qf(j)
			enddo
		enddo
		j1=id(3,i1)
		j2=id(3,i2)
		j3=id(3,i3)
		if (j1.le.mdof) rr(j1)=rr(j1)-q(1)*dtime
		if (j2.le.mdof) rr(j2)=rr(j2)-q(2)*dtime
		if (j3.le.mdof) rr(j3)=rr(j3)-q(3)*dtime
		goto 100
30		j1=id(3,i1)
		if (j1.le.mdof) rr(j1)=rr(j1)-(1.-beta)*aqw1(j)*dtime
!
100	continue
!
!	reinitiate matrix of flow to zero for second part of loading
!
	do i=1,3
		do j=1,100
			qw(i,j)=0.d0
		enddo
	enddo
!
	do i=1,3
		do j=1,nnw
			qw(i,j)=beta*aqw(j)
		enddo
	enddo
!
	do 110 i=1,nnw
		i1=infw(1,i)
		i2=infw(2,i)
		i3=infw(3,i)
		if (i2.eq.0.and.i3.eq.0) goto 40
		xelm(1)=x(1,i1)
		xelm(2)=x(1,i2)
		xelm(3)=x(1,i3)
		yelm(1)=x(2,i1)
		yelm(2)=x(2,i2)
		yelm(3)=x(2,i3)
		a=xelm(3)-2*xelm(2)+xelm(1)
		b=(xelm(3)-xelm(1))/2
		c=yelm(3)-2*yelm(2)+yelm(1)
		d=(yelm(3)-yelm(1))/2
		do k=1,3
			q(k)=0.d0
			qf(k)=0.d0
			do j=1,3	!3node*1[Pw]=3
				xmm(k,j)=0.d0
			enddo
		enddo
		qf(1)=qw(1,i)
		qf(2)=qw(2,i)
		qf(3)=qw(3,i)
		call mmatrice(xmm,xelm,yelm,1) ! kode=1; [Pw]
		do k=1,3
			do j=1,3
				q(k)=q(k)+xmm(k,j)*qf(j)
			enddo
		enddo
		j1=id(3,i1)
		j2=id(3,i2)
		j3=id(3,i3)
		if (j1.le.mdof) rr(j1)=rr(j1)-q(1)*dtime
		if (j2.le.mdof) rr(j2)=rr(j2)-q(2)*dtime
		if (j3.le.mdof) rr(j3)=rr(j3)-q(3)*dtime
		goto 110
40		j1=id(3,i1)
		if (j1.le.mdof) rr(j1)=rr(j1)-beta*aqw(j)*dtime
!
110	continue
!
200	continue
!
	return
	end
!
!
!****************************************************************************************
!		           FLOANODE subroutine ****(IT IS not MODifIED FOR 3 NODED ELEMENTS)*****
!	This subroutine calculates nodal air loads due to surface air flow.
!	this subroutine is called in this section: CHARGE2
!	variables used are:
!		rr		:
!		nna		:
!		ip		:
!		px		:
!****************************************************************************************
!
	subroutine floanode(x,id,rr)
!
	implicit double precision (a-h,o-z)
!
	common /a/ nnp,ifem,ibem,idyn,ieaq,nr
	common /h/ mdof,mdofn,ls,ls1,isolv
	common /c/ nload,ntime,itmax,itime,iter,iload,time,dtime,dtimei,dtimem,xinc
	common /sfn/ nnload,nnw,nna,nnw1,nna1
	common /snwm/ ip(3,100),px(7,100),infw(3,100),aqw(100),aqw1(100),infw1(3,100),&
				  infa(3,100),aqa(100),aqa1(100),infa1(3,100)
!
	dimension x(2,nnp),id(4,nnp),rr(mdof),q(3),qf(3),qa(3,100),xelm(3),yelm(3),xmm(12,12)
!
!
!
	if (nna.eq.0) goto 200
!
	beta=0.67
!
	do i=1,3
		do j=1,100
			qa(i,j)=0.d0
		enddo
	enddo
!
	do 10 i=1,3
		do 10 j=1,nna1
10			qa(i,j)=(1.-beta)*aqa1(j)
!
	do 100 i=1,nna1
		i1=infa1(1,i)
		i2=infa1(2,i)
		i3=infa1(3,i)
		if (i2.eq.0.and.i3.eq.0) goto 30
		xelm(1)=x(1,i1)
		xelm(2)=x(1,i2)
		xelm(3)=x(1,i3)
		yelm(1)=x(2,i1)
		yelm(2)=x(2,i2)
		yelm(3)=x(2,i3)
		a=xelm(3)-2*xelm(2)+xelm(1)
		b=(xelm(3)-xelm(1))/2
		c=yelm(3)-2*yelm(2)+yelm(1)
		d=(yelm(3)-yelm(1))/2
		do k=1,3
			q(k)=0.d0
			qf(k)=0.d0
			do j=1,3	!3node*1[Pw]=3
				xmm(k,j)=0.d0
			enddo
		enddo
		qf(1)=qa(1,i)
		qf(2)=qa(2,i)
		qf(3)=qa(3,i)
		call mmatrice(xmm,xelm,yelm,1) ! kode=1; [Pw]
		do k=1,3
			do j=1,3
				q(k)=q(k)+xmm(k,j)*qf(j)
			enddo
		enddo
		j1=id(4,i1)
		j2=id(4,i2)
		j3=id(4,i3)
		if (j1.le.mdof) rr(j1)=rr(j1)-q(1)*dtime
		if (j2.le.mdof) rr(j2)=rr(j2)-q(2)*dtime
		if (j3.le.mdof) rr(j3)=rr(j3)-q(3)*dtime
		goto 100
30		j1=id(3,i1)
		if (j1.le.mdof) rr(j1)=rr(j1)-(1.-beta)*aqa1(j)*dtime
!
100	continue
!
!	reinitiate matrix of flow to zero for second part of loading
!
	do i=1,3
		do j=1,100
			qa(i,j)=0.d0
		enddo
	enddo
!
	do i=1,3
		do j=1,nna
			qa(i,j)=beta*aqa(j)
		enddo
	enddo
!
	do 110 i=1,nna
		i1=infa1(1,i)
		i2=infa1(2,i)
		i3=infa1(3,i)
		if (i2.eq.0.and.i3.eq.0) goto 40
		xelm(1)=x(1,i1)
		xelm(2)=x(1,i2)
		xelm(3)=x(1,i3)
		yelm(1)=x(2,i1)
		yelm(2)=x(2,i2)
		yelm(3)=x(2,i3)
		a=xelm(3)-2*xelm(2)+xelm(1)
		b=(xelm(3)-xelm(1))/2
		c=yelm(3)-2*yelm(2)+yelm(1)
		d=(yelm(3)-yelm(1))/2
		do k=1,3
			q(k)=0.d0
			qf(k)=0.d0
			do j=1,3	!3node*1[Pw]=3
				xmm(k,j)=0.d0
			enddo
		enddo
		qf(1)=qa(1,i)
		qf(2)=qa(2,i)
		qf(3)=qa(3,i)
		call mmatrice(xmm,xelm,yelm,1) ! kode=1; [Pw]
		do k=1,3
			do j=1,3
				q(k)=q(k)+xmm(k,j)*qf(j)
			enddo
		enddo
		j1=id(4,i1)
		j2=id(4,i2)
		j3=id(4,i3)
		if (j1.le.mdof) rr(j1)=rr(j1)-q(1)*dtime
		if (j2.le.mdof) rr(j2)=rr(j2)-q(2)*dtime
		if (j3.le.mdof) rr(j3)=rr(j3)-q(3)*dtime
		goto 110
40		j1=id(4,i1)
		if (j1.le.mdof) rr(j1)=rr(j1)-beta*aqw(j)*dtime
!
110	continue
!
200	continue
!
	return
	end
!****************************************************************************************
!
!							           CHARGE3 subroutine
!
!	This subroutine calculates the sum of all external loads (rr vector).
!	This subroutine is called in: CONTROL
!	This subroutine calls: -
!
!	variables used are:
!
!		rr		: it is a vector which saves the mechanical loadings
!		totrr	: it is a vector which saves the sum of all mechanical loads (rr vector)
!		mdof	: maximum degree of freedom by considering the boundary conditions
!
!
!	INPUT	: rr, mdof
!
!
!	OUTPUT	: totrr
!
!
!****************************************************************************************
!
	subroutine charge3(rr,totrr)
!
	implicit double precision (a-h,o-z)
!
	common /h/ mdof,mdofn,ls,ls1,isolv
!
	dimension rr(mdof),totrr(mdof)
!
!
	do i=1,mdof
		totrr(i)=totrr(i)+rr(i)
	enddo
!
	return
	end
!****************************************************************************************
!
!							           CHARGE4 subroutine
!
!	This subroutine adds the external loads caused by bem areas & fem area to r.
!	This subroutine is called in: CONTROL
!	This subroutine calls: -
!
!	variables used are:
!
!		r		: all nodal forces (mechanical TOTRR & seismic RRR) are saved in this vector
!		totrr	: it is a vector which saves the sum of all mechanical loads (rr vector)
!		rrr		: all nodal forces (mechanical FRBED3 & seismic) caused by bem areas are
!				  transfered to the rrr force vector
!		rnn		:
!		rnn1	:
!
!
!	INPUT	:
!
!	OUTPUT	: r, rnn, rnn1
!
!
!****************************************************************************************
!
	subroutine charge4(r,totrr,rrr,rnn,rnn1)
!
 	implicit double precision (a-h,o-z)
!
	common /h/ mdof,mdofn,ls,ls1,isolv
!
	dimension r(mdof),totrr(mdof),rrr(mdof),rnn(mdof),rnn1(mdof)
!
	do i=1,mdof
		r(i)=r(i)+totrr(i)+rrr(i)
        rnn(i)=rnn(i)+totrr(i)
        rnn1(i)=rnn1(i)+totrr(i)
	enddo
!
	return
	end
!****************************************************************************************
!
!							           CHARGE5 subroutine
!
!	This subroutine applies the seismic loading on D3 nodes.
!	This subroutine is called in: BEM
!	This subroutine calls: -
!
!	variables used are:
!
!		ieaq	: earthquake code [0,1,2,3];
!							[0]: static case
!							[1]: acceleration registered by accelerogram??????????????
!							[2]: imposed displacement for finite number of nodes??????
!							[3]: incident wave displacement in x- & y-directions
!		uxeq	: (incident wave) displacement in x dir
!							uxeq (3,1:nnpbed3): in current time -> itime
!							uxeq (2,1:nnpbed3): in previous time -> itime-1
!							uxeq (1,1:nnpbed3): in two previous time -> itime-2
!		uyeq	: (incident wave) displacement in y dir
!							uyeq (3,1:nnpbed3): in current time -> itime
!							uyeq (2,1:nnpbed3): in previous time -> itime-1
!							uyeq (1,1:nnpbed3): in two previous time -> itime-2
!		ueaq	: if TETA-method is used, displacement field due to incident waves transforms
!				  to "Ueq(N)-[(teta-1)/teta]*Uex(N-1)"
!		feaq	: it is equal to M*G^(-1)*Ueaq
!		kbem		: number of degree of freedom per node in BE zone;
!								if IBEM=[1] -> KBEM=[2] [X,Y]
!								if IBEM=[2] -> KBEM=[3] [X,Y,Pw]
!								if IBEM=[3] -> KBEM=[4] [X,Y,Pw,Pa]
!		xgid3		: M*G^(-1)
!		gid3		: inverse matrice of gd3
!
!
!	INPUT	:
!
!
!	OUTPUT	:
!
!
!****************************************************************************************
!
	subroutine charge5(feaq,xgid3,ie3d3)
!
	implicit double precision (a-h,o-z)
!
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /b/ m8d,m8c,m8u,m8dep(10),m8cep(10),m8uep(10)
	common /chrg3/ nstime,neaq,nneaq(500),uxeaq,uyeaq,axeaq,ayeaq,uxeq(3,500),uyeq(3,500)
	common /stability/ anew1,anew2,wil,intb		! this common is used in DIM1
!
	dimension feaq(n12gh),xgid3(n12gh,n12gh),ueaq(n12gh),ie3d3(npbmax)
!
!
	n=kbem*nnpbed3
!
	do i=1,nnpbed3
		ueaq(kbem*(i-1)+1)=uxeq(3,i) -(wil-1.d0)/wil*uxeq(2,i)
        ueaq(kbem*(i-1)+2)=uyeq(3,i) -(wil-1.d0)/wil*uyeq(2,i)
		if (kbem.eq.3) ueaq(kbem*i)=0.d0
		if (kbem.eq.4) then
			ueaq(kbem*i-1)=0.d0
			ueaq(kbem*i)  =0.d0
		endif
	enddo
!
	do i=1,n
		feaq(i)=0.d0
		do j=1,n
			feaq(i) = feaq(i) + xgid3(i,j)*ueaq(j)
		enddo
	enddo
!
	return
	end
!****************************************************************************************
!
!							           CHARGE6 subroutine
!
!	This subroutine sums all nodal forces (internal & seismic) caused by bem areas and
!	transfers the sum to the rrr force vector.
!	This subroutine is called in: BEM
!	This subroutine calls: -
!
!	variables used are:
!
!		id		: degrees of freedom of each node from "1" to "mdof"; for the nodes with
!				  the boundary conditions, each node takes from "mdof+1" to "mdof+nbc"
!				  [X,Y]; which are free for the dry elm &
!				  [X,Y,Pw]; which are free for the saturated elm
!				  [X,Y,Pw,Pa]; which are free for the unsaturated elm
!		ie3d1	: ie3d1(1:nnpbed1,1:nbed): nodes' number which form the current D1 BEM area
!										   Numbering direction for external surface close
!										   (finite) domain is ANTICLOCKWISE
!				  ie3d1(nnpbed1(i)+1,1:nbed): it shows the type of just one material which
!											  is used:
!							[1]: dry boundary elm
!							[2]: saturated boundary elm
!							[3]: unsaturated boundary elm
!		ie3d2	: ie3d2(1:nnpbed2,1:nbed): nodes' number which form the current D2 BEM area
!										   Numbering direction for internal surface open
!										   (infinite) domain is CLOCKWISE
!				  ie3d2(nnpbed2(i)+1,1:nbed): it shows the type of just one material which
!											  is used:
!							[1]: dry boundary elm
!							[2]: saturated boundary elm
!							[3]: unsaturated boundary elm
!		ie3d3	: ie3d3(1:nnpbed3): nodes' number which form the current D3 BEM area
!										   Numbering direction for semi-infinite domain is
!										   ANTICLOCKWISE
!				  ie3d3(nnpbed3+1): it shows the type of just one material which
!											  is used:
!							[1]: dry boundary elm
!							[2]: saturated boundary elm
!							[3]: unsaturated boundary elm
!		feaq	: it is equal to M*G^(-1)*Ueaq
!		frbed1	: nodal internal forces in D1 BEM area (finite zone)
!		frbed2	: nodal internal forces in D2 BEM area (infinite zone)
!		frbed3	: nodal internal forces in D3 BEM area (semi-infinite zone)
!		rrr		: all nodal forces (mechanical & seismic) caused by bem areas are transfered
!				  to the rrr force vector
!		kbem		: number of degree of freedom per node in BE zone;
!								if IBEM=[1] -> KBEM=[2] [X,Y]
!								if IBEM=[2] -> KBEM=[3] [X,Y,Pw]
!								if IBEM=[3] -> KBEM=[4] [X,Y,Pw,Pa]
!		nbed1	: number of D1 BEM areas (finite areas), it can be at maximum 5 zones
!		nbed2	: number of D2 BEM areas (infinite areas), it can be at maximum 5 zones
!		nbed3	: number of D3 BEM areas (semi-infinite area), it can be at maximum 1 zone
!		nnpbed1	: number of nodes of each D1 BEM area; NNPBED1(1 to NBED1)
!		nnpbed2	: number of nodes of each D2 BEM area; NNPBED2(1 to NBED2)
!		nnpbed3	: number of nodes of one D3 BEM area
!		ide		: id matrix for BE zone which put in verctorial format
!		jdf		: global doF number of a node (1:nnpbed1) with the local dof number kbem
!		ndfe	: number of degree of freedom for BE zone
!		mdof	: maximum degree of freedom by considering the boundary conditions
!
!
!	INPUT	:
!
!
!	OUTPUT	:
!
!
!****************************************************************************************
!
	subroutine charge6(id,ie3d1,ie3d2,ie3d3,rrr,frbed1,frbed2,frbed3,feaq)
!
 	implicit double precision (a-h,o-z)
!
	common /a/ nnp,ifem,ibem,idyn,ieaq,nr
	common /a2/ nbed1,nbed2,nbed3,nnpbed1(5),nnpbed2(5),nnpbed3,nbeint,nnpen,&
				npbmax,n12gh,n3gh,kbem,kfsol,ne3u
	common /a3/ ne8d1,ne8c1,ne8u1,ne3u1,nload1,nbcx1,nbcy1,nbcw1,nbca1,m8d1,m8c1,m8u1,&
				nbed11,nbed21,nbed31,nbeint1,ne8dd,ne8cc,ne8uu
	common /b/ m8d,m8c,m8u,m8dep(10),m8cep(10),m8uep(10)
	common /h/ mdof,mdofn,ls,ls1,isolv
!
	dimension id(4,nnp),ie3d1(npbmax,nbed11),ie3d2(npbmax,nbed21),ie3d3(npbmax+1),rrr(mdof),&
			  frbed1(n12gh,nbed11),frbed2(n12gh,nbed21),frbed3(n12gh),ide(n12gh),jdf(4),&
			  feaq(n12gh)
!
!
!	---------------------------- bem area d1 ----------------------------
!
	if (nbed1.eq.0) goto 100
	do 90 i=1,nbed1
		ndfe=kbem*nnpbed1(i)
		do j=1,nnpbed1(i)
			ii=ie3d1(j,i)
			do k=1,kbem
				jdf(k)=kbem*j-kbem+k
				ide(jdf(k))=id(k,ii)
			enddo
		enddo
		do j=1,ndfe
			ii=ide(j)
			if (ii.le.mdof) rrr(ii)=rrr(ii)+frbed1(j,i)
		enddo
90	continue
!
!
!	---------------------------- bem area d2 ----------------------------
!
100	if (nbed2.eq.0) goto 200
	do 190 i=1,nbed2
		ndfe=kbem*nnpbed2(i)
		do j=1,nnpbed2(i)
			ii=ie3d2(j,i)
			do k=1,kbem
				jdf(k)=kbem*j-kbem+k
				ide(jdf(k))=id(k,ii)
			enddo
		enddo
		do j=1,ndfe
			ii=ide(j)
			if (ii.le.mdof) rrr(ii)=rrr(ii)+frbed2(j,i)
		enddo
190	continue
!
!
!	---------------------------- bem area d3 ----------------------------
!
200	if (nbed3.eq.0) goto 300
	ndfe=kbem*nnpbed3
	do j=1,nnpbed3
		ii=ie3d3(j)
		do k=1,kbem
			jdf(k)=kbem*j-kbem+k
			ide(jdf(k))=id(k,ii)
		enddo
	enddo
	do j=1,ndfe
		ii=ide(j)
		if (ii.le.mdof) then
			rrr(ii)=rrr(ii)+frbed3(j)
			if (ieaq.eq.3) rrr(ii)=rrr(ii)+feaq(j)
		endif
	enddo
!
300	return
	end

