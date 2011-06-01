      subroutine mexfunction(nlhs, plhs, nrhs, prhs)

%
%   Matlab wrapper for NEOART, to be called as:
%
%   >> [c1,c2,c3,c4] = neoart(ns,nc,M,T,zsp,den,ds,isel,rho{,RN,BN,q},ic,eparr{,eps,nreg,sigma,nleg,nenergy,ncof,neogeo,neofrc})
%
%   Inputs:
%	ns	number of species
%	nc	number of charge state per species, dim(1,ns)
%	M	species mass in kg, dim(1,ns)
%	T	species temperature in keV, dim(1,ns)
%	zsp	charge state of each component, dim(ns,max(nc))
%	den	density of each component, dim(ns,max(nc))
%	ds	first and second thermondynamic forces (pressure and temperature gradient) 
%		ds(i,j,1) = -d ln(p_i_j)/drho and ds(i,j,2) = -d ln(T_i)/drho
%		rho is the flux label actually used, ds>0 for a peaked profile
%	isel	select the geometry treatment
%		2 -> uses circular flux surfaces, user must specify RN, BN and q. The flux label is rho=(Rmax-Rmin)/2
%		4 -> reads the file neoart_geom.dat generated with CHEASE using NIDEAL=10. The flux label is rho=(Rmax-Rmin)/2/R0EXP, where R0EXP is specified in the CHEASE input file.
%	rho	flux surface label
%	RN	major radius [m] of the centre of the flux surface, for isel=2
%	BN	magnetic field [T] at RN, for isel=2
%	q	safety factor, for isel=2
%	ic	contribution for which the fluxes are calculated:
%			0 -> classical particle flux
%			1 -> banana plateau contribution
%			2 -> Pfirsch-Schlueter contribution
%			3 -> all together
%	eparr	the loop voltage [V] divided by 2pi
%	eps	accuracy [default eps=1e-5]
%	nreg	to force a regime in the viscosity calculation [default nreg=0]
%			0 -> banana and Pfirsch-Schlueter regimes are weighted
%			1 -> banana regime is forced
%			2 -> Pfirsch-Schlueter regime is forced
%	sigma	select the coupling terms kept in the Pfirsch-Schlueter equation
%			sigma(1)=0,1 -> coupling between heat eq. and the eq. for the third Laguerre harmonic
%			sigma(2)=0,1 -> coupling between the eq. for the third Laguerre harmonic and the eq. for the temperature pert.
%			sigma(3)=0,1 -> energy exchange between the different species
%			sigma(4)=0,1 -> coupling between the eq. for the third Laguerre harmonic and the eq. for the density pert.
%	nleg	number of Legendre polynomials in the expansion [default nleg=3]
%	nenergy	select or not the energy in the calculation of the viscosity [default nenergy=1]
%	ncof	select or not ion=electron collisions [default ncof=1]
%	neogeo	if 1 recalculates the geometry dep. parameters [default neogeo=1]
%	neofrc	if 0 recalculates the friction and viscosity matrixes [default neofrc=0]
%
%   The number of inputs must be 22 or 14 for isel=2 and  19 or 11 for isel=4
%
%	Outputs:
%	c1	particle flux, dim(ns,max(nc)) positive values for outward flux
%		c1(i,j) = <Gamma_ij . grad_rho>, with Gamma_ij in 1/m^2/s and <.> the FS average 
%	c2	heat flux, dim(ns,max(nc)) positive values for outward flux
%		c2(i,j) = <Q_ij . grad_rho>/T_i with Q_ij/T_i in 1/m^2/s
%	c3	parallel flow, dim(ns,max(nc))
%		c3(i,j) = q_ij*n_ij*<u_ij . B><B>/<B^2> with q_ij*n_ij*u_ij the current density of species i,j in A/m^2
%	c4	poloidal flow, dim(ns,max(nc))
%		c4(i,j) = <u_ij.grad_theta>/<B.grad_theta> in m/s/T with u_ij the flow of species i,j in m/s

