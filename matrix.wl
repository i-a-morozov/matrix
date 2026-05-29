(* Transport matrix computation *)
(* I.M., 2025-2026 *)

(* --------- RadFldPtcTrj --------- *)

ClearAll[RadFldPtcTrj] ;
Options[RadFldPtcTrj] = {
	Method -> {"FixedStep", Method -> {"StiffnessSwitching", Method -> {"ExplicitRungeKutta", Automatic}}}
} ;
RadFldPtcTrj::usage="RadFldPtcTrj[obj, E, {x0, dxdy0, z0, dzdy0}, {y0, y1}, np] -- compute transverse coordinates and its derivatives (angles/slopes) of a relativistic charged particle trajectory in 3D magnetic field produced by the object obj, using the NDSolve interface (see default method options). The particle energy is E [GeV], initial transverse coordinates and derivatives are {x0, dxdy0, z0, dzdy0}; the longitudinal coordinate y is varied from y0 to y1 in np steps. All positions are in millimeters and angles in radians." ;
RadFldPtcTrj[obj_, E_, {x0_, dxdy0_, z0_, dzdy0_}, {y0_, y1_}, np_, options:OptionsPattern[]] := Block[
	{alpha, bx, by, bz, system, x, xp, z, zp, xi, xpi, zi, zpi, y, solver, table, functions, positions},
	alpha = (0.299792458/E)/1000.0 ;
	bx[x_?NumericQ, y_?NumericQ, z_?NumericQ] := bx[x, y, z] = radFld[obj, "bx", {x, y, z}] ;
	by[x_?NumericQ, y_?NumericQ, z_?NumericQ] := by[x, y, z] = radFld[obj, "by", {x, y, z}] ;
	bz[x_?NumericQ, y_?NumericQ, z_?NumericQ] := bz[x, y, z] = radFld[obj, "bz", {x, y, z}] ;
	system = {
		x'[y] == xp[y],
		z'[y] == zp[y],
		xp'[y] == -alpha/Sqrt[1 + xp[y]^2 + zp[y]^2]*(zp[y]*by[x[y], y, z[y]] + (1 + xp[y]^2)*bz[x[y], y, z[y]] + xp[y]*zp[y]*bx[x[y], y, z[y]]),
		zp'[y] == +alpha/Sqrt[1 + xp[y]^2 + zp[y]^2]*(xp[y]*by[x[y], y, z[y]] + (1 + zp[y]^2)*bx[x[y], y, z[y]] + xp[y]*zp[y]*bz[x[y], y, z[y]]),
		x[y0] == xi,
		xp[y0] == xpi,
		z[y0] == zi,
		zp[y0] == zpi
	} ;
	solver = ParametricNDSolveValue[
		system,
		{x, xp, z, zp},
		{y, y0, y1},
		{xi, xpi, zi, zpi},
		Method -> OptionValue[Method],
		StartingStepSize -> (y1 - y0)/np,
		MaxStepSize  -> (y1 - y0)/np
	] ;
	functions = solver[x0, dxdy0, z0, dzdy0] ;
	positions = Subdivide[y0, y1, np - 1] ;
	Map[Flatten, Transpose[{positions, Map[Function[{position}, Through[functions[position]]], positions]}]]
] ;

(* --------- RadFldPtcTrjCnn --------- *)

ClearAll[RadFldPtcTrjCnn];
Options[RadFldPtcTrjCnn] = {
	"Delta" -> {0.1, 0.1},
	Method -> {"FixedStep", Method -> {"StiffnessSwitching", Method -> {"ExplicitRungeKutta", Automatic}}}
} ;
RadFldPtcTrjCnn::usage="RadFldPtcTrjCnn[obj, E, dE, {qx0, px0, qz0, pz0}, {y0, y1}, np, {extAx, extAy, extAz}] -- compute transverse canonical coordinates of a relativistic charged particle trajectory in 3D magnetic field produced by the object obj, using the NDSolve interface (see default options). The particle energy is E [GeV], initial transverse canonical coordinates are {qx0, px0, qz0, pz0}; the longitudinal coordinate y is varied from y0 to y1 in np steps. All positions are in millimeters." ;
RadFldPtcTrjCnn[obj_, E_, dE_, {qx0_, px0_, qz0_, pz0_}, {y0_, y1_}, np_, {extAx_, extAy_, extAz_}, options : OptionsPattern[]] := Block[
	{alpha, Ax, Ay, Az, ax, ay, az, daxdx, daydx, dazdx, daxdz, daydz, dazdz, hx, hz, flow, system, y, qx, qz, px, pz, qxi, qzi, pxi, pzi, solver, functions, positions},
	alpha = -(0.299792458/E)/1000.0 ;
	Ax[x_?NumericQ, y_?NumericQ, z_?NumericQ] := Ax[x, y, z] = radFld[obj, "ax", {x, y, z}] ;
	Ay[x_?NumericQ, y_?NumericQ, z_?NumericQ] := Ay[x, y, z] = radFld[obj, "ay", {x, y, z}] ;
	Az[x_?NumericQ, y_?NumericQ, z_?NumericQ] := Az[x, y, z] = radFld[obj, "az", {x, y, z}] ;
	ax[x_?NumericQ, y_?NumericQ, z_?NumericQ] := alpha*Ax[x, y, z] + extAx[x, y, z] ;
	ay[x_?NumericQ, y_?NumericQ, z_?NumericQ] := alpha*Ay[x, y, z] + extAy[x, y, z] ;
	az[x_?NumericQ, y_?NumericQ, z_?NumericQ] := alpha*Az[x, y, z] + extAz[x, y, z] ;
	daxdx[x_?NumericQ, y_?NumericQ, z_?NumericQ] := daxdx[x, y, z] = (ax[x + hx, y, z] - ax[x - hx, y, z])/(2 hx) ;
	daydx[x_?NumericQ, y_?NumericQ, z_?NumericQ] := daydx[x, y, z] = (ay[x + hx, y, z] - ay[x - hx, y, z])/(2 hx) ;
	dazdx[x_?NumericQ, y_?NumericQ, z_?NumericQ] := dazdx[x, y, z] = (az[x + hx, y, z] - az[x - hx, y, z])/(2 hx) ;
	daxdz[x_?NumericQ, y_?NumericQ, z_?NumericQ] := daxdz[x, y, z] = (ax[x, y, z + hz] - ax[x, y, z - hz])/(2 hz) ;
	daydz[x_?NumericQ, y_?NumericQ, z_?NumericQ] := daydz[x, y, z] = (ay[x, y, z + hz] - ay[x, y, z - hz])/(2 hz) ;
	dazdz[x_?NumericQ, y_?NumericQ, z_?NumericQ] := dazdz[x, y, z] = (az[x, y, z + hz] - az[x, y, z - hz])/(2 hz) ;
	{hx, hz} = OptionValue["Delta"] ;
	flow[qx_, qz_, px_, pz_, y_] := flow[qx, qz, px, pz, y] = Block[
		{axi, ayi, azi, pix, piz, sqrt, daxdxi, daydxi, dazdxi, daxdzi, daydzi, dazdzi, qxdot, qzdot, pxdot, pzdot},
		axi = ax[qx, y, qz] ;
		ayi = ay[qx, y, qz] ;
		azi = az[qx, y, qz] ;
		pix = px - axi ;
		piz = pz - azi ;
		sqrt = Sqrt[(1 + dE)^2 - pix^2 - piz^2] ;
		daxdxi = daxdx[qx, y, qz] ;
		daydxi = daydx[qx, y, qz] ;
		dazdxi = dazdx[qx, y, qz] ;
		daxdzi = daxdz[qx, y, qz] ;
		daydzi = daydz[qx, y, qz] ;
		dazdzi = dazdz[qx, y, qz] ;
		qxdot = pix/sqrt ;
		qzdot = piz/sqrt ;
		pxdot = daydxi + (pix*daxdxi + piz*dazdxi)/sqrt ;
		pzdot = daydzi + (pix*daxdzi + piz*dazdzi)/sqrt ;
		{qxdot, qzdot, pxdot, pzdot}
	] ;
	system = {
		{qx'[y], qz'[y], px'[y], pz'[y]} == flow[qx[y], qz[y], px[y], pz[y], y],
		qx[y0] == qxi,
		px[y0] == pxi,
		qz[y0] == qzi,
		pz[y0] == pzi
	} ;
	solver = ParametricNDSolveValue[
		system,
		{qx, px, qz, pz},
		{y, y0, y1},
		{qxi, pxi, qzi, pzi},
		Method -> OptionValue[Method],
		StartingStepSize -> (y1 - y0)/np,
		MaxStepSize -> (y1 - y0)/np
	] ;
	functions = solver[qx0, px0, qz0, pz0] ;
	positions = Subdivide[y0, y1, np - 1] ;
	Map[Flatten, Transpose[{positions, Map[Function[{position}, Through[functions[position]]], positions]}]]
] ;

(* --------- transport --------- *)

ClearAll[transport] ;
Options[transport] = {
	Method -> {"FixedStep", Method -> {"StiffnessSwitching", Method -> {"ExplicitRungeKutta", Automatic}}},
	"Thin" -> True
} ;
transport::usage = "transport[object, energy, delta, {start, stop}, steps, angles, solver][{qx, px, qz, pz}] -- track canonical initial condition (qx, px, qz, pz) with (qx, qz in m) through slope-based solver (x and z in mm)" ;
transport[
	object_,                  (* -- magnet object *)
	energy_,                  (* -- reference energy (GeV) *)
	delta_,                   (* -- energy deviation *)
	{start_, stop_},          (* -- integration start and stop positions relative to the object *)
	steps_,                   (* -- number of integration steps *)
	angles_:{0, 0, 0, 0},     (* -- angle kicks to be added to canonical momenta on entrance and exit (cxi, czi, cxf, czf) *)
	solver_:radFldPtcTrj,     (* -- solver (radFldPtcTrj or RadFldPtcTrj) *)
	options:OptionsPattern[]  (* -- options *)
][state_] := Block[
	{QX, PX, QZ, PZ, X, XP, Z, ZP, CXI, CZI, CXF, CZF},
	{QX, PX, QZ, PZ} = state ;
	{CXI, CZI, CXF, CZF} = angles ;
	If[
		OptionValue["Thin"],
		{QX, PX, QZ, PZ} = {QX - 1/2*(stop - start)/1000.0*PX/Sqrt[(1 + delta)^2 - PX^2 - PZ^2], PX, QZ - 1/2*(stop - start)/1000.0*PZ/Sqrt[(1 + delta)^2 - PX^2 - PZ^2], PZ} ;
	] ;
	{QX, PX, QZ, PZ} = {QX, PX + CXI, QZ, PZ + CZI} ;
    {X, XP, Z, ZP} = {QX, PX/Sqrt[(1 + delta)^2 - PX^2 - PZ^2], QZ, PZ/Sqrt[(1 + delta)^2 - PX^2 - PZ^2]} ;
    {X, XP, Z, ZP} = {1000.0*X, XP, 1000.0*Z, ZP} ;
    {X, XP, Z, ZP} = Rest[Last[solver[object, energy*(1 + delta), {X, XP, Z, ZP}, {start, stop}, steps, Sequence @@ FilterRules[options, Options[NDSolve]]]]] ;
    {X, XP, Z, ZP} = {X/1000.0, XP, Z/1000.0, ZP} ;
    {QX, PX, QZ, PZ} = {X, (1 + delta)*XP/Sqrt[1 + XP^2 + ZP^2], Z, (1 + delta)*ZP/Sqrt[1 + XP^2 + ZP^2]} ;
    {QX, PX, QZ, PZ} = {QX, PX + CXF, QZ, PZ + CZF} ;
	If[
		OptionValue["Thin"],
		{QX, PX, QZ, PZ} = {QX - 1/2*(stop - start)/1000.0*PX/Sqrt[(1 + delta)^2 - PX^2 - PZ^2], PX, QZ - 1/2*(stop - start)/1000.0*PZ/Sqrt[(1 + delta)^2 - PX^2 - PZ^2], PZ} ;
	] ;
    {QX, PX, QZ, PZ}
] ;

(* --------- amplitude --------- *)

ClearAll[amplitude] ;
Options[amplitude] = {"SamplesPerHarmonic" -> 32, "Samples" -> Automatic} ;
amplitude::usage = "amplitude[object, component, period, reference, harmonic, options] -- compute the amplitude of the n-th fourier harmonic of given periodic field component over one longitudinal period centered at given reference point using trapezoidal rule " ;
amplitude[                    (* -- amplitude *)
    object_,                  (* -- radia object *)
    component_,               (* -- field component *)
    period_,                  (* -- longitudinal period (mm) *)
    reference_,               (* -- reference point (mm) *)
    harmonic_,                (* -- harmonic number *)
    options:OptionsPattern[]  (* -- options *)
] := Block[
    {count, fields, weights, range, angles, cos, sin},
    count = If[OptionValue["Samples"] === Automatic, harmonic*OptionValue["SamplesPerHarmonic"], OptionValue["Samples"]] ;
    fields = radFldLst[object, component, reference - {0, period/2, 0}, reference + {0, period/2, 0}, count] ;
    weights = ConstantArray[1.0, count] ;
    weights[[+1]] = 0.5 ;
    weights[[-1]] = 0.5 ;
    range = Range[0, count - 1] ;
    angles = 2*Pi*harmonic*range/(count - 1) ;
    cos = (2/(count - 1))*Total[weights*fields*Cos[angles]] ;
    sin = (2/(count - 1))*Total[weights*fields*Sin[angles]] ;
    Sqrt[cos^2 + sin^2]
] ;

(* --------- Elleaume potential (periodic field) --------- *)

ClearAll[potential] ;
Options[potential] = Join[{"SamplesPerHarmonic" -> 32, "Samples" -> Automatic, NIntegrate -> False}, Options[NIntegrate]] ;
potential::usage = "potential[object, {x, z}, periods, harmonics, shift, options] -- compute Elleaume one-(super)period potential. By default sampled Fourier amplitudes are used; use NIntegrate -> True to compute the full period integral directly from shift - period/2 to shift + period/2 without harmonic expansion." ;
potential[                    (* -- potential (T^2 mm^3) *)
    object_,                  (* -- radia object *)
    {x_, z_},                 (* -- transverse evaluation point (mm) *)
    periods_,                 (* -- horizontal and vertical periods (mm), {ph, pv} = {n*p, p} or {p, n*p} and n*p -- super-period *)
    harmonics_,               (* -- list of harmonics *)
    shift_,                   (* -- longitudinal shift/position (mm) *)
    options:OptionsPattern[]  (* -- options *)
] := Block[{period, hx, hz, start, stop, integrate, field, primitive, value, bx, bz},
  period = Max[periods] ;
  {hx, hz} = Round[period/periods] ;
  If[
    OptionValue[NIntegrate],
    start = shift - period/2 ;
    stop = shift + period/2 ;
    integrate[integrand_, range_] := Apply[NIntegrate, Join[{integrand, range}, FilterRules[{options}, Options[NIntegrate]]]] ;
    field[component_, y_?NumericQ] := field[component, y] = radFld[object, component, {x, y, z}] ;
    primitive[component_, y_?NumericQ] := primitive[component, y] = integrate[field[component, s], {s, start, y}] ;
    value[component_] := value[component] = integrate[primitive[component, y], {y, start, stop}]/period ;
    integrate[(primitive["bx", y] - value["bx"])^2 + (primitive["bz", y] - value["bz"])^2, {y, start, stop}],
    bx = Table[amplitude[object, "bx", period, {x, shift, z}, hx*harmonic, Sequence @@ FilterRules[{options}, Options[amplitude]]], {harmonic, harmonics}] ;
    bz = Table[amplitude[object, "bz", period, {x, shift, z}, hz*harmonic, Sequence @@ FilterRules[{options}, Options[amplitude]]], {harmonic, harmonics}] ;
    0.5*period*(period/(2*Pi))^2*Total[(bx^2/hx^2 + bz^2/hz^2)/harmonics^2]
  ]
] ;

(* --------- horizontal slope kick (period) --------- *)

ClearAll[dxp] ;
Options[potential] = Join[{"SamplesPerHarmonic" -> 32, "Samples" -> Automatic, NIntegrate -> False}, Options[NIntegrate]] ;
dxp::usage = "dxp[object, {x, z}, periods, harmonics, shift, energy, delta, options] -- compute one-(super) period horizontal angle kick (murad) using central finite difference" ;
dxp[                          (* -- one-(super) period horizontal angle kick (murad) *)
    object_,                  (* -- radia object *)
    {x_, z_},                 (* -- transverse evaluation point (mm) *)
    periods_,                 (* -- horizontal and vertical periods (mm), {ph, pv} = {n*p, p} or {p, n*p} and n*p -- super-period *)
    harmonics_,               (* -- list of harmonics *)
    shift_,                   (* -- longitudinal shift/position (mm) *)
    energy_,                  (* -- reference energy (GeV) *)
    delta_,                   (* -- finite difference delta (mm) *)
    options:OptionsPattern[]  (* -- options *)
] := Block[{pa, pb},
    pa = potential[object, {x - delta/2, z}, periods, harmonics, shift, Sequence @@ FilterRules[{options}, Options[potential]]] ;
    pb = potential[object, {x + delta/2, z}, periods, harmonics, shift, Sequence @@ FilterRules[{options}, Options[potential]]] ;
    (pb - pa)/delta*0.5*(0.299792458/energy)^2
] ;

(* --------- vertical slope kick (period) --------- *)

ClearAll[dzp] ;
Options[potential] = Join[{"SamplesPerHarmonic" -> 32, "Samples" -> Automatic, NIntegrate -> False}, Options[NIntegrate]] ;
dzp::usage = "dzp[object, {x, z}, periods, harmonics, shift, energy, delta, options] -- compute one-(super) period vertical angle kick (murad) using central finite difference" ;
dzp[                          (* -- one-(super) period vertical angle kick (murad) *)
    object_,                  (* -- radia object *)
    {x_, z_},                 (* -- transverse evaluation point (mm) *)
    periods_,                 (* -- horizontal and vertical periods (mm), {ph, pv} = {n*p, p} or {p, n*p} and n*p -- super-period *)
    harmonics_,               (* -- list of harmonics *)
    shift_,                   (* -- longitudinal shift/position (mm) *)
    energy_,                  (* -- reference energy (GeV) *)
    delta_,                   (* -- finite difference delta (mm) *)
    options:OptionsPattern[]  (* -- options *)
] := Block[{pa, pb},
    pa = potential[object, {x, z - delta/2}, periods, harmonics, shift, Sequence @@ FilterRules[{options}, Options[potential]]] ;
    pb = potential[object, {x, z + delta/2}, periods, harmonics, shift, Sequence @@ FilterRules[{options}, Options[potential]]] ;
    (pb - pa)/delta*0.5*(0.299792458/energy)^2
] ;

(* --------- horizontal focusing strength (period) --------- *)

ClearAll[kx] ;
Options[potential] = Join[{"SamplesPerHarmonic" -> 32, "Samples" -> Automatic, NIntegrate -> False}, Options[NIntegrate]] ;
kx::usage = "kx[object, {x, z}, periods, harmonics, shift, energy, delta, options] -- compute horizontal focusing strength using central finite difference" ;
kx[                           (* -- horizontal focusing strength (1/m) *)
    object_,                  (* -- radia object *)
    {x_, z_},                 (* -- transverse evaluation point (mm) *)
    periods_,                 (* -- horizontal and vertical periods (mm), {ph, pv} = {n*p, p} or {p, n*p} and n*p -- super-period *)
    harmonics_,               (* -- list of harmonics *)
    shift_,                   (* -- longitudinal shift/position (mm) *)
    energy_,                  (* -- reference energy (GeV) *)
    delta_,                   (* -- finite difference delta (mm) *)
    options:OptionsPattern[]  (* -- options *)
] := Block[{pa, pb, pc},
    pa = potential[object, {x - delta/2, z}, periods, harmonics, shift, Sequence @@ FilterRules[{options}, Options[potential]]] ;
    pb = potential[object, {x, z}, periods, harmonics, shift, Sequence @@ FilterRules[{options}, Options[potential]]] ;
    pc = potential[object, {x + delta/2, z}, periods, harmonics, shift, Sequence @@ FilterRules[{options}, Options[potential]]] ;
    10.0^-3*4*(pa - 2*pb + pc)/delta^2*0.5*(0.299792458/energy)^2
] ;

(* --------- vertical focusing strength (period) --------- *)

ClearAll[kz] ;
Options[potential] = Join[{"SamplesPerHarmonic" -> 32, "Samples" -> Automatic, NIntegrate -> False}, Options[NIntegrate]] ;
kz::usage = "kz[object, {x, z}, periods, harmonics, shift, energy, delta, options] -- compute vertical focusing strength using central finite difference" ;
kz[                           (* -- vertical focusing strength (1/m) *)
    object_,                  (* -- radia object *)
    {x_, z_},                 (* -- transverse evaluation point (mm) *)
    periods_,                 (* -- horizontal and vertical periods (mm), {ph, pv} = {n*p, p} or {p, n*p} and n*p -- super-period *)
    harmonics_,               (* -- list of harmonics *)
    shift_,                   (* -- longitudinal shift/position (mm) *)
    energy_,                  (* -- reference energy (GeV) *)
    delta_,                   (* -- finite difference delta (mm) *)
    options:OptionsPattern[]  (* -- options *)
] := Block[{pa, pb, pc},
    pa = potential[object, {x, z - delta/2}, periods, harmonics, shift, Sequence @@ FilterRules[{options}, Options[potential]]] ;
    pb = potential[object, {x, z}, periods, harmonics, shift, Sequence @@ FilterRules[{options}, Options[potential]]] ;
    pc = potential[object, {x, z + delta/2}, periods, harmonics, shift, Sequence @@ FilterRules[{options}, Options[potential]]] ;
    10.0^-3*4*(pa - 2*pb + pc)/delta^2*0.5*(0.299792458/energy)^2
] ;

(* --------- dkd potential based tracking --------- *)

ClearAll[dkd] ;
Options[potential] = Join[{"SamplesPerHarmonic" -> 32, "Samples" -> Automatic, NIntegrate -> False}, Options[NIntegrate]] ;
dkd::usage = "dkd[object, energy, delta, periods, harmonics, shift, step, count, factors, options][{qx, px, qz, pz}] -- drift-kick-drift canonical tracking (period in mm, x and z in m)" ;
dkd[                          (* -- drift-kick-drift canonical tracking *)
	object_,                  (* -- radia object *)
	energy_,                  (* -- reference energy (GeV) *)
	delta_,                   (* -- energy delta *)
	periods_,                 (* -- horizontal and vertical periods (mm), {ph, pv} = {n*p, p} or {p, n*p} and n*p -- super-period *)
	harmonics_,               (* -- list of harmonics *)
	shift_,                   (* -- longitudinal shift/position (mm) *)
	step_,                    (* -- finite difference delta (mm) *)
	count_,                   (* -- total number of (super) periods *)
	factors_:{1.0, 1.0},      (* -- extra kick multiplication factors *)
    options:OptionsPattern[]  (* -- options *)
][state_] := Block[
    {FX, FZ, QX, PX, QZ, PZ, X, XP, Z, ZP, DL},
    {FX, FZ} = factors ;
    {QX, PX, QZ, PZ} = state ;
    DL = Max[periods]/2/10.0^3 ;
    {QX, PX, QZ, PZ} = {QX - count*DL*PX/Sqrt[(1 + delta)^2 - PX^2 - PZ^2], PX, QZ - count*DL*PZ/Sqrt[(1 + delta)^2 - PX^2 - PZ^2], PZ} ;
    Do[
		{QX, PX, QZ, PZ} = {QX + DL*PX/Sqrt[(1 + delta)^2 - PX^2 - PZ^2], PX, QZ + DL*PZ/Sqrt[(1 + delta)^2 - PX^2 - PZ^2], PZ} ;
		{X, XP, Z, ZP} = {QX, PX/Sqrt[(1 + delta)^2 - PX^2 - PZ^2], QZ, PZ/Sqrt[(1 + delta)^2 - PX^2 - PZ^2]} ;
		{X, XP, Z, ZP} = {
			X,
			XP - FX*dxp[object, 10.0^3*{X, Z}, periods, harmonics, shift, energy*(1 + delta), step, Sequence @@ FilterRules[{options}, Options[potential]]]*10^-6,
			Z,
			ZP - FZ*dzp[object, 10.0^3*{X, Z}, periods, harmonics, shift, energy*(1 + delta), step, Sequence @@ FilterRules[{options}, Options[potential]]]*10^-6
		} ;
		{QX, PX, QZ, PZ} = {X, (1 + delta)*XP/Sqrt[1 + XP^2 + ZP^2], Z, (1 + delta)*ZP/Sqrt[1 + XP^2 + ZP^2]} ;
		{QX, PX, QZ, PZ} = {QX + DL*PX/Sqrt[(1 + delta)^2 - PX^2 - PZ^2], PX, QZ + DL*PZ/Sqrt[(1 + delta)^2 - PX^2 - PZ^2], PZ},
		count
	] ;
    {QX, PX, QZ, PZ} = {QX - count*DL*PX/Sqrt[(1 + delta)^2 - PX^2 - PZ^2], PX, QZ - count*DL*PZ/Sqrt[(1 + delta)^2 - PX^2 - PZ^2], PZ} ;
	{QX, PX, QZ, PZ}
] ;

(* --------- explicit ID transport matrix (approximate) --------- *)

ClearAll[idtm] ;
idtm::usage = "idtm[{kx, kz}, {np, lp}, dp] -- compute id thin insertion exponent diagonal and corresponding transport matrix (second order in kx and kz)" ;
idtm[                         (* -- id thin insertion diagonal and transport matrix *)
	{kx_, kz_},               (* -- focusing strength (1/m) *)
	count_,                   (* -- total number of (super) periods *)
	period_,                  (* -- (super) period length (m) *)
	delta_                    (* -- energy delta *)
] := Block[
	{a, b, c, d, diagonal, matrix},
	a = (kx*count)/(1 + delta) - (kx^2*period*count*(-1 + count^2))/(6*(1 + delta)^3) ;
	b = (kx*period^2*count*(-1 + count^2))/(12*(1 + delta)^3) + (kx^2*period^3*count*(-1 + count^4))/(120*(1 + delta)^5) ;
	c = (kz*count)/(1 + delta) - (kz^2*period*count*(-1 + count^2))/(6*(1 + delta)^3); 
	d = (kz*period^2*count*(-1 + count^2))/(12*(1 + delta)^3) + (kz^2*period^3*count*(-1 + count^4))/(120*(1 + delta)^5) ;
	diagonal = {a, b, c, d} ;
	matrix = MatrixExp[{{0, 1, 0, 0}, {-1, 0, 0, 0}, {0, 0, 0, 1}, {0, 0, -1, 0}} . DiagonalMatrix[diagonal]] ;
	{diagonal, matrix}   
] ;

(* --------- symplectic identity matrix --------- *)

ClearAll[identity] ;
identity[dimension_] := KroneckerProduct[IdentityMatrix[dimension], {{0, 1}, {-1, 0}}] ;

(* --------- symplectify --------- *)

ClearAll[symplectify] ;
symplectify::usage = "symplectify[matrix] -- symplectify given matrix (symplectic projection)" ;
symplectify[matrix_] := Block[
  {size, dimension, E, S, V, W},
  size = Length[matrix] ;
  dimension = 1/2*size ;
  E = IdentityMatrix[size] ;
  S = identity[dimension] ;
  V = S.(E - matrix).LinearSolve[E + matrix, E] ;
  W = 1/2*(V + Transpose[V]) ;
  LinearSolve[S + W, S - W]
] ;

(* --------- transport matrix --------- *)

ClearAll[matrix] ;
Options[matrix] = {
	Method -> {"FixedStep", Method -> {"StiffnessSwitching", Method -> {"ExplicitRungeKutta", Automatic}}},
	"Thin" -> True
} ;
matrix::usage = "matrix[object, energy, delta, {start, stop}, steps, angles, epsilon, solver] -- compute transport matrix " ;
matrix[
	object_,                  (* -- magnet object *)
	energy_,                  (* -- reference energy (GeV) *)
	delta_,                   (* -- energy deviation *)
	{start_, stop_},          (* -- integration start and stop positions relative to the object *)
	steps_,                   (* -- number of integration steps *)
	angles_,                  (* -- angle kicks to be added to canonical momenta on entrance and exit (cxi, czi, cxf, czf) *)
	epsilon_,                 (* -- epsilon initial condition *)
	solver_:radFldPtcTrj,     (* -- solver (radFldPtcTrj or RadFldPtcTrj) *)
	options:OptionsPattern[]  (* -- options *)
] := Block[
	{initial, positive, negative},
	positive = transport[object, energy, delta, {start, stop}, steps, angles, solver, options] /@ (+ epsilon*IdentityMatrix[4]) ;
	negative = transport[object, energy, delta, {start, stop}, steps, angles, solver, options] /@ (- epsilon*IdentityMatrix[4]) ;
	Transpose[(positive - negative)/(2*epsilon)]
] ;

(* --------- transport matrix (parameterization) --------- *)

ClearAll[parameterize] ;
Options[parameterize] = {
	Method -> {"FixedStep", Method -> {"StiffnessSwitching", Method -> {"ExplicitRungeKutta", Automatic}}},
	"Thin" -> True
} ;
parameterize::usage = "parameterize[object, energy, delta, {start, stop}, steps, angles, epsilon, solver] -- parameterize transport matrix" ;
parameterize[
	object_,                  (* -- magnet object *)
	energy_,                  (* -- reference energy (GeV) *)
	delta_,                   (* -- energy deviation *)
	{start_, stop_},          (* -- integration start and stop positions relative to the object *)
	steps_,                   (* -- number of integration steps *)
	angles_,                  (* -- angle kicks to be added to canonical momenta on entrance and exit (cxi, czi, cxf, czf) *)
	epsilon_,                 (* -- epsilon initial condition *)
	solver_:radFldPtcTrj,     (* -- solver (radFldPtcTrj or RadFldPtcTrj) *)
	options:OptionsPattern[]  (* -- options *)
] := Block[
	{transport, symplectic, positive, negative, derivative, identity, A, B},
	transport = matrix[object, energy, 0.0, {start, stop}, steps, angles, epsilon, solver, options] ;
	symplectic = symplectify[transport] ;
	positive = matrix[object, energy, +delta, {start, stop}, steps, angles, epsilon, solver, options] ;
	negative = matrix[object, energy, -delta, {start, stop}, steps, angles, epsilon, solver, options] ;
	derivative = (positive - negative)/2/delta ;
	identity = {{0, 1, 0, 0}, {-1, 0, 0, 0}, {0, 0, 0, 1}, {0, 0, -1, 0}} ;
	A = Re[- identity . MatrixLog[symplectic]] ;
	B = - identity . MatrixExp[- identity . A] . derivative ;
	{transport, symplectic, A, B, 1/2 (B + Transpose[B])}
] ;

(* --------- kick map table generation & export (one period) --------- *)

ClearAll[ndxp];
Options[potential] = Join[{"SamplesPerHarmonic" -> 32, "Samples" -> Automatic, NIntegrate -> False}, Options[NIntegrate]] ;
ndxp::usage = "ndxp[object, {x, z}, periods, harmonics, shift, delta, options] -- compute the normalized one-(super)period horizontal slope kick from the Elleaume potential using a central finite difference. The transverse point {x, z}, periods, shift, and finite-difference delta are in millimeters; the returned value is not energy scaled." ;
ndxp[object_, {x_, z_}, periods_, harmonics_, shift_, delta_, options:OptionsPattern[]] := Module[{pa, pb},
    pa = potential[object, {x - delta/2, z}, periods, harmonics, shift, Sequence @@ FilterRules[{options}, Options[potential]]] ;
    pb = potential[object, {x + delta/2, z}, periods, harmonics, shift, Sequence @@ FilterRules[{options}, Options[potential]]] ;
    (pb - pa)/delta*0.5
] ;

ClearAll[ndzp];
Options[potential] = Join[{"SamplesPerHarmonic" -> 32, "Samples" -> Automatic, NIntegrate -> False}, Options[NIntegrate]] ;
ndzp::usage = "ndzp[object, {x, z}, periods, harmonics, shift, delta, options] -- compute the normalized one-(super)period vertical slope kick from the Elleaume potential using a central finite difference. The transverse point {x, z}, periods, shift, and finite-difference delta are in millimeters; the returned value is not energy scaled." ;
ndzp[object_, {x_, z_}, periods_, harmonics_, shift_, delta_, options : OptionsPattern[]] := Module[{pa, pb},
    pa = potential[object, {x, z - delta/2}, periods, harmonics, shift, Sequence @@ FilterRules[{options}, Options[potential]]] ;
    pb = potential[object, {x, z + delta/2}, periods, harmonics, shift, Sequence @@ FilterRules[{options}, Options[potential]]] ;
    (pb - pa)/delta*0.5
] ;

ClearAll[round] ;
round::usage = "round[x, digits] -- round x to the requested number of decimal digits. Use digits -> None to return N[x] without decimal rounding." ;
round[x_, digits_Integer?NonNegative] := N[Round[x*10^digits]/10^digits] ;
round[x_, None] := N[x] ;

ClearAll[table] ;
Options[table] = Join[{
	"XRange" -> {-20.0, 20.0},      (* -- (mm) *)
	"ZRange" -> {-20.0, 20.0},      (* -- (mm) *)
	"XStep" -> 0.5,                 (* -- (mm) *)
	"ZStep" -> 0.5,                 (* -- (mm) *)
	"Delta" -> 0.1,                 (* -- finite-difference step (mm) *)
	"Digits" -> None,               (* -- number of digits to keep *)
	"File" -> "table",              (* -- output file name *)
	"Point" -> {0.0, 0.0},          (* -- reference point transverse position (mm) *)
	"Energy" -> None,               (* -- None or energy value (GeV) *)
	"KickScales" -> 10.0^-6,        (* -- kick scale factors *)
	"KickSigns" -> {-1, -1},        (* -- kick signs (use {-1, -1} for AT and {1, 1} in WM) *)
	"Period" -> Automatic           (* -- period length (mm) *)
}, Options[potential]] ;
table::usage = "table[object, periods, harmonics, shift, options] -- generate a one-(super)period kick-map table on the configured transverse grid, export it as a MATLAB .mat file, and return the labeled data. XRange, ZRange, XStep, ZStep, Delta, Point, and Period are in millimeters; exported xtable, ytable, and Len are in meters. If Energy is None, kicks are normalized; otherwise kicks are scaled by (0.299792458/Energy)^2 and KickScales." ;
table[                        (* -- generate and export AT kick map table *)
    object_,                  (* -- radia object *)
	periods_,                 (* -- horizontal and vertical periods (mm), {ph, pv} = {n*p, p} or {p, n*p} and n*p -- super-period *)
	harmonics_,               (* -- list of harmonics *)
	shift_,                   (* -- longitudinal shift/position (mm) *)
	options: OptionsPattern[] (* -- option(s) *) 
] := Block[
	{xmin, xmax, dx, nptx, xGrid, xtable, zmin, zmax, dz, nptz, zGrid, ztable, delta, digits, point, energy, factor, scale, sx, sz, length, xp, zp, result, file},
	{xmin, xmax} = OptionValue["XRange"] ;
	{zmin, zmax} = OptionValue["ZRange"] ;
	dx = OptionValue["XStep"] ;
	dz = OptionValue["ZStep"] ;
	nptx = Round[(xmax - xmin)/dx] + 1 ;
	nptz = Round[(zmax - zmin)/dz] + 1 ;
	xGrid = N[Subdivide[xmin, xmax, nptx - 1]] ;
	zGrid = N[Subdivide[zmax, zmin, nptz - 1]] ;
	xtable = xGrid/1000.0 ;
	ztable = zGrid/1000.0 ;
	delta = OptionValue["Delta"] ;
	digits = OptionValue["Digits"] ;
	point = OptionValue["Point"] ;
	energy = OptionValue["Energy"] ;
	factor = If[energy === None, 1.0, (0.299792458/energy)^2] ;
	scale = factor*OptionValue["KickScales"];
	{sx, sz} = OptionValue["KickSigns"];
	length = OptionValue["Period"] ;
	length = If[length === Automatic, Max[periods], length] ;
	xp = Table[ndxp[object, point + {x, z}, periods, harmonics, shift, delta, Sequence @@ FilterRules[{options}, Options[ndxp]]], {z, zGrid}, {x, xGrid}] ;
	zp = Table[ndzp[object, point + {x, z}, periods, harmonics, shift, delta, Sequence @@ FilterRules[{options}, Options[ndzp]]], {z, zGrid}, {x, xGrid}] ;
	xp = Transpose[round[sx*scale*xp, digits]] ;
	zp = Transpose[round[sz*scale*zp, digits]] ;
	result = {
		"xkick" -> 0*xp,
		"ykick" -> 0*zp,
		"xkick1" -> xp,
		"ykick1" -> zp,
		"xtable" -> {xtable},
		"ytable" -> {ztable},
		"Len" -> {{N[length/1000.0]}}
	} ;
	file = OptionValue["File"] <> ".mat" ;
	Export[file, result, "LabeledData"] ;
	result
] ;

(* --------- kick map table tracking --------- *)

ClearAll[interpolate];
Options[interpolate] = {"InterpolationOrder" -> 3} ;
interpolate::usage = "interpolate[map, options] -- convert a labeled kick-map table to an association containing ListInterpolation functions for horizontal and vertical kicks, sorted transverse ranges, and map length. The input map is the labeled data returned by table or imported from the exported .mat file." ;
interpolate[map_, options : OptionsPattern[]] := Block[
	{table, xGrid, zGrid, xKick, zKick, xOrder, zOrder, order},
	table = Association[map] ;
	xGrid = N[Flatten[table["xtable"]]] ;
	zGrid = N[Flatten[table["ytable"]]] ;
	xKick = N[table["xkick1"]] ;
	zKick = N[table["ykick1"]] ;
	xOrder = Ordering[xGrid] ;
	zOrder = Ordering[zGrid] ;
	xGrid = xGrid[[xOrder]] ;
	zGrid = zGrid[[zOrder]] ;
	xKick = xKick[[xOrder, zOrder]] ;
	zKick = zKick[[xOrder, zOrder]] ;
	order = OptionValue["InterpolationOrder"] ;
	Association[
		"XKick" -> ListInterpolation[xKick, {xGrid, zGrid}, InterpolationOrder -> order],
		"YKick" -> ListInterpolation[zKick, {xGrid, zGrid}, InterpolationOrder -> order],
		"XRange" -> MinMax[xGrid],
		"ZRange" -> MinMax[zGrid],
		"Length" -> First[Flatten[table["Len"]]]
	]
] ;

(* --------- kick map dkd based tracking --------- *)

ClearAll[km] ;
Options[km] = {
	"Period" -> Automatic,          (* -- period length (m) *)
	"Energy" -> None,               (* -- None if map is already energy scaled or energy in GeV *)
	"KickScales" -> 1.0,            (* -- extra kick scale *)
	"KickSigns" -> {-1, -1},        (* -- optional sign flip during tracking *)
	"Delta" -> True,                (* -- use energy deviation *)
	"InterpolationOrder" -> 3       (* -- interpolation order *)
} ;
km::usage = "km[map, delta, count, factors, options][{qx, px, qz, pz}] -- track canonical coordinates through count drift-kick-drift periods using a kick-map table. Coordinates qx and qz are in meters; px and pz are canonical momenta. Period is in meters, or Automatic to use Len from the map. If Energy is None, the map is assumed already scaled; otherwise kicks are scaled by (0.299792458/(Energy*(1 + delta)))^2." ;
km[                          (* -- drift-kick-drift tracking using kick map *)
    map_,                    (* -- kick map table *)
    delta_,                  (* -- energy delta *)
    count_,                  (* -- total number of kicks/periods *)
    factors_:  {1.0, 1.0},   (* -- extra kick multiplication factors *)
    options:  OptionsPattern[]
][state_] := Block[
    {table, kx, kz, length, DL, FX, FZ, SX, SZ, QX, PX, QZ, PZ, X, XP, Z, ZP, energy, scale},
    table = interpolate[map, Sequence @@ FilterRules[{options}, Options[interpolate]]] ;
    kx = table["XKick"];
    kz = table["YKick"];
    length = OptionValue["Period"] ;
    length = If[length === Automatic, table["Length"], length] ;
    DL = length/2.0;
    {FX, FZ} = factors ;
    {SX, SZ} = OptionValue["KickSigns"] ;
    energy = OptionValue["Energy"] ;
    scale = OptionValue["KickScales"]*If[energy === None, If[OptionValue["Delta"], 1.0/(1.0 + delta)^2, 1.0], (0.299792458/(energy*(1.0 + delta)))^2] ;
    {QX, PX, QZ, PZ} = state;
    {QX, PX, QZ, PZ} = {QX - count*DL*PX/Sqrt[(1.0 + delta)^2 - PX^2 - PZ^2], PX, QZ - count*DL*PZ/Sqrt[(1.0 + delta)^2 - PX^2 - PZ^2], PZ} ;
    Do[
        {QX, PX, QZ, PZ} = {QX + DL*PX/Sqrt[(1.0 + delta)^2 - PX^2 - PZ^2], PX, QZ + DL*PZ/Sqrt[(1.0 + delta)^2 - PX^2 - PZ^2], PZ} ;
        {X, XP, Z, ZP} = {QX, PX/Sqrt[(1.0 + delta)^2 - PX^2 - PZ^2], QZ, PZ/Sqrt[(1.0 + delta)^2 - PX^2 - PZ^2]} ;
        {X, XP, Z, ZP} = {X, XP - FX*SX*scale*kx[X, Z], Z, ZP - FZ*SZ*scale*kz[X, Z]} ;
        {QX, PX, QZ, PZ} = {X, (1.0 + delta)*XP/Sqrt[1.0 + XP^2 + ZP^2], Z, (1.0 + delta)*ZP/Sqrt[1.0 + XP^2 + ZP^2]} ;
        {QX, PX, QZ, PZ} = {QX + DL*PX/Sqrt[(1.0 + delta)^2 - PX^2 - PZ^2], PX, QZ + DL*PZ/Sqrt[(1.0 + delta)^2 - PX^2 - PZ^2], PZ},
        count
    ] ;
    {QX, PX, QZ, PZ} = {QX - count*DL*PX/Sqrt[(1.0 + delta)^2 - PX^2 - PZ^2], PX, QZ - count*DL*PZ/Sqrt[(1.0 + delta)^2 - PX^2 - PZ^2], PZ} ;
    {QX, PX, QZ, PZ}
] ;