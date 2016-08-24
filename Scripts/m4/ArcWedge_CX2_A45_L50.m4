/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.3.0                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant/polyMesh";
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//Parametrized Orifice Geometry

//Run using:
//m4 -P blockMeshDict.m4 > blockMeshDict

//m4 definitions:
m4_changecom(//)m4_changequote([,])
m4_define(calc, [m4_esyscmd(perl -e 'use Math::Trig; printf ($1)')])
m4_define(VCOUNT, 0)
m4_define(vlabel, [[// ]Vertex $1 = VCOUNT m4_define($1, VCOUNT)m4_define([VCOUNT],
m4_incr(VCOUNT))])
//Mathematical constants:
m4_define(pi, 3.14159)
m4_define(rad, [calc($1*pi/180.0)])

//m4 spec:
m4_define(InletHeight, 3)
m4_define(OutletHeight, 3)
m4_define(ThroatHeight, 1)

m4_define(ChamferX, 0.2)
m4_define(ChamferAngle, rad(45.0))

m4_define(ThroatLength, 0.5)
m4_define(BackLength, 0.9)
m4_define(BackAngle, rad(60.0))
m4_define(OutletLength, 73.8)
m4_define(InletLength, calc(80-OutletLength))

m4_define(Boundary, 0.05)

m4_define(ChamferYCoor, calc(ThroatHeight+(ChamferX*tan(ChamferAngle))))
m4_define(ChamferXCoor, calc(-ChamferX))

m4_define(BackY, calc(ThroatHeight+(BackLength*tan(BackAngle))))
m4_define(BackX, calc(ThroatLength+BackLength))


//Arc E Inlet
m4_define(xE, calc(-InletLength))

//Arc F
m4_define(xFu, calc(-ChamferX))
m4_define(xFd, calc(-7.5*ChamferX))

//Arc G ChamferStart
m4_define(xGu, calc(-ChamferX))
m4_define(xGd, calc(-3.9*ChamferX))

//Arc H ThroatStart
m4_define(xHu, 0)
m4_define(xHd, calc(-1.25*ChamferX))

//Arc A ThroatEnd
m4_define(xAu, calc(ThroatLength))
m4_define(xAd, calc(1.3*ThroatLength))

//Arc B BackEnd
m4_define(xBu, calc(BackX))
m4_define(xBd, calc(1.53*BackX))

//Arc C
m4_define(xCu, calc(BackX))
m4_define(xCd, calc(1.85*BackX))

//Arc D Outlet
m4_define(xD, calc(OutletLength))

//plane z:
m4_define(z, calc(InletHeight*tan(rad(2.5))))
m4_define(zHA, calc(ThroatHeight*tan(rad(2.5))))
m4_define(zG, calc(ChamferYCoor*tan(rad(2.5))))
m4_define(zB, calc(BackY*tan(rad(2.5))))

m4_define(ArcL, 0.1)

//Angles Arc
m4_define(Ang, rad(45.0))
m4_define(AngG, calc((rad(270)-ChamferAngle)/1.8))
m4_define(AngH, calc((rad(180)+ChamferAngle)/1.5))
m4_define(AngA, calc((rad(180)+BackAngle)/1.5))
m4_define(AngB, calc((rad(270)-BackAngle)/1.5))

//Arc Points X
m4_define(ArcFx, calc(ChamferXCoor-(ArcL*sin(Ang))))
m4_define(ArcGx, calc(ChamferXCoor-(ArcL*sin(AngG))))
m4_define(ArcHx, calc(-ArcL*sin(AngH)))
m4_define(ArcAx, calc(ThroatLength+(ArcL*sin(AngA))))
m4_define(ArcBx, calc(BackX+(ArcL*sin(AngB))))
m4_define(ArcCx, calc(BackX+(ArcL*sin(Ang))))

//Arc Points Y
m4_define(ArcFy, calc(InletHeight-(ArcL*cos(Ang))))
m4_define(ArcGy, calc(ChamferYCoor+(ArcL*cos(AngG))))
m4_define(ArcHy, calc(ThroatHeight+(ArcL*cos(AngH))))
m4_define(ArcAy, calc(ThroatHeight+(ArcL*cos(AngA))))
m4_define(ArcBy, calc(BackY+(ArcL*cos(AngB))))
m4_define(ArcCy, calc(OutletHeight-(ArcL*cos(Ang))))

//Arc Points Z
m4_define(zArcF, calc(ArcFy*tan(rad(2.5))))
m4_define(zArcG, calc(ArcGy*tan(rad(2.5))))
m4_define(zArcH, calc(ArcHy*tan(rad(2.5))))
m4_define(zArcA, calc(ArcAy*tan(rad(2.5))))
m4_define(zArcB, calc(ArcBy*tan(rad(2.5))))
m4_define(zArcC, calc(ArcCy*tan(rad(2.5))))

/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.3.0                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant/polyMesh";
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

convertToMeters 0.001;

vertices
(
	// A
	(xAd 0 0) vlabel(A0)
	(xAu ThroatHeight -zHA) vlabel(A2)
	(xAu ThroatHeight zHA) vlabel(A3)

	// B
	(xBd 0 0) vlabel(B0)
	(xBu BackY -zB) vlabel(B2)
	(xBu BackY zB) vlabel(B3)

	// C
	(xCd 0 0) vlabel(C0)
	(xCu OutletHeight -z) vlabel(C2)
	(xCu OutletHeight z) vlabel(C3)

	// D
	(xD 0 0) vlabel(D0)
	(xD OutletHeight -z) vlabel(D2)
	(xD OutletHeight z) vlabel(D3)

	// E
	(xE 0 0) vlabel(E0)
	(xE InletHeight -z) vlabel(E2)
	(xE InletHeight z) vlabel(E3)

	// F
	(xFd 0 0) vlabel(F0)
	(xFu InletHeight -z) vlabel(F2)
	(xFu InletHeight z) vlabel(F3)

	// G
	(xGd 0 0) vlabel(G0)
	(xGu ChamferYCoor -zG) vlabel(G2)
	(xGu ChamferYCoor zG) vlabel(G3)

	// H
	(xHd 0 0) vlabel(H0)
	(xHu ThroatHeight -zHA) vlabel(H2)
	(xHu ThroatHeight zHA) vlabel(H3)
);

edges
(
	arc F0 F3 (ArcFx ArcFy zArcF)
	arc F0 F2 (ArcFx ArcFy -zArcF)
	arc G0 G3 (ArcGx ArcGy zArcG)
	arc G0 G2 (ArcGx ArcGy -zArcG)
	arc H0 H3 (ArcHx ArcHy zArcH)
	arc H0 H2 (ArcHx ArcHy -zArcH)
	arc A0 A3 (ArcAx ArcAy zArcA)
	arc A0 A2 (ArcAx ArcAy -zArcA)
	arc B0 B3 (ArcBx ArcBy zArcB)
	arc B0 B2 (ArcBx ArcBy -zArcB)
	arc C0 C3 (ArcCx ArcCy zArcC)
	arc C0 C2 (ArcCx ArcCy -zArcC)
);

blocks
(
	hex (A0 B0 B2 A2 A0 B0 B3 A3) (100 120 1)
	simpleGrading
	(
		1
		(
			(0.95 0.8 1)
			(0.05 0.2 0.25)

		)
		1
	)

	hex (B0 C0 C2 B2 B0 C0 C3 B3) (30 120 1)
	simpleGrading
	(
		1
		(
			(0.95 0.8 1)
			(0.05 0.2 0.25)

		)
		1
	)
	hex (C0 D0 D2 C2 C0 D0 D3 C3) (1000 120 1)
	simpleGrading
	(
		7
		(
			(0.95 0.8 1)
			(0.05 0.2 0.25)

		)
		1
	)
	hex (E0 F0 F2 E2 E0 F0 F3 E3) (160 120 1)
	simpleGrading
	(
		0.3
		(
			(0.95 0.8 1)
			(0.05 0.2 0.25)

		)
		1
	)
	hex (F0 G0 G2 F2 F0 G0 G3 F3) (100 120 1)
	simpleGrading
	(
		1
		(
			(0.95 0.8 1)
			(0.05 0.2 0.25)

		)
		1
	)
	hex (G0 H0 H2 G2 G0 H0 H3 G3) (40 120 1)
	simpleGrading
	(
		1
		(
			(0.95 0.8 1)
			(0.05 0.2 0.25)

		)
		1
	)
	hex (H0 A0 A2 H2 H0 A0 A3 H3) (120 120 1)
	simpleGrading
	(
		1
		(
			(0.95 0.8 1)
			(0.05 0.2 0.25)

		)
		1
	)
);

boundary
(
	inlet
	{
            type inlet;

            faces
            (
                (E3 E2 E0 E0)
            );
        }
        
	outlet
	{
            type outlet;

            faces
            (
                (D0 D2 D3 D0)
            );
        }
        
	fixedwall
	{
            type wall;
            
            faces
            (
                (E3 F3 F2 E2)
                (F3 G3 G2 F2)
                (G3 H3 H2 G2)
                (H3 A3 A2 H2)
                (A3 B3 B2 A2)
                (B3 C3 C2 B2)
                (C3 D3 D2 C2)
            );
        }
	wedge1
	{
            type wedge;
            
            faces
            (
                (E0 E2 F2 F0)
                (F0 F2 G2 G0)
                (G0 G2 H2 H0)
                (H0 H2 A2 A0)
                (A0 A2 B2 B0)
                (B0 B2 C2 C0)
                (C0 C2 D2 D0)
            );
            
        }
        
	wedge2
	{
            type wedge;
            
            faces
            (
                (F0 F3 E3 E0)
                (G0 G3 F3 F0)
                (H0 H3 G3 G0)
                (A0 A3 H3 H0)
                (B0 B3 A3 A0)
                (C0 C3 B3 B0)
                (D0 D3 C3 C0)
            );
        }
	axis
	{
            type empty;
            
            faces
            (
                (E0 E0 F0 F0)
                (F0 F0 G0 G0)
                (G0 G0 H0 H0)
                (H0 H0 A0 A0)
                (A0 A0 B0 B0)
                (B0 B0 C0 C0)
                (C0 C0 D0 D0)
            );
        }   
);

mergePatchPairs
(
);
// ************************************************************************* //
