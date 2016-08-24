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

m4_define(Chamfer, 0.2)

m4_define(ThroatLength, 1)
m4_define(BackLength, 0.9)
m4_define(BackAngle, rad(60.0))
m4_define(OutletLength, 73.8)
m4_define(InletLength, calc(80-OutletLength))

m4_define(Boundary, 0.05)

m4_define(ChamferY, calc(ThroatHeight+Chamfer))
m4_define(ChamferX, calc(-Chamfer))

m4_define(BackY, calc(ThroatHeight+(BackLength*tan(BackAngle))))
m4_define(BackX, calc(ThroatLength+BackLength))


//Arc E Inlet
m4_define(xE, calc(-InletLength))

//Arc F
m4_define(xFu, calc(-Chamfer))
m4_define(xFd, calc(-7.5*Chamfer))

//Arc G ChamferStart
m4_define(xGu, calc(-Chamfer))
m4_define(xGd, calc(-2.75*Chamfer))

//Arc H ThroatStart
m4_define(xHu, 0)
m4_define(xHd, calc(-1.25*Chamfer))

//Arc A ThroatEnd
m4_define(xAu, calc(ThroatLength))
m4_define(xAd, calc(1.25*ThroatLength))

//Arc B BackEnd
m4_define(xBu, calc(BackX))
m4_define(xBd, calc(1.25*BackX))

//Arc C
m4_define(xCu, calc(BackX))
m4_define(xCd, calc(1.43*BackX))

//Arc D Outlet
m4_define(xD, calc(OutletLength))

//plane z:
m4_define(z, 0.2)


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
	(xAd 0 -z) vlabel(A0)
	(xAd 0 z) vlabel(A1)
	(xAu ThroatHeight -z) vlabel(A2)
	(xAu ThroatHeight z) vlabel(A3)

	// B
	(xBd 0 -z) vlabel(B0)
	(xBd 0 z) vlabel(B1)
	(xBu BackY -z) vlabel(B2)
	(xBu BackY z) vlabel(B3)

	// C
	(xCd 0 -z) vlabel(C0)
	(xCd 0 z) vlabel(C1)
	(xCu OutletHeight -z) vlabel(C2)
	(xCu OutletHeight z) vlabel(C3)

	// D
	(xD 0 -z) vlabel(D0)
	(xD 0 z) vlabel(D1)
	(xD OutletHeight -z) vlabel(D2)
	(xD OutletHeight z) vlabel(D3)

	// E
	(xE 0 -z) vlabel(E0)
	(xE 0 z) vlabel(E1)
	(xE InletHeight -z) vlabel(E2)
	(xE InletHeight z) vlabel(E3)

	// F
	(xFd 0 -z) vlabel(F0)
	(xFd 0 z) vlabel(F1)
	(xFu InletHeight -z) vlabel(F2)
	(xFu InletHeight z) vlabel(F3)

	// G
	(xGd 0 -z) vlabel(G0)
	(xGd 0 z) vlabel(G1)
	(xGu ChamferY -z) vlabel(G2)
	(xGu ChamferY z) vlabel(G3)

	// H
	(xHd 0 -z) vlabel(H0)
	(xHd 0 z) vlabel(H1)
	(xHu ThroatHeight -z) vlabel(H2)
	(xHu ThroatHeight z) vlabel(H3)
);

edges
(
	arc F1 F3 (-1.4 0.6 z)
	arc F0 F2 (-1.4 0.6 -z)
	arc G1 G3 (-0.27 1.1 z)
	arc G0 G2 (-0.27 1.1 -z)
	arc H1 H3 (-0.15 0.6 z)
	arc H0 H2 (-0.15 0.6 -z)
	arc A1 A3 (1.15 0.6 z)
	arc A0 A2 (1.15 0.6 -z)
	arc B1 B3 (2.3 1 z)
	arc B0 B2 (2.3 1 -z)
	arc C1 C3 (2.6 1 z)
	arc C0 C2 (2.6 1 -z)
);

blocks
(
	hex (A0 B0 B2 A2 A1 B1 B3 A3) (100 100 10)
	simpleGrading
	(
		1
		(
			(0.95 0.8 1)
			(0.05 0.2 0.25)

		)
		1
	)

	hex (B0 C0 C2 B2 B1 C1 C3 B3) (30 100 10)
	simpleGrading
	(
		1
		(
			(0.95 0.8 1)
			(0.05 0.2 0.25)

		)
		1
	)
	hex (C0 D0 D2 C2 C1 D1 D3 C3) (900 100 10)
	simpleGrading
	(
		7
		(
			(0.95 0.8 1)
			(0.05 0.2 0.25)

		)
		1
	)
	hex (E0 F0 F2 E2 E1 F1 F3 E3) (100 100 10)
	simpleGrading
	(
		0.3
		(
			(0.95 0.8 1)
			(0.05 0.2 0.25)

		)
		1
	)
	hex (F0 G0 G2 F2 F1 G1 G3 F3) (100 100 10)
	simpleGrading
	(
		1
		(
			(0.95 0.8 1)
			(0.05 0.2 0.25)

		)
		1
	)
	hex (G0 H0 H2 G2 G1 H1 H3 G3) (40 100 10)
	simpleGrading
	(
		1
		(
			(0.95 0.8 1)
			(0.05 0.2 0.25)

		)
		1
	)
	hex (H0 A0 A2 H2 H1 A1 A3 H3) (70 100 10)
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
                (E1 E3 E2 E0)
            );
        }
        
	outlet
	{
            type outlet;

            faces
            (
                (D0 D2 D3 D1)
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
            type cyclic;
            
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
            
            neighbourPatch  wedge2;
        }
        
	wedge2
	{
            type cyclic;
            
            faces
            (
                (F1 F3 E3 E1)
                (G1 G3 F3 F1)
                (H1 H3 G3 G1)
                (A1 A3 H3 H1)
                (B1 B3 A3 A1)
                (C1 C3 B3 B1)
                (D1 D3 C3 C1)
            );
            neighbourPatch  wedge1;
        }
	axis
	{
            type symmetry;
            
            faces
            (
                (E1 E0 F0 F1)
                (F1 F0 G0 G1)
                (G1 G0 H0 H1)
                (H1 H0 A0 A1)
                (A1 A0 B0 B1)
                (B1 B0 C0 C1)
                (C1 C0 D0 D1)
            );
        }    
);

mergePatchPairs
(
);
// ************************************************************************* //
