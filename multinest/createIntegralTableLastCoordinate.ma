#!/usr/local/bin/MathematicaScript -script

M       = ToExpression[$ScriptCommandLine[[2]]];
nPoints = ToExpression[$ScriptCommandLine[[3]]];


xMin             = -1.25 * Sqrt[M];
xMax             = 0.;
rMin             = 0.;
rMax             = 1.25 * Sqrt[M];
integrationStart = -3 * Sqrt[M];

xStep = (xMax - xMin)/(nPoints - 1)
rStep = (rMax - rMin)/(nPoints - 1)

Print[M]
Print[0]
Print[""]

F[t0_, t1_, r_, N_] := NIntegrate[(Exp[-x^2]*(r^2+x^2)^N)/(Sqrt[Pi]*HypergeometricU[-N, 1/2 - N, r^2]),
                                                {x, t0, t1},
                                                AccuracyGoal -> 6,
                                                Method -> {"LocalAdaptive", "SymbolicProcessing" -> Automatic}
                                     ];

convertToString[n_] := ToString[StringReplace[ToString[StringReplace[ToString[n, StandardForm], {"`*^" :> "e", "`" :> ""}]], {" " :> ""}]]

For[ i = 1, i <= nPoints, i++,

	r = rMin + (i - 1)*rStep;
	xLast = integrationStart;
	v = 0.;

	For[ j = 1, j <= (nPoints - 1), j++,
		x = xMin + (j - 1)*xStep;
		If[ Abs[v-0.5] < 10^-7,
			v = 0.5;
		,
			v = v + F[ xLast, x, r, M ];
		];
		If[ v > 0.5,
			v = 0.5;
		];
		Print[convertToString[r] <> " " <> convertToString[v] <> " " <> convertToString[x]];
		xLast = x;
	];

	Print[convertToString[r] <> " 0.5 " <> convertToString[xMax]];

]
