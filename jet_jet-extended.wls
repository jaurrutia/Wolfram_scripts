

Jet[u_?NumericQ] :=  Blend[{{0, RGBColor[0, 0, 9/16]},
                            {1/9, Blue},
                            {23/63, Cyan},
                            {13/21, Yellow},
                            {47/63, Orange},
                            {55/63, Red},
                            {1, RGBColor[1/2, 0, 0]}
                         }, u] /; 0 <= u <= 1

JetExtendedColors = Table[If[i != 0 && i != 81,
                            {i/81, Jet[i/81]},
                            {i/81,  RGBColor[i, i, i]}
                    ], {i, 0, 81, 9}]

JetExtended[u_?NumericQ] := Blend[JetExtendedColors, u] /; 0 <= u <= 1
