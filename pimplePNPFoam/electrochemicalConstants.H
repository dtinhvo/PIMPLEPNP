const dimensionedScalar D // [C/mol] = [kg0⋅m2⋅s-1⋅K0⋅mol0⋅A0⋅cd0] 
        ( "binaryDiffusionCoefficient", 
            dimensionSet(0, 2, -1, 0, 0, 0, 0),
            1e-9);
const dimensionedScalar F // [C/mol] = [kg0⋅m0⋅s0⋅K0⋅mol-1⋅A1⋅cd0]
            ( "FaradayConstant", 
            dimensionSet(0, 0, 1, 0,-1, 1, 0),
            96485.33212331001); 
const dimensionedScalar e  // [C] = [kg0⋅m0⋅s1⋅K0⋅mol0⋅A1⋅cd0] 
            ("elementaryCharge",
            dimensionSet( 0, 0, 1, 0, 0, 1, 0 ),
            1.6021766208e-19);
const dimensionedScalar k_B  // [J/K] = [kg1⋅m2⋅s-2⋅K-1⋅mol0⋅A0⋅cd0] 
            ( "BoltzmannConstant", 
            dimensionSet( 1, 2, -2, -1, 0, 0, 0 ),
            1.38064852e-23);  
const dimensionedScalar T // [J/K] = [kg0⋅m0⋅s0⋅K1⋅mol0⋅A0⋅cd0]  
            ( "temperature", 
            dimensionSet( 0, 0, 0, 1, 0, 0, 0 ),
            300);  
const dimensionedScalar epsilon_0 = dimensionedScalar  
            ( "vacuumPermittivity", 
            dimensionSet(-1, -3, 4, 0, 0, 2, 0),
            8.8541878176e-12 ) ;  
