
fvScalarMatrix cKEqn 
(
    fvm::ddt(cK)
    + fvm::div( phi , cK, "div(phi,c)") // advective flux
==
    fvm::laplacian ( D , cK  )
    + fvm::div((D*e / (k_B * T ) )*fvc::snGrad(Psi_E) * mesh.magSf(), cK, "div(eMigFlux,cK)") 
);
cKEqn.relax();
cKEqn.solve(mesh.solutionDict().solver("c"));

fvScalarMatrix cClEqn 
(
    fvm::ddt(cCl)
    + fvm::div( phi , cCl, "div(phi,c)") // advective flux
==
    fvm::laplacian ( D , cCl  ) 
    + fvm::div((-1* D*e / (k_B * T ) )*fvc::snGrad(Psi_E) * mesh.magSf(), cCl, "div(eMigFlux,cCl)")
);
cClEqn.relax();
cClEqn.solve(mesh.solutionDict().solver("c"));
 
