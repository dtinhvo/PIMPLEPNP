// Time derivative matrix
fvVectorMatrix ddtUEqn(fvm::ddt(U));

volVectorField MaxwellForce
(
    - F*(cK - cCl ) * (fvc::grad(Psi_E) - E_field_ )()  
);
const_cast <dimensionSet&> (MaxwellForce.dimensions()) /=  dimMass/dimVolume ; 

// Convection-diffusion matrix
fvVectorMatrix HUEqn
(
    fvm::div(phi, U)
  + turbulence->divDevReff()
    ==
    MaxwellForce /1000 //  /rho : p is /rho, also divdevreff
);

if (pimple.momentumPredictor())
{
    // Solve momentum predictor
    solve(relax(ddtUEqn + HUEqn) == -fvc::grad(p));
}
else
{
    // Explicit update
    U = (ddtUEqn.H() + HUEqn.H() - fvc::grad(p))/(HUEqn.A() + ddtUEqn.A());
    U.correctBoundaryConditions();
}
