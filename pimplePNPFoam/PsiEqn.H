fvScalarMatrix PsiEqn
    (
        fvm::laplacian( 77.9 * epsilon_0, Psi_E ) // ! epsilon_r of "fluid", relative
        ==
        - (F  *  (cK - cCl ) )
    );
PsiEqn.relax();
PsiEqn.solve();  

            



