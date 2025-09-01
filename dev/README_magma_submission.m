This is the list of steps to do in order to submit/update the package in magma
It assumes that my own fork of the Magma-package repository has been cloned in
~/Magma-Maths-package
and that the AlgEt package is in
~/AlgEt

-) Make sure that licences are up-to-date. If needed update and run:
    magma ~/AlgEt/dev/magma_script_update_licence.txt
-) Run the following commands
    mkdir -p ~/Magma-Maths-package/package/Algebra/AlgEt/AlgEtQ/ && 
    cp ~/AlgEt/AlgEtQ/*.m ~/Magma-Maths-package/package/Algebra/AlgEt/AlgEtQ/
    cat ~/AlgEt/spec >> ~/Magma-Maths-package/package/Algebra/Algebra.spec
