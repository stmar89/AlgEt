This is the list of steps to do in order to submit/update the package in magma
It assumes that my own fork of the Magma-package repository has been cloned in
~/Magma-Maths-package
and that the AlgEt package is in
~/AlgEt

-) Make sure that licences are up-to-date. If needed update and run:
    magma ~/AlgEt/dev/magma_script_update_licence.txt
-) Run the following commands
    cp ~/AlgEt/AlgEtQ/*.m ~/Magma-Maths-package/package/Ring/AlgEtQ/
-) Make sure that 
    ~/Magma-Maths-package/package/Ring/AlgEtQ.spec
   contains all the .m files from
    ~/AlgEt/spec
-) Commit/push my fork of the magma repo.
-) Make a pull request to the original repo.
-) Send by email an explaination, with new doc file (if needed) and updated slow tests.

