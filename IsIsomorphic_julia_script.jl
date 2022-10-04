############################################################# 
# Stefano Marseglia, Utrecht University, s.marseglia@uu.nl
# http://www.staff.science.uu.nl/~marse004/
#############################################################
# Requires julia 1.7.0 or higher and Hecke
# http://www.thofma.com/Hecke.jl/dev/ 
# The keys step is the use if is_GLZ_conjugate, see
# http://www.thofma.com/Hecke.jl/dev/misc/conjugacy/
# which is based on https://arxiv.org/pdf/2202.03526.pdf
############################################################# 

    using Hecke;
    coeffs = eval(Meta.parse(open(f->read(f, String), ARGS[1])));
    n = Int(sqrt(length(coeffs[1])));
    candidates = [ matrix(ZZ,n,n,coeffs[1]) ];
    for i in 2:length(coeffs)
        push!(candidates,matrix(ZZ,n,n,coeffs[i]));
    end;
    # classes will contain the indices of the representatives, classes_mat the actual matrices
    classes = [1];
    classes_mat = [ matrix(ZZ,n,n,coeffs[1]) ];
    for i in 2:length(coeffs)
        I = matrix(ZZ,n,n,coeffs[i]);
        test = true;
        for j in 1:length(classes_mat)
            J = classes_mat[j];

            fl , T = is_GLZ_conjugate(I,J);
            if fl
                test = false;
                break; # break J
            end;
        end;
        if test
            push!(classes,i);
            push!(classes_mat,I);
        end;
    end;
    println(classes)
