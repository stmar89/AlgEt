############################################################# 
# Stefano Marseglia, stefano.marseglia89@gmail.com
# https://stmar89.github.io/index.html
#############################################################
# Requires julia 1.7.0 or higher and Hecke
# http://www.thofma.com/Hecke.jl/dev/ 
# The keys step is the use if is_GLZ_conjugate, see
# http://www.thofma.com/Hecke.jl/dev/misc/conjugacy/
# which is based on https://arxiv.org/pdf/2202.03526.pdf
############################################################# 

    using Hecke;
    file=ARGS[1];
    coeffs = eval(Meta.parse(open(f->read(f, String), file)));
    # we erase immediately the temporary file
    run(`rm $file`);
    n = Int(sqrt(length(coeffs[1])));
    candidates = [ matrix(ZZ,n,n,coeffs[1]) ];
    for i in 2:length(coeffs)
        push!(candidates,matrix(ZZ,n,n,coeffs[i]));
    end;
    # candidates now contains the input matrices.
    # classes will contain the indices of the representatives of the conjugacy classes of matrices from candidates
    # classes_mat will contain the actual matrices represnting the class.
    # We initialize classes and classes_mat with the first element of candidates.
    # for each other matrix in candidates we check if it is conjugate to something we have already found in 
    # classes_mat; if so, we record the index in candidates of the rep and the matrix T realizing the conjugation;
    # if not, we record the index of the current candidate and the indentity matrix.
    classes = [1];
    classes_mat = [ candidates[1] ];
    II = identity_matrix(ZZ,n);
    isom_mat = [II];
    isom_ind = [1];
    for i in 2:length(coeffs)
        I = candidates[i];
        test = true;
        for j in 1:length(classes_mat)
            J = classes_mat[j];
            fl , T = is_GLZ_conjugate(I,J);
            if fl
                @assert isone(abs(det(T))) && T * I == J * T;
                test = false;
                push!(isom_mat,T);
                push!(isom_ind,classes[j]);
                break; # break J
            end;
        end;
        if test
        # we got a new class: we store it.
            push!(classes,i);
            push!(classes_mat,I);
            push!(isom_mat,II);
            push!(isom_ind,i);
        end;
    end;
    @assert length(isom_mat) == length(candidates);
    for i in 1:length(candidates)
        j = isom_ind[i];    
        T = isom_mat[i];
        @assert T * candidates[i] == candidates[j] * T;
    end;
    println([classes,isom_ind,isom_mat])
