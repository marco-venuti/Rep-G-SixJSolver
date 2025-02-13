# Function to evaluate a cyclotomic expression
EvaluateCyclotomic := function(expr)
  local coeffs, n, realPart, k, angle, pi;
  #Define pi
  pi := 3.141592653589793;
  
  # Determine the order of the cyclotomic field (n)
  n := Conductor(expr);

  # Get the coefficients of the cyclotomic expression
  coeffs := CoeffsCyc(expr, n);

  # Initialize real and imaginary parts
  realPart := 0;

  # Loop through each coefficient to calculate the real and imaginary parts
  for k in [1..Length(coeffs)] do
    if coeffs[k] <> 0 then
      # Calculate the angle (2 * pi * (k-1) / n)
      angle := 2 * pi * (k - 1) / n;

      # Add the contribution to the real part
      realPart := realPart + coeffs[k] * Cos(angle);
    fi;
  od;

  return realPart;
end;


# Function to check if a matrix has a column with all zero entries
HasZeroColumn := function(mat)
    local n, m, col, i, isZeroColumn;

    # Get the number of rows and columns of the matrix
    n := Length(mat);
    m := Length(mat[1]);

    # Iterate through each column
    for col in [1..m] do
        # Assume the column is all zeros
        isZeroColumn := true;

        # Check each entry in the column
        for i in [1..n] do
            if Rat(EvaluateCyclotomic(mat[i][col])) <> 0 then
                isZeroColumn := false;
                break; # Exit loop early if a nonzero entry is found
            fi;
        od;

        # If an all-zero column is found, return true
        if isZeroColumn then
            return true;
        fi;
    od;

    # If no all-zero column is found, return false
    return false;
end;



NormalizeColumns := function(mat)
    
    local nrows, ncols, colNorm, i, j, resultMat, entry, entry_norm_squared, entry_string;
    
    # Ensure the input is a matrix
    if not IsMatrix(mat) then
        Error("Input must be a matrix.");
    fi;

    # Get the dimensions of the matrix
    nrows := DimensionsMat(mat)[1];
    ncols := DimensionsMat(mat)[2];

    # Initialize an empty matrix to store the result
    resultMat := NullMat(nrows, ncols);

    # Loop through each column
    for j in [1..ncols] do
        # Compute the norm of the j-th column
        colNorm := 0;
        for i in [1..nrows] do
            entry := mat[i][j];
            #Print("Here is the entry: ", entry, "\n");
            if not ImaginaryPart(entry) = 0 then
                #entry_norm_squared := AbsoluteValue(Norm(mat[i][j])); #this is the norm-squared of a non-real complex number
                entry_norm_squared := EvaluateCyclotomic(entry*ComplexConjugate(entry)); #this is the norm-squared of a non-real complex number
                
            else #(when the entry is real)
                #in GAP, Norm() method behaves differently for objects represented as linear combinations of roots of unity, even if real
                #Examples of this discrepancy include Norm(2)=2 but Norm(2i)=4; Similarly, Norm(E(8)+E(8)^7)=4, even though E(8)+E(8)^7 = Sqrt(2)
                #We must check if the entry contains roots of unity, and operate differently in each case
                
                entry_string := String(entry); #Solution: convert the entry to a string to check if E(n) is present or not
                if 'E' in entry_string then  
                    entry_norm_squared := EvaluateCyclotomic(entry)^2;
                else
                    entry_norm_squared:= entry^2; #the norm-squared for real numbers 
                fi;
            fi;
            #Print("Here is the entry norm: ", entry_norm_squared, "\n");
            colNorm := colNorm + entry_norm_squared;
        od;

        colNorm := (colNorm)^0.5;
        #Print("Here is the norm of the column: ", colNorm, "\n");
        
        
        for i in [1..nrows] do
            resultMat[i][j] := mat[i][j] / Rat(colNorm); #because division of complex numbers by floats is not supported, we use a rational approximation of the norm
        od;
    od;

    return resultMat;
end;
