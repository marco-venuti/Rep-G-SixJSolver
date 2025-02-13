#----------------------------------------------------------------------------------
#Title: A 6j symbol solver for the representation category of finite groups
#Author: Dhruv Bhat, BSc. Mathematics from New York University Abu Dhabi
#Created under the mentorship of Dr. Sachin Valera, Postdoctoral Associate at Center for Topological and Quantum Systems
#----------------------------------------------------------------------------------
#Notes to Reader
#This program has been written in GAP, which can be installed from https://www.gap-system.org/install/
#In GAP, array indexing starts at 1, not 0, unlike programming languages like Python.


#---------LOAD PACKAGES FOR CHAR TABLES, IRREPS AND UNITARISATION COMMANDS--------
if not IsPackageLoaded("ctbllib") then
    LoadPackage("ctbllib");
fi;

if not IsPackageLoaded("RepnDecomp") then
    LoadPackage("RepnDecomp");
fi;

Read("matrix_methods.g"); #contains functions to implement parts of Isao Sakata's Algorithms

#---------DEFINE GROUP, SORT IRREDUCIBLE CHARACTERS BY DEGREE----------
#Group Name
#G := SymmetricGroup(5); #(DEFINE GROUP NAME HERE)
G:=SmallGroup(6,1);

#Sorted List of Irreducible Characters
charTable := CharacterTable(G);
irreducibleCharacters := Irr(charTable);
sortedChars := SortedCharacters(charTable, irreducibleCharacters, "degree"); #irred. characters in ascending order of degree. trivial character is first by default

#-----------COMPUTE FUSION RULES---------------------------------------
Read("G_fusion_finder.g");  # This imports fusionRules and irrepNames


#------------COMPUTE UNITARISED IRREDUCIBLE REPRESENTATIONS------------------
#Computing Irreps in the sorted order of irreducible characters
irreps := IrreducibleRepresentationsDixon(G, sortedChars);

#Unitarise the representations
irreps_unitary := [];
iteration_num := 1;
for rho in irreps do
    if IsUnitaryRepresentation(rho) then
        Add(irreps_unitary, rho);
    else
        unitarized_rho:= UnitaryRepresentation(rho);
        Add(irreps_unitary, unitarized_rho.unitary_rep);
    fi;
od;


#-------------CREATE A RECORD OF UNITARISED IRREDUCIBLE REPRESENTATIONS AND THEIR RESPECTIVE LABELS FOR EASY ACCESS---------
# Define a "GAP record" (~ python dictionary) to store matrices by label
irrep_dictionary := rec();
for i in [1..Length(sortedChars)] do
    label := irrepNames[i];  # Get the label from the list
    representation := irreps_unitary[i];  # Define a matrix or any other value
    irrep_dictionary.(label) := representation;  # Assign the matrix to the record field
od;


#-------------------COMPUTE CLEBSCH-GORDON COEFFICIENTS OF THE GROUP-----------------------------
Print("We now compute unitary intertwiners from tensor product representations to their block diagonal form!", "\n");

u_intertwinerLIST := rec(); #a GAP "record" which works like a python dictionary. The keys will be ordered pairs of irreps and the values are the unitary intertwiners that block diagonalize them into a direct sum of irreps

fusionIndex := 1;
for i in [2..repCount] do
    for j in [i..repCount] do
        #We work with irreps_unitary[i] and irreps_unitary[j] in this loop
        identical_irreps := false;
        if i = j then #check if we are fusing identical irreps
            identical_irreps := true;
        fi;

        #-------------CREATE THE DECOMPOSITION OF THIS PAIR OF IRREPS---------------
        currentDecomp := fusionRules[fusionIndex]; #a list of fusion outcomes for current pair of irreps
        #Create the block diagonal decomposition representation
        block_diagonal_decomp := irrep_dictionary.(currentDecomp[1]); #set the first irrep into the block diagonal form
        for k in [1..Length(currentDecomp)] do
            if k > 1 then
                block_diagonal_decomp := DirectSumOfRepresentations([block_diagonal_decomp, irrep_dictionary.(currentDecomp[k])]);
            fi;
        od;

        #-----APPLY ISAO SAKATA'S TRICK TO FIND THE UNITARY INTERTWINER THAT RELATES THE TENSOR PRODUCT REP WITH THE BLOCK DIAGONAL REP CONSTRUCTED ABOVE---
        tensor_rep_dimension := DegreeOfRepresentation(irreps_unitary[i])*DegreeOfRepresentation(irreps_unitary[j]); #the dimension of the tensor product representation in this iteration of the loop

        #(1)- Find unitary intertwiner that block diagonalizes irreps_unitary[i]*irreps_unitary[j]
        sakata1 := NullMat(tensor_rep_dimension, tensor_rep_dimension, Integers);

        #Since we use a random matrix in this algorithm, we may need to change the choice of random matrix if sakata1 happens to have a 0 column
        iteration_count := 0; #a safety measure to prevent an infinite loop
        while HasZeroColumn(sakata1) and i<10 do
            randMatrix1 := RandomUnimodularMat(tensor_rep_dimension:domain:=[-10..10]);  # Random invertible square matrix of dimension matching the tensor product representation populated by integers
            for g in Elements(G) do
                sakata1:= sakata1 + KroneckerProduct(irreps_unitary[i](g), irreps_unitary[j](g))*randMatrix1*TransposedMat(ComplexConjugate(block_diagonal_decomp(g)));
            od;
        od;
        sakata1 := NormalizeColumns(sakata1); #Normalize every column of the matrix

        #Print("The sakata matrix for ", irrepNames[i], " \\tensor ", irrepNames[j], " is: ", "\n", sakata1, "\n");
        #Print("sakata1 is invertible/cols orthogonal?: ", IsDiagonalMat(TransposedMat(sakata1)*sakata1), ", and Deteterminant = ", EvaluateCyclotomic(DeterminantMat(sakata1)), "\n");

        #Label and store this matrix
        entry1 := Concatenation("U_",irrepNames[i], "*", irrepNames[j]); #the label for the intertwiner that block diagonalizes irrepNames[i] \tensor irrepNames[j]
        u_intertwinerLIST.(entry1) := sakata1; #add this intertwiner to the record (dictionary) of intertwiners with the appropriate label

        #(2)- Find unitary intertwiner that block diagonalizes irreps_unitary[j]*irreps_unitary[i]
        if not identical_irreps then #only if the two irreps are distinct
            sakata2 := NullMat(tensor_rep_dimension, tensor_rep_dimension, Integers);
            while HasZeroColumn(sakata2) and i<10 do
                randMatrix2 := RandomUnimodularMat(tensor_rep_dimension:domain:=[-10..10]);  # Random invertible square matrix of dimension matching the tensor product representation populated by integers
                for g in Elements(G) do
                    sakata2:= sakata2 + KroneckerProduct(irreps_unitary[j](g), irreps_unitary[i](g))*randMatrix2*TransposedMat(ComplexConjugate(block_diagonal_decomp(g)));
                od;
            od;
            sakata2 := NormalizeColumns(sakata2); #Normalize every column of the matrix
            #Print("The sakata matrix for ", irrepNames[j], " \\tensor ", irrepNames[i], " is: ", "\n", sakata2, "\n");
            #Print("sakata2 is invertible/cols orthogonal?: ", IsDiagonalMat(TransposedMat(sakata2)*sakata2), ", and Deteterminant = ", EvaluateCyclotomic(DeterminantMat(sakata2)), "\n");


            #Label and Store this matrix
            entry2 := Concatenation("U_", irrepNames[j], "*", irrepNames[i]); #the dictionary key for the intertwiner that block diagonalizes irrepNames[j] \tensor irrepNames[i]
            u_intertwinerLIST.(entry2) := sakata2; #add this intertwiner to the record (dictionary) of intertwiners with the appropriate label

        fi;

        #increment the index to work with the next fusion
        fusionIndex := fusionIndex + 1;

    od;

od;


#-------------------------------EXTRACT THE TRIVALENT VERTICES IN REP(G)-----------------------------
#This part of the code slices portions of each of the "Unitary Sakata Intertwiners" and assigns them to trivalent vertices
trivalent_vertices := rec();
fusionIndex := 1;

for i in [2..repCount] do
    for j in [i..repCount] do
        #We will define all (splitting) trivalent vertices pertaining to irreps_unitary[i] and irreps_unitary[j]
        #-----
        currentDecomp := fusionRules[fusionIndex]; #Get a list of the fusion outcome of irreps_unitary[i] tensor irreps_unitary[j]
        #-----
        #Check if we are fusing identical irreps
        identical_irreps := false;
        if i = j then
            identical_irreps := true;
        fi;
        #----

        #(1)- GET ALL (splitting) TRIVALENT VERTICES for the ORDERED PAIR (irreps_unitary[i],irreps_unitary[j])
        current_unitary_intertwiner := u_intertwinerLIST.(Concatenation("U_", irrepNames[i],"*",irrepNames[j])); #access the unitary intertwiner that block diagonalizes irreps_unitary[i] \tensor irreps_unitary[j] into a direct sum of irreps

        #Print("The BIG unitary block diagonalizing intertwiner for ", irrepNames[i], "*", irrepNames[j], " is: ", Float(current_unitary_intertwiner), "\n");
        u_size := Length(current_unitary_intertwiner); #this is the dimension of the current_unitary_intertwiner matrix

        #Create a variable to record the number of columns of current_unitary_intertwiner which have been extracted
        u_col_index := 1; #each time a trivalent vertex is extracted, this variable will be incremented by the dimension of the fusion outcome (the dimension of the irrep)

        #Iterate over each fusion outcome to create each trivalent vertex
        for k in [1..Length(currentDecomp)] do
            #Record the label and dimension of the fusion outcome in this iteration
            current_irrep_label := currentDecomp[k];
            irrep_dimension := DegreeOfRepresentation(irrep_dictionary.(current_irrep_label));

            Print("Extracting trivalent vertex for the fusion ", irrepNames[i], "*", irrepNames[j], " into ", current_irrep_label, " which is of dimension ", irrep_dimension, "...", "\n");

            #Create a label for the trivalent vertex for this triplet of irreps in the form "a*b_c"
            trivalent_vertex_label := Concatenation(irrepNames[i], "*", irrepNames[j], "_", current_irrep_label);

            #Pick the appropriate columns of the unitary intertwiner that block diagonalizes irrepNames[i] and irrepNames[j]
            if irrep_dimension = 1 then
                vertex_data := current_unitary_intertwiner{[1..u_size]}{[u_col_index]}; #this extracts only the column at index u_col_index (Note: GAP indexing starts at 1)
            else
                vertex_data := current_unitary_intertwiner{[1..u_size]}{[u_col_index..(u_col_index+irrep_dimension-1)]}; #this extracts all columns from col u_col_index and (u_col_index+irrep_dimension);
            fi;

            #Associate this vertex label with the vertex data in the trivalent vertex record (which works like a dictionary with key/value pairs)
            trivalent_vertices.(trivalent_vertex_label) := vertex_data;

            #Print(trivalent_vertex_label, " is ", Float(vertex_data), "\n");
            u_col_index := u_col_index + irrep_dimension ; #check if col_index was correct

        od;


        #(2)- GET ALL (splitting) TRIVALENT VERTICES for the ORDERED PAIR (irreps_unitary[j],irreps_unitary[i])
        if not identical_irreps then #as long as we have distinct irreps

            current_unitary_intertwiner := u_intertwinerLIST.(Concatenation("U_", irrepNames[j],"*",irrepNames[i])); ##access the unitary intertwiner that block diagonalizes irreps_unitary[j] \tensor irreps_unitary[i] into a direct sum of irreps
            #Print("The BIG unitary block diagonalizing intertwiner for ", irrepNames[j], "*", irrepNames[i], " is: ", Float(current_unitary_intertwiner), "\n");
            u_size := Length(current_unitary_intertwiner); #this is the dimension of the current_unitary_intertwiner matrix

            #Create a variable to record the number of columns of current_unitary_intertwiner which have been extracted
            u_col_index := 1; #each time a trivalent vertex is extracted, this variable will be incremented by the dimension of the fusion outcome (the dimension of the irrep)

            #Iterate over each fusion outcome to create each trivalent vertex
            for k in [1..Length(currentDecomp)] do
                #Record the label and dimension of the fusion outcome in this iteration
                current_irrep_label := currentDecomp[k];
                irrep_dimension := DegreeOfRepresentation(irrep_dictionary.(current_irrep_label));

                Print("Extracting trivalent vertex for the fusion ", irrepNames[j], "*", irrepNames[i], " into ", current_irrep_label, " which is of dimension ", irrep_dimension, "...", "\n");

                #Create a label for the trivalent vertex for this triplet of irreps in the form "a*b_c"
                trivalent_vertex_label := Concatenation(irrepNames[j], "*", irrepNames[i], "_", current_irrep_label);

                #Pick the appropriate columns of the unitary intertwiner that block diagonalizes irrepNames[j] and irrepNames[i]
                if irrep_dimension = 1 then
                    vertex_data := current_unitary_intertwiner{[1..u_size]}{[u_col_index]}; #this extracts only the column at index u_col_index (Note: GAP indexing starts at 1)
                else
                    vertex_data := current_unitary_intertwiner{[1..u_size]}{[u_col_index..(u_col_index+irrep_dimension-1)]}; #this extracts all columns from col u_col_index and (u_col_index+irrep_dimension);
                fi;

                #Associate this vertex label with the vertex data in the trivalent vertex record (which works like a dictionary with key/value pairs)
                trivalent_vertices.(trivalent_vertex_label) := vertex_data;

                #Print(trivalent_vertex_label, " is ", Float(vertex_data), "\n");
                u_col_index := u_col_index + irrep_dimension ; #check if col_index was correct

            od;

        fi;

        #Increment fusionIndex to work with the next fusion rule
        fusionIndex := fusionIndex + 1;

    od;
od;

#Also add all the trivial trivalent vertices (fusion of trivial irrep with irrep x into irrep x itself)
for i in [1..repCount] do
        degree_i := DegreeOfRepresentation(irreps_unitary[i]);
        trivial_label1 := Concatenation("R_I", "*", irrepNames[i], "_", irrepNames[i]); #Fusion of R_I with irrep i into irrep i
        trivial_label2 := Concatenation(irrepNames[i], "*", "R_I", "_", irrepNames[i]); #Fusion of irrep i with R_I into irrep i
        #Both the trivial trivalent vertices are simply the identity map on the representation space of irrep i
        trivalent_vertices.(trivial_label1) := IdentityMat(degree_i);
        trivalent_vertices.(trivial_label2) := IdentityMat(degree_i);
od;


#------------------------------------COMPUTE F SYMBOLS VIA THE GRAPHICAL CALCULUS OF A FUSION CATEGORY---------------------------------
#-----Create a record (dictionary) of fusion rules for easier access----
fusionRules_record := rec(); #dictionary of fusion rules in Rep(G). Keys are pairs of irrep labels in form a*b and values are lists of fusion outcomes of a and b
fusionIndex := 1; #to access each element of fusionRules sequentially
for i in [2..repCount] do
    for j in [i..repCount] do
        currentDecomp := fusionRules[fusionIndex];
        fusion_label1 := Concatenation(irrepNames[i], "*", irrepNames[j]);
        fusion_label2 := Concatenation(irrepNames[j], "*", irrepNames[i]);
        fusionRules_record.(fusion_label1) := currentDecomp;
        fusionRules_record.(fusion_label2) := currentDecomp;
        fusionIndex := fusionIndex + 1;
    od;
od;
#Also add trivial fusion rules into fusionRules_record for convenience
for i in [1..repCount] do
    trivial_rule_label1 := Concatenation("R_I", "*", irrepNames[i]);
    trivial_rule_label2 := Concatenation(irrepNames[i] ,"*", "R_I");
    fusionRules_record.(trivial_rule_label1):= [irrepNames[i]];
    fusionRules_record.(trivial_rule_label2):= [irrepNames[i]];
od;
#-------------------------------------------------------

#----Compute the Fusion Outcomes of Every Distinct Triplet of Irreps of the Group G----
triplet_fusion_outcomes := rec(); #a dictionary of fusion outcomes for every ordered triplet of irreps in Rep(G). Each key is an ordered tuple of non-trivial irreps and the corresponding value is a list of distinct fusion outcomes of the triplet

fusionIndex := 1;
for a in [2..repCount] do #first irrep label a
    for b in [a..repCount] do #second irrep label b
        currentDecomp := fusionRules[fusionIndex]; #contains the fusion outcomes of irrep a * irrep b

        for c in [2..repCount] do #third irrep label c (iterate over all other non-trivial irreps)
            #We will now find the list of all distinct fusion outcomes of a*b*c
            abc_distinct_outcomes := []; #initialize a list to store all distinct outcomes

            for k in [1..Length(currentDecomp)] do #iterate over all fusion outcomes of a and b
                ab_fusion_outcome := currentDecomp[k]; #a specific fusion outcome of a and b

                #Find all the fusion outcomes of irrep c with currentDecomp[k]
                if not ab_fusion_outcome = "R_I" then
                    k_c_outcomeList_label := Concatenation(ab_fusion_outcome, "*", irrepNames[c]); #the label of the list of fusion outcomes of currentDecomp[k] and irrep c
                    ab_k_c_outcomeList := fusionRules_record.(k_c_outcomeList_label); #A list of fusion outcomes of irrep c with currentDecomp[k] accessed from the fusionRules_record record (which works like a dictionary of fusion rules in Rep(G))
                else #if the current fusion outcome is the trivial irrep
                    ab_k_c_outcomeList := [irrepNames[c]]; #the fusion outcome is simply irrep c
                fi;
                for outcome in ab_k_c_outcomeList do
                    if not outcome in abc_distinct_outcomes then #if this is a new a*b*c fusion outcome
                        Add(abc_distinct_outcomes, outcome);
                    fi;
                od;
            od;

            #STORE THE DATA IN THE triplet_fusion_outcomes record (all orders of the triplet {a,b,c} are covered through this implementation)
            triplet_label1 := Concatenation(irrepNames[a], "*", irrepNames[b], "*", irrepNames[c]);
            triplet_fusion_outcomes.(triplet_label1) := abc_distinct_outcomes;
            Print("The outcomes of ", triplet_label1, " are: ", abc_distinct_outcomes, "\n");

            #If irrep a and b are distinct, then store the fusion rule for the swapped pair (b,a) to cover all permutations through which this fusion rule can be accessed
            if not a=b then
                triplet_label2 := Concatenation(irrepNames[b], "*", irrepNames[a], "*", irrepNames[c]);
                triplet_fusion_outcomes.(triplet_label2) := abc_distinct_outcomes;
                Print("The outcomes of ", triplet_label2, " are: ", abc_distinct_outcomes, "\n");

            fi;
        od;

        fusionIndex := fusionIndex + 1; #Increment to work with the next fusion rule
    od;
od;
#--------------------------------------------------

#--------Label and Compute F symbols--------------
Print("Computing the F symbols of Rep(G) through the graphical calculus!...", "\n");
F_symbol_record := rec(); #this is the dictionary of F symbols of Rep(G)
#A key in this record is the label of an F symbol. This label will be a list of 6 labels {a,b,c,e,f,d} which refers to entry (f,d) of the F-matrix [F^{abc}_e] (fusion of irreps a, b and c into e)
#Naturally, the value of each key will be the numerical value of the F symbol.

# We will now methodically LABEL and COMPUTE every non-trivial F symbol of Rep(G)
for a in [2..repCount] do #pick first irrep label irrepNames[a]
    for b in [2..repCount] do #pick second irrep label irrepNames[b]
        for c in [2..repCount] do #pick third irrep label irrepNames[c] (here iterate over all non-trivial irreps)
            #We are working with V^{abc} in this iteration of the loop
            #We want to compute F matrices for the fusion spaces V^{abc}_e for all valid e

            #First we find all valid e (i.e. access all the valid fusion outcomes of a,b and c via the triplet_fusion_record dictionary)
            triplet_label := Concatenation(irrepNames[a], "*", irrepNames[b], "*", irrepNames[c]);
            abc_distinct_outcomes := triplet_fusion_outcomes.(triplet_label); #access all the distinct fusion outcomes of a,b and c

            #Now, we iterate over all distinct fusion outcomes of a,b and c V^{abc}_e for all e
            for abc_outcome in abc_distinct_outcomes do
                #Now we are ready to work with V^{abc}_e for a particular, valid choice of e (e is the "abc_outcome" variable in the current iteration of the loop)
                #We will compute labels and entries of the F matrix F^{abc}_e, but for this we must first find the left and right hand basis of V^{abc}_e
                Print("Extracting the entries of the F matrix F^{", triplet_label, "}_", abc_outcome, "\n");

                #(1)- FIND THE LEFT AND RIGHT HAND FUSION BASIS FOR V^{abc}_e
                leftBasis := []; #will store fusion outcomes of a and b that produce e upon fusion with c (the value of abc_outcome in this loop)
                rightBasis := []; #will store fusion outcomes of b and c that produce e upon fusion with a (the value of abc_outcome in this loop)

                #COMPUTING LEFT BASIS- Fuse a and b first
                ab_outcomes := fusionRules_record.(Concatenation(irrepNames[a], "*", irrepNames[b])); #access the fusion outcomes of a*b
                for ab_outcome in ab_outcomes do #for an outcome d of a*b (d is the iteration variable ab_outcome in this loop)
                    d_c_label := Concatenation(ab_outcome, "*", irrepNames[c]); #the dictionary key for the list of fusion outcomes of d*c
                    dc_outcomeList := fusionRules_record.(d_c_label); #access the fusion outcomes of irrep d and c
                    if abc_outcome in dc_outcomeList then #if e (abc_outcome in this loop) appears amongst the fusion outcomes of d*c (where d is ab_outcome in this loop)
                        Add(leftBasis, ab_outcome); #add d (ab_outcome in this loop) to leftBasis (left hand basis list of V^{abc}_e)
                    fi;
                od;

                #COMPUTING RIGHT BASIS- Fuse a and b first
                bc_outcomes := fusionRules_record.(Concatenation(irrepNames[b], "*", irrepNames[c])); #access the fusion outcomes of a*b
                for bc_outcome in bc_outcomes do #for an outcome f of b*c (f is the iteration variable bc_outcome in this loop)
                    a_f_label := Concatenation(irrepNames[a], "*", bc_outcome); #the dictionary key for the list of the fusion outcomes of a*f
                    af_outcomeList := fusionRules_record.(a_f_label); #access the list of fusion outcomes of irrep a and f
                    if abc_outcome in af_outcomeList then #if e (abc_outcome in this loop) appears amongst the fusion outcomes of a*f (where f is bc_outcome in this loop)
                        Add(rightBasis, bc_outcome); #add f (bc_outcome in this loop) to rightBasis (right hand basis list of V^{abc}_e)
                    fi;
                od;

                #Print("The left hand basis: ", leftBasis, "\n");

                #Print("The right hand basis: ", rightBasis, "\n");


                #(2)-Compute all elements of the F matrix F^{abc}_e. The key for the entry (f,d) of F^{abc}_e will be represented as the string [a,b,c,e,f,d] (this order may look strange, but it becomes alphabetic in the binary tree representation)
                for label_d in leftBasis do
                    for label_f in rightBasis do
                        #Streamline the labels of all irreps
                        #We are computing the F matrix {a,b,c)
                        irrep_a := irrepNames[a];
                        irrep_b := irrepNames[b];
                        irrep_c := irrepNames[c];
                        irrep_e := abc_outcome; #a particular outcome of fusion of a,b and c
                        irrep_d := label_d; #a fusion outcome of a*b which can fuse with c to produce e
                        irrep_f := label_f; #a fusion outcome of b*c which can fuse with a to produce e
                        F_symbol_label := Concatenation(irrep_a, ",", irrep_b, ",", irrep_c, ",", irrep_e, ",", irrep_f ,",", irrep_d);

                        #We first compute the following binary tree (DIAGRAM 1):
                        #Note: Time flows bottom to top. If morphisms f1,f2 are represented by diagrams d1, d2 respectively then stacking d2 on d1 represents morhpism "f2 of f1"
                        #         |
                        #         e
                        #        / \
                        #       d   \
                        #      / \   \
                        #     a   b   c
                        c_degree := DegreeOfRepresentation(irreps_unitary[c]); #dimension of irrep c
                        c_identity := IdentityMat(c_degree); #identity map on the representation space of irrep c
                        ab_d_fusion_vertex := TransposedMat(ComplexConjugate(trivalent_vertices.(Concatenation(irrep_a, "*", irrep_b, "_", irrep_d)))); #Hermitian Conjugate of the Splitting Vertex d -> ab
                        dc_e_fusion_vertex := TransposedMat(ComplexConjugate(trivalent_vertices.(Concatenation(irrep_d, "*", irrep_c, "_", irrep_e)))); #Hermitian Conjugate of the Splitting Vertex e -> dc
                        diagram1 := dc_e_fusion_vertex*KroneckerProduct(ab_d_fusion_vertex, c_identity);

                        #Now we compute the following binary tree (DIAGRAM 2):
                        #         |
                        #         e
                        #        / \
                        #       /   f
                        #      /   / \
                        #     a   b   c
                        a_degree := DegreeOfRepresentation(irreps_unitary[a]); #dimension of irrep c
                        a_identity := IdentityMat(a_degree);#identity map on the representation space of irrep c
                        bc_f_fusion_vertex := TransposedMat(ComplexConjugate(trivalent_vertices.(Concatenation(irrep_b, "*", irrep_c, "_", irrep_f))));  #Hermitian Conjugate of the Splitting Vertex f -> bc
                        af_e_fusion_vertex := TransposedMat(ComplexConjugate(trivalent_vertices.(Concatenation(irrep_a, "*", irrep_f, "_", irrep_e))));  #Hermitian Conjugate of the Splitting Vertex e -> af
                        diagram2 := af_e_fusion_vertex*KroneckerProduct(a_identity, bc_f_fusion_vertex);


                        #The F symbol can be computed by stacking diagram 2 upon the hermitian conjugate of diagram 1
                        #This stacked diagram represents the identity map on the representation space of irrep e multiplied by [F^{abc}_e]_fd (proved via rules of graphical calculus in a unitary fusion category i.e. anyon diagrammatics)
                        #         |
                        #         e                           |
                        #        / \                          |
                        #       /   f                         |
                        #      /   / \                        |
                        #     a   b   c   = [F^{abc}_e]_fd    e
                        #      \ /   /                        |
                        #       d   /                         |
                        #        \ /                          |
                        #         e                           |
                        #         |
                        stacked_diagram := diagram2*TransposedMat(ComplexConjugate(diagram1));

                        #Thus, the F symbol is simply the trace of stacked_diagram divided by the dimension of irrep e
                        e_degree := DegreeOfRepresentation(irrep_dictionary.(irrep_e));
                        F_symbol_value := TraceMatrix(stacked_diagram)/e_degree;

                        #Store the information in the F_symbol_data matrix record/display it
                        F_symbol_record.(F_symbol_label) := F_symbol_value;
                        Print("The F symbol ", F_symbol_label, " is ", EvaluateCyclotomic(F_symbol_value), "\n");
                    od;
                od;
            od;
        od;
    od;
od;

#-----------------------------STORE INFO IN A TEXT FILE-------------------
