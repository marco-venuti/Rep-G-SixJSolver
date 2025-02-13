
#-------------------DEFINING GROUP, LOADING CHARACTER TABLE, OVERHEAD INFO---------------
LoadPackage("ctbllib");
#DEFINE THE GROUP HERE
G:=SmallGroup(6,1);
#G := SymmetricGroup(4);

#Find the character table and the irreducible characters
charTable := CharacterTable(G);
irreducibleCharacters := Irr(charTable);

#Validating the size of character table loaded
repCount:= Length(ConjugacyClasses(G)); #we expect as many irreps as the number of conj. classes
if not repCount = Length(irreducibleCharacters) then
    Print("WARNING: Character Table was not loaded correctly \n");
fi;

#----------------------------NAMING IRREDUCIBLE REPRESENTATIONS--------------------------------

#STEP 1: SORT the irreducible characters by degree in ascending order
sortedChars := SortedCharacters(charTable, irreducibleCharacters, "degree"); #trivial character is first by default


#STEP 2: Record degree of each character to decide whether to add secondary indexing later
irrepDegree := [];
for i in [1..Length(sortedChars)] do
    degree := Degree(sortedChars[i]);
    Add(irrepDegree, degree);
od;

#STEP 3: The label R_d will be ascribed to an irrep of dimension d. If there are more than one irreps of dimension d, a further (arbitrary) indexing R_d_j will be used to distinguish between them. R_I will refer to the trivial irrep.

irrepNames:=[]; #This list will contain the appropriate names of the irreps in the same order as sortedChars

#Iterate Over Each Irred. Character and Label it
countWithinDimension := 1; #maintains the count of irreps of the same degree
for i in [1..Length(sortedChars)] do

    #Labelling Trivial Irrep
    if i=1 then #the first irrep is always the trivial irrep
    	Add(irrepNames, "R_I"); #label of the trivial irrep

   	else
    	#Labelling General Irrep
    	irrepLabel := Concatenation("R_", String(irrepDegree[i])); #Label the irrep by the dimension

    	#Handle the indexing for irreps of the same degree
    	if Degree(sortedChars[i-1]) < irrepDegree[i] then
    		countWithinDimension := 1; #first irrep of this dimension found
    	else
    		countWithinDimension := countWithinDimension + 1; #increment the number of irreps of this dimension
    	fi;

    	#Add further indexing to irrep label if multiple irreps of the same degree exist
    	if countWithinDimension > 1 then
    		irrepLabel := Concatenation(irrepLabel, "_", String(countWithinDimension));
    	else
    		#If countWithinDimension is 1, but there are more to follow
    		if IsBound(sortedChars[i+1]) and Degree(sortedChars[i+1]) = irrepDegree[i] then
    			irrepLabel := Concatenation(irrepLabel, "_", String(countWithinDimension));
    		fi;
    	fi;

    	Add(irrepNames, irrepLabel); #Add this label to the list
    fi;
od;



# Display the irrep names for verification
Print("Irrep Labels: ");
for i in [1..Length(irreducibleCharacters)] do
    if i = Length(irreducibleCharacters) then
        Print(irrepNames[i], "\n");
    else
        Print(irrepNames[i], ", ");
    fi;
od;



#--------------------------------COMPUTING FUSION RULES----------------------------------------
# Initialize an empty list to store fusion rules
fusionRules := []; #Each entry will be a list of the products of fusion

for i in [2..Length(sortedChars)] do
    for j in [i..Length(sortedChars)] do
        chi := sortedChars[i];
        psi := sortedChars[j];

        # Compute the product of the characters
        product := chi * psi;

        # Decompose the product into irreducible characters
        decomposition := [];
        for k in [1..Length(sortedChars)] do
            phi := sortedChars[k];
            coefficient := ScalarProduct(charTable, product, phi);
            if coefficient > 0 then
                if coefficient = 1 then
                    Add(decomposition, irrepNames[k]);
                else
                    Add(decomposition, Concatenation(String(coefficient), irrepNames[k]));
                fi;
            fi;
        od;

        Add(fusionRules, decomposition);

    od;
od;

#---------------------------- PRINT THE FUSION RULES ---------------------------------
Print("The non-trivial fusion rules are: \n");
repCount := Length(sortedChars); #the total number of irreducible characters/representations of G
fusionIndex := 1; #an index for the fusionRules list
for i in [2..repCount] do
    for j in [i..repCount] do

        currentDecomp := fusionRules[fusionIndex]; #a list of fusion outcomes for current pair of irreps

        Print(irrepNames[i], " * ", irrepNames[j], " = ");
        for k in [1..Length(currentDecomp)] do
            if k =  Length(currentDecomp) then
                Print(currentDecomp[k], "\n");
            else
                Print(currentDecomp[k], " + ");
            fi;
        od;
        fusionIndex := fusionIndex + 1; #increment the index to print the next fusion rule
    od;

od;
