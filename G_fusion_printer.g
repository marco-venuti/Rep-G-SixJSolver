#-------------------DEFINING GROUP, LOADING CHARACTER TABLE, OVERHEAD INFO---------------
LoadPackage("ctbllib");
#DEFINE THE GROUP HERE
G := SmallGroup(3,1);

#Find the character table and the irreducible characters
charTable := CharacterTable(G);
irreducibleCharacters := Irr(charTable);

#Validating the size of character table loaded
repCount:= Length(ConjugacyClasses(G)); #we expect as many irreps as the number of conj. classes
if not repCount = Length(irreducibleCharacters) then
    Print("WARNING: Character Table was not loaded correctly \n");
fi;

#PRINT OVERHEAD INFORMATION
Print("This is the fusion ring data for ", G, ".", "\n");
Print("Note on Irrep Labels: The label R_d will be ascribed to an irrep of dimension d. If there are more than one irreps of dimension d, a further (arbitrary) indexing R_d_j will be used to distinguish between them. R_I refers to the trivial irrep. \n");
Print("This data will be stored in a text file named \"fusion_data.txt\" in the current directory. \n \n");

#----------------------------NAMING IRREDUCIBLE REPRESENTATIONS--------------------------------

#STEP 1: SORT the irreducible characters by degree
sortedChars := SortedCharacters(charTable, irreducibleCharacters, "degree"); #trivial character is first by default


#STEP 2: Record degree of each character with the same degree to create indexing later
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
    	Add(irrepNames, "{R_I}"); #label of the trivial irrep

   	else
    	#Labelling General Irrep
    	irrepLabel := Concatenation("{R_", String(irrepDegree[i]), "}"); #Label the irrep by the dimension

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



# Display the names for verification
Print("Irrep Labels: $\\{");
for i in [1..Length(irreducibleCharacters)] do
    if i = Length(irreducibleCharacters) then
        Print(irrepNames[i], "\\}$", "\n");
    else
        Print(irrepNames[i], ", ");
    fi;
od;



#--------------------------------COMPUTING FUSION RULES----------------------------------------
# Initialize an empty list to store fusion rules
fusionRules := [];

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

        # Store the fusion rule
        decompString := "";
        for irrep in [1..Length(decomposition)] do
            if irrep=Length(decomposition) then
                decompString := Concatenation(decompString, decomposition[irrep]);
            else
                decompString := Concatenation(decompString, decomposition[irrep], " + ");
            fi;
        od;
        #Print(decompString, "\n");
        Add(fusionRules, Concatenation(irrepNames[i], " \\times ", irrepNames[j], " = ", decompString));

    od;
od;

# Display fusion rules
Print("\\newline \n", "The non-trivial fusion rules are: \n");
for rule in fusionRules do
    Print("$$", rule, "$$", "\n");
od;



#-------------------------PRINTING DATA IN A TEXT FILE------------------------------------
# Open a file for writing the output data
filename:= "fusion_data.txt";;
stream := OutputTextFile(filename, false);

#OVERHEAD INFORMATION
if not repCount = Length(irreducibleCharacters) then
    PrintTo(stream, "WARNING: Character Table was not loaded correctly \n");
fi;
PrintTo(stream, "This is the fusion ring data for ", G, ".", "\n", "\n");
PrintTo(stream, "NOTE on Irrep Labels: The label R_d will be ascribed to an irrep of dimension d. If there are more than one irreps of dimension d, a further (arbitrary) indexing R_d_j will be used to distinguish between them. R_I refers to the trivial irrep. \n \n");

# IRREP LABELS
PrintTo(stream, "Irrep Labels: $\\{");
for i in [1..Length(irreducibleCharacters)] do
    if i = Length(irreducibleCharacters) then
        PrintTo(stream, irrepNames[i], "\\}$", "\n");
    else
        PrintTo(stream, irrepNames[i], ", ");
    fi;
od;

# FUSION DATA
PrintTo(stream, "\\newline \n", "The non-trivial fusion rules are: \n");
for rule in fusionRules do
    PrintTo(stream, "$$", rule, "$$", "\n");
od;

# Close the file stream
CloseStream(stream);
