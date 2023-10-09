#include <iostream>
#include <vector>
#include "PhenotypeIntegration.h"
using namespace std;

int main()
{
    PhenotypeIntegration phen;

    // Read user input, e.g., file with patirent's phenotypes HP:0001249, HP:0000003, ..."
    phen.read_user_HPOs("./Input/user_input_patient_phenotypes.csv"); 
    
    // Read external table with gene variants, fetching the GeneNames, e.g., table with some colums for GeneName and ACMG points
    phen.read_variant_table("./Input/gene_variants_table.csv"); 

    // Set up the HPO resource, precomputing Information Content for all HPO nodes in the tree, etc.
    // phen.setupHPO("./HPO/HPO_release20220414/"); 
    phen.setupHPO("./HPO/HPO_release20230127/"); 
    
    phen.set_buildOut("");
    string str = phen.get_buildOut();

    return 0;
}
