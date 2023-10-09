#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <math.h>
#include <algorithm>
#include "PhenotypeIntegration.h"
using namespace std;



void PhenotypeIntegration::checkQueryInHPOdata(vector<string> HPOSet)
{

    // Make sure there are no duplicated entries, as they would influence the calculation of semantic similarity
    sort(HPOSet.begin(), HPOSet.end());
    for (unsigned int i = 0; i < HPOSet.size(); i++)
    {
        for (unsigned int ii = i + 1; ii < HPOSet.size(); ii++)
        {
            if (HPOSet[i] == HPOSet[ii])
            {
                cerr<<"HPO query term present multiple times in set: " + HPOSet[i]<<endl;
                exit(111);
            }
        }

        // Also check that the user does not choose the root of a subontology with no gene associactions (e.g., HP:0032223 "Blood group")
        if (find(vUnderontGenTerm.begin(), vUnderontGenTerm.end(), (HPOSet[i])) != vUnderontGenTerm.end())
        {
            cerr<<"Sub-ontology root term chosen with no associations: " + HPOSet[i]<<endl;
            exit(222);
        }
    }

    for (unsigned int i = 0; i < HPOSet.size(); i++)
    {
        // cout << " " << HPOSet[i]<<endl;
        if ((mHPOandPars.find(HPOSet[i]) == mHPOandPars.end()) &&
            (mHPOandKids.find(HPOSet[i]) == mHPOandKids.end()))
        {
            // not found
            cerr<<"HPO query term not present in the HPO data released currently employed: " + HPOSet[i]<<endl;
            exit(333);
        }

        if (HPOSet[i] == "HP:0000001")
        {
            // not found
            cerr<<"HP:0000001 (HPO term 'All') is the overall root over the whole ontology, not a good idea to select this node."<<endl;
            exit(444);
        }
    }
}




void PhenotypeIntegration::readVarTab()
{
    // cout << endl << "Hello readVarTab" << endl;
    ifstream file(inpVarTab);
    if (file.is_open())
    {
        // Skip the header
        string header_genes_to_phenotype;
        getline(file, header_genes_to_phenotype);

        // Read data from the file
        string line;
        while (getline(file, line))
        {
            istringstream iss(line);

            string col1, col2, col3, col4, col5, col6, col7;
            double col8;
            if (getline(iss, col1, ';') &&
                getline(iss, col2, ';') &&
                getline(iss, col3, ';') && // assumes geneName is here
                getline(iss, col4, ';') &&
                getline(iss, col5, ';') &&
                getline(iss, col6, ';') &&
                getline(iss, col7, ';') &&
                (iss >> col8))  // assumes ACMG or other score of interest is here
            {
                vInpVarTabGenName.push_back(col3);
                vInpVarTabScore.push_back(col8);
            }
            else
            {
                cerr << "Failed to read line: " << line << endl;
            }
        }
        // no need to pop.back last entry here
        file.close();
    }
    else
    {
        cerr << "Unable to open file: " << inpVarTab << endl;
    }
}



void PhenotypeIntegration::readUserHPOs()
{
    // cout << endl << "Hello readUserHPOs" << endl;
    ifstream file(inpUserHPOs);
    if (file.is_open())
    {
        // Read data from the file
        string line;
        while (getline(file, line))
        {
            istringstream iss(line);

            string col1;
            if (getline(iss, col1, ' '))
            {
                vInpUserHPOs.push_back(col1);
            }
            else
            {
                cerr << "Failed to read line: " << line << endl;
            }
        }
        // no need to pop.back last entry here
        file.close();
    }
    else
    {
        cerr << "Unable to open file: " << inpUserHPOs << endl;
    }
}



void PhenotypeIntegration::readObo()
{
    // cout << endl << "Hello readObo" << endl;
    ifstream file(inpObo);
    if (file.is_open())
    {
        string str;
        while (!file.eof())
        {
            file >> str;
            vHPOs.push_back(str);
        }
        vHPOs.pop_back(); // this is needed to pop the overall last entry included twice
        file.close();
        frameObo();
    }
    else
    {
        cerr << "Unable to open file: " << inpObo << endl;
    }
}

void PhenotypeIntegration::frameObo()
{
    // cout << endl << "Hello frameObo" << endl;

    // 1. Make a map binary tree structure for efficient tree walking UP (nodes -> root)
    vector<string> tmp_v_parent;
    string tmp_node;
    int ii;

    for (unsigned int i = 0; i < vHPOs.size(); i++)
    {
        if (vHPOs[i] == "[Term]")
        { // catch the first [Term] flag (start)
            ii = i + 1;
            while (vHPOs[ii] != "[Term]" && vHPOs[ii] != "[Typedef]" && ii != vHPOs.size() - 1)
            { // while loop that moves forward till the end of the block
                // get the node name
                if (vHPOs[ii] == "id:")
                {
                    tmp_node = vHPOs[ii + 1];
                }
                // get the node parent
                if (vHPOs[ii] == "is_a:")
                {
                    tmp_v_parent.push_back(vHPOs[ii + 1]);
                }
                ii++;
            }
            mHPOandPars.insert(pair<string, vector<string>>(tmp_node, tmp_v_parent));
            tmp_v_parent.clear();
        } // END if [TERM]
    }

    // 2. Invert the above map for tree walking DOWN (nodes -> leaves)

    // Iterate through the original map and build the inverted map
    for (const auto &kv : mHPOandPars)
    {
        const string &key = kv.first;
        const vector<string> &values = kv.second;
        // Iterate through the vector of values
        for (const auto &value : values)
        {
            // Insert the inverted key-value pair into the inverted map
            mHPOandKids[value].push_back(key);
        }
    }
}

void PhenotypeIntegration::readGenToPhen()
{
    // cout << endl << "Hello readGenToPhen" << endl;
    ifstream file(inpGenToPhen);
    if (file.is_open())
    {
        // Skip the header
        string header_genes_to_phenotype;
        getline(file, header_genes_to_phenotype);

        // Read data from the file
        string line;
        while (getline(file, line))
        {
            istringstream iss(line);

            string col1, col2, col3, col4, col5, col6, col7, col8, col9;
            if (getline(iss, col1, '\t') && // entrez-gene-id
                getline(iss, col2, '\t') && // gene_symbol (some can have ncbi_gene_id but no gene_symbol ("-"))
                getline(iss, col3, '\t') && // HPO-Term-ID
                getline(iss, col4, '\t') && // HPO-Term-Name
                getline(iss, col5, '\t') &&
                getline(iss, col6, '\t') &&
                getline(iss, col7, '\t') &&
                getline(iss, col8, '\t') &&
                getline(iss, col9, '\t')) // disease-ID 
            {
                vGenToPhenGenID.push_back(col1); 
                vGenToPhenGenName.push_back(col2);
                vGenToPhenDisID.push_back(col9);
                vGenToPhenHPOterm.push_back(col3);
                replace(col4.begin(), col4.end(), ' ', '_');
                vGenToPhenHPOtermName.push_back(col4);
            }
            else
            {
                cerr << "Failed to read line: " << line << endl;
            }
        }
        // no need to pop.back last entry here
        file.close();
        frameGenToPhen();
    }
    else
    {
        cerr << "Unable to open file: " << inpGenToPhen << endl;
    }
}

void PhenotypeIntegration::frameGenToPhen()
{
    // cout << endl << "Hello frameGenToPhen" << endl;

    // Make a map binary tree structure for efficient grep/count genes mapped to HPO terms

    // Iterate through the HPO vector and insert its elements as keys and empty vectors as values in the map
    for (const auto &key : vGenToPhenHPOterm)
    {
        mHPOsGen[key] = vector<string>();
    }
    // Iterate through the GeneID vector and insert its elements
    for (unsigned int i = 0; i < vGenToPhenDisID.size(); i++)
    {
        mHPOsGen[vGenToPhenHPOterm[i]].push_back(vGenToPhenDisID[i]);
    }

    // Add a map that assosciates to the GeneName, the DiseaseOMIM(s), and the HPOs associated to the DiseaseOMIMs. These are the columns 2 and 9 of the HPO file "genes_to_phenotype.txt" (warning about HPO release 20230405 has 4 columns, before and after there are correctly 9 columns).
    // This way the system can receive only the gene name from the gene variant table, and it associates HPO terms corresponding to it according to a given HPO release, thereby avoiding the need for the input gene variant table or else to have genomic variants annotated with HPO terms.

    // Assuming that the HPO file "genes_to_phenotype.txt" is a list sorted by geneIDs/geneNames, i.e., that there are succeding blocks of the same gene one after the other, as this:
    // 8195 MKKS ORPHA:110 HP:0010747
    // 8195 MKKS ORPHA:110 HP:0000580
    // 8195 MKKS ORPHA:110 HP:0002230
    // ...
    // 8195 MKKS OMIM:236700 HP:0000145
    // 8195 MKKS OMIM:236700 HP:0001374
    // 8195 MKKS OMIM:236700 HP:0000113
    // ...
    // 253959 RALGAPA1 OMIM:618797 HP:0002059
    // 253959 RALGAPA1 OMIM:618797 HP:0000518
    // 253959 RALGAPA1 OMIM:618797 HP:0000007
    // then it associated GeneName-DiseaseID to the corresponding variable number of HPOs associated.

    vGenToPhenGenName.push_back("endOfTheVectorToAddForLastGeneNameBlock"); // This is only to make it read the very last block GeneID-Disease, otherwise in this loop-way comparison it would skip the last row/block.
    for (unsigned int i = 0; i < vGenToPhenGenName.size() - 1; i++)
    {
        // Conditional for same GeneName and same DiseaseID, then collect corresponding HPOs
        if ((vGenToPhenGenName[i] == vGenToPhenGenName[i + 1]) && (vGenToPhenDisID[i] == vGenToPhenDisID[i + 1]))
        {
            // cout<<" "<<vGenToPhenGenID[i];
            // cout<<" "<<vGenToPhenGenName[i];
            // cout<<" "<<vGenToPhenDisID[i];
            // cout<<" "<<vGenToPhenHPOterm[i]<<endl;
            vGenToPhenDisHPOs.push_back(vGenToPhenHPOterm[i]);
        }
        else
        {
            // read last HPO before the switch to next
            // cout<<" "<<vGenToPhenHPOterm[i];
            vGenToPhenDisHPOs.push_back(vGenToPhenHPOterm[i]);

            // push back also at the very end after HPO terms: GeneName, DiseaseID(OMIM/ORPHA), Disease Description with no space. They can be removed simply here below if not needed.
            vGenToPhenDisHPOs.push_back(vGenToPhenGenName[i]);
            vGenToPhenDisHPOs.push_back(vGenToPhenDisID[i]);
            vGenToPhenDisHPOs.push_back(vGenToPhenHPOtermName[i]);

            // all that goes in the 2-d vector, to further be associated to a geneName (it can be that a single geneName has multiple diseases associated, each with its HPOs)
            vvGenToPhenDisHPOs.push_back(vGenToPhenDisHPOs);
            vGenToPhenDisHPOs.clear();

        } // else END

        // finally put all stuff into map
        // cout<<"---> "<<vGenToPhenGenName[i]<<"-"<<vGenToPhenDisID[i]<<endl;
        // the map contains vectors in which the last elements are the DiseaseOMIM/ORPHA (also the GeneName, which is the key of the map) preceded by the corresponding HPO terms (DiseaseID and GeneName could be left out if necessary, but since they are unique and present once, they can be put at the end of the HPO term list).
        mGeneNamDiseHPOs[vGenToPhenGenName[i]] = vvGenToPhenDisHPOs;

        // vvGenToPhenDisHPOs.clear(); // not this way, othewise it does not accumulate different diseases-per-same-gene, but only the last one overwriting the previous ones, so:
        if ((vGenToPhenGenName[i] != vGenToPhenGenName[i + 1]))
        {
            vvGenToPhenDisHPOs.clear();
        }

    } // for END

    //     // Print the map
    //     for (const auto& pair : mGeneNamDiseHPOs) {
    //         const string& key = pair.first;
    //         const vector<vector<string>>& value = pair.second;

    //         cout << "Key: " << key << endl;

    //         for (const auto& innerVec : value) {
    //             for (const string& str : innerVec) {
    //                 cout << str << " ";
    //             }
    //             cout << endl;
    //         }
    //             cout << endl;
    //    }

    // sort and unique
    sort(vGenToPhenDisID.begin(), vGenToPhenDisID.end());
    vGenToPhenDisID.erase(unique(vGenToPhenDisID.begin(), vGenToPhenDisID.end()), vGenToPhenDisID.end());
    set_preCompICs();
}

void PhenotypeIntegration::preCompICs()
{
    // cout << "Hello preCompICs" << endl;

    cout << endl
         << "Precomputing IC values for all HPO terms" << endl
         << endl;

    // 0. Open an outputfile named with the release tag
    string name_outfile_IC;
    vector<pair<string, vector<double>>> v_fingerprint;

    // Extract the "release" part from inpObo to append to the IC file generated here
    // Split the path into parts using slash as delimiter
    // --------------------------------
    string path = inpObo;
    vector<string> parts;
    size_t start = 0, end = 0;
    while ((end = path.find("/", start)) != string::npos)
    {
        if (end != start)
        {
            string part = path.substr(start, end - start);
            parts.push_back(part);
        }
        start = end + 1;
    }
    if (end != start)
    {
        string part = path.substr(start);
        parts.push_back(part);
    }
    // Extract the desired part using its index
    int index = 2; // index of "documents"
    if (index >= 0 && index < parts.size())
    {
        // cout << parts[index] << endl;
        // cout << "release: IC_" << parts[index] << endl << endl;
        name_outfile_IC = "IC_" + parts[index] + ".dat";
    }
    // --------------------------------END

    // 1. get all sub-ontologies, e.g., for HPO_release20220414:
    // HP:0000005 HP:0000118 HP:0012823 HP:0032223 HP:0032443 HP:0040279

    map<string, vector<string>>::iterator it;
    it = mHPOandKids.find("HP:0000001"); // start from root and get the sub-ontology roots one level below
    // cout << (*it).first << endl;

    // Add the overall root to the list before dealing with the sub-ontologies
    // --------------------------------
    string term = (*it).first;
    cout << "Overall ontology root: " << term << ", number of genes associated (N) = ";
    treeWalkDown(term);
    set_uniqNodAllDesc();
    set_findHPOassGen(term);
    double N_genes_in_ontology = get_findHPOassGen();
    // // Change here too to make it like HPOSim:
    // N_genes_in_ontology = 8234;
    cout << N_genes_in_ontology << endl;
    vector<string> _str = get_treeWalkDown();
    // cout <<term << " ";

    vPathDown.clear();
    vPathUp.clear(); // clear also the other path vector otherwise interference!

    treeWalkDown(term);
    set_uniqNodAllDesc();
    set_findHPOassGen(term);
    double n = get_findHPOassGen();
    double IC = -log(n / N_genes_in_ontology);
    // // cout<< "IC = -ln(n/N) = -ln("<<n<<"/"<<N_genes_in_ontology<<") = ";
    // cout << IC << endl;

    vector<double> tmp_v_fingerprint{n, N_genes_in_ontology, IC};
    v_fingerprint.push_back(pair<string, vector<double>>(term, tmp_v_fingerprint));
    // --------------------------------END

    for (unsigned int j = 0; j < (*it).second.size(); j++)
    {
        string term = (*it).second[j];

        // 2. Count how many genes are associated to each sub-ontology
        cout << "Sub-ontology root: " << term << ", number of genes associated (N) = ";
        treeWalkDown(term);
        set_uniqNodAllDesc();
        set_findHPOassGen(term);
        double N_genes_in_subontology = get_findHPOassGen();
        // If needed to assume N being the overall HPO terms without bothering subdivision in sub ontologies:
        // N_genes_in_subontology = N_genes_in_ontology;
        // // Change this to make it equal to HPOSim (in which they assume N_genes_in_subontology (HP:0000118) is the number of DiseaseOMIM uniq, not the genes! (8234) instead of 4785):
        // N_genes_in_subontology = 8234;

        cout << N_genes_in_subontology << endl;

        // 3. Walk down and calculate for each HPO term the IC
        vector<string> _str = get_treeWalkDown();

        if (N_genes_in_subontology <= 1)
        { // if sub-ontology has 1 or less genes associated then skip because makes no sense
            cout << "No genes associated to sub-ontology of " << term << ", skip IC computation." << endl
                 << endl;

            // add the zero-genes-associated subontology to a vector
            vUnderontGenTerm.push_back(term);
        }
        else
        {
            cout << "Computing IC for all HPO terms..." << endl
                 << endl;
            for (unsigned int i = 0; i < _str.size(); i++)
            {
                // cout << _str[i] << " ";

                vPathDown.clear();
                vPathUp.clear();

                treeWalkDown(_str[i]);
                set_uniqNodAllDesc();
                set_findHPOassGen(_str[i]);
                double n = get_findHPOassGen();
                double IC = -log(n / N_genes_in_subontology);
                // // cout<< "IC = -ln(n/N) = -ln("<<n<<"/"<<N_genes_in_subontology<<") = ";
                // cout << IC << endl;

                vector<double> tmp_v_fingerprint{n, N_genes_in_subontology, IC};
                v_fingerprint.push_back(pair<string, vector<double>>(_str[i], tmp_v_fingerprint));
            }
        }

    } // end loop over sub-ontology roots

    // sort and unique
    sort(v_fingerprint.begin(), v_fingerprint.end());
    v_fingerprint.erase(unique(v_fingerprint.begin(), v_fingerprint.end()), v_fingerprint.end()); // erase the non-unique

    // open output file
    ofstream outfile_IC;
    outfile_IC.open(name_outfile_IC);

    // Store into a map for faster retrieval
    for (unsigned int i = 0; i < v_fingerprint.size(); i++)
    {
        mHPOsSign.insert(pair<string, vector<double>>(v_fingerprint[i].first, v_fingerprint[i].second));

        // write to outputfile
        outfile_IC << v_fingerprint[i].first << " ";
        outfile_IC << "IC = -ln(n/N) = -ln(" << v_fingerprint[i].second[0] << "/" << v_fingerprint[i].second[1] << ") = ";
        outfile_IC << v_fingerprint[i].second[2] << endl;
    }
    outfile_IC.close();
}



void PhenotypeIntegration::findHPOassGen(string qHPO)
{
    // cout << "Hello grepHPOforCorrespondingGene" << endl;

    // For all the HPO terms contained in a vector
    // (vPathDownSortUniq),
    // search all corresponding genes in "genes_to_phenotype.txt"
    // (which have been stored in a map HPO-GeneID for faster search)

    // The following is equivalent to something like this in bash:
    // grep -Ff list_HPO_terms_descendants genes_to_phenotype.txt | cut -f "$geneID" | sort | uniq | wc

    set_treeWalkDown(qHPO); // this generates vPathDownSortUniq
    map<string, vector<string>>::iterator it;
    vector<string> v_geneID;
    for (unsigned int i = 0; i < vPathDownSortUniq.size(); i++)
    {
        it = mHPOsGen.find(vPathDownSortUniq[i]);
        if (mHPOsGen.count(vPathDownSortUniq[i]) == 0)
        {
            // base case
        }
        else
        {
            for (unsigned int j = 0; j < (*it).second.size(); j++)
            {
                v_geneID.push_back((*it).second[j]);
            }
        }
    }
    // sort and unique
    sort(v_geneID.begin(), v_geneID.end());
    v_geneID.erase(unique(v_geneID.begin(), v_geneID.end()), v_geneID.end());
    // cout << "uniq = " << v_geneID.size()<< endl;
    // cout << "" << endl;
    mapGenToHPOkidz = v_geneID.size();
}

void PhenotypeIntegration::printTermPar()
{
    // cout << endl << "Hello printTermPar" << endl;
    map<string, vector<string>>::iterator it;
    // Tree strucure with node-parents
    cout << "Nodes -> Parents" << endl;
    for (it = mHPOandPars.begin(); it != mHPOandPars.end(); it++)
    {
        // Print current map content
        cout << (*it).first << " -> ";
        for (unsigned int j = 0; j < (*it).second.size(); j++)
        {
            cout << " " << (*it).second[j];
        }
        cout << "" << endl;
    }
}


void PhenotypeIntegration::printTermKids()
{
    // cout << endl << "Hello printTermKids" << endl;
    map<string, vector<string>>::iterator it;
    // Tree strucure with node-children
    cout << "Nodes -> Children" << endl;
    for (it = mHPOandKids.begin(); it != mHPOandKids.end(); it++)
    {
        // Print current map content
        cout << (*it).first << " -> ";
        for (unsigned int j = 0; j < (*it).second.size(); j++)
        {
            cout << " " << (*it).second[j];
        }
        cout << "" << endl;
    }
}


// Prints out delta IC between parent term and children, the value can then be mapped to the link in a visualization of the HPO tree
void PhenotypeIntegration::printInpNetForCytoscape()
{
    double IC1;
    double IC2;
    // cout << endl << "Hello printInpNetForCytoscape" << endl;
    map<string, vector<string>>::iterator it;
    for (it = mHPOandKids.begin(); it != mHPOandKids.end(); it++)
    {
        // Print current map content
        for (unsigned int j = 0; j < (*it).second.size(); j++)
        {
            // cout << (*it).first << " + " << (*it).second[j] << endl;
            IC1 = get_InformationContent((*it).first);
            IC2 = get_InformationContent((*it).second[j]);
            // cout << (*it).first << " " << IC2 << " + " << IC1 << " " << (*it).second[j] << endl;
            cout << (*it).first << " " << IC2 - IC1 << " " << (*it).second[j] << endl;
            // cout << endl;
        }
    }
}


void PhenotypeIntegration::queryTermForPar(string qHPO)
{
    // cout << endl << "Hello queryTermForPar" << endl;

    // Fast queries of node and returns parents:
    // cout << "Query (Node): " << qHPO << endl;
    // cout << "Parents: " << endl;
    map<string, vector<string>>::iterator it;
    it = mHPOandPars.find(qHPO);
    // cout << (*it).first << endl;
    for (unsigned int j = 0; j < (*it).second.size(); j++)
    {
        if (mHPOandPars.count(qHPO) == 0)
        {
            // base case: reached the bottom of the tree
            return;
        }
        cout << " " << (*it).second[j];
    }
    cout << "" << endl;
}


void PhenotypeIntegration::queryTermForKids(string qHPO)
{
    // cout << endl << "Hello queryTermForKids" << endl;

    // Fast queries of node and returns parents:
    // cout << "Query (Node): " << qHPO << endl;
    // cout << "Children: " << endl;
    map<string, vector<string>>::iterator it;
    it = mHPOandKids.find(qHPO);
    // cout << (*it).first << endl;
    for (unsigned int j = 0; j < (*it).second.size(); j++)
    {
        if (mHPOandKids.count(qHPO) == 0)
        {
            // base case: reached the bottom of the tree
            return;
        }
        cout << " " << (*it).second[j];
    }
    cout << "" << endl;
}


// Add a function that assosciated to the GeneName, the DiseaseOMIM(s), and the HPOs associated to the DiseaseOMIMs. These are the columns 2 and 9 of the HPO file "genes_to_phenotype.txt" (warning about release 20230405 has 4 columns, before and after there are correctly 9 columns).
void PhenotypeIntegration::getAssHPOs(string genNam)
{
    // cout << "Hello getAssHPOs" << endl;

    map<string, vector<vector<string>>>::iterator it;
    it = mGeneNamDiseHPOs.find(genNam);

    if (it != mGeneNamDiseHPOs.end()) {
        vvGenNamHPOsDiseaseID = (*it).second;

    // // // Printing the 2d vector, 1st D = DiseaseID ([0,]), 2nd D = each HPO connected to a set linked to one of the DiseaseID
    // for (unsigned int i = 0; i < (*it).second.size(); i++)
    // {
    // for (unsigned int ii = 0; ii < (*it).second[i].size(); ii++)
    // {
    //     cout<<" "<<(*it).second[i][ii];
    //     cout<<endl;
    // }
    //     cout<<endl;
    // }
    // cout<<endl;
    }
    else{ // It can be that the gene name in gene variant table has no correspondence in HPO file "genes_to_phenotype.txt", then assign an empty vector with no HPOs.
        cerr<<"Gene name not found in HPO: "<<genNam<<endl;
        vector<vector<string>> emptyVector;
        vvGenNamHPOsDiseaseID = emptyVector;
    }

}

// recursive function that walks UP till the root
void PhenotypeIntegration::treeWalkUp(string qHPO)
{
    if (mHPOandPars.count(qHPO) == 0)
    {
        // base case: reached the top of the tree
        return;
    }

    vector<string> &parents = mHPOandPars[qHPO];
    for (unsigned int j = 0; j < parents.size(); j++)
    {
        vPathUp.push_back(parents[j]);
        treeWalkUp(parents[j]);
    }
}


// recursive function that walks DOWN till the leaves
void PhenotypeIntegration::treeWalkDown(string qHPO)
{
    if (mHPOandKids.count(qHPO) == 0)
    {
        // base case: reached the bottom of the tree
        return;
    }

    vector<string> &children = mHPOandKids[qHPO];
    for (unsigned int j = 0; j < children.size(); j++)
    {
        vPathDown.push_back(children[j]);
        treeWalkDown(children[j]);
    }
}

string PhenotypeIntegration::findMICA(string qHPO1, string qHPO2)
{
    // cout << "Hello findMICA " << endl;

    // Compare two vectors containing the paths to the root,
    // and find the lowest common node
    // (i.e., the MICA = Most Informative Common Ancestor)

    string MICA;

    // Before calling recursive treeWalkUp: clear the path storage vector, and initiate it with the query for cases of direct ancestor MICA
    vPathUp.clear();
    vPathUp.push_back(qHPO1);
    treeWalkUp(qHPO1);
    vector<string> v1 = vPathUp;

    vPathUp.clear();
    vPathUp.push_back(qHPO2);
    treeWalkUp(qHPO2);
    vector<string> v2 = vPathUp;

    // easier if the order is reversed: from root (HP:1) to leaves
    reverse(v1.begin(), v1.end());
    reverse(v2.begin(), v2.end());

    // Check if a term is going to be present a second time, if yes prevent it (otherwise in some cases it can overwrite a MICA)
    for (unsigned int i = 0; i < v1.size(); i++)
    {
        for (unsigned int ii = i + 1; ii < v1.size(); ii++)
        {
            if (v1[i] == v1[ii])
            {
                v1.erase(v1.begin() + ii, v1.begin() + ii + 1);
            }
        }
    }

    for (unsigned int i = 0; i < v2.size(); i++)
    {
        for (unsigned int ii = i + 1; ii < v2.size(); ii++)
        {
            if (v2[i] == v2[ii])
            {
                v2.erase(v2.begin() + ii, v2.begin() + ii + 1);
            }
        }
    }

    // compares all the common elements from root to leaves, the last common one is the MICA
    for (unsigned int i = 0; i < v1.size(); i++)
    {
        for (unsigned int ii = 0; ii < v2.size(); ii++)
        {
            if (v1[i] == v2[ii]) // if in the common trunc...
            {
                MICA = v1[i];
                break;
            }
        }
    }
    return MICA;
}


void PhenotypeIntegration::uniqNodAllDesc()
{
    vPathDownSortUniq = vPathDown;
    // sort and unique
    sort(vPathDownSortUniq.begin(), vPathDownSortUniq.end());
    vPathDownSortUniq.erase(unique(vPathDownSortUniq.begin(), vPathDownSortUniq.end()), vPathDownSortUniq.end()); // erase the non-unique

    // for (unsigned int i = 0; i < vPathDownSortUniq.size(); i++)
    // {
    //   cout << " " << vPathDownSortUniq[i] << endl;
    // }
    // cout << "" << endl;

    // cout << "Number of nodes = "<< vPathDownSortUniq.size() << endl;
}

// Implementation of different combination functions all calculated here based on a similarity matrix, from which the different methods extract different things (max, max/average column/row, etc.)
void PhenotypeIntegration::calcSemSimSet(vector<string> HPOSet1, vector<string> HPOSet_2)
{
    // cout << "Hello calcSemSimSet" << endl;

    // HPOSet1: User terms (patient-annotated)
    // HPOSet_2: Disease terms (gene-annotated)

    // Define a similarity matrix user-terms x disease terms for both measures, later make all the operations like avg/max on columns/rows
    unsigned int rows = HPOSet1.size();
    unsigned int cols = HPOSet_2.size();
    vector<vector<double>> matrix_Resnik(rows, vector<double>(cols));
    vector<vector<double>> matrix_ERIC(rows, vector<double>(cols));

    for (unsigned int i = 0; i < HPOSet1.size(); i++)
    {
        string MICA;

        for (unsigned int ii = 0; ii < HPOSet_2.size(); ii++)
        {
            // here keep only the max for each HPOSet1

            MICA = findMICA(HPOSet1[i], HPOSet_2[ii]);

            // IC(MICA) = ResnikSim
            double ResnikSim = get_InformationContent(MICA); // Remember: Resnik(HPO1,HPO2) = IC( MICA(HPO1,HPO2) )

            // IC(t1)
            double IC1 = get_InformationContent(HPOSet1[i]);

            // IC(t2)
            double IC2 = get_InformationContent(HPOSet_2[ii]);

            double cutoff = min(IC1, IC2);
            double ERICSim = 0;
            if (2 * ResnikSim > cutoff)
            {
                ERICSim = 2 * ResnikSim - cutoff;
            }

            // cout<<HPOSet1[i]<<" "<<HPOSet_2[ii]<<endl;
            // cout << "Semantic Similarity ("<<HPOSet1[i]<<" - "<<HPOSet_2[ii]<<"). ";
            // cout<<"MICA: "<<MICA<<endl;
            // cout<<""<<endl;
            // cout << "IC1("<<HPOSet1[i]<<")="<<IC1<<", ";
            // cout << "IC2("<<HPOSet_2[ii]<<")="<<IC2<<". ";
            // cout<<""<<endl;
            // cout << "ResnikSim = " << ResnikSim <<". ";
            // cout<<""<<endl;
            // cout << "ERICSim = " << ERICSim <<endl;

            matrix_Resnik[i][ii] = ResnikSim;
            matrix_ERIC[i][ii] = ERICSim;
        }
    }

    // // display:
    // for (unsigned int i = 0; i < rows; i++) {
    //     for (unsigned int ii = 0; ii < cols; ii++) {
    //         cout << matrix_Resnik[i][ii] << " ";
    //     }
    //     cout << endl;
    // }

    // Iterate over the columns
    vector<double> v_col_sums_Resnik(cols, 0);
    vector<double> v_col_max_Resnik(cols, 0);
    for (unsigned int ii = 0; ii < cols; ii++)
    {
        double sum_col_Resnik = 0;
        double max_col_Resnik = 0;
        // Sum up the elements in the current column
        for (unsigned int i = 0; i < rows; i++)
        {
            sum_col_Resnik += matrix_Resnik[i][ii];
            if (max_col_Resnik < matrix_Resnik[i][ii])
            {
                max_col_Resnik = matrix_Resnik[i][ii];
            }
        }
        v_col_sums_Resnik[ii] = sum_col_Resnik;
        v_col_max_Resnik[ii] = max_col_Resnik;
    }
    // // display:
    // for (unsigned int ii = 0; ii < cols; ii++) {
    //     cout << "Sum of column " << ii << ": " << v_col_sums_Resnik[ii] << endl;
    //     cout << "Max of column " << ii << ": " << v_col_max_Resnik[ii] << endl;
    // }

    vector<double> v_col_sums_ERIC(cols, 0);
    vector<double> v_col_max_ERIC(cols, 0);
    for (unsigned int ii = 0; ii < cols; ii++)
    {
        double sum_col_ERIC = 0;
        double max_col_ERIC = 0;
        // Sum up the elements in the current column
        for (unsigned int i = 0; i < rows; i++)
        {
            sum_col_ERIC += matrix_ERIC[i][ii];
            if (max_col_ERIC < matrix_ERIC[i][ii])
            {
                max_col_ERIC = matrix_ERIC[i][ii];
            }
        }
        v_col_sums_ERIC[ii] = sum_col_ERIC;
        v_col_max_ERIC[ii] = max_col_ERIC;
    }
    // // display:
    // for (unsigned int ii = 0; ii < cols; ii++) {
    //     cout << "Sum of column " << ii << ": " << v_col_sums_ERIC[ii] << endl;
    //     cout << "Max of column " << ii << ": " << v_col_max_ERIC[ii] << endl;
    // }

    // Iterate over the rows
    vector<double> v_row_sums_Resnik(rows, 0);
    vector<double> v_row_max_Resnik(rows, 0);
    for (unsigned i = 0; i < rows; i++)
    {
        double sum_row_Resnik = 0;
        double max_row_Resnik = 0;
        // Sum up the elements in the current row
        for (unsigned int ii = 0; ii < cols; ii++)
        {
            sum_row_Resnik += matrix_Resnik[i][ii];
            if (max_row_Resnik < matrix_Resnik[i][ii])
            {
                max_row_Resnik = matrix_Resnik[i][ii];
            }
        }
        v_row_sums_Resnik[i] = sum_row_Resnik;
        v_row_max_Resnik[i] = max_row_Resnik;
    }
    // // display:
    // for (unsigned int i = 0; i < rows; i++) {
    //     cout << "Sum of row " << i << ": " << v_row_sums_Resnik[i] << endl;
    //     cout << "Max of row " << i << ": " << v_row_max_Resnik[i] << endl;
    // }

    vector<double> v_row_sums_ERIC(rows, 0);
    vector<double> v_row_max_ERIC(rows, 0);
    for (unsigned i = 0; i < rows; i++)
    {
        double sum_row_ERIC = 0;
        double max_row_ERIC = 0;
        // Sum up the elements in the current row
        for (unsigned int ii = 0; ii < cols; ii++)
        {
            sum_row_ERIC += matrix_ERIC[i][ii];
            if (max_row_ERIC < matrix_ERIC[i][ii])
            {
                max_row_ERIC = matrix_ERIC[i][ii];
            }
        }
        v_row_sums_ERIC[i] = sum_row_ERIC;
        v_row_max_ERIC[i] = max_row_ERIC;
    }
    // // display:
    // for (unsigned int i = 0; i < rows; i++) {
    //     cout << "Sum of row " << i << ": " << v_row_sums_ERIC[i] << endl;
    //     cout << "Max of row " << i << ": " << v_row_max_ERIC[i] << endl;
    // }

    // Less complicated combining scores:

    // 1. Max
    double max_overall_matrix_Resnik = 0;
    double max_overall_matrix_ERIC = 0;
    double sum_overall_Resnik = 0;
    double sum_overall_ERIC = 0;
    for (unsigned i = 0; i < rows; i++)
    {
        for (unsigned ii = 0; ii < cols; ii++)
        {
            sum_overall_Resnik += matrix_Resnik[i][ii];
            sum_overall_ERIC += matrix_ERIC[i][ii];
            // Update the maximum value if a larger value is found
            if (matrix_Resnik[i][ii] > max_overall_matrix_Resnik)
            {
                max_overall_matrix_Resnik = matrix_Resnik[i][ii];
            }
            if (matrix_ERIC[i][ii] > max_overall_matrix_ERIC)
            {
                max_overall_matrix_ERIC = matrix_ERIC[i][ii];
            }
        }
    }
    // cout<<endl<<"Max (Resnik): "<<max_overall_matrix_Resnik<<endl;
    // cout<<endl<<"Max (ERIC): "<<max_overall_matrix_ERIC<<endl;

    // Private variables
    maxRes = max_overall_matrix_Resnik;
    maxEric = max_overall_matrix_ERIC;

    // 2. Mean
    // cout<<endl<<"Mean (Resnik): "<<sum_overall_Resnik/(cols*rows)<<endl;
    // cout<<endl<<"Mean (ERIC): "<<sum_overall_ERIC/(cols*rows)<<endl;

    // Private variables
    meanRes = sum_overall_Resnik / (cols * rows);
    meanEric = sum_overall_ERIC / (cols * rows);

    // More complicated combining scores:

    double sum_rows_Resnik = 0;
    for (unsigned i = 0; i < rows; i++)
    {
        sum_rows_Resnik += v_row_max_Resnik[i];
    }
    double sum_cols_Resnik = 0;
    for (unsigned ii = 0; ii < cols; ii++)
    {
        sum_cols_Resnik += v_col_max_Resnik[ii];
    }

    double sum_rows_ERIC = 0;
    for (unsigned i = 0; i < rows; i++)
    {
        sum_rows_ERIC += v_row_max_ERIC[i];
    }
    double sum_cols_ERIC = 0;
    for (unsigned ii = 0; ii < cols; ii++)
    {
        sum_cols_ERIC += v_col_max_ERIC[ii];
    }

    // 3. BMA (Best Matching Average): Average of all maximum similarities on each row and column of the similarity matrix S.
    // Private variables
    BMARes = (sum_rows_Resnik + sum_cols_Resnik) / (rows + cols);
    BMAEric = (sum_rows_ERIC + sum_cols_ERIC) / (rows + cols);

    // cout<<endl<<"BMA (Resnik): "<<BMARes<<endl;
    // cout<<endl<<"BMA (ERIC): "<<BMAEric<<endl;

    // 4. BMS (Best Matching Sum): This should be as BMA but without denominator, however in ERIC/XRARE paper they just define as the first half, here is how they define it: "Our main question here is, how much the patient's phenotype set HPOSet1 (User terms) can be explained by the target gene corresponding to HPOSet_2 (Disease terms). Therefore, the similarity score between sets was defined as the sum of the best match (max) from HPOSet_2 for each phenotype term in HPOSet1." // (this is equivalent to my previous implementation):
    // Private variables
    BMSRes = (sum_rows_Resnik);
    BMSEric = (sum_rows_ERIC);

    // cout<<endl<<"BMS (Resnik) as in ERIC/XRARE: "<<BMSRes<<endl;
    // cout<<endl<<"BMS (ERIC) as in ERIC/XRARE: "<<BMSEric<<endl;

    // 5. FunSimAvg
    // Private variables
    FunSimAvRes = ((sum_rows_Resnik / rows) + (sum_cols_Resnik / cols)) / 2;
    FunSimAvEric = ((sum_rows_ERIC / rows) + (sum_cols_ERIC / cols)) / 2;
    // cout<<endl<<"FunSimAvg (Resnik): "<<FunSimAvRes<<endl;
    // cout<<endl<<"FunSimAvg (ERIC): "<<FunSimAvEric<<endl;

    // 6. FunSimMax
    FunSimMaxRes = max((sum_rows_Resnik / rows), (sum_cols_Resnik / cols));
    FunSimMaxEric = max((sum_rows_ERIC / rows), (sum_cols_ERIC / cols));
    // cout<<endl<<"FunSimMax (Resnik): "<<FunSimMaxRes<<endl;
    // cout<<endl<<"FunSimMax (ERIC): "<<FunSimMaxEric<<endl;
}

void PhenotypeIntegration::buildOut(string UTgenNam)
{
    ofstream outfile_combined_score_ERIC_ACMG;
    outfile_combined_score_ERIC_ACMG.open ("output_col_combined_score_ERIC_ACMG.dat");
    outfile_combined_score_ERIC_ACMG << "Gene-associated phenotypes" << endl; // header

    // cout << "Hello buildOut" << endl;

    // Particular way to make it accessible for unit tests: if UTgenNam is passed then the user input is replaces with single gene to be tested for HPO terms and Disease IDs extracted from "genes_to_phenotype.txt".
    // It the unit test sends over testing UTgenNam (& UnitTest_UserInput through a dedicated fucntion "set_buildOut_unit_test_inpUserHPOs()"), then those are used instead of conventional file input.
    if(UTgenNam!=""){
    // cout << "Hello Unit Test buildOut " <<UTgenNam<< endl;
    vInpVarTabGenName.clear();
    vInpVarTabScore.clear();
    vInpUserHPOs.clear();
    vInpVarTabGenName.push_back(UTgenNam);
    vInpVarTabScore.push_back(0); // for unit test assume ACMG score = 0 
    vInpUserHPOs=vUTInpUserHPOs;
    }


    for (unsigned int i = 0; i < vInpVarTabGenName.size(); i++) // screen over all GeneNames in variant table provided
    {
        string GeneNameInput = vInpVarTabGenName[i];
        double ACMGBayScore = vInpVarTabScore[i]; 
        set_getAssHPOs(GeneNameInput);
        vector<vector<string>> _vvGenNamHPOsDiseaseID = get_getAssHPOs();

        vector<tuple<double, string, string>> v_output;
        // vector < tuple<double, string, string, string> > v_output;
        for (unsigned int i = 0; i < _vvGenNamHPOsDiseaseID.size(); i++)
        {
            // Put in a vector only the HPO terms, leaving out the 3 final elements (GeneName, DiseaseID, Disease)
           vector<string> HPOSet; 
           copy(_vvGenNamHPOsDiseaseID[i].begin(), _vvGenNamHPOsDiseaseID[i].end() - 3, back_inserter(HPOSet));

            vector<double> _SemSim = get_calcSemSimSet(vInpUserHPOs, HPOSet);

            string DiseaseID = _vvGenNamHPOsDiseaseID[i][_vvGenNamHPOsDiseaseID[i].size() - 2];
            string DiseaseIDname = _vvGenNamHPOsDiseaseID[i][_vvGenNamHPOsDiseaseID[i].size() - 1];
            replace(DiseaseIDname.begin(), DiseaseIDname.end(), '_', ' ');

            // Make combination of ERIC + ACMG
            double combined_score_ERIC_ACMG = _SemSim[5] + ACMGBayScore; 
            
            v_output.push_back(make_tuple(combined_score_ERIC_ACMG, DiseaseID, DiseaseIDname)); // or if necessary to also add the (long) HPO list to the output (more for testing than for prod) (mod also below)
            // // --------
            // string separator = ", ";
            // string HPOlist = accumulate(vvGenNamHPOsDiseaseID[i].begin(), vvGenNamHPOsDiseaseID[i].end()-3, string(), [&separator](const string& a, const string& b) {return a.empty() ? b : a + separator + b;});
            // v_output.push_back(make_tuple(_SemSim[5],DiseaseID,DiseaseIDname,HPOlist));
            // // --------END
        }

        // Sort by decreasing SemSim score
        sort(v_output.begin(), v_output.end());
        reverse(v_output.begin(), v_output.end());



        outGenNamScoreDisID.clear();
        outGenNamScoreDisID = GeneNameInput + " { ";
        for (unsigned int i = 0; i < v_output.size(); i++)
        {
            // outfile_combined_score_ERIC_ACMG << get<0>(tuple) << " : " << get<1>(tuple) << " : " << get<2>(tuple) << endl;
            double score_ERIC_ACMG = get<0>(v_output[i]);
            string DiseaseID = get<1>(v_output[i]);
            string DiseaseIDname = get<2>(v_output[i]);
            outGenNamScoreDisID = outGenNamScoreDisID + to_string(score_ERIC_ACMG);
            outGenNamScoreDisID = outGenNamScoreDisID + " [";
            outGenNamScoreDisID = outGenNamScoreDisID + DiseaseID + " (";
            outGenNamScoreDisID = outGenNamScoreDisID + DiseaseIDname;


            // // also add the (long) HPO list
            // string HPOlist=get<3>(v_output[i]);
            // outGenNamScoreDisID=outGenNamScoreDisID+" (";
            // outGenNamScoreDisID=outGenNamScoreDisID+HPOlist ;
            // outGenNamScoreDisID=outGenNamScoreDisID+") ";

            outGenNamScoreDisID = outGenNamScoreDisID + ")] ";

            if (i != (_vvGenNamHPOsDiseaseID.size() - 1))
            {
                outGenNamScoreDisID = outGenNamScoreDisID + "| ";
            }
        }
        outGenNamScoreDisID = outGenNamScoreDisID + "}";

        // outfile_combined_score_ERIC_ACMG << vInpVarTabScore[i] << " ";
        outfile_combined_score_ERIC_ACMG << outGenNamScoreDisID << endl;
    }
    outfile_combined_score_ERIC_ACMG.close();
}