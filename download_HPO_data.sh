mkdir HPO
cd HPO

mkdir HPO_release20220414
cd HPO_release20220414
#wget https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2022-04-14/hp.obo
wget https://github.com/obophenotype/human-phenotype-ontology/raw/v2022-04-14/hp.obo
wget https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2022-04-14/genes_to_phenotype.txt
cd ..

mkdir HPO_release20230127
cd HPO_release20230127
wget https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2023-01-27/hp.obo
wget https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2023-01-27/genes_to_phenotype.txt
cd ..

cd ..
