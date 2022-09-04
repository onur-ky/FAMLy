mkdir data
mkdir data/cluster
mkdir data/anno
mkdir data/metadata
mkdir data/fasta
mkdir data/weights
mkdir data/eval

mkdir data/train
mkdir data/train/comp_train
mkdir data/train/groups_train
mkdir data/train/comp_test
mkdir data/train/groups_test

echo "Enter path to specieswise group data: "
read group_data_path
echo "Enter complementary species KEGG code: "
read comp_species
echo "Enter path to whole proteome of the complementary species: "
read comp_species_whole_proteome_path
echo "Enter path to the evaluation proteome: "
read eval_proteome_path


echo "Ordering data groupwise..."
python scripts/order_groupwise.py --path $group_data_path

echo "Annotating groupwise proteomes..."
for fpath in `ls ./data/fasta/*.fasta`
do
  fas.doAnno -i $fpath -o ./data/anno
done

rm -r ./data/anno/tmp/

echo "Extracting metadata..."
python scripts/preprocess/extract_metadata.py

echo "Clustering input files..."
for fpath in `ls ./data/fasta/*.fasta`
do
  fname="$(basename -- $fpath)"
  #echo $fname
  ./scripts/preprocess/cdhit/cd-hit -i $fpath -o "./data/cluster/$fname"
  rm "./data/cluster/$fname.clstr"
done

echo "Removing complementary data..."
python scripts/preprocess/remove_comp_data.py --remove 'sce'  # species code here

echo "Populating training data..."
python scripts/preprocess/populate_training_data.py --comp_gdata_path "$group_data_path$comp_species.fasta" --comp_whole_proteome_path $comp_species_whole_proteome_path

echo "Processing evaluation data..."
python scripts/preprocess/process_eval_data.py

echo "Preprocessing is complete."