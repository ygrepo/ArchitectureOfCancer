
# Check if a gene name is provided
if [ -z "$1" ]; then
    echo "Usage: $0 [gene name]"
    exit 1
fi

# Assign the first argument to a variable
GENE_NAME=$1

echo "GENE:" $GENE_NAME
mkdir -p "/hpc/users/greaty01/github/ArchitectureOfCancer/data/genes/${GENE_NAME}"
#find .. -type f -not -path "./qsub_results/*" | xargs grep -Rwl -e "$GENE_NAME" > qsub_results/"${GENE_NAME}_files_search.txt"
#grep -Rwl /sc/arion/projects/DiseaseGeneCell/Huang_lab_data/Annotated_Genebass/ -e "$GENE_NAME" > ./"${GENE_NAME}_files_search.txt"
find /sc/arion/projects/DiseaseGeneCell/Huang_lab_data/Annotated_Genebass/ -type f -name "*annotated.tsv" -exec grep -l -w -e "$GENE_NAME" {} + > ~/github/ArchitectureOfCancer/data/genes/${GENE_NAME}/"${GENE_NAME}_files_search.txt"
