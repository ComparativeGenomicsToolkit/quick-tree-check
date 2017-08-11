#!/bin/bash
# Re-generate the output TSV containing runtime & accuracy information.
set -o errexit
set -o nounset
set -o pipefail

DATA_DIR=/mnt/fasta
OUTPUT_FILE=output.tsv

python least_squares_fit.py --tsv testdata/tree.nh --files testdata/simCow.chr6 testdata/simDog.chr6 testdata/simHuman.chr6 testdata/simMouse.chr6 testdata/simRat.chr6 --labels simCow_chr6 simDog_chr6 simHuman_chr6 simMouse_chr6 simRat_chr6 > ${OUTPUT_FILE}
python least_squares_fit.py --tsv testdata/tree.nh --files testdata/simCow.chr6 testdata/simDog.chr6 testdata/simHuman.chr6 testdata/simMouse.chr6 testdata/simRat.chr6 --labels simCow_chr6 simDog_chr6 simHuman_chr6 simMouse_chr6 simRat_chr6 --method kmacs --noHeader >> ${OUTPUT_FILE}
python least_squares_fit.py --tsv testdata/tree.nh --files testdata/simCow.chr6 testdata/simDog.chr6 testdata/simHuman.chr6 testdata/simMouse.chr6 testdata/simRat.chr6 --labels simCow_chr6 simDog_chr6 simHuman_chr6 simMouse_chr6 simRat_chr6 --method spaced --noHeader >> ${OUTPUT_FILE}
python least_squares_fit.py --tsv ${DATA_DIR}/nematodes/tree.nh --files ${DATA_DIR}/nematodes/C_elegans.fa ${DATA_DIR}/nematodes/C_japonica.fa ${DATA_DIR}/nematodes/P_pacificus.fa --labels C_elegans C_japonica P_pacificus --noHeader >> ${OUTPUT_FILE}
python least_squares_fit.py --tsv ${DATA_DIR}/nematodes/tree.nh --files ${DATA_DIR}/nematodes/C_elegans.fa ${DATA_DIR}/nematodes/C_japonica.fa ${DATA_DIR}/nematodes/P_pacificus.fa --labels C_elegans C_japonica P_pacificus --method kmacs --noHeader >> ${OUTPUT_FILE}
python least_squares_fit.py --tsv ${DATA_DIR}/nematodes/tree.nh --files ${DATA_DIR}/nematodes/C_elegans.fa ${DATA_DIR}/nematodes/C_japonica.fa ${DATA_DIR}/nematodes/P_pacificus.fa --labels C_elegans C_japonica P_pacificus --method spaced --noHeader >> ${OUTPUT_FILE}
python least_squares_fit.py --tsv ${DATA_DIR}/apes/tree.nh --files ${DATA_DIR}/apes/hg38.fa ${DATA_DIR}/apes/panTro5.fa ${DATA_DIR}/apes/susie.fa ${DATA_DIR}/apes/ponAbe2.fa ${DATA_DIR}/apes/nomLeu3.fa --labels hg38 panTro5 susie ponAbe2 nomLeu3 --noHeader >> ${OUTPUT_FILE}
python least_squares_fit.py --tsv ${DATA_DIR}/apes/tree.nh --files ${DATA_DIR}/apes/hg38.fa ${DATA_DIR}/apes/panTro5.fa ${DATA_DIR}/apes/susie.fa ${DATA_DIR}/apes/ponAbe2.fa ${DATA_DIR}/apes/nomLeu3.fa --labels hg38 panTro5 susie ponAbe2 nomLeu3 --method kmacs --noHeader >> ${OUTPUT_FILE}
python least_squares_fit.py --tsv ${DATA_DIR}/apes/tree.nh --files ${DATA_DIR}/apes/hg38.fa ${DATA_DIR}/apes/panTro5.fa ${DATA_DIR}/apes/susie.fa ${DATA_DIR}/apes/ponAbe2.fa ${DATA_DIR}/apes/nomLeu3.fa --labels hg38 panTro5 susie ponAbe2 nomLeu3 --method spaced --noHeader >> ${OUTPUT_FILE}

# replace dataset column with friendlier names
sed -i 's^${DATA_DIR}/apes/tree.nh^apes^g' ${OUTPUT_FILE}
sed -i 's^${DATA_DIR}/nematodes/tree.nh^nematodes^g' ${OUTPUT_FILE}
sed -i 's^$testdata/tree.nh^simulatedMammals^g' ${OUTPUT_FILE}
