from subprocess import check_call

# http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgTfbsUniform/wgEncodeAwgTfbsBroadGm12878CtcfUniPk.narrowPeak.gz
label_index = 126
# label_index = 725
#label_index = 0

# check_call(['python3', 'gkm_datagen/split_pos_neg_for_gkm.py', 'chr1', 'train', str(label_index)])
# check_call(['python3', 'gkm_datagen/split_pos_neg_for_gkm.py', 'chr2', 'train', str(label_index)])
# check_call(['python3', 'gkm_datagen/split_pos_neg_for_gkm.py', 'chr3', 'train', str(label_index)])
# check_call(['python3', 'gkm_datagen/split_pos_neg_for_gkm.py', 'chr4', 'train', str(label_index)])
# check_call(['python3', 'gkm_datagen/split_pos_neg_for_gkm.py', 'chr5', 'train', str(label_index)])
# check_call(['python3', 'gkm_datagen/split_pos_neg_for_gkm.py', 'chr6', 'train', str(label_index)])
# check_call(['python3', 'gkm_datagen/split_pos_neg_for_gkm.py', 'chr7', 'train', str(label_index)])
# check_call(['python3', 'gkm_datagen/split_pos_neg_for_gkm.py', 'chr10', 'train', str(label_index)])
# check_call(['python3', 'gkm_datagen/split_pos_neg_for_gkm.py', 'chr11', 'train', str(label_index)])
# check_call(['python3', 'gkm_datagen/split_pos_neg_for_gkm.py', 'chr12', 'train', str(label_index)])
# check_call(['python3', 'gkm_datagen/split_pos_neg_for_gkm.py', 'chr13', 'train', str(label_index)])
# check_call(['python3', 'gkm_datagen/split_pos_neg_for_gkm.py', 'chr14', 'train', str(label_index)])
# check_call(['python3', 'gkm_datagen/split_pos_neg_for_gkm.py', 'chr15', 'train', str(label_index)])
# check_call(['python3', 'gkm_datagen/split_pos_neg_for_gkm.py', 'chr16', 'train', str(label_index)])
# check_call(['python3', 'gkm_datagen/split_pos_neg_for_gkm.py', 'chr17', 'train', str(label_index)])
# check_call(['python3', 'gkm_datagen/split_pos_neg_for_gkm.py', 'chr18', 'train', str(label_index)])
# check_call(['python3', 'gkm_datagen/split_pos_neg_for_gkm.py', 'chr19', 'train', str(label_index)])
# check_call(['python3', 'gkm_datagen/split_pos_neg_for_gkm.py', 'chr20', 'train', str(label_index)])
# check_call(['python3', 'gkm_datagen/split_pos_neg_for_gkm.py', 'chr21', 'train', str(label_index)])
# check_call(['python3', 'gkm_datagen/split_pos_neg_for_gkm.py', 'chr22', 'train', str(label_index)])

check_call(['python3', 'gkm_datagen/split_pos_neg_for_gkm.py', 'chr7', 'valid', str(label_index)])

check_call(['python3', 'gkm_datagen/split_pos_neg_for_gkm.py', 'chr8', 'test', str(label_index)])
check_call(['python3', 'gkm_datagen/split_pos_neg_for_gkm.py', 'chr9', 'test', str(label_index)])
