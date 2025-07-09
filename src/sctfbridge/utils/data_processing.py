import sys
import scanpy as sc
import episcanpy.api as epi
import numpy as np
import pandas as pd
import torch
import random
import pybedtools
import os
import anndata
import anndata as ad
from scipy.sparse import csr_matrix


def adata_multiomics_processing(adata_list: list,
                                output_path: str,
                                TF_list,
                                top_gene_num: int = 3000,
                                peak_percent: float = 0.01,
                                TF_keys=None
                                ):
    if TF_keys is None:
        TF_keys = ['ETS1', 'CTCF', 'STAT1', 'IRF1', 'RUNX1', 'SPI1']

    os.makedirs(f'{output_path}/', exist_ok=True)
    output_path = [f'{output_path}/RNA_filter.h5ad',
                   f'{output_path}/ATAC_filter.h5ad',
                   f'{output_path}/TF_filter.h5ad',
                   f'{output_path}/batch_info.csv']

    rna_adata = adata_list[0]

    # Find the indices of the TFs in rna_adata
    tf_indices = rna_adata.var_names.isin(TF_list)

    # Extract the TF data
    tf_data = rna_adata[:, tf_indices].copy()
    tf_mask = np.isinf(tf_data.X.data)
    tf_data.X.data[tf_mask] = 0

    tf_data = rna_adata[:, tf_indices].copy()

    # Step 2: Check for 'inf' values in the sparse matrix (using .data, which stores the non-zero values)
    tf_mask = np.isinf(tf_data.X.data)

    # Step 3: Replace 'inf' values with 0 in the sparse matrix
    tf_data.X.data[tf_mask] = 0

    sc.pp.normalize_total(tf_data, target_sum=1e4)
    sc.pp.log1p(tf_data)
    sc.pp.highly_variable_genes(tf_data, n_top_genes=128, subset=True)
    # tf_data.var.reset_index(drop=True, inplace=True)

    n = len(TF_keys)
    # 替换最后 n 个特征名称为 TF_list 中的对应值
    tf_data.var_names.values[-n:] = TF_keys

    # Find indices of genes that are NOT in TFName
    non_tf_indices = ~rna_adata.var_names.isin(TF_list)
    rna_adata = rna_adata[:, non_tf_indices]

    # 从adata.var中提取特征类型为ATAC的数据
    atac_adata = adata_list[1]

    # filtering rna data
    sc.pp.filter_cells(rna_adata, min_genes=200)  # 至少有200个基因在每个细胞中表达
    sc.pp.filter_genes(rna_adata, min_cells=1000)  # 每个细胞中至少有1000个UMI（唯一分子标记）

    mask = np.isinf(rna_adata.X.data)
    rna_adata.X.data[mask] = 0  # 或者你想要的其他替换值
    rna_adata.X.data[np.isnan(rna_adata.X.data)] = 0  # 填充NaN为0

    sc.pp.normalize_total(rna_adata, target_sum=1e4)
    sc.pp.log1p(rna_adata)
    sc.pp.highly_variable_genes(rna_adata, n_top_genes=top_gene_num, subset=True)

    # filtering atac data
    epi.pp.binarize(atac_adata)
    epi.pp.filter_features(atac_adata, min_cells=np.ceil(peak_percent * atac_adata.shape[0]))

    # Select the same cells
    shared_cells = np.intersect1d(rna_adata.obs_names, atac_adata.obs_names)
    rna_adata = rna_adata[shared_cells]
    atac_adata = atac_adata[shared_cells]

    # 保存批次信息
    if 'batch' in adata_list[0].obs.columns:
        batches_info = rna_adata.obs['batch'].values.tolist()
    else:
        batches_info = [0 for _ in range(rna_adata.X.shape[0])]

    one_hot_encoded_batches = pd.get_dummies(batches_info, prefix='batch', dtype=float)

    batch_info_df = pd.DataFrame({
        'batch': batches_info
    })

    batch_info_df['batch_encoded'] = one_hot_encoded_batches.values.tolist()

    # if '_index' in rna_adata.raw.var.columns:
    #     # Rename '_index' to something else
    #     rna_adata.raw.var.rename(columns={'_index': 'index'}, inplace=True)

    rna_adata.write(output_path[0])
    atac_adata.write(output_path[1])

    # if '_index' in tf_data.raw.var.columns:
    #     # Rename '_index' to something else
    #     tf_data.raw.var.rename(columns={'_index': 'index'}, inplace=True)

    tf_data.write(output_path[2])
    batch_info_df.to_csv(output_path[3], index=True)

    print('finish data pre-processing')


def multiomics_processing(adata_list: list,
                        output_path: str,
                        TF_list,  # Note: TF_list is now used to filter out TFs from the main rna_adata
                        top_gene_num: int = 3000,
                        peak_percent: float = 0.01,
                          TF_keys=None
                                ):
    """
    Processes multi-omic data based on the new logic:
    1. Extracts a broad set of TFs based on TF_list.
    2. Identifies the top 128 highly variable TFs from this pool.
    3. Ensures that all TFs from TF_keys are retained in the final TF set,
       regardless of their variability.
    """
    if TF_keys is None:
        TF_keys = ['ETS1', 'CTCF', 'STAT1', 'IRF1', 'RUNX1', 'SPI1']

    os.makedirs(f'{output_path}/', exist_ok=True)
    output_files = [f'{output_path}/RNA_filter.h5ad',
                    f'{output_path}/ATAC_filter.h5ad',
                    f'{output_path}/TF_filter.h5ad',
                    f'{output_path}/batch_info.csv']

    rna_adata = adata_list[0]

    # --- Start of Revised Logic for TF Data Creation ---

    # Step 1: 从 rna_adata 中提取 TF_list 定义的广泛TF，构建“TF候选池”
    tf_in_rna = [tf for tf in TF_list if tf in rna_adata.var_names]
    print(f"Found {len(tf_in_rna)} TFs from TF_list in the original RNA data.")

    if not tf_in_rna:
        raise ValueError("No TFs from TF_list were found in the rna_adata. Cannot proceed.")

    tf_adata = rna_adata[:, tf_in_rna].copy()

    # 对TF数据进行初步清理和预处理
    tf_adata.X.data[np.isinf(tf_adata.X.data)] = 0
    tf_adata.X.eliminate_zeros()
    sc.pp.normalize_total(tf_adata, target_sum=1e4)
    sc.pp.log1p(tf_adata)

    # Step 2: 在“TF候选池”中识别高变TF，但暂时不筛选 (subset=False)
    n_hvg_tf = 128
    # 确保请求的高变基因数量不超过总TF数
    if tf_adata.shape[1] < n_hvg_tf:
        n_hvg_tf = tf_adata.shape[1]
        print(f"Warning: Number of available TFs ({n_hvg_tf}) is less than 128. Using all available TFs as HVGs.")

    sc.pp.highly_variable_genes(tf_adata, n_top_genes=n_hvg_tf, subset=False)

    # 获取高变TF的列表
    highly_variable_tfs = tf_adata.var_names[tf_adata.var['highly_variable']].tolist()

    # Step 3: 合并高变TF和必须保留的TF_keys，生成最终要保留的TF列表
    # 首先，检查TF_keys中有哪些是在我们提取的tf_adata中确实存在的
    must_keep_tfs = [tf for tf in TF_keys if tf in tf_adata.var_names]
    missing_keys = [tf for tf in TF_keys if tf not in tf_adata.var_names]
    if missing_keys:
        print(f"Warning: The following essential TFs from TF_keys were not found in the initial TF_list pool: {missing_keys}")

    # 使用集合（set）来合并列表并自动去重
    final_tf_selection_set = set(highly_variable_tfs) | set(must_keep_tfs)
    final_tf_list = list(final_tf_selection_set)
    print(f"Final TF set includes {len(highly_variable_tfs)} HVG TFs and {len(must_keep_tfs)} essential TFs, for a total of {len(final_tf_list)} unique TFs.")

    # 根据最终的合并列表，对tf_adata进行切片，得到最终结果
    tf_data = tf_adata[:, final_tf_list].copy()

    # --- End of Revised Logic ---


    # --- 后续处理流程保持不变 ---
    # Find indices of all genes that are TFs (from the broader TF_list) and remove them from main RNA adata
    non_tf_indices = ~rna_adata.var_names.isin(TF_list)
    rna_adata = rna_adata[:, non_tf_indices]

    # Get ATAC data
    atac_adata = adata_list[1]

    # Filter RNA data
    sc.pp.filter_cells(rna_adata, min_genes=200)
    sc.pp.filter_genes(rna_adata, min_cells=1000)
    rna_adata.X.data[np.isinf(rna_adata.X.data)] = 0
    rna_adata.X.data[np.isnan(rna_adata.X.data)] = 0
    rna_adata.X.eliminate_zeros()
    sc.pp.normalize_total(rna_adata, target_sum=1e4)
    sc.pp.log1p(rna_adata)
    sc.pp.highly_variable_genes(rna_adata, n_top_genes=top_gene_num, subset=True)

    # Filter ATAC data
    epi.pp.binarize(atac_adata)
    epi.pp.filter_features(atac_adata, min_cells=np.ceil(peak_percent * atac_adata.shape[0]))

    # Find common cells across all modalities (RNA, ATAC, and TF)
    shared_cells = list(set(rna_adata.obs_names) & set(atac_adata.obs_names) & set(tf_data.obs_names))

    rna_adata = rna_adata[shared_cells, :].copy()
    atac_adata = atac_adata[shared_cells, :].copy()
    tf_data = tf_data[shared_cells, :].copy()

    # Save batch information
    if 'batch' in rna_adata.obs.columns:
        batches_info = rna_adata.obs['batch'].values
    else:
        batches_info = ['batch_0'] * rna_adata.shape[0]

    batch_info_df = pd.DataFrame({'batch': batches_info}, index=rna_adata.obs_names)

    # Save the processed data
    rna_adata.write(output_files[0])
    atac_adata.write(output_files[1])
    tf_data.write(output_files[2])
    batch_info_df.to_csv(output_files[3], index=True)

    print('Finish data pre-processing')



def five_fold_split_dataset(RNA_data, output_path):
    temp = [i for i in range(len(RNA_data.obs_names))]
    random.shuffle(temp)
    id_list = []
    test_count = int(0.2 * len(temp))
    validation_count = int(0.16 * len(temp))

    for i in range(5):
        test_id = temp[: test_count]
        validation_id = temp[test_count: test_count + validation_count]
        train_id = temp[test_count + validation_count:]
        temp.extend(test_id)
        temp = temp[test_count:]
        id_list.append([train_id, validation_id, test_id])
    df_fold = pd.DataFrame(id_list, columns=['train_id', 'validation_id', 'test_id'])
    df_fold.to_csv(output_path)


def preload_TF_binding(adata_path,
                       prior_path,
                       output_path):
    os.makedirs(output_path, exist_ok=True)
    # output_path = f'filter_data/{dataset_name}/TF_binding/'

    atac_data = anndata.read_h5ad(f'{adata_path}/ATAC_filter.h5ad')
    TF_data = anndata.read_h5ad(f'{adata_path}/TF_filter.h5ad')

    Match2 = pd.read_csv(prior_path + 'Match2.txt', sep='\t')
    Match2 = Match2.values
    motifWeight = pd.read_csv(prior_path + 'motifWeight.txt', index_col=0, sep='\t')

    TFName = TF_data.var.index.values
    Element_name = atac_data.var.index

    # 替换每个字符串中的第一个 '-' 为 ':' 只有BMMC用
    # Element_name = Element_name.str.replace('-', ':', 1)
    pd.DataFrame(Element_name).to_csv(f'{output_path}/Peaks.txt', header=None, index=None)
    extract_overlap_regions('hg38', prior_path, output_path, 'scTFBridge', f'{output_path}/Peaks.txt')
    print(Element_name)
    load_TFbinding(prior_path, motifWeight, Match2, TFName, Element_name, output_path)

    print('finish load TF_binding')


def load_motifbinding_chr(chrN, GRNdir, motifWeight, outdir):
    Motif_binding_temp = pd.read_csv(GRNdir + 'MotifTarget_Matrix_' + chrN + '.txt', sep='\t', index_col=0)
    REs = Motif_binding_temp.index
    march_hg19_Regrion = pd.read_csv(outdir + 'MotifTarget_hg19_hg38_' + chrN + '.txt', sep='\t', header=None)
    REoverlap = list(set(march_hg19_Regrion[1].values))
    Motif_binding_temp1 = Motif_binding_temp.loc[REoverlap]
    REs = Motif_binding_temp1.index
    Motif_binding_temp = np.zeros([march_hg19_Regrion.shape[0], Motif_binding_temp.shape[1]])
    Motif_binding_temp = Motif_binding_temp1.loc[march_hg19_Regrion[1].values].values
    Motif_binding_temp = pd.DataFrame(Motif_binding_temp, index=march_hg19_Regrion[0].values,
                                      columns=Motif_binding_temp1.columns)
    Motif_binding_temp1 = Motif_binding_temp.groupby(Motif_binding_temp.index).max()
    motifoverlap = list(set(Motif_binding_temp1.columns) & set(motifWeight.index))
    Motif_binding_temp1 = Motif_binding_temp1[motifoverlap]
    motifWeight = motifWeight.loc[Motif_binding_temp1.columns]
    Motif_binding = np.diag(1.0 / (motifWeight.T + 0.1)) * Motif_binding_temp1.values.T
    Motif_binding = np.log1p(Motif_binding)
    return Motif_binding_temp1


def load_TFbinding(GRNdir, motifWeight, Match2, TFName, Element_name, outdir):
    from tqdm import tqdm
    motif_binding = pd.DataFrame()
    chrall = ['chr' + str(i + 1) for i in range(22)]
    chrall.append('chrX')
    for chrN in tqdm(chrall):
        Motif_binding_temp1 = load_motifbinding_chr(chrN, GRNdir, motifWeight, outdir)
        motif_binding = pd.concat([motif_binding, Motif_binding_temp1], join='outer', axis=0)
    motif_binding = motif_binding.fillna(0)
    motif_binding = motif_binding.groupby(motif_binding.index).max()
    motifoverlap = list(set(motif_binding.columns) & set(motifWeight.index))
    Match2 = Match2[np.isin(Match2[:, 0], motifoverlap), :]
    TF_binding_temp = np.zeros((len(TFName), len(Element_name)))
    Motif_binding = np.zeros((motif_binding.shape[1], len(Element_name)))
    Element_name_idx = pd.DataFrame(range(len(Element_name)), index=Element_name)
    idx = Element_name_idx.loc[motif_binding.index][0].values
    Motif_binding = np.zeros((motif_binding.shape[1], len(Element_name)))
    Motif_binding[:, idx] = motif_binding.loc[Element_name[idx]].values.T
    Motif_binding = pd.DataFrame(Motif_binding, index=motif_binding.columns, columns=Element_name)
    Match2 = Match2[np.isin(Match2[:, 1], TFName), :]
    Motif_binding = Motif_binding.loc[Match2[:, 0]]
    Motif_binding.index = Match2[:, 1]
    TF_binding = Motif_binding.groupby(Motif_binding.index).sum()
    a = np.sum(TF_binding.values, axis=1)
    a[a == 0] = 1
    TF_binding_n = np.diag(1.0 / a) @ TF_binding.values
    TF_binding_n = pd.DataFrame(TF_binding_n.T, index=Element_name, columns=TF_binding.index)
    TF_binding = np.zeros((len(Element_name), len(TFName)))
    idx = np.isin(TFName, TF_binding_n.columns)
    TF_binding[:, idx] = TF_binding_n[TFName[idx]].values
    TF_binding = pd.DataFrame(TF_binding, index=Element_name, columns=TFName)
    TF_binding.to_csv(outdir + 'TF_binding.txt', sep='\t', index=None, header=None)


def extract_overlap_regions(genome, GRNdir, outdir, method, input_file_path):
    os.makedirs(outdir, exist_ok=True)
    input_file = input_file_path
    output_file = outdir + 'Region.bed'
    print(output_file)
    # Read the input file
    df = pd.read_csv(input_file, sep='\t', header=None)
    chromosomes = [item.split(':')[0] for item in df[0].values]
    # Drop the first row
    # Replace ':' and '-' with tabs
    df = df.replace({':': '\t', '-': '\t'}, regex=True)
    chrall = ['chr' + str(i + 1) for i in range(23)] + ['chrX']

    df = df[pd.DataFrame(chromosomes)[0].isin(chrall).values]
    print(df)

    df.to_csv(output_file, index=None, header=None)
    if method == 'scTFBridge':
        if genome == 'hg38':
            print(outdir + 'Region.bed')
            filepath = os.path.abspath(outdir + 'Region.bed')
            a = pybedtools.example_bedtool(filepath)
            bpath = os.path.abspath(GRNdir + 'hg38_hg19_pair.bed')
            b = pybedtools.example_bedtool(bpath)
            a_with_b = a.intersect(b, wa=True, wb=True)
            a_with_b.saveas(outdir + 'temp.bed')
            a_with_b = pd.read_csv(outdir + 'temp.bed', sep='\t', header=None)
            a_with_b[[6, 7, 8, 0, 1, 2]].to_csv(outdir + 'match_hg19_peak.bed', sep='\t', header=None, index=None)
        if genome == 'hg19':
            a = pybedtools.example_bedtool(outdir + 'Region.bed')
            b = pybedtools.example_bedtool(GRNdir + 'hg19_hg38_pair.bed')
            a_with_b = a.intersect(b, wa=True, wb=True)
            a_with_b.saveas(outdir + 'temp.bed')
            a_with_b = pd.read_csv(outdir + 'temp.bed', sep='\t', header=None)
            a_with_b[[6, 7, 8, 0, 1, 2]].to_csv(outdir + 'match_hg19_peak.bed', sep='\t', header=None, index=None)

        a = pybedtools.example_bedtool(os.path.abspath(outdir + 'match_hg19_peak.bed'))
        b = pybedtools.example_bedtool(os.path.abspath(GRNdir + 'RE_gene_corr_hg19.bed'))
        a_with_b = a.intersect(b, wa=True, wb=True)
        a_with_b.saveas(outdir + 'temp.bed')
        a_with_b = pd.read_csv(outdir + 'temp.bed', sep='\t', header=None)
        a_with_b = a_with_b[(a_with_b[1].values == a_with_b[7].values) & (a_with_b[2].values == a_with_b[8].values)]
        a_with_b_n = pd.DataFrame({
            'column1': a_with_b[0] + ':' + a_with_b[1].astype(str) + '-' + a_with_b[2].astype(str),
            'column2': a_with_b[3] + ':' + a_with_b[4].astype(str) + '-' + a_with_b[5].astype(str),
            'column3': a_with_b[9]})
        a_with_b_n = a_with_b_n.drop_duplicates()
        a_with_b_n.to_csv(outdir + 'hg19_Peak_hg19_gene_u.txt', sep='\t', header=None, index=None)
        chr_all = ['chr' + str(i + 1) for i in range(22)]
        chr_all.append('chrX')
        for chrtemp in chr_all:
            a = pybedtools.example_bedtool(os.path.abspath(outdir + 'match_hg19_peak.bed'))
            b = pybedtools.example_bedtool(os.path.abspath(GRNdir + 'MotifTarget_matrix_' + chrtemp + '.bed'))
            a_with_b = a.intersect(b, wa=True, wb=True)
            a_with_b.saveas(outdir + 'temp.bed')
            a_with_b = pd.read_csv(outdir + 'temp.bed', sep='\t', header=None)
            a_with_b = a_with_b[(a_with_b[1].values == a_with_b[7].values) & (a_with_b[2].values == a_with_b[8].values)]
            a_with_b_n = pd.DataFrame({
                'column1': a_with_b[3] + ':' + a_with_b[4].astype(str) + '-' + a_with_b[5].astype(str),
                'column2': a_with_b[6] + ':' + a_with_b[7].astype(str) + '-' + a_with_b[8].astype(str)})
            a_with_b_n = a_with_b_n.drop_duplicates()
            a_with_b_n.to_csv(outdir + 'MotifTarget_hg19_hg38_' + chrtemp + '.txt', sep='\t', header=None, index=None)
            a = pybedtools.example_bedtool(os.path.abspath(GRNdir + genome + '_Peaks_' + chrtemp + '.bed'))
            b = pybedtools.example_bedtool(os.path.abspath(outdir + 'Region.bed'))
            a_with_b = a.intersect(b, wa=True, wb=True)
            a_with_b.saveas(outdir + 'Region_overlap_' + chrtemp + '.bed')
    if method == 'baseline':
        chr_all = ['chr' + str(i + 1) for i in range(22)]
        chr_all.append('chrX')
        for chrtemp in chr_all:
            a = pybedtools.example_bedtool(GRNdir + genome + '_Peaks_' + chrtemp + '.bed')
            b = pybedtools.example_bedtool(outdir + 'Region.bed')
            a_with_b = a.intersect(b, wa=True, wb=True)
            a_with_b.saveas(outdir + 'Region_overlap_' + chrtemp + '.bed')


def merge_columns_in_bed_file(file_path, startcol):
    merged_values = []
    with open(file_path, 'r') as file:
        for line in file:
            columns = line.strip().split('\t')
            col1 = columns[-1 + startcol]
            col2 = columns[startcol]
            col3 = columns[1 + startcol]
            merged_value = f"{col1}:{col2}-{col3}"
            merged_values.append(merged_value)
    return merged_values


def merge_columns_in_bed_file2(file_path, startcol):
    merged_values = []
    with open(file_path, 'r') as file:
        for line in file:
            columns = line.strip().split('\t')
            col1 = columns[-1 + startcol]
            col2 = columns[startcol]
            col3 = columns[1 + startcol]
            merged_value = f"{col1}_{col2}_{col3}"
            merged_values.append(merged_value)
    return merged_values


def load_region(GRNdir, genome, chrN, outdir):
    O_overlap = merge_columns_in_bed_file(outdir + 'Region_overlap_' + chrN + '.bed', 1)
    N_overlap = merge_columns_in_bed_file(outdir + 'Region_overlap_' + chrN + '.bed', 4)
    O_overlap_u = list(set(O_overlap))
    N_overlap_u = list(set(N_overlap))
    #O_all=merge_columns_in_bed_file(GRNdir+'Peaks_'+chrN+'.bed',1)
    hg19_region = merge_columns_in_bed_file(GRNdir + 'hg19_Peaks_' + chrN + '.bed', 1)
    hg19_region = pd.DataFrame(range(len(hg19_region)), index=hg19_region)
    hg38_region = merge_columns_in_bed_file(GRNdir + 'hg38_Peaks_' + chrN + '.bed', 1)
    hg38_region = pd.DataFrame(range(len(hg38_region)), index=hg38_region)
    if genome == 'hg19':
        idx = hg19_region.loc[O_overlap_u][0].values
        O_overlap_u = hg38_region.index[idx].tolist()
        O_overlap_hg19_u = hg19_region.index[idx].tolist()
    if genome == 'hg38':
        idx = hg38_region.loc[O_overlap_u][0].values
        O_overlap_hg19_u = hg19_region.index[idx].tolist()
    return O_overlap, N_overlap, O_overlap_u, N_overlap_u, O_overlap_hg19_u


def load_RE_TG(GRNdir, chrN, O_overlap_u, O_overlap_hg19_u, O_overlap):
    #print('load prior RE-TG ...')
    from scipy.sparse import coo_matrix
    primary_s = pd.read_csv(GRNdir + 'Primary_RE_TG_' + chrN + '.txt', sep='\t')
    primary_s["RE"] = primary_s["RE"].apply(lambda x: x.split('_')[0] + ':' + x.split('_')[1] + '-' + x.split('_')[2])
    primary_s = primary_s[primary_s["RE"].isin(O_overlap_u)]
    TGset = primary_s["TG"].unique()
    REset = O_overlap_u
    # Create a dictionary mapping column names and row names to integer indices
    col_dict = {col: i for i, col in enumerate(TGset)}
    row_dict = {row: i for i, row in enumerate(REset)}
    # Map the column names and row names to integer indices in the DataFrame
    primary_s.loc[:, "col_index"] = primary_s["TG"].map(col_dict)
    primary_s.loc[:, "row_index"] = primary_s["RE"].map(row_dict)
    # Extract the column indices, row indices, and values from the DataFrame
    col_indices = primary_s["col_index"].tolist()
    row_indices = primary_s["row_index"].tolist()
    values = primary_s["score"].tolist()
    # Create the sparse matrix using coo_matrix
    sparse_S = coo_matrix((values, (row_indices, col_indices)))
    sparse_S.colnames = TGset
    sparse_S.rownames = REset
    array = sparse_S.toarray()
    O_overlap_u_df = pd.DataFrame(range(len(O_overlap_u)), index=O_overlap_u)
    hg19_38 = pd.DataFrame(O_overlap_u, index=O_overlap_hg19_u)
    array2 = np.zeros([len(O_overlap), array.shape[1]])
    index = O_overlap_u_df.loc[O_overlap][0].values
    array2 = array[index, :]
    array = pd.DataFrame(array2, index=O_overlap, columns=TGset)
    return array, TGset


def load_RE_TG_distance(GRNdir, chrN, O_overlap_hg19_u, O_overlap_u, O_overlap, TGoverlap):
    #print('load RE-TG distance for '+chrN+'...')
    from scipy.sparse import coo_matrix
    Dis = pd.read_csv(GRNdir + 'RE_TG_distance_' + chrN + '.txt', sep='\t', header=None)
    Dis.columns = ['RE', 'TG', 'dis']
    Dis["RE"] = Dis["RE"].apply(lambda x: x.split('_')[0] + ':' + x.split('_')[1] + '-' + x.split('_')[2])
    Dis = Dis[Dis["RE"].isin(O_overlap_hg19_u)]
    Dis = Dis[Dis['TG'].isin(TGoverlap)]
    col_dict = {col: i for i, col in enumerate(TGoverlap)}
    row_dict = {row: i for i, row in enumerate(O_overlap_hg19_u)}
    # Map the column names and row names to integer indices in the DataFrame
    Dis.loc[:, "col_index"] = Dis["TG"].map(col_dict)
    Dis.loc[:, "row_index"] = Dis["RE"].map(row_dict)
    col_indices = Dis["col_index"].tolist()
    row_indices = Dis["row_index"].tolist()
    values = Dis["dis"].tolist()
    # Create the sparse matrix using coo_matrix
    sparse_dis = coo_matrix((values, (row_indices, col_indices)), shape=(len(O_overlap_u), len(TGoverlap)))
    sparse_dis.colnames = TGoverlap
    sparse_dis.rownames = O_overlap_u
    sparse_dis = sparse_dis.tocsc()
    A = sparse_dis.multiply(1 / 25000)
    A.data += 0.5
    A.data = np.exp(-A.data)
    sparse_dis = A
    array = sparse_dis.toarray()
    O_overlap_u_df = pd.DataFrame(range(len(O_overlap_u)), index=O_overlap_u)
    hg19_38 = pd.DataFrame(O_overlap_u, index=O_overlap_hg19_u)
    array2 = np.zeros([len(O_overlap), array.shape[1]])
    index = O_overlap_u_df.loc[O_overlap][0].values
    array2 = array[index, :]
    array = pd.DataFrame(array2, index=O_overlap, columns=TGoverlap)
    return array
